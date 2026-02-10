#!/usr/bin/env Rscript

###############################################################################
# outlier_annotation.R
#
# Annotate outlier regions: extract sequences, BLAST per region against
# multiple DBs (from config file), map hits to GFF annotations (genes, GO).
# Output: overlapping genes with optional region metadata (Option A).
#
# Regions file (e.g. from identify_outliers outlier_regions*.csv) may contain
# extra columns (outlier_stat, n_windows, region_mean_coverage, etc.); when
# present they are preserved and included in the output so users can link
# genes back to region statistics.
#
# BLAST config: YAML or JSON with list of entries:
#   - db_path: path to BLAST DB (or FASTA for makeblastdb)
#   - name: label for this DB
#   - gff_path: path to GFF/GTF for mapping hits to genes
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Parse BLAST config (YAML or JSON)
read_blast_config <- function(path) {
  if (!file.exists(path)) stop("BLAST config file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yml", "yaml")) {
    if (!requireNamespace("yaml", quietly = TRUE))
      stop("Package 'yaml' required for YAML config. Install with: install.packages(\"yaml\")")
    cfg <- yaml::read_yaml(path)
  } else if (ext == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE))
      stop("Package 'jsonlite' required for JSON config. Install with: install.packages(\"jsonlite\")")
    cfg <- jsonlite::read_json(path)
  } else {
    stop("BLAST config must be YAML (.yml/.yaml) or JSON (.json)")
  }
  dbs <- cfg$databases %||% cfg
  if (!is.list(dbs)) dbs <- list(dbs)
  if (length(dbs) > 0 && !is.null(names(dbs)) && !"db_path" %in% names(dbs))
    dbs <- list(dbs)
  if (length(dbs) == 0) stop("BLAST config has no 'databases' or list of { db_path, name, gff_path }")
  dbs
}

# Read regions file (CSV/TSV from identify_outliers: chr, region_start, region_end, ...).
# Retains extra columns when present (e.g. outlier_stat, n_windows, region_mean_coverage,
# region_mean_mapping_quality, region_n_snps) for pass-through to output.
read_regions <- function(path) {
  delim <- if (grepl("\t", readLines(path, n = 1))) "\t" else ","
  d <- if (delim == "\t") read_tsv(path, show_col_types = FALSE) else read_csv(path, show_col_types = FALSE)
  d <- d %>% rename_all(tolower)
  start_col <- intersect(c("region_start", "start"), names(d))[1]
  end_col <- intersect(c("region_end", "end"), names(d))[1]
  chr_col <- intersect(c("chr", "chrom", "chromosome"), names(d))[1]
  if (is.na(chr_col) || is.na(start_col) || is.na(end_col))
    stop("Regions file must have chr/chrom, region_start/start, region_end/end")
  d <- d %>% rename(chr = !!chr_col, region_start = !!start_col, region_end = !!end_col)
  d$region_start <- as.numeric(d$region_start)
  d$region_end <- as.numeric(d$region_end)
  d %>% filter(!is.na(.data$chr), !is.na(.data$region_start), !is.na(.data$region_end))
}

# Write region FASTA using samtools faidx (reference must be indexed)
write_region_fasta <- function(regions, reference, out_fasta, samtools = "samtools") {
  dir.create(dirname(out_fasta), showWarnings = FALSE, recursive = TRUE)
  fai <- paste0(reference, ".fai")
  if (!file.exists(fai)) stop("Reference index not found: ", fai)
  unlink(out_fasta)
  for (i in seq_len(nrow(regions))) {
    r <- regions[i, ]
    reg <- paste0(r$chr, ":", r$region_start, "-", r$region_end)
    seq <- system2(samtools, c("faidx", reference, reg), stdout = TRUE, stderr = FALSE)
    cat(seq, sep = "\n", file = out_fasta, append = TRUE)
  }
  out_fasta
}

# Run BLAST (blastn) per region query; return path to BLAST output
run_blast <- function(query_fasta, db_path, out_tsv, blast_cmd = "blastn", nthread = 1) {
  args <- c("-query", query_fasta, "-db", db_path, "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue",
            "-num_threads", nthread, "-out", out_tsv)
  system2(blast_cmd, args)
  out_tsv
}

# Map BLAST hit IDs to GFF features (overlap by seqid and position)
# Simplified: assume BLAST subject id matches GFF seqid; we report genes overlapping hit positions
read_gff_genes <- function(gff_path) {
  if (!file.exists(gff_path)) return(NULL)
  d <- read_tsv(gff_path, comment = "#", col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"), show_col_types = FALSE, na = ".")
  d <- d %>% filter(.data$type %in% c("gene", "mRNA", "CDS"))
  d$start <- as.numeric(d$start)
  d$end <- as.numeric(d$end)
  # Extract gene_id/Name from attr if present
  if (nrow(d) > 0 && any(grepl("gene_id|Name|ID=", d$attr))) {
    d$gene_id <- gsub(".*gene_id=([^;]+).*", "\\1", d$attr)
    d$gene_id <- ifelse(d$gene_id == d$attr, gsub(".*Name=([^;]+).*", "\\1", d$attr), d$gene_id)
    d$gene_id <- ifelse(d$gene_id == d$attr, gsub(".*ID=([^;]+).*", "\\1", d$attr), d$gene_id)
  }
  if (!"gene_id" %in% names(d)) d$gene_id <- paste0(d$seqid, ":", d$start, "-", d$end)
  d
}

# Overlap BLAST hits (subject coords) with GFF genes
hits_to_genes <- function(blast_tsv, gff_path) {
  if (!file.exists(blast_tsv) || file.info(blast_tsv)$size == 0) return(NULL)
  b <- read_tsv(blast_tsv, col_names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue"), show_col_types = FALSE)
  if (nrow(b) == 0) return(NULL)
  gff <- read_gff_genes(gff_path)
  if (is.null(gff)) return(b %>% mutate(gene_id = NA_character_))
  b$hit_start <- pmin(b$sstart, b$send)
  b$hit_end <- pmax(b$sstart, b$send)
  genes <- gff %>% group_by(.data$seqid, .data$gene_id) %>% summarise(start = min(.data$start), end = max(.data$end), .groups = "drop")
  out <- b %>%
    left_join(genes, by = c("sseqid" = "seqid")) %>%
    filter(!is.na(.data$start), .data$hit_end >= .data$start, .data$hit_start <= .data$end) %>%
    distinct(.data$qseqid, .data$sseqid, .data$gene_id, .data$pident, .data$evalue)
  out
}

# Main
option_list <- list(
  make_option(c("--regions"), type = "character", default = NULL, help = "Regions CSV/TSV (from identify_outliers outlier_regions*.csv). Extra columns (outlier_stat, n_windows, region_mean_*, etc.) are preserved in output when present."),
  make_option(c("--reference"), type = "character", default = NULL, help = "Reference genome FASTA (indexed with samtools faidx)"),
  make_option(c("--blast-config"), type = "character", default = NULL, help = "BLAST config YAML/JSON (list of db_path, name, gff_path)"),
  make_option(c("--output-dir"), type = "character", default = ".", help = "Output directory"),
  make_option(c("--samtools"), type = "character", default = "samtools", help = "Path to samtools"),
  make_option(c("--blast-cmd"), type = "character", default = "blastn", help = "BLAST command (blastn/blastp/etc)"),
  make_option(c("--threads"), type = "integer", default = 1L, help = "Threads for BLAST"),
  make_option(c("--verbose"), action = "store_true", default = FALSE)
)
parser <- OptionParser(usage = "usage: %prog --regions FILE --reference FASTA --blast-config FILE [options]", option_list = option_list)
opts <- parse_args(parser)

if (is.null(opts$regions) || is.null(opts$reference) || is.null(opts$`blast-config`)) {
  stop("--regions, --reference, and --blast-config are required")
}

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)

regions <- read_regions(opts$regions)
if (nrow(regions) == 0) stop("No regions found in ", opts$regions)

config <- read_blast_config(opts$`blast-config`)
if (length(config) > 0 && !is.list(config[[1]])) config <- list(config)

# One FASTA per region (or combined); for simplicity one combined query FASTA with headers chr:start-end
query_fa <- file.path(opts$`output-dir`, "regions_query.fa")
regions$qseqid <- paste0(regions$chr, ":", regions$region_start, "-", regions$region_end)
write_region_fasta(regions, opts$reference, query_fa, opts$samtools)

all_genes <- list()
for (db in config) {
  db_path <- db$db_path %||% db[[1]]
  name <- db$name %||% db[[2]] %||% basename(db_path)
  gff_path <- db$gff_path %||% db$gff_path %||% db[[3]]
  out_blast <- file.path(opts$`output-dir`, paste0("blast_", gsub("[^A-Za-z0-9_]", "_", name), ".tsv"))
  run_blast(query_fa, db_path, out_blast, opts$`blast-cmd`, opts$threads)
  genes <- hits_to_genes(out_blast, gff_path)
  if (!is.null(genes)) {
    genes$db <- name
    all_genes[[length(all_genes) + 1L]] <- genes
  }
}

if (length(all_genes) > 0) {
  genes_out <- bind_rows(all_genes)
  # Option A: add region columns to each row (region_chr, region_start, region_end, and
  # any pass-through columns from the regions file when present)
  region_pass <- c("outlier_stat", "n_windows", "region_mean_coverage", "region_mean_mapping_quality", "region_n_snps")
  region_lookup <- regions %>%
    select("qseqid", "chr", "region_start", "region_end", any_of(region_pass)) %>%
    rename(region_chr = .data$chr)
  genes_out <- genes_out %>%
    left_join(region_lookup, by = "qseqid", multiple = "first")
  genes_file <- file.path(opts$`output-dir`, "outlier_regions_genes.csv")
  write_csv(genes_out, genes_file)
  message("Wrote ", genes_file, " (", nrow(genes_out), " hit-gene rows)")
} else {
  message("No BLAST hits or genes mapped.")
}

message("Done.")
