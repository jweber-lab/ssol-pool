#!/usr/bin/env Rscript

###############################################################################
# outlier_annotation.R
#
# Annotate outlier regions: extract sequences, BLAST per region against
# multiple DBs (from config file), map hits to GFF annotations, optional
# annotation TSV. Output: overlapping genes with region metadata.
#
# Regions file (e.g. from identify_outliers outlier_regions*.csv) may contain
# extra columns; when present they are preserved in the output.
#
# BLAST config: YAML or JSON with list of entries:
#   - db_path: path to BLAST DB (or FASTA for makeblastdb)
#   - name: label for this DB
#   - gff_path: path to GFF/GTF for mapping hits to genes (optional; warn if missing)
#   - annotation_tsv: optional TSV with gene_id column to join extra columns (e.g. product, GO)
#
# Optional --min-query-coverage FRAC (e.g. 0.8): only map hits to genes when the (qseqid,sseqid)
# pair has merged BLAST HSPs covering at least that fraction of the query region; reduces noise
# from short/local hits on subjects that do not align the full window.
#
# Outputs (filenames include a suffix from the regions file basename and, if set, --blast-task and --evalue):
#   outlier_regions_genes_<suffix>.csv: one row per (hit, overlapping gene); columns qseqid, sseqid,
#     gene_id, gene_name, product, gene_biotype, go_terms, ec_number, dbxref, strand,
#     feature_types, gff_source, pident, length, qstart, qend, sstart, send, evalue, db,
#     plus region pass-through (region_chr, region_start, region_end, etc.).
#   outlier_regions_genes_summary_<suffix>.csv: one row per (region, db); columns qseqid, db,
#     n_genes, gene_ids (comma-separated), best_evalue, products (semicolon-separated),
#     plus region pass-through.
# GFF/annotation fields are NA when not present; missing files produce warnings only.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
})

# Null-coalesce: return y when x is NULL (avoids infix for lint-clean script)
coalesce_null <- function(x, y) if (is.null(x)) y else x

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
  dbs <- coalesce_null(cfg$databases, cfg)
  if (!is.list(dbs)) dbs <- list(dbs)
  # If top-level is a single named list (one DB), wrap so we iterate over list of DBs
  if (length(dbs) > 0 && !is.null(names(dbs)) && !"db_path" %in% names(dbs))
    dbs <- list(dbs)
  if (length(dbs) == 0) stop("BLAST config has no 'databases' or list of { db_path, name, gff_path }")
  dbs
}

# Read regions file (CSV/TSV from identify_outliers: chr, region_start, region_end, ...).
# Retains extra columns for pass-through to output (outlier_stat, n_windows, region_*, etc.).
read_regions <- function(path) {
  first_line <- readLines(path, n = 1)
  if (length(first_line) == 0) stop("Regions file is empty: ", path)
  delim <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
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

# Write region FASTA using samtools faidx (reference must be indexed).
# Single samtools call with all regions for efficiency.
write_region_fasta <- function(regions, reference, out_fasta, samtools = "samtools") {
  dir.create(dirname(out_fasta), showWarnings = FALSE, recursive = TRUE)
  if (nrow(regions) == 0) {
    unlink(out_fasta)
    return(out_fasta)
  }
  fai <- paste0(reference, ".fai")
  if (!file.exists(fai)) stop("Reference index not found: ", fai)
  reg_vec <- paste0(regions$chr, ":", regions$region_start, "-", regions$region_end)
  seq <- system2(samtools, c("faidx", reference, reg_vec), stdout = TRUE, stderr = FALSE)
  writeLines(seq, out_fasta)
  out_fasta
}

# BLAST outfmt 6 column list (must be one argument; spaces cause split on some systems)
BLAST_OUTFMT_6 <- "6 qseqid sseqid pident length qstart qend sstart send evalue"

# Convert path to absolute and collapse . / .. (so BLAST gets a clean path; normalizePath can leave .. in place when path does not exist).
path_absolute <- function(path) {
  path <- trimws(as.character(path))
  if (is.na(path) || !nzchar(path)) return(path)
  if (!grepl("^[/~]", path)) path <- file.path(getwd(), path)
  path <- tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) path)
  # Collapse . and .. so BLAST sees a single canonical path (e.g. .../outliers/../../../genomes -> .../genomes)
  parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
  out <- character(0)
  for (p in parts) {
    if (p == "" || p == ".") next
    if (p == "..") { out <- head(out, -1L); next }
    out <- c(out, p)
  }
  if (length(out) == 0) return(path)
  paste(c("", out), collapse = "/")
}

# Run BLAST (blastn/blastp/etc) per region query; return path to BLAST output.
# Uses shell + quoting so -outfmt "6 qseqid ..." is passed as a single argument (avoids "Too many positional arguments: qseqid" on all BLAST variants).
# All paths are made absolute so the shell subprocess sees valid -query, -db, -out regardless of cwd.
run_blast <- function(query_fasta, db_path, out_tsv, blast_cmd = "blastn", nthread = 1, verbose = FALSE, blast_extra = character(0)) {
  db_path <- trimws(as.character(db_path))
  if (is.na(db_path) || !nzchar(db_path)) {
    message("BLAST skipped: db_path is missing or empty")
    return(out_tsv)
  }
  query_abs <- path_absolute(query_fasta)
  db_abs <- path_absolute(db_path)
  out_abs <- path_absolute(out_tsv)
  if (verbose) {
    message("  BLAST db_path: ", db_path)
    message("  BLAST db_path (absolute): ", db_abs)
    db_exists <- file.exists(db_abs)
    db_nin <- file.exists(paste0(db_abs, ".nin"))
    db_nsq <- file.exists(paste0(db_abs, ".nsq"))
    db_pin <- file.exists(paste0(db_abs, ".pin"))
    db_psq <- file.exists(paste0(db_abs, ".psq"))
    message("  DB path exists (file/dir): ", db_exists)
    message("  DB nucleotide index (.nin/.nsq): ", db_nin, " / ", db_nsq)
    message("  DB protein index (.pin/.psq): ", db_pin, " / ", db_psq)
    if (!db_exists && !db_nin && !db_pin)
      message("  (BLAST expects -db to be the DB prefix; ", db_abs, ".nin or .pin must exist)")
  }
  nthread <- as.character(as.integer(nthread))
  # NCBI BLAST manual (NBK569862): -outfmt must be one argument, e.g. -outfmt "6 qseqid sseqid ...".
  # Write command to a temp script and run with sh to avoid system2() mangling quotes/length when passing one long string (fixes "Too many positional arguments: qseqid" on some platforms).
  extra_str <- if (length(blast_extra) > 0) paste(paste(blast_extra, collapse = " ")) else ""
  cmd_line <- paste(
    blast_cmd,
    "-query", shQuote(query_abs),
    "-db", shQuote(db_abs),
    if (nzchar(extra_str)) extra_str,
    "-outfmt", shQuote(BLAST_OUTFMT_6),
    "-num_threads", nthread,
    "-out", shQuote(out_abs)
  )
  cmd_line <- paste(cmd_line[!is.na(cmd_line) & nzchar(cmd_line)], collapse = " ")
  if (verbose) message("  BLAST command: ", cmd_line)
  script <- tempfile(pattern = "blast_", fileext = ".sh")
  on.exit(unlink(script, force = TRUE))
  writeLines(cmd_line, script)
  stderr <- character(0)
  result <- system2("sh", script, stdout = FALSE, stderr = TRUE)
  exit <- attr(result, "status")
  if (is.null(exit)) exit <- 0L
  if (length(result) > 0) stderr <- result
  if (verbose || exit != 0L) {
    if (exit != 0L) message("  BLAST exit code: ", exit)
    if (length(stderr) > 0) message("  BLAST stderr: ", paste(stderr, collapse = "\n    "))
  }
  out_tsv
}

# Extract first capture group from GFF attr; returns NA_character_ where no match. Vectorized.
# Handles NA in attr: only indexes where we have a valid match (!is.na(m) & m > 0).
extract_attr_first <- function(attr, pattern) {
  out <- rep(NA_character_, length(attr))
  m <- regexpr(pattern, attr, perl = TRUE)
  hit <- !is.na(m) & (m > 0)
  if (!any(hit)) return(out)
  out[hit] <- gsub(pattern, "\\1", attr[hit], perl = TRUE)
  out
}

# Extract gene identifier from GFF attr (tries gene_id=, Name=, ID= in order). Vectorized.
# NA in attr remains NA (gsub/ifelse preserve NA).
extract_gff_gene_id <- function(attr) {
  out <- gsub(".*gene_id=([^;]+).*", "\\1", attr)
  out <- ifelse(is.na(out) | out == attr, gsub(".*Name=([^;]+).*", "\\1", attr), out)
  out <- ifelse(is.na(out) | out == attr, gsub(".*ID=([^;]+).*", "\\1", attr), out)
  out
}

# Extract optional GFF attr fields; missing or invalid attr yields NA. All vectorized.
# Returns list of vectors: gene_name, product, gene_biotype, go_terms, ec_number, dbxref.
extract_gff_attr_optional <- function(attr) {
  gene_name <- extract_attr_first(attr, "^(?:.*;)?(?:Name|gene_name)=([^;]+)")
  product <- extract_attr_first(attr, "^(?:.*;)?(?:product|description|Note)=([^;]+)")
  gene_biotype <- extract_attr_first(attr, "^(?:.*;)?(?:gene_biotype|biotype)=([^;]+)")
  ec_number <- extract_attr_first(attr, "^(?:.*;)?EC=([^;]+)")
  dbxref <- extract_attr_first(attr, "^(?:.*;)?Dbxref=([^;]+)")
  # GO: extract all GO: terms from Ontology_term= or Dbxref=, collapse with ";"
  go_raw <- extract_attr_first(attr, "^(?:.*;)?Ontology_term=([^;]+)")
  go_from_dbxref <- extract_attr_first(attr, "^(?:.*;)?Dbxref=([^;]+)")
  go_terms <- vapply(seq_along(attr), function(i) {
    a <- c(go_raw[i], go_from_dbxref[i])
    a <- a[!is.na(a) & nchar(a) > 0]
    if (length(a) == 0) return(NA_character_)
    gos <- unique(unlist(regmatches(a, gregexpr("GO:[0-9]+", a))))
    if (length(gos) == 0) return(NA_character_)
    paste(gos, collapse = ";")
  }, character(1))
  list(gene_name = gene_name, product = product, gene_biotype = gene_biotype, go_terms = go_terms, ec_number = ec_number, dbxref = dbxref)
}

# Read GFF/GTF and extract genes with optional fields (product, gene_name, strand, etc.).
# Returns NULL on missing file or read error (with warning). Graceful when attr lacks optional keys.
read_gff_genes <- function(gff_path, verbose = FALSE) {
  if (is.null(gff_path) || length(gff_path) == 0 || is.na(gff_path) || !nzchar(trimws(gff_path))) {
    if (verbose) message("GFF path missing or empty, skipping")
    return(NULL)
  }
  if (!file.exists(gff_path)) {
    if (verbose) message("GFF not found, skipping: ", gff_path)
    return(NULL)
  }
  d <- tryCatch({
    read_tsv(gff_path, comment = "#", col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"), show_col_types = FALSE, na = ".")
  }, error = function(e) {
    message("Warning: could not read GFF ", gff_path, ": ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(d) || nrow(d) == 0) return(d)
  d <- d %>% filter(.data$type %in% c("gene", "mRNA", "CDS"))
  if (nrow(d) == 0) return(d)
  d$start <- as.numeric(d$start)
  d$end <- as.numeric(d$end)
  d <- d %>% filter(!is.na(.data$start), !is.na(.data$end))
  if (nrow(d) == 0) return(d)
  has_attr <- !is.na(d$attr) & grepl("gene_id|Name|ID=", d$attr)
  if (any(has_attr)) {
    d$gene_id <- extract_gff_gene_id(d$attr)
    d$gene_id <- ifelse(is.na(d$gene_id) | !nzchar(trimws(d$gene_id)), paste0(d$seqid, ":", d$start, "-", d$end), d$gene_id)
  }
  if (!"gene_id" %in% names(d)) d$gene_id <- paste0(d$seqid, ":", d$start, "-", d$end)
  opt <- extract_gff_attr_optional(d$attr)
  d$gene_name <- opt$gene_name
  d$product <- opt$product
  d$gene_biotype <- opt$gene_biotype
  d$go_terms <- opt$go_terms
  d$ec_number <- opt$ec_number
  d$dbxref <- opt$dbxref
  # strand/type/source already in d
  d
}

# Normalize BLAST subject IDs so they can match GFF seqids (e.g. emb|CABIJS010000716.1| -> CABIJS010000716.1).
normalize_sseqid <- function(x) sub("^[^|]+\\|([^|]+)\\|?$", "\\1", x)

# Parse qseqid "chr:start-end" to get query length (region length).
query_length_from_qseqid <- function(qseqid) {
  start <- as.numeric(sub("^.*:([0-9]+)-[0-9]+$", "\\1", qseqid))
  end <- as.numeric(sub("^.*:[0-9]+-([0-9]+)$", "\\1", qseqid))
  pmax(1L, end - start + 1)
}

# Merge intervals [start, end] and return total covered length (handles overlapping and unordered).
merged_length <- function(start, end) {
  o <- order(start)
  start <- start[o]
  end <- end[o]
  n <- length(start)
  if (n == 0) return(0)
  total <- 0
  cur_end <- -Inf
  for (i in seq_len(n)) {
    if (start[i] > cur_end) {
      total <- total + (end[i] - start[i] + 1)
      cur_end <- end[i]
    } else if (end[i] > cur_end) {
      total <- total + (end[i] - cur_end)
      cur_end <- end[i]
    }
  }
  total
}

# Filter BLAST tibble to (qseqid, sseqid) pairs where query coverage >= min_coverage.
# b must have qseqid, sseqid, qstart, qend. regions must have qseqid and region_start, region_end (or we parse qseqid).
filter_by_query_coverage <- function(b, regions, min_coverage) {
  if (nrow(b) == 0 || min_coverage <= 0 || min_coverage > 1) return(b)
  # Query length: from regions if qseqid matches, else parse "chr:start-end"
  qlen <- if (!is.null(regions) && "qseqid" %in% names(regions) && "region_start" %in% names(regions)) {
    idx <- match(b$qseqid, regions$qseqid)
    (regions$region_end - regions$region_start + 1L)[idx]
  } else {
    query_length_from_qseqid(b$qseqid)
  }
  b$qlen <- qlen
  cov <- b %>%
    group_by(.data$qseqid, .data$sseqid) %>%
    summarise(covered = merged_length(.data$qstart, .data$qend), qlen = first(.data$qlen), .groups = "drop") %>%
    mutate(frac = .data$covered / .data$qlen)
  keep <- cov %>% filter(.data$frac >= min_coverage) %>% select("qseqid", "sseqid")
  b %>% select(-"qlen") %>% semi_join(keep, by = c("qseqid", "sseqid"))
}

# Overlap BLAST hits (subject coords) with GFF genes; keep BLAST cols and GFF annotation.
# blast_tsv: path to BLAST outfmt 6 TSV, or a tibble with same columns (e.g. after coverage filter).
hits_to_genes <- function(blast_tsv, gff_path, verbose = FALSE) {
  if (is.data.frame(blast_tsv)) {
    b <- blast_tsv
  } else {
    if (!file.exists(blast_tsv) || file.info(blast_tsv)$size == 0) return(NULL)
    b <- tryCatch({
      read_tsv(blast_tsv, col_names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue"), show_col_types = FALSE)
    }, error = function(e) {
      message("Warning: could not read BLAST output ", blast_tsv, ": ", conditionMessage(e))
      return(NULL)
    })
  }
  if (is.null(b) || nrow(b) == 0) return(b)
  b$sstart <- as.numeric(b$sstart)
  b$send <- as.numeric(b$send)
  b <- b %>% filter(!is.na(.data$sstart), !is.na(.data$send))
  if (nrow(b) == 0) return(b)
  # Normalize so NCBI-style sseqids (e.g. emb|ACCESSION|) match GFF seqid (ACCESSION)
  b$sseqid_for_join <- normalize_sseqid(b$sseqid)
  gff <- read_gff_genes(gff_path, verbose = verbose)
  if (is.null(gff)) {
    return(b %>% mutate(gene_id = NA_character_, gene_name = NA_character_, product = NA_character_, gene_biotype = NA_character_, go_terms = NA_character_, ec_number = NA_character_, dbxref = NA_character_, strand = NA_character_, feature_types = NA_character_, gff_source = NA_character_) %>% select(-any_of("sseqid_for_join")))
  }
  b$hit_start <- pmin(b$sstart, b$send)
  b$hit_end <- pmax(b$sstart, b$send)
  genes <- gff %>%
    group_by(.data$seqid, .data$gene_id) %>%
    summarise(
      start = min(.data$start, na.rm = TRUE),
      end = max(.data$end, na.rm = TRUE),
      gene_name = first(na.omit(.data$gene_name)),
      product = first(na.omit(.data$product)),
      gene_biotype = first(na.omit(.data$gene_biotype)),
      go_terms = first(na.omit(.data$go_terms)),
      ec_number = first(na.omit(.data$ec_number)),
      dbxref = first(na.omit(.data$dbxref)),
      strand = first(na.omit(.data$strand)),
      feature_types = paste(sort(unique(na.omit(.data$type))), collapse = ","),
      gff_source = first(na.omit(.data$source)),
      .groups = "drop"
    )
  # Overlap join: only (hit, gene) pairs that overlap in coordinates. Avoids Cartesian product
  # from joining on seqid only (which can produce billions of rows on large contigs).
  # Process per seqid, per hit (no full hits×genes matrix) to keep memory bounded.
  genes_by_seqid <- split(genes, genes$seqid)
  b_cols <- b %>% select(-"sseqid_for_join", -"hit_start", -"hit_end")
  out_list <- list()
  for (s in unique(b$sseqid_for_join)) {
    genes_s <- genes_by_seqid[[s]]
    if (is.null(genes_s) || nrow(genes_s) == 0) next
    rows_b <- which(b$sseqid_for_join == s)
    hits_s <- b[rows_b, ]
    for (j in seq_len(nrow(hits_s))) {
      idx <- which(genes_s$end >= hits_s$hit_start[j] & genes_s$start <= hits_s$hit_end[j])
      if (length(idx) == 0) next
      hit_part <- b_cols[rows_b[j], , drop = FALSE][rep(1L, length(idx)), , drop = FALSE]
      gene_part <- genes_s[idx, setdiff(names(genes_s), "seqid"), drop = FALSE]
      out_list[[length(out_list) + 1L]] <- dplyr::bind_cols(hit_part, gene_part)
    }
  }
  if (length(out_list) == 0) {
    return(tibble(
      qseqid = character(), sseqid = character(), gene_id = character(), gene_name = character(),
      product = character(), gene_biotype = character(), go_terms = character(), ec_number = character(),
      dbxref = character(), strand = character(), feature_types = character(), gff_source = character(),
      pident = numeric(), length = numeric(), qstart = numeric(), qend = numeric(),
      sstart = numeric(), send = numeric(), evalue = numeric()
    ))
  }
  out <- dplyr::bind_rows(out_list)
  out <- out %>%
    distinct(.data$qseqid, .data$sseqid, .data$gene_id, .keep_all = TRUE) %>%
    select("qseqid", "sseqid", "gene_id", "gene_name", "product", "gene_biotype", "go_terms", "ec_number", "dbxref", "strand", "feature_types", "gff_source", "pident", "length", "qstart", "qend", "sstart", "send", "evalue")
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
  make_option(c("--threads"), type = "integer", default = 1L, help = "Threads per BLAST run"),
  make_option(c("--parallel-dbs"), type = "integer", default = 1L, help = "Run this many BLAST DBs in parallel (Unix/macOS only; 1 = sequential)"),
  make_option(c("--blast-task"), type = "character", default = NA_character_, help = "BLAST -task (e.g. blastn for sensitive search; default megablast can miss divergent hits)"),
  make_option(c("--evalue"), type = "character", default = NA_character_, help = "BLAST -evalue threshold (default 10)"),
  make_option(c("--min-query-coverage"), type = "numeric", default = NA_real_, metavar = "FRAC",
    help = "Minimum fraction of query region covered by BLAST hits per (qseqid,sseqid). Only hits from subject sequences that align to at least this fraction of the query region are mapped to genes (e.g. 0.8 = 80%%). Omit or NA to keep all hits (backward compatible)."),
  make_option(c("--verbose"), action = "store_true", default = FALSE)
)
parser <- OptionParser(usage = "usage: %prog --regions FILE --reference FASTA --blast-config FILE [options]", option_list = option_list)
opts <- parse_args(parser)

if (is.null(opts$regions) || is.null(opts$reference) || is.null(opts$`blast-config`)) {
  stop("--regions, --reference, and --blast-config are required")
}
# In verbose mode, print warnings as they occur; at end we always dump any stored warnings to the log
if (opts$verbose) options(warn = 1)

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)

# Output filename suffix from regions file and optional run parameters (so outputs are distinguishable).
output_suffix <- tools::file_path_sans_ext(basename(opts$regions))
output_suffix <- sub("^outlier_regions_?", "", output_suffix)
if (!nzchar(output_suffix)) output_suffix <- "regions"
if (!is.na(opts$`blast-task`) && nzchar(trimws(opts$`blast-task`)))
  output_suffix <- paste0(output_suffix, "_task_", gsub("[^A-Za-z0-9_-]", "_", trimws(opts$`blast-task`)))
if (!is.na(opts$evalue) && nzchar(trimws(opts$evalue)))
  output_suffix <- paste0(output_suffix, "_evalue_", gsub("[^A-Za-z0-9.+-]", "_", trimws(opts$evalue)))

regions <- read_regions(opts$regions)
if (nrow(regions) == 0) stop("No regions found in ", opts$regions)

config <- read_blast_config(opts$`blast-config`)
if (length(config) > 0 && !is.list(config[[1]])) config <- list(config)

# One FASTA per region (or combined); for simplicity one combined query FASTA with headers chr:start-end
query_fa <- file.path(opts$`output-dir`, "regions_query.fa")
regions$qseqid <- paste0(regions$chr, ":", regions$region_start, "-", regions$region_end)
invisible(write_region_fasta(regions, opts$reference, query_fa, opts$samtools))
if (opts$verbose) message("Query FASTA: ", normalizePath(query_fa, mustWork = FALSE))

# Per-DB work: run BLAST, map hits to genes, optional annotation join. Returns genes tibble or NULL.
process_one_db <- function(db) {
  # Prefer explicit db_path (YAML order can make db[[1]] the name, not the path)
  db_path <- coalesce_null(db$db_path, db[["db_path"]])
  db_path <- if (is.null(db_path)) NULL else trimws(as.character(db_path))
  name <- coalesce_null(coalesce_null(db$name, db[[2]]), if (nzchar(db_path)) basename(db_path) else "unknown")
  if (is.null(db_path) || !nzchar(db_path)) {
    message("Skipping DB ", name, ": db_path missing or empty in config")
    return(NULL)
  }
  gff_path <- coalesce_null(db$gff_path, db[[3]])
  annotation_tsv <- coalesce_null(db$annotation_tsv, db$annotation_table)
  out_blast <- file.path(opts$`output-dir`, paste0("blast_", gsub("[^A-Za-z0-9_]", "_", name), ".tsv"))
  if (opts$verbose) message("Running BLAST for DB: ", name)
  blast_extra <- character(0)
  if (!is.na(opts$`blast-task`) && nzchar(trimws(opts$`blast-task`))) blast_extra <- c(blast_extra, "-task", trimws(opts$`blast-task`))
  if (!is.na(opts$evalue) && nzchar(trimws(opts$evalue))) blast_extra <- c(blast_extra, "-evalue", trimws(opts$evalue))
  run_blast(query_fa, db_path, out_blast, opts$`blast-cmd`, opts$threads, verbose = opts$verbose, blast_extra = blast_extra)
  min_cov <- opts$`min-query-coverage`
  use_coverage_filter <- is.numeric(min_cov) && is.finite(min_cov) && min_cov > 0 && min_cov <= 1
  if (use_coverage_filter) {
    if (!file.exists(out_blast) || file.info(out_blast)$size == 0) {
      genes <- NULL
    } else {
      b <- read_tsv(out_blast, col_names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue"), show_col_types = FALSE)
      b <- filter_by_query_coverage(b, regions, min_cov)
      genes <- hits_to_genes(b, gff_path, verbose = opts$verbose)
    }
  } else {
    genes <- hits_to_genes(out_blast, gff_path, verbose = opts$verbose)
  }
  if (is.null(genes)) return(NULL)
  genes$db <- name
  anno_ok <- !is.null(annotation_tsv) && length(annotation_tsv) > 0 && !is.na(annotation_tsv) && nzchar(trimws(annotation_tsv))
  if (anno_ok) {
    if (!file.exists(annotation_tsv)) {
      message("Warning: optional annotation file not found, skipping: ", annotation_tsv)
    } else {
      anno <- tryCatch({
        read_tsv(annotation_tsv, show_col_types = FALSE) %>% rename_all(tolower)
      }, error = function(e) {
        message("Warning: could not read annotation file ", annotation_tsv, ": ", conditionMessage(e))
        NULL
      })
      if (!is.null(anno) && nrow(anno) > 0 && "gene_id" %in% names(anno)) {
        idx <- match(genes$gene_id, anno$gene_id)
        common <- setdiff(intersect(names(genes), names(anno)), "gene_id")
        for (col in common) {
          fill <- anno[[col]][idx]
          genes[[col]] <- coalesce(genes[[col]], as.character(fill))
        }
        extra <- setdiff(names(anno), names(genes))
        if (length(extra) > 0) {
          genes <- genes %>%
            left_join(anno %>% distinct(.data$gene_id, .keep_all = TRUE) %>% select("gene_id", any_of(extra)), by = "gene_id", multiple = "first")
        }
      } else if (!is.null(anno) && (!"gene_id" %in% names(anno))) {
        if (opts$verbose) message("Warning: annotation file has no gene_id column: ", annotation_tsv)
      }
    }
  }
  genes
}

n_parallel <- max(1L, min(as.integer(opts$`parallel-dbs`), length(config)))
use_parallel <- n_parallel > 1L && .Platform$OS.type == "unix"
if (use_parallel && opts$verbose) message("Running ", n_parallel, " BLAST DBs in parallel")

if (use_parallel) {
  all_genes <- parallel::mclapply(config, process_one_db, mc.cores = n_parallel)
  # When a worker errors, its element may be an error object, not a data frame; keep only valid results.
  all_genes <- all_genes[vapply(all_genes, is.data.frame, logical(1))]
  n_failed <- length(config) - length(all_genes)
  if (n_failed > 0)
    message("Warning: ", n_failed, " DB(s) failed in parallel (re-run with --parallel-dbs 1 to see the error).")
} else {
  if (n_parallel > 1L && .Platform$OS.type != "unix")
    message("Note: --parallel-dbs > 1 is only supported on Unix/macOS; running sequentially")
  all_genes <- list()
  for (db in config) {
    genes <- process_one_db(db)
    if (!is.null(genes) && is.data.frame(genes)) all_genes[[length(all_genes) + 1L]] <- genes
  }
}

if (length(all_genes) > 0) {
  genes_out <- bind_rows(all_genes)
  # Pass through all region columns (chr -> region_chr; region_start, region_end; plus any extra)
  region_extra <- setdiff(names(regions), c("qseqid", "chr", "region_start", "region_end"))
  region_lookup <- regions %>%
    select("qseqid", "chr", "region_start", "region_end", any_of(region_extra)) %>%
    rename(region_chr = chr)
  genes_out <- genes_out %>%
    left_join(region_lookup, by = "qseqid", multiple = "first")
  genes_file <- file.path(opts$`output-dir`, paste0("outlier_regions_genes_", output_suffix, ".csv"))
  write_csv(genes_out, genes_file)
  message("Wrote ", genes_file, " (", nrow(genes_out), " hit-gene rows)")

  # Summary: one row per (region, db) with comma-separated gene_ids, n_genes, best_evalue, products
  summary_out <- genes_out %>%
    group_by(.data$qseqid, .data$db) %>%
    summarise(
      n_genes = dplyr::n_distinct(.data$gene_id),
      gene_ids = paste(sort(unique(na.omit(.data$gene_id))), collapse = ","),
      best_evalue = min(as.numeric(.data$evalue), na.rm = TRUE),
      products = paste(sort(unique(na.omit(.data$product))), collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(best_evalue = if_else(is.finite(.data$best_evalue), .data$best_evalue, NA_real_)) %>%
    left_join(region_lookup, by = "qseqid", multiple = "first")
  summary_file <- file.path(opts$`output-dir`, paste0("outlier_regions_genes_summary_", output_suffix, ".csv"))
  write_csv(summary_out, summary_file)
  message("Wrote ", summary_file, " (", nrow(summary_out), " region-DB rows)")

  if (nrow(genes_out) == 0)
    message("Tip: for divergent taxa, try --blast-task blastn (more sensitive than default megablast) and/or inspect the raw BLAST TSV files in ", opts$`output-dir`)
} else {
  message("No BLAST hits or genes mapped.")
  message("Tip: try --blast-task blastn for more sensitive search; check raw BLAST TSV files in ", opts$`output-dir`)
}

w <- warnings()
if (length(w) > 0) {
  message("Warnings (", length(w), "):")
  print(w)
}

message("Done.")
