#!/usr/bin/env Rscript

###############################################################################
# gene_collation.R
#
# Reorganize outlier_annotation.sh output by assembly GFF gene, instead of by
# variant or region. The GFF is the gene universe; annotation CSVs supply the
# outlier regions, csq sites, and optional BLAST homologs (joined by region).
#
# Inputs: --gff plus either --annotation-dir (discovers outlier_csq_*.csv and
# optional outlier_genes_*.csv) or explicit --sites / --genes / --summary /
# --blast-genes. --suffix picks one annotation run when several exist.
#
# Outputs (suffix from the annotation files, or "genes"):
#   outlier_gene_summary_<suffix>.csv  One row per GFF gene that overlaps an
#     outlier region (or every gene with --include-non-outlier). GFF fields,
#     n_regions, qseqids, outlier_stat, site counts / coding_effects, and BLAST
#     homologs collapsed across overlapping regions (other-assembly gene_ids).
#   outlier_gene_regions_<suffix>.csv  One row per (gene, region).
#   outlier_gene_sites_<suffix>.csv    One row per (gene, site): sites whose
#     position falls in the gene and/or whose csq `genes` field names it.
#     qseqids are semicolon-collapsed when a site sits in several regions.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
  library(tidyr)
})

coalesce_null <- function(x, y) if (is.null(x)) y else x

extract_attr_first <- function(attr, pattern) {
  out <- rep(NA_character_, length(attr))
  ok <- !is.na(attr)
  if (!any(ok)) return(out)
  m <- regmatches(attr[ok], regexec(pattern, attr[ok], perl = TRUE))
  out[ok] <- vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  out
}

extract_gff_gene_id <- function(attr) {
  out <- gsub(".*gene_id=([^;]+).*", "\\1", attr)
  out <- ifelse(is.na(out) | out == attr, gsub(".*Name=([^;]+).*", "\\1", attr), out)
  out <- ifelse(is.na(out) | out == attr, gsub(".*ID=([^;]+).*", "\\1", attr), out)
  out
}

extract_gff_attr_optional <- function(attr) {
  gene_name <- extract_attr_first(attr, "^(?:.*;)?(?:Name|gene_name)=([^;]+)")
  product <- extract_attr_first(attr, "^(?:.*;)?(?:product|description|Note)=([^;]+)")
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
  list(gene_name = gene_name, product = product, go_terms = go_terms)
}

read_assembly_genes <- function(gff_path, verbose = FALSE) {
  if (is.null(gff_path) || !file.exists(gff_path)) stop("GFF not found: ", gff_path)
  d <- tryCatch({
    read_tsv(gff_path, comment = "#",
             col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
             show_col_types = FALSE, na = ".")
  }, error = function(e) stop("Could not read GFF ", gff_path, ": ", conditionMessage(e)))
  if (nrow(d) == 0) stop("GFF is empty: ", gff_path)
  d <- d %>% filter(.data$type %in% c("gene", "mRNA", "CDS"))
  if (nrow(d) == 0) stop("No gene/mRNA/CDS features in ", gff_path)
  d$start <- as.numeric(d$start)
  d$end <- as.numeric(d$end)
  d <- d %>% filter(!is.na(.data$start), !is.na(.data$end))
  d$gene_id <- extract_gff_gene_id(d$attr)
  d$gene_id <- ifelse(is.na(d$gene_id) | !nzchar(trimws(d$gene_id)),
                      paste0(d$seqid, ":", d$start, "-", d$end), d$gene_id)
  opt <- extract_gff_attr_optional(d$attr)
  d$gene_name <- opt$gene_name
  d$product <- opt$product
  d$go_terms <- opt$go_terms
  id_attr <- extract_attr_first(d$attr, "^(?:.*;)?ID=([^;]+)")
  d$csq_key <- dplyr::coalesce(extract_attr_first(d$attr, "^(?:.*;)?Name=([^;]+)"),
                               sub("^gene:", "", id_attr))
  g <- if (any(d$type == "gene")) d %>% filter(.data$type == "gene") else d
  first_non_na <- function(x) {
    x <- x[!is.na(x) & !(is.character(x) & !nzchar(x))]
    if (length(x) == 0) return(NA)
    x[[1]]
  }
  out <- g %>%
    group_by(.data$seqid, .data$gene_id) %>%
    summarise(
      gene_start = min(.data$start), gene_end = max(.data$end),
      csq_key = first_non_na(.data$csq_key),
      gene_name = first_non_na(.data$gene_name),
      product = first_non_na(.data$product),
      go_terms = first_non_na(.data$go_terms),
      strand = first_non_na(.data$strand),
      .groups = "drop"
    ) %>%
    mutate(chr = .data$seqid, gene_length = .data$gene_end - .data$gene_start + 1L)
  if (verbose) message("GFF genes: ", nrow(out), " from ", gff_path)
  out
}

read_table_auto <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  first_line <- readLines(path, n = 1)
  if (length(first_line) == 0) return(NULL)
  delim <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
  if (delim == "\t") read_tsv(path, show_col_types = FALSE) else read_csv(path, show_col_types = FALSE)
}

find_col <- function(d, candidates) {
  nms <- names(d)
  low <- tolower(nms)
  for (cnd in candidates) {
    hit <- which(low == tolower(cnd))
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

rename_known <- function(d, mapping) {
  for (i in seq_along(mapping)) {
    src <- find_col(d, mapping[[i]])
    dest <- names(mapping)[i]
    if (!is.na(src) && src != dest && !dest %in% names(d)) names(d)[names(d) == src] <- dest
  }
  d
}

split_semicolon <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(character(0))
  v <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  v[nzchar(v)]
}

join_unique <- function(x) {
  v <- sort(unique(split_semicolon(x)))
  if (length(v) == 0) NA_character_ else paste(v, collapse = ";")
}

effect_counts_string <- function(x) {
  x <- x[!is.na(x) & x != "none" & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  tb <- sort(table(x), decreasing = TRUE)
  paste0(names(tb), ":", as.integer(tb), collapse = ";")
}

safe_min <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

safe_sum <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else sum(x)
}

annotation_suffix <- function(path) {
  bn <- tools::file_path_sans_ext(basename(path))
  sub("^outlier_csq_(sites|genes|summary)_", "", bn, perl = TRUE)
}

pick_one <- function(files, kind, suffix, verbose) {
  if (length(files) == 0) return(NULL)
  if (!is.null(suffix) && nzchar(suffix)) {
    files <- files[grepl(suffix, basename(files), fixed = TRUE)]
    if (length(files) == 0) return(NULL)
  }
  if (length(files) > 1L) {
    stop("Several ", kind, " files found. Pass --suffix to choose one:\n  ",
         paste(basename(files), collapse = "\n  "))
  }
  if (verbose) message("Using ", kind, ": ", files)
  files
}

discover_annotation <- function(dir, suffix, verbose) {
  if (is.null(dir) || !dir.exists(dir)) return(list())
  list(
    sites = pick_one(list.files(dir, pattern = "^outlier_csq_sites_.*\\.(csv|tsv)$", full.names = TRUE),
                     "sites", suffix, verbose),
    genes = pick_one(list.files(dir, pattern = "^outlier_csq_genes_.*\\.(csv|tsv)$", full.names = TRUE),
                     "csq genes", suffix, verbose),
    summary = pick_one(list.files(dir, pattern = "^outlier_csq_summary_.*\\.(csv|tsv)$", full.names = TRUE),
                       "summary", suffix, verbose),
    blast = pick_one(
      grep("outlier_genes_summary_", list.files(dir, pattern = "^outlier_genes_.*\\.(csv|tsv)$", full.names = TRUE),
           value = TRUE, invert = TRUE),
      "BLAST genes", suffix, verbose)
  )
}

# Interval overlap: left rows vs right intervals, by chromosome.
overlap_hits <- function(left_chr, left_start, left_end, right) {
  if (length(left_chr) == 0 || nrow(right) == 0)
    return(tibble(left_row = integer(), region_row = integer()))
  hits <- list()
  by_chr <- split(seq_along(left_chr), left_chr)
  for (i in seq_len(nrow(right))) {
    idx <- by_chr[[as.character(right$chr[i])]]
    if (is.null(idx)) next
    keep <- left_end[idx] >= right$region_start[i] & left_start[idx] <= right$region_end[i]
    if (any(keep)) hits[[length(hits) + 1L]] <- tibble(left_row = idx[keep], region_row = i)
  }
  bind_rows(hits)
}

gene_keys_row <- function(gene_id, csq_key, gene_name) {
  v <- unique(c(gene_id, sub("^gene:", "", gene_id), csq_key, gene_name))
  v[!is.na(v) & nzchar(v)]
}

# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--gff"), type = "character", default = NULL,
              help = "Assembly GFF3 (the one given to variant_csq.sh / outlier_annotation.sh)"),
  make_option(c("--annotation-dir"), type = "character", default = NULL,
              help = "Directory of outlier_annotation.sh output (outlier_csq_*.csv, optional outlier_genes_*.csv)"),
  make_option(c("--suffix"), type = "character", default = NULL,
              help = "Pick annotation files whose names contain this suffix (required when several runs share a directory)"),
  make_option(c("--sites"), type = "character", default = NULL, help = "outlier_csq_sites_*.csv"),
  make_option(c("--genes"), type = "character", default = NULL, help = "outlier_csq_genes_*.csv"),
  make_option(c("--summary"), type = "character", default = NULL, help = "outlier_csq_summary_*.csv"),
  make_option(c("--blast-genes"), type = "character", default = NULL, help = "outlier_genes_*.csv (BLAST hits)"),
  make_option(c("--include-non-outlier"), action = "store_true", default = FALSE,
              help = "Keep GFF genes that overlap no outlier region"),
  make_option(c("--output-dir"), type = "character", default = ".", help = "Output directory [default: %default]"),
  make_option(c("--verbose"), action = "store_true", default = FALSE)
)
opts <- parse_args(OptionParser(
  usage = "usage: %prog --gff FILE --annotation-dir DIR [options]",
  option_list = option_list
))
if (is.null(opts$gff)) stop("--gff is required")
if (is.null(opts$`annotation-dir`) && is.null(opts$sites) && is.null(opts$genes) && is.null(opts$summary))
  stop("Give --annotation-dir or at least one of --sites, --genes, --summary")

found <- discover_annotation(opts$`annotation-dir`, opts$suffix, opts$verbose)
sites_path <- coalesce_null(opts$sites, found$sites)
genes_path <- coalesce_null(opts$genes, found$genes)
summary_path <- coalesce_null(opts$summary, found$summary)
blast_path <- coalesce_null(opts$`blast-genes`, found$blast)
if (is.null(sites_path) && is.null(genes_path) && is.null(summary_path))
  stop("No outlier_csq_sites/genes/summary files found")

suffix <- opts$suffix
if (is.null(suffix) || !nzchar(suffix)) {
  src <- coalesce_null(summary_path, coalesce_null(sites_path, genes_path))
  suffix <- annotation_suffix(src)
}
if (!nzchar(suffix)) suffix <- "genes"

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)
asm <- read_assembly_genes(opts$gff, verbose = opts$verbose)

sites <- read_table_auto(sites_path)
genes_ann <- read_table_auto(genes_path)
summary_ann <- read_table_auto(summary_path)
blast <- read_table_auto(blast_path)

sites <- if (is.null(sites)) sites else rename_known(sites, list(
  qseqid = "qseqid", region_chr = c("region_chr", "chr"),
  region_start = "region_start", region_end = "region_end",
  chr = "chr", pos = "pos"
))
genes_ann <- if (is.null(genes_ann)) genes_ann else rename_known(genes_ann, list(
  qseqid = "qseqid", gene = c("gene", "gene_id", "csq_key"),
  region_chr = c("region_chr", "chr"), region_start = "region_start", region_end = "region_end"
))
summary_ann <- if (is.null(summary_ann)) summary_ann else rename_known(summary_ann, list(
  qseqid = "qseqid", region_chr = c("region_chr", "chr"),
  region_start = "region_start", region_end = "region_end"
))
blast <- if (is.null(blast)) blast else rename_known(blast, list(qseqid = "qseqid"))

# Regions: prefer summary (one row per region, with pass-through stats); else unique from sites/genes.
regions <- NULL
if (!is.null(summary_ann) && nrow(summary_ann) > 0) {
  regions <- summary_ann
} else {
  src <- bind_rows(
    if (!is.null(sites) && nrow(sites) > 0) sites %>% select(any_of(c("qseqid", "region_chr", "region_start", "region_end", "outlier_stat", "n_windows"))) else NULL,
    if (!is.null(genes_ann) && nrow(genes_ann) > 0) genes_ann %>% select(any_of(c("qseqid", "region_chr", "region_start", "region_end", "outlier_stat", "n_windows"))) else NULL
  )
  if (!is.null(src) && nrow(src) > 0) {
    regions <- src %>%
      mutate(qseqid = if ("qseqid" %in% names(.)) .data$qseqid else paste0(.data$region_chr, ":", .data$region_start, "-", .data$region_end)) %>%
      distinct(.data$qseqid, .keep_all = TRUE)
  }
}
if (is.null(regions) || nrow(regions) == 0) stop("No outlier regions in the annotation files")
if (!"region_chr" %in% names(regions) && "chr" %in% names(regions)) regions$region_chr <- regions$chr
if (!all(c("region_chr", "region_start", "region_end") %in% names(regions)))
  stop("Annotation regions need region_chr, region_start, region_end (or qseqid plus those columns)")
regions$region_start <- as.numeric(regions$region_start)
regions$region_end <- as.numeric(regions$region_end)
if (!"qseqid" %in% names(regions)) {
  regions$qseqid <- paste0(regions$region_chr, ":", regions$region_start, "-", regions$region_end)
} else {
  miss <- is.na(regions$qseqid) | !nzchar(as.character(regions$qseqid))
  if (any(miss))
    regions$qseqid[miss] <- paste0(regions$region_chr[miss], ":", regions$region_start[miss], "-", regions$region_end[miss])
}
regions$chr <- regions$region_chr
if (opts$verbose) message("Regions: ", nrow(regions))

# GFF genes overlapping each region.
olap <- overlap_hits(asm$chr, asm$gene_start, asm$gene_end, regions)
gene_regions <- if (nrow(olap) == 0) {
  tibble(gene_id = character(), qseqid = character())
} else {
  bind_cols(
    asm[olap$left_row, c("gene_id", "csq_key", "gene_name", "product", "go_terms", "strand",
                         "chr", "gene_start", "gene_end", "gene_length"), drop = FALSE],
    regions[olap$region_row, , drop = FALSE] %>% select(-any_of(c("chr")))
  ) %>% mutate(overlaps_region = TRUE)
}

# Csq-named genes that missed the GFF overlap (name mismatch): add from outlier_csq_genes.
if (!is.null(genes_ann) && nrow(genes_ann) > 0 && "gene" %in% names(genes_ann)) {
  key_map <- lapply(seq_len(nrow(asm)), function(i) gene_keys_row(asm$gene_id[i], asm$csq_key[i], asm$gene_name[i]))
  names(key_map) <- asm$gene_id
  key_to_id <- unlist(lapply(names(key_map), function(id) setNames(rep(id, length(key_map[[id]])), key_map[[id]])))
  genes_ann$gene_id_mapped <- unname(key_to_id[as.character(genes_ann$gene)])
  extra <- genes_ann %>% filter(is.na(.data$gene_id_mapped))
  matched <- genes_ann %>% filter(!is.na(.data$gene_id_mapped))
  if (nrow(matched) > 0) {
    have <- paste(gene_regions$gene_id, gene_regions$qseqid, sep = "\t")
    add <- matched %>%
      transmute(gene_id = .data$gene_id_mapped, qseqid = .data$qseqid) %>%
      filter(!paste(.data$gene_id, .data$qseqid, sep = "\t") %in% have)
    if (nrow(add) > 0) {
      add_asm <- asm %>% semi_join(add, by = "gene_id")
      add_reg <- regions %>% semi_join(add, by = "qseqid")
      add_rows <- add %>%
        left_join(add_asm, by = "gene_id") %>%
        left_join(add_reg %>% select(-any_of("chr")), by = "qseqid") %>%
        mutate(overlaps_region = TRUE)
      gene_regions <- bind_rows(gene_regions, add_rows)
    }
  }
  if (nrow(extra) > 0 && opts$verbose)
    message("Csq gene names not in GFF (dropped from gene tables): ",
            paste(unique(extra$gene), collapse = ", "))
}

# Sites -> genes by position and/or csq name.
gene_sites <- tibble()
if (!is.null(sites) && nrow(sites) > 0 && "pos" %in% names(sites)) {
  sites$pos <- as.numeric(sites$pos)
  if (!"chr" %in% names(sites) && "region_chr" %in% names(sites)) sites$chr <- sites$region_chr
  site_dummy_end <- sites$pos
  pos_hits <- overlap_hits(sites$chr, sites$pos, site_dummy_end,
                           asm %>% transmute(chr = .data$chr, region_start = .data$gene_start, region_end = .data$gene_end))
  pos_pairs <- if (nrow(pos_hits) == 0) tibble(site_row = integer(), gene_id = character(), assigned_by = character()) else
    tibble(site_row = pos_hits$left_row, gene_id = asm$gene_id[pos_hits$region_row], assigned_by = "position")

  name_pairs <- tibble(site_row = integer(), gene_id = character(), assigned_by = character())
  if ("genes" %in% names(sites)) {
    key_to_id <- unlist(lapply(seq_len(nrow(asm)), function(i) {
      k <- gene_keys_row(asm$gene_id[i], asm$csq_key[i], asm$gene_name[i])
      setNames(rep(asm$gene_id[i], length(k)), k)
    }))
    for (i in seq_len(nrow(sites))) {
      names_i <- split_semicolon(sites$genes[i])
      ids <- unique(unname(key_to_id[names_i]))
      ids <- ids[!is.na(ids)]
      if (length(ids)) name_pairs <- bind_rows(name_pairs, tibble(site_row = i, gene_id = ids, assigned_by = "csq"))
    }
  }
  pairs <- bind_rows(pos_pairs, name_pairs) %>%
    group_by(.data$site_row, .data$gene_id) %>%
    summarise(assigned_by = paste(sort(unique(.data$assigned_by)), collapse = "+"), .groups = "drop")
  if (nrow(pairs) > 0) {
    site_keep <- setdiff(names(sites), c("gene_id", "gene_name", "csq_key", "product", "go_terms",
                                         "strand", "gene_start", "gene_end", "gene_length", "assigned_by"))
    gene_sites <- bind_cols(
      asm[match(pairs$gene_id, asm$gene_id), c("gene_id", "csq_key", "gene_name", "product", "go_terms",
                                               "strand", "chr", "gene_start", "gene_end", "gene_length"), drop = FALSE] %>%
        rename(gene_chr = chr),
      sites[pairs$site_row, site_keep, drop = FALSE],
      tibble(assigned_by = pairs$assigned_by)
    )
    # One row per (gene, site); collapse regions.
    site_key <- intersect(c("chr", "pos", "ref", "alt"), names(gene_sites))
    grp <- c("gene_id", site_key)
    list_cols <- intersect(c("qseqid", "outlier_stat"), names(gene_sites))
    other <- setdiff(names(gene_sites), c(grp, list_cols, "assigned_by",
                                          "region_chr", "region_start", "region_end", "n_windows"))
    gene_sites <- gene_sites %>%
      group_by(across(all_of(grp))) %>%
      summarise(
        qseqids = if ("qseqid" %in% names(gene_sites)) join_unique(.data$qseqid) else NA_character_,
        outlier_stat = if ("outlier_stat" %in% names(gene_sites)) join_unique(.data$outlier_stat) else NA_character_,
        assigned_by = paste(sort(unique(.data$assigned_by)), collapse = "+"),
        across(any_of(other), ~ {
          v <- .[!is.na(.)]
          if (is.character(v)) v <- v[nzchar(v)]
          if (length(v) == 0) NA else v[[1]]
        }),
        .groups = "drop"
      )
  }
}

# Per-(gene, region) site counts from gene_sites qseqids.
if (nrow(gene_sites) > 0 && "qseqids" %in% names(gene_sites) && nrow(gene_regions) > 0) {
  gs_long <- gene_sites %>%
    mutate(qseqid = strsplit(coalesce(.data$qseqids, ""), ";", fixed = TRUE)) %>%
    unnest(cols = "qseqid") %>%
    mutate(qseqid = trimws(.data$qseqid)) %>%
    filter(nzchar(.data$qseqid))
  per <- gs_long %>%
    group_by(.data$gene_id, .data$qseqid) %>%
    summarise(
      n_sites = dplyr::n(),
      n_sites_snp = if ("variant_type" %in% names(gs_long)) sum(.data$variant_type == "snp", na.rm = TRUE) else NA_integer_,
      n_sites_indel = if ("variant_type" %in% names(gs_long)) sum(.data$variant_type == "indel", na.rm = TRUE) else NA_integer_,
      n_protein_altering = if ("genes_protein_altering" %in% names(gs_long))
        sum(!is.na(.data$genes_protein_altering) & nzchar(.data$genes_protein_altering)) else 0L,
      coding_effects = if ("coding_effect" %in% names(gs_long)) effect_counts_string(.data$coding_effect) else NA_character_,
      .groups = "drop"
    )
  gene_regions <- gene_regions %>%
    select(-any_of(c("n_sites", "n_sites_snp", "n_sites_indel", "n_protein_altering", "coding_effects"))) %>%
    left_join(per, by = c("gene_id", "qseqid"))
}
if (nrow(gene_regions) > 0) {
  for (col in c("n_sites", "n_sites_snp", "n_sites_indel", "n_protein_altering")) {
    if (col %in% names(gene_regions)) gene_regions[[col]] <- coalesce(as.integer(gene_regions[[col]]), 0L)
  }
}

# BLAST homologs attach by qseqid (other assembly), not by gene_id.
blast_by_region <- NULL
if (!is.null(blast) && nrow(blast) > 0 && "qseqid" %in% names(blast)) {
  blast_by_region <- blast %>%
    group_by(.data$qseqid) %>%
    summarise(
      n_blast_hits = dplyr::n(),
      blast_dbs = if ("db" %in% names(blast)) join_unique(.data$db) else NA_character_,
      blast_gene_ids = if ("gene_id" %in% names(blast)) join_unique(.data$gene_id) else NA_character_,
      blast_products = if ("product" %in% names(blast)) join_unique(.data$product) else NA_character_,
      best_evalue = if ("evalue" %in% names(blast)) safe_min(as.numeric(.data$evalue)) else NA_real_,
      .groups = "drop"
    )
  gene_regions <- gene_regions %>% left_join(blast_by_region, by = "qseqid")
}

# One row per gene.
if (opts$`include-non-outlier`) {
  gene_ids <- asm$gene_id
} else {
  gene_ids <- unique(c(gene_regions$gene_id, if (nrow(gene_sites) > 0) gene_sites$gene_id else character(0)))
}
gene_summary <- asm %>% filter(.data$gene_id %in% gene_ids)
if (nrow(gene_regions) > 0) {
  reg_agg <- gene_regions %>%
    group_by(.data$gene_id) %>%
    summarise(
      n_regions = dplyr::n_distinct(.data$qseqid),
      qseqids = join_unique(.data$qseqid),
      outlier_stat = if ("outlier_stat" %in% names(gene_regions)) join_unique(.data$outlier_stat) else NA_character_,
      n_windows = if ("n_windows" %in% names(gene_regions)) safe_sum(as.numeric(.data$n_windows)) else NA_real_,
      n_blast_hits = if ("n_blast_hits" %in% names(gene_regions)) safe_sum(.data$n_blast_hits) else NA_real_,
      blast_dbs = if ("blast_dbs" %in% names(gene_regions)) join_unique(.data$blast_dbs) else NA_character_,
      blast_gene_ids = if ("blast_gene_ids" %in% names(gene_regions)) join_unique(.data$blast_gene_ids) else NA_character_,
      blast_products = if ("blast_products" %in% names(gene_regions)) join_unique(.data$blast_products) else NA_character_,
      best_evalue = if ("best_evalue" %in% names(gene_regions)) safe_min(.data$best_evalue) else NA_real_,
      .groups = "drop"
    )
  gene_summary <- gene_summary %>% left_join(reg_agg, by = "gene_id")
}
if (nrow(gene_sites) > 0) {
  site_agg <- gene_sites %>%
    group_by(.data$gene_id) %>%
    summarise(
      n_sites = dplyr::n(),
      n_sites_snp = if ("variant_type" %in% names(gene_sites)) sum(.data$variant_type == "snp", na.rm = TRUE) else NA_integer_,
      n_sites_indel = if ("variant_type" %in% names(gene_sites)) sum(.data$variant_type == "indel", na.rm = TRUE) else NA_integer_,
      n_protein_altering = if ("genes_protein_altering" %in% names(gene_sites))
        sum(!is.na(.data$genes_protein_altering) & nzchar(.data$genes_protein_altering)) else 0L,
      coding_effects = if ("coding_effect" %in% names(gene_sites)) effect_counts_string(.data$coding_effect) else NA_character_,
      .groups = "drop"
    )
  gene_summary <- gene_summary %>% left_join(site_agg, by = "gene_id")
}
if (!"n_regions" %in% names(gene_summary)) gene_summary$n_regions <- 0L
if (!"n_sites" %in% names(gene_summary)) gene_summary$n_sites <- 0L
gene_summary$n_regions <- coalesce(as.integer(gene_summary$n_regions), 0L)
gene_summary$n_sites <- coalesce(as.integer(gene_summary$n_sites), 0L)
gene_summary$in_outlier <- gene_summary$n_regions > 0 | gene_summary$n_sites > 0
if ("n_protein_altering" %in% names(gene_summary))
  gene_summary$n_protein_altering <- coalesce(as.integer(gene_summary$n_protein_altering), 0L)

lead <- c("gene_id", "gene_name", "csq_key", "product", "go_terms", "strand",
          "chr", "gene_start", "gene_end", "gene_length", "in_outlier",
          "n_regions", "n_sites", "n_sites_snp", "n_sites_indel", "n_protein_altering",
          "coding_effects", "outlier_stat", "n_windows", "qseqids")
list_last <- intersect(c("blast_gene_ids", "blast_products", "blast_dbs", "qseqids"), names(gene_summary))
gene_summary <- gene_summary %>%
  select(-any_of("seqid")) %>%
  select(any_of(lead), setdiff(names(gene_summary), c(lead, list_last, "seqid")), any_of(list_last)) %>%
  arrange(desc(.data$n_protein_altering), desc(.data$n_sites), .data$chr, .data$gene_start)

if (nrow(gene_regions) > 0) {
  gr_lead <- c("gene_id", "gene_name", "csq_key", "chr", "gene_start", "gene_end",
               "qseqid", "region_chr", "region_start", "region_end", "outlier_stat",
               "n_windows", "overlaps_region", "n_sites", "n_protein_altering", "coding_effects")
  gene_regions <- gene_regions %>%
    select(any_of(gr_lead), setdiff(names(gene_regions), gr_lead)) %>%
    arrange(.data$chr, .data$gene_start, .data$region_start)
}

if (nrow(gene_sites) > 0) {
  gs_lead <- c("gene_id", "gene_name", "csq_key", "gene_chr", "gene_start", "gene_end",
               "assigned_by", "chr", "pos", "ref", "alt", "variant_type",
               "region_category", "coding_effect", "genes", "genes_protein_altering",
               "qseqids", "outlier_stat")
  # gene_chr may be named chr from asm join after collapse
  if (!"gene_chr" %in% names(gene_sites) && "chr" %in% names(gene_sites)) {
    # site chr is the variant chromosome; gene interval is gene_start/end on the same chr
  }
  gene_sites <- gene_sites %>%
    select(any_of(gs_lead), setdiff(names(gene_sites), gs_lead)) %>%
    arrange(.data$gene_id, .data$chr, .data$pos)
}

out_summary <- file.path(opts$`output-dir`, paste0("outlier_gene_summary_", suffix, ".csv"))
out_regions <- file.path(opts$`output-dir`, paste0("outlier_gene_regions_", suffix, ".csv"))
out_sites <- file.path(opts$`output-dir`, paste0("outlier_gene_sites_", suffix, ".csv"))
write_csv(gene_summary, out_summary, na = "")
write_csv(gene_regions, out_regions, na = "")
write_csv(gene_sites, out_sites, na = "")
message("Wrote ", out_summary, " (", nrow(gene_summary), " genes)")
message("Wrote ", out_regions, " (", nrow(gene_regions), " gene-region rows)")
message("Wrote ", out_sites, " (", nrow(gene_sites), " gene-site rows)")
message("Done.")
