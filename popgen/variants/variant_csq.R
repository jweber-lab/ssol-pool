#!/usr/bin/env Rscript

###############################################################################
# variant_csq.R
#
# Parse bcftools csq BCSQ INFO, keep raw string, add csq_consequences (verbatim
# tokens), region_category + coding_effect (mapped for plotting), genes and
# genes_protein_altering (gene names from BCSQ), variant_type (snp/indel/other),
# optional long table (one row per BCSQ sub-record).
#
# Sibling: variant_csq.sh (runs bcftools csq + query, then this script).
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
  library(tibble)
})

#' Split all BCSQ strings into one table of records (comma-separated; one per transcript/ALT).
#' Each record: consequence|gene|transcript|biotype|strand|...; the consequence
#' field can join several terms with "&" (e.g. stop_gained&frameshift).
#' Records starting with "@" point to a compound record at another position.
#' Returns list(rec = records with .row index, term = terms with .rec index).
parse_bcsq <- function(bcsq) {
  b <- as.character(bcsq)
  has <- which(!is.na(b) & nzchar(b) & !b %in% c(".", "NA"))
  parts <- strsplit(b[has], ",", fixed = TRUE)
  row <- rep(has, lengths(parts))
  piece <- trimws(unlist(parts, use.names = FALSE))
  keep <- nzchar(piece) & !startsWith(piece, "@")
  row <- row[keep]; piece <- piece[keep]
  fld <- strsplit(piece, "|", fixed = TRUE)
  field <- function(k) vapply(fld, `[`, character(1L), k)
  cons <- trimws(sub("^[[:space:]]*\\*", "", field(1L))) # leading * = downstream of stop
  gene <- field(2L)
  ok <- !is.na(cons) & nzchar(cons)
  rec <- tibble(
    .row = row[ok], consequence = cons[ok],
    gene = ifelse(!is.na(gene[ok]) & nzchar(gene[ok]), gene[ok], NA_character_),
    transcript = field(3L)[ok], strand = field(5L)[ok]
  )
  tl <- strsplit(rec$consequence, "&", fixed = TRUE)
  term <- tibble(.rec = rep(seq_len(nrow(rec)), lengths(tl)), term = trimws(unlist(tl, use.names = FALSE)))
  term <- term[nzchar(term$term), ]
  list(rec = rec, term = term)
}

# coding_effect values that change the protein (genes_protein_altering)
PROTEIN_ALTERING <- c("missense", "nonsense", "frameshift", "inframe_indel", "other_coding")

variant_type_of <- function(ref, alt) {
  alt1 <- sub(",.*", "", alt)
  case_when(
    is.na(ref) | is.na(alt1) ~ NA_character_,
    nchar(ref) == 1L & nchar(alt1) == 1L ~ "snp",
    nchar(ref) != nchar(alt1) ~ "indel",
    TRUE ~ "other"
  )
}

#' Per group (integer id in 1..n): unique non-empty values in first-seen order, joined by ";"
join_unique_by <- function(group, x, n) {
  out <- rep(NA_character_, n)
  ok <- !is.na(x) & nzchar(x)
  group <- group[ok]; x <- x[ok]
  u <- !duplicated(paste(group, x, sep = "\r"))
  if (!any(u)) return(out)
  s <- split(x[u], group[u])
  out[as.integer(names(s))] <- vapply(s, paste, character(1L), collapse = ";")
  out
}

#' Per group (integer id in 1..n): mapped categories of the most severe term (first on ties)
pick_by <- function(group, region, effect, severity, n) {
  out_region <- rep("none", n); out_effect <- rep("none", n)
  o <- order(group, -severity, seq_along(group))
  first <- o[!duplicated(group[o])]
  out_region[group[first]] <- region[first]
  out_effect[group[first]] <- effect[first]
  list(region_category = out_region, coding_effect = out_effect)
}

# Default mapping: bcftools csq consequence token -> region_category, coding_effect, severity (higher = more severe for wide-row pick)
default_csq_map <- tribble(
  ~consequence, ~region_category, ~coding_effect, ~severity,
  "missense", "cds", "missense", 85L,
  "missense_variant", "cds", "missense", 85L,
  "synonymous", "cds", "synonymous", 55L,
  "synonymous_variant", "cds", "synonymous", 55L,
  "stop_gained", "cds", "nonsense", 100L,
  "stop_lost", "cds", "other_coding", 78L,
  "stop_retained", "cds", "synonymous", 52L,
  "start_lost", "cds", "other_coding", 77L,
  "frameshift", "cds", "frameshift", 95L,
  "frameshift_variant", "cds", "frameshift", 95L,
  "inframe_insertion", "cds", "inframe_indel", 70L,
  "inframe_deletion", "cds", "inframe_indel", 70L,
  "inframe_altering", "cds", "inframe_indel", 70L,
  "feature_elongation", "cds", "other_coding", 75L,
  "coding_sequence", "cds", "other_coding", 60L,
  "start_retained", "cds", "synonymous", 52L,
  "non_coding", "ncRNA", "not_coding", 33L,
  "splice_acceptor", "splice", "not_coding", 62L,
  "splice_donor", "splice", "not_coding", 62L,
  "splice_region", "splice", "not_coding", 58L,
  "5_prime_UTR_variant", "five_prime_utr", "not_coding", 42L,
  "5_prime_utr", "five_prime_utr", "not_coding", 42L,
  "3_prime_UTR_variant", "three_prime_utr", "not_coding", 41L,
  "3_prime_utr", "three_prime_utr", "not_coding", 41L,
  "intron_variant", "intron", "not_coding", 35L,
  "intron", "intron", "not_coding", 35L,
  "intergenic_variant", "intergenic", "not_coding", 12L,
  "intergenic", "intergenic", "not_coding", 12L,
  "upstream_gene_variant", "upstream", "not_coding", 15L,
  "downstream_gene_variant", "downstream", "not_coding", 14L,
  "non_coding_transcript_exon_variant", "ncRNA", "not_coding", 33L,
  "non_coding_transcript_variant", "ncRNA", "not_coding", 32L,
  "mature_miRNA_variant", "ncRNA", "not_coding", 31L,
  "NMD_transcript_variant", "cds", "other_coding", 60L
)

#' Map terms to region_category / coding_effect / severity (unknown tokens -> other / unknown / 5)
map_terms <- function(term, cq_map) {
  m <- match(term, cq_map$consequence)
  list(
    region = ifelse(is.na(m), "other", as.character(cq_map$region_category[m])),
    effect = ifelse(is.na(m), "unknown", as.character(cq_map$coding_effect[m])),
    severity = ifelse(is.na(m), 5L, as.integer(cq_map$severity[m]))
  )
}

sanitize_colname <- function(s) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", s))
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) x <- "na"
  x
}

main <- function() {
  option_list <- list(
    make_option(c("--input-tsv"), type = "character", help = "TSV from bcftools query (must include BCSQ or INFO/BCSQ column)"),
    make_option(c("--output-tsv"), type = "character", help = "Annotated wide site TSV output path"),
    make_option(c("--output-long-tsv"), type = "character", default = NULL, help = "Optional long-format BCSQ rows (one per sub-record)"),
    make_option(c("--map-tsv"), type = "character", default = NULL, help = "Optional TSV with columns consequence, region_category, coding_effect, severity (extends/overrides defaults)")
  )
  opts <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opts$`input-tsv`) || is.null(opts$`output-tsv`)) {
    stop("--input-tsv and --output-tsv are required")
  }
  path <- opts$`input-tsv`
  if (!file.exists(path)) stop("Input not found: ", path)

  d <- read_tab(path) %>% rename_all(tolower)
  bcol <- intersect(c("bcsq", "info/bcsq", "info_bcsq"), names(d))[1L]
  if (is.na(bcol)) stop("No BCSQ column found (expected BCSQ or INFO/BCSQ)")
  names(d)[names(d) == bcol] <- "bcsq"

  cq_map <- default_csq_map
  if (!is.null(opts$`map-tsv`) && file.exists(opts$`map-tsv`)) {
    extra <- readr::read_tsv(opts$`map-tsv`, show_col_types = FALSE)
    req <- c("consequence", "region_category", "coding_effect", "severity")
    if (!all(req %in% names(extra))) stop("--map-tsv must have columns: ", paste(req, collapse = ", "))
    cq_map <- bind_rows(cq_map, extra %>% distinct(.data$consequence, .keep_all = TRUE)) %>%
      group_by(.data$consequence) %>% slice_max(.data$severity, n = 1L, with_ties = FALSE) %>% ungroup()
  }

  n <- nrow(d)
  p <- parse_bcsq(d$bcsq)
  rec <- p$rec; term <- p$term
  mt <- map_terms(term$term, cq_map)
  term_row <- rec$.row[term$.rec]

  # Site: most severe term over all records; record: most severe of its own terms
  site_pick <- pick_by(term_row, mt$region, mt$effect, mt$severity, n)
  rec_pick <- pick_by(term$.rec, mt$region, mt$effect, mt$severity, nrow(rec))
  altering <- rec_pick$coding_effect %in% PROTEIN_ALTERING

  out <- bind_cols(d, tibble(
    region_category = site_pick$region_category,
    coding_effect = site_pick$coding_effect,
    csq_consequences = join_unique_by(term_row, term$term, n),
    genes = join_unique_by(rec$.row, rec$gene, n),
    genes_protein_altering = join_unique_by(rec$.row[altering], rec$gene[altering], n)
  ))
  if (all(c("ref", "alt") %in% names(out))) out$variant_type <- variant_type_of(as.character(out$ref), as.character(out$alt))
  readr::write_tsv(out, opts$`output-tsv`, na = "")
  message("Wrote ", opts$`output-tsv`)

  if (!is.null(opts$`output-long-tsv`)) {
    col <- function(nm, f = as.character, na = NA_character_) if (nm %in% names(d)) suppressWarnings(f(d[[nm]])) else rep(na, n)
    site <- tibble(
      chr = col("chr"), pos = col("pos", as.numeric, NA_real_), ref = col("ref"), alt = col("alt"),
      qual = col("qual", as.numeric, NA_real_),
      filter = if ("filter" %in% names(d)) col("filter") else col("filt"),
      dp = col("dp", as.numeric, NA_real_)
    )
    # one row per BCSQ record; sites without records get a single "none" row
    none <- setdiff(seq_len(n), rec$.row)
    long_all <- bind_rows(
      tibble(.row = rec$.row, bcsq = as.character(d$bcsq[rec$.row]), consequence = rec$consequence, gene = rec$gene,
             transcript = rec$transcript, strand = rec$strand,
             region_category = rec_pick$region_category, coding_effect = rec_pick$coding_effect),
      tibble(.row = none, bcsq = NA_character_, consequence = NA_character_, gene = NA_character_,
             transcript = NA_character_, strand = NA_character_, region_category = "none", coding_effect = "none")
    )
    long_all <- long_all[order(long_all$.row), ]
    long_all <- bind_cols(site[long_all$.row, ], long_all %>% select(-".row"))
    readr::write_tsv(long_all, opts$`output-long-tsv`, na = "")
    message("Wrote ", opts$`output-long-tsv`)
  }
  invisible(NULL)
}

read_tab <- function(path, delim = NULL, na = c("", "NA", "nan", "NaN", "NAN")) {
  first <- readLines(path, n = 1L)
  if (is.null(delim)) {
    delim <- if (grepl("\t", first) && !grepl(",", first)) "\t" else if (grepl(",", first)) "," else "\t"
  }
  if (delim == "\t") readr::read_tsv(path, na = na, show_col_types = FALSE)
  else readr::read_csv(path, na = na, show_col_types = FALSE)
}

if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(e) {
      message(conditionMessage(e))
      quit(status = 1L)
    }
  )
}

