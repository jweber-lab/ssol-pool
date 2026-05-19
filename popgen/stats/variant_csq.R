#!/usr/bin/env Rscript

###############################################################################
# variant_csq.R
#
# Parse bcftools csq BCSQ INFO, keep raw string, add csq_consequences (verbatim
# tokens), region_category + coding_effect (mapped for plotting), optional long
# table (one row per BCSQ sub-record).
#
# Sibling: variant_csq.sh (runs bcftools csq + query, then this script).
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
  library(tibble)
})

#' Split BCSQ into records (comma between ALTs; comma can also separate
#' consequences — bcftools uses comma to separate consequence entries).
#' Each entry: consequence|gene|strand-aware fields...
parse_bcsq_string <- function(bcsq) {
  if (length(bcsq) != 1L) return(list())
  s <- bcsq[[1L]]
  if (is.na(s) || !nzchar(s) || s %in% c(".", "NA")) return(list())
  # Multiple consequences: comma-separated; each part has pipes
  parts <- strsplit(s, ",", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  out <- list()
  for (p in parts) {
    fld <- strsplit(p, "|", fixed = TRUE)[[1L]]
    if (length(fld) < 1L) next
    cons <- sub("^[[:space:]]*\\*", "", fld[[1L]]) # leading * = downstream of stop
    cons <- trimws(cons)
    if (!nzchar(cons)) next
    out[[length(out) + 1L]] <- list(
      consequence = cons,
      raw_piece = p,
      fields = fld
    )
  }
  out
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

map_consequence <- function(tokens, cq_map) {
  if (length(tokens) == 0L) {
    return(c(region_category = "none", coding_effect = "none", severity = 0L))
  }
  rows <- lapply(tokens, function(tok) {
    m <- cq_map %>% filter(.data$consequence == tok)
    if (nrow(m) >= 1L) {
      m[1L, ]
    } else {
      tibble(
        consequence = tok, region_category = "other", coding_effect = "unknown",
        severity = 5L
      )
    }
  })
  rows <- bind_rows(rows)
  win <- rows %>% slice_max(.data$severity, n = 1L, with_ties = FALSE)
  c(
    region_category = as.character(win$region_category[1L]),
    coding_effect = as.character(win$coding_effect[1L]),
    severity = as.integer(win$severity[1L])
  )
}

#' Build long table from parsed rows + original variant context
bcsq_to_long <- function(chr, pos, ref, alt, qual, filter_val, dp, bcsq_raw, records, cq_map) {
  if (length(records) == 0L) {
    return(tibble(
      chr = chr, pos = pos, ref = ref, alt = alt, qual = qual, filter = filter_val, dp = dp,
      bcsq = NA_character_, consequence = NA_character_, gene = NA_character_,
      transcript = NA_character_, strand = NA_character_,
      region_category = "none", coding_effect = "none"
    ))
  }
  parse_one <- function(rec) {
    fld <- rec$fields
    gene <- if (length(fld) >= 2L) fld[[2L]] else NA_character_
    # transcript id often 3rd; biotype 4th; strand 5th — layout from bcftools csq
    transcript <- if (length(fld) >= 3L) fld[[3L]] else NA_character_
    strand <- if (length(fld) >= 5L) fld[[5L]] else NA_character_
    tok <- rec$consequence
    m <- cq_map %>% filter(.data$consequence == tok)
    if (nrow(m) == 0L) {
      rc <- "other"
      ce <- "unknown"
    } else {
      rc <- m$region_category[1L]
      ce <- m$coding_effect[1L]
    }
    tibble(
      chr = chr, pos = pos, ref = ref, alt = alt, qual = qual, filter = filter_val, dp = dp,
      bcsq = bcsq_raw, consequence = tok, gene = gene, transcript = transcript, strand = strand,
      region_category = rc, coding_effect = ce
    )
  }
  bind_rows(lapply(records, parse_one))
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

  wide_rows <- vector("list", nrow(d))
  long_parts <- vector("list", nrow(d))

  for (i in seq_len(nrow(d))) {
    recs <- parse_bcsq_string(d$bcsq[i])
    toks <- vapply(recs, function(z) z$consequence, character(1L))
    toks <- unique(toks[nzchar(toks)])
    picked <- map_consequence(toks, cq_map)
    csq_consequences <- if (length(toks) == 0L) NA_character_ else paste(toks, collapse = ";")
    wide_rows[[i]] <- tibble(
      region_category = picked[["region_category"]],
      coding_effect = picked[["coding_effect"]],
      csq_consequences = csq_consequences
    )
    chr_i <- if ("chr" %in% names(d)) as.character(d$chr[i]) else NA_character_
    pos_i <- if ("pos" %in% names(d)) as.numeric(d$pos[i]) else NA_real_
    ref_i <- if ("ref" %in% names(d)) as.character(d$ref[i]) else NA_character_
    alt_i <- if ("alt" %in% names(d)) as.character(d$alt[i]) else NA_character_
    qual_i <- if ("qual" %in% names(d)) as.numeric(d$qual[i]) else NA_real_
    filt_i <- if ("filter" %in% names(d)) as.character(d$filter[i]) else if ("filt" %in% names(d)) as.character(d$filt[i]) else NA_character_
    dp_i <- if ("dp" %in% names(d)) as.numeric(d$dp[i]) else NA_real_
    long_parts[[i]] <- bcsq_to_long(chr_i, pos_i, ref_i, alt_i, qual_i, filt_i, dp_i, d$bcsq[i], recs, cq_map)
  }

  wide_add <- bind_rows(wide_rows)
  out <- bind_cols(d, wide_add)
  readr::write_tsv(out, opts$`output-tsv`, na = "")
  message("Wrote ", opts$`output-tsv`)

  if (!is.null(opts$`output-long-tsv`)) {
    long_all <- bind_rows(long_parts)
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

