#!/usr/bin/env Rscript

###############################################################################
# gene_report.R
#
# Collect genes from one or more gene tables into one output directory.
# A gene whose folder is already there is left in place, so later lists add
# only the genes that are new. Two summaries in the output directory are
# rewritten every run from every gene in the directory:
#
#   gene_summary_brief.csv  id, chromosome, strand, coordinates, variable SNP
#                           and indel counts, and the stats that pass the
#                           0.99 and 0.999 high and low tails
#   gene_summary_full.csv   those fields, plus per-stat pass counts, the
#                           maximum overlapping value and quantile, and
#                           variant-consequence counts
#
# Stat values and quantiles come from the single-site and windowed HDF5 files,
# not from the outlier table that supplied the gene. Each new gene folder holds
# the genomic sequence, the spliced CDS when a GFF is given, the windows and
# sites from the gene table, and the proteome blastx hits.
###############################################################################

GENE_COLLATION_FUNCS_ONLY <- TRUE
cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
if (length(script_path) == 0 || !nzchar(script_path)) script_path <- "gene_report.R"
source(file.path(dirname(normalizePath(script_path)), "..", "outliers", "gene_collation.R"), local = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
})

sibling_table <- function(genes_path, from, to) {
  bn <- basename(genes_path)
  if (!grepl(from, bn, fixed = TRUE)) return(NULL)
  cand <- file.path(dirname(genes_path), sub(from, to, bn, fixed = TRUE))
  if (file.exists(cand)) cand else NULL
}

safe_gene_dir <- function(gene_id, used) {
  s <- gsub("[^A-Za-z0-9._-]+", "_", gene_id)
  s <- gsub("^_+|_+$", "", s)
  if (!nzchar(s)) s <- "gene"
  base <- s
  n <- 2L
  while (s %in% used) {
    s <- paste0(base, "_", n)
    n <- n + 1L
  }
  s
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

scale_from_name <- function(bn) {
  if (!grepl("_w[0-9]+", bn)) return("single")
  w <- sub(".*_w([0-9]+).*", "\\1", bn)
  if (!grepl("_s[0-9]+", bn)) return(paste0("w", w))
  paste0("w", w, "_s", sub(".*_s([0-9]+).*", "\\1", bn))
}

strip_sample <- function(sample_name) {
  vapply(as.character(sample_name), function(s) {
    if (grepl("^diversity_w[0-9]+_s[0-9]+_sample", s, ignore.case = TRUE)) {
      stripped <- sub("^diversity_w[0-9]+_s[0-9]+_sample", "", s, ignore.case = TRUE)
      if (nzchar(trimws(stripped))) return(stripped)
    }
    s
  }, character(1), USE.NAMES = FALSE)
}

normalize_wide_stat <- function(col) {
  if (grepl("\\.(fst|pbe|qibe|lrpi|pi|theta|tajima_d)$", col)) return(col)
  if (!grepl("^(fst|pbe|qibe|lrpi|pi|theta|tajima_d)_.+", col)) return(col)
  stat <- sub("^(fst|pbe|qibe|lrpi|pi|theta|tajima_d)_.+$", "\\1", col)
  rest <- sub("^(fst|pbe|qibe|lrpi|pi|theta|tajima_d)_", "", col)
  paste0(gsub("_", ":", rest, fixed = TRUE), ".", stat)
}

empty_stats <- function() {
  data.frame(
    scale = character(), stat = character(),
    n_high_0.99 = integer(), n_high_0.999 = integer(),
    n_low_0.99 = integer(), n_low_0.999 = integer(),
    max_value = numeric(), max_quantile = numeric(),
    stringsAsFactors = FALSE
  )
}

append_scores <- function(bucket, gene_id, scale, stat, value, quant) {
  q <- as.numeric(quant)
  v <- as.numeric(value)
  row <- data.frame(
    scale = scale, stat = stat,
    n_high_0.99 = sum(q >= 0.99, na.rm = TRUE),
    n_high_0.999 = sum(q >= 0.999, na.rm = TRUE),
    n_low_0.99 = sum(q <= 0.01, na.rm = TRUE),
    n_low_0.999 = sum(q <= 0.001, na.rm = TRUE),
    max_value = safe_max(v), max_quantile = safe_max(q),
    stringsAsFactors = FALSE
  )
  prev <- bucket[[gene_id]]
  bucket[[gene_id]] <- if (is.null(prev)) row else rbind(prev, row)
}

consolidate_stats <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(empty_stats())
  df %>%
    group_by(.data$scale, .data$stat) %>%
    summarise(
      n_high_0.99 = sum(.data$n_high_0.99),
      n_high_0.999 = sum(.data$n_high_0.999),
      n_low_0.99 = sum(.data$n_low_0.99),
      n_low_0.999 = sum(.data$n_low_0.999),
      max_value = safe_max(.data$max_value),
      max_quantile = safe_max(.data$max_quantile),
      .groups = "drop"
    ) %>%
    as.data.frame()
}

score_rows <- function(bucket, genes, genes_by_chr, chr, start, end, value, quant, stat, scale) {
  by_chr <- split(seq_along(chr), chr)
  for (ch in names(genes_by_chr)) {
    rows <- by_chr[[ch]]
    if (is.null(rows)) next
    st <- start[rows]
    en <- end[rows]
    vv <- value[rows]
    qq <- quant[rows]
    ss <- stat[rows]
    for (gi in genes_by_chr[[ch]]) {
      hit <- !is.na(en) & !is.na(st) & en >= genes$gene_start[[gi]] & st <= genes$gene_end[[gi]]
      if (!any(hit)) next
      parts <- split(seq_len(sum(hit)), ss[hit])
      for (nm in names(parts)) {
        ii <- which(hit)[parts[[nm]]]
        append_scores(bucket, genes$gene_id[[gi]], scale, nm, vv[ii], qq[ii])
      }
    }
  }
}

pick_h5_group <- function(h5, scale) {
  groups <- names(h5)
  if (identical(scale, "single") && "sites" %in% groups) return("sites")
  if ("windows" %in% groups) return("windows")
  if ("sites" %in% groups) return("sites")
  NULL
}

h5_col <- function(grp, name, idx) {
  if (!name %in% names(grp)) return(NULL)
  grp[[name]][][idx]
}

score_h5_file <- function(path, genes, genes_by_chr, bucket) {
  scale <- scale_from_name(basename(path))
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  grp_name <- pick_h5_group(h5, scale)
  if (is.null(grp_name)) {
    message("Skipping ", basename(path), ": no /windows or /sites group")
    return(invisible(NULL))
  }
  message("Reading ", basename(path), " (", scale, ")")
  grp <- h5[[grp_name]]
  gnames <- names(grp)
  if (!"chr" %in% gnames) return(invisible(NULL))
  chr_all <- as.character(grp[["chr"]][])
  n <- length(chr_all)
  start_all <- if ("start" %in% gnames) as.numeric(grp[["start"]][]) else if ("pos" %in% gnames) as.numeric(grp[["pos"]][]) else NULL
  if (is.null(start_all)) return(invisible(NULL))
  end_all <- if ("end" %in% gnames) as.numeric(grp[["end"]][]) else start_all
  keep <- !is.na(chr_all) & chr_all %in% genes$chr & !is.na(start_all) & !is.na(end_all)
  if (!any(keep)) return(invisible(NULL))
  idx <- which(keep)
  chr <- chr_all[idx]
  start <- start_all[idx]
  end <- end_all[idx]
  quant_or_na <- function(qcol) {
    if (qcol %in% gnames) as.numeric(h5_col(grp, qcol, idx)) else rep(NA_real_, length(idx))
  }

  if ("sample" %in% gnames && any(c("pi", "theta", "tajima_d") %in% gnames)) {
    sample <- strip_sample(h5_col(grp, "sample", idx))
    for (col in c("pi", "theta", "tajima_d")) {
      if (!col %in% gnames) next
      score_rows(bucket, genes, genes_by_chr, chr, start, end,
                 as.numeric(h5_col(grp, col, idx)), quant_or_na(paste0(col, "_quantile")),
                 paste0(sample, ".", col), scale)
    }
    return(invisible(NULL))
  }

  branch <- intersect(c("fst", "pbe", "qibe", "lrpi"), gnames)
  if ("pop1" %in% gnames && length(branch) > 0) {
    pop_vals <- list(pop1 = as.character(h5_col(grp, "pop1", idx)))
    for (pc in c("pop2", "pop3", "pop4")) {
      if (pc %in% gnames) pop_vals[[pc]] <- as.character(h5_col(grp, pc, idx))
    }
    for (col in branch) {
      np <- switch(col, fst = 2L, lrpi = 2L, pbe = 3L, qibe = 4L, 2L)
      if (length(pop_vals) < np) next
      id <- do.call(paste, c(pop_vals[seq_len(np)], list(sep = ":")))
      score_rows(bucket, genes, genes_by_chr, chr, start, end,
                 as.numeric(h5_col(grp, col, idx)), quant_or_na(paste0(col, "_quantile")),
                 paste0(id, ".", col), scale)
    }
    return(invisible(NULL))
  }

  wide_as_stored <- basename(path) == "single_position.h5"
  qnames <- grep("_quantile$", gnames, value = TRUE)
  for (qcol in qnames) {
    stat_col <- sub("_quantile$", "", qcol)
    if (!stat_col %in% gnames || grepl("_rank$", stat_col)) next
    value <- suppressWarnings(as.numeric(h5_col(grp, stat_col, idx)))
    if (length(value) == 0 || all(is.na(value))) next
    stat_name <- if (wide_as_stored) stat_col else normalize_wide_stat(stat_col)
    score_rows(bucket, genes, genes_by_chr, chr, start, end, value, quant_or_na(qcol),
               rep(stat_name, length(idx)), scale)
  }
  invisible(NULL)
}

stat_h5_files <- function(hdf5_dir) {
  files <- list.files(hdf5_dir, pattern = "\\.h5$", full.names = TRUE)
  if (length(files) == 0) return(character(0))
  bn <- basename(files)
  keep <- grepl("^(diversity|fst|pbe|qibe|lrpi).+\\.h5$", bn) | bn == "single_position.h5"
  if ("single_position.h5" %in% bn) {
    keep <- keep & !bn %in% c("fst_single.h5", "pbe_single.h5", "qibe_single.h5", "lrpi_single.h5")
  }
  files[keep]
}

variant_h5_path <- function(hdf5_dir) {
  merged <- file.path(hdf5_dir, "single_position.h5")
  if (file.exists(merged)) return(merged)
  sites <- file.path(hdf5_dir, "variants.h5")
  if (file.exists(sites)) return(sites)
  NULL
}

count_split_field <- function(x) {
  parts <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  parts <- parts[!is.na(parts) & nzchar(parts) & parts != "none"]
  if (length(parts) == 0) return(NA_character_)
  tb <- sort(table(parts), decreasing = TRUE)
  paste0(names(tb), ":", as.integer(tb), collapse = ";")
}

variant_counts <- function(genes, d) {
  out <- lapply(seq_len(nrow(genes)), function(i) {
    empty <- list(n_snps = NA_integer_, n_indels = NA_integer_, n_sites = NA_integer_,
                  n_protein_altering = NA_integer_, coding_effects = NA_character_,
                  region_categories = NA_character_, csq_consequences = NA_character_)
    if (is.null(d)) return(empty)
    known_zero <- function() {
      empty$n_sites <- 0L
      if ("variant_type" %in% names(d)) {
        empty$n_snps <- 0L
        empty$n_indels <- 0L
      }
      if ("genes_protein_altering" %in% names(d)) empty$n_protein_altering <- 0L
      empty
    }
    if (nrow(d) == 0) return(known_zero())
    hit <- d$chr == genes$chr[[i]] & d$pos >= genes$gene_start[[i]] & d$pos <= genes$gene_end[[i]]
    hit <- hit & !is.na(hit)
    sub <- d[hit, , drop = FALSE]
    if (nrow(sub) == 0) return(known_zero())
    empty$n_sites <- nrow(sub)
    if ("variant_type" %in% names(sub)) {
      empty$n_snps <- sum(sub$variant_type == "snp", na.rm = TRUE)
      empty$n_indels <- sum(sub$variant_type == "indel", na.rm = TRUE)
    }
    if ("genes_protein_altering" %in% names(sub)) {
      empty$n_protein_altering <- sum(!is.na(sub$genes_protein_altering) & nzchar(sub$genes_protein_altering))
    }
    if ("coding_effect" %in% names(sub)) empty$coding_effects <- effect_counts_string(sub$coding_effect)
    if ("region_category" %in% names(sub)) empty$region_categories <- count_split_field(sub$region_category)
    if ("csq_consequences" %in% names(sub)) empty$csq_consequences <- count_split_field(sub$csq_consequences)
    empty
  })
  bind_rows(out)
}

read_variant_h5 <- function(path, genes) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  grp_name <- if ("sites" %in% names(h5)) "sites" else if ("windows" %in% names(h5)) "windows" else NULL
  if (is.null(grp_name)) return(NULL)
  grp <- h5[[grp_name]]
  gnames <- names(grp)
  if (!"chr" %in% gnames) return(NULL)
  if (!any(c("variant_type", "coding_effect", "region_category", "csq_consequences", "genes_protein_altering") %in% gnames))
    return(NULL)
  pos <- if ("pos" %in% gnames) as.numeric(grp[["pos"]][]) else if ("start" %in% gnames) as.numeric(grp[["start"]][]) else NULL
  if (is.null(pos)) return(NULL)
  chr <- as.character(grp[["chr"]][])
  keep <- !is.na(chr) & chr %in% genes$chr & !is.na(pos)
  if (!any(keep)) keep <- rep(FALSE, length(chr))
  d <- data.frame(chr = chr[keep], pos = pos[keep], stringsAsFactors = FALSE)
  for (col in c("variant_type", "coding_effect", "region_category", "csq_consequences", "genes_protein_altering")) {
    if (col %in% gnames) d[[col]] <- as.character(grp[[col]][][keep])
  }
  d
}

h5_stamp <- function(hdf5_dir) {
  files <- list.files(hdf5_dir, pattern = "\\.h5$", full.names = TRUE)
  if (length(files) == 0) return(character(0))
  info <- file.info(files)
  stamp <- setNames(as.character(info$mtime), basename(files))
  stamp[order(names(stamp))]
}

first_col <- function(row, candidates) {
  for (col in candidates) {
    if (col %in% names(row) && length(row[[col]]) > 0 && !is.na(row[[col]][[1]]) && nzchar(as.character(row[[col]][[1]]))) {
      return(row[[col]][[1]])
    }
  }
  NA_character_
}

int_col <- function(row, candidates) {
  for (col in candidates) {
    if (!col %in% names(row)) next
    v <- suppressWarnings(as.integer(row[[col]][[1]]))
    if (length(v) == 1 && !is.na(v)) return(v)
  }
  NA_integer_
}

sites_for_counts <- function(sites_sub) {
  if (is.null(sites_sub) || nrow(sites_sub) == 0 || !"chr" %in% names(sites_sub)) return(NULL)
  d <- sites_sub
  if (!"pos" %in% names(d)) {
    if (!"start" %in% names(d)) return(NULL)
    d$pos <- d$start
  }
  d$chr <- as.character(d$chr)
  d$pos <- as.numeric(d$pos)
  d
}

meta_from_sources <- function(row, variant_row, sites_sub) {
  meta <- data.frame(
    gene_id = as.character(row$gene_id[[1]]),
    chr = as.character(row$chr[[1]]),
    strand = as.character(first_col(row, "strand")),
    gene_start = as.integer(row$gene_start[[1]]),
    gene_end = as.integer(row$gene_end[[1]]),
    n_snps = int_col(row, c("n_sites_snp", "n_snps")),
    n_indels = int_col(row, c("n_sites_indel", "n_indels")),
    n_sites = int_col(row, c("n_sites")),
    n_protein_altering = int_col(row, c("n_protein_altering")),
    coding_effects = as.character(first_col(row, "coding_effects")),
    region_categories = as.character(first_col(row, c("region_categories", "region_category"))),
    csq_consequences = as.character(first_col(row, c("csq_consequences"))),
    stringsAsFactors = FALSE
  )
  sites_ready <- sites_for_counts(sites_sub)
  if (!is.null(sites_ready)) {
    counted <- variant_counts(
      data.frame(chr = meta$chr, gene_start = meta$gene_start, gene_end = meta$gene_end, stringsAsFactors = FALSE),
      sites_ready
    )
    for (col in names(counted)) {
      if (!is.na(counted[[col]][[1]])) meta[[col]] <- counted[[col]][[1]]
    }
  }
  if (!is.null(variant_row)) {
    for (col in c("n_snps", "n_indels", "n_sites", "n_protein_altering", "coding_effects", "region_categories", "csq_consequences")) {
      if (!is.null(variant_row[[col]]) && length(variant_row[[col]]) > 0 && !is.na(variant_row[[col]][[1]])) {
        meta[[col]] <- variant_row[[col]][[1]]
      }
    }
  }
  meta
}

pass_list <- function(stats, col) {
  if (is.null(stats) || nrow(stats) == 0) return(NA_character_)
  keys <- paste(stats$scale, stats$stat, sep = "|")
  keep <- stats[[col]] > 0
  keep[is.na(keep)] <- FALSE
  keys <- sort(keys[keep])
  if (length(keys) == 0) NA_character_ else paste(keys, collapse = ";")
}

stat_measures <- c("n_high_0.99", "n_high_0.999", "n_low_0.99", "n_low_0.999", "max_value", "max_quantile")

write_summaries <- function(caches, out_dir) {
  if (length(caches) == 0) stop("No genes to summarize")
  metas <- bind_rows(lapply(caches, `[[`, "meta"))
  metas <- metas[order(metas$chr, metas$gene_start, metas$gene_id), , drop = FALSE]
  by_id <- setNames(caches, vapply(caches, function(z) z$meta$gene_id[[1]], ""))
  brief <- data.frame(
    gene_id = metas$gene_id, chr = metas$chr, strand = metas$strand,
    gene_start = metas$gene_start, gene_end = metas$gene_end,
    n_snps = metas$n_snps, n_indels = metas$n_indels,
    stats_high_0.99 = vapply(metas$gene_id, function(id) pass_list(by_id[[id]]$stats, "n_high_0.99"), ""),
    stats_high_0.999 = vapply(metas$gene_id, function(id) pass_list(by_id[[id]]$stats, "n_high_0.999"), ""),
    stats_low_0.99 = vapply(metas$gene_id, function(id) pass_list(by_id[[id]]$stats, "n_low_0.99"), ""),
    stats_low_0.999 = vapply(metas$gene_id, function(id) pass_list(by_id[[id]]$stats, "n_low_0.999"), ""),
    stringsAsFactors = FALSE
  )
  brief_path <- file.path(out_dir, "gene_summary_brief.csv")
  write_csv(brief, brief_path, na = "")
  message("Wrote ", brief_path, " (", nrow(brief), " genes)")

  keys <- unique(unlist(lapply(caches, function(z) {
    if (is.null(z$stats) || nrow(z$stats) == 0) character(0) else paste(z$stats$scale, z$stats$stat, sep = "|")
  }), use.names = FALSE))
  keys <- sort(keys)
  wide_names <- unlist(lapply(keys, function(k) paste(k, stat_measures, sep = "|")), use.names = FALSE)
  wide <- matrix(NA_real_, nrow = nrow(metas), ncol = length(wide_names),
                 dimnames = list(metas$gene_id, wide_names))
  for (i in seq_len(nrow(metas))) {
    st <- by_id[[metas$gene_id[[i]]]]$stats
    if (is.null(st) || nrow(st) == 0) next
    for (r in seq_len(nrow(st))) {
      key <- paste(st$scale[[r]], st$stat[[r]], sep = "|")
      for (m in stat_measures) wide[i, paste(key, m, sep = "|")] <- st[[m]][[r]]
    }
  }
  id_cols <- metas[, c("gene_id", "chr", "strand", "gene_start", "gene_end", "n_snps", "n_indels",
                       "n_sites", "n_protein_altering", "coding_effects", "region_categories", "csq_consequences"),
                   drop = FALSE]
  full <- if (length(wide_names) == 0) {
    id_cols
  } else {
    cbind(id_cols, as.data.frame(wide, optional = TRUE), stringsAsFactors = FALSE)
  }
  full_path <- file.path(out_dir, "gene_summary_full.csv")
  write_csv(full, full_path, na = "")
  message("Wrote ", full_path, " (", nrow(full), " genes, ", length(keys), " stats)")
  invisible(list(brief = brief, full = full))
}

existing_gene_folders <- function(out_dir) {
  if (!dir.exists(out_dir)) return(list(ids = character(0), folders = character(0)))
  ids <- character(0)
  folders <- character(0)
  for (d in list.dirs(out_dir, recursive = FALSE, full.names = TRUE)) {
    bn <- basename(d)
    if (bn %in% c("log", "stat_cache")) next
    id <- NA_character_
    id_file <- file.path(d, "gene_id.txt")
    row_file <- file.path(d, "gene_summary_row.csv")
    if (file.exists(id_file)) {
      id <- readLines(id_file, n = 1L, warn = FALSE)
    } else if (file.exists(row_file)) {
      row <- tryCatch(read_csv(row_file, n_max = 1L, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(row) && "gene_id" %in% names(row)) id <- as.character(row$gene_id[[1]])
    }
    if (is.na(id) || !nzchar(id)) next
    ids <- c(ids, id)
    folders <- c(folders, bn)
  }
  if (any(duplicated(ids))) message("Repeated gene ids in the output directory; keeping the first folder")
  keep <- !duplicated(ids)
  list(ids = ids[keep], folders = setNames(folders[keep], ids[keep]))
}

load_cache <- function(path) {
  obj <- readRDS(path)
  if (is.null(obj$stats)) obj$stats <- empty_stats()
  obj$stats <- as.data.frame(obj$stats, stringsAsFactors = FALSE)
  obj
}

load_caches <- function(cache_dir) {
  files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)
  files <- files[basename(files) != "hdf5_stamp.rds"]
  lapply(files, load_cache)
}

write_span_fasta <- function(chr, start, end, header, reference, samtools, out_fa) {
  reg <- paste0(chr, ":", fmt_bp(start), "-", fmt_bp(end))
  seq <- system2(samtools, c("faidx", reference, reg), stdout = TRUE, stderr = FALSE)
  if (!is.null(attr(seq, "status")) && attr(seq, "status") != 0L) return(FALSE)
  hdr <- which(startsWith(seq, ">"))
  if (length(hdr) != 1L) return(FALSE)
  seq[hdr] <- paste0(">", header, " ", reg)
  writeLines(seq, out_fa)
  TRUE
}

annotate_hit_table <- function(b, qmap, lookup) {
  if (nrow(b) == 0) return(b)
  parsed <- parse_protein_title(b$sseqid, b$stitle)
  gene <- product <- go <- character(nrow(b))
  for (j in seq_len(nrow(b))) {
    rec <- lookup_proteome(lookup, b$sseqid[[j]], parsed$accession[[j]])
    gene[[j]] <- fill_from_lookup(parsed$gene[[j]], if (is.null(rec)) NA_character_ else rec$gene)
    product[[j]] <- fill_from_lookup(parsed$product[[j]], if (is.null(rec)) NA_character_ else rec$product)
    go[[j]] <- fill_from_lookup(NA_character_, if (is.null(rec)) NA_character_ else rec$go)
  }
  b$gene_id <- qmap$gene_id[match(b$qseqid, qmap$qseqid)]
  b$gene <- gene
  b$product <- product
  b$go <- go
  b$accession <- parsed$accession
  b
}

blast_proteome_details <- function(genes, gff_path, reference, config_path, out_dir, folders, samtools, blastx_cmd,
                                   threads, evalue, max_target_seqs, max_hsps) {
  config <- read_blast_config(config_path)
  query_fa <- file.path(out_dir, "genes_cds_query.fa")
  qmap <- write_cds_fasta(genes, gff_path, reference, samtools, query_fa)
  if (is.null(qmap) || nrow(qmap) == 0) {
    message("Proteome blastx: no CDS features for these genes")
    return(invisible(NULL))
  }
  message("Proteome query: ", nrow(qmap), " spliced CDS")
  for (i in seq_along(config)) {
    db <- config[[i]]
    db_path <- db_field(db, "db_path")
    name <- safe_db_name(coalesce_null(db_field(db, "name"), paste0("proteome", i)), i)
    if (name == "proteome") name <- "proteome_db"
    if (is.null(db_path)) {
      message("Skipping proteome ", name, ": db_path missing")
      next
    }
    db_files <- paste0(db_path, c(".phr", ".pin", ".psq", ".pdb", ".pal"))
    if (!any(file.exists(db_files))) {
      message("Skipping proteome ", name, ": BLAST database not found at ", db_path)
      next
    }
    out_tsv <- file.path(out_dir, paste0("blast_genes_", name, ".tsv"))
    message("blastx ", nrow(qmap), " CDS against ", name, " (max_target_seqs=", search_target_seqs(max_target_seqs),
            ", max_hsps=", max_hsps, ", evalue=", evalue, ")")
    ok <- run_gene_blast(query_fa, db_path, out_tsv, blastx_cmd, threads, evalue, max_target_seqs, max_hsps,
                         blast_task = NULL,
                         outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle")
    empty_hits <- tibble(
      qseqid = character(), sseqid = character(), pident = numeric(), length = integer(),
      qstart = integer(), qend = integer(), sstart = integer(), send = integer(),
      evalue = character(), bitscore = numeric(), stitle = character(),
      gene_id = character(), gene = character(), product = character(), go = character(),
      accession = character()
    )
    hits <- if (!ok || !file.exists(out_tsv) || file.info(out_tsv)$size == 0) {
      empty_hits
    } else {
      read_tsv(out_tsv,
               col_names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"),
               show_col_types = FALSE)
    }
    lookup <- load_proteome_lookup(db_field(db, "annotation_tsv", "annotation_table"),
                                   db_field(db, "id_map", "id_map_tsv"), name)
    if (nrow(hits) > 0) {
      hits$evalue_num <- suppressWarnings(as.numeric(hits$evalue))
      hits$bitscore <- suppressWarnings(as.numeric(hits$bitscore))
      hits <- hits[order(hits$evalue_num, -hits$bitscore), , drop = FALSE]
      hits <- annotate_hit_table(hits, qmap, lookup)
    }
    for (j in seq_len(nrow(qmap))) {
      gid <- qmap$gene_id[[j]]
      dest <- file.path(folders[[gid]], paste0("blast_", name, ".csv"))
      sub <- hits[hits$gene_id == gid, , drop = FALSE]
      if (nrow(sub) > 0) {
        sub$qseqid <- gid
        sub <- sub %>% select("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
                              "evalue", "bitscore", "accession", "gene", "product", "go", "stitle")
      }
      write_csv(sub, dest, na = "")
    }
    missed <- setdiff(genes$gene_id, qmap$gene_id)
    for (gid in missed) {
      write_csv(empty_hits[0, c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
                                "evalue", "bitscore", "accession", "gene", "product", "go", "stitle")],
                file.path(folders[[gid]], paste0("blast_", name, ".csv")), na = "")
    }
  }
  invisible(NULL)
}

option_list <- list(
  make_option(c("--genes"), type = "character", default = NULL,
              help = "Gene table (gene_id, chr, gene_start, gene_end). Repeat to combine lists. Genes already in --output-dir are kept"),
  make_option(c("--hdf5-dir"), type = "character", default = NULL,
              help = "Directory of collated HDF5 files. Every diversity, fst, pbe, qibe, lrpi, and single_position file is scored"),
  make_option(c("--regions"), type = "character", default = NULL,
              help = "outlier_gene_regions_*.csv. Default: each summary's sibling regions file, copied into new gene folders"),
  make_option(c("--sites"), type = "character", default = NULL,
              help = "outlier_gene_sites_*.csv. Default: each summary's sibling sites file, copied into new gene folders"),
  make_option(c("--output-dir"), type = "character", default = "gene_report", help = "Output directory [default: %default]"),
  make_option(c("--reference"), type = "character", default = NULL, help = "Assembly FASTA (indexed). Writes gene.fa and, with --gff, cds.fa"),
  make_option(c("--gff"), type = "character", default = NULL, help = "Assembly GFF3, required to extract CDS and to run proteome blastx"),
  make_option(c("--proteome-config"), type = "character", default = NULL,
              help = "YAML/JSON protein databases. blastx of the spliced CDS; every returned HSP is written per gene"),
  make_option(c("--samtools"), type = "character", default = "samtools", help = "samtools [default: %default]"),
  make_option(c("--blastx-cmd"), type = "character", default = "blastx", help = "blastx program [default: %default]"),
  make_option(c("--threads"), type = "integer", default = 4L, help = "Threads per blastx [default: %default]"),
  make_option(c("--evalue"), type = "double", default = 1e-5, help = "blastx -evalue [default: %default]"),
  make_option(c("--max-target-seqs"), type = "integer", default = 25L,
              help = "blastx -max_target_seqs [default: %default]. Values below 5 are raised so the best hit is not dropped"),
  make_option(c("--max-hsps"), type = "integer", default = 5L, help = "blastx -max_hsps [default: %default]"),
  make_option(c("--verbose"), action = "store_true", default = FALSE)
)
raw_args <- commandArgs(trailingOnly = TRUE)
gene_files <- character(0)
keep_arg <- rep(TRUE, length(raw_args))
i <- 1L
while (i <= length(raw_args)) {
  arg <- raw_args[[i]]
  if (arg == "--genes" || startsWith(arg, "--genes=")) {
    if (arg == "--genes") {
      if (i == length(raw_args)) stop("--genes needs a file")
      gene_files <- c(gene_files, raw_args[[i + 1L]])
      keep_arg[i] <- FALSE
      keep_arg[i + 1L] <- FALSE
      i <- i + 2L
    } else {
      gene_files <- c(gene_files, sub("^--genes=", "", arg))
      keep_arg[i] <- FALSE
      i <- i + 1L
    }
  } else {
    i <- i + 1L
  }
}
opts <- parse_args(OptionParser(
  usage = "usage: %prog --genes outlier_gene_summary.csv [--genes another.csv] [--hdf5-dir DIR] [options]",
  option_list = option_list
), args = raw_args[keep_arg])

if (length(gene_files) == 0) stop("--genes is required")
missing_genes <- gene_files[!file.exists(gene_files)]
if (length(missing_genes) > 0) stop("Gene table not found: ", paste(missing_genes, collapse = ", "))
if (!is.null(opts$threads) && opts$threads < 1L) stop("--threads must be >= 1")
if (!is.null(opts$`max-target-seqs`) && opts$`max-target-seqs` < 1L) stop("--max-target-seqs must be >= 1")
if (!is.null(opts$`max-hsps`) && opts$`max-hsps` < 1L) stop("--max-hsps must be >= 1")
if (!is.null(opts$`proteome-config`) && (is.null(opts$reference) || is.null(opts$gff)))
  stop("--reference and --gff are required with --proteome-config")
if (!is.null(opts$reference) && !file.exists(paste0(opts$reference, ".fai")))
  stop("Reference index not found: ", opts$reference, ".fai")
if (!is.null(opts$`hdf5-dir`) && !dir.exists(opts$`hdf5-dir`))
  stop("HDF5 directory not found: ", opts$`hdf5-dir`)
if (!is.null(opts$`hdf5-dir`) && !requireNamespace("hdf5r", quietly = TRUE))
  stop("Package hdf5r required for --hdf5-dir. Install with: install.packages(\"hdf5r\")")

read_gene_file <- function(path) {
  d <- read_table_auto(path)
  if (is.null(d) || nrow(d) == 0) stop("Gene table is empty: ", path)
  for (col in c("gene_id", "chr", "gene_start", "gene_end")) {
    if (!col %in% names(d)) stop("Gene table needs a ", col, " column: ", path)
  }
  d$gene_id <- as.character(d$gene_id)
  d$chr <- as.character(d$chr)
  d$gene_start <- suppressWarnings(as.integer(d$gene_start))
  d$gene_end <- suppressWarnings(as.integer(d$gene_end))
  d$.source <- path
  d[!is.na(d$gene_id) & nzchar(d$gene_id), , drop = FALSE]
}

genes <- bind_rows(lapply(gene_files, read_gene_file))
n_in <- nrow(genes)
genes <- genes[!duplicated(genes$gene_id), , drop = FALSE]
if (nrow(genes) < n_in) message("Kept the first row for ", n_in - nrow(genes), " repeated gene ids")
bad <- is.na(genes$chr) | !nzchar(genes$chr) | is.na(genes$gene_start) | is.na(genes$gene_end)
if (any(bad)) {
  message("Skipping ", sum(bad), " genes with missing coordinates")
  genes <- genes[!bad, , drop = FALSE]
}

out_dir <- opts$`output-dir`
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
existing <- existing_gene_folders(out_dir)
new_genes <- genes[!genes$gene_id %in% existing$ids, , drop = FALSE]
n_have <- sum(genes$gene_id %in% existing$ids)
message(n_have, " gene", if (n_have == 1L) "" else "s", " already present; adding ", nrow(new_genes))

read_bound <- function(paths, label) {
  if (length(paths) == 0) return(NULL)
  tabs <- lapply(paths, function(p) {
    message(label, ": ", p)
    d <- read_table_auto(p)
    if (!is.null(d) && "gene_id" %in% names(d)) d$gene_id <- as.character(d$gene_id)
    d
  })
  bind_rows(tabs)
}
region_paths <- if (!is.null(opts$regions)) opts$regions else {
  unlist(lapply(gene_files, function(p) sibling_table(p, "outlier_gene_summary_", "outlier_gene_regions_")), use.names = FALSE)
}
site_paths <- if (!is.null(opts$sites)) opts$sites else {
  unlist(lapply(gene_files, function(p) sibling_table(p, "outlier_gene_summary_", "outlier_gene_sites_")), use.names = FALSE)
}
if (is.null(region_paths)) region_paths <- character(0)
if (is.null(site_paths)) site_paths <- character(0)
region_paths <- unique(region_paths[!is.na(region_paths) & nzchar(region_paths)])
site_paths <- unique(site_paths[!is.na(site_paths) & nzchar(site_paths)])
regions <- read_bound(region_paths, "Regions")
sites <- read_bound(site_paths, "Sites")

new_folders <- character(0)
used <- unname(existing$folders)
if (length(used) == 0) used <- character(0)
for (i in seq_len(nrow(new_genes))) {
  gid <- new_genes$gene_id[[i]]
  folder_name <- safe_gene_dir(gid, used)
  used <- c(used, folder_name)
  dest <- file.path(out_dir, folder_name)
  new_folders[gid] <- dest
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  writeLines(gid, file.path(dest, "gene_id.txt"))
  row <- new_genes[i, setdiff(names(new_genes), ".source"), drop = FALSE]
  write_csv(row, file.path(dest, "gene_summary_row.csv"), na = "")
  if (!is.null(regions) && "gene_id" %in% names(regions)) {
    write_csv(regions[regions$gene_id == gid, , drop = FALSE], file.path(dest, "windows.csv"), na = "")
  }
  if (!is.null(sites) && "gene_id" %in% names(sites)) {
    write_csv(sites[sites$gene_id == gid, , drop = FALSE], file.path(dest, "sites.csv"), na = "")
  }
  if (!is.null(opts$reference)) {
    ok <- write_span_fasta(new_genes$chr[[i]], new_genes$gene_start[[i]], new_genes$gene_end[[i]], gid,
                           opts$reference, opts$samtools, file.path(dest, "gene.fa"))
    if (!ok) message("Could not extract sequence for ", gid)
  }
}

if (nrow(new_genes) > 0 && !is.null(opts$gff) && !is.null(opts$reference)) {
  cds_fa <- file.path(out_dir, "genes_cds_query.fa")
  qmap <- write_cds_fasta(new_genes, opts$gff, opts$reference, opts$samtools, cds_fa)
  if (!is.null(qmap) && nrow(qmap) > 0 && file.exists(cds_fa)) {
    lines <- readLines(cds_fa)
    hdr <- which(startsWith(lines, ">"))
    for (j in seq_len(nrow(qmap))) {
      end <- if (j == nrow(qmap)) length(lines) else hdr[j + 1L] - 1L
      chunk <- lines[hdr[j]:end]
      chunk[1] <- paste0(">", qmap$gene_id[[j]])
      writeLines(chunk, file.path(new_folders[[qmap$gene_id[[j]]]], "cds.fa"))
    }
  }
}

if (nrow(new_genes) > 0 && !is.null(opts$`proteome-config`)) {
  blast_proteome_details(new_genes, opts$gff, opts$reference, opts$`proteome-config`, out_dir, new_folders,
                         opts$samtools, opts$`blastx-cmd`, opts$threads, opts$evalue,
                         opts$`max-target-seqs`, opts$`max-hsps`)
}

cache_dir <- file.path(out_dir, "stat_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
cached <- load_caches(cache_dir)
cached_ids <- if (length(cached) == 0) character(0) else vapply(cached, function(z) as.character(z$meta$gene_id[[1]]), "")
hdf5_dir <- opts$`hdf5-dir`
stamp_path <- file.path(cache_dir, "hdf5_stamp.rds")
old_stamp <- if (!is.null(hdf5_dir) && file.exists(stamp_path)) readRDS(stamp_path) else NULL
new_stamp <- if (!is.null(hdf5_dir)) h5_stamp(hdf5_dir) else NULL
rescore <- !is.null(hdf5_dir) && !identical(old_stamp, new_stamp)

folder_of <- existing$folders
if (length(new_folders) > 0) folder_of[names(new_folders)] <- new_folders
want <- character(0)
for (gid in c(genes$gene_id, existing$ids)) {
  if (gid %in% want) next
  if (rescore || !gid %in% cached_ids) want <- c(want, gid)
}

score_parts <- list()
for (gid in want) {
  if (gid %in% genes$gene_id) {
    score_parts[[gid]] <- genes[genes$gene_id == gid, , drop = FALSE]
    next
  }
  row_file <- file.path(out_dir, folder_of[[gid]], "gene_summary_row.csv")
  if (!file.exists(row_file)) {
    message("No coordinates for ", gid, "; its previous summary is unchanged")
    next
  }
  row <- read_csv(row_file, show_col_types = FALSE)
  if (!all(c("gene_id", "chr", "gene_start", "gene_end") %in% names(row))) {
    message("No coordinates for ", gid, "; its previous summary is unchanged")
    next
  }
  row$gene_id <- as.character(row$gene_id)
  row$chr <- as.character(row$chr)
  row$gene_start <- suppressWarnings(as.integer(row$gene_start))
  row$gene_end <- suppressWarnings(as.integer(row$gene_end))
  if (is.na(row$chr[[1]]) || is.na(row$gene_start[[1]]) || is.na(row$gene_end[[1]])) {
    message("No coordinates for ", gid, "; its previous summary is unchanged")
    next
  }
  score_parts[[gid]] <- row
}
score_genes <- if (length(score_parts) == 0) genes[0, , drop = FALSE] else bind_rows(score_parts)

if (nrow(score_genes) > 0) {
  n_score <- nrow(score_genes)
  if (rescore) message("Scoring ", n_score, " gene", if (n_score == 1L) "" else "s", " against HDF5")
  else message("Scoring ", n_score, " gene", if (n_score == 1L) "" else "s", " not yet in the stat cache")
  bucket <- new.env(parent = emptyenv())
  if (is.null(hdf5_dir)) {
    message("No --hdf5-dir. Stat columns are left empty until a later run passes one.")
  } else {
    files <- stat_h5_files(hdf5_dir)
    if (length(files) == 0) message("No diversity, fst, pbe, qibe, lrpi, or single_position HDF5 files in ", hdf5_dir)
    genes_by_chr <- split(seq_len(nrow(score_genes)), score_genes$chr)
    for (f in files) score_h5_file(f, score_genes, genes_by_chr, bucket)
  }
  variant_rows <- NULL
  vpath <- if (!is.null(hdf5_dir)) variant_h5_path(hdf5_dir) else NULL
  if (!is.null(vpath)) {
    message("Variant calls: ", basename(vpath))
    variant_rows <- variant_counts(score_genes, read_variant_h5(vpath, score_genes))
  }
  for (i in seq_len(nrow(score_genes))) {
    gid <- score_genes$gene_id[[i]]
    sites_sub <- if (!is.null(sites) && "gene_id" %in% names(sites)) sites[sites$gene_id == gid, , drop = FALSE] else NULL
    variant_row <- if (is.null(variant_rows)) NULL else variant_rows[i, , drop = FALSE]
    meta <- meta_from_sources(score_genes[i, , drop = FALSE], variant_row, sites_sub)
    stats_df <- consolidate_stats(bucket[[gid]])
    bn <- folder_of[[gid]]
    if (is.null(bn) || !nzchar(bn)) bn <- safe_gene_dir(gid, used)
    saveRDS(list(meta = meta, stats = stats_df), file.path(cache_dir, paste0(basename(bn), ".rds")))
  }
}

if (!is.null(hdf5_dir) && length(setdiff(want, score_genes$gene_id)) == 0)
  saveRDS(new_stamp, stamp_path)

caches <- load_caches(cache_dir)
if (length(caches) == 0) stop("No genes to summarize")
write_summaries(caches, out_dir)
message("Done.")
