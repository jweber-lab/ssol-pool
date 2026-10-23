#!/usr/bin/env Rscript

###############################################################################
# variant_csq_windows.R
#
# Aggregate annotated variant sites (variant_csq.R output) into per-window counts.
# Windows: the same (chr, start, end) as masked diversity / FST when --template-dir
# (or --template-tsv) points at those tables; otherwise grenedalf's interval grid from
# the reference .fai (--reference with --window-size / --step-size). Templates are used
# even when --reference is also given, so chromosomes grenedalf omitted (no data) and
# all-NA windows dropped by --no-nan-windows stay aligned for collate. --filter-mask-
# total-bed / --filter-mask-total-fasta exclude sites in the same mask as the stats.
#
# Outputs variant_features_w{W}_s{S}.tsv with n_sites, n_sites_<type> (snp/indel),
# n_csq_*, n_region_*, n_effect_*, and genes / genes_protein_altering (semicolon-separated).
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
  library(tibble)
})

read_tab <- function(path, na = c("", "NA", "nan", "NaN", "NAN")) {
  first <- readLines(path, n = 1L)
  delim <- if (grepl("\t", first) && !grepl(",", first)) "\t" else if (grepl(",", first)) "," else "\t"
  if (delim == "\t") readr::read_tsv(path, na = na, show_col_types = FALSE)
  else readr::read_csv(path, na = na, show_col_types = FALSE)
}

parse_window_step_name <- function(b) {
  m <- regmatches(b, regexec("w([0-9]+)_s([0-9]+)", b, perl = TRUE))[[1L]]
  if (length(m) >= 3L) {
    list(window_size = as.numeric(m[2L]), step_size = as.numeric(m[3L]), scale = paste0("w", m[2L], "_s", m[3L]))
  } else {
    list(window_size = NA_real_, step_size = NA_real_, scale = NA_character_)
  }
}

sanitize_suffix <- function(x) {
  s <- gsub("[^A-Za-z0-9]+", "_", x)
  s <- gsub("^_+|_+$", "", s)
  if (!nzchar(s)) s <- "na"
  tolower(s)
}

prepare_sites <- function(sites) {
  d <- sites %>% rename_all(tolower)
  if (!all(c("chr", "pos") %in% names(d))) stop("sites TSV must have chr and pos")
  d$chr <- as.character(d$chr)
  d$pos <- as.numeric(d$pos)
  if (!"variant_type" %in% names(d) && all(c("ref", "alt") %in% names(d))) {
    ref <- as.character(d$ref)
    alt <- sub(",.*", "", as.character(d$alt))
    d$variant_type <- ifelse(nchar(ref) == 1L & nchar(alt) == 1L, "snp", ifelse(nchar(ref) != nchar(alt), "indel", "other"))
  }
  for (gc in c("genes", "genes_protein_altering", "csq_consequences", "region_category", "coding_effect")) {
    if (!gc %in% names(d)) d[[gc]] <- NA_character_
    d[[gc]] <- as.character(d[[gc]])
  }
  d
}

#' grenedalf interval grid: start = 1 + k * step; the last window is the first whose end reaches
#' the chromosome length (end not clipped)
grid_windows <- function(fai, window_size, step_size) {
  n_win <- ifelse(fai$length > window_size, ceiling((fai$length - window_size) / step_size) + 1, 1)
  k <- sequence(n_win) - 1
  start <- 1 + k * step_size
  tibble(chr = rep(fai$chr, n_win), start = start, end = start + window_size - 1)
}

extract_windows_template <- function(path) {
  x <- read_tab(path) %>% rename_all(tolower)
  chr_col <- intersect(c("chr", "chromosome", "chrom"), names(x))[1L]
  if (is.na(chr_col)) stop("template missing chr: ", path)
  if (!all(c("start", "end") %in% names(x))) stop("template missing start/end: ", path)
  x %>%
    rename(chr = !!chr_col) %>%
    mutate(chr = as.character(.data$chr), start = as.numeric(.data$start), end = as.numeric(.data$end)) %>%
    distinct(.data$chr, .data$start, .data$end) %>%
    filter(!is.na(.data$chr), !is.na(.data$start), !is.na(.data$end))
}

#' Union of (chr, start, end) from stats TSVs (diversity / FST / PBE).
extract_windows_templates <- function(files) {
  distinct(bind_rows(lapply(files, extract_windows_template)), .data$chr, .data$start, .data$end)
}

#' Prefer the windows written by masked diversity/FST; fall back to the .fai grid.
windows_for_scale <- function(scale, z, fai) {
  tpl <- if (length(z$files) > 0L) extract_windows_templates(z$files) else NULL
  grid <- if (!is.null(fai)) grid_windows(fai, z$w, z$s) else NULL
  if (!is.null(tpl) && !is.null(grid)) {
    off <- anti_join(tpl, grid, by = c("chr", "start", "end"))
    if (nrow(off) > 0L) {
      i <- 1L
      message("WARNING ", scale, ": ", nrow(off), " template window(s) are off the reference grid (e.g. ",
              off$chr[i], ":", off$start[i], "-", off$end[i], "). Re-run stats with the same --reference-genome.")
    }
    extra <- nrow(anti_join(grid, tpl, by = c("chr", "start", "end")))
    message(scale, ": using ", nrow(tpl), " windows from stats templates",
            if (extra > 0L) paste0(" (omitted ", extra, " reference-grid window(s) not in those tables; typically chromosomes with no data, or all-NA windows dropped by --no-nan-windows)") else "")
    return(tpl)
  }
  if (!is.null(tpl)) {
    message(scale, ": using ", nrow(tpl), " windows from stats templates (no --reference to check the grid)")
    return(tpl)
  }
  if (!is.null(grid)) return(grid)
  stop("No windows for ", scale, ": give --reference or a template TSV with start/end")
}

read_mask_bed <- function(path) {
  raw <- readLines(path, warn = FALSE)
  raw <- raw[!grepl("^(#|track|browser)", raw, ignore.case = TRUE)]
  raw <- raw[nzchar(trimws(raw))]
  if (length(raw) == 0L) return(tibble(chr = character(), start0 = numeric(), end = numeric()))
  parts <- strsplit(raw, "\t", fixed = TRUE)
  n3 <- vapply(parts, function(p) length(p) >= 3L, logical(1L))
  parts <- parts[n3]
  tibble(
    chr = vapply(parts, `[[`, character(1L), 1L),
    start0 = as.numeric(vapply(parts, `[[`, character(1L), 2L)),
    end = as.numeric(vapply(parts, `[[`, character(1L), 3L))
  ) %>%
    filter(!is.na(.data$start0), !is.na(.data$end), .data$end > .data$start0)
}

#' vcftools-style mask FASTA: N/n/X/x/1 are masked (0-based half-open intervals).
read_mask_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0L) return(tibble(chr = character(), start0 = numeric(), end = numeric()))
  hdr <- grepl("^>", lines)
  if (!any(hdr)) stop("Mask FASTA has no '>' headers: ", path)
  idx <- which(hdr)
  chr <- sub("^>(\\S+).*", "\\1", lines[idx])
  end_line <- c(idx[-1L] - 1L, length(lines))
  out_chr <- character(); out_s <- numeric(); out_e <- numeric()
  for (i in seq_along(idx)) {
    from <- idx[i] + 1L
    to <- end_line[i]
    if (from > to) next
    seq <- paste(lines[from:to], collapse = "")
    hits <- gregexpr("[NnXx1]+", seq, perl = TRUE)[[1L]]
    if (hits[1L] == -1L) next
    len <- attr(hits, "match.length")
    out_chr <- c(out_chr, rep(chr[i], length(hits)))
    out_s <- c(out_s, hits - 1L)
    out_e <- c(out_e, hits - 1L + len)
  }
  tibble(chr = out_chr, start0 = as.numeric(out_s), end = as.numeric(out_e))
}

merge_intervals <- function(start0, end) {
  if (length(start0) == 0L) return(list(start0 = start0, end = end))
  o <- order(start0, end)
  start0 <- start0[o]
  end <- end[o]
  cs <- start0[1L]
  ce <- end[1L]
  out_s <- numeric()
  out_e <- numeric()
  for (i in seq_along(start0)) {
    if (start0[i] <= ce) {
      ce <- max(ce, end[i])
    } else {
      out_s <- c(out_s, cs)
      out_e <- c(out_e, ce)
      cs <- start0[i]
      ce <- end[i]
    }
  }
  list(start0 = c(out_s, cs), end = c(out_e, ce))
}

#' 1-based pos is masked when it falls in a 0-based half-open BED interval.
sites_in_mask <- function(sites, bed) {
  masked <- logical(nrow(sites))
  if (nrow(sites) == 0L || nrow(bed) == 0L) return(masked)
  bi <- split(seq_len(nrow(bed)), bed$chr)
  si <- split(seq_len(nrow(sites)), sites$chr)
  for (chr in intersect(names(si), names(bi))) {
    iv <- merge_intervals(bed$start0[bi[[chr]]], bed$end[bi[[chr]]])
    pos <- sites$pos[si[[chr]]]
    i <- findInterval(pos - 1, iv$start0)
    hit <- logical(length(pos))
    ok <- !is.na(i) & i > 0L
    hit[ok] <- (pos[ok] - 1) < iv$end[i[ok]]
    masked[si[[chr]]] <- hit
  }
  masked
}

#' (site, window) index pairs for sites inside windows. Windows of one chromosome must have
#' non-decreasing starts and ends (true for any fixed-width grid), so each site's windows form
#' one contiguous run found with findInterval.
site_window_hits <- function(sites, windows) {
  site_i <- integer(0); win_i <- integer(0)
  wi_by_chr <- split(seq_len(nrow(windows)), windows$chr)
  si_by_chr <- split(seq_len(nrow(sites)), sites$chr)
  for (chr in intersect(names(si_by_chr), names(wi_by_chr))) {
    wi <- wi_by_chr[[chr]]
    wi <- wi[order(windows$start[wi], windows$end[wi])]
    st <- windows$start[wi]; en <- windows$end[wi]
    if (is.unsorted(en)) stop("Windows on ", chr, " have decreasing ends; window widths must be constant")
    si <- si_by_chr[[chr]]
    si <- si[!is.na(sites$pos[si])]
    pos <- sites$pos[si]
    lo <- findInterval(pos, en, left.open = TRUE) + 1L   # first window with end >= pos
    hi <- findInterval(pos, st)                          # last window with start <= pos
    n <- pmax(hi - lo + 1L, 0L)
    site_i <- c(site_i, rep(si, n))
    win_i <- c(win_i, wi[rep(lo, n) + sequence(n) - 1L])
  }
  list(site = site_i, win = win_i)
}

#' Per window: sorted union of semicolon-separated gene lists over its sites
genes_by_window <- function(genes, hits, n_win) {
  out <- rep(NA_character_, n_win)
  keep <- !is.na(genes[hits$site]) & nzchar(genes[hits$site])
  if (!any(keep)) return(out)
  g <- strsplit(genes[hits$site[keep]], ";", fixed = TRUE)
  win <- rep(hits$win[keep], lengths(g))
  g <- trimws(unlist(g, use.names = FALSE))
  ok <- nzchar(g)
  pairs <- unique(data.frame(win = win[ok], g = g[ok], stringsAsFactors = FALSE))
  pairs <- pairs[order(pairs$win, pairs$g), ]
  s <- split(pairs$g, pairs$win)
  out[as.integer(names(s))] <- vapply(s, paste, character(1L), collapse = ";")
  out
}

aggregate_hits <- function(sites, windows) {
  if (nrow(windows) == 0L) return(tibble())
  n_win <- nrow(windows)
  hits <- site_window_hits(sites, windows)
  count_where <- function(flag) tabulate(hits$win[flag[hits$site]], nbins = n_win)

  out <- windows
  out$n_sites <- tabulate(hits$win, nbins = n_win)

  type_lev <- sort(unique(na.omit(sites$variant_type)))
  for (vt in type_lev) out[[paste0("n_sites_", sanitize_suffix(vt))]] <- as.numeric(count_where(sites$variant_type %in% vt))

  # csq tokens: one count per site carrying the token
  tok <- strsplit(ifelse(is.na(sites$csq_consequences), "", sites$csq_consequences), ";", fixed = TRUE)
  tok_site <- rep(seq_len(nrow(sites)), lengths(tok))
  tok <- trimws(unlist(tok, use.names = FALSE))
  keep <- nzchar(tok)
  tok_site <- tok_site[keep]; tok <- tok[keep]
  for (tk in sort(unique(tok))) {
    has <- logical(nrow(sites)); has[tok_site[tok == tk]] <- TRUE
    out[[paste0("n_csq_", sanitize_suffix(tk))]] <- as.numeric(count_where(has))
  }
  for (rc in sort(unique(na.omit(sites$region_category)))) {
    out[[paste0("n_region_", sanitize_suffix(rc))]] <- as.numeric(count_where(sites$region_category %in% rc))
  }
  for (ec in sort(unique(na.omit(sites$coding_effect)))) {
    out[[paste0("n_effect_", sanitize_suffix(ec))]] <- as.numeric(count_where(sites$coding_effect %in% ec))
  }
  out$genes <- genes_by_window(sites$genes, hits, n_win)
  out$genes_protein_altering <- genes_by_window(sites$genes_protein_altering, hits, n_win)
  out %>% arrange(.data$chr, .data$start, .data$end)
}

read_fai <- function(path) {
  fai <- if (grepl("\\.fai$", path)) path else paste0(path, ".fai")
  if (!file.exists(fai)) stop("Reference index not found: ", fai, " (run samtools faidx)")
  x <- utils::read.delim(fai, header = FALSE, colClasses = c("character", "numeric", rep("NULL", 3)))
  tibble(chr = x[[1L]], length = x[[2L]])
}

split_num <- function(x) if (is.null(x)) numeric(0) else as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))

main <- function() {
  option_list <- list(
    make_option(c("--sites-tsv"), type = "character", help = "Annotated variants TSV(s) from variant_csq.R; comma-separate to combine (e.g. SNP and indel outputs)"),
    make_option(c("--output-dir"), type = "character", default = ".", help = "Directory for variant_features_wW_sS.tsv"),
    make_option(c("--reference"), type = "character", default = NULL, help = "Reference FASTA (with .fai) or .fai: windows are built on the grenedalf interval grid from chromosome lengths"),
    make_option(c("--window-size"), type = "character", default = NULL, help = "Window size(s), comma-separated (with --reference)"),
    make_option(c("--step-size"), type = "character", default = NULL, help = "Step size(s), comma-separated and paired with --window-size [default: window/2]"),
    make_option(c("--template-dir"), type = "character", default = NULL, help = "Directory with windowed diversity / FST / PBE tables (scales from file names; their (chr,start,end) are used even with --reference, so windows match masked stats)"),
    make_option(c("--template-tsv"), type = "character", default = NULL, help = "Single template TSV (if only one scale)"),
    make_option(c("--filter-mask-total-bed"), type = "character", default = NULL, help = "Same BED mask as calculate_pi_theta.sh / calculate_fst.sh; sites inside it are omitted from window counts"),
    make_option(c("--filter-mask-total-fasta"), type = "character", default = NULL, help = "vcftools-style mask FASTA (N/n/X/x/1 masked); use instead of --filter-mask-total-bed"),
    make_option(c("--no-nan-windows"), action = "store_true", default = TRUE, help = "When building from --reference, omit windows with n_sites=0 [default]"),
    make_option(c("--keep-nan-windows"), action = "store_true", default = FALSE, help = "Write empty reference-grid windows")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  sites_paths <- if (is.null(opt$`sites-tsv`)) character(0) else trimws(strsplit(opt$`sites-tsv`, ",", fixed = TRUE)[[1L]])
  if (length(sites_paths) == 0L || !all(file.exists(sites_paths))) {
    stop("--sites-tsv must exist: ", paste(sites_paths[!file.exists(sites_paths)], collapse = ", "))
  }

  # Scales: explicit sizes, else template file names
  ws <- split_num(opt$`window-size`)
  ss <- split_num(opt$`step-size`)
  if (length(ss) == 0L) ss <- ws / 2
  if (length(ws) != length(ss)) stop("--window-size and --step-size need the same number of values")
  if (any(ws <= 0 | ss <= 0 | ws %% 1 != 0 | ss %% 1 != 0)) stop("window and step sizes must be positive integers")
  if (!is.null(opt$`filter-mask-total-bed`) && !is.null(opt$`filter-mask-total-fasta`)) {
    stop("Use only one of --filter-mask-total-bed and --filter-mask-total-fasta")
  }
  mask_bed <- NULL
  if (!is.null(opt$`filter-mask-total-bed`)) {
    if (!file.exists(opt$`filter-mask-total-bed`)) stop("Mask BED not found: ", opt$`filter-mask-total-bed`)
    mask_bed <- read_mask_bed(opt$`filter-mask-total-bed`)
  } else if (!is.null(opt$`filter-mask-total-fasta`)) {
    if (!file.exists(opt$`filter-mask-total-fasta`)) stop("Mask FASTA not found: ", opt$`filter-mask-total-fasta`)
    mask_bed <- read_mask_fasta(opt$`filter-mask-total-fasta`)
  }

  template_files <- character(0)
  if (!is.null(opt$`template-tsv`)) template_files <- opt$`template-tsv`
  if (!is.null(opt$`template-dir`)) {
    if (!dir.exists(opt$`template-dir`)) stop("Template directory not found: ", opt$`template-dir`)
    template_files <- c(template_files, list.files(
      opt$`template-dir`,
      pattern = ".*(fst|diversity|pi|theta|pbe).*w[0-9]+_s[0-9]+.*\\.(tsv|csv)$",
      full.names = TRUE, recursive = TRUE, ignore.case = TRUE
    ))
    template_files <- template_files[!grepl("(^|/)\\.", basename(template_files))]
  }
  template_files <- unique(template_files)
  scales <- list()
  for (i in seq_along(ws)) scales[[sprintf("w%d_s%d", as.integer(ws[i]), as.integer(ss[i]))]] <- list(w = ws[i], s = ss[i], files = character(0))
  unnamed <- character(0)
  for (f in template_files) {
    inf <- parse_window_step_name(basename(f))
    if (is.na(inf$scale)) {
      unnamed <- c(unnamed, f)
      next
    }
    if (is.null(scales[[inf$scale]])) scales[[inf$scale]] <- list(w = inf$window_size, s = inf$step_size, files = character(0))
    scales[[inf$scale]]$files <- unique(c(scales[[inf$scale]]$files, f))
  }
  if (length(unnamed) > 0L) {
    if (length(scales) == 1L) {
      scales[[1L]]$files <- unique(c(scales[[1L]]$files, unnamed))
    } else {
      message("Ignoring template file(s) with no w*_s* in the name (pass a single --window-size/--step-size to use them): ",
              paste(basename(unnamed), collapse = ", "))
    }
  }
  if (length(scales) == 0L) stop("Give --reference with --window-size [--step-size], or --template-dir / --template-tsv with w*_s* file names")
  fai <- if (!is.null(opt$reference)) read_fai(opt$reference) else NULL
  if (is.null(fai)) {
    no_tpl <- names(scales)[vapply(scales, function(z) length(z$files) == 0L, logical(1L))]
    if (length(no_tpl) > 0L) stop("--reference is required for scales without template files: ", paste(no_tpl, collapse = ", "))
  }

  sites <- prepare_sites(bind_rows(lapply(sites_paths, function(p) read_tab(p) %>% mutate(across(everything(), as.character)))) %>% readr::type_convert(col_types = readr::cols()))
  if (!is.null(fai)) {
    unknown <- setdiff(unique(sites$chr), fai$chr)
    if (length(unknown) > 0L) message("Sites on ", length(unknown), " sequence(s) not in the reference index are skipped (e.g. ", paste(head(unknown, 3), collapse = ", "), ")")
    sites <- sites %>% filter(.data$chr %in% fai$chr)
  }
  if (!is.null(mask_bed) && nrow(mask_bed) > 0L) {
    n_before <- nrow(sites)
    sites <- sites[!sites_in_mask(sites, mask_bed), , drop = FALSE]
    n_drop <- n_before - nrow(sites)
    if (n_drop > 0L) message("Mask: omitted ", n_drop, " site(s) overlapping --filter-mask-total-bed/fasta (", nrow(sites), " remaining)")
  }

  no_nan_windows <- isTRUE(opt$`no-nan-windows`) && !isTRUE(opt$`keep-nan-windows`)

  dir.create(opt$`output-dir`, showWarnings = FALSE, recursive = TRUE)
  for (sc in names(scales)) {
    z <- scales[[sc]]
    used_template <- length(z$files) > 0L
    win <- windows_for_scale(sc, z, fai)
    agg <- aggregate_hits(sites, win)
    # Match grenedalf --no-nan-windows on the reference grid. Template windows
    # already follow the stats tables (including their --no-nan-windows gaps).
    if (no_nan_windows && !used_template && nrow(agg) > 0L && "n_sites" %in% names(agg)) {
      keep <- !is.na(agg$n_sites) & agg$n_sites > 0
      n_drop <- sum(!keep)
      if (n_drop > 0L) {
        message(sc, ": --no-nan-windows dropped ", n_drop, " empty window(s) (n_sites=0)")
      }
      agg <- agg[keep, , drop = FALSE]
    }
    out_name <- file.path(opt$`output-dir`, sprintf("variant_features_%s.tsv", sc))
    readr::write_tsv(agg, out_name, na = "")
    message("Wrote ", out_name, " (rows=", nrow(agg), ", scale=", sc, ")")
  }
  invisible(NULL)
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
