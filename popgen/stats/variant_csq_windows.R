#!/usr/bin/env Rscript

###############################################################################
# variant_csq_windows.R
#
# Aggregate annotated variant sites (variant_csq.R output) into per-window
# counts using template grenedalf windows from FST or diversity TSV.
#
# Outputs variant_features_w{W}_s{S}.tsv with n_sites, n_csq_*, n_region_*, n_effect_*.
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
  if ("csq_consequences" %in% names(d)) {
    d$toks <- strsplit(ifelse(is.na(d$csq_consequences), "", d$csq_consequences), ";", fixed = TRUE)
    d$toks <- lapply(d$toks, function(z) unique(trimws(z[nzchar(trimws(z))])))
  } else {
    d$toks <- replicate(nrow(d), character(0L), simplify = FALSE)
  }
  if (!"region_category" %in% names(d)) d$region_category <- NA_character_
  if (!"coding_effect" %in% names(d)) d$coding_effect <- NA_character_
  d$region_category <- as.character(d$region_category)
  d$coding_effect <- as.character(d$coding_effect)
  d
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

aggregate_hits <- function(sites, windows) {
  if (nrow(windows) == 0L) {
    return(tibble())
  }
  all_toks <- unique(unlist(sites$toks))
  all_toks <- all_toks[nzchar(all_toks)]
  all_toks <- sort(all_toks)

  rc_lev <- sort(unique(sites$region_category))
  rc_lev <- rc_lev[!is.na(rc_lev)]
  ec_lev <- sort(unique(sites$coding_effect))
  ec_lev <- ec_lev[!is.na(ec_lev)]

  # Join sites to all windows on chr, then position overlap
  hits <- sites %>%
    inner_join(windows, by = "chr", suffix = c("", ".w")) %>%
    filter(.data$pos >= .data$start, .data$pos <= .data$end)

  by_win <- if (nrow(hits) > 0L) hits %>% group_by(.data$chr, .data$start, .data$end) else NULL

  first <- if (nrow(hits) > 0L) {
    by_win %>% summarise(n_sites = dplyr::n(), .groups = "drop")
  } else {
    tibble(chr = character(0), start = numeric(0), end = numeric(0), n_sites = integer(0))
  }

  if (nrow(hits) > 0L) {
    for (tk in all_toks) {
      cn <- paste0("n_csq_", sanitize_suffix(tk))
      cnt <- by_win %>%
        summarise(
          !!cn := sum(vapply(.data$toks, function(z) tk %in% z, logical(1L))),
          .groups = "drop"
        )
      first <- first %>% left_join(cnt, by = c("chr", "start", "end"))
    }

    for (rc in rc_lev) {
      cn <- paste0("n_region_", sanitize_suffix(rc))
      cnt <- by_win %>%
        summarise(!!cn := sum(.data$region_category == rc, na.rm = TRUE), .groups = "drop")
      first <- first %>% left_join(cnt, by = c("chr", "start", "end"))
    }

    for (ec in ec_lev) {
      cn <- paste0("n_effect_", sanitize_suffix(ec))
      cnt <- by_win %>%
        summarise(!!cn := sum(.data$coding_effect == ec, na.rm = TRUE), .groups = "drop")
      first <- first %>% left_join(cnt, by = c("chr", "start", "end"))
    }
  }

  # Ensure every template window appears; fill NAs with 0
  out <- windows %>%
    left_join(first, by = c("chr", "start", "end")) %>%
    mutate(n_sites = if_else(is.na(.data$n_sites), 0L, as.integer(.data$n_sites)))

  num_cols <- setdiff(names(out), c("chr", "start", "end"))
  for (nm in num_cols) {
    v <- out[[nm]]
    v[is.na(v)] <- 0
    out[[nm]] <- if (nm == "n_sites") as.integer(v) else as.numeric(v)
  }

  # Fill missing counter columns with 0 when no hits at all
  for (tk in all_toks) {
    cn <- paste0("n_csq_", sanitize_suffix(tk))
    if (!cn %in% names(out)) out[[cn]] <- 0L
  }
  for (rc in rc_lev) {
    cn <- paste0("n_region_", sanitize_suffix(rc))
    if (!cn %in% names(out)) out[[cn]] <- 0L
  }
  for (ec in ec_lev) {
    cn <- paste0("n_effect_", sanitize_suffix(ec))
    if (!cn %in% names(out)) out[[cn]] <- 0L
  }

  out %>% arrange(.data$chr, .data$start, .data$end)
}

main <- function() {
  option_list <- list(
    make_option(c("--sites-tsv"), type = "character", help = "Annotated variants TSV (variant_csq.R output)"),
    make_option(c("--output-dir"), type = "character", default = ".", help = "Directory for variant_features_wW_sS.tsv"),
    make_option(c("--template-dir"), type = "character", default = NULL, help = "Directory with *fst*w*_s*.tsv or similar window TSVs"),
    make_option(c("--template-tsv"), type = "character", default = NULL, help = "Single template TSV (if only one scale)")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$`sites-tsv`) || !file.exists(opt$`sites-tsv`)) {
    stop("--sites-tsv must exist")
  }
  if (is.null(opt$`template-dir`) && is.null(opt$`template-tsv`)) {
    stop("Provide --template-dir or --template-tsv")
  }

  sites <- prepare_sites(read_tab(opt$`sites-tsv`))

  template_files <- if (!is.null(opt$`template-tsv`)) {
    c(opt$`template-tsv`)
  } else {
    list.files(
      opt$`template-dir`,
      pattern = ".*(fst|diversity|pi|theta).*w[0-9]+_s[0-9]+.*\\.(tsv|csv)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }
  template_files <- setdiff(template_files, NA_character_)
  if (length(template_files) == 0L) {
    stop("No template files found; check --template-dir pattern or pass --template-tsv")
  }

  by_path <- lapply(template_files, function(f) {
    inf <- parse_window_step_name(basename(f))
    list(path = f, scale = inf$scale, window_size = inf$window_size, step_size = inf$step_size)
  })
  scales <- unique(vapply(by_path, function(z) z$scale, character(1L)))
  scales <- scales[!is.na(scales)]
  if (length(scales) == 0L) {
    stop("Could not infer w*_s* scale from template file names. Rename files to include w1000_s500 style suffix.")
  }

  dir.create(opt$`output-dir`, showWarnings = FALSE, recursive = TRUE)

  for (sc in scales) {
    tpls <- Filter(function(z) identical(z$scale, sc), by_path)
    tpl_paths <- vapply(tpls, function(z) z$path, character(1L))
    # Union windows from all templates of this scale
    win <- NULL
    for (tp in tpl_paths) {
      wi <- extract_windows_template(tp)
      win <- if (is.null(win)) wi else bind_rows(win, wi)
    }
    win <- win %>% distinct(.data$chr, .data$start, .data$end)

    agg <- aggregate_hits(sites, win)

    ## Parse W and S from scale
    m <- regmatches(sc, regexec("^w([0-9]+)_s([0-9]+)$", sc))[[1L]]
    out_name <- file.path(opt$`output-dir`, sprintf("variant_features_%s.tsv", sc))
    if (length(m) >= 3L) {
      out_name <- file.path(opt$`output-dir`, sprintf("variant_features_w%s_s%s.tsv", m[2L], m[3L]))
    }
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
