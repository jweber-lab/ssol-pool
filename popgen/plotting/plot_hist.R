#!/usr/bin/env Rscript

###############################################################################
# plot_hist.R
#
# Genome-wide (or subset) distributions for a chosen statistic.
# By default only a histogram is drawn; optional panels: kernel density, ECDF, diff.
# Histogram bin breaks and --q-lines quantiles use the pooled (all-group) distribution
# so facets/overlays are comparable. --compare A,B filters to two groups and adds a
# diff panel: proportion in bin for A minus proportion for B on those same breaks.
# --hist-mode controls histogram bar heights (count vs density), not the optional density panel.
#
# Inputs mirror plot_region.R: reads from TSV dirs and/or HDF5 (collate output).
#
# Requires: ggplot2, dplyr, tidyr, readr, optparse, patchwork
# Optional: hdf5r when using --hdf5-dir
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)
})
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required. Install with: install.packages(\"patchwork\")")
}
library(patchwork)

initial_args <- commandArgs(trailingOnly = FALSE)
file_arg <- initial_args[substr(initial_args, 1, 7) == "--file="]
script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg)) else "."
plot_common_path <- file.path(script_dir, "plot_common.R")
if (!file.exists(plot_common_path)) plot_common_path <- file.path(getwd(), "plot_common.R")
if (!file.exists(plot_common_path)) stop("plot_common.R not found in script dir or getwd()")
source(plot_common_path)

option_list <- list(
  make_option("--chromosome", type = "character", default = NULL, help = "Subset to chromosome (optional)"),
  make_option("--region", type = "character", default = NULL, help = "Subset to region CHR:START-END (optional)"),
  make_option("--diversity-dir", type = "character", default = NULL, help = "Directory with diversity TSVs"),
  make_option("--fst-dir", type = "character", default = NULL, help = "Directory with FST TSVs"),
  make_option("--pbe-dir", type = "character", default = NULL, help = "Directory with PBE TSVs"),
  make_option("--seq-qual-dir", type = "character", default = NULL, help = "Directory with seq_qual TSVs"),
  make_option("--hdf5-dir", type = "character", default = NULL, help = "Collate output dir with .h5 files"),
  make_option("--window-size", type = "numeric", default = NULL, help = "Window size to use (default: first available)"),
  make_option("--step-size", type = "numeric", default = NULL, help = "Step size (optional)"),
  make_option("--output-dir", type = "character", default = ".", help = "Output directory"),
  make_option("--file-prefix", type = "character", default = "", help = "Prefix for output filename"),
  make_option("--stat", type = "character", default = "pi",
    help = "Statistic: coverage, mapq, n_snps, pi, theta, tajima_d, fst, pbe [default: %default]"),
  make_option("--y-value", type = "character", default = "value",
    help = "For diversity/FST/PBE: value, rank, or quantile [default: %default]"),
  make_option("--transform", type = "character", default = "none",
    help = "Transform on the statistic axis: none, log, or asinh [default: %default]"),
  make_option("--bins", type = "integer", default = 60L, help = "Histogram bins [default: %default]"),
  make_option("--panels", type = "character", default = "histogram",
    help = paste0("Comma-separated figure panels: histogram, density, ecdf, diff ",
      "(diff requires --compare; default: histogram only; density = kernel curve, not --hist-mode).")),
  make_option("--hist-mode", type = "character", default = "count",
    help = "Histogram bar heights: count or after_stat(density) for bars only [default: %default]"),
  make_option("--overlay", action = "store_true", default = FALSE,
    help = "Overlay multiple groups in a single panel (color/fill = group)"),
  make_option("--max-groups", type = "integer", default = 30L,
    help = "Maximum number of groups to plot (to avoid unreadable legends) [default: %default]"),
  make_option("--q-lines", type = "character", default = "0.5,0.95,0.99",
    help = "Comma-separated quantiles of the pooled sample to mark as vertical lines [default: %default]"),
  make_option("--compare", type = "character", default = NULL,
    help = "Two group names 'A,B' (first minus second in diff panel); limits data to those groups and adds diff panel"),
  make_option("--width", type = "numeric", default = 12, help = "Figure width (in)"),
  make_option("--height", type = "numeric", default = 8, help = "Figure height (in)"),
  make_option("--dpi", type = "numeric", default = NULL, help = "DPI for PNG (default: 300)"),
  make_option("--plot-format", type = "character", default = "png", help = "png, pdf, svg, both, or all [default: %default]"),
  make_option("--verbose", action = "store_true", default = FALSE, help = "Verbose")
)
opts <- parse_args(OptionParser(option_list = option_list, usage = "usage: %prog [options]"))

compare_pair <- NULL
if (!is.null(opts$compare) && nzchar(opts$compare)) {
  parts <- trimws(strsplit(opts$compare, ",")[[1]])
  if (length(parts) != 2L || !nzchar(parts[1L]) || !nzchar(parts[2L])) {
    stop("--compare must be exactly two comma-separated group names, e.g. Echo,Cheney")
  }
  compare_pair <- parts
}

valid_stats <- c("coverage", "mapq", "n_snps", "pi", "theta", "tajima_d", "fst", "pbe")
if (!opts$stat %in% valid_stats) stop("--stat must be one of: ", paste(valid_stats, collapse = ", "))
if (!opts$`y-value` %in% c("value", "rank", "quantile")) stop("--y-value must be one of: value, rank, quantile")
if (!opts$transform %in% c("none", "log", "asinh")) stop("--transform must be one of: none, log, asinh")
if (!opts$`hist-mode` %in% c("count", "density")) stop("--hist-mode must be one of: count, density")
panels_req <- trimws(tolower(strsplit(opts$panels, ",")[[1]]))
panels_req <- panels_req[nzchar(panels_req)]
valid_panels <- c("histogram", "density", "ecdf", "diff")
bad_p <- setdiff(panels_req, valid_panels)
if (length(bad_p) > 0) stop("Unknown --panels: ", paste(bad_p, collapse = ", "), ". Valid: ", paste(valid_panels, collapse = ", "))
if (length(panels_req) == 0) stop("--panels must list at least one of: ", paste(valid_panels, collapse = ", "))
panels_req <- panels_req[!duplicated(panels_req)]
if (!is.null(compare_pair) && !"diff" %in% panels_req) panels_req <- c(panels_req, "diff")
panels_req <- panels_req[!duplicated(panels_req)]
if ("diff" %in% panels_req && is.null(compare_pair)) stop("Panel 'diff' requires --compare 'A,B'")

format_parts <- trimws(tolower(strsplit(opts$`plot-format`, ",")[[1]]))
if ("both" %in% format_parts) format_parts <- c(setdiff(format_parts, "both"), "png", "pdf")
if ("all" %in% format_parts) format_parts <- c(setdiff(format_parts, "all"), "png", "pdf", "svg")
format_parts <- unique(format_parts)

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Parse region subset (optional)
# ---------------------------------------------------------------------------
chr_region <- opts$chromosome
start_region <- NA_real_
end_region <- NA_real_
if (!is.null(opts$region) && nzchar(opts$region)) {
  parts <- strsplit(trimws(opts$region), ":")[[1]]
  if (length(parts) >= 1) chr_region <- trimws(parts[1])
  if (length(parts) >= 2) {
    range_parts <- strsplit(trimws(parts[2]), "-")[[1]]
    if (length(range_parts) >= 2) {
      start_region <- as.numeric(trimws(range_parts[1]))
      end_region <- as.numeric(trimws(range_parts[2]))
    }
  }
}

add_pos <- function(d) {
  if ("start" %in% names(d) && "end" %in% names(d))
    d <- d %>% mutate(pos = (.data$start + .data$end) / 2)
  else if (!"pos" %in% names(d) && "position" %in% names(d))
    d <- d %>% rename(pos = "position")
  else if (!"pos" %in% names(d) && "start" %in% names(d))
    d <- d %>% mutate(pos = .data$start)
  d
}

filter_region <- function(d) {
  if (!is.null(chr_region) && nzchar(chr_region)) d <- d %>% filter(.data$chr == chr_region)
  if (!is.na(start_region) && "end" %in% names(d)) d <- d %>% filter(.data$end >= start_region)
  if (!is.na(end_region) && "start" %in% names(d)) d <- d %>% filter(.data$start <= end_region)
  d
}

apply_transform <- function(values, tfm) {
  if (tfm == "log") {
    min_v <- min(values, na.rm = TRUE)
    if (min_v <= 0) log1p(pmax(0, values)) else log(values)
  } else if (tfm == "asinh") {
    med <- median(values, na.rm = TRUE)
    s <- sd(values, na.rm = TRUE)
    if (is.na(s) || s == 0) s <- 1
    asinh((values - med) / s)
  } else {
    values
  }
}

build_x_scale <- function(original_values, tfm) {
  if (tfm == "none") return(NULL)
  orig_range <- range(original_values, na.rm = TRUE)
  if (any(!is.finite(orig_range))) return(NULL)

  fmt_label <- function(x) {
    sapply(x, function(v) {
      if (!is.finite(v) || is.na(v)) return("")
      av <- abs(v)
      if (av == 0) "0"
      else if (av >= 1e4 || av < 1e-3) formatC(v, format = "g", digits = 3)
      else {
        s <- formatC(v, format = "f", digits = 4)
        sub("0+$", "", sub("\\.$", "", s))
      }
    })
  }

  if (tfm == "log") {
    min_v <- min(original_values, na.rm = TRUE)
    use_log1p <- (min_v <= 0)
    if (!use_log1p) {
      # Data x is log(value): space ticks evenly in log-space (not linear pretty on originals),
      # otherwise marks cluster toward high x.
      ov <- original_values[is.finite(original_values) & original_values > 0]
      if (length(ov) < 1L) return(NULL)
      lr <- range(log(ov), na.rm = TRUE)
      if (!all(is.finite(lr))) return(NULL)
      x_breaks <- pretty(lr, n = 7)
      x_breaks <- x_breaks[is.finite(x_breaks)]
      orig_ticks <- exp(x_breaks)
      ok <- is.finite(orig_ticks) & orig_ticks > 0
      x_breaks <- x_breaks[ok]
      orig_ticks <- orig_ticks[ok]
      if (length(x_breaks) < 1L) return(NULL)
      scale_x_continuous(breaks = x_breaks, labels = fmt_label(orig_ticks))
    } else {
      xv <- log1p(pmax(0, original_values[is.finite(original_values)]))
      if (length(xv) < 1L) return(NULL)
      xr <- range(xv, na.rm = TRUE)
      if (!all(is.finite(xr))) return(NULL)
      x_breaks <- pretty(xr, n = 7)
      x_breaks <- x_breaks[is.finite(x_breaks) & x_breaks >= 0]
      if (length(x_breaks) < 1L) return(NULL)
      orig_lab <- pmax(0, expm1(x_breaks))
      scale_x_continuous(breaks = x_breaks, labels = fmt_label(orig_lab))
    }
  } else if (tfm == "asinh") {
    med <- median(original_values, na.rm = TRUE)
    s <- sd(original_values, na.rm = TRUE)
    if (is.na(s) || s == 0) s <- 1
    fwd <- function(x) asinh((x - med) / s)
    orig_ticks <- pretty(orig_range, n = 7)
    x_breaks <- fwd(orig_ticks)
    scale_x_continuous(breaks = x_breaks, labels = fmt_label(orig_ticks))
  } else {
    NULL
  }
}

# ---------------------------------------------------------------------------
# Load from HDF5 (preferred) with TSV fallbacks like plot_region.R
# ---------------------------------------------------------------------------
div_data <- NULL
fst_data <- NULL
pbe_data <- NULL
cov_data <- NULL
mapq_data <- NULL
window_size_use <- opts$`window-size`
step_size_use <- opts$`step-size`

if (!is.null(opts$`hdf5-dir`) && dir.exists(opts$`hdf5-dir`)) {
  hdir <- opts$`hdf5-dir`

  div_h5 <- list.files(hdir, pattern = "^diversity_.*\\.h5$", full.names = TRUE)
  if (length(div_h5) > 0) {
    if (!is.null(window_size_use)) div_h5 <- div_h5[grepl(paste0("_w", window_size_use, "_"), div_h5, fixed = TRUE)]
    if (length(div_h5) > 0 && !is.null(step_size_use)) div_h5 <- div_h5[grepl(paste0("_s", step_size_use, "\\."), div_h5, fixed = TRUE)]
    if (length(div_h5) > 0) {
      div_h5 <- div_h5[1]
      ws <- parse_window_step_from_filename(basename(div_h5))
      if (is.null(window_size_use) || is.na(window_size_use)) window_size_use <- ws$window_size
      d <- read_h5_windows(div_h5, "windows")
      if (!is.null(d) && nrow(d) > 0) {
        d <- add_pos(d) %>% filter_region()
        if (nrow(d) > 0) {
          div_data <- d
          for (stat in c("pi", "theta", "tajima_d")) {
            y_col <- switch(opts$`y-value`, rank = paste0(stat, "_rank"), quantile = paste0(stat, "_quantile"), stat)
            if (y_col %in% names(div_data) && y_col != stat) div_data[[stat]] <- div_data[[y_col]]
          }
          if ("mean_coverage" %in% names(div_data) && is.null(cov_data)) {
            cov_data <- div_data %>%
              select(any_of(c("chr", "start", "end", "sample", "mean_coverage", "pos"))) %>%
              filter(!is.na(mean_coverage))
          }
          if ("mean_mapping_quality" %in% names(div_data) && is.null(mapq_data)) {
            mapq_data <- div_data %>%
              select(any_of(c("chr", "start", "end", "sample", "mean_mapping_quality", "pos"))) %>%
              filter(!is.na(mean_mapping_quality))
          }
        }
      }
    }
  }

  fst_h5 <- list.files(hdir, pattern = "^fst_.*\\.h5$", full.names = TRUE)
  if (length(fst_h5) > 0) {
    if (!is.null(window_size_use) && !is.na(window_size_use)) fst_h5 <- fst_h5[grepl(paste0("_w", window_size_use), fst_h5, fixed = TRUE)]
    if (length(fst_h5) > 0) {
      fst_h5 <- fst_h5[1]
      d <- read_h5_windows(fst_h5, "windows")
      if (is.null(d)) d <- read_h5_windows(fst_h5, "sites")
      if (!is.null(d) && nrow(d) > 0) {
        if (!"sample_pair" %in% names(d) && all(c("pop1", "pop2") %in% names(d))) d <- d %>% mutate(sample_pair = paste(pop1, pop2, sep = ":"))
        d <- add_pos(d) %>% filter_region()
        if (nrow(d) > 0) {
          fst_data <- d
          y_col <- switch(opts$`y-value`, rank = "fst_rank", quantile = "fst_quantile", "fst")
          if (y_col %in% names(fst_data) && y_col != "fst") fst_data[["fst"]] <- fst_data[[y_col]]
        }
      }
    }
  }

  pbe_h5 <- list.files(hdir, pattern = "^pbe_.*\\.h5$", full.names = TRUE)
  if (length(pbe_h5) > 0) {
    if (!is.null(window_size_use) && !is.na(window_size_use)) pbe_h5 <- pbe_h5[grepl(paste0("_w", window_size_use), pbe_h5, fixed = TRUE)]
    if (length(pbe_h5) > 0) {
      pbe_h5 <- pbe_h5[1]
      d <- read_h5_windows(pbe_h5, "windows")
      if (is.null(d)) d <- read_h5_windows(pbe_h5, "sites")
      if (!is.null(d) && nrow(d) > 0) {
        if (!"trio" %in% names(d) && all(c("pop1", "pop2", "pop3") %in% names(d))) d <- d %>% mutate(trio = paste(pop1, pop2, pop3, sep = ":"))
        d <- add_pos(d) %>% filter_region()
        if (nrow(d) > 0) {
          pbe_data <- d
          y_col <- switch(opts$`y-value`, rank = "pbe_rank", quantile = "pbe_quantile", "pbe")
          if (y_col %in% names(pbe_data) && y_col != "pbe") pbe_data[["pbe"]] <- pbe_data[[y_col]]
        }
      }
    }
  }
}

if ((is.null(cov_data) || is.null(mapq_data)) && !is.null(opts$`seq-qual-dir`) && dir.exists(opts$`seq-qual-dir`)) {
  sq_files <- list.files(opts$`seq-qual-dir`, pattern = "seq_qual.*\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (length(sq_files) > 0) {
    sq_list <- lapply(sq_files, function(f) {
      x <- read_tsv(f, show_col_types = FALSE)
      if ("chromosome" %in% names(x) && !"chr" %in% names(x)) x <- x %>% rename(chr = chromosome)
      add_pos(x)
    })
    sq_combined <- bind_rows(sq_list) %>% filter_region()
    if (is.null(cov_data) && "mean_coverage" %in% names(sq_combined)) cov_data <- sq_combined %>% filter(!is.na(mean_coverage))
    if (is.null(mapq_data) && "mean_mapping_quality" %in% names(sq_combined)) mapq_data <- sq_combined %>% filter(!is.na(mean_mapping_quality))
  }
}

# Basic TSV fallbacks (single file) for compatibility with plot_region.R
if (is.null(div_data) && !is.null(opts$`diversity-dir`) && dir.exists(opts$`diversity-dir`)) {
  div_tsv <- list.files(opts$`diversity-dir`, pattern = "diversity.*\\.(tsv|csv)$", full.names = TRUE)
  if (length(div_tsv) > 0) {
    d <- tryCatch(read_tsv(div_tsv[1], show_col_types = FALSE), error = function(e) NULL)
    if (is.null(d)) d <- tryCatch(read_csv(div_tsv[1], show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(d)) {
      if ("chromosome" %in% names(d) && !"chr" %in% names(d)) d <- d %>% rename(chr = chromosome)
      if ("position" %in% names(d) && !"pos" %in% names(d)) d <- d %>% rename(pos = position)
      if (!"sample" %in% names(d)) d <- d %>% mutate(sample = "sample")
      div_data <- add_pos(d) %>% filter_region()
    }
  }
}
if (is.null(fst_data) && !is.null(opts$`fst-dir`) && dir.exists(opts$`fst-dir`)) {
  fst_tsv <- list.files(opts$`fst-dir`, pattern = "fst.*\\.(tsv|csv)$", full.names = TRUE)
  if (length(fst_tsv) > 0) {
    d <- tryCatch(read_tsv(fst_tsv[1], show_col_types = FALSE), error = function(e) NULL)
    if (is.null(d)) d <- tryCatch(read_csv(fst_tsv[1], show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(d)) {
      if ("chromosome" %in% names(d) && !"chr" %in% names(d)) d <- d %>% rename(chr = chromosome)
      if ("position" %in% names(d) && !"pos" %in% names(d)) d <- d %>% rename(pos = position)
      fst_col <- names(d)[grepl("\\.fst$|^fst$", names(d), ignore.case = TRUE)][1]
      if (!is.na(fst_col)) {
        d <- d %>% rename(fst = all_of(fst_col))
        if (!"sample_pair" %in% names(d)) d <- d %>% mutate(sample_pair = sub("\\.fst$", "", fst_col))
        fst_data <- add_pos(d) %>% filter_region()
      }
    }
  }
}
if (is.null(pbe_data) && !is.null(opts$`pbe-dir`) && dir.exists(opts$`pbe-dir`)) {
  pbe_tsv <- list.files(opts$`pbe-dir`, pattern = "pbe.*\\.(tsv|csv)$", full.names = TRUE)
  if (length(pbe_tsv) > 0) {
    d <- tryCatch(read_tsv(pbe_tsv[1], show_col_types = FALSE), error = function(e) NULL)
    if (is.null(d)) d <- tryCatch(read_csv(pbe_tsv[1], show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(d)) {
      if ("chromosome" %in% names(d) && !"chr" %in% names(d)) d <- d %>% rename(chr = chromosome)
      if ("position" %in% names(d) && !"pos" %in% names(d)) d <- d %>% rename(pos = position)
      pbe_col <- names(d)[grepl("^pbe$", names(d), ignore.case = TRUE)][1]
      if (!is.na(pbe_col)) {
        d <- d %>% rename(pbe = all_of(pbe_col))
        if (!"trio" %in% names(d)) d <- d %>% mutate(trio = "trio")
        pbe_data <- add_pos(d) %>% filter_region()
      }
    }
  }
}

if (opts$verbose) {
  message("[plot_hist] subset chr=", ifelse(is.null(chr_region), "ALL", chr_region),
    " start=", ifelse(is.na(start_region), "NA", start_region),
    " end=", ifelse(is.na(end_region), "NA", end_region))
}

# ---------------------------------------------------------------------------
# Build plotting dataframe for requested stat
# ---------------------------------------------------------------------------
plot_df <- NULL
group_col <- NULL
value_col <- NULL
legend_title <- NULL
base_label <- NULL

if (opts$stat == "coverage") {
  if (is.null(cov_data) || nrow(cov_data) == 0) stop("No coverage data found (need seq_qual or diversity HDF5 with mean_coverage).")
  if (!"sample" %in% names(cov_data)) cov_data <- cov_data %>% mutate(sample = "sample")
  plot_df <- cov_data %>% transmute(chr, pos, group = .data$sample, value = .data$mean_coverage) %>% filter(!is.na(value))
  legend_title <- "Sample"
  base_label <- "Mean coverage"
} else if (opts$stat == "mapq") {
  if (is.null(mapq_data) || nrow(mapq_data) == 0) stop("No MAPQ data found (need seq_qual or diversity HDF5 with mean_mapping_quality).")
  if (!"sample" %in% names(mapq_data)) mapq_data <- mapq_data %>% mutate(sample = "sample")
  plot_df <- mapq_data %>% transmute(chr, pos, group = .data$sample, value = .data$mean_mapping_quality) %>% filter(!is.na(value))
  legend_title <- "Sample"
  base_label <- "Mean MAPQ"
} else if (opts$stat == "n_snps") {
  if (is.null(div_data) || nrow(div_data) == 0 || !"n_snps" %in% names(div_data)) stop("No n_snps data found (need diversity input with n_snps).")
  if (!"sample" %in% names(div_data)) div_data <- div_data %>% mutate(sample = "sample")
  plot_df <- div_data %>% transmute(chr, pos, group = .data$sample, value = .data$n_snps) %>% filter(!is.na(value))
  legend_title <- "Sample"
  base_label <- "Number of SNPs"
} else if (opts$stat %in% c("pi", "theta", "tajima_d")) {
  if (is.null(div_data) || nrow(div_data) == 0 || !opts$stat %in% names(div_data)) stop("No diversity data found for stat: ", opts$stat)
  if (!"sample" %in% names(div_data)) div_data <- div_data %>% mutate(sample = "sample")
  plot_df <- div_data %>% transmute(chr, pos, group = .data$sample, value = .data[[opts$stat]]) %>% filter(!is.na(value))
  legend_title <- "Sample"
  base_label <- switch(opts$stat, pi = "\u03c0", theta = "\u03b8", tajima_d = "Tajima's D")
  if (opts$`y-value` == "rank") base_label <- paste(base_label, "rank")
  if (opts$`y-value` == "quantile") base_label <- paste(base_label, "quantile")
} else if (opts$stat == "fst") {
  if (is.null(fst_data) || nrow(fst_data) == 0 || !"fst" %in% names(fst_data)) stop("No FST data found.")
  if (!"sample_pair" %in% names(fst_data) && all(c("pop1", "pop2") %in% names(fst_data))) fst_data <- fst_data %>% mutate(sample_pair = paste(pop1, pop2, sep = ":"))
  plot_df <- fst_data %>% transmute(chr, pos, group = .data$sample_pair, value = .data$fst) %>% filter(!is.na(value))
  legend_title <- "Pair"
  base_label <- if (opts$`y-value` == "rank") "FST rank" else if (opts$`y-value` == "quantile") "FST quantile" else "FST"
} else if (opts$stat == "pbe") {
  if (is.null(pbe_data) || nrow(pbe_data) == 0 || !"pbe" %in% names(pbe_data)) stop("No PBE data found.")
  if (!"trio" %in% names(pbe_data) && all(c("pop1", "pop2", "pop3") %in% names(pbe_data))) pbe_data <- pbe_data %>% mutate(trio = paste(pop1, pop2, pop3, sep = ":"))
  plot_df <- pbe_data %>% transmute(chr, pos, group = .data$trio, value = .data$pbe) %>% filter(!is.na(value))
  legend_title <- "Trio"
  base_label <- if (opts$`y-value` == "rank") "PBE rank" else if (opts$`y-value` == "quantile") "PBE quantile" else "PBE"
}

if (is.null(plot_df) || nrow(plot_df) == 0) stop("No data to plot after subsetting and NA filtering.")
plot_df$group <- as.character(plot_df$group)

groups <- sort(unique(plot_df$group))
if (length(groups) > opts$`max-groups`) {
  stop("Too many groups to plot (", length(groups), " > --max-groups ", opts$`max-groups`, "). Use subsetting inputs or increase --max-groups.")
}

plot_df <- plot_df %>% mutate(x_val = apply_transform(.data$value, opts$transform))

# Pooled x range -> shared histogram bin breaks for every group / facet / overlay
xv_all <- plot_df$x_val[is.finite(plot_df$x_val)]
if (length(xv_all) < 1L) stop("No finite x values to plot (after transform).")
rng_x <- range(xv_all, na.rm = TRUE)
if (!is.finite(rng_x[1]) || !is.finite(rng_x[2])) stop("Non-finite range for transformed x values.")
if (rng_x[1] == rng_x[2]) {
  eps <- if (rng_x[1] == 0) 0.5 else max(abs(rng_x[1]) * 1e-6, .Machine$double.eps)
  rng_x <- rng_x + c(-eps, eps)
}
nb <- as.integer(opts$bins)
if (is.na(nb) || nb < 1L) stop("--bins must be a positive integer")
global_breaks <- seq(rng_x[1], rng_x[2], length.out = nb + 1L)

# Quantiles of the pooled distribution (raw --stat values), drawn at transformed positions
q_probs <- suppressWarnings(as.numeric(trimws(strsplit(opts$`q-lines`, ",")[[1]])))
q_probs <- q_probs[is.finite(q_probs) & q_probs >= 0 & q_probs <= 1]
q_probs <- unique(q_probs)
q_df <- if (length(q_probs) > 0L) {
  tibble(
    q_prob = q_probs,
    q_value = as.numeric(stats::quantile(plot_df$value, probs = q_probs, na.rm = TRUE))
  ) %>%
    mutate(q_x = apply_transform(.data$q_value, opts$transform))
} else {
  tibble(q_prob = numeric(0), q_value = numeric(0), q_x = numeric(0))
}

# Original-scale values for x-axis tick mapping (pooled; unchanged by --compare filter)
orig_vals_pooled_for_x_scale <- plot_df$value

compare_suffix <- ""
if (!is.null(compare_pair)) {
  miss <- setdiff(compare_pair, unique(plot_df$group))
  if (length(miss) > 0L) {
    stop("--compare name(s) not in data: ", paste(miss, collapse = ", "), ". Available: ", paste(sort(unique(plot_df$group)), collapse = ", "))
  }
  compare_suffix <- paste0(
    "_compare_",
    gsub("[^A-Za-z0-9]+", "_", compare_pair[1L]),
    "_vs_",
    gsub("[^A-Za-z0-9]+", "_", compare_pair[2L])
  )
  plot_df <- plot_df %>% filter(.data$group %in% compare_pair)
  plot_df$group <- factor(plot_df$group, levels = compare_pair)
  groups <- compare_pair
}

if (nrow(plot_df) == 0L) stop("No rows left after --compare filter.")

plot_df$group <- as.character(plot_df$group)

pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(groups)), groups)

# Bin proportion difference (compare_pair[1] - compare_pair[2]) on global breaks
diff_df <- NULL
bar_width_diff <- mean(diff(global_breaks), na.rm = TRUE) * 0.92
if (!is.finite(bar_width_diff) || bar_width_diff <= 0) bar_width_diff <- 0.01
if ("diff" %in% panels_req) {
  g1 <- compare_pair[1L]
  g2 <- compare_pair[2L]
  x1 <- plot_df$x_val[plot_df$group == g1]
  x2 <- plot_df$x_val[plot_df$group == g2]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  h1 <- graphics::hist(x1, breaks = global_breaks, plot = FALSE, include.lowest = TRUE)
  h2 <- graphics::hist(x2, breaks = global_breaks, plot = FALSE, include.lowest = TRUE)
  n1 <- sum(h1$counts)
  n2 <- sum(h2$counts)
  p1 <- if (n1 > 0L) h1$counts / n1 else rep(0, length(h1$counts))
  p2 <- if (n2 > 0L) h2$counts / n2 else rep(0, length(h2$counts))
  mids <- (h1$breaks[-length(h1$breaks)] + h1$breaks[-1]) / 2
  diff_df <- tibble(mid_x = mids, diff_prop = as.numeric(p1 - p2))
}

theme_panel <- theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
  theme(panel.grid.minor = element_blank())

xlab <- base_label
if (opts$transform == "log") xlab <- paste0(xlab, " (log-scaled axis labels)")
if (opts$transform == "asinh") xlab <- paste0(xlab, " (asinh-scaled axis labels)")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
hist_aes <- aes(x = .data$x_val)
dens_aes <- aes(x = .data$x_val)
ecdf_aes <- aes(x = .data$x_val)
if (opts$overlay) {
  hist_aes <- aes(x = .data$x_val, fill = .data$group)
  dens_aes <- aes(x = .data$x_val, color = .data$group)
  ecdf_aes <- aes(x = .data$x_val, color = .data$group)
}

p_hist <- ggplot(plot_df, hist_aes) +
  geom_histogram(
    aes(y = if (opts$`hist-mode` == "density") after_stat(density) else after_stat(count)),
    breaks = global_breaks,
    alpha = if (opts$overlay) 0.4 else 0.85,
    position = "identity",
    color = "grey30",
    linewidth = 0.15,
    na.rm = TRUE
  ) +
  labs(x = xlab, y = if (opts$`hist-mode` == "density") "Density" else "Count", title = "Histogram") +
  theme_panel
if (opts$overlay) p_hist <- p_hist + scale_fill_manual(name = legend_title, values = pal, drop = FALSE)
if (!opts$overlay) p_hist <- p_hist + facet_wrap(~ group, scales = "free_y")

p_dens <- ggplot(plot_df, dens_aes) +
  geom_density(linewidth = 0.8, alpha = 0.3, na.rm = TRUE) +
  labs(x = xlab, y = "Density", title = "Density") +
  theme_panel
if (opts$overlay) p_dens <- p_dens + scale_color_manual(name = legend_title, values = pal, drop = FALSE)
if (!opts$overlay) p_dens <- p_dens + facet_wrap(~ group, scales = "free_y")

p_ecdf <- ggplot(plot_df, ecdf_aes) +
  stat_ecdf(linewidth = 0.8, alpha = 0.85, na.rm = TRUE) +
  labs(x = xlab, y = "ECDF", title = "ECDF") +
  theme_panel
if (opts$overlay) p_ecdf <- p_ecdf + scale_color_manual(name = legend_title, values = pal, drop = FALSE)
if (!opts$overlay) p_ecdf <- p_ecdf + facet_wrap(~ group)

# Quantile lines (pooled distribution): same x-intercepts in every facet/overlay
if (nrow(q_df) > 0) {
  if ("histogram" %in% panels_req) {
    p_hist <- p_hist + geom_vline(data = q_df, aes(xintercept = .data$q_x), linewidth = 0.4, alpha = 0.85, color = "grey20")
  }
  if ("density" %in% panels_req) {
    p_dens <- p_dens + geom_vline(data = q_df, aes(xintercept = .data$q_x), linewidth = 0.4, alpha = 0.85, color = "grey20")
  }
}

xs <- build_x_scale(orig_vals_pooled_for_x_scale, opts$transform)
if (!is.null(xs)) {
  if ("histogram" %in% panels_req) p_hist <- p_hist + xs
  if ("density" %in% panels_req) p_dens <- p_dens + xs
  if ("ecdf" %in% panels_req) p_ecdf <- p_ecdf + xs
}

p_diff <- NULL
if ("diff" %in% panels_req) {
  p_diff <- ggplot(diff_df, aes(x = .data$mid_x, y = .data$diff_prop)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, color = "grey55") +
    geom_col(width = bar_width_diff, fill = "steelblue4", alpha = 0.88, color = "grey25", linewidth = 0.12) +
    labs(
      x = xlab,
      y = "\u0394 bin proportion (1st \u2212 2nd in --compare)",
      title = paste0("Bin proportion difference (", compare_pair[1L], " \u2212 ", compare_pair[2L], ")")
    ) +
    theme_panel
  if (!is.null(xs)) p_diff <- p_diff + xs
}

plot_list <- list()
for (pname in panels_req) {
  plot_list[[length(plot_list) + 1L]] <- switch(pname,
    histogram = p_hist, density = p_dens, ecdf = p_ecdf, diff = p_diff)
}
n_p <- length(plot_list)
# Match previous 3-panel layout: total height opts$height for three stacked panels; one panel uses full opts$height
out_height <- if (n_p == 1L) opts$height else opts$height * n_p / 3
combined <- if (n_p == 1L) {
  p1 <- plot_list[[1L]]
  if (opts$overlay) p1 <- p1 + theme(legend.position = "right") else p1 <- p1 + theme(legend.position = "none")
  p1
} else {
  wrap_plots(plot_list, ncol = 1L) + plot_layout(guides = "collect") &
    theme(legend.position = if (opts$overlay) "right" else "none")
}

chr_safe <- if (!is.null(chr_region) && nzchar(chr_region)) gsub("[^A-Za-z0-9]", "_", chr_region) else "All"
pos_suffix <- if (!is.na(start_region) && !is.na(end_region)) paste0("_", chr_safe, "_", start_region, "_", end_region) else paste0("_", chr_safe)
stat_safe <- gsub("[^A-Za-z0-9]", "_", opts$stat)
overlay_suffix <- if (opts$overlay) "_overlay" else ""
tfm_suffix <- if (opts$transform != "none") paste0("_", opts$transform, "trans") else ""
yv_suffix <- if (opts$`y-value` != "value" && opts$stat %in% c("pi", "theta", "tajima_d", "fst", "pbe")) paste0("_", opts$`y-value`) else ""
panels_suffix <- if (length(panels_req) == 1L && panels_req[1] == "histogram") "" else paste0("_", paste(panels_req, collapse = "-"))
base_name <- paste0(opts$`file-prefix`, "hist_", stat_safe, pos_suffix, overlay_suffix, compare_suffix, tfm_suffix, yv_suffix, panels_suffix)

dpi_use <- if (!is.null(opts$dpi) && !is.na(opts$dpi)) opts$dpi else PLOT_DPI
for (fmt in format_parts) {
  out_path <- file.path(opts$`output-dir`, paste0(base_name, ".", fmt))
  if (fmt == "png") ggsave(out_path, combined, width = opts$width, height = out_height, dpi = dpi_use)
  else ggsave(out_path, combined, width = opts$width, height = out_height)
  message("Saved: ", out_path)
}

message("Done.")

