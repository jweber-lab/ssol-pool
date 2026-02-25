#!/usr/bin/env Rscript

###############################################################################
# plot_region.R
#
# One stacked figure for a genomic region: coverage, π, θ, Tajima's D, FST,
# PBE panels with shared x-axis and shared color key.
# Reads from TSV dirs and/or HDF5 (collate output).
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
  make_option("--chromosome", type = "character", default = NULL,
    help = "Chromosome to plot (full)"),
  make_option("--region", type = "character", default = NULL,
    help = "Region CHR:START-END"),
  make_option("--diversity-dir", type = "character", default = NULL,
    help = "Directory with diversity TSVs"),
  make_option("--fst-dir", type = "character", default = NULL,
    help = "Directory with FST TSVs"),
  make_option("--pbe-dir", type = "character", default = NULL,
    help = "Directory with PBE TSVs"),
  make_option("--seq-qual-dir", type = "character", default = NULL,
    help = "Directory with seq_qual TSVs"),
  make_option("--hdf5-dir", type = "character", default = NULL,
    help = "Collate output dir with .h5 files"),
  make_option("--window-size", type = "numeric", default = NULL,
    help = "Window size to use"),
  make_option("--step-size", type = "numeric", default = NULL,
    help = "Step size (optional)"),
  make_option("--reference-genome", type = "character", default = NULL,
    help = "Reference genome FASTA"),
  make_option("--output-dir", type = "character", default = "."),
  make_option("--file-prefix", type = "character", default = ""),
  make_option("--y-value", type = "character", default = "value",
    help = "Default y-axis for all stats: value, rank, or quantile [default: value]"),
  make_option("--statistics", type = "character", default = NULL,
    help = paste0("Comma-separated panels to include, in order. ",
      "Valid: coverage, pi, theta, tajima_d, fst, pbe. ",
      "Default: all that have data")),
  make_option("--transform", type = "character", default = "",
    help = paste0("Per-stat transforms as STAT:TRANSFORM pairs, comma-separated. ",
      "TRANSFORM is none, log, or asinh. Example: coverage:log,pi:none,fst:asinh. ",
      "Stats not listed use 'none'.")),
  make_option("--plot-style", type = "character", default = "line",
    help = "line or line_points [default: line]"),
  make_option("--width", type = "numeric", default = 12),
  make_option("--height", type = "numeric", default = 10),
  make_option("--dpi", type = "numeric", default = NULL),
  make_option("--plot-format", type = "character", default = "png")
)
opts <- parse_args(OptionParser(option_list = option_list, usage = "usage: %prog [options]"))

# ---------------------------------------------------------------------------
# Parse region
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
if (is.null(chr_region) || !nzchar(chr_region))
  stop("--chromosome or --region CHR:START-END is required")

# ---------------------------------------------------------------------------
# Read reference genome .fai (suppress "New names" by giving proper col_names)
# ---------------------------------------------------------------------------
chr_lengths <- tibble(chr = character(), length = numeric())
if (!is.null(opts$`reference-genome`) && nzchar(opts$`reference-genome`)) {
  fai <- paste0(opts$`reference-genome`, ".fai")
  if (file.exists(fai)) {
    chr_lengths <- read_tsv(fai,
      col_names = c("chr", "length", "offset", "linebases", "linewidth"),
      col_types = "cnnnn") %>%
      select(chr, length)
  }
}

dpi_use <- if (!is.null(opts$dpi) && !is.na(opts$dpi)) opts$dpi else PLOT_DPI

# ---------------------------------------------------------------------------
# Parse --statistics (panel order)
# ---------------------------------------------------------------------------
valid_stats <- c("coverage", "mapq", "pi", "theta", "tajima_d", "fst", "pbe")
requested_stats <- NULL
if (!is.null(opts$statistics) && nzchar(opts$statistics)) {
  requested_stats <- trimws(strsplit(opts$statistics, ",")[[1]])
  bad <- setdiff(requested_stats, valid_stats)
  if (length(bad) > 0) stop("Unknown --statistics: ", paste(bad, collapse = ", "),
    ". Valid: ", paste(valid_stats, collapse = ", "))
}

# ---------------------------------------------------------------------------
# Parse --transform (per-stat transforms)
# ---------------------------------------------------------------------------
stat_transforms <- setNames(rep("none", length(valid_stats)), valid_stats)
if (nzchar(opts$transform)) {
  pairs <- trimws(strsplit(opts$transform, ",")[[1]])
  for (p in pairs) {
    kv <- trimws(strsplit(p, ":")[[1]])
    if (length(kv) != 2) stop("Bad --transform entry: '", p, "'. Expected STAT:TRANSFORM")
    if (!kv[1] %in% valid_stats) stop("Unknown stat in --transform: ", kv[1])
    if (!kv[2] %in% c("none", "log", "asinh")) stop("Unknown transform: ", kv[2], ". Use none, log, or asinh")
    stat_transforms[kv[1]] <- kv[2]
  }
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
filter_region <- function(d) {
  d <- d %>% filter(.data$chr == chr_region)
  if (!is.na(start_region) && "end" %in% names(d))
    d <- d %>% filter(.data$end >= start_region)
  if (!is.na(end_region) && "start" %in% names(d))
    d <- d %>% filter(.data$start <= end_region)
  d
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

apply_transform <- function(values, tfm) {
  if (tfm == "log") {
    min_v <- min(values, na.rm = TRUE)
    if (min_v <= 0) log1p(pmax(0, values))
    else log(values)
  } else if (tfm == "asinh") {
    med <- median(values, na.rm = TRUE)
    s <- sd(values, na.rm = TRUE)
    if (is.na(s) || s == 0) s <- 1
    asinh((values - med) / s)
  } else {
    values
  }
}

# Build scale_y_continuous with ticks at nice original-scale values,
# positioned at their transformed coordinates, labeled with the original values.
build_y_scale <- function(original_values, tfm) {
  if (tfm == "none") return(NULL)

  orig_range <- range(original_values, na.rm = TRUE)
  if (any(!is.finite(orig_range))) return(NULL)

  min_v <- min(original_values, na.rm = TRUE)

  if (tfm == "log") {
    use_log1p <- (min_v <= 0)
    fwd <- if (use_log1p) function(x) log1p(pmax(0, x)) else log
    # Nice breaks in original space
    if (min_v > 0 && orig_range[2] / orig_range[1] > 10) {
      # Span multiple orders of magnitude: use log-spaced breaks
      log_lo <- floor(log10(max(orig_range[1], 1e-10)))
      log_hi <- ceiling(log10(orig_range[2]))
      orig_ticks <- 10^seq(log_lo, log_hi)
      orig_ticks <- orig_ticks[orig_ticks >= orig_range[1] * 0.9 & orig_ticks <= orig_range[2] * 1.1]
      if (length(orig_ticks) < 3) orig_ticks <- pretty(orig_range, n = 6)
    } else {
      orig_ticks <- pretty(orig_range, n = 6)
    }
    orig_ticks <- orig_ticks[orig_ticks >= 0 | !use_log1p]
    if (length(orig_ticks) < 2) orig_ticks <- pretty(orig_range, n = 6)
    y_breaks <- fwd(orig_ticks)

  } else if (tfm == "asinh") {
    med <- median(original_values, na.rm = TRUE)
    s <- sd(original_values, na.rm = TRUE)
    if (is.na(s) || s == 0) s <- 1
    fwd <- function(x) asinh((x - med) / s)
    orig_ticks <- pretty(orig_range, n = 6)
    y_breaks <- fwd(orig_ticks)

  } else {
    return(NULL)
  }

  # Format tick labels concisely
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

  scale_y_continuous(breaks = y_breaks, labels = fmt_label(orig_ticks))
}

y_label <- function(base_label, tfm) {
  # Keep original stat name; the inverse-labeled ticks make the transform clear
  base_label
}

# ---------------------------------------------------------------------------
# Load data from HDF5
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

  # --- Diversity HDF5 ---
  div_h5 <- list.files(hdir, pattern = "^diversity_.*\\.h5$", full.names = TRUE)
  if (length(div_h5) > 0) {
    if (!is.null(window_size_use))
      div_h5 <- div_h5[grepl(paste0("_w", window_size_use, "_"), div_h5, fixed = TRUE)]
    if (length(div_h5) > 0 && !is.null(step_size_use))
      div_h5 <- div_h5[grepl(paste0("_s", step_size_use, "\\."), div_h5, fixed = TRUE)]
    if (length(div_h5) > 0) {
      div_h5 <- div_h5[1]
      ws <- parse_window_step_from_filename(basename(div_h5))
      if (is.null(window_size_use)) window_size_use <- ws$window_size
      d <- read_h5_windows(div_h5, "windows")
      if (!is.null(d) && nrow(d) > 0) {
        d <- add_pos(d) %>% filter_region()
        if (nrow(d) > 0) {
          div_data <- d
          # Handle y-value mapping for each diversity stat
          for (stat in c("pi", "theta", "tajima_d")) {
            y_col <- switch(opts$`y-value`,
              rank = paste0(stat, "_rank"),
              quantile = paste0(stat, "_quantile"),
              stat)
            if (y_col %in% names(div_data) && y_col != stat)
              div_data[[stat]] <- div_data[[y_col]]
          }
          # Coverage from diversity HDF5 if present
          if ("mean_coverage" %in% names(div_data) && is.null(cov_data)) {
            cov_data <- div_data %>%
              select(any_of(c("chr", "start", "end", "sample", "mean_coverage", "pos"))) %>%
              filter(!is.na(mean_coverage))
          }
          # Mapping quality from diversity HDF5 if present
          if ("mean_mapping_quality" %in% names(div_data) && is.null(mapq_data)) {
            mapq_data <- div_data %>%
              select(any_of(c("chr", "start", "end", "sample", "mean_mapping_quality", "pos"))) %>%
              filter(!is.na(mean_mapping_quality))
          }
        }
      }
    }
  }

  # --- FST HDF5 ---
  fst_h5 <- list.files(hdir, pattern = "^fst_.*\\.h5$", full.names = TRUE)
  if (length(fst_h5) > 0) {
    if (!is.null(window_size_use))
      fst_h5 <- fst_h5[grepl(paste0("_w", window_size_use), fst_h5, fixed = TRUE)]
    if (length(fst_h5) > 0) {
      fst_h5 <- fst_h5[1]
      d <- read_h5_windows(fst_h5, "windows")
      if (is.null(d)) d <- read_h5_windows(fst_h5, "sites")
      if (!is.null(d) && nrow(d) > 0) {
        if (!"sample_pair" %in% names(d) && all(c("pop1", "pop2") %in% names(d)))
          d <- d %>% mutate(sample_pair = paste(pop1, pop2, sep = ":"))
        d <- add_pos(d) %>% filter_region()
        if (nrow(d) > 0) {
          fst_data <- d
          y_col <- switch(opts$`y-value`, rank = "fst_rank", quantile = "fst_quantile", "fst")
          if (y_col %in% names(fst_data) && y_col != "fst")
            fst_data[["fst"]] <- fst_data[[y_col]]
        }
      }
    }
  }

  # --- PBE HDF5 ---
  pbe_h5 <- list.files(hdir, pattern = "^pbe_.*\\.h5$", full.names = TRUE)
  if (length(pbe_h5) > 0) {
    if (!is.null(window_size_use))
      pbe_h5 <- pbe_h5[grepl(paste0("_w", window_size_use), pbe_h5, fixed = TRUE)]
    if (length(pbe_h5) > 0) {
      pbe_h5 <- pbe_h5[1]
      d <- read_h5_windows(pbe_h5, "windows")
      if (is.null(d)) d <- read_h5_windows(pbe_h5, "sites")
      if (!is.null(d) && nrow(d) > 0) {
        if (!"trio" %in% names(d) && all(c("pop1", "pop2", "pop3") %in% names(d)))
          d <- d %>% mutate(trio = paste(pop1, pop2, pop3, sep = ":"))
        d <- add_pos(d) %>% filter_region()
        if (nrow(d) > 0) {
          pbe_data <- d
          y_col <- switch(opts$`y-value`, rank = "pbe_rank", quantile = "pbe_quantile", "pbe")
          if (y_col %in% names(pbe_data) && y_col != "pbe")
            pbe_data[["pbe"]] <- pbe_data[[y_col]]
        }
      }
    }
  }
}

# ---------------------------------------------------------------------------
# TSV fallbacks
# ---------------------------------------------------------------------------
if ((is.null(cov_data) || is.null(mapq_data)) && !is.null(opts$`seq-qual-dir`) && dir.exists(opts$`seq-qual-dir`)) {
  sq_files <- list.files(opts$`seq-qual-dir`, pattern = "seq_qual.*\\.tsv$",
    full.names = TRUE, recursive = TRUE)
  if (length(sq_files) > 0) {
    sq_list <- lapply(sq_files, function(f) {
      x <- read_tsv(f, show_col_types = FALSE)
      if ("chromosome" %in% names(x) && !"chr" %in% names(x)) x <- x %>% rename(chr = chromosome)
      if (any(c("mean_coverage", "mean_mapping_quality") %in% names(x))) add_pos(x) else NULL
    })
    sq_list <- sq_list[!sapply(sq_list, is.null)]
    if (length(sq_list) > 0) {
      sq_combined <- bind_rows(sq_list) %>% filter_region()
      if (nrow(sq_combined) > 0 && !"pos" %in% names(sq_combined)) sq_combined <- add_pos(sq_combined)
      if (is.null(cov_data) && "mean_coverage" %in% names(sq_combined))
        cov_data <- sq_combined %>% filter(!is.na(mean_coverage))
      if (is.null(mapq_data) && "mean_mapping_quality" %in% names(sq_combined))
        mapq_data <- sq_combined %>% filter(!is.na(mean_mapping_quality))
    }
  }
}

if (is.null(div_data) && !is.null(opts$`diversity-dir`) && dir.exists(opts$`diversity-dir`)) {
  div_tsv <- list.files(opts$`diversity-dir`, pattern = "diversity.*\\.tsv$", full.names = TRUE)
  if (length(div_tsv) > 0) {
    d <- read_tsv(div_tsv[1], show_col_types = FALSE)
    if ("chromosome" %in% names(d) && !"chr" %in% names(d)) d <- d %>% rename(chr = chromosome)
    if ("position" %in% names(d) && !"pos" %in% names(d)) d <- d %>% rename(pos = position)
    if ("pi" %in% names(d)) {
      if (!"sample" %in% names(d)) d <- d %>% mutate(sample = "sample")
      div_data <- add_pos(d) %>% filter_region()
    }
  }
}

if (is.null(fst_data) && !is.null(opts$`fst-dir`) && dir.exists(opts$`fst-dir`)) {
  fst_tsv <- list.files(opts$`fst-dir`, pattern = "fst.*\\.tsv$", full.names = TRUE)
  if (length(fst_tsv) > 0) {
    d <- read_tsv(fst_tsv[1], show_col_types = FALSE)
    if ("chromosome" %in% names(d) && !"chr" %in% names(d)) d <- d %>% rename(chr = chromosome)
    if ("position" %in% names(d) && !"pos" %in% names(d)) d <- d %>% rename(pos = position)
    fst_col <- names(d)[grepl("\\.fst$", names(d), ignore.case = TRUE)][1]
    if (!is.na(fst_col)) {
      d <- d %>% rename(fst = all_of(fst_col))
      if (!"sample_pair" %in% names(d)) d <- d %>% mutate(sample_pair = sub("\\.fst$", "", fst_col))
      fst_data <- add_pos(d) %>% filter_region()
    }
  }
}

if (is.null(pbe_data) && !is.null(opts$`pbe-dir`) && dir.exists(opts$`pbe-dir`)) {
  pbe_tsv <- list.files(opts$`pbe-dir`, pattern = "pbe.*\\.tsv$", full.names = TRUE)
  if (length(pbe_tsv) > 0) {
    d <- read_tsv(pbe_tsv[1], show_col_types = FALSE)
    if ("chromosome" %in% names(d) && !"chr" %in% names(d)) d <- d %>% rename(chr = chromosome)
    if ("position" %in% names(d) && !"pos" %in% names(d)) d <- d %>% rename(pos = position)
    pbe_col <- names(d)[grepl("pbe", names(d), ignore.case = TRUE)][1]
    if (!is.na(pbe_col)) {
      d <- d %>% rename(pbe = all_of(pbe_col))
      if (!"trio" %in% names(d)) d <- d %>% mutate(trio = "trio")
      pbe_data <- add_pos(d) %>% filter_region()
    }
  }
}

# ---------------------------------------------------------------------------
# Determine which panels have data
# ---------------------------------------------------------------------------
has_data <- c(
  coverage  = !is.null(cov_data) && nrow(cov_data) > 0 && "mean_coverage" %in% names(cov_data),
  mapq      = !is.null(mapq_data) && nrow(mapq_data) > 0 && "mean_mapping_quality" %in% names(mapq_data),
  pi        = !is.null(div_data) && nrow(div_data) > 0 && "pi" %in% names(div_data),
  theta     = !is.null(div_data) && nrow(div_data) > 0 && "theta" %in% names(div_data),
  tajima_d  = !is.null(div_data) && nrow(div_data) > 0 && "tajima_d" %in% names(div_data),
  fst       = !is.null(fst_data) && nrow(fst_data) > 0 && "fst" %in% names(fst_data),
  pbe       = !is.null(pbe_data) && nrow(pbe_data) > 0 && "pbe" %in% names(pbe_data)
)

if (is.null(requested_stats)) {
  panels_to_plot <- valid_stats[has_data[valid_stats]]
} else {
  panels_to_plot <- requested_stats[has_data[requested_stats]]
  skipped <- requested_stats[!has_data[requested_stats]]
  if (length(skipped) > 0) message("No data for: ", paste(skipped, collapse = ", "))
}

if (length(panels_to_plot) == 0)
  stop("No data found for the specified region and inputs. Check --chromosome/--region and input dirs.")

# ---------------------------------------------------------------------------
# Shared theme (no x-axis title; added only to bottom panel)
# ---------------------------------------------------------------------------
theme_panel <- theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
  theme(panel.grid.minor = element_blank())

use_points <- (opts$`plot-style` == "line_points")

# ---------------------------------------------------------------------------
# Build one panel per stat
# ---------------------------------------------------------------------------
panels <- list()

for (stat_name in panels_to_plot) {
  tfm <- stat_transforms[stat_name]

  if (stat_name == "coverage") {
    if (!"sample" %in% names(cov_data)) cov_data <- cov_data %>% mutate(sample = "sample")
    samples <- sort(unique(cov_data$sample))
    pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(samples)), samples)
    orig_vals <- cov_data$mean_coverage
    cov_data <- cov_data %>% mutate(y_val = apply_transform(mean_coverage, tfm))
    p <- ggplot(cov_data, aes(x = pos, y = y_val, color = sample, group = sample)) +
      geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
      scale_color_manual(name = "Sample", values = pal, drop = FALSE) +
      labs(y = y_label("Mean coverage", tfm)) + theme_panel
    ys <- build_y_scale(orig_vals, tfm); if (!is.null(ys)) p <- p + ys
    if (use_points) p <- p + geom_point(size = 1, alpha = 0.7, na.rm = TRUE)
    panels[[length(panels) + 1]] <- p

  } else if (stat_name == "mapq") {
    if (!"sample" %in% names(mapq_data)) mapq_data <- mapq_data %>% mutate(sample = "sample")
    samples <- sort(unique(mapq_data$sample))
    pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(samples)), samples)
    orig_vals <- mapq_data$mean_mapping_quality
    mapq_data <- mapq_data %>% mutate(y_val = apply_transform(mean_mapping_quality, tfm))
    p <- ggplot(mapq_data, aes(x = pos, y = y_val, color = sample, group = sample)) +
      geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
      scale_color_manual(name = "Sample", values = pal, drop = FALSE) +
      labs(y = y_label("Mean MAPQ", tfm)) + theme_panel
    ys <- build_y_scale(orig_vals, tfm); if (!is.null(ys)) p <- p + ys
    if (use_points) p <- p + geom_point(size = 1, alpha = 0.7, na.rm = TRUE)
    panels[[length(panels) + 1]] <- p

  } else if (stat_name %in% c("pi", "theta", "tajima_d")) {
    if (!"sample" %in% names(div_data)) div_data <- div_data %>% mutate(sample = "sample")
    if (!stat_name %in% names(div_data)) next
    samples <- sort(unique(div_data$sample))
    pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(samples)), samples)
    orig_vals <- div_data[[stat_name]]
    plot_df <- div_data %>%
      filter(!is.na(.data[[stat_name]])) %>%
      mutate(y_val = apply_transform(.data[[stat_name]], tfm))
    base_lbl <- switch(stat_name, pi = "\u03c0", theta = "\u03b8", tajima_d = "Tajima's D")
    if (opts$`y-value` == "rank") base_lbl <- paste(base_lbl, "rank")
    else if (opts$`y-value` == "quantile") base_lbl <- paste(base_lbl, "quantile")
    p <- ggplot(plot_df, aes(x = pos, y = y_val, color = sample, group = sample)) +
      geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
      scale_color_manual(name = "Sample", values = pal, drop = FALSE) +
      labs(y = y_label(base_lbl, tfm)) + theme_panel
    ys <- build_y_scale(orig_vals[!is.na(orig_vals)], tfm); if (!is.null(ys)) p <- p + ys
    if (use_points) p <- p + geom_point(size = 1, alpha = 0.7, na.rm = TRUE)
    panels[[length(panels) + 1]] <- p

  } else if (stat_name == "fst") {
    if (!"sample_pair" %in% names(fst_data)) next
    pairs <- sort(unique(fst_data$sample_pair))
    pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(pairs)), pairs)
    orig_vals <- fst_data$fst
    fst_data <- fst_data %>% mutate(y_val = apply_transform(fst, tfm))
    base_lbl <- if (opts$`y-value` == "rank") "FST rank" else if (opts$`y-value` == "quantile") "FST quantile" else "FST"
    p <- ggplot(fst_data, aes(x = pos, y = y_val, color = sample_pair, group = sample_pair)) +
      geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
      scale_color_manual(name = "Pair", values = pal, drop = FALSE) +
      labs(y = y_label(base_lbl, tfm)) + theme_panel
    ys <- build_y_scale(orig_vals[!is.na(orig_vals)], tfm); if (!is.null(ys)) p <- p + ys
    if (use_points) p <- p + geom_point(size = 1, alpha = 0.7, na.rm = TRUE)
    panels[[length(panels) + 1]] <- p

  } else if (stat_name == "pbe") {
    if (!"trio" %in% names(pbe_data)) next
    trios <- sort(unique(pbe_data$trio))
    pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(trios)), trios)
    orig_vals <- pbe_data$pbe
    pbe_data <- pbe_data %>% mutate(y_val = apply_transform(pbe, tfm))
    base_lbl <- if (opts$`y-value` == "rank") "PBE rank" else if (opts$`y-value` == "quantile") "PBE quantile" else "PBE"
    p <- ggplot(pbe_data, aes(x = pos, y = y_val, color = trio, group = trio)) +
      geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
      scale_color_manual(name = "Trio", values = pal, drop = FALSE) +
      labs(y = y_label(base_lbl, tfm)) + theme_panel
    ys <- build_y_scale(orig_vals[!is.na(orig_vals)], tfm); if (!is.null(ys)) p <- p + ys
    if (use_points) p <- p + geom_point(size = 1, alpha = 0.7, na.rm = TRUE)
    panels[[length(panels) + 1]] <- p
  }
}

if (length(panels) == 0)
  stop("No panels could be built. Check data availability for requested statistics.")

# ---------------------------------------------------------------------------
# X-axis: label bottom panel with scaffold; hide title + ticks on upper panels
# ---------------------------------------------------------------------------
x_label <- paste0("Position on ", chr_region, " (bp)")
n_panels <- length(panels)
for (i in seq_len(n_panels)) {
  if (i < n_panels) {
    panels[[i]] <- panels[[i]] +
      labs(x = NULL) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    panels[[i]] <- panels[[i]] + labs(x = x_label)
  }
}

# ---------------------------------------------------------------------------
# Combine with patchwork
# ---------------------------------------------------------------------------
combined <- wrap_plots(panels, ncol = 1, heights = rep(1, n_panels)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "right",
    plot.tag = element_text(size = PLOT_BASE_SIZE + 2, face = "bold"),
    plot.tag.position = "topleft")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
chr_safe <- gsub("[^A-Za-z0-9]", "_", chr_region)
pos_suffix <- if (!is.na(start_region) && !is.na(end_region))
  paste0("_", chr_safe, "_", start_region, "_", end_region) else paste0("_", chr_safe)
base_name <- paste0(opts$`file-prefix`, "region_plot", pos_suffix)

format_parts <- trimws(tolower(strsplit(opts$`plot-format`, ",")[[1]]))
if ("both" %in% format_parts) format_parts <- c(setdiff(format_parts, "both"), "png", "pdf")
if ("all" %in% format_parts) format_parts <- c(setdiff(format_parts, "all"), "png", "pdf", "svg")
for (fmt in unique(format_parts)) {
  out_path <- file.path(opts$`output-dir`, paste0(base_name, ".", fmt))
  if (fmt == "png") {
    ggsave(out_path, combined, width = opts$width, height = opts$height, dpi = dpi_use)
  } else {
    ggsave(out_path, combined, width = opts$width, height = opts$height)
  }
  message("Saved: ", out_path)
}

message("Done.")
