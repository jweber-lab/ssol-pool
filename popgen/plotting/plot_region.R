#!/usr/bin/env Rscript

###############################################################################
# plot_region.R
#
# One stacked figure for a genomic region: coverage, π, FST, PBE panels with
# shared x-axis and shared color key for sample/pair/trio.
# Reads from TSV dirs and/or HDF5 (collate output).
#
# Requires: ggplot2, dplyr, tidyr, readr, optparse, patchwork (for combining).
# Optional: hdf5r when using --hdf5-dir.
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

# Source plot_common (read_h5_windows, PLOT_PALETTE_QUALITATIVE, PLOT_DPI, etc.)
initial_args <- commandArgs(trailingOnly = FALSE)
file_arg <- initial_args[substr(initial_args, 1, 7) == "--file="]
script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg)) else "."
plot_common_path <- file.path(script_dir, "plot_common.R")
if (!file.exists(plot_common_path)) plot_common_path <- file.path(getwd(), "plot_common.R")
if (!file.exists(plot_common_path)) stop("plot_common.R not found in script dir or getwd()")
source(plot_common_path)

option_list <- list(
  make_option(c("--chromosome"), type = "character", default = NULL, help = "Chromosome to plot (full)"),
  make_option(c("--region"), type = "character", default = NULL, help = "Region CHR:START-END"),
  make_option(c("--diversity-dir"), type = "character", default = NULL),
  make_option(c("--fst-dir"), type = "character", default = NULL),
  make_option(c("--pbe-dir"), type = "character", default = NULL),
  make_option(c("--seq-qual-dir"), type = "character", default = NULL),
  make_option(c("--hdf5-dir"), type = "character", default = NULL),
  make_option(c("--window-size"), type = "numeric", default = NULL),
  make_option(c("--step-size"), type = "numeric", default = NULL),
  make_option(c("--reference-genome"), type = "character", default = NULL),
  make_option(c("--output-dir"), type = "character", default = "."),
  make_option(c("--file-prefix"), type = "character", default = ""),
  make_option(c("--y-value"), type = "character", default = "value", help = "value, rank, or quantile"),
  make_option(c("--width"), type = "numeric", default = 12),
  make_option(c("--height"), type = "numeric", default = 10),
  make_option(c("--dpi"), type = "numeric", default = NULL),
  make_option(c("--plot-format"), type = "character", default = "png")
)
opts <- parse_args(OptionParser(option_list = option_list, usage = "usage: %prog [options]"))

# Parse region
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
if (is.null(chr_region) || !nzchar(chr_region)) stop("--chromosome or --region CHR:START-END is required")

# Chromosome lengths (optional)
chr_lengths <- tibble(chr = character(), length = numeric())
if (!is.null(opts$`reference-genome`) && nzchar(opts$`reference-genome`)) {
  fai <- paste0(opts$`reference-genome`, ".fai")
  if (file.exists(fai)) {
    chr_lengths <- readr::read_tsv(fai, col_names = c("chr", "length", NA, NA, NA), show_col_types = FALSE)
  }
}

dpi_use <- if (!is.null(opts$dpi) && !is.na(opts$dpi)) opts$dpi else PLOT_DPI

# Filter to region helper
filter_region <- function(d, chr_col = "chr", start_col = "start", end_col = "end") {
  d <- d %>% filter(.data[[chr_col]] == chr_region)
  if (!is.na(start_region) && start_col %in% names(d))
    d <- d %>% filter(.data[[end_col]] >= start_region)
  if (!is.na(end_region) && end_col %in% names(d))
    d <- d %>% filter(.data[[start_col]] <= end_region)
  d
}

# Shared x: use mid position or start for alignment
add_pos <- function(d) {
  if ("start" %in% names(d) && "end" %in% names(d))
    d <- d %>% mutate(pos = (.data$start + .data$end) / 2)
  else if ("pos" %in% names(d)) { }
  else if ("position" %in% names(d)) d <- d %>% rename(pos = .data$position)
  else if ("start" %in% names(d)) d <- d %>% mutate(pos = .data$start)
  d
}

# --- Load data: HDF5 path takes precedence when provided ---
div_data <- NULL
fst_data <- NULL
pbe_data <- NULL
cov_data <- NULL
window_size_use <- opts$`window-size`
step_size_use <- opts$`step-size`

if (!is.null(opts$`hdf5-dir`) && dir.exists(opts$`hdf5-dir`)) {
  hdir <- opts$`hdf5-dir`
  # Diversity
  div_h5 <- list.files(hdir, pattern = "^diversity_.*\\.h5$", full.names = TRUE)
  if (length(div_h5) > 0) {
    # Match window/step if specified
    if (!is.null(window_size_use)) {
      div_h5 <- div_h5[grepl(paste0("_w", window_size_use, "_"), div_h5, fixed = TRUE)]
      if (length(div_h5) > 0 && !is.null(step_size_use))
        div_h5 <- div_h5[grepl(paste0("_s", step_size_use, "\\."), div_h5, fixed = TRUE)]
    }
    if (length(div_h5) > 0) {
      div_h5 <- div_h5[1]
      ws <- parse_window_step_from_filename(basename(div_h5))
      if (is.null(window_size_use)) window_size_use <- ws$window_size
      d <- read_h5_windows(div_h5, "windows")
      if (!is.null(d) && nrow(d) > 0) {
        d <- d %>% mutate(window_size = ws$window_size, step_size = ws$step_size)
        div_data <- add_pos(d) %>% filter_region()
        y_col <- switch(opts$`y-value`, rank = "pi_rank", quantile = "pi_quantile", "pi")
        if (y_col %in% names(div_data)) div_data <- div_data %>% rename(pi = .data[[y_col]])
        if ("mean_coverage" %in% names(div_data) && is.null(cov_data)) {
          cov_data <- div_data %>% select(.data$chr, .data$start, .data$end, .data$sample, .data$mean_coverage, .data$pos) %>%
            filter(!is.na(.data$mean_coverage))
        }
      }
    }
  }
  # FST
  fst_h5 <- list.files(hdir, pattern = "^fst_.*\\.h5$", full.names = TRUE)
  if (length(fst_h5) > 0) {
    if (!is.null(window_size_use)) fst_h5 <- fst_h5[grepl(paste0("_w", window_size_use), fst_h5, fixed = TRUE)]
    if (length(fst_h5) > 0) {
      fst_h5 <- fst_h5[1]
      d <- read_h5_windows(fst_h5, "windows")
      if (is.null(d) && "sites" %in% names(tryCatch(hdf5r::H5File$new(fst_h5, "r"), error = function(e) NULL)))
        d <- read_h5_windows(fst_h5, "sites")
      if (!is.null(d) && nrow(d) > 0) {
        if (!"sample_pair" %in% names(d) && all(c("pop1", "pop2") %in% names(d)))
          d <- d %>% mutate(sample_pair = paste(.data$pop1, .data$pop2, sep = ":"))
        fst_data <- add_pos(d) %>% filter_region()
        y_col <- switch(opts$`y-value`, rank = "fst_rank", quantile = "fst_quantile", "fst")
        if (y_col %in% names(fst_data)) fst_data <- fst_data %>% rename(fst = .data[[y_col]])
      }
    }
  }
  # PBE
  pbe_h5 <- list.files(hdir, pattern = "^pbe_.*\\.h5$", full.names = TRUE)
  if (length(pbe_h5) > 0) {
    if (!is.null(window_size_use)) pbe_h5 <- pbe_h5[grepl(paste0("_w", window_size_use), pbe_h5, fixed = TRUE)]
    if (length(pbe_h5) > 0) {
      pbe_h5 <- pbe_h5[1]
      d <- read_h5_windows(pbe_h5, "windows")
      if (is.null(d)) d <- read_h5_windows(pbe_h5, "sites")
      if (!is.null(d) && nrow(d) > 0) {
        if (!"trio" %in% names(d) && all(c("pop1", "pop2", "pop3") %in% names(d)))
          d <- d %>% mutate(trio = paste(.data$pop1, .data$pop2, .data$pop3, sep = ":"))
        pbe_data <- add_pos(d) %>% filter_region()
        y_col <- switch(opts$`y-value`, rank = "pbe_rank", quantile = "pbe_quantile", "pbe")
        if (y_col %in% names(pbe_data)) pbe_data <- pbe_data %>% rename(pbe = .data[[y_col]])
      }
    }
  }
}

# TSV fallbacks
if (is.null(cov_data) && !is.null(opts$`seq-qual-dir`) && dir.exists(opts$`seq-qual-dir`)) {
  sq_files <- list.files(opts$`seq-qual-dir`, pattern = "seq_qual.*\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (length(sq_files) > 0) {
    sq_list <- lapply(sq_files, function(f) {
      x <- readr::read_tsv(f, show_col_types = FALSE)
      x <- x %>% rename(chr = any_of(c("chromosome", "chr")), start = any_of("start"), end = any_of("end"))
      if ("mean_coverage" %in% names(x)) add_pos(x) else NULL
    })
    sq_list <- sq_list[!sapply(sq_list, is.null)]
    if (length(sq_list) > 0) {
      cov_data <- bind_rows(sq_list) %>% filter_region()
      if (nrow(cov_data) > 0 && !"pos" %in% names(cov_data)) cov_data <- add_pos(cov_data)
    }
  }
}

if (is.null(div_data) && !is.null(opts$`diversity-dir`) && dir.exists(opts$`diversity-dir`)) {
  div_tsv <- list.files(opts$`diversity-dir`, pattern = "diversity.*\\.tsv$", full.names = TRUE)
  if (length(div_tsv) > 0) {
    d <- readr::read_tsv(div_tsv[1], show_col_types = FALSE)
    d <- d %>% rename(chr = any_of(c("chromosome", "chr")), pos = any_of(c("position", "pos", "start")))
    if ("pi" %in% names(d)) {
      if (!"sample" %in% names(d)) d <- d %>% mutate(sample = "sample")
      div_data <- add_pos(d) %>% filter_region()
    }
  }
}

if (is.null(fst_data) && !is.null(opts$`fst-dir`) && dir.exists(opts$`fst-dir`)) {
  fst_tsv <- list.files(opts$`fst-dir`, pattern = "fst.*\\.tsv$", full.names = TRUE)
  if (length(fst_tsv) > 0) {
    d <- readr::read_tsv(fst_tsv[1], show_col_types = FALSE)
    d <- d %>% rename(chr = any_of(c("chromosome", "chr")), pos = any_of(c("position", "pos", "start")))
    fst_col <- names(d)[grepl("\\.fst$", names(d), ignore.case = TRUE)][1]
    if (!is.na(fst_col)) {
      d <- d %>% rename(fst = .data[[fst_col]])
      if (!"sample_pair" %in% names(d)) d <- d %>% mutate(sample_pair = sub("\\.fst$", "", fst_col))
      fst_data <- add_pos(d) %>% filter_region()
    }
  }
}

if (is.null(pbe_data) && !is.null(opts$`pbe-dir`) && dir.exists(opts$`pbe-dir`)) {
  pbe_tsv <- list.files(opts$`pbe-dir`, pattern = "pbe.*\\.tsv$", full.names = TRUE)
  if (length(pbe_tsv) > 0) {
    d <- readr::read_tsv(pbe_tsv[1], show_col_types = FALSE)
    d <- d %>% rename(chr = any_of(c("chromosome", "chr")), pos = any_of(c("position", "pos", "start")))
    pbe_col <- names(d)[grepl("pbe", names(d), ignore.case = TRUE)][1]
    if (!is.na(pbe_col)) {
      d <- d %>% rename(pbe = .data[[pbe_col]])
      if (!"trio" %in% names(d)) d <- d %>% mutate(trio = "trio")
      pbe_data <- add_pos(d) %>% filter_region()
    }
  }
}

# Build panels (only include if data present)
panels <- list()

theme_panel <- theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank())

# Coverage
if (!is.null(cov_data) && nrow(cov_data) > 0 && "sample" %in% names(cov_data)) {
  samples <- sort(unique(cov_data$sample))
  pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(samples)), samples)
  p_cov <- ggplot(cov_data, aes(x = .data$pos, y = .data$mean_coverage, color = .data$sample, group = .data$sample)) +
    geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
    scale_color_manual(name = "Sample", values = pal, drop = FALSE) +
    labs(y = "Mean coverage") +
    theme_panel
  panels[[length(panels) + 1]] <- p_cov
}

# π
if (!is.null(div_data) && nrow(div_data) > 0 && "pi" %in% names(div_data)) {
  if (!"sample" %in% names(div_data)) div_data <- div_data %>% mutate(sample = "sample")
  samples <- sort(unique(div_data$sample))
  pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(samples)), samples)
  p_pi <- ggplot(div_data, aes(x = .data$pos, y = .data$pi, color = .data$sample, group = .data$sample)) +
    geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
    scale_color_manual(name = "Sample", values = pal, drop = FALSE) +
    labs(y = if (opts$`y-value` == "rank") "π rank" else if (opts$`y-value` == "quantile") "π quantile" else "π") +
    theme_panel
  panels[[length(panels) + 1]] <- p_pi
}

# FST
if (!is.null(fst_data) && nrow(fst_data) > 0 && "fst" %in% names(fst_data)) {
  pairs <- sort(unique(fst_data$sample_pair))
  pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(pairs)), pairs)
  p_fst <- ggplot(fst_data, aes(x = .data$pos, y = .data$fst, color = .data$sample_pair, group = .data$sample_pair)) +
    geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
    scale_color_manual(name = "Pair", values = pal, drop = FALSE) +
    labs(y = if (opts$`y-value` == "rank") "FST rank" else if (opts$`y-value` == "quantile") "FST quantile" else "FST") +
    theme_panel
  panels[[length(panels) + 1]] <- p_fst
}

# PBE
if (!is.null(pbe_data) && nrow(pbe_data) > 0 && "pbe" %in% names(pbe_data)) {
  trios <- sort(unique(pbe_data$trio))
  pal <- setNames(rep(PLOT_PALETTE_QUALITATIVE, length.out = length(trios)), trios)
  p_pbe <- ggplot(pbe_data, aes(x = .data$pos, y = .data$pbe, color = .data$trio, group = .data$trio)) +
    geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
    scale_color_manual(name = "Trio", values = pal, drop = FALSE) +
    labs(y = if (opts$`y-value` == "rank") "PBE rank" else if (opts$`y-value` == "quantile") "PBE quantile" else "PBE") +
    theme_panel
  panels[[length(panels) + 1]] <- p_pbe
}

if (length(panels) == 0) {
  stop("No data found for the specified region and inputs. Check --chromosome/--region and input dirs.")
}

# Add x-axis label to bottom panel only
panels[[length(panels)]] <- panels[[length(panels)]] + labs(x = "Position (bp)")
# Combine: stack vertically, collect legends, panel labels (a), (b), ...
combined <- wrap_plots(panels, ncol = 1, heights = rep(1, length(panels))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "right", plot.tag = element_text(size = PLOT_BASE_SIZE + 2, face = "bold"), plot.tag.position = "topleft")

# Position label
chr_safe <- gsub("[^A-Za-z0-9]", "_", chr_region)
pos_suffix <- if (!is.na(start_region) && !is.na(end_region)) paste0("_", chr_safe, "_", start_region, "_", end_region) else paste0("_", chr_safe)
base_name <- paste0(opts$`file-prefix`, "region_plot", pos_suffix)

format_parts <- trimws(tolower(strsplit(opts$`plot-format`, ",")[[1]]))
if ("both" %in% format_parts) format_parts <- c(setdiff(format_parts, "both"), "png", "pdf")
if ("all" %in% format_parts) format_parts <- c(setdiff(format_parts, "all"), "png", "pdf", "svg")
for (fmt in unique(format_parts)) {
  if (fmt == "png") {
    path <- file.path(opts$`output-dir`, paste0(base_name, ".png"))
    ggsave(path, combined, width = opts$width, height = opts$height, dpi = dpi_use)
    message("Saved: ", path)
  } else if (fmt == "pdf") {
    path <- file.path(opts$`output-dir`, paste0(base_name, ".pdf"))
    ggsave(path, combined, width = opts$width, height = opts$height)
    message("Saved: ", path)
  } else if (fmt == "svg") {
    path <- file.path(opts$`output-dir`, paste0(base_name, ".svg"))
    ggsave(path, combined, width = opts$width, height = opts$height)
    message("Saved: ", path)
  }
}

message("Done.")
