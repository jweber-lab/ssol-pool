#!/usr/bin/env Rscript

###############################################################################
# plot_outlier_csq.R
#
# Consequence composition of sites in outlier regions, per outlier statistic, against the
# genome-wide background. Inputs come from outliers/outlier_annotation.sh --variant-csq:
#   outlier_csq_sites_<suffix>.csv       (one row per region x site, with outlier_stat)
#   outlier_csq_background_<suffix>.csv  (genome-wide counts; all_sites and non_outlier)
#
# A region's outlier_stat may list several statistics ("fst:pbe:pi"); its sites count toward
# each of them. Each site is counted once per statistic even if it falls in several regions.
#
# Outputs (in --output-dir, <cat> = --category):
#   outlier_csq_composition_<cat>.png  Stacked proportions per statistic + background bars,
#                                      faceted by variant type; n sites above each bar
#   outlier_csq_enrichment_<cat>.png   log2 odds ratio (outlier-statistic sites vs non-outlier
#                                      sites) per category, 95% CI; filled points = BH q < --alpha
#   outlier_csq_enrichment_<cat>.csv   Counts, odds ratio, CI, Fisher p, BH q per stat x type x category
###############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

initial_args <- commandArgs(trailingOnly = FALSE)
file_arg <- initial_args[substr(initial_args, 1, 7) == "--file="]
script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg)) else "."
plot_common_path <- file.path(script_dir, "plot_common.R")
if (!file.exists(plot_common_path)) plot_common_path <- file.path(getwd(), "plot_common.R")
if (!file.exists(plot_common_path)) stop("plot_common.R not found in script dir or getwd()")
source(plot_common_path)

option_list <- list(
  make_option("--sites", type = "character", default = NULL,
    help = "outlier_csq_sites_*.csv from outlier_annotation.sh (required unless --annotation-dir)"),
  make_option("--background", type = "character", default = NULL,
    help = "outlier_csq_background_*.csv [default: same suffix as --sites]"),
  make_option("--annotation-dir", type = "character", default = NULL,
    help = "outlier_annotation output dir; uses its single outlier_csq_sites_*.csv"),
  make_option("--category", type = "character", default = "coding_effect",
    help = "coding_effect or region_category [default: %default]"),
  make_option("--variant-type", type = "character", default = "split",
    help = "split (facet SNP/indel), all (pool), snp, or indel [default: %default]"),
  make_option("--drop-none", action = "store_true", default = FALSE,
    help = "Drop sites with no consequence record (category 'none') before computing proportions"),
  make_option("--min-sites", type = "integer", default = 5L,
    help = "Minimum sites in a statistic group for its enrichment row [default: %default]"),
  make_option("--alpha", type = "double", default = 0.05, help = "BH q threshold for filled points [default: %default]"),
  make_option(c("--output-dir", "-o"), type = "character", default = "./", help = "Output directory"),
  make_option("--file-prefix", type = "character", default = "", help = "Prefix for output file names"),
  make_option("--width", type = "double", default = 10, help = "Width in inches [default: %default]"),
  make_option("--height", type = "double", default = 6, help = "Height in inches [default: %default]"),
  make_option("--dpi", type = "integer", default = PLOT_DPI, help = "PNG DPI [default: %default]"),
  make_option("--plot-format", type = "character", default = "png", help = "png, pdf, svg, both, or all [default: %default]")
)
opts <- parse_args(OptionParser(option_list = option_list))

category <- match.arg(opts$category, c("coding_effect", "region_category"))
vt_mode <- match.arg(opts$`variant-type`, c("split", "all", "snp", "indel"))

sites_path <- opts$sites
if (is.null(sites_path) && !is.null(opts$`annotation-dir`)) {
  cand <- list.files(opts$`annotation-dir`, pattern = "^outlier_csq_sites_.*\\.csv$", full.names = TRUE)
  if (length(cand) != 1) stop("Expected one outlier_csq_sites_*.csv in ", opts$`annotation-dir`, ", found ", length(cand), "; use --sites")
  sites_path <- cand
}
if (is.null(sites_path) || !file.exists(sites_path)) stop("--sites (or --annotation-dir) is required")
bg_path <- opts$background
if (is.null(bg_path)) bg_path <- file.path(dirname(sites_path), sub("^outlier_csq_sites_", "outlier_csq_background_", basename(sites_path)))
if (!file.exists(bg_path)) stop("Background counts not found: ", bg_path, " (from outlier_annotation.sh --variant-csq)")

sites <- read_csv(sites_path, col_types = cols(.default = col_character()), show_col_types = FALSE)
bg <- read_csv(bg_path, show_col_types = FALSE, col_types = cols(n = col_double(), .default = col_character()))
if (!"outlier_stat" %in% names(sites)) stop("Sites table has no outlier_stat column (regions file from identify_outliers needed)")
if (!category %in% names(sites)) stop("Sites table has no ", category, " column")

stat_sites <- sites %>%
  mutate(stat = strsplit(coalesce(.data$outlier_stat, "unknown"), ":", fixed = TRUE)) %>%
  unnest("stat") %>%
  distinct(.data$stat, .data$chr, .data$pos, .data$ref, .data$alt, .keep_all = TRUE) %>%
  transmute(group = .data$stat, variant_type = .data$variant_type, cat = coalesce(.data[[category]], "none"))
bg_counts <- bg %>%
  transmute(set = .data$set, variant_type = .data$variant_type, cat = coalesce(.data[[category]], "none"), n = .data$n)

apply_vt <- function(d) {
  if (vt_mode == "all") d$variant_type <- "all sites"
  else if (vt_mode %in% c("snp", "indel")) d <- d %>% filter(.data$variant_type == vt_mode)
  d
}
stat_sites <- apply_vt(stat_sites)
bg_counts <- apply_vt(bg_counts)
if (opts$`drop-none`) {
  stat_sites <- stat_sites %>% filter(.data$cat != "none")
  bg_counts <- bg_counts %>% filter(.data$cat != "none")
}
if (nrow(stat_sites) == 0) stop("No outlier sites left to plot")

stat_counts <- stat_sites %>% count(.data$group, .data$variant_type, .data$cat, name = "n")
bg_labels <- c(all_sites = "all sites", non_outlier = "non-outlier sites")
comp <- bind_rows(
  stat_counts,
  bg_counts %>% mutate(group = bg_labels[.data$set]) %>% group_by(.data$group, .data$variant_type, .data$cat) %>%
    summarise(n = sum(.data$n), .groups = "drop")
) %>%
  group_by(.data$group, .data$variant_type) %>%
  mutate(total = sum(.data$n), prop = .data$n / .data$total) %>%
  ungroup()
stat_levels <- sort(unique(stat_counts$group))
comp$group <- factor(comp$group, levels = c(stat_levels, unname(bg_labels)))
comp$cat <- csq_factor(as.character(comp$cat), category)
fill_vals <- csq_fill_values(levels(comp$cat))

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)
save_plot <- function(p, stem, width = opts$width, height = opts$height) {
  for (fmt in plot_formats(opts$`plot-format`)) {
    out <- file.path(opts$`output-dir`, paste0(opts$`file-prefix`, stem, ".", fmt))
    if (fmt == "png") ggsave(out, p, width = width, height = height, dpi = opts$dpi) else ggsave(out, p, width = width, height = height)
    message("Saved: ", out)
  }
}

totals <- comp %>% distinct(.data$group, .data$variant_type, .data$total)
cat_label <- if (category == "coding_effect") "Coding effect" else "Region"
p_comp <- ggplot(comp, aes(x = .data$group, y = .data$prop, fill = .data$cat)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  geom_text(data = totals, aes(x = .data$group, y = 1.02, label = paste0("n=", .data$total)),
            inherit.aes = FALSE, size = 3, vjust = 0) +
  facet_wrap(~ variant_type, nrow = 1) +
  scale_fill_manual(values = fill_vals, name = cat_label) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1.08), expand = c(0, 0)) +
  labs(x = "Outlier statistic", y = "Share of sites",
       title = paste0("Consequences of sites in outlier regions (", cat_label, ")")) +
  theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major.x = element_blank())
save_plot(p_comp, paste0("outlier_csq_composition_", category))

# Enrichment: sites of one outlier statistic vs non-outlier sites, per category (2x2 Fisher; Haldane-corrected
# odds ratio with Woolf 95% CI for display so zero cells stay finite).
ref <- bg_counts %>% filter(.data$set == "non_outlier") %>%
  group_by(.data$variant_type, .data$cat) %>% summarise(n_ref = sum(.data$n), .groups = "drop")
grid <- stat_counts %>% select("group", "variant_type") %>% distinct() %>%
  inner_join(bind_rows(stat_counts %>% select("variant_type", "cat"), ref %>% select("variant_type", "cat")) %>% distinct(),
             by = "variant_type", relationship = "many-to-many")
enr <- grid %>%
  left_join(stat_counts, by = c("group", "variant_type", "cat")) %>%
  left_join(ref, by = c("variant_type", "cat")) %>%
  mutate(n = coalesce(.data$n, 0L), n_ref = coalesce(.data$n_ref, 0)) %>%
  group_by(.data$group, .data$variant_type) %>%
  mutate(group_total = sum(.data$n)) %>%
  ungroup() %>%
  left_join(ref %>% group_by(.data$variant_type) %>% summarise(ref_total = sum(.data$n_ref), .groups = "drop"),
            by = "variant_type") %>%
  mutate(ref_total = coalesce(.data$ref_total, 0)) %>%
  filter(.data$group_total >= opts$`min-sites`, .data$ref_total > 0)
if (nrow(enr) > 0) {
  a <- enr$n; b <- enr$group_total - enr$n; c <- enr$n_ref; d <- enr$ref_total - enr$n_ref
  enr$p_value <- mapply(function(a, b, c, d) fisher.test(matrix(c(a, b, c, d), 2))$p.value, a, b, c, d)
  log_or <- log((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)))
  se <- sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5))
  enr <- enr %>%
    mutate(prop_outlier = .data$n / .data$group_total, prop_non_outlier = .data$n_ref / .data$ref_total,
           odds_ratio = exp(log_or), ci_low = exp(log_or - 1.96 * se), ci_high = exp(log_or + 1.96 * se)) %>%
    group_by(.data$variant_type) %>%
    mutate(q_value = p.adjust(.data$p_value, method = "BH")) %>%
    ungroup() %>%
    rename(outlier_stat = "group", !!category := "cat", n_sites = "n", n_sites_stat = "group_total",
           n_non_outlier = "n_ref", n_non_outlier_total = "ref_total") %>%
    arrange(.data$variant_type, .data$outlier_stat, csq_factor(.data[[category]], category))
  enr_file <- file.path(opts$`output-dir`, paste0(opts$`file-prefix`, "outlier_csq_enrichment_", category, ".csv"))
  write_csv(enr, enr_file)
  message("Wrote ", enr_file)

  pe <- enr %>% mutate(cat = csq_factor(.data[[category]], category), sig = factor(.data$q_value < opts$alpha, levels = c(FALSE, TRUE)))
  dodge <- position_dodge(width = 0.6)
  p_enr <- ggplot(pe, aes(x = .data$cat, y = log2(.data$odds_ratio), colour = .data$outlier_stat)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = log2(.data$ci_low), ymax = log2(.data$ci_high)), width = 0.2, position = dodge) +
    geom_point(aes(shape = .data$sig), size = 2.4, position = dodge, fill = "white") +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 16), name = paste0("BH q < ", opts$alpha), drop = FALSE) +
    scale_colour_manual(values = rep(PLOT_PALETTE_QUALITATIVE, length.out = length(stat_levels)), name = "Outlier statistic") +
    facet_wrap(~ variant_type, nrow = 1, scales = "free_x") +
    labs(x = cat_label, y = "log2 odds ratio (outlier vs non-outlier sites)",
         title = "Consequence enrichment in outlier regions") +
    theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p_enr, paste0("outlier_csq_enrichment_", category))
} else {
  message("No statistic group has >= ", opts$`min-sites`, " sites; enrichment skipped")
}
message("Done.")
