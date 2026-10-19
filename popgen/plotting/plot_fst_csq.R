#!/usr/bin/env Rscript

###############################################################################
# plot_fst_csq.R
#
# Per-site FST by consequence class, from collate output (--hdf5-dir):
#   single_position.h5 (collate --single-position-merged), or variants.h5 + fst_single.h5.
# FST exists only on SNP rows; indels are skipped. Consequence columns (coding_effect,
# region_category) come from collate --variant-csq.
#
# Optional --regions (identify_outliers outlier_regions_*.csv) splits sites into inside vs
# outside outlier regions (regions whose outlier_stat includes --region-stat).
#
# Outputs (in --output-dir, <cat> = --category):
#   fst_csq_<cat>.png            Violin + box of per-site FST per category, one panel per pair
#                                (fill = inside/outside outlier regions with --regions)
#   fst_csq_<cat>.csv            Per pair x category (x region set): n, mean, median, q95, Wilcoxon
#                                p vs --reference-category and (with --regions) inside vs outside, BH q
#   fst_csq_manhattan_<cat>.png  Per-site FST along the genome, coloured by category, outlier regions
#                                shaded (with --manhattan)
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
  make_option("--hdf5-dir", type = "character", default = NULL,
    help = "Collate output dir with single_position.h5, or variants.h5 + fst_single.h5 (required)"),
  make_option("--category", type = "character", default = "coding_effect",
    help = "coding_effect or region_category [default: %default]"),
  make_option("--y-value", type = "character", default = "value", help = "value or quantile [default: %default]"),
  make_option("--sample-pairs", type = "character", default = "", help = "Comma-separated pairs (Pop1:Pop2) [default: all]"),
  make_option("--chromosome", type = "character", default = NULL, help = "Restrict to one chromosome"),
  make_option("--regions", type = "character", default = NULL,
    help = "outlier_regions_*.csv from identify_outliers: compare sites inside vs outside outlier regions"),
  make_option("--region-stat", type = "character", default = "fst",
    help = "Use regions whose outlier_stat includes this statistic; 'any' for all regions [default: %default]"),
  make_option("--reference-category", type = "character", default = "synonymous",
    help = "Category each class is tested against (Wilcoxon) [default: %default; falls back to 'none']"),
  make_option("--min-sites", type = "integer", default = 3L, help = "Minimum sites to draw/test a category [default: %default]"),
  make_option("--drop-none", action = "store_true", default = FALSE, help = "Drop sites with no consequence record ('none')"),
  make_option("--manhattan", action = "store_true", default = FALSE, help = "Also plot per-site FST along the genome"),
  make_option(c("--output-dir", "-o"), type = "character", default = "./", help = "Output directory"),
  make_option("--file-prefix", type = "character", default = "", help = "Prefix for output file names"),
  make_option("--width", type = "double", default = 12, help = "Width in inches [default: %default]"),
  make_option("--height", type = "double", default = 7, help = "Height in inches [default: %default]"),
  make_option("--dpi", type = "integer", default = PLOT_DPI, help = "PNG DPI [default: %default]"),
  make_option("--plot-format", type = "character", default = "png", help = "png, pdf, svg, both, or all [default: %default]")
)
opts <- parse_args(OptionParser(option_list = option_list))
if (is.null(opts$`hdf5-dir`)) stop("--hdf5-dir is required")
category <- match.arg(opts$category, c("coding_effect", "region_category"))
y_value <- match.arg(opts$`y-value`, c("value", "quantile"))

sites <- read_h5_site_table(opts$`hdf5-dir`)
if (is.null(sites) || nrow(sites) == 0)
  stop("No single_position.h5 or variants.h5 in ", opts$`hdf5-dir`, " (collate --variant-tsv-dir/--variant-csq)")
if (!category %in% names(sites)) stop("No ", category, " column in collated sites; run collate with --variant-csq")
fst_cols <- grep("\\.fst$", names(sites), value = TRUE)
if (length(fst_cols) == 0) stop("No per-site FST columns (<pop1>:<pop2>.fst); collate needs fst_single input from --fst-dir")
if (!"variant_type" %in% names(sites)) sites$variant_type <- "snp"
sites <- sites %>% filter(.data$variant_type == "snp")
if (!is.null(opts$chromosome)) sites <- sites %>% filter(.data$chr == opts$chromosome)

pairs <- sub("\\.fst$", "", fst_cols)
if (nzchar(opts$`sample-pairs`)) {
  want <- trimws(strsplit(opts$`sample-pairs`, ",", fixed = TRUE)[[1]])
  flip <- vapply(strsplit(want, ":", fixed = TRUE), function(p) paste(rev(p), collapse = ":"), character(1))
  keep <- pairs %in% want | pairs %in% flip
  if (!any(keep)) stop("None of --sample-pairs found; available: ", paste(pairs, collapse = ", "))
  pairs <- pairs[keep]
}
value_suffix <- if (y_value == "quantile") ".fst_quantile" else ".fst"
missing_q <- setdiff(paste0(pairs, value_suffix), names(sites))
if (length(missing_q) > 0) stop("Columns not found: ", paste(missing_q, collapse = ", "))

long <- bind_rows(lapply(pairs, function(p) {
  tibble(chr = sites$chr, pos = sites$pos, pair = p, y = as.numeric(sites[[paste0(p, value_suffix)]]),
         cat = as.character(sites[[category]]))
})) %>% filter(!is.na(.data$y))
long$cat[is.na(long$cat) | !nzchar(long$cat)] <- "none"
if (opts$`drop-none`) long <- long %>% filter(.data$cat != "none")
if (nrow(long) == 0) stop("No SNP sites with FST to plot")

regions <- NULL
if (!is.null(opts$regions)) {
  regions <- read_csv(opts$regions, show_col_types = FALSE) %>% rename_with(tolower)
  sc <- intersect(c("region_start", "start"), names(regions))[1]
  ec <- intersect(c("region_end", "end"), names(regions))[1]
  if (is.na(sc) || is.na(ec) || !"chr" %in% names(regions)) stop("--regions needs chr, region_start, region_end")
  regions <- regions %>% transmute(chr = as.character(.data$chr), start = as.numeric(.data[[sc]]), end = as.numeric(.data[[ec]]),
                                   outlier_stat = if ("outlier_stat" %in% names(regions)) .data$outlier_stat else NA_character_)
  if (opts$`region-stat` != "any") {
    regions <- regions %>% filter(vapply(strsplit(coalesce(.data$outlier_stat, ""), ":", fixed = TRUE),
                                         function(s) opts$`region-stat` %in% s, logical(1)))
  }
  message(nrow(regions), " outlier regions used (", opts$`region-stat`, ")")
  inside <- rep(FALSE, nrow(long))
  for (i in seq_len(nrow(regions))) {
    inside <- inside | (long$chr == regions$chr[i] & long$pos >= regions$start[i] & long$pos <= regions$end[i])
  }
  long$set <- factor(ifelse(inside, "outlier region", "outside"), levels = c("outside", "outlier region"))
} else {
  long$set <- factor("all sites")
}

long$cat <- csq_factor(long$cat, category)
cat_n <- long %>% count(.data$pair, .data$cat, .data$set)
plot_d <- long %>% semi_join(cat_n %>% filter(.data$n >= opts$`min-sites`), by = c("pair", "cat", "set"))

ref_cat <- opts$`reference-category`
if (!ref_cat %in% long$cat) ref_cat <- "none"
wtest <- function(x, y) {
  if (length(x) < opts$`min-sites` || length(y) < opts$`min-sites`) return(NA_real_)
  suppressWarnings(wilcox.test(x, y)$p.value)
}
stats_tab <- long %>%
  group_by(.data$pair, .data$set, .data$cat) %>%
  summarise(n = dplyr::n(), mean = mean(.data$y), median = median(.data$y),
            q95 = unname(quantile(.data$y, 0.95)), .groups = "drop")
stats_tab$p_vs_reference <- mapply(function(p, s, cc) {
  if (cc == ref_cat) return(NA_real_)
  d <- long[long$pair == p & long$set == s, ]
  wtest(d$y[d$cat == cc], d$y[d$cat == ref_cat])
}, stats_tab$pair, as.character(stats_tab$set), as.character(stats_tab$cat))
stats_tab$q_vs_reference <- p.adjust(stats_tab$p_vs_reference, method = "BH")
if (!is.null(regions)) {
  io <- long %>% group_by(.data$pair, .data$cat) %>%
    summarise(p_inside_vs_outside = wtest(.data$y[.data$set == "outlier region"], .data$y[.data$set == "outside"]),
              .groups = "drop") %>%
    mutate(q_inside_vs_outside = p.adjust(.data$p_inside_vs_outside, method = "BH"))
  stats_tab <- stats_tab %>% left_join(io, by = c("pair", "cat"))
}
stats_tab <- stats_tab %>% mutate(reference_category = ref_cat) %>% rename(!!category := "cat", region_set = "set") %>%
  arrange(.data$pair, .data$region_set, .data[[category]])

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)
save_plot <- function(p, stem, width = opts$width, height = opts$height) {
  for (fmt in plot_formats(opts$`plot-format`)) {
    out <- file.path(opts$`output-dir`, paste0(opts$`file-prefix`, stem, ".", fmt))
    if (fmt == "png") ggsave(out, p, width = width, height = height, dpi = opts$dpi) else ggsave(out, p, width = width, height = height)
    message("Saved: ", out)
  }
}
tab_file <- file.path(opts$`output-dir`, paste0(opts$`file-prefix`, "fst_csq_", category, ".csv"))
write_csv(stats_tab, tab_file)
message("Wrote ", tab_file)

cat_label <- if (category == "coding_effect") "Coding effect" else "Region"
y_label <- if (y_value == "quantile") "Per-site FST quantile" else "Per-site FST"
fill_vals <- csq_fill_values(levels(long$cat))
dodge <- position_dodge(width = 0.85, preserve = "single")
n_lab <- cat_n %>% filter(.data$n >= opts$`min-sites`)
y_floor <- min(plot_d$y, na.rm = TRUE)
y_span <- diff(range(plot_d$y, na.rm = TRUE))
if (!is.finite(y_span) || y_span == 0) y_span <- 1
if (is.null(regions)) {
  p <- ggplot(plot_d, aes(x = .data$cat, y = .data$y, fill = .data$cat)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.6) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
    geom_text(data = n_lab, aes(x = .data$cat, y = y_floor - 0.06 * y_span, label = .data$n), inherit.aes = FALSE, size = 2.6) +
    scale_fill_manual(values = fill_vals, guide = "none")
} else {
  p <- ggplot(plot_d, aes(x = .data$cat, y = .data$y, fill = .data$set)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.6, position = dodge) +
    geom_boxplot(width = 0.18, outlier.shape = NA, position = dodge) +
    geom_text(data = n_lab, aes(x = .data$cat, y = y_floor - 0.06 * y_span, label = .data$n, group = .data$set),
              inherit.aes = FALSE, size = 2.4, position = position_dodge(width = 0.85)) +
    scale_fill_manual(values = c(outside = "#BBBBBB", `outlier region` = "#D55E00"), name = NULL)
}
p <- p +
  facet_wrap(~ pair) +
  labs(x = cat_label, y = y_label, title = paste0(y_label, " by consequence (SNPs)"),
       subtitle = paste0("Numbers: sites per class; tests vs '", ref_cat, "' in ", basename(tab_file))) +
  theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
save_plot(p, paste0("fst_csq_", category))

if (opts$manhattan) {
  chr_len <- long %>% group_by(.data$chr) %>% summarise(len = max(.data$pos), .groups = "drop") %>%
    arrange(desc(.data$len)) %>% mutate(offset = cumsum(c(0, head(.data$len, -1))))
  md <- long %>% inner_join(chr_len, by = "chr") %>% mutate(x = .data$pos + .data$offset,
                                                            background = .data$cat %in% c("none", "not_coding", "intergenic", "intron"))
  stripes <- chr_len %>% mutate(xmin = .data$offset, xmax = .data$offset + .data$len,
                                shade = rep(c("a", "b"), length.out = dplyr::n()))
  pm <- ggplot() +
    geom_rect(data = stripes, aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf, fill = .data$shade),
              alpha = 0.4, show.legend = FALSE) +
    scale_fill_manual(values = c(a = "white", b = "grey90"))
  if (!is.null(regions) && nrow(regions) > 0) {
    rr <- regions %>% inner_join(chr_len, by = "chr")
    pm <- pm + geom_rect(data = rr, aes(xmin = .data$start + .data$offset, xmax = .data$end + .data$offset, ymin = -Inf, ymax = Inf),
                         fill = "#D55E00", alpha = 0.12)
  }
  pm <- pm +
    geom_point(data = md %>% filter(.data$background), aes(x = .data$x, y = .data$y), colour = "grey70", size = 0.4) +
    geom_point(data = md %>% filter(!.data$background), aes(x = .data$x, y = .data$y, colour = .data$cat), size = 1.1) +
    scale_colour_manual(values = fill_vals, name = cat_label, drop = TRUE) +
    scale_x_continuous(breaks = chr_len$offset + chr_len$len / 2,
                       labels = if (nrow(chr_len) <= 30) chr_len$chr else NULL, expand = c(0.005, 0)) +
    facet_wrap(~ pair, ncol = 1) +
    labs(x = "Chromosome", y = y_label,
         title = paste0(y_label, " along the genome; coloured sites have a coding consequence")) +
    theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(pm, paste0("fst_csq_manhattan_", category), height = max(opts$height, 2.2 * length(pairs)))
}
message("Done.")
