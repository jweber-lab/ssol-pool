###############################################################################
# plot_common.R
#
# Shared helpers for popgen plotting scripts: HDF5 reading, theme constants,
# and color palettes. Sourced by plot_diversity.R, plot_fst.R, plot_pbe.R,
# and plot_region.R.
#
# Requires: hdf5r (when using read_h5_windows)
###############################################################################

# Default theme and DPI (used by Section 3 figure standards)
PLOT_BASE_SIZE <- 11
PLOT_DPI <- 300

# Colorblind-friendly qualitative palette (Okabe-Ito inspired; 8 colors)
# Use for sample / pair / trio in overlay and plot_region
PLOT_PALETTE_QUALITATIVE <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

# Consequence categories from variant_csq.R (coding_effect, region_category), most to least severe;
# "none" = no BCSQ record (intergenic or not annotated). Unknown values are appended.
CSQ_EFFECT_LEVELS <- c("nonsense", "frameshift", "missense", "inframe_indel", "other_coding",
                       "synonymous", "not_coding", "none")
CSQ_REGION_LEVELS <- c("cds", "splice", "five_prime_utr", "three_prime_utr", "ncRNA", "intron",
                       "upstream", "downstream", "intergenic", "none")
CSQ_EFFECT_COLORS <- c(nonsense = "#7F0000", frameshift = "#D55E00", missense = "#E69F00",
                       inframe_indel = "#CC79A7", other_coding = "#F0E442", synonymous = "#56B4E9",
                       not_coding = "#009E73", none = "#BBBBBB")

#' Factor with consequence levels in severity order (unknown values last).
#' @noRd
csq_factor <- function(x, category = c("coding_effect", "region_category")) {
  category <- match.arg(category)
  base <- if (category == "coding_effect") CSQ_EFFECT_LEVELS else CSQ_REGION_LEVELS
  x[is.na(x) | !nzchar(x)] <- "none"
  factor(x, levels = c(base[base %in% x], sort(setdiff(unique(x), base))))
}

#' Fill scale for consequence categories (fixed colors for coding_effect, palette otherwise).
#' @noRd
csq_fill_values <- function(levels) {
  cols <- CSQ_EFFECT_COLORS[levels]
  miss <- is.na(cols)
  if (any(miss)) {
    pal <- grDevices::hcl.colors(max(3L, sum(miss)), "Dark 3")
    cols[miss] <- pal[seq_len(sum(miss))]
  }
  stats::setNames(unname(cols), levels)
}

#' Expand plot-format option ("png,pdf", "both", "all") into file extensions.
#' @noRd
plot_formats <- function(fmt) {
  f <- trimws(tolower(strsplit(fmt, ",", fixed = TRUE)[[1]]))
  if ("both" %in% f) f <- c(setdiff(f, "both"), "png", "pdf")
  if ("all" %in% f) f <- c(setdiff(f, "all"), "png", "pdf", "svg")
  unique(f)
}

#' Read an HDF5 group into a tibble (one column per dataset).
#'
#' @param path Path to .h5 file.
#' @param group_name Group to read (e.g. "windows" or "sites").
#' @return Tibble with same column names as group datasets; character columns
#'   for chr/sample/pop1/pop2/pop3 and string datasets (genes, ref, alt, bcsq,
#'   coding_effect, <sample>_AD, ...), numeric for rest. Returns NULL on error.
#' @noRd
read_h5_windows <- function(path, group_name = "windows") {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required for HDF5 input. Install with: install.packages(\"hdf5r\")")
  }
  if (!file.exists(path)) return(NULL)
  h5 <- tryCatch(hdf5r::H5File$new(path, mode = "r"), error = function(e) NULL)
  if (is.null(h5)) return(NULL)
  on.exit(h5$close_all(), add = TRUE)
  if (!group_name %in% names(h5)) return(NULL)
  grp <- h5[[group_name]]
  nms <- names(grp)
  if (length(nms) == 0) return(NULL)
  n <- length(grp[[nms[1L]]][])
  char_cols <- c("chr", "chromosome", "sample", "pop1", "pop2", "pop3")
  out <- lapply(nms, function(nm) {
    v <- grp[[nm]][]
    if (nm %in% char_cols || is.character(v)) as.character(v) else as.numeric(v)
  })
  names(out) <- nms
  tibble::as_tibble(out, .name_repair = "minimal")
}

#' Drop per-site collate files (fst_single.h5, pbe_single.h5) when windowed files are present,
#' so window plots never mix in or default to single-position values.
#' @noRd
prefer_window_h5 <- function(files) {
  windowed <- files[!grepl("_single\\.h5$", basename(files))]
  if (length(windowed) > 0) windowed else files
}

#' Parse window_size and step_size from HDF5 or TSV filename (e.g. diversity_w10000_s5000.h5).
#' @return list with window_size, step_size (numeric), or NA if not matched.
#' @noRd
parse_window_step_from_filename <- function(basename_file) {
  m <- regmatches(basename_file, regexpr("w(\\d+)_s(\\d+)", basename_file, ignore.case = TRUE))
  if (length(m) == 0) return(list(window_size = NA_real_, step_size = NA_real_))
  parts <- strsplit(m, "_")[[1]]
  list(
    window_size = as.numeric(sub("w", "", parts[1], ignore.case = TRUE)),
    step_size = as.numeric(sub("s", "", parts[2], ignore.case = TRUE))
  )
}

#' Read companion *_summary.tsv for an HDF5 file (same dir, same base name).
#' @param h5_path Full path to .h5 file.
#' collate writes percentile columns as q001/q005/q095/q099 (percent, zero-padded);
#' they are renamed to q1/q5/q95/q99 to match the plotting scripts' summaries.
#' @return Tibble or NULL if file not found / read error.
#' @noRd
read_summary_tsv_for_h5 <- function(h5_path) {
  summary_path <- sub("\\.h5$", "_summary.tsv", h5_path)
  if (!file.exists(summary_path)) return(NULL)
  d <- tryCatch(
    readr::read_tsv(summary_path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(d)) return(NULL)
  qcols <- grep("^q0[0-9]+$", names(d), value = TRUE)
  names(d)[match(qcols, names(d))] <- paste0("q", as.integer(sub("^q", "", qcols)))
  d
}

#' Read collated per-site variant tables from a collate output directory.
#'
#' Uses single_position.h5 (--single-position-merged; FST/PBE already wide on SNP rows)
#' when present, otherwise variants.h5 with fst_single.h5 pivoted to wide
#' (<pop1>:<pop2>.fst, .fst_quantile) and joined onto SNP rows by chr, pos.
#' @return Tibble with chr, pos, variant_type, consequence/gene columns and FST columns, or NULL.
#' @noRd
read_h5_site_table <- function(hdf5_dir) {
  sp <- file.path(hdf5_dir, "single_position.h5")
  if (file.exists(sp)) return(read_h5_windows(sp, "sites"))
  vp <- file.path(hdf5_dir, "variants.h5")
  if (!file.exists(vp)) return(NULL)
  v <- read_h5_windows(vp, "sites")
  if (is.null(v)) return(NULL)
  if (!"variant_type" %in% names(v)) v$variant_type <- "snp"
  fp <- file.path(hdf5_dir, "fst_single.h5")
  f <- if (file.exists(fp)) read_h5_windows(fp, "sites") else NULL
  if (is.null(f) || !all(c("pop1", "pop2", "fst") %in% names(f))) return(v)
  f <- f %>%
    dplyr::mutate(pair = paste(.data$pop1, .data$pop2, sep = ":")) %>%
    dplyr::select("chr", "pos", "pair", "fst", dplyr::any_of("fst_quantile")) %>%
    dplyr::distinct(.data$chr, .data$pos, .data$pair, .keep_all = TRUE)
  value_cols <- intersect(c("fst", "fst_quantile"), names(f))
  fw <- tidyr::pivot_wider(f, names_from = "pair", values_from = dplyr::all_of(value_cols),
                           names_glue = "{pair}.{.value}")
  snp <- v$variant_type == "snp"
  out <- dplyr::bind_rows(
    dplyr::left_join(v[snp, , drop = FALSE], fw, by = c("chr", "pos")),
    v[!snp, , drop = FALSE]
  )
  out %>% dplyr::arrange(.data$chr, .data$pos)
}
