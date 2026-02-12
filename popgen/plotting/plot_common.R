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

#' Read an HDF5 group into a tibble (one column per dataset).
#'
#' @param path Path to .h5 file.
#' @param group_name Group to read (e.g. "windows" or "sites").
#' @return Tibble with same column names as group datasets; character columns
#'   for chr/sample/pop1/pop2/pop3, numeric for rest. Returns NULL on error.
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
    if (nm %in% char_cols) as.character(v) else as.numeric(v)
  })
  names(out) <- nms
  d <- as.data.frame(out, stringsAsFactors = FALSE)
  dplyr::as_tibble(d)
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
#' @return Tibble or NULL if file not found / read error.
#' @noRd
read_summary_tsv_for_h5 <- function(h5_path) {
  summary_path <- sub("\\.h5$", "_summary.tsv", h5_path)
  if (!file.exists(summary_path)) return(NULL)
  tryCatch(
    readr::read_tsv(summary_path, show_col_types = FALSE),
    error = function(e) NULL
  )
}
