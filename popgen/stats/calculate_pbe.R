#!/usr/bin/env Rscript

###############################################################################
# calculate_pbe.R
#
# Calculate Population Branch Statistic (PBS) and Population Branch Excess
# (PBE) from FST TSV/CSV output from grenedalf/calculate_fst.sh.
#
# PBS1 = (T12 + T13 - T23) / 2, where Tij = -ln(1 - Fst_ij)
# PBE1 = PBS1 - (T23 * median(PBS1)) / median(T23)
# where medians are whole-genome or chromosome-specific (--median-scope).
#
# Structure: input, calculation, and output are separated.
# Output format: csv, tsv, or hdf5 (requires hdf5r when using hdf5).
# Option --drop-all-na: drop rows where any of fst_12/fst_13/fst_23 is NA; skip output if no rows remain.
###############################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(optparse)
    library(readr)
    library(tibble)
})

# -----------------------------------------------------------------------------
# INPUT
# -----------------------------------------------------------------------------

#' Detect delimiter of a tabular file from first line
detect_delimiter <- function(path) {
    first <- readLines(path, n = 1L)
    if (grepl("\t", first) && !grepl(",", first)) return("\t")
    if (grepl(",", first)) return(",")
    "\t"
}

#' Find FST column names for three population pairs.
#' Grenedalf uses columns ending in .fst or _fst containing both sample names.
find_fst_columns <- function(header_names, pop1, pop2, pop3) {
    is_fst_col <- function(name, a, b) {
        n <- tolower(name)
        if (!grepl("(\\.fst|_fst)$", n, ignore.case = TRUE)) return(FALSE)
        (grepl(a, name, fixed = TRUE) && grepl(b, name, fixed = TRUE))
    }
    col_12 <- col_13 <- col_23 <- NULL
    for (i in seq_along(header_names)) {
        nm <- header_names[i]
        if (is.null(col_12) && is_fst_col(nm, pop1, pop2)) col_12 <- nm
        if (is.null(col_13) && is_fst_col(nm, pop1, pop3)) col_13 <- nm
        if (is.null(col_23) && is_fst_col(nm, pop2, pop3)) col_23 <- nm
    }
    list(col_12 = col_12, col_13 = col_13, col_23 = col_23)
}

#' Read FST file and return a tibble with chr, pos, fst_12, fst_13, fst_23.
#' Handles csv/tsv and standard NA/nan string representations.
read_fst_for_pbs <- function(path, pop1, pop2, pop3) {
    delim <- detect_delimiter(path)
    na_vec <- c("", "NA", "na", "nan", "NaN", "NAN", "nan", "NA")
    if (delim == "\t") {
        raw <- readr::read_tsv(path, na = na_vec, show_col_types = FALSE)
    } else {
        raw <- readr::read_csv(path, na = na_vec, show_col_types = FALSE)
    }
    header_names <- colnames(raw)
    cols <- find_fst_columns(header_names, pop1, pop2, pop3)
    if (is.null(cols$col_12) || is.null(cols$col_13) || is.null(cols$col_23)) {
        stop("Could not find FST columns for all three pairs. Header: ",
             paste(header_names, collapse = " | "))
    }
    chr_col <- names(raw)[1L]
    pos_col <- names(raw)[2L]
    if (is.null(chr_col) || is.null(pos_col)) {
        chr_col <- header_names[1L]
        pos_col <- header_names[2L]
    }
    tibble(
        chr = raw[[chr_col]],
        pos = as.numeric(raw[[pos_col]]),
        fst_12 = as.numeric(raw[[cols$col_12]]),
        fst_13 = as.numeric(raw[[cols$col_13]]),
        fst_23 = as.numeric(raw[[cols$col_23]])
    )
}

# -----------------------------------------------------------------------------
# CALCULATION (vectorized)
# -----------------------------------------------------------------------------

#' Branch length Tij = -log(1 - Fst_ij). Fst in [0,1); invalid -> NA.
fst_to_branch <- function(fst) {
    fst <- as.numeric(fst)
    out <- rep(NA_real_, length(fst))
    ok <- is.finite(fst) & fst >= 0 & fst < 1
    out[ok & fst == 0] <- 0
    ok_pos <- ok & fst > 0
    out[ok_pos] <- -log(1 - fst[ok_pos])
    out[!ok & is.finite(fst) & fst >= 1] <- 10  # cap as in original
    out
}

#' PBS = (T12 + T13 - T23) / 2, vectorized.
compute_pbs <- function(fst_12, fst_13, fst_23) {
    T12 <- fst_to_branch(fst_12)
    T13 <- fst_to_branch(fst_13)
    T23 <- fst_to_branch(fst_23)
    (T12 + T13 - T23) / 2
}

#' PBE1 = PBS1 - (T23 * median(PBS1)) / median(T23)
#' med_PBS and med_T23 must be vectors of length length(pbs) (genome: repeated scalar; chromosome: per-chr medians).
compute_pbe1 <- function(pbs, T23, med_PBS, med_T23) {
    out <- rep(NA_real_, length(pbs))
    ok <- is.finite(pbs) & is.finite(T23) & is.finite(med_PBS) & is.finite(med_T23) & (med_T23 > 0)
    out[ok] <- pbs[ok] - (T23[ok] * med_PBS[ok]) / med_T23[ok]
    out
}

# -----------------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------------

#' Write results as CSV or TSV
write_table_format <- function(dat, path, format = c("tsv", "csv")) {
    format <- match.arg(format)
    if (format == "tsv") {
        readr::write_tsv(dat, path)
    } else {
        readr::write_csv(dat, path)
    }
}

#' Write results as HDF5 (requires hdf5r). Stores one group "pbs_pbe" with
#' dataset "table" as a compound layout or as separate vectors for simplicity.
write_hdf5 <- function(dat, path) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Output format 'hdf5' requires package 'hdf5r'. Install with: install.packages(\"hdf5r\")")
    }
    h5 <- hdf5r::H5File$new(path, mode = "w")
    on.exit(h5$close_all(), add = TRUE)
    # Write as CSV string inside HDF5 for simplicity (keeps columns together)
    # Alternatively write numeric matrices; we do both for convenience.
    grp <- h5$create_group("pbs_pbe")
    grp[["chr"]] <- as.character(dat$chr)
    grp[["pos"]] <- dat$pos
    grp[["fst_12"]] <- dat$fst_12
    grp[["fst_13"]] <- dat$fst_13
    grp[["fst_23"]] <- dat$fst_23
    grp[["T12"]] <- dat$T12
    grp[["T13"]] <- dat$T13
    grp[["T23"]] <- dat$T23
    grp[["PBS"]] <- dat$PBS
    grp[["PBE"]] <- dat$PBE
    # Also write one table as CSV blob for recovery
    tmp <- tempfile(fileext = ".csv")
    readr::write_csv(dat, tmp)
    tbl_txt <- paste(readLines(tmp), collapse = "\n")
    grp[["table_csv"]] <- tbl_txt
    unlink(tmp)
    invisible()
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

option_list <- list(
    make_option(c("--fst-file"), type = "character", default = NULL,
                help = "FST TSV/CSV file (from grenedalf/calculate_fst.sh)", metavar = "FILE"),
    make_option(c("--pop1-name"), type = "character", default = NULL,
                help = "Population 1 (target) name", metavar = "NAME"),
    make_option(c("--pop2-name"), type = "character", default = NULL,
                help = "Population 2 name", metavar = "NAME"),
    make_option(c("--pop3-name"), type = "character", default = NULL,
                help = "Population 3 name", metavar = "NAME"),
    make_option(c("--output-dir"), type = "character", default = ".",
                help = "Output directory [default: %default]", metavar = "DIR"),
    make_option(c("--output-format"), type = "character", default = "tsv",
                help = "Output format: tsv, csv, or hdf5 [default: %default]", metavar = "FMT"),
    make_option(c("--file-prefix"), type = "character", default = "",
                help = "Optional prefix for output files", metavar = "PREFIX"),
    make_option(c("--median-scope"), type = "character", default = "genome",
                help = "Medians for PBE: genome or chromosome [default: %default]", metavar = "SCOPE"),
    make_option(c("--drop-all-na"), action = "store_true", default = FALSE,
                help = "Drop rows where any of fst_12/fst_13/fst_23 is NA (keep only rows with all three non-NA); skip output if no rows remain"),
    make_option(c("--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

parser <- OptionParser(usage = "usage: %prog --fst-file FILE --pop1-name A --pop2-name B --pop3-name C [options]",
                       option_list = option_list)
opts <- parse_args(parser)

if (is.null(opts$`fst-file`) || is.null(opts$`pop1-name`) ||
    is.null(opts$`pop2-name`) || is.null(opts$`pop3-name`)) {
    stop("Must specify --fst-file, --pop1-name, --pop2-name, --pop3-name")
}

fmt <- tolower(opts$`output-format`)
if (!fmt %in% c("tsv", "csv", "hdf5")) {
    stop("--output-format must be one of: tsv, csv, hdf5")
}

median_scope <- tolower(opts$`median-scope`)
if (!median_scope %in% c("genome", "chromosome")) {
    stop("--median-scope must be one of: genome, chromosome")
}

if (!dir.exists(opts$`output-dir`)) {
    dir.create(opts$`output-dir`, recursive = TRUE)
}

# ---- INPUT ----
if (opts$verbose) cat("Reading FST file:", opts$`fst-file`, "\n")
dat <- read_fst_for_pbs(opts$`fst-file`, opts$`pop1-name`, opts$`pop2-name`, opts$`pop3-name`)

if (opts$`drop-all-na`) {
    n_before <- nrow(dat)
    dat <- dat %>% filter(!is.na(fst_12) & !is.na(fst_13) & !is.na(fst_23))
    n_dropped <- n_before - nrow(dat)
    message("Dropped ", n_dropped, " rows with NA in any of fst_12/fst_13/fst_23")
}

if (nrow(dat) == 0) {
    warning("No rows remaining after dropping NA; skipping output.")
} else {
# ---- CALCULATION ----
dat$T12 <- fst_to_branch(dat$fst_12)
dat$T13 <- fst_to_branch(dat$fst_13)
dat$T23 <- fst_to_branch(dat$fst_23)
dat$PBS <- compute_pbs(dat$fst_12, dat$fst_13, dat$fst_23)

# PBE1 = PBS1 - (T23 * median(PBS1)) / median(T23); medians genome-wide or per-chromosome
if (median_scope == "genome") {
    med_PBS <- rep(median(dat$PBS, na.rm = TRUE), nrow(dat))
    med_T23 <- rep(median(dat$T23, na.rm = TRUE), nrow(dat))
    dat$PBE <- compute_pbe1(dat$PBS, dat$T23, med_PBS, med_T23)
} else {
    by_chr <- dat %>% group_by(chr) %>% summarise(med_PBS = median(PBS, na.rm = TRUE), med_T23 = median(T23, na.rm = TRUE), .groups = "drop")
    dat <- dat %>% left_join(by_chr, by = "chr")
    dat$PBE <- compute_pbe1(dat$PBS, dat$T23, dat$med_PBS, dat$med_T23)
    dat$med_PBS <- NULL
    dat$med_T23 <- NULL
}

# ---- OUTPUT ----
prefix <- if (nzchar(opts$`file-prefix`)) paste0(opts$`file-prefix`, "_") else ""

if (fmt == "hdf5") {
    out_path <- file.path(opts$`output-dir`, paste0(prefix, "pbs_pbe.h5"))
    write_hdf5(dat, out_path)
    cat("Written HDF5:", out_path, "\n")
} else {
    # Use display names for CSV/TSV
    out_names <- as.data.frame(dat)
    names(out_names)[names(out_names) == "fst_12"] <- paste0("Fst_", opts$`pop1-name`, "_", opts$`pop2-name`)
    names(out_names)[names(out_names) == "fst_13"] <- paste0("Fst_", opts$`pop1-name`, "_", opts$`pop3-name`)
    names(out_names)[names(out_names) == "fst_23"] <- paste0("Fst_", opts$`pop2-name`, "_", opts$`pop3-name`)
    pbs_path <- file.path(opts$`output-dir`, paste0(prefix, "pbs.", fmt))
    pbe_path <- file.path(opts$`output-dir`, paste0(prefix, "pbe.", fmt))
    pbs_cols <- c("chr", "pos", paste0("Fst_", opts$`pop1-name`, "_", opts$`pop2-name`),
                  paste0("Fst_", opts$`pop1-name`, "_", opts$`pop3-name`),
                  paste0("Fst_", opts$`pop2-name`, "_", opts$`pop3-name`),
                  "T12", "T13", "T23", "PBS")
    pbe_cols <- c(pbs_cols, "PBE")
    write_table_format(out_names[, pbs_cols], pbs_path, format = fmt)
    write_table_format(out_names[, pbe_cols], pbe_path, format = fmt)
    cat("Written PBS:", pbs_path, "\n")
    cat("Written PBE:", pbe_path, "\n")
}
}
