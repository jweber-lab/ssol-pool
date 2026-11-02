#!/usr/bin/env Rscript

###############################################################################
# calculate_lrpi.R
#
# Log ratio of nucleotide diversity between two populations:
#   LRPI = ln(pi_pop1 / pi_pop2)
# from per-sample diversity tables written by calculate_pi_theta.sh.
#
# pi must be finite and > 0 in both populations; otherwise LRPI is NA.
# Swapping the two populations negates the value.
#
# The shell calculate_lrpi.sh discovers files and pairs. This script joins
# one pair of tables on the shared window grid (chr, start, end).
###############################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(optparse)
    library(readr)
    library(tibble)
})

detect_delimiter <- function(path) {
    first <- readLines(path, n = 1L)
    if (grepl("\t", first) && !grepl(",", first)) return("\t")
    if (grepl(",", first)) return(",")
    "\t"
}

first_named <- function(d, candidates) {
    lo <- tolower(names(d))
    hit <- match(candidates, lo)
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0L) return(NA_character_)
    names(d)[hit[1L]]
}

#' Grenedalf writes sample.theta_pi. Also accept theta_pi or pi.
find_pi_column <- function(header_names, sample) {
    lo <- tolower(header_names)
    prefer <- which(lo == tolower(paste0(sample, ".theta_pi")))
    if (length(prefer) == 1L) return(header_names[prefer])
    dotted <- grep("\\.theta_pi$", header_names)
    if (length(dotted) == 1L) return(header_names[dotted])
    exact <- which(lo %in% c("theta_pi", "pi"))
    if (length(exact) >= 1L) return(header_names[exact[1L]])
    NULL
}

read_pi_table <- function(path, sample) {
    delim <- detect_delimiter(path)
    na_vec <- c("", "NA", "na", "nan", "NaN", "NAN", ".")
    raw <- if (delim == "\t") {
        readr::read_tsv(path, na = na_vec, show_col_types = FALSE)
    } else {
        readr::read_csv(path, na = na_vec, show_col_types = FALSE)
    }
    pi_col <- find_pi_column(colnames(raw), sample)
    if (is.null(pi_col)) {
        stop("No pi column in ", path, " (expected ", sample, ".theta_pi, or theta_pi). Header: ",
             paste(colnames(raw), collapse = " | "))
    }
    chr_col <- first_named(raw, c("chr", "chromosome", "chrom"))
    start_col <- first_named(raw, c("start", "begin"))
    end_col <- first_named(raw, c("end", "stop"))
    pos_col <- first_named(raw, c("position", "pos"))
    if (is.na(pos_col)) pos_col <- start_col
    if (is.na(chr_col) || is.na(pos_col)) {
        stop("Diversity file needs a chromosome column and pos or start: ", path)
    }
    out <- tibble(
        chr = as.character(raw[[chr_col]]),
        pos = as.numeric(raw[[pos_col]]),
        pi = as.numeric(raw[[pi_col]])
    )
    if (!is.na(start_col) && !is.na(end_col)) {
        out$start <- as.numeric(raw[[start_col]])
        out$end <- as.numeric(raw[[end_col]])
    }
    out
}

#' ln(pi_a / pi_b). Non-positive or non-finite pi is NA.
compute_lrpi <- function(pi_a, pi_b) {
    out <- rep(NA_real_, length(pi_a))
    ok <- is.finite(pi_a) & is.finite(pi_b) & pi_a > 0 & pi_b > 0
    out[ok] <- log(pi_a[ok] / pi_b[ok])
    out
}

write_table_format <- function(dat, path, format = c("tsv", "csv")) {
    format <- match.arg(format)
    if (format == "tsv") readr::write_tsv(dat, path) else readr::write_csv(dat, path)
}

write_hdf5 <- function(dat, path) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Output format 'hdf5' requires package 'hdf5r'. Install with: install.packages(\"hdf5r\")")
    }
    h5 <- hdf5r::H5File$new(path, mode = "w")
    on.exit(h5$close_all(), add = TRUE)
    grp <- h5$create_group("lrpi")
    grp[["chr"]] <- as.character(dat$chr)
    grp[["pos"]] <- dat$pos
    if ("start" %in% names(dat)) grp[["start"]] <- dat$start
    if ("end" %in% names(dat)) grp[["end"]] <- dat$end
    grp[["pop1"]] <- as.character(dat$pop1)
    grp[["pop2"]] <- as.character(dat$pop2)
    grp[["pi_pop1"]] <- dat$pi_pop1
    grp[["pi_pop2"]] <- dat$pi_pop2
    grp[["lrpi"]] <- dat$lrpi
    invisible()
}

option_list <- list(
    make_option(c("--diversity-a"), type = "character", default = NULL,
                help = "Diversity TSV/CSV for population 1 (numerator)", metavar = "FILE"),
    make_option(c("--diversity-b"), type = "character", default = NULL,
                help = "Diversity TSV/CSV for population 2 (denominator)", metavar = "FILE"),
    make_option(c("--pop1-name"), type = "character", default = NULL,
                help = "Population 1 name (numerator of pi_A/pi_B)", metavar = "NAME"),
    make_option(c("--pop2-name"), type = "character", default = NULL,
                help = "Population 2 name (denominator)", metavar = "NAME"),
    make_option(c("--output-dir"), type = "character", default = ".",
                help = "Output directory [default: %default]", metavar = "DIR"),
    make_option(c("--output-format"), type = "character", default = "tsv",
                help = "Output format: tsv, csv, or hdf5 [default: %default]", metavar = "FMT"),
    make_option(c("--file-prefix"), type = "character", default = "",
                help = "Optional prefix for the output file", metavar = "PREFIX"),
    make_option(c("--drop-all-na"), action = "store_true", default = FALSE,
                help = "Drop rows where LRPI is NA"),
    make_option(c("--verbose"), action = "store_true", default = FALSE,
                help = "Print row counts and how many pi values are missing or non-positive")
)

parser <- OptionParser(
    usage = "usage: %prog --diversity-a FILE --diversity-b FILE --pop1-name A --pop2-name B [options]",
    option_list = option_list
)
opts <- parse_args(parser)

if (is.null(opts$`diversity-a`) || is.null(opts$`diversity-b`) ||
    is.null(opts$`pop1-name`) || is.null(opts$`pop2-name`)) {
    stop("Must specify --diversity-a, --diversity-b, --pop1-name, and --pop2-name")
}
if (identical(opts$`pop1-name`, opts$`pop2-name`)) {
    stop("pop1 and pop2 must be different populations")
}
fmt <- tolower(opts$`output-format`)
if (!fmt %in% c("tsv", "csv", "hdf5")) stop("--output-format must be one of: tsv, csv, hdf5")
if (!dir.exists(opts$`output-dir`)) dir.create(opts$`output-dir`, recursive = TRUE)

a <- read_pi_table(opts$`diversity-a`, opts$`pop1-name`)
b <- read_pi_table(opts$`diversity-b`, opts$`pop2-name`)
keys <- if (all(c("start", "end") %in% names(a)) && all(c("start", "end") %in% names(b))) {
    c("chr", "start", "end")
} else {
    c("chr", "pos")
}
if (anyDuplicated(a %>% select(all_of(keys))) > 0L) stop("Duplicate windows in ", opts$`diversity-a`, " on ", paste(keys, collapse = ", "))
if (anyDuplicated(b %>% select(all_of(keys))) > 0L) stop("Duplicate windows in ", opts$`diversity-b`, " on ", paste(keys, collapse = ", "))

a_only <- a %>% anti_join(b, by = keys)
b_only <- b %>% anti_join(a, by = keys)
dat <- a %>%
    rename(pi_pop1 = "pi") %>%
    inner_join(b %>% select(all_of(c(keys, "pi"))) %>% rename(pi_pop2 = "pi"), by = keys)
dat$lrpi <- compute_lrpi(dat$pi_pop1, dat$pi_pop2)
dat$pop1 <- opts$`pop1-name`
dat$pop2 <- opts$`pop2-name`

if (opts$verbose) {
    message("Joined on ", paste(keys, collapse = ", "),
            ": rows=", nrow(dat),
            " only_in_", opts$`pop1-name`, "=", nrow(a_only),
            " only_in_", opts$`pop2-name`, "=", nrow(b_only))
    message("LRPI NA=", sum(is.na(dat$lrpi)),
            " (pi_pop1 non-positive=", sum(is.finite(dat$pi_pop1) & dat$pi_pop1 <= 0),
            ", pi_pop2 non-positive=", sum(is.finite(dat$pi_pop2) & dat$pi_pop2 <= 0), ")")
}

if (opts$`drop-all-na`) {
    n_before <- nrow(dat)
    dat <- dat %>% filter(!is.na(.data$lrpi))
    message("Dropped ", n_before - nrow(dat), " rows with NA LRPI")
}

if (nrow(dat) == 0L) {
    warning("No rows remaining; skipping output.")
} else {
    lead <- c("chr", "pos", intersect(c("start", "end"), names(dat)), "pop1", "pop2", "pi_pop1", "pi_pop2", "lrpi")
    dat <- dat %>% select(all_of(lead))
    prefix <- if (nzchar(opts$`file-prefix`)) paste0(opts$`file-prefix`, "_") else ""
    if (fmt == "hdf5") {
        out_path <- file.path(opts$`output-dir`, paste0(prefix, "lrpi.h5"))
        write_hdf5(dat, out_path)
        cat("Written LRPI:", out_path, "\n")
    } else {
        out_path <- file.path(opts$`output-dir`, paste0(prefix, "lrpi.", fmt))
        write_table_format(dat, out_path, format = fmt)
        cat("Written LRPI:", out_path, "\n")
    }
}
