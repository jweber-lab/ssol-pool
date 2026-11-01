#!/usr/bin/env Rscript

###############################################################################
# calculate_qibe.R
#
# Calculate Quartet Internal Branch Statistic (QIBS) and Quartet Internal
# Branch Excess (QIBE) from FST TSV/CSV output from calculate_fst.sh.
#
# Topology ((1,2),(3,4)):
#   QIBS = (T13 + T14 + T23 + T24) / 4 - (T12 + T34) / 2
#   QEBS = T12 + T34
#   QIBE = QIBS - (QEBS / median(QEBS)) * median(QIBS)
# where Tij = -ln(1 - Fst_ij), and medians are whole-genome or per chromosome
# (--median-scope). A wrong grouping can make QIBS negative.
#
# Output format: csv, tsv, or hdf5 (requires hdf5r when using hdf5).
# Option --drop-all-na: drop rows where any of the six FST pairs is NA.
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

detect_delimiter <- function(path) {
    first <- readLines(path, n = 1L)
    if (grepl("\t", first) && !grepl(",", first)) return("\t")
    if (grepl(",", first)) return(",")
    "\t"
}

# Grenedalf columns end in .fst or _fst and contain both sample names
# (e.g. Echo:Kjer.fst). A colon, underscore, or dot pair separator keeps a
# short name from matching inside a longer one (Ann vs Anne).
is_fst_col <- function(name, a, b) {
    if (!grepl("(\\.fst|_fst)$", name, ignore.case = TRUE)) return(FALSE)
    core <- sub("(\\.fst|_fst)$", "", name, ignore.case = TRUE)
    if (grepl(":", core, fixed = TRUE)) {
        parts <- strsplit(core, ":", fixed = TRUE)[[1L]]
        return(length(parts) == 2L && setequal(parts, c(a, b)))
    }
    core == paste(a, b, sep = "_") || core == paste(b, a, sep = "_")
}

find_pair_col <- function(header_names, a, b) {
    for (nm in header_names) if (is_fst_col(nm, a, b)) return(nm)
    NULL
}

first_named <- function(d, candidates) {
    lo <- tolower(names(d))
    hit <- match(candidates, lo)
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0L) return(NA_character_)
    names(d)[hit[1L]]
}

#' Read FST and return chr, pos, optional start/end, and the six quartet pairs.
read_fst_for_qibs <- function(path, pop1, pop2, pop3, pop4) {
    delim <- detect_delimiter(path)
    na_vec <- c("", "NA", "na", "nan", "NaN", "NAN")
    if (delim == "\t") {
        raw <- readr::read_tsv(path, na = na_vec, show_col_types = FALSE)
    } else {
        raw <- readr::read_csv(path, na = na_vec, show_col_types = FALSE)
    }
    header_names <- colnames(raw)
    pair_pops <- list(
        fst_12 = c(pop1, pop2), fst_34 = c(pop3, pop4),
        fst_13 = c(pop1, pop3), fst_14 = c(pop1, pop4),
        fst_23 = c(pop2, pop3), fst_24 = c(pop2, pop4)
    )
    cols <- lapply(pair_pops, function(ab) find_pair_col(header_names, ab[1], ab[2]))
    missing <- names(cols)[vapply(cols, is.null, logical(1))]
    if (length(missing) > 0L) {
        stop("Could not find FST columns for: ", paste(missing, collapse = ", "),
             ". Header: ", paste(header_names, collapse = " | "))
    }
    chr_col <- first_named(raw, c("chr", "chromosome", "chrom"))
    start_col <- first_named(raw, c("start", "begin"))
    end_col <- first_named(raw, c("end", "stop"))
    pos_col <- first_named(raw, c("position", "pos"))
    if (is.na(pos_col)) pos_col <- start_col
    if (is.na(chr_col) || is.na(pos_col)) {
        stop("FST file needs a chromosome column and pos or start. Header: ",
             paste(header_names, collapse = " | "))
    }
    out <- tibble(
        chr = as.character(raw[[chr_col]]),
        pos = as.numeric(raw[[pos_col]]),
        fst_12 = as.numeric(raw[[cols$fst_12]]),
        fst_13 = as.numeric(raw[[cols$fst_13]]),
        fst_14 = as.numeric(raw[[cols$fst_14]]),
        fst_23 = as.numeric(raw[[cols$fst_23]]),
        fst_24 = as.numeric(raw[[cols$fst_24]]),
        fst_34 = as.numeric(raw[[cols$fst_34]])
    )
    if (!is.na(start_col) && !is.na(end_col)) {
        out$start <- as.numeric(raw[[start_col]])
        out$end <- as.numeric(raw[[end_col]])
        out <- out %>% relocate("start", "end", .after = "pos")
    }
    out
}

# -----------------------------------------------------------------------------
# CALCULATION (vectorized)
# -----------------------------------------------------------------------------

#' Tij = -log(1 - Fst_ij). Fst in [0, 1); Fst >= 1 is capped at 10, as in calculate_pbe.R.
fst_to_branch <- function(fst) {
    fst <- as.numeric(fst)
    out <- rep(NA_real_, length(fst))
    ok <- is.finite(fst) & fst >= 0 & fst < 1
    out[ok & fst == 0] <- 0
    ok_pos <- ok & fst > 0
    out[ok_pos] <- -log(1 - fst[ok_pos])
    out[!ok & is.finite(fst) & fst >= 1] <- 10
    out
}

#' QIBS = (T13 + T14 + T23 + T24) / 4 - (T12 + T34) / 2. QEBS = T12 + T34.
compute_qibs_qebs <- function(fst_12, fst_13, fst_14, fst_23, fst_24, fst_34) {
    T12 <- fst_to_branch(fst_12)
    T13 <- fst_to_branch(fst_13)
    T14 <- fst_to_branch(fst_14)
    T23 <- fst_to_branch(fst_23)
    T24 <- fst_to_branch(fst_24)
    T34 <- fst_to_branch(fst_34)
    list(
        T12 = T12, T13 = T13, T14 = T14, T23 = T23, T24 = T24, T34 = T34,
        QEBS = T12 + T34,
        QIBS = (T13 + T14 + T23 + T24) / 4 - (T12 + T34) / 2
    )
}

#' QIBE = QIBS - (QEBS * median(QIBS)) / median(QEBS).
#' med_QIBS and med_QEBS have length length(qibs).
compute_qibe <- function(qibs, qebs, med_QIBS, med_QEBS) {
    out <- rep(NA_real_, length(qibs))
    ok <- is.finite(qibs) & is.finite(qebs) & is.finite(med_QIBS) & is.finite(med_QEBS) & (med_QEBS > 0)
    out[ok] <- qibs[ok] - (qebs[ok] * med_QIBS[ok]) / med_QEBS[ok]
    out
}

# -----------------------------------------------------------------------------
# OUTPUT
# -----------------------------------------------------------------------------

write_table_format <- function(dat, path, format = c("tsv", "csv")) {
    format <- match.arg(format)
    if (format == "tsv") {
        readr::write_tsv(dat, path)
    } else {
        readr::write_csv(dat, path)
    }
}

write_hdf5 <- function(dat, path) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Output format 'hdf5' requires package 'hdf5r'. Install with: install.packages(\"hdf5r\")")
    }
    h5 <- hdf5r::H5File$new(path, mode = "w")
    on.exit(h5$close_all(), add = TRUE)
    grp <- h5$create_group("qibs_qibe")
    grp[["chr"]] <- as.character(dat$chr)
    grp[["pos"]] <- dat$pos
    if ("start" %in% names(dat)) grp[["start"]] <- dat$start
    if ("end" %in% names(dat)) grp[["end"]] <- dat$end
    for (nm in c("pop1", "pop2", "pop3", "pop4")) {
        if (nm %in% names(dat)) grp[[nm]] <- as.character(dat[[nm]])
    }
    for (nm in c("fst_12", "fst_13", "fst_14", "fst_23", "fst_24", "fst_34",
                 "T12", "T13", "T14", "T23", "T24", "T34", "QEBS", "QIBS", "QIBE")) {
        grp[[nm]] <- dat[[nm]]
    }
    tmp <- tempfile(fileext = ".csv")
    readr::write_csv(dat, tmp)
    grp[["table_csv"]] <- paste(readLines(tmp), collapse = "\n")
    unlink(tmp)
    invisible()
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

option_list <- list(
    make_option(c("--fst-file"), type = "character", default = NULL,
                help = "FST TSV/CSV file (from calculate_fst.sh)", metavar = "FILE"),
    make_option(c("--pop1-name"), type = "character", default = NULL,
                help = "Population 1 (sister of population 2)", metavar = "NAME"),
    make_option(c("--pop2-name"), type = "character", default = NULL,
                help = "Population 2 (sister of population 1)", metavar = "NAME"),
    make_option(c("--pop3-name"), type = "character", default = NULL,
                help = "Population 3 (sister of population 4)", metavar = "NAME"),
    make_option(c("--pop4-name"), type = "character", default = NULL,
                help = "Population 4 (sister of population 3)", metavar = "NAME"),
    make_option(c("--output-dir"), type = "character", default = ".",
                help = "Output directory [default: %default]", metavar = "DIR"),
    make_option(c("--output-format"), type = "character", default = "tsv",
                help = "Output format: tsv, csv, or hdf5 [default: %default]", metavar = "FMT"),
    make_option(c("--file-prefix"), type = "character", default = "",
                help = "Optional prefix for output files", metavar = "PREFIX"),
    make_option(c("--median-scope"), type = "character", default = "genome",
                help = "Medians for QIBE: genome or chromosome [default: %default]", metavar = "SCOPE"),
    make_option(c("--drop-all-na"), action = "store_true", default = FALSE,
                help = "Drop rows where any of the six FST pairs is NA; skip output if none remain"),
    make_option(c("--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

parser <- OptionParser(
    usage = "usage: %prog --fst-file FILE --pop1-name A --pop2-name B --pop3-name C --pop4-name D [options]",
    option_list = option_list)
opts <- parse_args(parser)

if (is.null(opts$`fst-file`) || is.null(opts$`pop1-name`) || is.null(opts$`pop2-name`) ||
    is.null(opts$`pop3-name`) || is.null(opts$`pop4-name`)) {
    stop("Must specify --fst-file, --pop1-name, --pop2-name, --pop3-name, --pop4-name")
}

pops <- c(opts$`pop1-name`, opts$`pop2-name`, opts$`pop3-name`, opts$`pop4-name`)
if (anyDuplicated(pops)) stop("Quartet populations must be four different names: ", paste(pops, collapse = ", "))

fmt <- tolower(opts$`output-format`)
if (!fmt %in% c("tsv", "csv", "hdf5")) stop("--output-format must be one of: tsv, csv, hdf5")

median_scope <- tolower(opts$`median-scope`)
if (!median_scope %in% c("genome", "chromosome")) stop("--median-scope must be one of: genome, chromosome")

if (!dir.exists(opts$`output-dir`)) dir.create(opts$`output-dir`, recursive = TRUE)

fst_cols <- c("fst_12", "fst_13", "fst_14", "fst_23", "fst_24", "fst_34")
dat <- read_fst_for_qibs(opts$`fst-file`, pops[1], pops[2], pops[3], pops[4])
if (opts$verbose) {
    message("FST file: ", opts$`fst-file`)
    message("Rows: ", nrow(dat), "; columns: ", paste(names(dat), collapse = ", "))
    message("NA counts: ", paste(sprintf("%s=%d", fst_cols, vapply(dat[fst_cols], function(v) sum(is.na(v)), integer(1))), collapse = " "))
}

if (opts$`drop-all-na`) {
    n_before <- nrow(dat)
    keep <- !Reduce(`|`, lapply(dat[fst_cols], is.na))
    dat <- dat[keep, , drop = FALSE]
    message("Dropped ", n_before - nrow(dat), " rows with NA in any of the six FST pairs")
}

if (nrow(dat) == 0) {
    warning("No rows remaining after dropping NA; skipping output.")
} else {
    branches <- compute_qibs_qebs(dat$fst_12, dat$fst_13, dat$fst_14, dat$fst_23, dat$fst_24, dat$fst_34)
    dat$T12 <- branches$T12
    dat$T13 <- branches$T13
    dat$T14 <- branches$T14
    dat$T23 <- branches$T23
    dat$T24 <- branches$T24
    dat$T34 <- branches$T34
    dat$QEBS <- branches$QEBS
    dat$QIBS <- branches$QIBS
    dat$pop1 <- pops[1]
    dat$pop2 <- pops[2]
    dat$pop3 <- pops[3]
    dat$pop4 <- pops[4]

    if (median_scope == "genome") {
        med_QIBS <- rep(median(dat$QIBS, na.rm = TRUE), nrow(dat))
        med_QEBS <- rep(median(dat$QEBS, na.rm = TRUE), nrow(dat))
        if (opts$verbose) message("Genome medians: QIBS=", signif(med_QIBS[1], 6), " QEBS=", signif(med_QEBS[1], 6))
        dat$QIBE <- compute_qibe(dat$QIBS, dat$QEBS, med_QIBS, med_QEBS)
    } else {
        if (opts$verbose) message("Medians computed per chromosome")
        by_chr <- dat %>% group_by(.data$chr) %>%
            summarise(med_QIBS = median(.data$QIBS, na.rm = TRUE),
                      med_QEBS = median(.data$QEBS, na.rm = TRUE), .groups = "drop")
        dat <- dat %>% left_join(by_chr, by = "chr")
        dat$QIBE <- compute_qibe(dat$QIBS, dat$QEBS, dat$med_QIBS, dat$med_QEBS)
        dat$med_QIBS <- NULL
        dat$med_QEBS <- NULL
    }

    prefix <- if (nzchar(opts$`file-prefix`)) paste0(opts$`file-prefix`, "_") else ""
    if (fmt == "hdf5") {
        out_path <- file.path(opts$`output-dir`, paste0(prefix, "qibs_qibe.h5"))
        write_hdf5(dat, out_path)
        cat("Written HDF5:", out_path, "\n")
    } else {
        out_names <- as.data.frame(dat)
        pair_labels <- c(
            fst_12 = paste0("Fst_", pops[1], "_", pops[2]),
            fst_13 = paste0("Fst_", pops[1], "_", pops[3]),
            fst_14 = paste0("Fst_", pops[1], "_", pops[4]),
            fst_23 = paste0("Fst_", pops[2], "_", pops[3]),
            fst_24 = paste0("Fst_", pops[2], "_", pops[4]),
            fst_34 = paste0("Fst_", pops[3], "_", pops[4])
        )
        for (nm in names(pair_labels)) names(out_names)[names(out_names) == nm] <- pair_labels[[nm]]
        lead <- c("chr", "pos", intersect(c("start", "end"), names(out_names)), "pop1", "pop2", "pop3", "pop4")
        qibs_cols <- c(lead, unname(pair_labels), "T12", "T13", "T14", "T23", "T24", "T34", "QEBS", "QIBS")
        qibe_cols <- c(qibs_cols, "QIBE")
        qibs_path <- file.path(opts$`output-dir`, paste0(prefix, "qibs.", fmt))
        qibe_path <- file.path(opts$`output-dir`, paste0(prefix, "qibe.", fmt))
        write_table_format(out_names[, qibs_cols], qibs_path, format = fmt)
        write_table_format(out_names[, qibe_cols], qibe_path, format = fmt)
        cat("Written QIBS:", qibs_path, "\n")
        cat("Written QIBE:", qibe_path, "\n")
    }
}
