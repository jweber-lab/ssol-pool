#!/usr/bin/env Rscript

###############################################################################
# identify_outliers.R
#
# Purpose: Identify outlier candidate windows from diversity (pi, theta, tajima_d),
# FST, and PBE statistics based on quantile thresholds. Use --statistics to select
# which stats to process (pi, theta, tajima_d, pbe, fst). Outputs summary
# statistics and outlier windows/regions to CSV.
# Usage: Rscript identify_outliers.R --hdf5-dir DIR --output-dir DIR
#        --high-quantile Q --low-quantile Q [--statistics pi,theta,tajima_d] [options]
#
# When using --hdf5-dir, expected files are collate output: diversity_w*_s*.h5,
# fst_w*.h5 or fst_single.h5, pbe_*.h5, with group /windows or /sites and columns
# as documented in collate.R (e.g. chr, start, end, sample/pop1,pop2, stat values,
# *_rank, *_quantile, and where present: n_snps, n_total, mean_coverage, mean_mapping_quality).
#
# Output files (see identify_outliers.sh for option summary):
#   outlier_windows_*.csv   One row per (chr, start, end) that is an outlier in at least one statistic.
#     Columns: chr, start, end; window_size; optional mean_coverage, mean_mapping_quality, n_snps;
#     outlier_stat (colon-separated stats, e.g. pi:fst); outlier_direction (high, low, or high:low);
#     optional window_type (seed/expanded when seed-expand mode); then value columns
#     samplename.statname (e.g. Cheney.pi) and samplename.statname_quantile (window quantile 0-1);
#     sample1:sample2.fst, sample1:sample2:sample3.pbe and their _quantile columns.
#   outlier_regions_*.csv   One row per merged/expanded region when --merge-distance or
#     seed-expand is used. Filename includes _across_samples or _within_samples and
#     _seed_expand or _merge_only. Columns: chr, region_start, region_end, n_windows,
#     outlier_stat; optional region_mean_coverage, region_mean_mapping_quality, region_n_snps;
#     per-stat summaries (for each value column): region_mean_<name> (overlap-corrected),
#     region_max_<name>, region_max_<name>_quantile (quantile of the most extreme window).
#
# How regions are built (see also --merge-across-samples):
#   (1) Merge-only path (no --expand-*): merge_outlier_windows_to_regions() groups
#       outlier windows by chr and by gap; when gap between consecutive (sorted by start)
#       windows > merge_distance a new region starts. region_start/end = min(start)/max(end)
#       of that group. With --merge-across-samples true, all samples are merged together;
#       with false, merging is done separately per sample (or per sample_pair for FST/PBE).
#   (2) Seed-expand path (--expand-high-quantile etc.): For each (stat, sample) we call
#       seed_expand_regions(); then with --merge-across-samples true, overlapping/nearby
#       regions from different samples are merged into one set; with false, we keep one
#       row per distinct (chr, region_start, region_end) from any sample (so regions can
#       overlap when different samples produced different extents). n_windows = number of
#       *rows* in the group or overlapping tagged windows (see below), not region_span/step.
#   n_windows: In merge-only path it is the count of outlier (window x sample) rows in
#       the merged group. In seed-expand path, add_region_aggregates() sets n_windows =
#       number of rows in the tagged windows (seed + expanded) that overlap the region;
#       a window can appear as both "seed" and "expanded", so the same interval may be
#       counted twice. So n_windows is not "unique window positions in [region_start,
#       region_end]" and is not region_span/step (e.g. 458 for 136500 bp at step 500
#       is expected when counting tagged rows, possibly with double-counting).
#
# Summary statistics are not written here; use plot/other scripts for that.
# Coordinates: 1-based inclusive where applicable.
###############################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(optparse)
    library(readr)
    library(tidyr)
})

# Parse command-line arguments
option_list <- list(
    make_option(c("--fst-dir"), type="character", default=NULL,
                help="Directory containing FST TSV files (from calculate_fst.sh)", metavar="DIR"),
    make_option(c("--diversity-dir"), type="character", default=NULL,
                help="Directory containing diversity TSV files (from calculate_pi_theta.sh)", metavar="DIR"),
    make_option(c("--output-dir"), type="character", default=".",
                help="Output directory for results [default: %default]", metavar="DIR"),
    make_option(c("--reference-genome"), type="character", default=NULL,
                help="Reference genome FASTA file (optional, for chromosome lengths)", metavar="FILE"),
    make_option(c("--high-quantile"), type="numeric", default=NULL,
                help="High quantile threshold (e.g., 0.99 for top 1%) [required]", metavar="NUMBER"),
    make_option(c("--low-quantile"), type="numeric", default=NULL,
                help="Low quantile threshold (e.g., 0.01 for bottom 1%) [required]", metavar="NUMBER"),
    make_option(c("--seed-high-quantile"), type="numeric", default=NULL,
                help="Strict high quantile for seed windows (seed-then-expand mode; requires --expand-high-quantile)", metavar="NUMBER"),
    make_option(c("--seed-low-quantile"), type="numeric", default=NULL,
                help="Strict low quantile for seed windows (seed-then-expand mode; requires --expand-low-quantile)", metavar="NUMBER"),
    make_option(c("--expand-high-quantile"), type="numeric", default=NULL,
                help="Soft high quantile for expanding regions (seed-then-expand mode)", metavar="NUMBER"),
    make_option(c("--expand-low-quantile"), type="numeric", default=NULL,
                help="Soft low quantile for expanding regions (seed-then-expand mode)", metavar="NUMBER"),
    make_option(c("--top-n-chromosomes"), type="numeric", default=NULL,
                help="Process only the N longest chromosomes/scaffolds (optional)", metavar="N"),
    make_option(c("--min-chromosome-length"), type="numeric", default=NULL,
                help="Process only chromosomes/scaffolds longer than N bp (optional)", metavar="N"),
    make_option(c("--chromosome"), type="character", default=NULL,
                help="Restrict to a single chromosome/scaffold (optional)", metavar="CHR"),
    make_option(c("--region-start"), type="numeric", default=NULL,
                help="Restrict to interval start on --chromosome (requires --chromosome)", metavar="N"),
    make_option(c("--region-end"), type="numeric", default=NULL,
                help="Restrict to interval end on --chromosome (requires --chromosome)", metavar="N"),
    make_option(c("--window-size"), type="numeric", default=NULL,
                help="Filter to specific window size (default: all window sizes)", metavar="NUMBER"),
    make_option(c("--statistics"), type="character", default="pi,theta,tajima_d",
                help="Comma-separated list: pi, theta, tajima_d, pbe, fst (default: pi,theta,tajima_d). Diversity stats use diversity data; fst and pbe use FST/PBE data when provided.", metavar="LIST"),
    make_option(c("--sample-pairs"), type="character", default="",
                help="Comma-separated list of FST pairs to process (default: all)", metavar="LIST"),
    make_option(c("--top-n-extreme"), type="numeric", default=NULL,
                help="Return only the N most extreme values in each quantile (optional)", metavar="N"),
    make_option(c("--hdf5-dir"), type="character", default=NULL,
                help="Directory containing collated HDF5 (diversity_w*.h5, fst_w*.h5); used instead of --diversity-dir/--fst-dir when set", metavar="DIR"),
    make_option(c("--merge-distance"), type="character", default="auto",
                help="Merge nearby outlier windows into regions: integer bp, or 'auto' = max(2*window_size, 2000) [default: auto]", metavar="NUMBER"),
    make_option(c("--merge-across-samples"), type="character", default="true",
                help="When building regions: merge windows/regions across samples (true) or only within each sample (false). Affects both merge-only and seed-expand modes [default: true]", metavar="true|false"),
    make_option(c("--min-depth"), type="numeric", default=NULL,
                help="Minimum mean coverage (depth) for a window to be a seed or in expansion [optional]", metavar="NUMBER"),
    make_option(c("--max-depth"), type="numeric", default=NULL,
                help="Maximum mean coverage (depth) for a window to be a seed or in expansion [optional]", metavar="NUMBER"),
    make_option(c("--min-mapping-quality"), type="numeric", default=NULL,
                help="Minimum mean mapping quality for a window to be a seed or in expansion [optional]", metavar="NUMBER"),
    make_option(c("--max-mapping-quality"), type="numeric", default=NULL,
                help="Maximum mean mapping quality for a window to be a seed or in expansion [optional]", metavar="NUMBER"),
    make_option(c("--min-snps"), type="numeric", default=NULL,
                help="Minimum SNPs per window for a window to be a seed or in expansion [optional]", metavar="NUMBER"),
    make_option(c("--max-snps"), type="numeric", default=NULL,
                help="Maximum SNPs per window for a window to be a seed or in expansion [optional]", metavar="NUMBER"),
    make_option(c("--min-region-snps"), type="numeric", default=NULL,
                help="Minimum total SNPs in a region to output the region [optional]", metavar="NUMBER"),
    make_option(c("--max-region-snps"), type="numeric", default=NULL,
                help="Maximum total SNPs in a region to output the region [optional]", metavar="NUMBER"),
    make_option(c("--min-total-sites"), type="numeric", default=NULL,
                help="Minimum total sites per window (n_total: passed + invariant) for seed/expansion [optional]", metavar="NUMBER"),
    make_option(c("--max-total-sites"), type="numeric", default=NULL,
                help="Maximum total sites per window (n_total) for seed/expansion [optional]", metavar="NUMBER"),
    make_option(c("--verbose"), action="store_true", default=FALSE,
                help="Enable verbose output for debugging", metavar="FLAG")
)

opt_parser <- OptionParser(option_list=option_list, usage="usage: %prog [options]")
opts <- parse_args(opt_parser)
# Coerce --merge-across-samples to logical (default TRUE)
merge_across_val <- tolower(trimws(opts$`merge-across-samples`))
opts$merge_across_samples <- !identical(merge_across_val, "false") && !identical(merge_across_val, "0")

# Validate arguments
if (is.null(opts$`fst-dir`) && is.null(opts$`diversity-dir`) && is.null(opts$`hdf5-dir`)) {
    stop("At least one of --fst-dir, --diversity-dir, or --hdf5-dir must be specified")
}

expand_mode <- !is.null(opts$`expand-high-quantile`) && !is.null(opts$`expand-low-quantile`)
if (expand_mode) {
    if (is.null(opts$`seed-high-quantile`) || is.null(opts$`seed-low-quantile`)) {
        stop("Seed-then-expand requires --seed-high-quantile and --seed-low-quantile when --expand-high-quantile and --expand-low-quantile are set")
    }
    if (opts$`seed-high-quantile` < opts$`expand-high-quantile`) {
        stop("--seed-high-quantile must be >= --expand-high-quantile")
    }
    if (opts$`seed-low-quantile` > opts$`expand-low-quantile`) {
        stop("--seed-low-quantile must be <= --expand-low-quantile")
    }
    if (is.null(opts$`high-quantile`)) opts$`high-quantile` <- opts$`expand-high-quantile`
    if (is.null(opts$`low-quantile`)) opts$`low-quantile` <- opts$`expand-low-quantile`
} else {
    if (is.null(opts$`high-quantile`)) {
        stop("--high-quantile is required (or use seed/expand quantiles for seed-then-expand mode)")
    }
    if (is.null(opts$`low-quantile`)) {
        stop("--low-quantile is required (or use seed/expand quantiles for seed-then-expand mode)")
    }
}

if (opts$`high-quantile` < 0 || opts$`high-quantile` > 1) {
    stop("--high-quantile must be between 0 and 1")
}
if (opts$`low-quantile` < 0 || opts$`low-quantile` > 1) {
    stop("--low-quantile must be between 0 and 1")
}

if ((!is.null(opts$`region-start`) || !is.null(opts$`region-end`)) && is.null(opts$chromosome)) {
    stop("--region-start and --region-end require --chromosome")
}

# Parse statistics list
statistics <- strsplit(opts$statistics, ",")[[1]]
statistics <- trimws(statistics)
valid_stats <- c("pi", "theta", "tajima_d", "pbe", "fst")
if (!all(statistics %in% valid_stats)) {
    stop("--statistics must be one or more of: pi, theta, tajima_d, pbe, fst")
}

# Create output directory if it doesn't exist
if (!dir.exists(opts$`output-dir`)) {
    dir.create(opts$`output-dir`, recursive=TRUE)
}

# Function to read FASTA index (.fai) to get chromosome lengths
read_fai <- function(fai_file) {
    if (!file.exists(fai_file)) {
        return(NULL)
    }
    fai <- read.table(fai_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(fai) <- c("chr", "length", "offset", "linebases", "linewidth")
    return(fai[, c("chr", "length")])
}

# Function to get chromosome lengths
get_chr_lengths <- function(data, reference_genome=NULL) {
    if (!is.null(reference_genome)) {
        # Try to read .fai file
        fai_file <- paste0(reference_genome, ".fai")
        if (file.exists(fai_file)) {
            fai <- read_fai(fai_file)
            if (!is.null(fai)) {
                if (opts$verbose) cat("Using chromosome lengths from FASTA index:", fai_file, "\n")
                return(fai)
            }
        }
        # Try .dict file
        dict_file <- sub("\\.fa$", ".dict", reference_genome)
        dict_file <- sub("\\.fasta$", ".dict", dict_file)
        if (file.exists(dict_file)) {
            # Parse SAM dict format (simplified)
            dict_lines <- readLines(dict_file)
            chr_lines <- dict_lines[grepl("^@SQ", dict_lines)]
            if (length(chr_lines) > 0) {
                chr_info <- data.frame(
                    chr = gsub(".*SN:([^\\t]+).*", "\\1", chr_lines),
                    length = as.numeric(gsub(".*LN:([^\\t]+).*", "\\1", chr_lines)),
                    stringsAsFactors=FALSE
                )
                if (opts$verbose) cat("Using chromosome lengths from SAM dict:", dict_file, "\n")
                return(chr_info)
            }
        }
    }
    
    # Fall back to data-inferred chromosome lengths
    if (opts$verbose) cat("Inferring chromosome lengths from data...\n")
    chr_lengths <- data %>%
        group_by(.data$chr) %>%
        summarise(length = max(.data$pos, na.rm=TRUE), .groups="drop") %>%
        filter(!is.na(.data$length), .data$length > 0) %>%
        arrange(desc(.data$length))  # Order by length, longest first
    
    if (opts$verbose) cat("Found", nrow(chr_lengths), "chromosomes (ordered by length, longest first)\n")
    return(chr_lengths)
}

# Strip diversity subdir-style prefix from sample name (e.g. diversity_w1000_s500_sampleCheney -> Cheney)
# Vectorized: safe to use in mutate(sample = normalize_diversity_sample_from_h5(.data$sample))
normalize_diversity_sample_from_h5 <- function(sample_name) {
    vapply(as.character(sample_name), function(one) {
        s <- one
        if (grepl("^diversity_w[0-9]+_s[0-9]+_sample", s, ignore.case = TRUE)) {
            s <- sub("^diversity_w[0-9]+_s[0-9]+_sample", "", s, ignore.case = TRUE)
        }
        if (nchar(trimws(s)) == 0 || s == one) {
            return(normalize_sample_name(one))
        }
        normalize_sample_name(s)
    }, character(1), USE.NAMES = FALSE)
}

# Helper function to normalize sample names (from plot_fst.R)
normalize_sample_name <- function(sample_name) {
    normalized <- tolower(sample_name)
    # Remove .dedup suffix
    normalized <- sub("\\.dedup$", "", normalized, ignore.case=TRUE)
    # Remove _all_seq suffix
    normalized <- sub("_all_seq$", "", normalized, ignore.case=TRUE)
    # Remove _pool_s[0-9]+ pattern (e.g., _pool_s3)
    normalized <- sub("_pool_s[0-9]+.*$", "", normalized, ignore.case=TRUE)
    # Remove _pool suffix
    normalized <- sub("_pool$", "", normalized, ignore.case=TRUE)
    
    # Try to extract shorter sample name (e.g., "Cheney" from "cheney_pool_s3")
    # If there's an underscore, try to get the first meaningful part
    if (grepl("_", normalized)) {
        parts <- strsplit(normalized, "_")[[1]]
        # Take the first part if it looks like a name (starts with letter)
        if (length(parts) > 0 && grepl("^[a-z]", parts[1], ignore.case=TRUE)) {
            normalized <- parts[1]
        }
    }
    
    # Capitalize first letter for consistency
    if (nchar(normalized) > 0) {
        normalized <- paste0(toupper(substring(normalized, 1, 1)), substring(normalized, 2))
    }
    
    return(normalized)
}

# Helper function to extract window size from filename
extract_window_size_from_filename <- function(filename) {
    window_match <- regmatches(filename, regexpr("w(\\d+)_s(\\d+)", filename))
    if (length(window_match) > 0) {
        window_parts <- strsplit(window_match, "_")[[1]]
        window_size <- as.numeric(sub("w", "", window_parts[1]))
        return(window_size)
    }
    # Try simpler pattern (just w followed by digits)
    window_match <- regmatches(filename, regexpr("w(\\d+)", filename))
    if (length(window_match) > 0) {
        return(as.numeric(sub("w", "", window_match)))
    }
    return(NA)
}

# Read diversity data from collated HDF5 (diversity_w*.h5 with /windows group)
read_diversity_from_h5 <- function(hdf5_dir, selected_window_size = NULL) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Package hdf5r required for --hdf5-dir. Install with: install.packages(\"hdf5r\")")
    }
    files <- list.files(hdf5_dir, pattern = "^diversity.*\\.h5$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    out <- list()
    for (f in files) {
        win_size <- NA
        m <- regmatches(basename(f), regexpr("w([0-9]+)", basename(f), ignore.case = TRUE))
        if (length(m)) win_size <- as.numeric(sub("w([0-9]+).*", "\\1", m))
        if (!is.null(selected_window_size) && !is.na(win_size) && win_size != selected_window_size) next
        h5 <- hdf5r::H5File$new(f, mode = "r")
        on.exit(h5$close_all(), add = TRUE)
        if (!"windows" %in% names(h5)) { h5$close_all(); next }
        grp <- h5[["windows"]]
        d <- data.frame(
            chr = grp[["chr"]][],
            start = grp[["start"]][],
            end = grp[["end"]][],
            sample = grp[["sample"]][],
            stringsAsFactors = FALSE
        )
        d$pos <- d$start
        d$window_size <- win_size
        for (col in c("pi", "theta", "tajima_d")) {
            if (col %in% names(grp)) d[[col]] <- grp[[col]][]
            qcol <- paste0(col, "_quantile")
            if (qcol %in% names(grp)) d[[qcol]] <- grp[[qcol]][]
        }
        for (col in c("mean_coverage", "mean_mapping_quality", "n_snps", "n_total")) {
            if (col %in% names(grp)) d[[col]] <- grp[[col]][]
        }
        h5$close_all()
        on.exit()
        out[[length(out) + 1L]] <- d
    }
    if (length(out) == 0) return(NULL)
    combined_data <- bind_rows(out)
    combined_data <- combined_data %>% mutate(sample = normalize_diversity_sample_from_h5(.data$sample))
    combined_data %>% filter(!is.na(.data$chr), !is.na(.data$pos))
}

# Read FST data from collated HDF5 (fst_w*.h5 or fst_single.h5). Group /windows or /sites.
# Supports wide format (fst_Sample1_Sample2) and long format (pop1, pop2, fst).
read_fst_from_h5 <- function(hdf5_dir, selected_window_size = NULL, selected_pairs = NULL) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Package hdf5r required for --hdf5-dir. Install with: install.packages(\"hdf5r\")")
    }
    files <- list.files(hdf5_dir, pattern = "^fst.*\\.h5$", full.names = TRUE)
    if (length(files) == 0) return(list(data = NULL, pairs = list()))
    out <- list()
    pair_names <- character(0)
    for (f in files) {
        win_size <- NA
        m <- regmatches(basename(f), regexpr("w([0-9]+)", basename(f), ignore.case = TRUE))
        if (length(m)) win_size <- as.numeric(sub("w([0-9]+).*", "\\1", m))
        if (!is.null(selected_window_size) && !is.na(win_size) && win_size != selected_window_size) next
        h5 <- hdf5r::H5File$new(f, mode = "r")
        on.exit(h5$close_all(), add = TRUE)
        grp_name <- if ("windows" %in% names(h5)) "windows" else if ("sites" %in% names(h5)) "sites" else NULL
        if (is.null(grp_name)) { h5$close_all(); next }
        grp <- h5[[grp_name]]
        grp_names <- names(grp)
        # Long format: pop1, pop2, fst, fst_quantile
        if ("pop1" %in% grp_names && "pop2" %in% grp_names && "fst" %in% grp_names) {
            d <- data.frame(
                chr = grp[["chr"]][],
                start = grp[["start"]][],
                end = grp[["end"]][],
                pos = grp[["start"]][],
                window_size = win_size,
                fst = grp[["fst"]][],
                sample_pair = paste0(grp[["pop1"]][], ":", grp[["pop2"]][]),
                stringsAsFactors = FALSE
            )
            if ("fst_quantile" %in% grp_names) d$fst_quantile <- grp[["fst_quantile"]][]
            if ("n_snps" %in% grp_names) d$n_snps <- grp[["n_snps"]][]
            if ("step_size" %in% grp_names) d$step_size <- grp[["step_size"]][]
            pair_names <- c(pair_names, unique(d$sample_pair))
            if (!is.null(selected_pairs) && length(selected_pairs) > 0) {
                d <- d %>% filter(.data$sample_pair %in% selected_pairs)
                if (nrow(d) == 0) next
            }
            out[[length(out) + 1L]] <- d
            h5$close_all()
            on.exit()
            next
        }
        # Wide format: fst_Sample1_Sample2, fst_Sample1_Sample2_quantile
        fst_cols <- grep("^fst_.+_.+$", grp_names, value = TRUE)
        fst_cols <- fst_cols[!grepl("_rank$|_quantile$", fst_cols)]
        for (fc in fst_cols) {
            pair_name <- sub("^fst_", "", fc)
            pair_name <- gsub("_", ":", pair_name)
            pair_names <- c(pair_names, pair_name)
            if (!is.null(selected_pairs) && length(selected_pairs) > 0 &&
                !pair_name %in% selected_pairs) next
            d <- data.frame(
                chr = grp[["chr"]][],
                start = grp[["start"]][],
                end = grp[["end"]][],
                pos = grp[["start"]][],
                window_size = win_size,
                fst = grp[[fc]][],
                sample_pair = pair_name,
                stringsAsFactors = FALSE
            )
            qcol <- paste0(fc, "_quantile")
            if (qcol %in% grp_names) d$fst_quantile <- grp[[qcol]][]
            out[[length(out) + 1L]] <- d
        }
        h5$close_all()
        on.exit()
    }
    if (length(out) == 0) return(list(data = NULL, pairs = list()))
    unique_pairs <- unique(pair_names)
    pairs_named <- setNames(as.list(rep("fst", length(unique_pairs))), unique_pairs)
    list(data = bind_rows(out), pairs = pairs_named)
}

# Read PBE data from collated HDF5 (pbe_w*.h5 or pbe_single.h5). Group /windows or /sites.
read_pbe_from_h5 <- function(hdf5_dir, selected_window_size = NULL) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Package hdf5r required for PBE from HDF5. Install with: install.packages(\"hdf5r\")")
    }
    files <- list.files(hdf5_dir, pattern = "^pbe.*\\.h5$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    out <- list()
    for (f in files) {
        grp_name <- NULL
        h5 <- hdf5r::H5File$new(f, mode = "r")
        on.exit(h5$close_all(), add = TRUE)
        if ("windows" %in% names(h5)) grp_name <- "windows" else if ("sites" %in% names(h5)) grp_name <- "sites"
        if (is.null(grp_name)) { h5$close_all(); next }
        grp <- h5[[grp_name]]
        win_size <- NA
        m <- regmatches(basename(f), regexpr("w([0-9]+)", basename(f), ignore.case = TRUE))
        if (length(m)) win_size <- as.numeric(sub("w([0-9]+).*", "\\1", m))
        if (!is.null(selected_window_size) && !is.na(win_size) && win_size != selected_window_size) { h5$close_all(); next }
        if (!"pbe" %in% names(grp) || !"pop1" %in% names(grp)) { h5$close_all(); next }
        d <- data.frame(
            chr = grp[["chr"]][],
            start = grp[["start"]][],
            end = grp[["end"]][],
            pos = grp[["start"]][],
            window_size = win_size,
            pop1 = as.character(grp[["pop1"]][]),
            pop2 = as.character(grp[["pop2"]][]),
            pop3 = as.character(grp[["pop3"]][]),
            pbe = as.numeric(grp[["pbe"]][]),
            stringsAsFactors = FALSE
        )
        if ("pbe_quantile" %in% names(grp)) d$pbe_quantile <- grp[["pbe_quantile"]][]
        if ("pbe_rank" %in% names(grp)) d$pbe_rank <- grp[["pbe_rank"]][]
        for (col in c("mean_coverage", "mean_mapping_quality", "n_snps")) {
            if (col %in% names(grp)) d[[col]] <- grp[[col]][]
        }
        d$trio_id <- paste0(d$pop1, ":", d$pop2, ":", d$pop3)
        out[[length(out) + 1L]] <- d
        h5$close_all()
        on.exit()
    }
    if (length(out) == 0) return(NULL)
    bind_rows(out) %>% filter(!is.na(.data$chr), !is.na(.data$pos))
}

# Check that every used data source has columns required by set filter options; stop() with clear message if not.
check_quality_columns_and_abort <- function(data_sources, opts, has_diversity, has_fst, has_pbe) {
    need_coverage <- !is.null(opts$`min-depth`) || !is.null(opts$`max-depth`)
    need_mapq <- !is.null(opts$`min-mapping-quality`) || !is.null(opts$`max-mapping-quality`)
    need_n_snps <- !is.null(opts$`min-snps`) || !is.null(opts$`max-snps`) ||
        !is.null(opts$`min-region-snps`) || !is.null(opts$`max-region-snps`)
    if (!need_coverage && !need_mapq && !need_n_snps) return(invisible(NULL))
    check <- function(df, label) {
        if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
        if (need_coverage && !"mean_coverage" %in% colnames(df)) {
            stop("Filter --min-depth/--max-depth requires column 'mean_coverage' in all used data. ",
                 "Provide collated HDF5 with diversity (or data that includes mean_coverage) and re-run, or omit the depth filter.")
        }
        if (need_mapq && !"mean_mapping_quality" %in% colnames(df)) {
            stop("Filter --min-mapping-quality/--max-mapping-quality requires column 'mean_mapping_quality' in all used data. ",
                 "Provide collated HDF5 with diversity (or data that includes mean_mapping_quality) and re-run, or omit the mapping quality filter.")
        }
        if (need_n_snps && !"n_snps" %in% colnames(df)) {
            stop("Filter --min-snps/--max-snps/--min-region-snps/--max-region-snps requires column 'n_snps' in all used data. ",
                 "Provide collated HDF5 with diversity (or data that includes n_snps) and re-run, or omit the SNP filter(s).")
        }
    }
    if (has_diversity) check(data_sources$diversity, "diversity")
    if (has_fst) check(data_sources$fst, "FST")
    if (has_pbe) check(data_sources$pbe, "PBE")
    invisible(NULL)
}

# Attach region_mean_coverage, region_mean_mapping_quality, region_n_snps to regions by aggregating overlapping windows.
add_region_aggregates <- function(regions, windows) {
    if (is.null(regions) || nrow(regions) == 0) return(regions)
    if (!"start" %in% colnames(windows)) windows$start <- windows$pos
    if (!"end" %in% colnames(windows)) windows$end <- windows$pos
    out_list <- list()
    for (i in seq_len(nrow(regions))) {
        r <- regions[i, ]
        w <- windows %>%
            filter(.data$chr == r$chr, .data$start <= r$region_end, .data$end >= r$region_start)
        r$n_windows <- nrow(w)
        if ("mean_coverage" %in% colnames(w)) r$region_mean_coverage <- mean(w$mean_coverage, na.rm = TRUE)
        if ("mean_mapping_quality" %in% colnames(w)) r$region_mean_mapping_quality <- mean(w$mean_mapping_quality, na.rm = TRUE)
        if ("n_snps" %in% colnames(w)) r$region_n_snps <- sum(w$n_snps, na.rm = TRUE)
        out_list[[i]] <- r
    }
    bind_rows(out_list)
}

# For each (region, value_col): overlap-corrected mean stat, max stat, quantile of the window with max.
# windows_by_value_col: named list of tibbles with chr, start, end, value, and optional quantile column.
# Weights: 1 / (number of windows overlapping this window); then weighted mean = sum(value*weight)/sum(weight).
add_region_stat_summaries <- function(regions, windows_by_value_col) {
    if (is.null(regions) || nrow(regions) == 0) return(regions)
    if (length(windows_by_value_col) == 0) return(regions)
    for (vc in names(windows_by_value_col)) {
        w <- windows_by_value_col[[vc]]
        if (is.null(w) || nrow(w) == 0 || !"value" %in% colnames(w)) next
        qcol <- intersect(c("quantile", paste0(gsub(".*\\.", "", vc), "_quantile")), colnames(w))[1L]
        if (length(qcol) == 0 || is.na(qcol)) qcol <- NULL
        safe_name <- gsub("[:.]", "_", vc)
        mean_col <- paste0("region_mean_", safe_name)
        max_col <- paste0("region_max_", safe_name)
        max_quant_col <- paste0("region_max_", safe_name, "_quantile")
        out_means <- numeric(nrow(regions))
        out_maxs <- numeric(nrow(regions))
        out_max_quants <- rep(NA_real_, nrow(regions))
        for (i in seq_len(nrow(regions))) {
            r <- regions[i, ]
            sub <- w %>%
                filter(.data$chr == r$chr, .data$start <= r$region_end, .data$end >= r$region_start) %>%
                filter(!is.na(.data$value))
            if (nrow(sub) == 0) {
                out_means[i] <- NA_real_
                out_maxs[i] <- NA_real_
                next
            }
            # Overlap count: window i overlaps j when start[i] <= end[j] and end[i] >= start[j]; vectorized
            n_overlap <- colSums(outer(sub$start, sub$end, "<=") & outer(sub$end, sub$start, ">="))
            weight <- 1 / pmax(n_overlap, 1)
            out_means[i] <- sum(sub$value * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE)
            idx_max <- which.max(sub$value)
            out_maxs[i] <- sub$value[idx_max]
            if (!is.null(qcol) && qcol %in% colnames(sub)) out_max_quants[i] <- sub[[qcol]][idx_max]
        }
        regions[[mean_col]] <- out_means
        regions[[max_col]] <- out_maxs
        regions[[max_quant_col]] <- out_max_quants
    }
    regions
}

# Build list of window tibbles (chr, start, end, value, quantile) keyed by value_col for add_region_stat_summaries.
# value_cols: character vector of names like "Cheney.pi", "S1:S2.fst", "S1:S2:S3.pbe".
build_windows_by_value_col <- function(value_cols, diversity_data, fst_data, pbe_data) {
    out <- list()
    if (length(value_cols) == 0) return(out)
    for (vc in value_cols) {
        if (grepl("\\.pbe$", vc)) {
            trio <- sub("\\.pbe$", "", vc)
            if (is.null(pbe_data) || nrow(pbe_data) == 0) next
            if ("trio_id" %in% colnames(pbe_data)) {
                sub <- pbe_data %>% filter(.data$trio_id == trio)
            } else {
                trio_parts <- strsplit(trio, ":", fixed = TRUE)[[1]]
                if (length(trio_parts) != 3) next
                sub <- pbe_data %>%
                    filter(.data$pop1 == trio_parts[1], .data$pop2 == trio_parts[2], .data$pop3 == trio_parts[3])
            }
            sub <- sub %>%
                transmute(chr = .data$chr, start = .data$start, end = .data$end, value = .data$pbe,
                         quantile = if ("pbe_quantile" %in% colnames(pbe_data)) .data$pbe_quantile else NA_real_)
            if (nrow(sub) > 0) out[[vc]] <- sub
        } else if (grepl("\\.fst$", vc)) {
            pair <- sub("\\.fst$", "", vc)
            if (is.null(fst_data) || nrow(fst_data) == 0) next
            sub <- fst_data %>% filter(.data$sample_pair == pair) %>%
                transmute(chr = .data$chr, start = .data$start, end = .data$end, value = .data$fst,
                         quantile = if ("fst_quantile" %in% colnames(fst_data)) .data$fst_quantile else NA_real_)
            if (nrow(sub) > 0) out[[vc]] <- sub
        } else {
            # diversity: "Cheney.pi", "Echo.theta" (last segment is stat, rest is sample)
            parts <- strsplit(vc, ".", fixed = TRUE)[[1]]
            if (length(parts) < 2) next
            stat_col <- parts[length(parts)]
            sample_name <- paste(parts[-length(parts)], collapse = ".")
            if (is.null(diversity_data) || nrow(diversity_data) == 0 || !stat_col %in% colnames(diversity_data)) next
            qcol <- paste0(stat_col, "_quantile")
            sub <- diversity_data %>%
                filter(.data$sample == sample_name) %>%
                transmute(chr = .data$chr, start = .data$start, end = .data$end,
                         value = !!sym(stat_col),
                         quantile = if (qcol %in% colnames(diversity_data)) !!sym(qcol) else NA_real_)
            if (nrow(sub) > 0) out[[vc]] <- sub
        }
    }
    out
}

# Standardize outlier table: add quantile from source column, rename value column, select out_cols + extra_cols.
# Use when building outlier_results so all branches produce the same column set.
standardize_outlier_columns <- function(outliers, value_col, quantile_col, out_cols, extra_cols) {
    outliers %>%
        mutate(quantile = if (quantile_col %in% colnames(outliers)) !!sym(quantile_col) else NA_real_) %>%
        rename(value = !!sym(value_col)) %>%
        select(any_of(c(out_cols, extra_cols)))
}

# Apply region-level depth, mapping quality, and region SNP filters from opts. Returns filtered regions df.
filter_regions_by_opts <- function(regions, opts) {
    if (is.null(regions) || nrow(regions) == 0) return(regions)
    if (!is.null(opts$`min-depth`) && "region_mean_coverage" %in% colnames(regions)) {
        regions <- regions %>% filter(.data$region_mean_coverage >= opts$`min-depth` | is.na(.data$region_mean_coverage))
    }
    if (!is.null(opts$`max-depth`) && "region_mean_coverage" %in% colnames(regions)) {
        regions <- regions %>% filter(.data$region_mean_coverage <= opts$`max-depth` | is.na(.data$region_mean_coverage))
    }
    if (!is.null(opts$`min-mapping-quality`) && "region_mean_mapping_quality" %in% colnames(regions)) {
        regions <- regions %>% filter(.data$region_mean_mapping_quality >= opts$`min-mapping-quality` | is.na(.data$region_mean_mapping_quality))
    }
    if (!is.null(opts$`max-mapping-quality`) && "region_mean_mapping_quality" %in% colnames(regions)) {
        regions <- regions %>% filter(.data$region_mean_mapping_quality <= opts$`max-mapping-quality` | is.na(.data$region_mean_mapping_quality))
    }
    if (!is.null(opts$`min-region-snps`) && "region_n_snps" %in% colnames(regions)) {
        regions <- regions %>% filter(.data$region_n_snps >= opts$`min-region-snps` | is.na(.data$region_n_snps))
    }
    if (!is.null(opts$`max-region-snps`) && "region_n_snps" %in% colnames(regions)) {
        regions <- regions %>% filter(.data$region_n_snps <= opts$`max-region-snps` | is.na(.data$region_n_snps))
    }
    regions
}

# Apply optional depth, mapping quality, and per-window SNP filters. Returns TRUE for rows that pass.
# If a threshold is set and the column is missing, the row fails (excluded).
passes_quality_filters <- function(df, opts) {
    if (nrow(df) == 0) return(logical(0))
    pass <- rep(TRUE, nrow(df))
    if (!is.null(opts$`min-depth`)) {
        if ("mean_coverage" %in% colnames(df)) {
            pass <- pass & (df$mean_coverage >= opts$`min-depth` | is.na(df$mean_coverage))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    if (!is.null(opts$`max-depth`)) {
        if ("mean_coverage" %in% colnames(df)) {
            pass <- pass & (df$mean_coverage <= opts$`max-depth` | is.na(df$mean_coverage))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    if (!is.null(opts$`min-mapping-quality`)) {
        if ("mean_mapping_quality" %in% colnames(df)) {
            pass <- pass & (df$mean_mapping_quality >= opts$`min-mapping-quality` | is.na(df$mean_mapping_quality))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    if (!is.null(opts$`max-mapping-quality`)) {
        if ("mean_mapping_quality" %in% colnames(df)) {
            pass <- pass & (df$mean_mapping_quality <= opts$`max-mapping-quality` | is.na(df$mean_mapping_quality))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    if (!is.null(opts$`min-snps`)) {
        if ("n_snps" %in% colnames(df)) {
            pass <- pass & (df$n_snps >= opts$`min-snps` | is.na(df$n_snps))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    if (!is.null(opts$`max-snps`)) {
        if ("n_snps" %in% colnames(df)) {
            pass <- pass & (df$n_snps <= opts$`max-snps` | is.na(df$n_snps))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    if (!is.null(opts$`min-total-sites`)) {
        if ("n_total" %in% colnames(df)) {
            pass <- pass & (df$n_total >= opts$`min-total-sites` | is.na(df$n_total))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    if (!is.null(opts$`max-total-sites`)) {
        if ("n_total" %in% colnames(df)) {
            pass <- pass & (df$n_total <= opts$`max-total-sites` | is.na(df$n_total))
        } else {
            pass <- rep(FALSE, nrow(df))
        }
    }
    pass
}

# Merge overlapping regions (data frame with chr, start, end or region_start, region_end). Merges when gap <= max_gap.
merge_regions <- function(regions_df, max_gap = 0) {
    if (is.null(regions_df) || nrow(regions_df) == 0) return(NULL)
    nms <- names(regions_df)
    start_col <- intersect(c("region_start", "start"), nms)[1L]
    end_col <- intersect(c("region_end", "end"), nms)[1L]
    if (is.na(start_col) || is.na(end_col)) return(regions_df)
    by_chr <- regions_df %>%
        rename(region_start = !!sym(start_col), region_end = !!sym(end_col)) %>%
        arrange(.data$chr, .data$region_start)
    out_list <- list()
    for (ch in unique(by_chr$chr)) {
        sub <- by_chr %>% filter(.data$chr == ch)
        gap <- c(Inf, sub$region_start[-1] - sub$region_end[-nrow(sub)])
        sub$grp <- cumsum(gap > max_gap)
        reg <- sub %>%
            group_by(.data$grp) %>%
            summarise(chr = .data$chr[1], region_start = min(.data$region_start, na.rm = TRUE), region_end = max(.data$region_end, na.rm = TRUE), .groups = "drop") %>%
            select(-"grp")
        out_list[[length(out_list) + 1L]] <- reg
    }
    bind_rows(out_list)
}

# Seed-then-expand: merge seed windows to regions, expand each region with expand windows within merge_distance, merge overlapping expanded regions. Returns list(regions = data frame of expanded regions, windows = tagged windows with window_type "seed" or "expanded").
seed_expand_regions <- function(seed_windows, expand_windows, merge_distance, id_cols = NULL, verbose = FALSE) {
    if (is.null(seed_windows) || nrow(seed_windows) == 0) return(list(regions = NULL, windows = seed_windows %>% mutate(window_type = "seed")))
    seed_regions <- merge_outlier_windows_to_regions(seed_windows, merge_distance)
    if (is.null(seed_regions) || nrow(seed_regions) == 0) return(list(regions = NULL, windows = seed_windows %>% mutate(window_type = "seed")))
    start_col <- if ("region_start" %in% names(seed_regions)) "region_start" else "start"
    end_col <- if ("region_end" %in% names(seed_regions)) "region_end" else "end"
    seed_regions <- seed_regions %>% rename(region_start = !!sym(start_col), region_end = !!sym(end_col))
    expanded_list <- list()
    for (i in seq_len(nrow(seed_regions))) {
        r <- seed_regions[i, ]
        ch <- r$chr
        rstart <- r$region_start - merge_distance
        rend <- r$region_end + merge_distance
        exp_sub <- expand_windows %>%
            filter(.data$chr == ch, .data$start <= rend, .data$end >= rstart)
        if (nrow(exp_sub) == 0) {
            expanded_list[[length(expanded_list) + 1L]] <- data.frame(chr = ch, region_start = r$region_start, region_end = r$region_end, stringsAsFactors = FALSE)
        } else {
            expanded_list[[length(expanded_list) + 1L]] <- data.frame(chr = ch, region_start = min(c(r$region_start, exp_sub$start), na.rm = TRUE), region_end = max(c(r$region_end, exp_sub$end), na.rm = TRUE), stringsAsFactors = FALSE)
        }
    }
    expanded <- bind_rows(expanded_list)
    expanded <- merge_regions(expanded, max_gap = 0)
    seed_windows$window_type <- "seed"
    n_expand_before <- nrow(expand_windows)
    n_regions <- nrow(expanded)
    expand_used <- expand_windows %>%
        inner_join(expanded %>% select("chr", "region_start", "region_end"), by = "chr", relationship = "many-to-many") %>%
        filter(.data$start <= .data$region_end, .data$end >= .data$region_start) %>%
        select(all_of(names(expand_windows))) %>%
        distinct() %>%
        mutate(window_type = "expanded")
    if (verbose) cat("  seed_expand: expand_windows", n_expand_before, "x expanded regions", n_regions, "(join by chr, many-to-many) -> tagged expanded windows", nrow(expand_used), "\n")
    all_tagged <- bind_rows(seed_windows, expand_used)
    expanded <- add_region_aggregates(expanded, all_tagged)
    if (!is.null(id_cols)) {
        for (c in id_cols) if (c %in% names(seed_regions) && !c %in% names(expanded)) expanded[[c]] <- seed_regions[[c]][1]
    }
    list(regions = expanded, windows = all_tagged)
}

# Safe min/max: return NA when all x are NA (avoids min(..., na.rm=TRUE) returning Inf/-Inf).
safe_min <- function(x) { if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE) }
safe_max <- function(x) { if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE) }

# Merge outlier windows into regions (gap <= merge_distance bp). Optionally compute region-level depth/mapq/n_snps summaries.
merge_outlier_windows_to_regions <- function(all_outliers, merge_distance) {
    if (!"start" %in% colnames(all_outliers)) all_outliers$start <- all_outliers$pos
    if (!"end" %in% colnames(all_outliers)) all_outliers$end <- all_outliers$pos
    id_col <- if ("sample" %in% colnames(all_outliers)) "sample" else "sample_pair"
    if (!id_col %in% colnames(all_outliers)) id_col <- NULL
    by_chr <- all_outliers %>%
        filter(!is.na(.data$chr), .data$chr != "genome") %>%
        arrange(.data$chr, .data$start)
    if (nrow(by_chr) == 0) return(NULL)
    regions_list <- list()
    for (ch in unique(by_chr$chr)) {
        sub <- by_chr %>% filter(.data$chr == ch) %>% arrange(.data$start)
        gap <- c(Inf, sub$start[-1] - sub$end[-nrow(sub)])
        sub$region_id <- cumsum(gap > merge_distance)
        sum_exprs <- list(
            chr = quote(.data$chr[1]),
            region_start = quote(min(.data$start, na.rm = TRUE)),
            region_end = quote(max(.data$end, na.rm = TRUE)),
            n_windows = quote(n())
        )
        if (!is.null(id_col)) {
            sum_exprs$n_populations <- quote(n_distinct(.data[[id_col]], na.rm = TRUE))
        }
        if ("mean_coverage" %in% colnames(sub)) {
            sum_exprs$region_mean_coverage <- quote(mean(.data$mean_coverage, na.rm = TRUE))
            sum_exprs$region_min_coverage <- quote(safe_min(.data$mean_coverage))
            sum_exprs$region_max_coverage <- quote(safe_max(.data$mean_coverage))
        }
        if ("mean_mapping_quality" %in% colnames(sub)) {
            sum_exprs$region_mean_mapping_quality <- quote(mean(.data$mean_mapping_quality, na.rm = TRUE))
            sum_exprs$region_min_mapping_quality <- quote(safe_min(.data$mean_mapping_quality))
            sum_exprs$region_max_mapping_quality <- quote(safe_max(.data$mean_mapping_quality))
        }
        if ("n_snps" %in% colnames(sub)) {
            sum_exprs$region_n_snps <- quote(sum(.data$n_snps, na.rm = TRUE))
            sum_exprs$region_mean_snps_per_window <- quote(mean(.data$n_snps, na.rm = TRUE))
        }
        reg <- sub %>% group_by(.data$region_id) %>%
            summarise(!!!sum_exprs, .groups = "drop") %>%
            select(-"region_id")
        regions_list[[length(regions_list) + 1L]] <- reg
    }
    bind_rows(regions_list)
}

# Helper function to peek at column names in a file (read just header)
peek_column_names <- function(filepath) {
    # Read just first line to get column names
    first_line <- readLines(filepath, n=1)
    has_tabs <- grepl("\t", first_line)
    has_commas <- grepl(",", first_line)
    
    if (has_tabs && (!has_commas || length(strsplit(first_line, "\t")[[1]]) > length(strsplit(first_line, ",")[[1]]))) {
        # Tab-separated
        cols <- strsplit(first_line, "\t")[[1]]
    } else if (has_commas) {
        # Comma-separated
        cols <- strsplit(first_line, ",")[[1]]
    } else {
        # Default to tab
        cols <- strsplit(first_line, "\t")[[1]]
    }
    
    return(tolower(trimws(cols)))
}

# Function to read and process diversity files
read_diversity_files <- function(input_dir, selected_window_size=NULL) {
    if (is.null(input_dir)) {
        return(list())
    }
    
    # Find all diversity files
    tsv_files <- list.files(input_dir, pattern=".*diversity.*\\.tsv$", full.names=TRUE, recursive=TRUE)
    csv_files <- list.files(input_dir, pattern=".*diversity.*\\.csv$", full.names=TRUE, recursive=TRUE)
    all_files <- c(tsv_files, csv_files)
    
    if (length(all_files) == 0) {
        if (opts$verbose) cat("No diversity files found in:", input_dir, "\n")
        return(list())
    }
    
    if (opts$verbose) {
        cat("Found", length(all_files), "diversity file(s) (", length(tsv_files), " TSV, ", length(csv_files), " CSV)\n", sep="")
    }
    
    # Filter files by window size if specified
    if (!is.null(selected_window_size)) {
        file_window_sizes <- sapply(all_files, function(f) extract_window_size_from_filename(basename(f)))
        all_files <- all_files[is.na(file_window_sizes) | file_window_sizes == selected_window_size]
        if (opts$verbose) {
            cat("Filtered to", length(all_files), "file(s) matching window size", selected_window_size, "\n")
        }
        if (length(all_files) == 0) {
            if (opts$verbose) cat("No files found matching window size", selected_window_size, "\n")
            return(list())
        }
    }
    
    # Read and combine all diversity files
    all_data <- list()
    for (file in all_files) {
        if (opts$verbose) {
            cat("\n=== Reading:", basename(file), "===\n")
        }
        
        # Detect separator
        first_line <- readLines(file, n=1)
        has_tabs <- grepl("\t", first_line)
        has_commas <- grepl(",", first_line)
        
        if (has_tabs && (!has_commas || length(strsplit(first_line, "\t")[[1]]) > length(strsplit(first_line, ",")[[1]]))) {
            data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
            detected_sep <- "tab"
        } else if (has_commas) {
            data <- read_csv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
            detected_sep <- "comma"
        } else {
            data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
            detected_sep <- "tab"
        }
        
        # Standardize column names to lowercase
        data <- data %>%
            rename_all(tolower)
        
        # Convert base numeric columns
        data <- data %>%
            mutate(across(any_of(c("start", "end", "position", "pos")), as.numeric))
        
        # Convert statistics columns to numeric
        stat_cols <- grep("\\.(theta_pi|theta_watterson|tajimas?_d)$", colnames(data), ignore.case=FALSE, value=TRUE)
        if (length(stat_cols) > 0) {
            data <- data %>%
                mutate(across(all_of(stat_cols), as.numeric))
        }
        
        # Handle chromosome column
        if (!"chr" %in% colnames(data)) {
            if ("chromosome" %in% colnames(data)) {
                data <- data %>% mutate(chr = as.character(chromosome))
            } else if ("chrom" %in% colnames(data)) {
                data <- data %>% mutate(chr = as.character(chrom))
            } else if ("contig" %in% colnames(data)) {
                data <- data %>% mutate(chr = as.character(contig))
            } else {
                chr_col <- grep("^(chr|chrom|contig|scaffold)", colnames(data), ignore.case=FALSE, value=TRUE)
                if (length(chr_col) > 0) {
                    data <- data %>% mutate(chr = as.character(.data[[chr_col[1]]]))
                } else {
                    cat("  ERROR: Could not find chromosome column\n")
                    next
                }
            }
        } else {
            data <- data %>% mutate(chr = as.character(chr))
        }
        
        # Handle position column
        if (!"pos" %in% colnames(data)) {
            if ("start" %in% colnames(data)) {
                data <- data %>% mutate(pos = as.numeric(start))
            } else if ("position" %in% colnames(data)) {
                data <- data %>% mutate(pos = as.numeric(position))
            } else if ("end" %in% colnames(data)) {
                data <- data %>% mutate(pos = as.numeric(end))
            } else {
                pos_col <- grep("^(pos|start|end|coordinate)", colnames(data), ignore.case=FALSE, value=TRUE)
                if (length(pos_col) > 0) {
                    data <- data %>% mutate(pos = as.numeric(.data[[pos_col[1]]]))
                } else {
                    cat("  ERROR: Could not find position column\n")
                    next
                }
            }
        } else {
            data <- data %>% mutate(pos = as.numeric(pos))
        }
        
        # Extract window size from filename
        filename <- basename(file)
        window_match <- regmatches(filename, regexpr("w(\\d+)_s(\\d+)", filename))
        if (length(window_match) > 0) {
            window_parts <- strsplit(window_match, "_")[[1]]
            window_size <- as.numeric(sub("w", "", window_parts[1]))
            step_size <- as.numeric(sub("s", "", window_parts[2]))
            data$window_size <- window_size
            data$step_size <- step_size
        } else {
            data$window_size <- NA
            data$step_size <- NA
        }
        
        # Extract sample name
        file_dir <- normalizePath(dirname(file), winslash="/", mustWork=FALSE)
        input_dir_norm <- normalizePath(input_dir, winslash="/", mustWork=FALSE)
        dir_parts <- strsplit(file_dir, "/")[[1]]
        input_dir_parts <- strsplit(input_dir_norm, "/")[[1]]
        
        sample_name <- "unknown"
        if (length(dir_parts) > length(input_dir_parts)) {
            sample_idx <- length(input_dir_parts) + 1
            if (sample_idx <= length(dir_parts)) {
                sample_name <- dir_parts[sample_idx]
            }
        } else {
            filename_no_ext <- tools::file_path_sans_ext(basename(file))
            if (grepl("^([A-Za-z0-9_]+)_diversity", filename_no_ext, ignore.case=TRUE)) {
                sample_name <- sub("^([A-Za-z0-9_]+)_diversity.*", "\\1", filename_no_ext, ignore.case=TRUE)
            }
        }
        
        # Extract diversity statistics
        theta_pi_cols <- grep("\\.theta_pi$", colnames(data), ignore.case=FALSE, value=TRUE)
        if (length(theta_pi_cols) > 0) {
            theta_pi_col <- theta_pi_cols[1]
            data <- data %>%
                mutate(pi = as.numeric(.data[[theta_pi_col]]))
            
            if (sample_name == "unknown") {
                sample_from_col <- sub("\\.theta_pi$", "", theta_pi_col, ignore.case=FALSE)
                sample_from_col <- sub("_all_seq\\.dedup$", "", sample_from_col, ignore.case=TRUE)
                sample_from_col <- sub("\\.dedup$", "", sample_from_col, ignore.case=TRUE)
                sample_from_col <- sub("_all_seq$", "", sample_from_col, ignore.case=TRUE)
                if (grepl("^([a-z0-9_]+)_pool_s[0-9]+", sample_from_col, ignore.case=TRUE)) {
                    sample_name <- sub("^([a-z0-9_]+)_pool_s[0-9]+.*", "\\1", sample_from_col, ignore.case=TRUE)
                } else {
                    sample_name <- sample_from_col
                }
            }
        }
        
        theta_watterson_cols <- grep("\\.theta_watterson$", colnames(data), ignore.case=FALSE, value=TRUE)
        if (length(theta_watterson_cols) > 0) {
            theta_watterson_col <- theta_watterson_cols[1]
            data <- data %>%
                mutate(theta = as.numeric(.data[[theta_watterson_col]]))
        }
        
        tajima_cols <- grep("\\.tajimas?_d$", colnames(data), ignore.case=FALSE, value=TRUE)
        if (length(tajima_cols) > 0) {
            tajima_col <- tajima_cols[1]
            data <- data %>%
                mutate(tajima_d = as.numeric(.data[[tajima_col]]))
        } else if ("tajima-d" %in% colnames(data)) {
            data <- data %>% mutate(tajima_d = as.numeric(`tajima-d`))
        } else if ("tajima_d" %in% colnames(data)) {
            data <- data %>% mutate(tajima_d = as.numeric(tajima_d))
        } else if ("tajimas_d" %in% colnames(data)) {
            data <- data %>% mutate(tajima_d = as.numeric(tajimas_d))
        }
        
        data <- data %>% mutate(sample = sample_name)
        all_data[[file]] <- data
    }
    
    if (length(all_data) == 0) {
        return(NULL)
    }
    
    # Combine all data
    combined_data <- bind_rows(all_data)
    
    # Filter out rows with missing chr or pos
    combined_data <- combined_data %>%
        filter(!is.na(.data$chr), !is.na(.data$pos))
    
    return(combined_data)
}

# Function to read and process FST files
read_fst_files <- function(input_dir, selected_window_size=NULL, selected_pairs=NULL) {
    if (is.null(input_dir)) {
        return(list())
    }
    
    # Find all FST files
    tsv_files <- list.files(input_dir, pattern=".*fst.*\\.tsv$", full.names=TRUE, recursive=TRUE)
    csv_files <- list.files(input_dir, pattern=".*fst.*\\.csv$", full.names=TRUE, recursive=TRUE)
    all_files <- c(tsv_files, csv_files)
    
    if (length(all_files) == 0) {
        if (opts$verbose) cat("No FST files found in:", input_dir, "\n")
        return(list())
    }
    
    if (opts$verbose) {
        cat("Found", length(all_files), "FST file(s) (", length(tsv_files), " TSV, ", length(csv_files), " CSV)\n", sep="")
    }
    
    # Filter files by window size if specified
    if (!is.null(selected_window_size)) {
        file_window_sizes <- sapply(all_files, function(f) extract_window_size_from_filename(basename(f)))
        all_files <- all_files[is.na(file_window_sizes) | file_window_sizes == selected_window_size]
        if (opts$verbose) {
            cat("Filtered to", length(all_files), "file(s) matching window size", selected_window_size, "\n")
        }
        if (length(all_files) == 0) {
            if (opts$verbose) cat("No files found matching window size", selected_window_size, "\n")
            return(list())
        }
    }
    
    # Filter files by sample pairs if specified (peek at headers)
    selected_pairs_normalized <- NULL
    if (!is.null(selected_pairs) && length(selected_pairs) > 0) {
        selected_pairs_normalized <- sapply(selected_pairs, function(p) {
            parts <- strsplit(p, ":")[[1]]
            if (length(parts) == 2) {
                sample1 <- normalize_sample_name(parts[1])
                sample2 <- normalize_sample_name(parts[2])
                return(paste(sample1, sample2, sep=":"))
            }
            return(p)
        })
        
        if (opts$verbose) {
            cat("Filtering files by sample pairs:", paste(selected_pairs_normalized, collapse=", "), "\n")
        }
        
        # Check each file for matching pairs
        files_to_keep <- c()
        for (file in all_files) {
            cols <- peek_column_names(file)
            fst_cols <- grep("\\.fst$", cols, ignore.case=FALSE, value=TRUE)
            
            has_pair <- FALSE
            for (col in fst_cols) {
                # Try to parse pair from column name
                pair_match <- sub("\\.fst$", "", col, ignore.case=TRUE)
                if (grepl(":", pair_match)) {
                    pair_parts <- strsplit(pair_match, ":")[[1]]
                    if (length(pair_parts) == 2) {
                        sample1 <- normalize_sample_name(pair_parts[1])
                        sample2 <- normalize_sample_name(pair_parts[2])
                        pair_name <- paste(sample1, sample2, sep=":")
                        if (pair_name %in% selected_pairs_normalized) {
                            has_pair <- TRUE
                            break
                        }
                    }
                }
            }
            
            if (has_pair) {
                files_to_keep <- c(files_to_keep, file)
            }
        }
        
        all_files <- files_to_keep
        if (opts$verbose) {
            cat("Filtered to", length(all_files), "file(s) containing requested sample pairs\n")
        }
        if (length(all_files) == 0) {
            if (opts$verbose) cat("No files found containing requested sample pairs\n")
            return(list())
        }
    }
    
    # Read and combine all FST files
    all_data <- list()
    for (file in all_files) {
        if (opts$verbose) {
            cat("\n=== Reading:", basename(file), "===\n")
        }
        
        # Detect separator
        first_line <- readLines(file, n=1)
        has_tabs <- grepl("\t", first_line)
        has_commas <- grepl(",", first_line)
        
        if (has_tabs && (!has_commas || length(strsplit(first_line, "\t")[[1]]) > length(strsplit(first_line, ",")[[1]]))) {
            data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        } else if (has_commas) {
            data <- read_csv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        } else {
            data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        }
        
        # Standardize column names to lowercase
        data <- data %>%
            rename_all(tolower)
        
        # Convert base numeric columns
        data <- data %>%
            mutate(across(any_of(c("start", "end", "position", "pos")), as.numeric))
        
        # Convert FST columns to numeric
        fst_stat_cols <- grep("\\.fst$", colnames(data), ignore.case=FALSE, value=TRUE)
        if (length(fst_stat_cols) > 0) {
            data <- data %>%
                mutate(across(all_of(fst_stat_cols), as.numeric))
        }
        
        # Handle chromosome column
        if (!"chr" %in% colnames(data)) {
            if ("chromosome" %in% colnames(data)) {
                data <- data %>% mutate(chr = as.character(chromosome))
            } else if ("chrom" %in% colnames(data)) {
                data <- data %>% mutate(chr = as.character(chrom))
            } else if ("contig" %in% colnames(data)) {
                data <- data %>% mutate(chr = as.character(contig))
            } else {
                chr_col <- grep("^(chr|chrom|contig|scaffold)", colnames(data), ignore.case=FALSE, value=TRUE)
                if (length(chr_col) > 0) {
                    data <- data %>% mutate(chr = as.character(.data[[chr_col[1]]]))
                } else {
                    cat("  ERROR: Could not find chromosome column\n")
                    next
                }
            }
        } else {
            data <- data %>% mutate(chr = as.character(chr))
        }
        
        # Handle position column
        if (!"pos" %in% colnames(data)) {
            if ("start" %in% colnames(data)) {
                data <- data %>% mutate(pos = as.numeric(start))
            } else if ("position" %in% colnames(data)) {
                data <- data %>% mutate(pos = as.numeric(position))
            } else if ("end" %in% colnames(data)) {
                data <- data %>% mutate(pos = as.numeric(end))
            } else {
                pos_col <- grep("^(pos|start|end|coordinate)", colnames(data), ignore.case=FALSE, value=TRUE)
                if (length(pos_col) > 0) {
                    data <- data %>% mutate(pos = as.numeric(.data[[pos_col[1]]]))
                } else {
                    cat("  ERROR: Could not find position column\n")
                    next
                }
            }
        } else {
            data <- data %>% mutate(pos = as.numeric(pos))
        }
        
        # Extract window size from filename
        filename <- basename(file)
        window_match <- regmatches(filename, regexpr("w(\\d+)_s(\\d+)", filename))
        if (length(window_match) > 0) {
            window_parts <- strsplit(window_match, "_")[[1]]
            window_size <- as.numeric(sub("w", "", window_parts[1]))
            step_size <- as.numeric(sub("s", "", window_parts[2]))
            data$window_size <- window_size
            data$step_size <- step_size
        } else {
            data$window_size <- NA
            data$step_size <- NA
        }
        
        all_data[[file]] <- data
    }
    
    if (length(all_data) == 0) {
        return(list())
    }
    
    # Combine all data
    combined_data <- bind_rows(all_data)
    
    # Filter out rows with missing chr or pos
    combined_data <- combined_data %>%
        filter(!is.na(.data$chr), !is.na(.data$pos))
    
    # Parse FST pairs from column names
    fst_cols <- grep("\\.fst$", colnames(combined_data), ignore.case=FALSE, value=TRUE)
    if (length(fst_cols) == 0) {
        if (opts$verbose) cat("No FST columns found in combined data\n")
        return(list())
    }
    
    # First pass: collect all mappings (normalized pair -> list of columns)
    fst_pairs_raw <- list()
    for (col in fst_cols) {
        if (grepl(":", col)) {
            pair_match <- sub("\\.fst$", "", col, ignore.case=TRUE)
            pair_parts <- strsplit(pair_match, ":")[[1]]
            if (length(pair_parts) == 2) {
                sample1_raw <- pair_parts[1]
                sample2_raw <- pair_parts[2]
                sample1 <- normalize_sample_name(sample1_raw)
                sample2 <- normalize_sample_name(sample2_raw)
                pair_name <- paste(sample1, sample2, sep=":")
                if (!pair_name %in% names(fst_pairs_raw)) {
                    fst_pairs_raw[[pair_name]] <- list()
                }
                fst_pairs_raw[[pair_name]] <- c(fst_pairs_raw[[pair_name]], col)
            }
        } else {
            pair_match <- sub("\\.fst$", "", col, ignore.case=TRUE)
            pair_match <- sub("_fst$", "", pair_match, ignore.case=TRUE)
            parts <- strsplit(pair_match, "_")[[1]]
            if (length(parts) >= 2) {
                sample1_raw <- paste(parts[1:(length(parts)-1)], collapse="_")
                sample2_raw <- parts[length(parts)]
                sample1 <- normalize_sample_name(sample1_raw)
                sample2 <- normalize_sample_name(sample2_raw)
                pair_name <- paste(sample1, sample2, sep=":")
                if (!pair_name %in% names(fst_pairs_raw)) {
                    fst_pairs_raw[[pair_name]] <- list()
                }
                fst_pairs_raw[[pair_name]] <- c(fst_pairs_raw[[pair_name]], col)
            }
        }
    }
    
    # Second pass: for each normalized pair, select the best column
    fst_pairs <- list()
    for (pair_name in names(fst_pairs_raw)) {
        cols_for_pair <- fst_pairs_raw[[pair_name]]
        if (length(cols_for_pair) == 1) {
            fst_pairs[[pair_name]] <- cols_for_pair[1]
        } else {
            best_col <- NULL
            best_count <- -1
            for (col in cols_for_pair) {
                if (col %in% colnames(combined_data)) {
                    non_na_count <- sum(!is.na(combined_data[[col]]))
                    if (non_na_count > best_count) {
                        best_count <- non_na_count
                        best_col <- col
                    }
                }
            }
            if (!is.null(best_col)) {
                fst_pairs[[pair_name]] <- best_col
            } else {
                fst_pairs[[pair_name]] <- cols_for_pair[1]
            }
        }
    }
    
    if (length(fst_pairs) == 0) {
        if (opts$verbose) cat("Could not parse sample pairs from FST column names\n")
        return(list())
    }
    
    # Filter pairs if specified
    if (!is.null(selected_pairs_normalized) && length(selected_pairs_normalized) > 0) {
        fst_pairs <- fst_pairs[names(fst_pairs) %in% selected_pairs_normalized]
        if (length(fst_pairs) == 0) {
            if (opts$verbose) cat("No pairs found matching selected pairs\n")
            return(list())
        }
    }
    
    return(list(data=combined_data, pairs=fst_pairs))
}

# Read data
diversity_data <- NULL
fst_data_list <- list()

diversity_data <- NULL
fst_data_list <- list(data = NULL, pairs = character(0))

pbe_data <- NULL
if (!is.null(opts$`hdf5-dir`) && dir.exists(opts$`hdf5-dir`)) {
    if (opts$verbose) cat("Reading from HDF5 directory:", opts$`hdf5-dir`, "\n")
    diversity_data <- read_diversity_from_h5(opts$`hdf5-dir`, opts$`window-size`)
    selected_pairs <- NULL
    if (opts$`sample-pairs` != "") {
        selected_pairs <- strsplit(opts$`sample-pairs`, ",")[[1]]
        selected_pairs <- trimws(selected_pairs)
    }
    if ("fst" %in% statistics) {
        fst_data_list <- read_fst_from_h5(opts$`hdf5-dir`, opts$`window-size`, selected_pairs)
    }
    if ("pbe" %in% statistics) {
        pbe_data <- read_pbe_from_h5(opts$`hdf5-dir`, opts$`window-size`)
        if (!is.null(pbe_data) && nrow(pbe_data) > 0 && !is.null(diversity_data)) {
            sum_exprs <- list()
            if ("mean_coverage" %in% colnames(diversity_data)) sum_exprs$mean_coverage <- quote(mean(mean_coverage, na.rm = TRUE))
            if ("mean_mapping_quality" %in% colnames(diversity_data)) sum_exprs$mean_mapping_quality <- quote(mean(mean_mapping_quality, na.rm = TRUE))
            if ("n_snps" %in% colnames(diversity_data)) sum_exprs$n_snps <- quote(first(n_snps, default = NA_real_))
            if (length(sum_exprs) > 0) {
                div_win <- diversity_data %>%
                    group_by(.data$chr, .data$start, .data$end) %>%
                    summarise(!!!sum_exprs, .groups = "drop")
                pbe_data <- pbe_data %>%
                    left_join(div_win, by = c("chr", "start", "end"))
            }
        }
    }
}

if (is.null(diversity_data) && !is.null(opts$`diversity-dir`)) {
    diversity_data <- read_diversity_files(opts$`diversity-dir`, opts$`window-size`)
}

if ("fst" %in% statistics && (length(fst_data_list) == 0 || is.null(fst_data_list$data)) && !is.null(opts$`fst-dir`)) {
    selected_pairs <- NULL
    if (opts$`sample-pairs` != "") {
        selected_pairs <- strsplit(opts$`sample-pairs`, ",")[[1]]
        selected_pairs <- trimws(selected_pairs)
    }
    fst_data_list <- read_fst_files(opts$`fst-dir`, opts$`window-size`, selected_pairs)
}

# Require at least one requested statistic to have data
has_diversity <- !is.null(diversity_data) && any(c("pi", "theta", "tajima_d") %in% statistics)
has_fst <- length(fst_data_list) > 0 && !is.null(fst_data_list$data) && length(fst_data_list$pairs) > 0 && ("fst" %in% statistics)
has_pbe <- !is.null(pbe_data) && nrow(pbe_data) > 0 && ("pbe" %in% statistics)
if (!has_diversity && !has_fst && !has_pbe) {
    stop("No data found for requested statistics. Provide diversity data for pi/theta/tajima_d, FST data for fst, and/or PBE data for pbe.")
}

# Attach depth, mapping quality, n_snps to FST (and later PBE) from diversity when available
if (!is.null(diversity_data) && length(fst_data_list) > 0 && !is.null(fst_data_list$data)) {
    sum_exprs <- list()
    if ("mean_coverage" %in% colnames(diversity_data)) sum_exprs$mean_coverage <- quote(mean(mean_coverage, na.rm = TRUE))
    if ("mean_mapping_quality" %in% colnames(diversity_data)) sum_exprs$mean_mapping_quality <- quote(mean(mean_mapping_quality, na.rm = TRUE))
    if ("n_snps" %in% colnames(diversity_data)) sum_exprs$n_snps <- quote(first(n_snps, default = NA_real_))
    if (length(sum_exprs) > 0) {
        div_win <- diversity_data %>%
            group_by(.data$chr, .data$start, .data$end) %>%
            summarise(!!!sum_exprs, .groups = "drop")
        fst_data_list$data <- fst_data_list$data %>%
            left_join(div_win, by = c("chr", "start", "end"))
    }
}

# Abort before any computation if a filter is set but required column is missing in used data
check_quality_columns_and_abort(
    list(diversity = diversity_data, fst = if (length(fst_data_list) > 0) fst_data_list$data else NULL, pbe = pbe_data),
    opts, has_diversity, has_fst, has_pbe
)

# Combine diversity and FST data for chromosome length calculation
all_data_for_chr_lengths <- list()
if (!is.null(diversity_data)) {
    all_data_for_chr_lengths <- c(all_data_for_chr_lengths, list(diversity_data))
}
if (length(fst_data_list) > 0 && !is.null(fst_data_list$data)) {
    all_data_for_chr_lengths <- c(all_data_for_chr_lengths, list(fst_data_list$data))
}

if (length(all_data_for_chr_lengths) > 0) {
    combined_for_chr <- bind_rows(all_data_for_chr_lengths)
} else {
    stop("No data available for chromosome length calculation")
}

# Get chromosome lengths
chr_lengths <- get_chr_lengths(combined_for_chr, opts$`reference-genome`)

# Filter chromosomes by length if requested
if (!is.null(opts$`top-n-chromosomes`)) {
    top_n <- as.integer(opts$`top-n-chromosomes`)
    if (top_n > 0 && top_n <= nrow(chr_lengths)) {
        chr_lengths <- chr_lengths %>% head(top_n)
        if (opts$verbose) {
            cat("\nFiltering to top", top_n, "longest chromosomes\n")
            cat("Selected chromosomes:", paste(chr_lengths$chr, collapse=", "), "\n")
        }
    } else {
        if (opts$verbose) cat("Warning: --top-n-chromosomes (", top_n, ") is invalid, using all chromosomes\n", sep="")
    }
}

if (!is.null(opts$`min-chromosome-length`)) {
    min_length <- as.numeric(opts$`min-chromosome-length`)
    if (min_length > 0) {
        chr_lengths <- chr_lengths %>% filter(.data$length >= min_length)
        if (opts$verbose) {
            cat("\nFiltering to chromosomes >= ", min_length, " bp\n", sep="")
            cat("Selected chromosomes:", paste(chr_lengths$chr, collapse=", "), "\n")
            cat("Number of chromosomes:", nrow(chr_lengths), "\n")
        }
    } else {
        if (opts$verbose) cat("Warning: --min-chromosome-length (", min_length, ") is invalid, using all chromosomes\n", sep="")
    }
}

# Restrict to single chromosome and/or interval when requested
region_start <- opts$`region-start`
region_end <- opts$`region-end`
region_filename_suffix <- ""
if (!is.null(opts$chromosome) && nchar(trimws(opts$chromosome)) > 0) {
    req_chr <- trimws(opts$chromosome)
    if (!req_chr %in% chr_lengths$chr) {
        stop("--chromosome '", req_chr, "' not found in data (or filtered out by --top-n-chromosomes/--min-chromosome-length)")
    }
    chr_lengths <- chr_lengths %>% filter(.data$chr == req_chr)
    if (!is.null(region_start) || !is.null(region_end)) {
        region_start <- if (!is.null(region_start)) as.numeric(region_start) else NA_real_
        region_end <- if (!is.null(region_end)) as.numeric(region_end) else NA_real_
        if (is.na(region_start)) region_start <- 1
        if (is.na(region_end)) region_end <- chr_lengths$length[1L]
        region_filename_suffix <- paste0("_region_", gsub("[^A-Za-z0-9]", "_", req_chr), "_", as.integer(region_start), "_", as.integer(region_end))
        if (opts$verbose) cat("Restricting to interval [", region_start, ",", region_end, "] on ", req_chr, "\n", sep = "")
    } else {
        region_filename_suffix <- paste0("_chr_", gsub("[^A-Za-z0-9]", "_", req_chr))
        if (opts$verbose) cat("Restricting to chromosome:", req_chr, "\n")
    }
}

# Filter data to selected chromosomes
selected_chrs <- chr_lengths$chr
if (!is.null(diversity_data)) {
    diversity_data <- diversity_data %>%
        filter(.data$chr %in% selected_chrs)
    if (opts$verbose) {
        cat("Diversity data after chromosome filter: ", nrow(diversity_data), " rows, chr(s): ", paste(selected_chrs, collapse = ", "), "\n", sep = "")
        if (nrow(diversity_data) > 0 && "pi" %in% colnames(diversity_data)) {
            qcol <- "pi_quantile"
            if (qcol %in% colnames(diversity_data)) {
                q <- diversity_data[[qcol]]
                q <- q[!is.na(q)]
                cat("  pi_quantile: min=", if (length(q)) min(q) else NA, ", max=", if (length(q)) max(q) else NA, ", n_non_NA=", length(q), "\n", sep = "")
                if (length(q) > 0) {
                    cat("  Windows with pi_quantile >= 0.99: ", sum(q >= 0.99), ", <= 0.01: ", sum(q <= 0.01), "\n", sep = "")
                }
            }
            cat("  Samples: ", paste(unique(diversity_data$sample), collapse = ", "), "\n", sep = "")
            if ("window_size" %in% colnames(diversity_data)) cat("  Window sizes: ", paste(sort(unique(diversity_data$window_size[!is.na(diversity_data$window_size)])), collapse = ", "), "\n", sep = "")
        }
    }
}

if (length(fst_data_list) > 0 && !is.null(fst_data_list$data)) {
    fst_data_list$data <- fst_data_list$data %>%
        filter(.data$chr %in% selected_chrs)
}
if (!is.null(pbe_data)) {
    pbe_data <- pbe_data %>% filter(.data$chr %in% selected_chrs)
}

# Further restrict to [region_start, region_end] when set (windows overlapping the interval)
if (!is.null(opts$chromosome) && (is.numeric(region_start) && is.numeric(region_end))) {
    overlap_filter <- function(df) {
        if (!"start" %in% colnames(df) || !"end" %in% colnames(df)) return(df)
        df %>% filter(.data$start <= region_end, .data$end >= region_start)
    }
    if (!is.null(diversity_data)) diversity_data <- overlap_filter(diversity_data)
    if (length(fst_data_list) > 0 && !is.null(fst_data_list$data)) fst_data_list$data <- overlap_filter(fst_data_list$data)
    if (!is.null(pbe_data)) pbe_data <- overlap_filter(pbe_data)
}

# Whether any depth/mapq/SNP quality options are set (for filtering seeds/expansion)
has_quality_opts <- !all(vapply(list(opts$`min-depth`, opts$`max-depth`, opts$`min-mapping-quality`, opts$`max-mapping-quality`, opts$`min-snps`, opts$`max-snps`), is.null, logical(1L)))

# Now process statistics and identify outliers
outlier_results <- list()
expanded_region_results <- list()  # when expand_mode, holds expanded regions per key

# Process diversity statistics
if (is.null(diversity_data)) {
    if (opts$verbose) cat("No diversity data (diversity_data is NULL); skipping diversity statistics.\n")
}
if (!is.null(diversity_data)) {
    for (stat in statistics) {
        stat_col <- switch(stat,
            "pi" = "pi",
            "theta" = "theta",
            "tajima_d" = "tajima_d"
        )
        
        if (!stat_col %in% colnames(diversity_data)) {
            if (opts$verbose) cat("Warning: Column '", stat_col, "' not found in diversity data, skipping\n", sep="")
            next
        }
        
        if (all(is.na(diversity_data[[stat_col]]))) {
            if (opts$verbose) cat("Warning: Column '", stat_col, "' contains only NA values, skipping\n", sep="")
            next
        }
        
        # Filter out NA values for this statistic
        stat_data <- diversity_data %>%
            filter(!is.na(.data[[stat_col]]))
        
        if (nrow(stat_data) == 0) {
            if (opts$verbose) cat("Warning: No valid data for statistic '", stat, "', skipping\n", sep="")
            next
        }
        if (opts$verbose) cat("Statistic '", stat, "': ", nrow(stat_data), " rows, samples: ", paste(unique(stat_data$sample), collapse = ", "), "\n", sep = "")
        
        # Identify outlier windows
        # Calculate quantile thresholds per sample and per window size
        for (sample_name in unique(stat_data$sample)) {
            sample_data <- stat_data %>%
                filter(.data$sample == sample_name)
            
            if (nrow(sample_data) == 0) next
            
            # Group by window size if available
            if ("window_size" %in% colnames(sample_data) && !all(is.na(sample_data$window_size))) {
                window_sizes <- unique(sample_data$window_size[!is.na(sample_data$window_size)])
            } else {
                window_sizes <- c(NA)
            }
            
            all_high_outliers <- list()
            all_low_outliers <- list()
            all_tagged <- list()
            all_expanded_regions <- list()
            
            for (ws in window_sizes) {
                if (is.na(ws)) {
                    ws_data <- sample_data %>% filter(is.na(.data$window_size))
                } else {
                    ws_data <- sample_data %>% filter(.data$window_size == ws)
                }
                
                if (nrow(ws_data) == 0) {
                    if (opts$verbose) cat("  [", stat, "] sample=", sample_name, " window_size=", ws, " n=0 (skipping)\n", sep = "")
                    next
                }
                
                merge_dist_ws <- opts$`merge-distance`
                if (!is.null(merge_dist_ws) && (merge_dist_ws == "auto" || merge_dist_ws == "0")) {
                    merge_dist_ws <- if (is.na(ws)) 2000 else max(2 * ws, 2000)
                } else if (!is.null(merge_dist_ws)) {
                    merge_dist_ws <- as.numeric(merge_dist_ws)
                    if (is.na(merge_dist_ws) || merge_dist_ws < 0) merge_dist_ws <- 2000
                }
                
                quantile_col <- paste0(stat_col, "_quantile")
                use_quantile <- quantile_col %in% colnames(ws_data) && !all(is.na(ws_data[[quantile_col]]))
                if (opts$verbose) {
                    cat("  [", stat, "] sample=", sample_name, " window_size=", ws, " n=", nrow(ws_data), " use_quantile_column=", use_quantile, sep = "")
                    if (use_quantile && quantile_col %in% colnames(ws_data)) {
                        qv <- ws_data[[quantile_col]]
                        qv <- qv[!is.na(qv)]
                        if (length(qv) > 0) cat(" quantile_range=[", min(qv), ",", max(qv), "]", sep = "")
                    }
                    cat("\n")
                }
                
                if (expand_mode) {
                    if (use_quantile) {
                        seed_high <- ws_data %>% filter(.data[[quantile_col]] >= opts$`seed-high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`seed-high-quantile`)
                        seed_low <- ws_data %>% filter(.data[[quantile_col]] <= opts$`seed-low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`seed-low-quantile`)
                        expand_high <- ws_data %>% filter(.data[[quantile_col]] >= opts$`expand-high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`expand-high-quantile`)
                        expand_low <- ws_data %>% filter(.data[[quantile_col]] <= opts$`expand-low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`expand-low-quantile`)
                        if (opts$verbose) cat("    seed_high(>=", opts$`seed-high-quantile`, ")=", nrow(seed_high), " seed_low(<=", opts$`seed-low-quantile`, ")=", nrow(seed_low), " expand_high(>=", opts$`expand-high-quantile`, ")=", nrow(expand_high), " expand_low(<=", opts$`expand-low-quantile`, ")=", nrow(expand_low), "\n", sep = "")
                    } else {
                        sh <- quantile(ws_data[[stat_col]], opts$`seed-high-quantile`, na.rm = TRUE)
                        sl <- quantile(ws_data[[stat_col]], opts$`seed-low-quantile`, na.rm = TRUE)
                        eh <- quantile(ws_data[[stat_col]], opts$`expand-high-quantile`, na.rm = TRUE)
                        el <- quantile(ws_data[[stat_col]], opts$`expand-low-quantile`, na.rm = TRUE)
                        if (opts$verbose) cat("    Using value thresholds (no quantile column): seed_high>=", sh, " seed_low<=", sl, " expand_high>=", eh, " expand_low<=", el, "\n", sep = "")
                        seed_high <- ws_data %>% filter(.data[[stat_col]] >= sh) %>% mutate(quantile_type = "high", quantile_value = opts$`seed-high-quantile`)
                        seed_low <- ws_data %>% filter(.data[[stat_col]] <= sl) %>% mutate(quantile_type = "low", quantile_value = opts$`seed-low-quantile`)
                        expand_high <- ws_data %>% filter(.data[[stat_col]] >= eh) %>% mutate(quantile_type = "high", quantile_value = opts$`expand-high-quantile`)
                        expand_low <- ws_data %>% filter(.data[[stat_col]] <= el) %>% mutate(quantile_type = "low", quantile_value = opts$`expand-low-quantile`)
                    }
                    if (has_quality_opts) {
                        seed_high <- seed_high %>% filter(passes_quality_filters(seed_high, opts))
                        seed_low <- seed_low %>% filter(passes_quality_filters(seed_low, opts))
                        expand_high <- expand_high %>% filter(passes_quality_filters(expand_high, opts))
                        expand_low <- expand_low %>% filter(passes_quality_filters(expand_low, opts))
                        if (opts$verbose) cat("    After quality filters: seed_high=", nrow(seed_high), " seed_low=", nrow(seed_low), " expand_high=", nrow(expand_high), " expand_low=", nrow(expand_low), "\n", sep = "")
                    }
                    if (!is.null(opts$`top-n-extreme`) && opts$`top-n-extreme` > 0) {
                        n_extreme <- as.integer(opts$`top-n-extreme`)
                        if (nrow(seed_high) > n_extreme) seed_high <- seed_high %>% arrange(desc(.data[[stat_col]])) %>% head(n_extreme)
                        if (nrow(seed_low) > n_extreme) seed_low <- seed_low %>% arrange(.data[[stat_col]]) %>% head(n_extreme)
                    }
                    seed_windows <- bind_rows(if (nrow(seed_high) > 0) seed_high else NULL, if (nrow(seed_low) > 0) seed_low else NULL)
                    expand_windows <- bind_rows(if (nrow(expand_high) > 0) expand_high else NULL, if (nrow(expand_low) > 0) expand_low else NULL)
                    if (opts$verbose) cat("    seed_windows=", nrow(seed_windows), " expand_windows=", nrow(expand_windows), " merge_dist_ws=", merge_dist_ws, "\n", sep = "")
                    if (nrow(seed_windows) > 0 && !is.null(merge_dist_ws) && merge_dist_ws > 0) {
                        res <- seed_expand_regions(seed_windows, expand_windows, merge_dist_ws, verbose = opts$verbose)
                        if (nrow(res$windows) > 0) {
                            all_tagged[[length(all_tagged) + 1L]] <- res$windows
                            if (!is.null(res$regions) && nrow(res$regions) > 0) {
                                res$regions$statistic <- stat
                                res$regions$sample <- sample_name
                                res$regions$sample_pair <- NA_character_
                                all_expanded_regions[[length(all_expanded_regions) + 1L]] <- res$regions
                            }
                        }
                    } else if (nrow(seed_windows) > 0) {
                        all_tagged[[length(all_tagged) + 1L]] <- seed_windows %>% mutate(window_type = "seed")
                    }
                } else {
                    if (use_quantile) {
                        high_outliers <- ws_data %>%
                            filter(.data[[quantile_col]] >= opts$`high-quantile`) %>%
                            mutate(quantile_type = "high", quantile_value = opts$`high-quantile`)
                        low_outliers <- ws_data %>%
                            filter(.data[[quantile_col]] <= opts$`low-quantile`) %>%
                            mutate(quantile_type = "low", quantile_value = opts$`low-quantile`)
                    } else {
                        high_threshold <- quantile(ws_data[[stat_col]], opts$`high-quantile`, na.rm=TRUE)
                        low_threshold <- quantile(ws_data[[stat_col]], opts$`low-quantile`, na.rm=TRUE)
                        high_outliers <- ws_data %>% filter(.data[[stat_col]] >= high_threshold) %>% mutate(quantile_type = "high", quantile_value = opts$`high-quantile`)
                        low_outliers <- ws_data %>% filter(.data[[stat_col]] <= low_threshold) %>% mutate(quantile_type = "low", quantile_value = opts$`low-quantile`)
                    }
                    if (has_quality_opts) {
                        high_outliers <- high_outliers %>% filter(passes_quality_filters(high_outliers, opts))
                        low_outliers <- low_outliers %>% filter(passes_quality_filters(low_outliers, opts))
                    }
                    if (!is.null(opts$`top-n-extreme`) && opts$`top-n-extreme` > 0) {
                        n_extreme <- as.integer(opts$`top-n-extreme`)
                        if (nrow(high_outliers) > n_extreme) high_outliers <- high_outliers %>% arrange(desc(.data[[stat_col]])) %>% head(n_extreme)
                        if (nrow(low_outliers) > n_extreme) low_outliers <- low_outliers %>% arrange(.data[[stat_col]]) %>% head(n_extreme)
                    }
                    if (nrow(high_outliers) > 0) all_high_outliers[[length(all_high_outliers) + 1]] <- high_outliers
                    if (nrow(low_outliers) > 0) all_low_outliers[[length(all_low_outliers) + 1]] <- low_outliers
                }
            }
            
            if (expand_mode && length(all_tagged) > 0) {
                outliers <- bind_rows(all_tagged) %>%
                    mutate(statistic = stat, sample_pair = NA_character_)
                if (!"start" %in% colnames(outliers)) outliers$start <- NA_real_
                if (!"end" %in% colnames(outliers)) outliers$end <- NA_real_
                out_cols <- c("statistic", "sample", "sample_pair", "window_size", "chr", "start", "end", "pos", "value", "quantile", "quantile_type", "quantile_value", "window_type")
                extra_cols <- intersect(c("mean_coverage", "mean_mapping_quality", "n_snps"), colnames(outliers))
                outliers <- standardize_outlier_columns(outliers, stat_col, paste0(stat_col, "_quantile"), out_cols, extra_cols)
                outlier_results[[paste0("diversity_", stat, "_", sample_name)]] <- outliers
                if (length(all_expanded_regions) > 0) {
                    expanded_region_results[[paste0("diversity_", stat, "_", sample_name)]] <- bind_rows(all_expanded_regions)
                }
            } else if (length(all_high_outliers) > 0 || length(all_low_outliers) > 0) {
                outliers <- bind_rows(
                    if (length(all_high_outliers) > 0) bind_rows(all_high_outliers) else NULL,
                    if (length(all_low_outliers) > 0) bind_rows(all_low_outliers) else NULL
                ) %>%
                    mutate(statistic = stat, sample_pair = NA_character_)
                if (!"start" %in% colnames(outliers)) outliers$start <- NA_real_
                if (!"end" %in% colnames(outliers)) outliers$end <- NA_real_
                out_cols <- c("statistic", "sample", "sample_pair", "window_size", "chr", "start", "end", "pos", "value", "quantile", "quantile_type", "quantile_value")
                extra_cols <- intersect(c("mean_coverage", "mean_mapping_quality", "n_snps"), colnames(outliers))
                outliers <- standardize_outlier_columns(outliers, stat_col, paste0(stat_col, "_quantile"), out_cols, extra_cols)
                outlier_results[[paste0("diversity_", stat, "_", sample_name)]] <- outliers
            }
        }
    }
}

# Process FST statistics
if (length(fst_data_list) > 0 && !is.null(fst_data_list$data) && length(fst_data_list$pairs) > 0) {
    for (pair_name in names(fst_data_list$pairs)) {
        fst_col <- fst_data_list$pairs[[pair_name]]
        
        if (!fst_col %in% colnames(fst_data_list$data)) {
            if (opts$verbose) cat("Warning: Column '", fst_col, "' not found, skipping pair ", pair_name, "\n", sep="")
            next
        }
        
        # Extract FST values for this pair (filter by sample_pair when data is long-format, e.g. from HDF5)
        pair_data <- fst_data_list$data %>%
            filter(!is.na(.data[[fst_col]]))
        if ("sample_pair" %in% colnames(fst_data_list$data)) {
            pair_data <- pair_data %>% filter(.data$sample_pair == pair_name)
        }
        pair_data <- pair_data %>%
            mutate(
                fst = .data[[fst_col]],
                sample_pair = pair_name
            )
        
        if (nrow(pair_data) == 0) {
            if (opts$verbose) cat("Warning: No valid FST data for pair '", pair_name, "', skipping\n", sep="")
            next
        }
        
        # Identify outlier windows
        if ("window_size" %in% colnames(pair_data) && !all(is.na(pair_data$window_size))) {
            window_sizes <- unique(pair_data$window_size[!is.na(pair_data$window_size)])
        } else {
            window_sizes <- c(NA)
        }
        
        all_high_outliers <- list()
        all_low_outliers <- list()
        all_tagged_fst <- list()
        all_expanded_regions_fst <- list()
        
        for (ws in window_sizes) {
            if (is.na(ws)) {
                ws_data <- pair_data %>% filter(is.na(.data$window_size))
            } else {
                ws_data <- pair_data %>% filter(.data$window_size == ws)
            }
            if (nrow(ws_data) == 0) next
            
            merge_dist_ws <- opts$`merge-distance`
            if (!is.null(merge_dist_ws) && (merge_dist_ws == "auto" || merge_dist_ws == "0")) {
                merge_dist_ws <- if (is.na(ws)) 2000 else max(2 * ws, 2000)
            } else if (!is.null(merge_dist_ws)) {
                merge_dist_ws <- as.numeric(merge_dist_ws)
                if (is.na(merge_dist_ws) || merge_dist_ws < 0) merge_dist_ws <- 2000
            }
            
            use_fst_quantile <- "fst_quantile" %in% colnames(ws_data) && !all(is.na(ws_data$fst_quantile))
            
            if (expand_mode) {
                if (use_fst_quantile) {
                    seed_high <- ws_data %>% filter(.data$fst_quantile >= opts$`seed-high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`seed-high-quantile`)
                    seed_low <- ws_data %>% filter(.data$fst_quantile <= opts$`seed-low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`seed-low-quantile`)
                    expand_high <- ws_data %>% filter(.data$fst_quantile >= opts$`expand-high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`expand-high-quantile`)
                    expand_low <- ws_data %>% filter(.data$fst_quantile <= opts$`expand-low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`expand-low-quantile`)
                } else {
                    sh <- quantile(ws_data$fst, opts$`seed-high-quantile`, na.rm = TRUE)
                    sl <- quantile(ws_data$fst, opts$`seed-low-quantile`, na.rm = TRUE)
                    eh <- quantile(ws_data$fst, opts$`expand-high-quantile`, na.rm = TRUE)
                    el <- quantile(ws_data$fst, opts$`expand-low-quantile`, na.rm = TRUE)
                    seed_high <- ws_data %>% filter(.data$fst >= sh) %>% mutate(quantile_type = "high", quantile_value = opts$`seed-high-quantile`)
                    seed_low <- ws_data %>% filter(.data$fst <= sl) %>% mutate(quantile_type = "low", quantile_value = opts$`seed-low-quantile`)
                    expand_high <- ws_data %>% filter(.data$fst >= eh) %>% mutate(quantile_type = "high", quantile_value = opts$`expand-high-quantile`)
                    expand_low <- ws_data %>% filter(.data$fst <= el) %>% mutate(quantile_type = "low", quantile_value = opts$`expand-low-quantile`)
                }
                if (has_quality_opts) {
                    seed_high <- seed_high %>% filter(passes_quality_filters(seed_high, opts))
                    seed_low <- seed_low %>% filter(passes_quality_filters(seed_low, opts))
                    expand_high <- expand_high %>% filter(passes_quality_filters(expand_high, opts))
                    expand_low <- expand_low %>% filter(passes_quality_filters(expand_low, opts))
                }
                if (!is.null(opts$`top-n-extreme`) && opts$`top-n-extreme` > 0) {
                    n_extreme <- as.integer(opts$`top-n-extreme`)
                    if (nrow(seed_high) > n_extreme) seed_high <- seed_high %>% arrange(desc(.data$fst)) %>% head(n_extreme)
                    if (nrow(seed_low) > n_extreme) seed_low <- seed_low %>% arrange(.data$fst) %>% head(n_extreme)
                }
                seed_windows <- bind_rows(if (nrow(seed_high) > 0) seed_high else NULL, if (nrow(seed_low) > 0) seed_low else NULL)
                expand_windows <- bind_rows(if (nrow(expand_high) > 0) expand_high else NULL, if (nrow(expand_low) > 0) expand_low else NULL)
                if (nrow(seed_windows) > 0 && !is.null(merge_dist_ws) && merge_dist_ws > 0) {
                    res <- seed_expand_regions(seed_windows, expand_windows, merge_dist_ws, verbose = opts$verbose)
                    if (nrow(res$windows) > 0) {
                        all_tagged_fst[[length(all_tagged_fst) + 1L]] <- res$windows
                        if (!is.null(res$regions) && nrow(res$regions) > 0) {
                            res$regions$statistic <- "fst"
                            res$regions$sample <- NA_character_
                            res$regions$sample_pair <- pair_name
                            all_expanded_regions_fst[[length(all_expanded_regions_fst) + 1L]] <- res$regions
                        }
                    }
                } else if (nrow(seed_windows) > 0) {
                    all_tagged_fst[[length(all_tagged_fst) + 1L]] <- seed_windows %>% mutate(window_type = "seed")
                }
            } else {
                if (use_fst_quantile) {
                    high_outliers <- ws_data %>% filter(.data$fst_quantile >= opts$`high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`high-quantile`)
                    low_outliers <- ws_data %>% filter(.data$fst_quantile <= opts$`low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`low-quantile`)
                } else {
                    high_threshold <- quantile(ws_data$fst, opts$`high-quantile`, na.rm=TRUE)
                    low_threshold <- quantile(ws_data$fst, opts$`low-quantile`, na.rm=TRUE)
                    high_outliers <- ws_data %>% filter(.data$fst >= high_threshold) %>% mutate(quantile_type = "high", quantile_value = opts$`high-quantile`)
                    low_outliers <- ws_data %>% filter(.data$fst <= low_threshold) %>% mutate(quantile_type = "low", quantile_value = opts$`low-quantile`)
                }
                if (has_quality_opts) {
                    high_outliers <- high_outliers %>% filter(passes_quality_filters(high_outliers, opts))
                    low_outliers <- low_outliers %>% filter(passes_quality_filters(low_outliers, opts))
                }
                if (!is.null(opts$`top-n-extreme`) && opts$`top-n-extreme` > 0) {
                    n_extreme <- as.integer(opts$`top-n-extreme`)
                    if (nrow(high_outliers) > n_extreme) high_outliers <- high_outliers %>% arrange(desc(.data$fst)) %>% head(n_extreme)
                    if (nrow(low_outliers) > n_extreme) low_outliers <- low_outliers %>% arrange(.data$fst) %>% head(n_extreme)
                }
                if (nrow(high_outliers) > 0) all_high_outliers[[length(all_high_outliers) + 1]] <- high_outliers
                if (nrow(low_outliers) > 0) all_low_outliers[[length(all_low_outliers) + 1]] <- low_outliers
            }
        }
        
        if (expand_mode && length(all_tagged_fst) > 0) {
            outliers <- bind_rows(all_tagged_fst) %>% mutate(statistic = "fst", sample = NA_character_)
            if (!"start" %in% colnames(outliers)) outliers$start <- NA_real_
            if (!"end" %in% colnames(outliers)) outliers$end <- NA_real_
            out_cols_fst <- c("statistic", "sample", "sample_pair", "window_size", "chr", "start", "end", "pos", "value", "quantile", "quantile_type", "quantile_value", "window_type")
            extra_cols_fst <- intersect(c("mean_coverage", "mean_mapping_quality", "n_snps"), colnames(outliers))
            outliers <- standardize_outlier_columns(outliers, "fst", "fst_quantile", out_cols_fst, extra_cols_fst)
            outlier_results[[paste0("fst_", pair_name)]] <- outliers
            if (length(all_expanded_regions_fst) > 0) {
                expanded_region_results[[paste0("fst_", pair_name)]] <- bind_rows(all_expanded_regions_fst)
            }
        } else if (length(all_high_outliers) > 0 || length(all_low_outliers) > 0) {
            outliers <- bind_rows(
                if (length(all_high_outliers) > 0) bind_rows(all_high_outliers) else NULL,
                if (length(all_low_outliers) > 0) bind_rows(all_low_outliers) else NULL
            ) %>% mutate(statistic = "fst", sample = NA_character_)
            if (!"start" %in% colnames(outliers)) outliers$start <- NA_real_
            if (!"end" %in% colnames(outliers)) outliers$end <- NA_real_
            out_cols_fst <- c("statistic", "sample", "sample_pair", "window_size", "chr", "start", "end", "pos", "value", "quantile", "quantile_type", "quantile_value")
            extra_cols_fst <- intersect(c("mean_coverage", "mean_mapping_quality", "n_snps"), colnames(outliers))
            outliers <- standardize_outlier_columns(outliers, "fst", "fst_quantile", out_cols_fst, extra_cols_fst)
            outlier_results[[paste0("fst_", pair_name)]] <- outliers
        }
    }
}

# Process PBE statistics (from HDF5 only)
if (!is.null(pbe_data) && nrow(pbe_data) > 0 && "pbe" %in% statistics) {
    trios <- pbe_data %>% distinct(pop1, pop2, pop3) %>% mutate(trio_id = paste0(.data$pop1, ":", .data$pop2, ":", .data$pop3))
    for (i in seq_len(nrow(trios))) {
        trio_id <- trios$trio_id[i]
        pair_data <- pbe_data %>%
            filter(.data$pop1 == trios$pop1[i], .data$pop2 == trios$pop2[i], .data$pop3 == trios$pop3[i], !is.na(.data$pbe))
        if (nrow(pair_data) == 0) next
        if ("window_size" %in% colnames(pair_data) && !all(is.na(pair_data$window_size))) {
            window_sizes <- unique(pair_data$window_size[!is.na(pair_data$window_size)])
        } else {
            window_sizes <- c(NA)
        }
        all_high_outliers <- list()
        all_low_outliers <- list()
        all_tagged_pbe <- list()
        all_expanded_regions_pbe <- list()
        for (ws in window_sizes) {
            if (is.na(ws)) {
                ws_data <- pair_data %>% filter(is.na(.data$window_size))
            } else {
                ws_data <- pair_data %>% filter(.data$window_size == ws)
            }
            if (nrow(ws_data) == 0) next
            merge_dist_ws <- opts$`merge-distance`
            if (!is.null(merge_dist_ws) && (merge_dist_ws == "auto" || merge_dist_ws == "0")) {
                merge_dist_ws <- if (is.na(ws)) 2000 else max(2 * ws, 2000)
            } else if (!is.null(merge_dist_ws)) {
                merge_dist_ws <- as.numeric(merge_dist_ws)
                if (is.na(merge_dist_ws) || merge_dist_ws < 0) merge_dist_ws <- 2000
            }
            use_pbe_quantile <- "pbe_quantile" %in% colnames(ws_data) && !all(is.na(ws_data$pbe_quantile))
            if (expand_mode) {
                if (use_pbe_quantile) {
                    seed_high <- ws_data %>% filter(.data$pbe_quantile >= opts$`seed-high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`seed-high-quantile`)
                    seed_low <- ws_data %>% filter(.data$pbe_quantile <= opts$`seed-low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`seed-low-quantile`)
                    expand_high <- ws_data %>% filter(.data$pbe_quantile >= opts$`expand-high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`expand-high-quantile`)
                    expand_low <- ws_data %>% filter(.data$pbe_quantile <= opts$`expand-low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`expand-low-quantile`)
                } else {
                    sh <- quantile(ws_data$pbe, opts$`seed-high-quantile`, na.rm = TRUE)
                    sl <- quantile(ws_data$pbe, opts$`seed-low-quantile`, na.rm = TRUE)
                    eh <- quantile(ws_data$pbe, opts$`expand-high-quantile`, na.rm = TRUE)
                    el <- quantile(ws_data$pbe, opts$`expand-low-quantile`, na.rm = TRUE)
                    seed_high <- ws_data %>% filter(.data$pbe >= sh) %>% mutate(quantile_type = "high", quantile_value = opts$`seed-high-quantile`)
                    seed_low <- ws_data %>% filter(.data$pbe <= sl) %>% mutate(quantile_type = "low", quantile_value = opts$`seed-low-quantile`)
                    expand_high <- ws_data %>% filter(.data$pbe >= eh) %>% mutate(quantile_type = "high", quantile_value = opts$`expand-high-quantile`)
                    expand_low <- ws_data %>% filter(.data$pbe <= el) %>% mutate(quantile_type = "low", quantile_value = opts$`expand-low-quantile`)
                }
                if (has_quality_opts) {
                    seed_high <- seed_high %>% filter(passes_quality_filters(seed_high, opts))
                    seed_low <- seed_low %>% filter(passes_quality_filters(seed_low, opts))
                    expand_high <- expand_high %>% filter(passes_quality_filters(expand_high, opts))
                    expand_low <- expand_low %>% filter(passes_quality_filters(expand_low, opts))
                }
                if (!is.null(opts$`top-n-extreme`) && opts$`top-n-extreme` > 0) {
                    n_extreme <- as.integer(opts$`top-n-extreme`)
                    if (nrow(seed_high) > n_extreme) seed_high <- seed_high %>% arrange(desc(.data$pbe)) %>% head(n_extreme)
                    if (nrow(seed_low) > n_extreme) seed_low <- seed_low %>% arrange(.data$pbe) %>% head(n_extreme)
                }
                seed_windows <- bind_rows(if (nrow(seed_high) > 0) seed_high else NULL, if (nrow(seed_low) > 0) seed_low else NULL)
                expand_windows <- bind_rows(if (nrow(expand_high) > 0) expand_high else NULL, if (nrow(expand_low) > 0) expand_low else NULL)
                if (nrow(seed_windows) > 0 && !is.null(merge_dist_ws) && merge_dist_ws > 0) {
                    res <- seed_expand_regions(seed_windows, expand_windows, merge_dist_ws, verbose = opts$verbose)
                    if (nrow(res$windows) > 0) {
                        all_tagged_pbe[[length(all_tagged_pbe) + 1L]] <- res$windows
                        if (!is.null(res$regions) && nrow(res$regions) > 0) {
                            res$regions$statistic <- "pbe"
                            res$regions$sample <- NA_character_
                            res$regions$sample_pair <- trio_id
                            all_expanded_regions_pbe[[length(all_expanded_regions_pbe) + 1L]] <- res$regions
                        }
                    }
                } else if (nrow(seed_windows) > 0) {
                    all_tagged_pbe[[length(all_tagged_pbe) + 1L]] <- seed_windows %>% mutate(window_type = "seed")
                }
            } else {
                if (use_pbe_quantile) {
                    high_outliers <- ws_data %>% filter(.data$pbe_quantile >= opts$`high-quantile`) %>% mutate(quantile_type = "high", quantile_value = opts$`high-quantile`)
                    low_outliers <- ws_data %>% filter(.data$pbe_quantile <= opts$`low-quantile`) %>% mutate(quantile_type = "low", quantile_value = opts$`low-quantile`)
                } else {
                    high_threshold <- quantile(ws_data$pbe, opts$`high-quantile`, na.rm = TRUE)
                    low_threshold <- quantile(ws_data$pbe, opts$`low-quantile`, na.rm = TRUE)
                    high_outliers <- ws_data %>% filter(.data$pbe >= high_threshold) %>% mutate(quantile_type = "high", quantile_value = opts$`high-quantile`)
                    low_outliers <- ws_data %>% filter(.data$pbe <= low_threshold) %>% mutate(quantile_type = "low", quantile_value = opts$`low-quantile`)
                }
                if (has_quality_opts) {
                    high_outliers <- high_outliers %>% filter(passes_quality_filters(high_outliers, opts))
                    low_outliers <- low_outliers %>% filter(passes_quality_filters(low_outliers, opts))
                }
                if (!is.null(opts$`top-n-extreme`) && opts$`top-n-extreme` > 0) {
                    n_extreme <- as.integer(opts$`top-n-extreme`)
                    if (nrow(high_outliers) > n_extreme) high_outliers <- high_outliers %>% arrange(desc(.data$pbe)) %>% head(n_extreme)
                    if (nrow(low_outliers) > n_extreme) low_outliers <- low_outliers %>% arrange(.data$pbe) %>% head(n_extreme)
                }
                if (nrow(high_outliers) > 0) all_high_outliers[[length(all_high_outliers) + 1]] <- high_outliers
                if (nrow(low_outliers) > 0) all_low_outliers[[length(all_low_outliers) + 1]] <- low_outliers
            }
        }
        if (expand_mode && length(all_tagged_pbe) > 0) {
            outliers <- bind_rows(all_tagged_pbe) %>% mutate(statistic = "pbe", sample = NA_character_, sample_pair = trio_id)
            if (!"start" %in% colnames(outliers)) outliers$start <- NA_real_
            if (!"end" %in% colnames(outliers)) outliers$end <- NA_real_
            out_cols_pbe <- c("statistic", "sample", "sample_pair", "window_size", "chr", "start", "end", "pos", "value", "quantile", "quantile_type", "quantile_value", "window_type")
            extra_cols_pbe <- intersect(c("mean_coverage", "mean_mapping_quality", "n_snps"), colnames(outliers))
            outliers <- standardize_outlier_columns(outliers, "pbe", "pbe_quantile", out_cols_pbe, extra_cols_pbe)
            outlier_results[[paste0("pbe_", trio_id)]] <- outliers
            if (length(all_expanded_regions_pbe) > 0) {
                expanded_region_results[[paste0("pbe_", trio_id)]] <- bind_rows(all_expanded_regions_pbe)
            }
        } else if (length(all_high_outliers) > 0 || length(all_low_outliers) > 0) {
            outliers <- bind_rows(
                if (length(all_high_outliers) > 0) bind_rows(all_high_outliers) else NULL,
                if (length(all_low_outliers) > 0) bind_rows(all_low_outliers) else NULL
            ) %>% mutate(statistic = "pbe", sample = NA_character_, sample_pair = trio_id)
            if (!"start" %in% colnames(outliers)) outliers$start <- NA_real_
            if (!"end" %in% colnames(outliers)) outliers$end <- NA_real_
            out_cols_pbe <- c("statistic", "sample", "sample_pair", "window_size", "chr", "start", "end", "pos", "value", "quantile", "quantile_type", "quantile_value")
            extra_cols_pbe <- intersect(c("mean_coverage", "mean_mapping_quality", "n_snps"), colnames(outliers))
            outliers <- standardize_outlier_columns(outliers, "pbe", "pbe_quantile", out_cols_pbe, extra_cols_pbe)
            outlier_results[[paste0("pbe_", trio_id)]] <- outliers
        }
    }
}

# Combine all outlier windows and build one wide table (chr, start, end, samplename.stat, ..., outlier_stat)
if (length(outlier_results) > 0) {
    all_outliers <- bind_rows(outlier_results)
    if (!"start" %in% colnames(all_outliers)) all_outliers$start <- all_outliers$pos
    if (!"end" %in% colnames(all_outliers)) all_outliers$end <- all_outliers$pos
    all_outliers$value_col <- ifelse(
        all_outliers$statistic %in% c("pi", "theta", "tajima_d"),
        paste0(ifelse(is.na(all_outliers$sample), "unknown", all_outliers$sample), ".", all_outliers$statistic),
        ifelse(all_outliers$statistic == "fst",
               paste0(ifelse(is.na(all_outliers$sample_pair), "unknown", all_outliers$sample_pair), ".fst"),
               paste0(ifelse(is.na(all_outliers$sample_pair), "unknown", all_outliers$sample_pair), ".pbe"))
    )
    qtype_col <- if ("quantile_type" %in% colnames(all_outliers)) "quantile_type" else NULL
    stats_per_window <- all_outliers %>%
        group_by(.data$chr, .data$start, .data$end) %>%
        summarise(
            outlier_stat = paste(sort(unique(.data$statistic)), collapse = ":"),
            outlier_direction = if (!is.null(qtype_col)) {
                dir <- paste(sort(unique(.data$quantile_type[!is.na(.data$quantile_type)])), collapse = ":")
                ifelse(nchar(dir) == 0, NA_character_, dir)
            } else NA_character_,
            .groups = "drop"
        )
    wide_long <- all_outliers %>% select("chr", "start", "end", "value_col", "value")
    wide_long <- wide_long %>% distinct(chr, start, end, value_col, .keep_all = TRUE)
    wide_outliers <- wide_long %>%
        pivot_wider(names_from = "value_col", values_from = "value") %>%
        left_join(stats_per_window, by = c("chr", "start", "end"))
    if ("quantile" %in% colnames(all_outliers)) {
        wide_quantile <- all_outliers %>%
            select("chr", "start", "end", "value_col", "quantile") %>%
            distinct(chr, start, end, value_col, .keep_all = TRUE) %>%
            pivot_wider(names_from = "value_col", values_from = "quantile", names_glue = "{.name}_quantile")
        wide_outliers <- wide_outliers %>% left_join(wide_quantile, by = c("chr", "start", "end"))
    }
    window_meta <- all_outliers %>%
        group_by(.data$chr, .data$start, .data$end) %>%
        slice(1L) %>%
        ungroup() %>%
        select("chr", "start", "end", any_of(c("window_size", "mean_coverage", "mean_mapping_quality", "n_snps", "window_type")))
    wide_outliers <- wide_outliers %>% left_join(window_meta, by = c("chr", "start", "end"))
    value_cols <- sort(setdiff(names(wide_outliers), c("chr", "start", "end", "outlier_stat", "outlier_direction", "window_size", "mean_coverage", "mean_mapping_quality", "n_snps", "window_type")))
    meta_cols <- c("window_size", "mean_coverage", "mean_mapping_quality", "n_snps", "outlier_stat", "outlier_direction", "window_type")
    meta_cols <- intersect(meta_cols, names(wide_outliers))
    wide_outliers <- wide_outliers %>% select("chr", "start", "end", any_of(meta_cols), all_of(value_cols))

    unique_windows_outliers <- if ("window_size" %in% colnames(all_outliers)) {
        sort(unique(all_outliers$window_size[!is.na(all_outliers$window_size)]))
    } else {
        numeric(0)
    }
    unique_chrs_outliers <- unique(all_outliers$chr[!is.na(all_outliers$chr) & all_outliers$chr != "genome"])
    if (!is.null(opts$`window-size`)) {
        window_suffix <- paste0("_w", opts$`window-size`)
    } else if (length(unique_windows_outliers) > 0) {
        window_suffix <- paste0("_w", paste(unique_windows_outliers, collapse = "-"))
    } else {
        window_suffix <- ""
    }
    if (length(unique_chrs_outliers) > 0) {
        chr_order <- match(unique_chrs_outliers, chr_lengths$chr)
        unique_chrs_outliers <- unique_chrs_outliers[order(chr_order)]
        chr_list <- if (length(unique_chrs_outliers) <= 3) {
            paste(gsub("[^A-Za-z0-9]", "_", unique_chrs_outliers), collapse = "-")
        } else {
            paste0(gsub("[^A-Za-z0-9]", "_", unique_chrs_outliers[1]), "-and", length(unique_chrs_outliers) - 1L, "more")
        }
        chr_suffix <- paste0("_chr", chr_list)
    } else {
        chr_suffix <- "_chrAll"
    }
    if (!is.null(opts$`top-n-chromosomes`)) {
        chr_suffix <- paste0(chr_suffix, "_top", opts$`top-n-chromosomes`, "chr")
    } else if (!is.null(opts$`min-chromosome-length`)) {
        chr_suffix <- paste0(chr_suffix, "_min", opts$`min-chromosome-length`, "bp")
    }
    quantile_suffix <- paste0("_high", opts$`high-quantile`, "_low", opts$`low-quantile`)
    if (!is.null(opts$`top-n-extreme`) && opts$`top-n-extreme` > 0) {
        quantile_suffix <- paste0(quantile_suffix, "_top", opts$`top-n-extreme`)
    }
    merge_suffix <- if (opts$merge_across_samples) "_across_samples" else "_within_samples"
    value_cols_for_regions <- setdiff(names(wide_outliers), c("chr", "start", "end", "outlier_stat", "outlier_direction", "window_size", "mean_coverage", "mean_mapping_quality", "n_snps", "window_type"))
    value_cols_for_regions <- sort(value_cols_for_regions[!grepl("_quantile$", value_cols_for_regions)])
    windows_by_value_col <- build_windows_by_value_col(
        value_cols_for_regions,
        diversity_data,
        if (length(fst_data_list) > 0 && !is.null(fst_data_list$data)) fst_data_list$data else NULL,
        pbe_data
    )

    outliers_file <- file.path(opts$`output-dir`, paste0("outlier_windows", window_suffix, chr_suffix, quantile_suffix, region_filename_suffix, ".csv"))
    write_csv(wide_outliers, outliers_file)
    cat("Saved outlier windows to:", outliers_file, "\n")
    cat("Total outlier windows:", nrow(wide_outliers), "\n")

    if (expand_mode && length(expanded_region_results) > 0) {
        all_expanded <- bind_rows(expanded_region_results)
        if (opts$verbose) cat("Seed-expand: all_expanded rows =", nrow(all_expanded), "\n")
        if (opts$merge_across_samples) {
            # Merge overlapping/nearby regions from different samples into one set of regions.
            merge_dist <- opts$`merge-distance`
            if (merge_dist == "auto") {
                wins <- unique(all_outliers$window_size[!is.na(all_outliers$window_size)])
                merge_dist <- if (length(wins) > 0) max(2 * max(wins), 2000) else 2000
            } else {
                merge_dist <- as.numeric(merge_dist)
                if (is.na(merge_dist) || merge_dist < 0) merge_dist <- 2000
            }
            merged_bounds <- merge_regions(all_expanded %>% select("chr", "region_start", "region_end"), max_gap = merge_dist)
            if (opts$verbose) cat("  merged_bounds (merged regions) =", nrow(merged_bounds), "\n")
            # Reattach outlier_stat and other cols from overlapping all_expanded rows (join by chr is many-to-many; filter keeps overlapping bounds).
            overlap_join <- all_expanded %>%
                inner_join(merged_bounds %>% rename(region_start_m = "region_start", region_end_m = "region_end"), by = "chr", relationship = "many-to-many") %>%
                filter(.data$region_start <= .data$region_end_m, .data$region_end >= .data$region_start_m)
            if (opts$verbose) cat("  overlap_join rows (after many-to-many by chr + overlap filter) =", nrow(overlap_join), "\n")
            regions <- overlap_join %>%
                group_by(.data$chr, .data$region_start_m, .data$region_end_m) %>%
                summarise(outlier_stat = paste(sort(unique(.data$statistic)), collapse = ":"), .groups = "drop") %>%
                rename(region_start = "region_start_m", region_end = "region_end_m")
            other_cols <- overlap_join %>%
                group_by(.data$chr, .data$region_start_m, .data$region_end_m) %>%
                slice(1L) %>%
                select("chr", region_start = "region_start_m", region_end = "region_end_m", any_of(c("n_windows", "region_mean_coverage", "region_mean_mapping_quality", "region_n_snps")))
            regions <- regions %>% left_join(other_cols, by = c("chr", "region_start", "region_end"))
        } else {
            # One row per (chr, region_start, region_end) from any sample; no cross-sample merging.
            regions <- all_expanded %>%
                group_by(.data$chr, .data$region_start, .data$region_end) %>%
                summarise(outlier_stat = paste(sort(unique(.data$statistic)), collapse = ":"), .groups = "drop")
            other_cols <- all_expanded %>%
                group_by(.data$chr, .data$region_start, .data$region_end) %>%
                slice(1L) %>%
                select("chr", "region_start", "region_end", any_of(c("n_windows", "region_mean_coverage", "region_mean_mapping_quality", "region_n_snps")))
            regions <- regions %>% left_join(other_cols, by = c("chr", "region_start", "region_end"))
        }
        regions <- add_region_stat_summaries(regions, windows_by_value_col)
        regions <- filter_regions_by_opts(regions, opts)
        if (nrow(regions) > 0) {
            regions_file <- file.path(opts$`output-dir`, paste0("outlier_regions", window_suffix, chr_suffix, quantile_suffix, merge_suffix, "_seed_expand", region_filename_suffix, ".csv"))
            write_csv(regions, regions_file)
            cat("Saved outlier regions to:", regions_file, "\n")
            cat("Total regions:", nrow(regions), "\n")
        }
    } else {
        merge_dist <- opts$`merge-distance`
        if (!is.null(merge_dist) && (merge_dist != "0" && merge_dist != 0)) {
            if (merge_dist == "auto") {
                wins <- unique(all_outliers$window_size[!is.na(all_outliers$window_size)])
                merge_dist <- if (length(wins) > 0) max(2 * max(wins), 2000) else 2000
                if (opts$verbose) cat("Merge distance (auto):", merge_dist, "bp\n")
            } else {
                merge_dist <- as.numeric(merge_dist)
                if (is.na(merge_dist) || merge_dist < 0) merge_dist <- 2000
            }
            if (opts$merge_across_samples) {
                regions <- merge_outlier_windows_to_regions(all_outliers, merge_dist)
            } else {
                # Merge within each sample (or sample_pair for FST/PBE).
                merge_by <- if ("sample" %in% colnames(all_outliers)) "sample" else "sample_pair"
                if (!merge_by %in% colnames(all_outliers)) merge_by <- NULL
                if (!is.null(merge_by)) {
                    regions_list <- all_outliers %>%
                        mutate(merge_by = coalesce(!!sym(merge_by), "unknown")) %>%
                        group_by(.data$merge_by) %>%
                        group_map(~ merge_outlier_windows_to_regions(.x %>% select(-"merge_by"), merge_dist), .keep = TRUE)
                    regions_list <- Filter(Negate(is.null), regions_list)
                    regions <- if (length(regions_list) > 0) bind_rows(regions_list) else NULL
                } else {
                    regions <- merge_outlier_windows_to_regions(all_outliers, merge_dist)
                }
            }
            if (!is.null(regions) && nrow(regions) > 0) {
                region_stats <- all_outliers %>%
                    filter(!is.na(.data$chr)) %>%
                    inner_join(regions %>% select("chr", "region_start", "region_end"), by = "chr", relationship = "many-to-many") %>%
                    filter(.data$start <= .data$region_end, .data$end >= .data$region_start) %>%
                    group_by(.data$chr, .data$region_start, .data$region_end) %>%
                    summarise(outlier_stat = paste(sort(unique(.data$statistic)), collapse = ":"), .groups = "drop")
                regions <- regions %>% left_join(region_stats, by = c("chr", "region_start", "region_end"))
                regions <- add_region_stat_summaries(regions, windows_by_value_col)
                regions <- filter_regions_by_opts(regions, opts)
                if (nrow(regions) > 0) {
                    regions_file <- file.path(opts$`output-dir`, paste0("outlier_regions", window_suffix, chr_suffix, quantile_suffix, merge_suffix, "_merge_only", region_filename_suffix, ".csv"))
                    write_csv(regions, regions_file)
                    cat("Saved outlier regions to:", regions_file, "\n")
                    cat("Total regions:", nrow(regions), "\n")
                }
            }
        }
    }
} else {
    cat("Warning: No outlier windows were identified\n")
}

if (opts$verbose) cat("\nDone!\n")
