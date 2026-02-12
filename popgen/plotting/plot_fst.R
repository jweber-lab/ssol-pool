#!/usr/bin/env Rscript

###############################################################################
# plot_fst.R
# 
# Create ggplot2 plots of FST statistics from grenedalf FST output files.
#
# Features:
# - Chromosomal positions on x-axis with alternating white/grey stripes
# - Paneling by window size, sample pair, or both
# - Horizontal reference lines (median, mean, 95th percentile)
# - Statistics calculated genome-wide and per-chromosome
#
# Author: Based on workflow by JW
# Usage: See README.md or run with --help
###############################################################################

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(optparse)
    library(readr)
    library(tidyr)
})

initial_args <- commandArgs(trailingOnly = FALSE)
file_arg <- initial_args[substr(initial_args, 1, 7) == "--file="]
script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg)) else "."
plot_common_path <- file.path(script_dir, "plot_common.R")
if (!file.exists(plot_common_path)) plot_common_path <- file.path(getwd(), "plot_common.R")
if (file.exists(plot_common_path)) source(plot_common_path)

# Parse command-line arguments
option_list <- list(
    make_option(c("--input-dir"), type="character", default=NULL,
                help="Directory containing FST TSV or HDF5 files", metavar="DIR"),
    make_option(c("--input-format"), type="character", default="auto",
                help="Input format: tsv, hdf5, or auto [default: %default]", metavar="FORMAT"),
    make_option(c("--y-value"), type="character", default="value",
                help="Y-axis: value, rank, or quantile [default: %default]", metavar="VALUE"),
    make_option(c("--plot-style"), type="character", default="line",
                help="Plot style: line or line_points [default: %default]", metavar="STYLE"),
    make_option(c("--overlay-pairs"), action="store_true", default=FALSE,
                help="Plot all pairs in one panel per window (color = pair); no facet by pair", metavar="FLAG"),
    make_option(c("--output-dir"), type="character", default=".",
                help="Output directory for plots [default: %default]", metavar="DIR"),
    make_option(c("--reference-genome"), type="character", default=NULL,
                help="Reference genome FASTA file (optional, for chromosome lengths)", metavar="FILE"),
    make_option(c("--panel-by"), type="character", default="both",
                help="Paneling option: window, pair, both, or none [default: %default]",
                metavar="OPTION"),
    make_option(c("--plot-format"), type="character", default="png",
                help="Plot format: png, pdf, svg, or any combination (e.g., both, all) [default: %default]", metavar="FORMAT"),
    make_option(c("--width"), type="numeric", default=12,
                help="Plot width in inches [default: %default]", metavar="NUMBER"),
    make_option(c("--height"), type="numeric", default=8,
                help="Plot height in inches [default: %default]", metavar="NUMBER"),
    make_option(c("--dpi"), type="numeric", default=NULL,
                help="DPI for PNG (default: 300)", metavar="NUMBER"),
    make_option(c("--file-prefix"), type="character", default="",
                help="Prefix for output files [default: no prefix]", metavar="PREFIX"),
    make_option(c("--sample-pairs"), type="character", default="",
                help="Comma-separated list of sample pairs to plot (default: all pairs)", metavar="LIST"),
    make_option(c("--chromosome"), type="character", default=NULL,
                help="Single chromosome/scaffold to plot (default: all chromosomes)", metavar="CHR"),
    make_option(c("--window-size"), type="numeric", default=NULL,
                help="Single window size to plot (default: all window sizes)", metavar="NUMBER"),
    make_option(c("--top-n-chromosomes"), type="numeric", default=NULL,
                help="Plot only the N longest chromosomes/scaffolds (optional)", metavar="N"),
    make_option(c("--min-chromosome-length"), type="numeric", default=NULL,
                help="Plot only chromosomes/scaffolds longer than N bp (optional)", metavar="N"),
    make_option(c("--verbose"), action="store_true", default=FALSE,
                help="Enable verbose output for debugging", metavar="FLAG"),
    make_option(c("--transform"), type="character", default="none",
                help="Y-axis transformation: none, asinh (centered on median), or log [default: %default]", metavar="TRANSFORM"),
    make_option(c("--asinh-scale"), type="numeric", default=NULL,
                help="Scale factor for asinh transformation (default: use global standard deviation)", metavar="NUMBER")
)

opt_parser <- OptionParser(option_list=option_list, usage="usage: %prog [options]")
opts <- parse_args(opt_parser)

# Validate arguments
if (is.null(opts$`input-dir`)) {
    stop("--input-dir is required")
}

if (!opts$`panel-by` %in% c("window", "pair", "both", "none")) {
    stop("--panel-by must be one of: window, pair, both, none")
}

if (!opts$`input-format` %in% c("tsv", "hdf5", "auto")) {
    stop("--input-format must be one of: tsv, hdf5, auto")
}
if (!opts$`y-value` %in% c("value", "rank", "quantile")) {
    stop("--y-value must be one of: value, rank, quantile")
}

if (!opts$`plot-style` %in% c("line", "line_points")) {
    stop("--plot-style must be one of: line, line_points")
}

if (!opts$transform %in% c("asinh", "log", "none")) {
    stop("--transform must be one of: asinh, log, none")
}

valid_formats <- c("png", "pdf", "svg", "both", "all")
format_parts <- strsplit(opts$`plot-format`, ",")[[1]]
format_parts <- trimws(format_parts)
if (!all(format_parts %in% valid_formats)) {
    stop("--plot-format must be one or more of: png, pdf, svg, both, all (comma-separated)")
}
# Expand "both" to "png,pdf" and "all" to "png,pdf,svg"
if ("both" %in% format_parts) {
    format_parts <- setdiff(format_parts, "both")
    format_parts <- c(format_parts, "png", "pdf")
}
if ("all" %in% format_parts) {
    format_parts <- setdiff(format_parts, "all")
    format_parts <- c(format_parts, "png", "pdf", "svg")
}
format_parts <- unique(format_parts)

# Create output directory if it doesn't exist
if (!dir.exists(opts$`output-dir`)) {
    dir.create(opts$`output-dir`, recursive=TRUE)
}

input_dir <- opts$`input-dir`
input_format_used <- opts$`input-format`
if (input_format_used == "auto") {
    h5_candidates <- list.files(input_dir, pattern = "^fst_.*\\.h5$", full.names = TRUE, recursive = FALSE)
    input_format_used <- if (length(h5_candidates) > 0) "hdf5" else "tsv"
}
if (input_format_used == "hdf5" && !requireNamespace("hdf5r", quietly = TRUE)) {
    stop("--input-format hdf5 requires package 'hdf5r'. Install with: install.packages(\"hdf5r\")")
}
summary_tsv <- NULL

# Function to read FASTA index (.fai) to get chromosome lengths
read_fai <- function(fai_file) {
    if (!file.exists(fai_file)) {
        return(NULL)
    }
    fai <- read.table(fai_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    colnames(fai) <- c("chr", "length", "offset", "linebases", "linewidth")
    return(fai[, c("chr", "length")])
}

# Function to get chromosome lengths using tidyverse
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
    cat("Inferring chromosome lengths from data...\n")
    chr_lengths <- data %>%
        group_by(.data$chr) %>%
        summarise(length = max(.data$pos, na.rm=TRUE), .groups="drop") %>%
        filter(!is.na(.data$length), .data$length > 0) %>%
        arrange(desc(.data$length))  # Order by length, longest first
    
    cat("Found", nrow(chr_lengths), "chromosomes (ordered by length, longest first)\n")
    return(chr_lengths)
}

# Function to calculate cumulative positions using tidyverse
add_cumulative_positions <- function(data, chr_lengths) {
    # Use the order from chr_lengths (already sorted by length, longest first)
    # Match chromosomes that exist in both data and chr_lengths
    chr_order <- chr_lengths$chr[chr_lengths$chr %in% unique(data$chr)]
    
    # Calculate cumulative lengths using tidyverse
    chr_lengths_ordered <- chr_lengths %>%
        filter(.data$chr %in% chr_order) %>%
        slice(match(chr_order, .data$chr)) %>%
        mutate(cum_start = c(0, cumsum(.data$length)[-n()]))
    
    # Merge with data
    data <- data %>%
        left_join(chr_lengths_ordered, by="chr") %>%
        mutate(cum_pos = .data$cum_start + .data$pos)
    
    return(data)
}

# Function to create chromosome stripe rectangles using tidyverse
create_chr_stripes <- function(chr_lengths) {
    # Use the order from chr_lengths (already sorted by length, longest first)
    chr_order <- chr_lengths$chr
    
    stripes <- chr_lengths %>%
        slice(match(chr_order, .data$chr)) %>%
        mutate(
            cum_start = c(0, cumsum(.data$length)[-n()]),
            cum_end = cumsum(.data$length),
            fill = rep(c("white", "grey90"), length.out=n())
        ) %>%
        select(xmin = cum_start, xmax = cum_end, fill)
    
    return(stripes)
}

# Function to normalize sample names (remove common suffixes)
normalize_sample_name <- function(sample_name) {
    # Convert to lowercase for consistency
    normalized <- tolower(sample_name)
    
    # Remove common suffixes (order matters - remove longer patterns first)
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

# Find all FST files and read (TSV or HDF5)
selected_window_size <- opts$`window-size`
selected_pairs <- NULL
if (opts$`sample-pairs` != "") {
    selected_pairs <- strsplit(opts$`sample-pairs`, ",")[[1]]
    selected_pairs <- trimws(selected_pairs)
    selected_pairs_normalized <- sapply(selected_pairs, function(p) {
        parts <- strsplit(p, ":")[[1]]
        if (length(parts) == 2) {
            sample1 <- normalize_sample_name(parts[1])
            sample2 <- normalize_sample_name(parts[2])
            return(paste(sample1, sample2, sep=":"))
        }
        return(p)
    })
}

combined_data_long_format <- FALSE
if (input_format_used == "hdf5") {
    all_h5 <- list.files(input_dir, pattern = "^fst_.*\\.h5$", full.names = TRUE, recursive = FALSE)
    if (length(all_h5) == 0) stop("No FST HDF5 files found in: ", input_dir)
    if (!is.null(selected_window_size)) {
        file_win_sizes <- vapply(all_h5, function(f) { ps <- parse_window_step_from_filename(basename(f)); ps$window_size }, numeric(1))
        all_h5 <- all_h5[!is.na(file_win_sizes) & file_win_sizes == selected_window_size]
        if (length(all_h5) == 0) stop("No HDF5 files matching window size ", selected_window_size)
    }
    all_data_list <- list()
    summary_list <- list()
    for (path in all_h5) {
        grp_name <- if (grepl("single", basename(path), ignore.case = TRUE)) "sites" else "windows"
        d <- read_h5_windows(path, grp_name)
        if (is.null(d) || nrow(d) == 0) next
        ps <- parse_window_step_from_filename(basename(path))
        d$window_size <- ps$window_size
        d$step_size <- ps$step_size
        if (!"pos" %in% names(d) && "start" %in% names(d)) d$pos <- as.numeric(d$start)
        else if (!"pos" %in% names(d) && "end" %in% names(d)) d$pos <- as.numeric(d$end)
        else if (!"pos" %in% names(d) && "start" %in% names(d) && "end" %in% names(d)) d$pos <- (as.numeric(d$start) + as.numeric(d$end)) / 2
        if ("pop1" %in% names(d) && "pop2" %in% names(d)) d$sample_pair <- paste(d$pop1, d$pop2, sep = ":")
        all_data_list[[path]] <- d
        sum_t <- read_summary_tsv_for_h5(path)
        if (!is.null(sum_t)) summary_list[[path]] <- sum_t
    }
    combined_data <- bind_rows(all_data_list)
    combined_data_long_format <- TRUE
    if (length(summary_list) > 0) summary_tsv <- bind_rows(summary_list)
    if (!is.null(selected_pairs_normalized)) {
        combined_data <- combined_data %>% filter(.data$sample_pair %in% selected_pairs_normalized)
        if (nrow(combined_data) == 0) stop("None of the specified sample pairs found in HDF5 data")
    }
    if (opts$verbose) cat("Read", length(all_h5), "HDF5 file(s); combined", nrow(combined_data), "rows (long format)\n")
} else {
    tsv_files <- list.files(input_dir, pattern = ".*fst.*\\.tsv$", full.names = TRUE, recursive = TRUE)
    csv_files <- list.files(input_dir, pattern = ".*fst.*\\.csv$", full.names = TRUE, recursive = TRUE)
    all_files <- c(tsv_files, csv_files)
    if (length(all_files) == 0) stop("No FST files found in: ", input_dir)
    if (!is.null(selected_window_size)) {
        file_window_sizes <- sapply(all_files, function(f) extract_window_size_from_filename(basename(f)))
        has_single_label <- grepl("single", basename(all_files), ignore.case = TRUE)
        all_files <- all_files[!has_single_label & (is.na(file_window_sizes) | file_window_sizes == selected_window_size)]
        if (length(all_files) == 0) stop("No files matching window size ", selected_window_size)
    }
    if (!is.null(selected_pairs) && length(selected_pairs_normalized) > 0) {
        if (opts$verbose) cat("Filtering files by sample pairs:", paste(selected_pairs_normalized, collapse = ", "), "\n")
        files_to_keep <- character(0)
        for (file in all_files) {
        col_names <- peek_column_names(file)
        fst_cols <- col_names[grepl("\\.fst$", col_names, ignore.case=TRUE)]
        
        # Check if any FST column matches a selected pair
        file_has_pair <- FALSE
        for (col in fst_cols) {
            # Parse pair from column name
            if (grepl(":", col)) {
                pair_match <- sub("\\.fst$", "", col, ignore.case=TRUE)
                pair_parts <- strsplit(pair_match, ":")[[1]]
                if (length(pair_parts) == 2) {
                    sample1 <- normalize_sample_name(pair_parts[1])
                    sample2 <- normalize_sample_name(pair_parts[2])
                    pair_name <- paste(sample1, sample2, sep=":")
                    if (pair_name %in% selected_pairs_normalized) {
                        file_has_pair <- TRUE
                        break
                    }
                }
            } else {
                # Try underscore-separated format
                pair_match <- sub("\\.fst$", "", col, ignore.case=TRUE)
                pair_match <- sub("_fst$", "", pair_match, ignore.case=TRUE)
                parts <- strsplit(pair_match, "_")[[1]]
                if (length(parts) >= 2) {
                    sample1_raw <- paste(parts[1:(length(parts)-1)], collapse="_")
                    sample2_raw <- parts[length(parts)]
                    sample1 <- normalize_sample_name(sample1_raw)
                    sample2 <- normalize_sample_name(sample2_raw)
                    pair_name <- paste(sample1, sample2, sep=":")
                    if (pair_name %in% selected_pairs_normalized) {
                        file_has_pair <- TRUE
                        break
                    }
                }
            }
        }
        
        if (file_has_pair) {
            files_to_keep <- c(files_to_keep, file)
        }
    }
    all_files <- files_to_keep
    if (length(all_files) == 0) stop("No files found containing the specified sample pairs")
    if (opts$verbose) cat("Filtered to", length(all_files), "file(s) containing requested sample pairs\n")
    }
    all_data <- list()
    for (file in all_files) {
    if (opts$verbose) {
        cat("\n=== Reading:", basename(file), "===\n")
    }
    
    # Detect separator by reading first line
    first_line <- readLines(file, n=1)
    has_tabs <- grepl("\t", first_line)
    has_commas <- grepl(",", first_line)
    
    # Determine separator and file type
    if (has_tabs && (!has_commas || length(strsplit(first_line, "\t")[[1]]) > length(strsplit(first_line, ",")[[1]]))) {
        # Tab-separated (more tabs than commas, or only tabs)
        data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        detected_sep <- "tab"
        if (opts$verbose) cat("  Detected separator: tab\n")
    } else if (has_commas) {
        # Comma-separated
        data <- read_csv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        detected_sep <- "comma"
        if (opts$verbose) cat("  Detected separator: comma\n")
    } else {
        # Default to TSV if neither detected (shouldn't happen, but fallback)
        data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        detected_sep <- "tab"
        if (opts$verbose) cat("  Detected separator: tab (fallback)\n")
    }
    
    # Warn if extension doesn't match separator
    file_ext <- tolower(tools::file_ext(file))
    if ((detected_sep == "tab" && file_ext == "csv") || (detected_sep == "comma" && file_ext == "tsv")) {
        if (opts$verbose) cat("  Warning: File has", file_ext, "extension but", detected_sep, "separator\n")
    }
    
    if (opts$verbose) {
        cat("  Original columns:", length(colnames(data)), "\n")
        cat("  Column names:", paste(head(colnames(data), 10), collapse=", "), if(length(colnames(data)) > 10) "..." else "", "\n")
    }
    
    # Standardize column names to lowercase
    data <- data %>%
        rename_all(tolower)
    
    if (opts$verbose) {
        cat("  After tolower - columns:", paste(head(colnames(data), 10), collapse=", "), if(length(colnames(data)) > 10) "..." else "", "\n")
    }
    
    # Convert base numeric columns (start, end, position, pos) to numeric using mutate
    # NOTE: Do NOT convert chrom/chromosome/chr to numeric - they are character identifiers!
    data <- data %>%
        mutate(across(any_of(c("start", "end", "position", "pos")), as.numeric))
    
    # Convert FST columns to numeric (sample-prefixed columns ending in .fst)
    fst_stat_cols <- grep("\\.fst$", colnames(data), ignore.case=FALSE, value=TRUE)
    if (length(fst_stat_cols) > 0) {
        if (opts$verbose) {
            cat("  Found FST columns:", paste(head(fst_stat_cols, 5), collapse=", "), if(length(fst_stat_cols) > 5) "..." else "", "\n")
        }
        data <- data %>%
            mutate(across(all_of(fst_stat_cols), as.numeric))
    } else {
        if (opts$verbose) cat("  Warning: No FST columns found matching pattern *.fst\n")
    }
    
    # Handle different column name variations - ensure chr and pos exist using tidyverse
    # Check for chromosome column (various names) - MUST be character, not numeric!
    if (!"chr" %in% colnames(data)) {
        if ("chromosome" %in% colnames(data)) {
            data <- data %>% mutate(chr = as.character(chromosome))
            if (opts$verbose) cat("  Mapped 'chromosome' -> 'chr' (as character)\n")
        } else if ("chrom" %in% colnames(data)) {
            data <- data %>% mutate(chr = as.character(chrom))
            if (opts$verbose) cat("  Mapped 'chrom' -> 'chr' (as character)\n")
        } else if ("contig" %in% colnames(data)) {
            data <- data %>% mutate(chr = as.character(contig))
            if (opts$verbose) cat("  Mapped 'contig' -> 'chr' (as character)\n")
        } else {
            # Try to find any column that might be chromosome
            chr_col <- grep("^(chr|chrom|contig|scaffold)", colnames(data), ignore.case=FALSE, value=TRUE)
            if (length(chr_col) > 0) {
                data <- data %>% mutate(chr = as.character(.data[[chr_col[1]]]))
                if (opts$verbose) cat("  Mapped '", chr_col[1], "' -> 'chr' (as character)\n", sep="")
            } else {
                cat("  ERROR: Could not find chromosome column\n")
                cat("  Available columns:", paste(colnames(data), collapse=", "), "\n")
            }
        }
    } else {
        # Ensure chr is character type (in case it was read as something else)
        data <- data %>% mutate(chr = as.character(chr))
        if (opts$verbose) cat("  Found 'chr' column (ensured character type)\n")
    }
    
    # Check for position column (various names) - convert to numeric
    # Prefer 'start' over 'end' for window-based data
    if (!"pos" %in% colnames(data)) {
        if ("start" %in% colnames(data)) {
            data <- data %>% mutate(pos = as.numeric(start))
            if (opts$verbose) cat("  Mapped 'start' -> 'pos' (as numeric)\n")
        } else if ("position" %in% colnames(data)) {
            data <- data %>% mutate(pos = as.numeric(position))
            if (opts$verbose) cat("  Mapped 'position' -> 'pos' (as numeric)\n")
        } else if ("end" %in% colnames(data)) {
            data <- data %>% mutate(pos = as.numeric(end))
            if (opts$verbose) cat("  Mapped 'end' -> 'pos' (as numeric)\n")
        } else {
            # Try to find any column that might be position
            pos_col <- grep("^(pos|start|end|coordinate)", colnames(data), ignore.case=FALSE, value=TRUE)
            if (length(pos_col) > 0) {
                data <- data %>% mutate(pos = as.numeric(.data[[pos_col[1]]]))
                if (opts$verbose) cat("  Mapped '", pos_col[1], "' -> 'pos' (as numeric)\n", sep="")
            } else {
                cat("  ERROR: Could not find position column\n")
                cat("  Available columns:", paste(colnames(data), collapse=", "), "\n")
            }
        }
    } else {
        # Ensure pos is numeric type
        data <- data %>% mutate(pos = as.numeric(pos))
        if (opts$verbose) cat("  Found 'pos' column (ensured numeric type)\n")
    }
    
    # Extract window size from filename if present
    filename <- basename(file)
    window_match <- regmatches(filename, regexpr("w(\\d+)_s(\\d+)", filename))
    if (length(window_match) > 0) {
        window_parts <- strsplit(window_match, "_")[[1]]
        window_size <- as.numeric(sub("w", "", window_parts[1]))
        step_size <- as.numeric(sub("s", "", window_parts[2]))
        data <- data %>% mutate(window_size = window_size, step_size = step_size)
    } else {
        data <- data %>% mutate(window_size = NA_real_, step_size = NA_real_)
    }
    
    if (opts$verbose) cat("  Rows:", nrow(data), "\n")
    all_data[[file]] <- data
}

    if (opts$verbose) cat("\n=== Combining data from", length(all_data), "file(s) ===\n")
    combined_data <- bind_rows(all_data)
}

if (opts$verbose) {
    cat("Combined data dimensions:", nrow(combined_data), "rows,", ncol(combined_data), "columns\n")
    cat("Combined columns:", paste(head(colnames(combined_data), 15), collapse=", "), if(length(colnames(combined_data)) > 15) "..." else "", "\n")
}

# Check if required columns exist
if (!"chr" %in% colnames(combined_data)) {
    stop("Error: 'chr' column not found in any input files. Available columns: ", paste(colnames(combined_data), collapse=", "))
}
if (!"pos" %in% colnames(combined_data)) {
    stop("Error: 'pos' column not found in any input files. Available columns: ", paste(colnames(combined_data), collapse=", "))
}

if (opts$verbose) {
    cat("Found required columns: chr, pos\n")
    
    # Debug: Check chr and pos columns before filtering
    cat("Before filtering:\n")
    cat("  chr column type:", class(combined_data$chr), "\n")
    cat("  pos column type:", class(combined_data$pos), "\n")
    cat("  chr NA count:", sum(is.na(combined_data$chr)), "out of", nrow(combined_data), "\n")
    cat("  pos NA count:", sum(is.na(combined_data$pos)), "out of", nrow(combined_data), "\n")
    cat("  chr unique values (first 10):", paste(head(unique(combined_data$chr), 10), collapse=", "), "\n")
    cat("  pos range:", min(combined_data$pos, na.rm=TRUE), "to", max(combined_data$pos, na.rm=TRUE), "\n")

# Remove rows with missing chromosome or position
    cat("\nFiltering rows with missing chr or pos...\n")
    cat("Before filtering:", nrow(combined_data), "rows\n")
}
combined_data <- combined_data %>%
    filter(!is.na(.data$chr), !is.na(.data$pos))
if (opts$verbose) {
    cat("After filtering:", nrow(combined_data), "rows\n")
}

if (nrow(combined_data) == 0) {
    stop("Error: All rows were filtered out! Check chr and pos columns.")
}

# Parse FST columns (or use long-format pairs when from HDF5)
if (combined_data_long_format) {
    y_col <- switch(opts$`y-value`, rank = "fst_rank", quantile = "fst_quantile", "fst")
    if (!y_col %in% colnames(combined_data)) stop("--y-value ", opts$`y-value`, " requested but column '", y_col, "' not found in HDF5.")
    pair_names <- unique(combined_data$sample_pair[!is.na(combined_data$sample_pair)])
    fst_pairs <- setNames(rep(y_col, length(pair_names)), pair_names)
    if (!is.null(selected_pairs_normalized)) fst_pairs <- fst_pairs[names(fst_pairs) %in% selected_pairs_normalized]
    if (length(fst_pairs) == 0) stop("None of the specified sample pairs found in data")
} else {
fst_cols <- colnames(combined_data)[grepl("\\.fst$", colnames(combined_data), ignore.case=TRUE)]
if (length(fst_cols) == 0) {
    stop("No FST columns found in input files. Expected columns ending in .fst")
}

if (opts$verbose) {
    cat("Found", length(fst_cols), "FST column(s)\n")
    cat("FST columns:", paste(fst_cols, collapse=", "), "\n")
}

# Parse sample pairs from column names
# First pass: collect all mappings (normalized pair -> list of columns)
fst_pairs_raw <- list()
for (col in fst_cols) {
    # Try colon-separated format first
    if (grepl(":", col)) {
        pair_match <- sub("\\.fst$", "", col, ignore.case=TRUE)
        pair_parts <- strsplit(pair_match, ":")[[1]]
        if (length(pair_parts) == 2) {
            sample1_raw <- pair_parts[1]
            sample2_raw <- pair_parts[2]
            # Normalize sample names
            sample1 <- normalize_sample_name(sample1_raw)
            sample2 <- normalize_sample_name(sample2_raw)
            pair_name <- paste(sample1, sample2, sep=":")
            if (!pair_name %in% names(fst_pairs_raw)) {
                fst_pairs_raw[[pair_name]] <- list()
            }
            fst_pairs_raw[[pair_name]] <- c(fst_pairs_raw[[pair_name]], col)
        }
    } else {
        # Try underscore-separated format
        pair_match <- sub("\\.fst$", "", col, ignore.case=TRUE)
        pair_match <- sub("_fst$", "", pair_match, ignore.case=TRUE)
        parts <- strsplit(pair_match, "_")[[1]]
        if (length(parts) >= 2) {
            sample1_raw <- paste(parts[1:(length(parts)-1)], collapse="_")
            sample2_raw <- parts[length(parts)]
            # Normalize sample names
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

# Second pass: for each normalized pair, select the best column (one with most non-NA data)
fst_pairs <- list()
for (pair_name in names(fst_pairs_raw)) {
    cols_for_pair <- fst_pairs_raw[[pair_name]]
    if (length(cols_for_pair) == 1) {
        # Only one column, use it
        fst_pairs[[pair_name]] <- cols_for_pair[1]
    } else {
        # Multiple columns map to same normalized pair
        # Choose the one with the most non-NA data
        if (opts$verbose) {
            cat("  Multiple columns for pair '", pair_name, "': ", paste(cols_for_pair, collapse=", "), "\n", sep="")
        }
        best_col <- NULL
        best_count <- -1
        for (col in cols_for_pair) {
            if (col %in% colnames(combined_data)) {
                non_na_count <- sum(!is.na(combined_data[[col]]))
                if (opts$verbose) {
                    cat("    ", col, ": ", non_na_count, " non-NA values\n", sep="")
                }
                if (non_na_count > best_count) {
                    best_count <- non_na_count
                    best_col <- col
                }
            }
        }
        if (!is.null(best_col)) {
            fst_pairs[[pair_name]] <- best_col
            if (opts$verbose && length(cols_for_pair) > 1) {
                cat("    Selected: ", best_col, " (", best_count, " non-NA values)\n", sep="")
            }
        } else {
            # Fallback: use first column
            fst_pairs[[pair_name]] <- cols_for_pair[1]
            if (opts$verbose) {
                cat("    Warning: No valid column found, using first: ", cols_for_pair[1], "\n", sep="")
            }
        }
    }
}

if (length(fst_pairs) == 0) {
    stop("Could not parse sample pairs from FST column names")
}

if (opts$verbose) {
    cat("\nNormalized sample pairs:\n")
    for (pair_name in names(fst_pairs)) {
        cat("  ", pair_name, " -> ", fst_pairs[[pair_name]], "\n", sep="")
    }
}
}

# Filter pairs if specified (use normalized pair names)
if (!is.null(selected_pairs)) {
    # selected_pairs_normalized was created during file filtering
    # If it doesn't exist (no file filtering happened), create it now
    if (!exists("selected_pairs_normalized")) {
        selected_pairs_normalized <- sapply(selected_pairs, function(p) {
            parts <- strsplit(p, ":")[[1]]
            if (length(parts) == 2) {
                sample1 <- normalize_sample_name(parts[1])
                sample2 <- normalize_sample_name(parts[2])
                return(paste(sample1, sample2, sep=":"))
            }
            return(p)
        })
    }
    fst_pairs <- fst_pairs[names(fst_pairs) %in% selected_pairs_normalized]
    if (length(fst_pairs) == 0) {
        stop("None of the specified sample pairs were found in the data")
    }
}

if (opts$verbose) {
    cat("Processing", length(fst_pairs), "sample pair(s)\n")
}

# Filter to specific chromosome if requested
selected_chr <- opts$chromosome
if (!is.null(selected_chr)) {
    if (opts$verbose) cat("\nFiltering to chromosome:", selected_chr, "\n")
    if (!selected_chr %in% combined_data$chr) {
        stop("Error: Chromosome '", selected_chr, "' not found in data. Available chromosomes: ", 
             paste(unique(combined_data$chr), collapse=", "))
    }
    combined_data <- combined_data %>%
        filter(.data$chr == selected_chr)
    if (opts$verbose) cat("Filtered data:", nrow(combined_data), "rows\n")
}

# Filter to specific window size if requested (already filtered files, but double-check data)
if (!is.null(selected_window_size)) {
    if (opts$verbose) {
        cat("\nFiltering to window size:", selected_window_size, "\n")
    }
    if ("window_size" %in% colnames(combined_data)) {
        available_windows <- unique(combined_data$window_size[!is.na(combined_data$window_size)])
        if (!selected_window_size %in% available_windows) {
            stop("Error: Window size '", selected_window_size, "' not found in data. Available window sizes: ", 
                 paste(available_windows, collapse=", "))
        }
        combined_data <- combined_data %>%
            filter(is.na(.data$window_size) | .data$window_size == selected_window_size)
        if (opts$verbose) {
            cat("Filtered data:", nrow(combined_data), "rows\n")
        }
    } else {
        if (opts$verbose) {
            cat("Warning: No window_size column found, ignoring --window-size filter\n")
        }
    }
}

# Get chromosome lengths
chr_lengths <- get_chr_lengths(combined_data, opts$`reference-genome`)

# Filter chromosomes by length if requested
if (!is.null(opts$`top-n-chromosomes`)) {
    top_n <- as.integer(opts$`top-n-chromosomes`)
    if (top_n > 0 && top_n <= nrow(chr_lengths)) {
        # chr_lengths is already sorted by length (longest first)
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

# Filter combined_data to only include selected chromosomes
if (nrow(chr_lengths) > 0) {
    selected_chrs <- chr_lengths$chr
    combined_data <- combined_data %>%
        filter(.data$chr %in% selected_chrs)
    if (opts$verbose) cat("Filtered data to", length(selected_chrs), "chromosome(s):", nrow(combined_data), "rows\n")
} else {
    if (opts$verbose) cat("Warning: No chromosomes remain after filtering\n")
}

# Add cumulative positions
combined_data <- add_cumulative_positions(combined_data, chr_lengths)

# Create chromosome stripes
chr_stripes <- create_chr_stripes(chr_lengths)

# Process each sample pair (or overlay all pairs in one panel)
stats_results <- list()

overlay_pairs <- opts$`overlay-pairs` && combined_data_long_format && "sample_pair" %in% colnames(combined_data)
if (overlay_pairs && length(fst_pairs) == 0) overlay_pairs <- FALSE
if (overlay_pairs) {
    y_col_overlay <- fst_pairs[[1L]]
    overlay_data <- combined_data %>%
        rename(fst = .data[[y_col_overlay]]) %>%
        filter(!is.na(.data$fst))
    if (nrow(overlay_data) == 0) stop("No FST data for overlay")
    transform_this <- if (opts$`y-value` %in% c("rank", "quantile")) "none" else opts$transform
    if (transform_this == "asinh") {
        global_median <- median(overlay_data$fst, na.rm = TRUE)
        scale_factor <- if (is.null(opts$`asinh-scale`)) sd(overlay_data$fst, na.rm = TRUE) else opts$`asinh-scale`
        overlay_data <- overlay_data %>% mutate(fst = asinh((.data$fst - global_median) / scale_factor))
    } else if (transform_this == "log") {
        min_val <- min(overlay_data$fst, na.rm = TRUE)
        if (min_val <= 0) overlay_data <- overlay_data %>% mutate(fst = log1p(pmax(0, .data$fst + 1e-10)))
        else overlay_data <- overlay_data %>% mutate(fst = log(.data$fst))
    }
    pair_names_sorted <- sort(unique(overlay_data$sample_pair))
    p <- ggplot(overlay_data, aes(x = cum_pos, y = fst, color = .data$sample_pair, group = .data$sample_pair)) +
        scale_color_manual(name = "Pair", values = setNames(
            rep(PLOT_PALETTE_QUALITATIVE, length.out = length(pair_names_sorted)),
            pair_names_sorted), drop = FALSE) +
        geom_rect(data = chr_stripes, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.3) +
        scale_fill_identity() +
        geom_line(alpha = 0.7, linewidth = 0.5, na.rm = TRUE) +
        labs(x = "Genomic Position (cumulative)",
            y = if (transform_this == "asinh") "asinh((FST - median) / SD)" else if (transform_this == "log") "log(FST)" else if (opts$`y-value` == "rank") "FST rank" else if (opts$`y-value` == "quantile") "FST quantile" else "FST (Fixation Index)",
            title = "FST (all pairs)") +
        theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    if (opts$`plot-style` == "line_points") p <- p + geom_point(size = 1, alpha = 0.7, na.rm = TRUE)
    has_window_size <- "window_size" %in% colnames(overlay_data) && !all(is.na(overlay_data$window_size))
    unique_windows <- if (has_window_size) unique(overlay_data$window_size[!is.na(overlay_data$window_size)]) else NULL
    multiple_windows <- has_window_size && length(unique_windows) > 1
    if (multiple_windows) p <- p + facet_wrap(~ window_size, scales = "fixed", ncol = 1)
    chr_labels <- chr_lengths$chr[chr_lengths$chr %in% unique(overlay_data$chr)]
    chr_centers <- chr_stripes %>% mutate(chr_idx = row_number()) %>%
        filter(chr_idx <= length(chr_labels)) %>% mutate(center = (.data$xmin + .data$xmax) / 2) %>% select(center)
    if (nrow(chr_centers) == length(chr_labels) && nrow(chr_centers) > 0) {
        p <- p + scale_x_continuous(breaks = chr_centers$center, labels = chr_labels, expand = c(0.02, 0))
    }
    out_name <- paste0(opts$`file-prefix`, "fst_overlay")
    if (opts$`y-value` != "value") out_name <- paste0(out_name, "_", opts$`y-value`)
    for (fmt in strsplit(opts$`plot-format`, ",")[[1]]) {
        fmt <- trimws(tolower(fmt))
        if (fmt %in% c("png", "pdf", "svg")) {
            ext <- if (fmt == "png") ".png" else if (fmt == "pdf") ".pdf" else ".svg"
            ggsave(file.path(opts$`output-dir`, paste0(out_name, ext)), plot = p, width = opts$width, height = opts$height, dpi = PLOT_DPI)
        }
    }
} else {

for (pair_name in names(fst_pairs)) {
    fst_col <- fst_pairs[[pair_name]]
    if (opts$verbose) {
        cat("\nProcessing sample pair:", pair_name, "\n")
    }
    
    # Extract FST values for this pair using tidyverse
    # Check if fst_col exists in combined_data
    if (!fst_col %in% colnames(combined_data)) {
        cat("Warning: Column '", fst_col, "' not found in combined_data for pair '", pair_name, "', skipping\n", sep="")
        if (opts$verbose) {
            cat("  Available columns:", paste(head(colnames(combined_data), 20), collapse=", "), if(length(colnames(combined_data)) > 20) "..." else "", "\n")
            cat("  Looking for column:", fst_col, "\n")
            cat("  Total rows in combined_data:", nrow(combined_data), "\n")
        }
        next
    }
    
    # Check how many non-NA values exist in this column after filtering
    if (opts$verbose) {
        non_na_count <- sum(!is.na(combined_data[[fst_col]]))
        cat("  Column '", fst_col, "' has ", non_na_count, " non-NA values in filtered combined_data\n", sep="")
        if (non_na_count == 0) {
            cat("  ERROR: Column exists but has no non-NA values after filtering!\n")
            cat("  This suggests the column was selected based on unfiltered data but doesn't exist in filtered data\n")
        }
    }
    
    # Select columns that exist (cum_pos might not exist if add_cumulative_positions failed)
    cols_to_select <- c("chr", "pos", "window_size", fst_col)
    if ("cum_pos" %in% colnames(combined_data)) {
        cols_to_select <- c("chr", "pos", "cum_pos", "window_size", fst_col)
    }
    if (combined_data_long_format && "sample_pair" %in% colnames(combined_data)) {
        pair_data <- combined_data %>%
            filter(.data$sample_pair == pair_name) %>%
            select(any_of(c(cols_to_select, "sample_pair"))) %>%
            rename(fst = all_of(fst_col)) %>%
            filter(!is.na(.data$fst))
    } else {
    pair_data <- combined_data %>%
        select(all_of(cols_to_select)) %>%
        rename(fst = all_of(fst_col)) %>%
        filter(!is.na(.data$fst))
    }
    
    # Add cum_pos if it's missing (shouldn't happen, but safety check)
    if (!"cum_pos" %in% colnames(pair_data)) {
        if (opts$verbose) {
            cat("  Warning: cum_pos missing, recalculating...\n")
        }
        pair_data <- add_cumulative_positions(pair_data, chr_lengths)
    }
    
    # Sort by cum_pos to ensure line connects properly (critical when combining multiple files)
    pair_data <- pair_data %>%
        arrange(.data$cum_pos) %>%
        mutate(sample_pair = pair_name)
    
    if (opts$verbose) {
        cat("  Data sorted by cum_pos for proper line connection\n")
        # Check for any duplicate positions that might cause issues
        dup_positions <- sum(duplicated(pair_data$cum_pos))
        if (dup_positions > 0) {
            cat("  Warning: Found", dup_positions, "duplicate cum_pos values (may cause line artifacts)\n")
        }
    }
    
    if (nrow(pair_data) == 0) {
        cat("Warning: No valid FST data for pair '", pair_name, "', skipping\n", sep="")
        if (opts$verbose) {
            cat("  Total rows in combined_data:", nrow(combined_data), "\n")
            cat("  Non-NA values in column '", fst_col, "': ", sum(!is.na(combined_data[[fst_col]])), "\n", sep="")
        }
        next
    }
    
    if (opts$verbose) {
        cat("  Extracted", nrow(pair_data), "rows with valid FST data\n")
        cat("  FST range:", min(pair_data$fst, na.rm=TRUE), "to", max(pair_data$fst, na.rm=TRUE), "\n")
    }
    
    # Store original data for statistics calculation
    pair_data_original <- pair_data
    transform_this <- if (opts$`y-value` %in% c("rank", "quantile")) "none" else opts$transform

    # Calculate global median and standard deviation for asinh transformation (if needed)
    # Also store log shift if needed for axis label inverse transformation
    global_median <- NULL
    scale_factor <- NULL
    log_shift <- NULL
    if (transform_this == "asinh") {
        global_median <- median(pair_data_original$fst, na.rm=TRUE)
        global_sd <- sd(pair_data_original$fst, na.rm=TRUE)
        
        # Use provided scale or default to standard deviation
        scale_factor <- opts$`asinh-scale`
        if (is.null(scale_factor)) {
            scale_factor <- global_sd
            if (opts$verbose) cat("Applied asinh transformation: center=", global_median, ", scale=", scale_factor, " (global SD)\n", sep="")
        } else {
            if (opts$verbose) cat("Applied asinh transformation: center=", global_median, ", scale=", scale_factor, " (user-specified)\n", sep="")
        }
    }
    
    # Apply transformation if requested (for plotting)
    if (transform_this != "none") {
        if (transform_this == "asinh") {
            # Apply asinh transformation: asinh((x - global_median) / scale_factor)
            pair_data <- pair_data %>%
                mutate(fst = asinh((.data$fst - global_median) / scale_factor))
        } else if (transform_this == "log") {
            # Check for zeros or negatives
            min_val <- min(pair_data_original$fst, na.rm=TRUE)
            if (min_val <= 0) {
                # Use log1p to handle zeros/negatives: log1p(x) = log(1 + x)
                # Shift by minimum value to ensure all positive
                log_shift <- abs(min_val) + 1
                pair_data <- pair_data %>%
                    mutate(fst = log1p(.data$fst + log_shift))
                if (opts$verbose) cat("Applied log1p transformation with shift:", log_shift, "\n")
            } else {
                log_shift <- 0  # No shift needed
                pair_data <- pair_data %>%
                    mutate(fst = log(.data$fst))
                if (opts$verbose) cat("Applied log transformation\n")
            }
        }
    }
    
    # Calculate statistics from ORIGINAL data (before transformation)
    # This ensures reference lines are calculated correctly
    if (opts$`panel-by` == "none" || opts$`panel-by` == "window") {
        # Genome-wide and per-chromosome, optionally by window
        if (opts$`panel-by` == "window" && "window_size" %in% colnames(pair_data_original) && !all(is.na(pair_data_original$window_size))) {
            stats_summary <- pair_data_original %>%
                group_by(.data$sample_pair, .data$window_size, .data$chr) %>%
                summarise(
                    median = median(.data$fst, na.rm=TRUE),
                    mean = mean(.data$fst, na.rm=TRUE),
                    sd = sd(.data$fst, na.rm=TRUE),
                    q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    pair_data_original %>%
                        group_by(.data$sample_pair, .data$window_size) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$fst, na.rm=TRUE),
                            mean = mean(.data$fst, na.rm=TRUE),
                            sd = sd(.data$fst, na.rm=TRUE),
                            q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        } else {
            stats_summary <- pair_data_original %>%
                group_by(.data$sample_pair, .data$chr) %>%
                summarise(
                    median = median(.data$fst, na.rm=TRUE),
                    mean = mean(.data$fst, na.rm=TRUE),
                    sd = sd(.data$fst, na.rm=TRUE),
                    q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    pair_data_original %>%
                        group_by(.data$sample_pair) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$fst, na.rm=TRUE),
                            mean = mean(.data$fst, na.rm=TRUE),
                            sd = sd(.data$fst, na.rm=TRUE),
                            q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
            if (!"window_size" %in% colnames(stats_summary)) {
                stats_summary$window_size <- NA
            }
        }
    } else {
        # Similar logic for other paneling options
        grouping_vars <- c("sample_pair")
        if ("window_size" %in% colnames(pair_data_original) && !all(is.na(pair_data_original$window_size))) {
            # If paneling includes window dimension (both or pair with multiple windows), include it
            if (opts$`panel-by` %in% c("both", "pair")) {
            grouping_vars <- c("sample_pair", "window_size")
            }
        }
        
        # Build group_by expression dynamically
        if (length(grouping_vars) == 1) {
            stats_summary <- pair_data_original %>%
                group_by(.data[[grouping_vars[1]]], .data$chr) %>%
                summarise(
                    median = median(.data$fst, na.rm=TRUE),
                    mean = mean(.data$fst, na.rm=TRUE),
                    sd = sd(.data$fst, na.rm=TRUE),
                    q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    pair_data_original %>%
                        group_by(.data[[grouping_vars[1]]]) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$fst, na.rm=TRUE),
                            mean = mean(.data$fst, na.rm=TRUE),
                            sd = sd(.data$fst, na.rm=TRUE),
                            q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        } else {
            stats_summary <- pair_data_original %>%
                group_by(.data[[grouping_vars[1]]], .data[[grouping_vars[2]]], .data$chr) %>%
                summarise(
                    median = median(.data$fst, na.rm=TRUE),
                    mean = mean(.data$fst, na.rm=TRUE),
                    sd = sd(.data$fst, na.rm=TRUE),
                    q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    pair_data_original %>%
                        group_by(.data[[grouping_vars[1]]], .data[[grouping_vars[2]]]) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$fst, na.rm=TRUE),
                            mean = mean(.data$fst, na.rm=TRUE),
                            sd = sd(.data$fst, na.rm=TRUE),
                            q5 = quantile(.data$fst, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$fst, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$fst, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$fst, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$fst, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$fst, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        }
        
        if (!"window_size" %in% colnames(stats_summary)) {
            stats_summary$window_size <- NA
        }
    }
    
    stats_results[[pair_name]] <- stats_summary
    
    # Prepare reference lines data for legend, panel-aware (pair/window_size if present)
    # First, determine which panel columns will actually be used for faceting
    # This must match the faceting logic below
    has_window_size <- "window_size" %in% colnames(pair_data) && !all(is.na(pair_data$window_size))
    if (has_window_size) {
        unique_windows <- unique(pair_data$window_size[!is.na(pair_data$window_size)])
        multiple_windows <- length(unique_windows) > 1
    } else {
        multiple_windows <- FALSE
    }
    unique_pairs <- unique(pair_data$sample_pair)
    multiple_pairs <- length(unique_pairs) > 1
    
    # Determine which panel columns will be used for faceting based on panel-by option
    facet_panel_cols <- character(0)
    if (opts$`panel-by` == "window" && has_window_size) {
        facet_panel_cols <- c("window_size")
        if (multiple_pairs) {
            facet_panel_cols <- c(facet_panel_cols, "sample_pair")
        }
    } else if (opts$`panel-by` == "pair" && multiple_pairs) {
        facet_panel_cols <- c("sample_pair")
        if (has_window_size && multiple_windows) {
            facet_panel_cols <- c(facet_panel_cols, "window_size")
        }
    } else if (opts$`panel-by` == "both" && has_window_size) {
        facet_panel_cols <- c("sample_pair", "window_size")
    }
    # If panel-by == "none", no panel columns
    
    # Only include panel columns that will actually be used for faceting
    panel_cols <- intersect(facet_panel_cols, colnames(stats_summary))
    
    # Apply transformation to reference lines if needed
    if (transform_this == "asinh") {
        # For asinh, reference lines need to be transformed using the global median and scale
        # The median itself becomes asinh((median - global_median) / scale_factor)
        # Other stats are transformed as asinh((stat - global_median) / scale_factor)
        stats_summary <- stats_summary %>%
            mutate(
                median_transformed = asinh((median - global_median) / scale_factor),
                mean_transformed = asinh((mean - global_median) / scale_factor),
                q5_transformed = asinh((q5 - global_median) / scale_factor),
                q1_transformed = asinh((q1 - global_median) / scale_factor),
                q0.2_transformed = asinh((q0.2 - global_median) / scale_factor),
                q95_transformed = asinh((q95 - global_median) / scale_factor),
                q99_transformed = asinh((q99 - global_median) / scale_factor),
                q99.8_transformed = asinh((q99.8 - global_median) / scale_factor)
            )
    } else if (transform_this == "log") {
        # For log transformation, apply the same transformation as to the data
        # Use original data to determine shift
        min_val <- min(pair_data_original$fst, na.rm=TRUE)
        if (min_val <= 0) {
            shift <- abs(min_val) + 1
            stats_summary <- stats_summary %>%
                mutate(
                    median_transformed = log1p(median + shift),
                    mean_transformed = log1p(mean + shift),
                    q5_transformed = log1p(q5 + shift),
                    q1_transformed = log1p(q1 + shift),
                    q0.2_transformed = log1p(q0.2 + shift),
                    q95_transformed = log1p(q95 + shift),
                    q99_transformed = log1p(q99 + shift),
                    q99.8_transformed = log1p(q99.8 + shift)
                )
        } else {
            stats_summary <- stats_summary %>%
                mutate(
                    median_transformed = log(median),
                    mean_transformed = log(mean),
                    q5_transformed = log(q5),
                    q1_transformed = log(q1),
                    q0.2_transformed = log(q0.2),
                    q95_transformed = log(q95),
                    q99_transformed = log(q99),
                    q99.8_transformed = log(q99.8)
                )
        }
    } else {
        # No transformation
        stats_summary <- stats_summary %>%
            mutate(
                median_transformed = median,
                mean_transformed = mean,
                q5_transformed = q5,
                q1_transformed = q1,
                q0.2_transformed = q0.2,
                q95_transformed = q95,
                q99_transformed = q99,
                q99.8_transformed = q99.8
            )
    }
    
    if (opts$verbose) {
        cat("  Panel structure for reference lines:\n")
        cat("    has_window_size:", has_window_size, "\n")
        cat("    multiple_windows:", multiple_windows, "\n")
        cat("    multiple_pairs:", multiple_pairs, "\n")
        cat("    panel-by option:", opts$`panel-by`, "\n")
        cat("    facet_panel_cols:", paste(facet_panel_cols, collapse=", "), "\n")
        cat("    panel_cols (after intersect):", paste(panel_cols, collapse=", "), "\n")
    }
    
    # Get unique panel combinations from pair_data to ensure ref_lines matches
    if (length(panel_cols) > 0) {
        pair_panels <- pair_data %>%
            select(all_of(panel_cols)) %>%
            distinct()
        if (opts$verbose) {
            cat("    Unique panel combinations in pair_data:", nrow(pair_panels), "\n")
            if (nrow(pair_panels) > 0) {
                print(pair_panels)
            }
        }
    } else {
        pair_panels <- data.frame()
        if (opts$verbose) {
            cat("    No panel columns, using empty data frame\n")
        }
    }
    
    if (!is.null(selected_chr)) {
        # Use stats from the selected chromosome
        chr_stats <- stats_summary %>% filter(.data$chr == selected_chr)
        if (opts$verbose) {
            cat("  Using stats from selected chromosome:", selected_chr, "\n")
            cat("    chr_stats rows:", nrow(chr_stats), "\n")
        }
        if (nrow(chr_stats) > 0) {
            # Ensure we have stats for all panel combinations in pair_data
            if (nrow(pair_panels) > 0 && length(panel_cols) > 0) {
                if (opts$verbose) {
                    cat("    Merging chr_stats with pair_panels by:", paste(panel_cols, collapse=", "), "\n")
                    cat("    chr_stats before merge:", nrow(chr_stats), "rows\n")
                }
                # Merge with pair_panels to ensure all combinations are present
                chr_stats <- pair_panels %>%
                    left_join(chr_stats, by=panel_cols)
                if (opts$verbose) {
                    cat("    chr_stats after merge:", nrow(chr_stats), "rows\n")
                }
                # Fill missing stats with NA (will be filtered out later)
            }
            
            ref_lines <- chr_stats %>%
                mutate(
                    Median = median_transformed,
                    Mean   = mean_transformed,
                    `5th percentile` = q5_transformed,
                    `1st percentile` = q1_transformed,
                    `0.2th percentile` = q0.2_transformed,
                    `95th percentile` = q95_transformed,
                    `99th percentile` = q99_transformed,
                    `99.8th percentile` = q99.8_transformed
                ) %>%
                select(all_of(panel_cols), Median, Mean, `5th percentile`, `1st percentile`, `0.2th percentile`, `95th percentile`, `99th percentile`, `99.8th percentile`) %>%
                tidyr::pivot_longer(
                    cols = c(Median, Mean, `5th percentile`, `1st percentile`, `0.2th percentile`, `95th percentile`, `99th percentile`, `99.8th percentile`),
                    names_to = "label",
                    values_to = "yintercept"
                ) %>%
                filter(!is.na(.data$yintercept)) %>%
                mutate(label = factor(label, levels=c("0.2th percentile","1st percentile","5th percentile","Median","Mean","95th percentile","99th percentile","99.8th percentile"))) %>%
                mutate(color = case_when(
                    label == "Median" ~ "blue",
                    label == "Mean" ~ "red",
                    label == "5th percentile" ~ "lightblue",
                    label == "1st percentile" ~ "cyan",
                    label == "0.2th percentile" ~ "lightcyan",
                    label == "95th percentile" ~ "orange",
                    label == "99th percentile" ~ "purple",
                    label == "99.8th percentile" ~ "darkgreen",
                    TRUE ~ "black"
                ))
            if (opts$verbose) {
                cat("    ref_lines rows:", nrow(ref_lines), "\n")
                cat("    ref_lines columns:", paste(colnames(ref_lines), collapse=", "), "\n")
                if (nrow(ref_lines) > 0 && length(panel_cols) > 0) {
                    cat("    ref_lines panel combinations:\n")
                    print(ref_lines %>% select(all_of(panel_cols), label) %>% distinct())
                }
            }
        } else {
            if (opts$verbose) {
                cat("    No chromosome stats found, creating empty ref_lines\n")
            }
            ref_lines <- data.frame(
                yintercept = numeric(0),
                label = character(0),
                color = character(0),
                stringsAsFactors = FALSE
            )
            if (length(panel_cols) > 0) {
                for (col in panel_cols) {
                    ref_lines[[col]] <- character(0)
                }
            }
        }
    } else {
        # Use genome-wide stats
        genome_stats <- stats_summary %>% filter(.data$chr == "genome")
        if (opts$verbose) {
            cat("  Using genome-wide stats\n")
            cat("    genome_stats rows:", nrow(genome_stats), "\n")
        }
        if (nrow(genome_stats) > 0) {
            # Ensure we have stats for all panel combinations in pair_data
            if (nrow(pair_panels) > 0 && length(panel_cols) > 0) {
                if (opts$verbose) {
                    cat("    Merging genome_stats with pair_panels by:", paste(panel_cols, collapse=", "), "\n")
                    cat("    genome_stats before merge:", nrow(genome_stats), "rows\n")
                }
                # Merge with pair_panels to ensure all combinations are present
                genome_stats <- pair_panels %>%
                    left_join(genome_stats, by=panel_cols)
                if (opts$verbose) {
                    cat("    genome_stats after merge:", nrow(genome_stats), "rows\n")
                }
                # Fill missing stats with NA (will be filtered out later)
            }
            
            ref_lines <- genome_stats %>%
                mutate(
                    Median = median_transformed,
                    Mean   = mean_transformed,
                    `5th percentile` = q5_transformed,
                    `1st percentile` = q1_transformed,
                    `0.2th percentile` = q0.2_transformed,
                    `95th percentile` = q95_transformed,
                    `99th percentile` = q99_transformed,
                    `99.8th percentile` = q99.8_transformed
                ) %>%
                select(all_of(panel_cols), Median, Mean, `5th percentile`, `1st percentile`, `0.2th percentile`, `95th percentile`, `99th percentile`, `99.8th percentile`) %>%
                tidyr::pivot_longer(
                    cols = c(Median, Mean, `5th percentile`, `1st percentile`, `0.2th percentile`, `95th percentile`, `99th percentile`, `99.8th percentile`),
                    names_to = "label",
                    values_to = "yintercept"
                ) %>%
                filter(!is.na(.data$yintercept)) %>%
                mutate(label = factor(label, levels=c("0.2th percentile","1st percentile","5th percentile","Median","Mean","95th percentile","99th percentile","99.8th percentile"))) %>%
                mutate(color = case_when(
                    label == "Median" ~ "blue",
                    label == "Mean" ~ "red",
                    label == "5th percentile" ~ "lightblue",
                    label == "1st percentile" ~ "cyan",
                    label == "0.2th percentile" ~ "lightcyan",
                    label == "95th percentile" ~ "orange",
                    label == "99th percentile" ~ "purple",
                    label == "99.8th percentile" ~ "darkgreen",
                    TRUE ~ "black"
                ))
            if (opts$verbose) {
                cat("    ref_lines rows:", nrow(ref_lines), "\n")
                cat("    ref_lines columns:", paste(colnames(ref_lines), collapse=", "), "\n")
                if (nrow(ref_lines) > 0 && length(panel_cols) > 0) {
                    cat("    ref_lines panel combinations:\n")
                    print(ref_lines %>% select(all_of(panel_cols), label) %>% distinct())
                }
            }
        } else {
            if (opts$verbose) {
                cat("    No genome-wide stats found, creating empty ref_lines\n")
            }
            ref_lines <- data.frame(
                yintercept = numeric(0),
                label = character(0),
                color = character(0),
                stringsAsFactors = FALSE
            )
            if (length(panel_cols) > 0) {
                for (col in panel_cols) {
                    ref_lines[[col]] <- character(0)
                }
            }
        }
    }
    
    if (opts$verbose) {
        cat("  Final ref_lines structure:\n")
        cat("    Rows:", nrow(ref_lines), "\n")
        cat("    Columns:", paste(colnames(ref_lines), collapse=", "), "\n")
        if (nrow(ref_lines) > 0) {
            cat("    Unique labels:", paste(unique(ref_lines$label), collapse=", "), "\n")
        }
    }
    
    # Verify we have data to plot
    if (opts$verbose) {
        cat("  Preparing plot with", nrow(pair_data), "data points\n")
        cat("  cum_pos range:", if("cum_pos" %in% colnames(pair_data)) paste(range(pair_data$cum_pos, na.rm=TRUE), collapse=" to ") else "MISSING", "\n")
        cat("  fst range:", paste(range(pair_data$fst, na.rm=TRUE), collapse=" to "), "\n")
        cat("  Sample of first 5 rows:\n")
        print(head(pair_data %>% select(chr, pos, cum_pos, fst, sample_pair), 5))
    }
    
    if (!"cum_pos" %in% colnames(pair_data)) {
        cat("Error: cum_pos column missing for pair '", pair_name, "'. Cannot create plot.\n", sep="")
        next
    }
    
    if (nrow(pair_data) == 0) {
        cat("Error: No data points to plot for pair '", pair_name, "'\n", sep="")
        next
    }
    
    # Check if all FST values are the same (which would make the line invisible)
    fst_unique <- length(unique(pair_data$fst[!is.na(pair_data$fst)]))
    if (fst_unique == 1) {
        cat("Warning: All FST values are identical (", unique(pair_data$fst[!is.na(pair_data$fst)])[1], ") for pair '", pair_name, "'\n", sep="")
    }
    
    # Create plot
    # Ensure data is sorted by cum_pos (critical for line plotting)
    pair_data <- pair_data %>% arrange(.data$cum_pos)
    
    # Create plot (per-pair mode; overlay-pairs uses a separate code path above)
    p <- ggplot(pair_data, aes(x=cum_pos, y=fst)) +
        # Add chromosome stripes
        geom_rect(data=chr_stripes, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=fill),
                  inherit.aes=FALSE, alpha=0.3) +
        scale_fill_identity() +
        # Add data - use group=1 to ensure single continuous line across all data
        geom_line(alpha=0.7, linewidth=0.5, group=1, na.rm=TRUE) +
        labs(
            x="Genomic Position (cumulative)",
            y=if (transform_this == "asinh") {
                if (is.null(opts$`asinh-scale`)) {
                    "asinh((FST - median) / SD)"
                } else {
                    paste0("asinh((FST - ", round(global_median, 4), ") / ", round(scale_factor, 4), ")")
                }
            } else if (transform_this == "log") {
                "log(FST)"
            } else if (opts$`y-value` == "rank") {
                "FST rank"
            } else if (opts$`y-value` == "quantile") {
                "FST quantile"
            } else {
                "FST (Fixation Index)"
            },
            title=paste("FST:", pair_name)
        ) +
        theme_bw(base_size = PLOT_BASE_SIZE, base_family = "sans") +
        theme(
            panel.grid.minor=element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1)
        )
    if (opts$`plot-style` == "line_points") {
        p <- p + geom_point(size = 1, alpha = 0.7, na.rm = TRUE)
    }

    # Add y-axis scale with inverse transformation for labels (if transformation is applied)
    # This labels the axis ticks with original (untransformed) values while plotting transformed data
    # Helper function to format numbers with scientific notation when appropriate
    format_tick_label <- function(x) {
        # Format each number individually
        formatted <- sapply(x, function(val) {
            if (val == 0) {
                return("0")
            }
            # Use scientific notation for very small (< 0.001) or very large (> 1000) numbers
            # or if the fixed format would be too long
            abs_val <- abs(val)
            if (abs_val < 1e-3 || abs_val > 1e3) {
                # Use 'e' format for scientific notation (e.g., 2e-4, 1.5e2)
                return(formatC(val, format="e", digits=2))
            } else {
                # Use fixed format for moderate values
                # Check if fixed format would be reasonable length
                fixed_fmt <- formatC(val, format="f", digits=4)
                if (nchar(fixed_fmt) > 8) {
                    # Too long, use scientific
                    return(formatC(val, format="e", digits=2))
                } else {
                    # Remove trailing zeros
                    return(sub("\\.?0+$", "", fixed_fmt))
                }
            }
        })
        return(formatted)
    }
    
    if (transform_this == "asinh") {
        # Get range of original (untransformed) values
        orig_range <- range(pair_data_original$fst, na.rm=TRUE)
        
        # Generate nice tick values in original space, ensuring 0 is included
        # Expand range slightly to ensure we include 0 if it's near the range
        expanded_min <- min(orig_range[1], 0)
        expanded_max <- max(orig_range[2], 0)
        
        # Generate nice breaks in original space
        orig_ticks <- pretty(c(expanded_min, expanded_max), n=8)
        
        # Ensure 0 is included (add it if not already present)
        if (!0 %in% orig_ticks) {
            orig_ticks <- sort(c(orig_ticks, 0))
        }
        
        # Transform these original values to get tick positions in transformed space
        # asinh transformation: y = asinh((x - median) / scale)
        y_breaks <- asinh((orig_ticks - global_median) / scale_factor)
        
        # Format labels with scientific notation when appropriate
        y_labels_formatted <- format_tick_label(orig_ticks)
        
        p <- p + scale_y_continuous(
            breaks = y_breaks,
            labels = y_labels_formatted
        )
    } else if (transform_this == "log") {
        # Get range of original (untransformed) values
        orig_range <- range(pair_data_original$fst, na.rm=TRUE)
        
        # For log transformation, we need to handle the shift
        if (!is.null(log_shift) && log_shift > 0) {
            # Data was shifted: log1p(x + shift)
            # Original range is already in the shifted space
            expanded_min <- min(orig_range[1], 0)
            expanded_max <- max(orig_range[2], 0)
        } else {
            # No shift: log(x), so values must be positive
            expanded_min <- max(orig_range[1], 1e-10)  # Avoid log(0)
            expanded_max <- orig_range[2]
        }
        
        # Generate nice breaks in original space
        orig_ticks <- pretty(c(expanded_min, expanded_max), n=8)
        
        # Ensure 0 is included if possible (only if shift allows it)
        if (!is.null(log_shift) && log_shift > 0) {
            # With shift, we can include 0
            if (!0 %in% orig_ticks) {
                orig_ticks <- sort(c(orig_ticks, 0))
            }
        } else {
            # For log without shift, we can't have 0, but ensure we have a small positive value
            # pretty() should already give us reasonable values, but ensure minimum is reasonable
            if (min(orig_ticks) > 1e-6) {
                # If all ticks are > 1e-6, add a small one near the minimum
                orig_ticks <- sort(c(orig_ticks, max(1e-10, expanded_min * 0.1)))
            }
        }
        
        # Transform these original values to get tick positions in transformed space
        if (!is.null(log_shift) && log_shift > 0) {
            # log1p(x + shift): y = log1p(x + shift)
            y_breaks <- log1p(orig_ticks + log_shift)
        } else {
            # log(x): y = log(x), but ensure x > 0
            orig_ticks <- pmax(orig_ticks, 1e-10)
            y_breaks <- log(orig_ticks)
        }
        
        # Format labels with scientific notation when appropriate
        y_labels_formatted <- format_tick_label(orig_ticks)
        
        p <- p + scale_y_continuous(
            breaks = y_breaks,
            labels = y_labels_formatted
        )
    }
    
    # Add reference lines (genome-wide) with legend
    # Only use labels that are actually present in ref_lines (not all factor levels)
    if (nrow(ref_lines) > 0) {
        # Get unique labels that are actually present (not all factor levels)
        present_labels <- unique(ref_lines$label)
        present_labels <- present_labels[!is.na(present_labels)]
        
        if (opts$verbose) {
            cat("  Creating reference lines with", length(present_labels), "labels:", paste(present_labels, collapse=", "), "\n")
        }
        
        # Create color and linetype mappings only for present labels
        color_map <- setNames(ref_lines$color, ref_lines$label)[as.character(present_labels)]
        linetype_map <- setNames(rep("dashed", length(present_labels)), as.character(present_labels))
        
        if (opts$verbose) {
            cat("  color_map length:", length(color_map), "\n")
            cat("  linetype_map length:", length(linetype_map), "\n")
            cat("  present_labels length:", length(present_labels), "\n")
        }
        
        p <- p +
            geom_hline(data=ref_lines, 
                       aes(yintercept=.data$yintercept, linetype=.data$label, color=.data$label),
                       linewidth=0.8, show.legend=TRUE) +
            scale_linetype_manual(
                name="Reference",
                values=linetype_map,
                breaks=as.character(present_labels),
                guide=guide_legend(override.aes=list(color=color_map))
            ) +
            scale_color_manual(
                name="Reference",
                values=color_map,
                breaks=as.character(present_labels),
                guide="legend"
            )
    } else {
        if (opts$verbose) {
            cat("  No reference lines to add (ref_lines is empty)\n")
        }
    }
    
    # Add paneling
    # Check if we have multiple window sizes and multiple pairs
    has_window_size <- "window_size" %in% colnames(pair_data) && !all(is.na(pair_data$window_size))
    if (has_window_size) {
        unique_windows <- unique(pair_data$window_size[!is.na(pair_data$window_size)])
        multiple_windows <- length(unique_windows) > 1
    } else {
        multiple_windows <- FALSE
    }
    unique_pairs <- unique(pair_data$sample_pair)
    multiple_pairs <- length(unique_pairs) > 1
    
    if (opts$`panel-by` == "window" && has_window_size) {
        # Panel by window size vertically
        if (multiple_pairs) {
            # If multiple pairs, add horizontal paneling by pair
            p <- p + facet_grid(window_size ~ sample_pair, scales="fixed")
        } else {
            # Single pair, just panel by window size
            p <- p + facet_wrap(~ window_size, scales="fixed", ncol=1)
        }
    } else if (opts$`panel-by` == "pair" && multiple_pairs) {
        # Panel by pair vertically
        if (has_window_size && multiple_windows) {
            # If multiple window sizes, add horizontal paneling by window size
            p <- p + facet_grid(sample_pair ~ window_size, scales="fixed")
        } else {
            # Single window size or no window size, just panel by pair
            p <- p + facet_wrap(~ sample_pair, scales="fixed", ncol=1)
        }
    } else if (opts$`panel-by` == "both" && has_window_size) {
        # Explicitly panel by both: pairs as rows, window sizes as columns
        p <- p + facet_grid(sample_pair ~ window_size, scales="fixed")
    }
    # If panel-by == "none", no paneling is added
    
    # Add chromosome boundary labels (ordered by length, longest first)
    # Get chromosome labels in the correct order (longest first) - only those in data
    chr_labels <- chr_lengths$chr[chr_lengths$chr %in% unique(pair_data$chr)]
    
    # Filter chr_stripes to only include chromosomes present in the data
    # Match by checking which chromosomes from chr_lengths are in the data
    chr_centers <- chr_stripes %>%
        # Create a mapping: find which stripe corresponds to which chromosome
        # chr_stripes are in the same order as chr_lengths (longest first)
        mutate(chr_idx = row_number()) %>%
        filter(chr_idx <= length(chr_labels)) %>%
        mutate(center = (.data$xmin + .data$xmax) / 2) %>%
        select(center)
    
    # Ensure breaks and labels have the same length
    if (nrow(chr_centers) == length(chr_labels) && nrow(chr_centers) > 0) {
    p <- p + scale_x_continuous(
        breaks=chr_centers$center,
            labels=chr_labels,
        expand=expansion(mult=0.01)
    )
    } else if (length(chr_labels) > 0) {
        # Fallback: use subset if lengths don't match
        min_len <- min(nrow(chr_centers), length(chr_labels))
        if (min_len > 0) {
            p <- p + scale_x_continuous(
                breaks=chr_centers$center[1:min_len],
                labels=chr_labels[1:min_len],
        expand=expansion(mult=0.01)
    )
        } else {
            p <- p + scale_x_continuous(expand=expansion(mult=0.01))
        }
    } else {
        # No chromosome breaks if no data
        p <- p + scale_x_continuous(expand=expansion(mult=0.01))
    }
    
    # Save plot
    prefix <- if (opts$`file-prefix` != "") paste0(opts$`file-prefix`, "_") else ""
    pair_safe <- gsub(":", "_", pair_name)
    pair_safe <- gsub("[^A-Za-z0-9_]", "_", pair_safe)
    # Truncate pair name if too long to avoid filename length issues
    if (nchar(pair_safe) > 80) {
        pair_safe <- substr(pair_safe, 1, 80)
    }
    
    # Build filename suffixes with actual values separated by dashes
    # Window size suffix - always include
    if (!is.null(selected_window_size)) {
        window_suffix <- paste0("_w", selected_window_size)
    } else if (has_window_size) {
        # List all window sizes separated by dashes
        unique_windows <- sort(unique(pair_data$window_size[!is.na(pair_data$window_size)]))
        if (length(unique_windows) > 0) {
            window_suffix <- paste0("_w", paste(unique_windows, collapse="-"))
        } else {
            window_suffix <- ""
        }
    } else {
        window_suffix <- ""
    }
    
    # Chromosome suffix - always include
    if (!is.null(selected_chr)) {
        chr_safe <- gsub("[^A-Za-z0-9]", "_", selected_chr)
        chr_suffix <- paste0("_chr", chr_safe)
    } else {
        # List all chromosomes separated by dashes
        unique_chrs <- unique(pair_data$chr[!is.na(pair_data$chr)])
        unique_chrs <- unique_chrs[unique_chrs != "genome"]  # Exclude genome-wide
        if (length(unique_chrs) > 0) {
            # Sort by order in chr_lengths (longest first)
            chr_order <- match(unique_chrs, chr_lengths$chr)
            unique_chrs <- unique_chrs[order(chr_order)]
            # If more than 3, just show first one and count
            if (length(unique_chrs) <= 3) {
                chr_list <- paste(gsub("[^A-Za-z0-9]", "_", unique_chrs), collapse="-")
            } else {
                chr_list <- paste0(gsub("[^A-Za-z0-9]", "_", unique_chrs[1]), "-and", length(unique_chrs)-1, "more")
            }
            chr_suffix <- paste0("_chr", chr_list)
        } else {
            chr_suffix <- "_chrAll"
        }
    }
    
    # Transform suffix - add if transform applied
    transform_suffix <- ""
    if (transform_this != "none") {
        if (transform_this == "asinh") transform_suffix <- "_asinhtrans"
        else if (transform_this == "log") transform_suffix <- "_logtrans"
    }
    if (opts$`y-value` == "rank") transform_suffix <- paste0(transform_suffix, "_rank")
    if (opts$`y-value` == "quantile") transform_suffix <- paste0(transform_suffix, "_quantile")

    if ("png" %in% format_parts) {
        png_file <- file.path(opts$`output-dir`, paste0(prefix, "fst_", pair_safe, window_suffix, chr_suffix, transform_suffix, ".png"))
        ggsave(png_file, p, width=opts$width, height=opts$height, dpi=if (!is.null(opts$dpi)) opts$dpi else PLOT_DPI)
        cat("Saved:", png_file, "\n")
    }
    
    if ("pdf" %in% format_parts) {
        pdf_file <- file.path(opts$`output-dir`, paste0(prefix, "fst_", pair_safe, window_suffix, chr_suffix, transform_suffix, ".pdf"))
        ggsave(pdf_file, p, width=opts$width, height=opts$height)
        cat("Saved:", pdf_file, "\n")
    }
    
    if ("svg" %in% format_parts) {
        svg_file <- file.path(opts$`output-dir`, paste0(prefix, "fst_", pair_safe, window_suffix, chr_suffix, transform_suffix, ".svg"))
        ggsave(svg_file, p, width=opts$width, height=opts$height, device="svg")
        cat("Saved:", svg_file, "\n")
    }
}
}

# Combine and save statistics
if (length(stats_results) > 0) {
    all_stats <- bind_rows(stats_results) %>%
        select(sample_pair, window_size, chr, median, mean, sd,
               q5, q1, q0.2, q95, q99, q99.8) %>%
        arrange(sample_pair, window_size, chr)
    
    prefix <- if (opts$`file-prefix` != "") paste0(opts$`file-prefix`, "_") else ""
    
    # Build filename suffixes for statistics file (same logic as plot files)
    # Get unique values from all_stats
    unique_windows_stats <- if ("window_size" %in% colnames(all_stats)) {
        sort(unique(all_stats$window_size[!is.na(all_stats$window_size)]))
    } else {
        numeric(0)
    }
    
    unique_chrs_stats <- if ("chr" %in% colnames(all_stats)) {
        unique(all_stats$chr[!is.na(all_stats$chr) & all_stats$chr != "genome"])
    } else {
        character(0)
    }
    
    unique_pairs_stats <- if ("sample_pair" %in% colnames(all_stats)) {
        sort(unique(all_stats$sample_pair[!is.na(all_stats$sample_pair)]))
    } else {
        character(0)
    }
    
    # Window size suffix
    if (!is.null(selected_window_size)) {
        window_suffix <- paste0("_w", selected_window_size)
    } else if (length(unique_windows_stats) > 0) {
        window_suffix <- paste0("_w", paste(unique_windows_stats, collapse="-"))
    } else {
        window_suffix <- ""
    }
    
    # Chromosome suffix
    if (!is.null(selected_chr)) {
        chr_safe <- gsub("[^A-Za-z0-9]", "_", selected_chr)
        chr_suffix <- paste0("_chr", chr_safe)
    } else if (length(unique_chrs_stats) > 0) {
        # Sort by order in chr_lengths
        chr_order <- match(unique_chrs_stats, chr_lengths$chr)
        unique_chrs_stats <- unique_chrs_stats[order(chr_order)]
        # If more than 3, just show first one and count
        if (length(unique_chrs_stats) <= 3) {
            chr_list <- paste(gsub("[^A-Za-z0-9]", "_", unique_chrs_stats), collapse="-")
        } else {
            chr_list <- paste0(gsub("[^A-Za-z0-9]", "_", unique_chrs_stats[1]), "-and", length(unique_chrs_stats)-1, "more")
        }
        chr_suffix <- paste0("_chr", chr_list)
    } else {
        chr_suffix <- "_chrAll"
    }
    
    # Pair suffix (for FST, we have pairs instead of samples)
    if (length(unique_pairs_stats) > 1) {
        # Clean pair names
        pair_names_clean <- gsub(":", "_", unique_pairs_stats)
        pair_names_clean <- gsub("[^A-Za-z0-9_]", "_", pair_names_clean)
        
        # If more than 3, just show first one and count
        if (length(pair_names_clean) > 3) {
            pair_list <- paste0(pair_names_clean[1], "-and", length(pair_names_clean)-1, "more")
            pair_suffix <- paste0("_pairs", pair_list)
        } else {
            pair_list <- paste(pair_names_clean, collapse="-")
            pair_suffix <- paste0("_pairs", pair_list)
        }
    } else if (length(unique_pairs_stats) == 1) {
        pair_safe <- gsub(":", "_", unique_pairs_stats[1])
        pair_safe <- gsub("[^A-Za-z0-9_]", "_", pair_safe)
        # Truncate single pair name if too long
        if (nchar(pair_safe) > 50) {
            pair_safe <- substr(pair_safe, 1, 50)
        }
        pair_suffix <- paste0("_pair", pair_safe)
    } else {
        pair_suffix <- ""
    }
    
    # Transform suffix for statistics file
    transform_suffix_stats <- ""
    if (opts$transform != "none") {
        if (opts$transform == "asinh") {
            transform_suffix_stats <- "_asinhtrans"
        } else if (opts$transform == "log") {
            transform_suffix_stats <- "_logtrans"
        }
    }
    
    stats_file <- file.path(opts$`output-dir`, paste0(prefix, "fst_statistics", window_suffix, chr_suffix, pair_suffix, transform_suffix_stats, ".csv"))
    write_csv(all_stats, stats_file)
    cat("\nSaved statistics to:", stats_file, "\n")
} else {
    cat("\nWarning: No statistics were calculated\n")
}

cat("\nDone!\n")
