#!/usr/bin/env Rscript

###############################################################################
# plot_diversity.R
#
# Create ggplot2 plots of diversity statistics (π, θ, Tajima's D) from
# grenedalf diversity output files.
#
# Features:
# - Chromosomal positions on x-axis with alternating white/grey stripes
# - Paneling by window size, sample name, or both
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

# Parse command-line arguments
option_list <- list(
    make_option(c("--input-dir"), type="character", default=NULL,
                help="Directory containing diversity TSV files", metavar="DIR"),
    make_option(c("--output-dir"), type="character", default=".",
                help="Output directory for plots [default: %default]", metavar="DIR"),
    make_option(c("--reference-genome"), type="character", default=NULL,
                help="Reference genome FASTA file (optional, for chromosome lengths)", metavar="FILE"),
    make_option(c("--panel-by"), type="character", default="both",
                help="Paneling option: window, sample, both, or none [default: %default]",
                metavar="OPTION"),
    make_option(c("--statistics"), type="character", default="pi,theta,tajima_d",
                help="Comma-separated list of statistics to plot [default: %default]",
                metavar="LIST"),
    make_option(c("--plot-format"), type="character", default="png",
                help="Plot format: png, pdf, svg, or any combination (e.g., both, all) [default: %default]", metavar="FORMAT"),
    make_option(c("--width"), type="numeric", default=12,
                help="Plot width in inches [default: %default]", metavar="NUMBER"),
    make_option(c("--height"), type="numeric", default=8,
                help="Plot height in inches [default: %default]", metavar="NUMBER"),
    make_option(c("--file-prefix"), type="character", default="",
                help="Prefix for output files [default: no prefix]", metavar="PREFIX"),
    make_option(c("--chromosome"), type="character", default=NULL,
                help="Single chromosome/scaffold to plot (default: all chromosomes)", metavar="CHR"),
    make_option(c("--window-size"), type="numeric", default=NULL,
                help="Single window size to plot (default: all window sizes)", metavar="NUMBER"),
    make_option(c("--top-n-chromosomes"), type="numeric", default=NULL,
                help="Plot only the N longest chromosomes/scaffolds (optional)", metavar="N"),
    make_option(c("--min-chromosome-length"), type="numeric", default=NULL,
                help="Plot only chromosomes/scaffolds longer than N bp (optional)", metavar="N"),
    make_option(c("--transform"), type="character", default="none",
                help="Y-axis transformation: none, asinh (centered on median), or log [default: %default]", metavar="TRANSFORM"),
    make_option(c("--asinh-scale"), type="numeric", default=NULL,
                help="Scale factor for asinh transformation (default: use global standard deviation)", metavar="NUMBER"),
    make_option(c("--verbose"), action="store_true", default=FALSE,
                help="Enable verbose output for debugging", metavar="FLAG")
)

opt_parser <- OptionParser(option_list=option_list, usage="usage: %prog [options]")
opts <- parse_args(opt_parser)

# Validate arguments
if (is.null(opts$`input-dir`)) {
    stop("--input-dir is required")
}

if (!opts$`panel-by` %in% c("window", "sample", "both", "none")) {
    stop("--panel-by must be one of: window, sample, both, none")
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

# Parse statistics list
statistics <- strsplit(opts$statistics, ",")[[1]]
statistics <- trimws(statistics)
valid_stats <- c("pi", "theta", "tajima_d")
if (!all(statistics %in% valid_stats)) {
    stop("Invalid statistics. Valid options: ", paste(valid_stats, collapse=", "))
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

# Function to calculate cumulative positions
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

# Helper function to extract window size from filename
extract_window_size_from_filename <- function(filename) {
    window_match <- regmatches(filename, regexpr("w(\\d+)_s(\\d+)", filename))
    if (length(window_match) > 0) {
        window_parts <- strsplit(window_match, "_")[[1]]
        window_size <- as.numeric(sub("w", "", window_parts[1]))
        return(window_size)
    }
    return(NA)
}

# Helper function to extract sample name from file path/filename
extract_sample_from_filepath <- function(filepath, input_dir) {
    # Try directory structure first
    file_dir <- normalizePath(dirname(filepath), winslash="/", mustWork=FALSE)
    input_dir_norm <- normalizePath(input_dir, winslash="/", mustWork=FALSE)
    
    dir_parts <- strsplit(file_dir, "/")[[1]]
    input_dir_parts <- strsplit(input_dir_norm, "/")[[1]]
    
    if (length(dir_parts) > length(input_dir_parts)) {
        # File is in a subdirectory of input_dir
        sample_idx <- length(input_dir_parts) + 1
        if (sample_idx <= length(dir_parts)) {
            return(dir_parts[sample_idx])
        }
    }
    
    # Try filename
    filename_no_ext <- tools::file_path_sans_ext(basename(filepath))
    if (grepl("^([A-Za-z0-9_]+)_diversity", filename_no_ext, ignore.case=TRUE)) {
        return(sub("^([A-Za-z0-9_]+)_diversity.*", "\\1", filename_no_ext, ignore.case=TRUE))
    }
    
    return(NA)
}

# Find all diversity files (TSV or CSV)
input_dir <- opts$`input-dir`
tsv_files <- list.files(input_dir, pattern=".*diversity.*\\.tsv$", full.names=TRUE, recursive=TRUE)
csv_files <- list.files(input_dir, pattern=".*diversity.*\\.csv$", full.names=TRUE, recursive=TRUE)
all_files <- c(tsv_files, csv_files)

if (length(all_files) == 0) {
    stop("No diversity files found in: ", input_dir, " (expected files matching pattern: *diversity*.tsv or *diversity*.csv)")
}

if (opts$verbose) {
    cat("Found", length(all_files), "diversity file(s) (", length(tsv_files), " TSV, ", length(csv_files), " CSV)\n", sep="")
}

# Filter files by window size if specified (before reading)
selected_window_size <- opts$`window-size`
if (!is.null(selected_window_size)) {
    file_window_sizes <- sapply(all_files, extract_window_size_from_filename)
    # Check for "single" label in filenames (single-position windows, 1bp windows at each SNP)
    has_single_label <- grepl("single", basename(all_files), ignore.case=TRUE)
    # Keep files that:
    # - Do not have "single" label, AND
    # - Have no window size indication in filename (is.na), OR match the selected window size
    # Exclude files with "single" label or different window sizes
    all_files <- all_files[!has_single_label & (is.na(file_window_sizes) | file_window_sizes == selected_window_size)]
    if (opts$verbose) {
        cat("Filtered to", length(all_files), "file(s) matching window size", selected_window_size, " (excluding *single* and different window sizes)\n")
    }
    if (length(all_files) == 0) {
        stop("No files found matching window size ", selected_window_size)
    }
}

# Read and combine all diversity files using tidyverse
all_data <- list()
for (file in all_files) {
    cat("\n=== Reading:", basename(file), "===\n")
    
    # Detect separator by reading first line
    first_line <- readLines(file, n=1)
    has_tabs <- grepl("\t", first_line)
    has_commas <- grepl(",", first_line)
    
    # Determine separator and file type
    if (has_tabs && (!has_commas || length(strsplit(first_line, "\t")[[1]]) > length(strsplit(first_line, ",")[[1]]))) {
        # Tab-separated (more tabs than commas, or only tabs)
        data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        detected_sep <- "tab"
        cat("  Detected separator: tab\n")
    } else if (has_commas) {
        # Comma-separated
        data <- read_csv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        detected_sep <- "comma"
        cat("  Detected separator: comma\n")
    } else {
        # Default to TSV if neither detected (shouldn't happen, but fallback)
        data <- read_tsv(file, col_types=cols(.default="c"), show_col_types=FALSE, progress=FALSE)
        detected_sep <- "tab"
        cat("  Detected separator: tab (fallback)\n")
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
    
    # Convert statistics columns to numeric (sample-prefixed columns)
    stat_cols <- grep("\\.(theta_pi|theta_watterson|tajimas?_d)$", colnames(data), ignore.case=FALSE, value=TRUE)
    if (length(stat_cols) > 0) {
        if (opts$verbose) cat("  Found statistics columns:", paste(stat_cols, collapse=", "), "\n")
        data <- data %>%
            mutate(across(all_of(stat_cols), as.numeric))
    } else {
        cat("  Warning: No statistics columns found matching pattern *.theta_pi, *.theta_watterson, *.tajimas_d\n")
        if (opts$verbose) cat("  All columns:", paste(colnames(data), collapse=", "), "\n")
    }
    
    # Handle different column name variations - ensure chr and pos exist using tidyverse
    # Check for chromosome column (various names) - MUST be character, not numeric!
    if (!"chr" %in% colnames(data)) {
        if ("chromosome" %in% colnames(data)) {
            data <- data %>% mutate(chr = as.character(chromosome))
            if (opts$verbose) cat("  Mapped 'chromosome' -> 'chr' (as character)\n")
        } else if ("chrom" %in% colnames(data)) {
            # chrom might have been read as character, ensure it stays character
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
    
    # Extract sample name and window size from filename first
    filename <- basename(file)
    
    # Try to extract window size from filename (e.g., diversity_w10000_s5000.csv)
    window_match <- regmatches(filename, regexpr("w(\\d+)_s(\\d+)", filename))
    if (length(window_match) > 0) {
        window_parts <- strsplit(window_match, "_")[[1]]
        window_size <- as.numeric(sub("w", "", window_parts[1]))
        step_size <- as.numeric(sub("s", "", window_parts[2]))
        data$window_size <- window_size
        data$step_size <- step_size
    } else {
        # Try to extract from directory structure or use default
        data$window_size <- NA
        data$step_size <- NA
    }
    
    # Try to extract sample name from directory or filename
    # Check if file is in a sample-specific directory
    # Normalize path separators (handle both / and \)
    file_dir <- normalizePath(dirname(file), winslash="/", mustWork=FALSE)
    input_dir_norm <- normalizePath(input_dir, winslash="/", mustWork=FALSE)
    
    dir_parts <- strsplit(file_dir, "/")[[1]]
    input_dir_parts <- strsplit(input_dir_norm, "/")[[1]]
    
    sample_name <- "unknown"
    if (length(dir_parts) > length(input_dir_parts)) {
        # File is in a subdirectory of input_dir
        # Sample name should be the first subdirectory level
        sample_idx <- length(input_dir_parts) + 1
        if (sample_idx <= length(dir_parts)) {
            sample_name <- dir_parts[sample_idx]
        }
    } else {
        # File is directly in input_dir, try to extract from filename
        filename_no_ext <- tools::file_path_sans_ext(basename(file))
        # Try common patterns like "sample_diversity..." or "diversity_sample..."
        if (grepl("^([A-Za-z0-9_]+)_diversity", filename_no_ext, ignore.case=TRUE)) {
            sample_name <- sub("^([A-Za-z0-9_]+)_diversity.*", "\\1", filename_no_ext, ignore.case=TRUE)
        }
    }
    
    # Extract diversity statistics from sample-prefixed columns using tidyverse
    # Look for columns matching patterns like: sample.theta_pi, sample.theta_watterson, sample.tajimas_d
    # Grenedalf format: sample_name.theta_pi, sample_name.theta_watterson, sample_name.tajimas_d
    
    # Find theta_pi columns (π / pi)
    theta_pi_cols <- grep("\\.theta_pi$", colnames(data), ignore.case=FALSE, value=TRUE)
    if (opts$verbose) cat("  Searching for theta_pi columns...\n")
    if (length(theta_pi_cols) > 0) {
        if (opts$verbose) cat("    Found theta_pi columns:", paste(theta_pi_cols, collapse=", "), "\n")
        # Use the first theta_pi column found
        theta_pi_col <- theta_pi_cols[1]
        data <- data %>%
            mutate(pi = as.numeric(.data[[theta_pi_col]]))
        if (opts$verbose) {
            cat("    Extracted pi from column:", theta_pi_col, "\n")
            cat("    Sample pi values (first 5):", paste(head(data$pi, 5), collapse=", "), "\n")
            cat("    Non-NA pi values:", sum(!is.na(data$pi)), "out of", nrow(data), "\n")
        }
        
        # Extract sample name from column name if not already set from directory
        if (sample_name == "unknown") {
            sample_from_col <- sub("\\.theta_pi$", "", theta_pi_col, ignore.case=FALSE)
            # Remove common suffixes
            sample_from_col <- sub("_all_seq\\.dedup$", "", sample_from_col, ignore.case=TRUE)
            sample_from_col <- sub("\\.dedup$", "", sample_from_col, ignore.case=TRUE)
            sample_from_col <- sub("_all_seq$", "", sample_from_col, ignore.case=TRUE)
            # Try to extract shorter sample name (e.g., "Cheney" from "cheney_pool_s3")
            if (grepl("^([a-z0-9_]+)_pool_s[0-9]+", sample_from_col, ignore.case=TRUE)) {
                sample_name <- sub("^([a-z0-9_]+)_pool_s[0-9]+.*", "\\1", sample_from_col, ignore.case=TRUE)
            } else {
                sample_name <- sample_from_col
            }
            if (opts$verbose) cat("    Extracted sample name from column:", sample_name, "\n")
        }
    } else {
        if (opts$verbose) cat("    WARNING: No theta_pi columns found\n")
    }
    
    # Find theta_watterson columns (θ / theta)
    theta_watterson_cols <- grep("\\.theta_watterson$", colnames(data), ignore.case=FALSE, value=TRUE)
    if (opts$verbose) cat("  Searching for theta_watterson columns...\n")
    if (length(theta_watterson_cols) > 0) {
        if (opts$verbose) cat("    Found theta_watterson columns:", paste(theta_watterson_cols, collapse=", "), "\n")
        theta_watterson_col <- theta_watterson_cols[1]
        data <- data %>%
            mutate(theta = as.numeric(.data[[theta_watterson_col]]))
        if (opts$verbose) {
            cat("    Extracted theta from column:", theta_watterson_col, "\n")
            cat("    Non-NA theta values:", sum(!is.na(data$theta)), "out of", nrow(data), "\n")
        }
    } else {
        if (opts$verbose) cat("    WARNING: No theta_watterson columns found\n")
    }
    
    # Find tajimas_d or tajima_d columns (Tajima's D)
    # Note: grenedalf outputs "tajimas_d" (with 's') in the format: sample.tajimas_d
    tajima_cols <- grep("\\.tajimas?_d$", colnames(data), ignore.case=FALSE, value=TRUE)
    if (opts$verbose) cat("  Searching for tajima columns...\n")
    if (length(tajima_cols) > 0) {
        if (opts$verbose) cat("    Found tajima columns:", paste(tajima_cols, collapse=", "), "\n")
        tajima_col <- tajima_cols[1]
        data <- data %>%
            mutate(tajima_d = as.numeric(.data[[tajima_col]]))
        if (opts$verbose) {
            cat("    Extracted tajima_d from column:", tajima_col, "\n")
            cat("    Non-NA tajima_d values:", sum(!is.na(data$tajima_d)), "out of", nrow(data), "\n")
        }
    } else if ("tajima-d" %in% colnames(data)) {
        data <- data %>% mutate(tajima_d = as.numeric(`tajima-d`))
        if (opts$verbose) cat("    Found 'tajima-d' column\n")
    } else if ("tajima_d" %in% colnames(data)) {
        data <- data %>% mutate(tajima_d = as.numeric(tajima_d))
        if (opts$verbose) cat("    Found 'tajima_d' column\n")
    } else if ("tajimas_d" %in% colnames(data)) {
        data <- data %>% mutate(tajima_d = as.numeric(tajimas_d))
        if (opts$verbose) cat("    Found 'tajimas_d' column\n")
    } else {
        if (opts$verbose) cat("    WARNING: No tajima columns found\n")
    }
    
    # Set sample name
    data <- data %>% mutate(sample = sample_name)
    if (opts$verbose) {
        cat("  Sample name:", sample_name, "\n")
        cat("  Window size:", if(exists("window_size")) window_size else "NA", "\n")
        cat("  Rows:", nrow(data), "\n")
    }
    
    all_data[[file]] <- data
}

# Combine all data using tidyverse
if (opts$verbose) {
    cat("\n=== Combining data from", length(all_data), "file(s) ===\n")
}
combined_data <- bind_rows(all_data)

if (opts$verbose) {
    cat("Combined data dimensions:", nrow(combined_data), "rows,", ncol(combined_data), "columns\n")
    cat("Combined columns:", paste(colnames(combined_data), collapse=", "), "\n")
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
}

# Check statistics columns
stat_cols_found <- c()
if ("pi" %in% colnames(combined_data)) {
    stat_cols_found <- c(stat_cols_found, "pi")
    if (opts$verbose) {
        cat("Found 'pi' column with", sum(!is.na(combined_data$pi)), "non-NA values out of", nrow(combined_data), "\n")
        cat("  pi type:", class(combined_data$pi), "\n")
        cat("  pi sample values (first 5):", paste(head(combined_data$pi[!is.na(combined_data$pi)], 5), collapse=", "), "\n")
    }
}
if ("theta" %in% colnames(combined_data)) {
    stat_cols_found <- c(stat_cols_found, "theta")
    if (opts$verbose) cat("Found 'theta' column with", sum(!is.na(combined_data$theta)), "non-NA values out of", nrow(combined_data), "\n")
}
if ("tajima_d" %in% colnames(combined_data)) {
    stat_cols_found <- c(stat_cols_found, "tajima_d")
    if (opts$verbose) cat("Found 'tajima_d' column with", sum(!is.na(combined_data$tajima_d)), "non-NA values out of", nrow(combined_data), "\n")
}

if (length(stat_cols_found) == 0) {
    cat("ERROR: No statistics columns found in combined data!\n")
    cat("Available columns:", paste(colnames(combined_data), collapse=", "), "\n")
}

# Remove rows with missing chromosome or position
if (opts$verbose) {
    cat("\nFiltering rows with missing chr or pos...\n")
    cat("Before filtering:", nrow(combined_data), "rows\n")
}
combined_data <- combined_data %>%
    filter(!is.na(.data$chr), !is.na(.data$pos))
if (opts$verbose) cat("After filtering:", nrow(combined_data), "rows\n")

if (nrow(combined_data) == 0) {
    stop("Error: All rows were filtered out! Check chr and pos columns.")
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

# Filter to specific window size if requested
selected_window_size <- opts$`window-size`
if (!is.null(selected_window_size)) {
    if (opts$verbose) cat("\nFiltering to window size:", selected_window_size, "\n")
    if ("window_size" %in% colnames(combined_data)) {
        available_windows <- unique(combined_data$window_size[!is.na(combined_data$window_size)])
        if (!selected_window_size %in% available_windows) {
            stop("Error: Window size '", selected_window_size, "' not found in data. Available window sizes: ", 
                 paste(available_windows, collapse=", "))
        }
        combined_data <- combined_data %>%
            filter(is.na(.data$window_size) | .data$window_size == selected_window_size)
        if (opts$verbose) cat("Filtered data:", nrow(combined_data), "rows\n")
    } else {
        if (opts$verbose) cat("Warning: No window_size column found, ignoring --window-size filter\n")
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

# Calculate statistics for each statistic type
stats_results <- list()

for (stat in statistics) {
    if (opts$verbose) cat("\nProcessing statistic:", stat, "\n")
    
    # Map statistic name to column name
    stat_col <- switch(stat,
        "pi" = "pi",
        "theta" = "theta",
        "tajima_d" = "tajima_d"
    )
    
    # Check if column exists
    if (!stat_col %in% colnames(combined_data)) {
        cat("Warning: Column '", stat_col, "' not found in combined data\n", sep="")
        cat("  Available columns:", paste(colnames(combined_data), collapse=", "), "\n")
        cat("  Skipping statistic:", stat, "\n")
        next
    }
    
    # Check if column has any non-NA values
    if (all(is.na(combined_data[[stat_col]]))) {
        cat("Warning: Column '", stat_col, "' contains only NA values, skipping\n", sep="")
        next
    }
    
    # Filter out NA values for this statistic
    stat_data <- combined_data %>%
        filter(!is.na(.data[[stat_col]]))
    
    if (nrow(stat_data) == 0) {
        cat("Warning: No valid data for statistic '", stat, "', skipping\n", sep="")
        next
    }
    
    # Store original data for statistics calculation
    stat_data_original <- stat_data
    
    # Calculate global median and standard deviation for asinh transformation (if needed)
    # Also store log shift if needed for axis label inverse transformation
    global_median <- NULL
    scale_factor <- NULL
    log_shift <- NULL
    if (opts$transform == "asinh") {
        global_median <- median(stat_data_original[[stat_col]], na.rm=TRUE)
        global_sd <- sd(stat_data_original[[stat_col]], na.rm=TRUE)
        
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
    if (opts$transform != "none") {
        if (opts$transform == "asinh") {
            # Apply asinh transformation: asinh((x - global_median) / scale_factor)
            stat_data <- stat_data %>%
                mutate(!!stat_col := asinh((.data[[stat_col]] - global_median) / scale_factor))
        } else if (opts$transform == "log") {
            # Check for zeros or negatives
            min_val <- min(stat_data_original[[stat_col]], na.rm=TRUE)
            if (min_val <= 0) {
                # Use log1p to handle zeros/negatives: log1p(x) = log(1 + x)
                # Shift by minimum value to ensure all positive
                log_shift <- abs(min_val) + 1
                stat_data <- stat_data %>%
                    mutate(!!stat_col := log1p(.data[[stat_col]] + log_shift))
                if (opts$verbose) cat("Applied log1p transformation with shift:", log_shift, "\n")
            } else {
                log_shift <- 0  # No shift needed
                stat_data <- stat_data %>%
                    mutate(!!stat_col := log(.data[[stat_col]]))
                if (opts$verbose) cat("Applied log transformation\n")
            }
        }
    }
    
    # Calculate statistics from ORIGINAL data (before transformation)
    # This ensures reference lines are calculated correctly
    if (opts$`panel-by` == "none" || opts$`panel-by` == "window") {
        # Genome-wide and per-chromosome, optionally by window
        if (opts$`panel-by` == "window" && "window_size" %in% colnames(stat_data_original) && !all(is.na(stat_data_original$window_size))) {
            stats_summary <- stat_data_original %>%
                group_by(.data$sample, .data$window_size, .data$chr) %>%
                summarise(
                    median = median(.data[[stat_col]], na.rm=TRUE),
                    mean = mean(.data[[stat_col]], na.rm=TRUE),
                    sd = sd(.data[[stat_col]], na.rm=TRUE),
                    q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                    q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                    q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                    q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    stat_data_original %>%
                        group_by(.data$sample, .data$window_size) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data[[stat_col]], na.rm=TRUE),
                            mean = mean(.data[[stat_col]], na.rm=TRUE),
                            sd = sd(.data[[stat_col]], na.rm=TRUE),
                            q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                            q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                            q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                            q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        } else {
            stats_summary <- stat_data_original %>%
                group_by(.data$sample, .data$chr) %>%
                summarise(
                    median = median(.data[[stat_col]], na.rm=TRUE),
                    mean = mean(.data[[stat_col]], na.rm=TRUE),
                    sd = sd(.data[[stat_col]], na.rm=TRUE),
                    q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                    q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                    q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                    q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    stat_data_original %>%
                        group_by(.data$sample) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data[[stat_col]], na.rm=TRUE),
                            mean = mean(.data[[stat_col]], na.rm=TRUE),
                            sd = sd(.data[[stat_col]], na.rm=TRUE),
                            q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                            q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                            q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                            q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
            if (!"window_size" %in% colnames(stats_summary)) {
                stats_summary$window_size <- NA
            }
        }
    } else {
        # Similar logic for other paneling options
        grouping_vars <- c("sample")
        if ("window_size" %in% colnames(stat_data_original) && !all(is.na(stat_data_original$window_size))) {
            # If paneling includes window dimension (both or sample with multiple windows), include it
            if (opts$`panel-by` %in% c("both", "sample")) {
            grouping_vars <- c("sample", "window_size")
            }
        }
        
        # Build group_by expression dynamically
        if (length(grouping_vars) == 1) {
            stats_summary <- stat_data_original %>%
                group_by(.data[[grouping_vars[1]]], .data$chr) %>%
                summarise(
                    median = median(.data[[stat_col]], na.rm=TRUE),
                    mean = mean(.data[[stat_col]], na.rm=TRUE),
                    sd = sd(.data[[stat_col]], na.rm=TRUE),
                    q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                    q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                    q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                    q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    stat_data_original %>%
                        group_by(.data[[grouping_vars[1]]]) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data[[stat_col]], na.rm=TRUE),
                            mean = mean(.data[[stat_col]], na.rm=TRUE),
                            sd = sd(.data[[stat_col]], na.rm=TRUE),
                            q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                            q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                            q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                            q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        } else {
            stats_summary <- stat_data_original %>%
                group_by(.data[[grouping_vars[1]]], .data[[grouping_vars[2]]], .data$chr) %>%
                summarise(
                    median = median(.data[[stat_col]], na.rm=TRUE),
                    mean = mean(.data[[stat_col]], na.rm=TRUE),
                    sd = sd(.data[[stat_col]], na.rm=TRUE),
                    q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                    q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                    q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                    q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    stat_data_original %>%
                        group_by(.data[[grouping_vars[1]]], .data[[grouping_vars[2]]]) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data[[stat_col]], na.rm=TRUE),
                            mean = mean(.data[[stat_col]], na.rm=TRUE),
                            sd = sd(.data[[stat_col]], na.rm=TRUE),
                            q5 = quantile(.data[[stat_col]], 0.05, na.rm=TRUE),
                            q1 = quantile(.data[[stat_col]], 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data[[stat_col]], 0.002, na.rm=TRUE),
                            q95 = quantile(.data[[stat_col]], 0.95, na.rm=TRUE),
                            q99 = quantile(.data[[stat_col]], 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data[[stat_col]], 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        }
        
        if (!"window_size" %in% colnames(stats_summary)) {
            stats_summary$window_size <- NA
        }
    }
    
    # Add statistic name
    stats_summary$statistic <- stat
    
    stats_results[[stat]] <- stats_summary
    
    # Prepare reference lines data for legend, panel-aware (sample/window_size if present)
    panel_cols <- intersect(c("sample", "window_size"), colnames(stats_summary))
    
    # Apply transformation to reference lines if needed
    if (opts$transform == "asinh") {
        # For asinh, reference lines need to be transformed using the global median and scale
        # The median itself becomes 0 after asinh((median - global_median) / scale_factor), so median = 0
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
    } else if (opts$transform == "log") {
        # For log transformation, apply the same transformation as to the data
        # Use original data to determine shift
        min_val <- min(stat_data_original[[stat_col]], na.rm=TRUE)
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
    
    if (!is.null(selected_chr)) {
        # Use stats from the selected chromosome
        chr_stats <- stats_summary %>% filter(.data$chr == selected_chr)
        if (nrow(chr_stats) > 0) {
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
        } else {
            ref_lines <- data.frame(
                yintercept = numeric(0),
                label = character(0),
                color = character(0),
                stringsAsFactors = FALSE
            )
        }
    } else {
        # Use genome-wide stats
        genome_stats <- stats_summary %>% filter(.data$chr == "genome")
        if (nrow(genome_stats) > 0) {
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
        } else {
            ref_lines <- data.frame(
                yintercept = numeric(0),
                label = character(0),
                color = character(0),
                stringsAsFactors = FALSE
            )
        }
    }
    
    # Create plot
    p <- ggplot(stat_data, aes(x=cum_pos, y=.data[[stat_col]])) +
        # Add chromosome stripes
        geom_rect(data=chr_stripes, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=fill),
                  inherit.aes=FALSE, alpha=0.3) +
        scale_fill_identity() +
        # Add data points
        geom_line(alpha=0.7, linewidth=0.5) +
        labs(
            x="Genomic Position (cumulative)",
            y=if (opts$transform == "asinh") {
                if (is.null(opts$`asinh-scale`)) {
                    switch(stat,
                        "pi" = expression(asinh((pi - median) / SD)),
                        "theta" = expression(asinh((theta - median) / SD)),
                        "tajima_d" = "asinh((Tajima's D - median) / SD)"
                    )
                } else {
                    switch(stat,
                        "pi" = paste0("asinh((pi - ", round(global_median, 4), ") / ", round(scale_factor, 4), ")"),
                        "theta" = paste0("asinh((theta - ", round(global_median, 4), ") / ", round(scale_factor, 4), ")"),
                        "tajima_d" = paste0("asinh((Tajima's D - ", round(global_median, 4), ") / ", round(scale_factor, 4), ")")
                    )
                }
            } else if (opts$transform == "log") {
                switch(stat,
                    "pi" = expression(log(pi)),
                    "theta" = expression(log(theta)),
                    "tajima_d" = "log(Tajima's D)"
                )
            } else {
                switch(stat,
                    "pi" = expression(pi ~ "(Nucleotide Diversity)"),
                    "theta" = expression(theta ~ "(Watterson's Theta)"),
                    "tajima_d" = "Tajima's D"
                )
            },
            title=paste("Diversity Statistic:", switch(stat,
                "pi" = expression(pi),
                "theta" = expression(theta),
                "tajima_d" = "Tajima's D"
            ))
        ) +
        theme_bw() +
        theme(
            panel.grid.minor=element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1)
        )
    
    # Add reference lines (genome-wide) with legend
    # Only use labels that are actually present in ref_lines (not all factor levels)
    if (nrow(ref_lines) > 0) {
        # Get unique labels that are actually present (not all factor levels)
        present_labels <- unique(ref_lines$label)
        present_labels <- present_labels[!is.na(present_labels)]
        
        # Create color and linetype mappings only for present labels
        color_map <- setNames(ref_lines$color, ref_lines$label)[as.character(present_labels)]
        linetype_map <- setNames(rep("dashed", length(present_labels)), as.character(present_labels))
        
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
    }
    
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
    
    # Add y-axis scale with inverse transformation for labels (if transformation is applied)
    # This labels the axis ticks with original (untransformed) values while plotting transformed data
    if (opts$transform == "asinh") {
        # Get range of original (untransformed) values
        orig_range <- range(stat_data_original[[stat_col]], na.rm=TRUE)
        
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
    } else if (opts$transform == "log") {
        # Get range of original (untransformed) values
        orig_range <- range(stat_data_original[[stat_col]], na.rm=TRUE)
        
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
    
    # Add paneling
    # Check if we have multiple window sizes and multiple samples
    has_window_size <- "window_size" %in% colnames(stat_data) && !all(is.na(stat_data$window_size))
    if (has_window_size) {
        unique_windows <- unique(stat_data$window_size[!is.na(stat_data$window_size)])
        multiple_windows <- length(unique_windows) > 1
    } else {
        multiple_windows <- FALSE
    }
    unique_samples <- unique(stat_data$sample)
    multiple_samples <- length(unique_samples) > 1
    
    if (opts$`panel-by` == "window" && has_window_size) {
        # Panel by window size vertically
        if (multiple_samples) {
            # If multiple samples, add horizontal paneling by sample
            p <- p + facet_grid(window_size ~ sample, scales="fixed")
        } else {
            # Single sample, just panel by window size
            p <- p + facet_wrap(~ window_size, scales="fixed", ncol=1)
        }
    } else if (opts$`panel-by` == "sample") {
        # Panel by sample vertically
        if (has_window_size && multiple_windows) {
            # If multiple window sizes, add horizontal paneling by window size
            p <- p + facet_grid(sample ~ window_size, scales="fixed")
        } else {
            # Single window size or no window size, just panel by sample
            p <- p + facet_wrap(~ sample, scales="fixed", ncol=1)
        }
    } else if (opts$`panel-by` == "both" && has_window_size) {
        # Explicitly panel by both: samples as rows, window sizes as columns
        p <- p + facet_grid(sample ~ window_size, scales="fixed")
    }
    # If panel-by == "none", no paneling is added
    
    # Add chromosome boundary labels (ordered by length, longest first)
    # Get chromosome labels in the correct order (longest first) - only those in data
    chr_labels <- chr_lengths$chr[chr_lengths$chr %in% unique(stat_data$chr)]
    
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
    stat_label <- switch(stat,
        "pi" = "pi",
        "theta" = "theta",
        "tajima_d" = "tajima_d"
    )
    
    # Build filename suffixes with actual values separated by dashes
    # Window size suffix - always include
    if (!is.null(selected_window_size)) {
        window_suffix <- paste0("_w", selected_window_size)
    } else if (has_window_size) {
        # List all window sizes separated by dashes
        unique_windows <- sort(unique(stat_data$window_size[!is.na(stat_data$window_size)]))
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
        unique_chrs <- unique(stat_data$chr[!is.na(stat_data$chr)])
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
    
    # Sample suffix - always include
    if (multiple_samples) {
        # List all sample names separated by dashes
        unique_samples_list <- sort(unique(stat_data$sample[!is.na(stat_data$sample)]))
        # If more than 3, just show first one and count
        if (length(unique_samples_list) > 3) {
            sample_list <- paste0(gsub("[^A-Za-z0-9]", "_", unique_samples_list[1]), "-and", length(unique_samples_list)-1, "more")
            sample_suffix <- paste0("_samples", sample_list)
        } else {
            sample_list <- paste(gsub("[^A-Za-z0-9]", "_", unique_samples_list), collapse="-")
            sample_suffix <- paste0("_samples", sample_list)
        }
    } else if (length(unique_samples) == 1) {
        # Single sample - include its name
        sample_safe <- gsub("[^A-Za-z0-9]", "_", unique_samples[1])
        # Truncate single sample name if too long
        if (nchar(sample_safe) > 80) {
            sample_safe <- substr(sample_safe, 1, 80)
        }
        sample_suffix <- paste0("_sample", sample_safe)
    } else {
        sample_suffix <- ""
    }
    
    # Transform suffix - add if transform is not "none"
    transform_suffix <- ""
    if (opts$transform != "none") {
        if (opts$transform == "asinh") {
            transform_suffix <- "_asinhtrans"
        } else if (opts$transform == "log") {
            transform_suffix <- "_logtrans"
        }
    }
    
    if ("png" %in% format_parts) {
        png_file <- file.path(opts$`output-dir`, paste0(prefix, stat_label, window_suffix, chr_suffix, sample_suffix, transform_suffix, ".png"))
        ggsave(png_file, p, width=opts$width, height=opts$height, dpi=300)
        cat("Saved:", png_file, "\n")
    }
    
    if ("pdf" %in% format_parts) {
        pdf_file <- file.path(opts$`output-dir`, paste0(prefix, stat_label, window_suffix, chr_suffix, sample_suffix, transform_suffix, ".pdf"))
        ggsave(pdf_file, p, width=opts$width, height=opts$height)
        cat("Saved:", pdf_file, "\n")
    }
    
    if ("svg" %in% format_parts) {
        svg_file <- file.path(opts$`output-dir`, paste0(prefix, stat_label, window_suffix, chr_suffix, sample_suffix, transform_suffix, ".svg"))
        ggsave(svg_file, p, width=opts$width, height=opts$height, device="svg")
        cat("Saved:", svg_file, "\n")
    }
}

# Combine and save statistics
if (length(stats_results) > 0) {
    all_stats <- bind_rows(stats_results) %>%
        select(statistic, sample, window_size, chr, median, mean, sd,
               q5, q1, q0.2, q95, q99, q99.8) %>%
        arrange(statistic, sample, window_size, chr)
    
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
    
    unique_samples_stats <- if ("sample" %in% colnames(all_stats)) {
        sort(unique(all_stats$sample[!is.na(all_stats$sample)]))
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
    
    # Sample suffix
    if (length(unique_samples_stats) > 1) {
        # If more than 3, just show first one and count
        if (length(unique_samples_stats) > 3) {
            sample_list <- paste0(gsub("[^A-Za-z0-9]", "_", unique_samples_stats[1]), "-and", length(unique_samples_stats)-1, "more")
            sample_suffix <- paste0("_samples", sample_list)
        } else {
            sample_list <- paste(gsub("[^A-Za-z0-9]", "_", unique_samples_stats), collapse="-")
            sample_suffix <- paste0("_samples", sample_list)
        }
    } else if (length(unique_samples_stats) == 1) {
        sample_safe <- gsub("[^A-Za-z0-9]", "_", unique_samples_stats[1])
        # Truncate single sample name if too long
        if (nchar(sample_safe) > 50) {
            sample_safe <- substr(sample_safe, 1, 50)
        }
        sample_suffix <- paste0("_sample", sample_safe)
    } else {
        sample_suffix <- ""
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
    
    stats_file <- file.path(opts$`output-dir`, paste0(prefix, "diversity_statistics", window_suffix, chr_suffix, sample_suffix, transform_suffix_stats, ".csv"))
    write_csv(all_stats, stats_file)
    cat("Saved statistics to:", stats_file, "\n")
} else {
    if (opts$verbose) cat("\nWarning: No statistics were calculated\n")
}

if (opts$verbose) cat("\nDone!\n")
