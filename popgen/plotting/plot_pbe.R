#!/usr/bin/env Rscript

###############################################################################
# plot_pbe.R
# 
# Create ggplot2 plots of PBE statistics from calculate_pbe.sh output files.
#
# Features:
# - Chromosomal positions on x-axis with alternating white/grey stripes
# - Paneling by window size, sample trio, or both
# - Horizontal reference lines (median, mean, 95th percentile)
# - Statistics calculated genome-wide and per-chromosome
# - Handles symmetric ordering: pop1:pop2:pop3 is equivalent to pop1:pop3:pop2
#
# Author: Based on plot_fst.R
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
                help="Directory containing PBE TSV files", metavar="DIR"),
    make_option(c("--output-dir"), type="character", default=".",
                help="Output directory for plots [default: %default]", metavar="DIR"),
    make_option(c("--reference-genome"), type="character", default=NULL,
                help="Reference genome FASTA file (optional, for chromosome lengths)", metavar="FILE"),
    make_option(c("--panel-by"), type="character", default="both",
                help="Paneling option: window, trio, both, or none [default: %default]",
                metavar="OPTION"),
    make_option(c("--plot-format"), type="character", default="png",
                help="Plot format: png, pdf, svg, or any combination (e.g., both, all) [default: %default]", metavar="FORMAT"),
    make_option(c("--width"), type="numeric", default=12,
                help="Plot width in inches [default: %default]", metavar="NUMBER"),
    make_option(c("--height"), type="numeric", default=8,
                help="Plot height in inches [default: %default]", metavar="NUMBER"),
    make_option(c("--file-prefix"), type="character", default="",
                help="Prefix for output files [default: no prefix]", metavar="PREFIX"),
    make_option(c("--sample-trios"), type="character", default="",
                help="Comma-separated list of sample trios to plot (format: pop1:pop2:pop3, default: all trios)", metavar="LIST"),
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

if (!opts$`panel-by` %in% c("window", "trio", "both", "none")) {
    stop("--panel-by must be one of: window, trio, both, none")
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

# Find all PBE files (TSV or CSV)
input_dir <- opts$`input-dir`
tsv_files <- list.files(input_dir, pattern=".*pbe.*\\.tsv$", full.names=TRUE, recursive=TRUE)
csv_files <- list.files(input_dir, pattern=".*pbe.*\\.csv$", full.names=TRUE, recursive=TRUE)
all_files <- c(tsv_files, csv_files)

if (length(all_files) == 0) {
    stop("No PBE files found in: ", input_dir, " (expected files matching pattern: *pbe*.tsv or *pbe*.csv)")
}

if (opts$verbose) {
    cat("Found", length(all_files), "PBE file(s) (", length(tsv_files), " TSV, ", length(csv_files), " CSV)\n", sep="")
}

# Filter files by window size if specified (before reading)
selected_window_size <- opts$`window-size`
if (!is.null(selected_window_size)) {
    file_window_sizes <- sapply(all_files, function(f) extract_window_size_from_filename(basename(f)))
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

# Helper function to extract trio names from filename
extract_trio_from_filename <- function(filename) {
    # Try pattern: w1000_s500_Echo_Kjer_Cheney_pbe.tsv
    # Extract the part between window/step and _pbe
    base <- basename(filename)
    # Remove extension
    base <- sub("\\.(tsv|csv)$", "", base, ignore.case=TRUE)
    # Remove _pbe suffix if present
    base <- sub("_pbe$", "", base, ignore.case=TRUE)
    # Try to extract pop names after window pattern
    # Pattern: w\d+_s\d+_(.+)
    trio_match <- regmatches(base, regexpr("w\\d+_s\\d+_(.+)", base))
    if (length(trio_match) > 0) {
        trio_part <- sub("w\\d+_s\\d+_", "", trio_match)
        parts <- strsplit(trio_part, "_")[[1]]
        if (length(parts) >= 3) {
            return(paste(parts[1:3], collapse=":"))
        }
    }
    # Try pattern without window: Echo_Kjer_Cheney_pbe.tsv
    trio_match <- regmatches(base, regexpr("^([A-Za-z0-9_]+)_([A-Za-z0-9_]+)_([A-Za-z0-9_]+)$", base))
    if (length(trio_match) > 0) {
        parts <- strsplit(trio_match, "_")[[1]]
        if (length(parts) >= 3) {
            return(paste(parts[1:3], collapse=":"))
        }
    }
    return(NULL)
}

# Helper function to normalize trio name (handle symmetric ordering)
normalize_trio_name <- function(trio_str) {
    parts <- strsplit(trio_str, ":")[[1]]
    if (length(parts) != 3) return(trio_str)
    pop1 <- normalize_sample_name(parts[1])
    pop2 <- normalize_sample_name(parts[2])
    pop3 <- normalize_sample_name(parts[3])
    # PBE is symmetric with respect to pop2 and pop3, so normalize order
    if (pop2 > pop3) {
        return(paste(pop1, pop3, pop2, sep=":"))
    }
    return(paste(pop1, pop2, pop3, sep=":"))
}

# Filter files by sample trios if specified (peek at headers and filenames)
selected_trios <- NULL
if (opts$`sample-trios` != "") {
    selected_trios <- strsplit(opts$`sample-trios`, ",")[[1]]
    selected_trios <- trimws(selected_trios)
    # Normalize selected trios (handle symmetric ordering)
    selected_trios_normalized <- sapply(selected_trios, function(t) {
        parts <- strsplit(t, ":")[[1]]
        if (length(parts) == 3) {
            return(normalize_trio_name(t))
        }
        return(t)
    })
    
    if (opts$verbose) {
        cat("Filtering files by sample trios:", paste(selected_trios_normalized, collapse=", "), "\n")
    }
    
    # Peek at each file's column names and filename to see if it contains the requested trios
    files_to_keep <- character(0)
    for (file in all_files) {
        file_has_trio <- FALSE
        
        # First, try to extract trio from filename
        filename_trio <- extract_trio_from_filename(file)
        if (!is.null(filename_trio)) {
            filename_trio_normalized <- normalize_trio_name(filename_trio)
            # Check if this matches any selected trio (or its symmetric variant)
            for (sel_trio in selected_trios_normalized) {
                sel_parts <- strsplit(sel_trio, ":")[[1]]
                filename_parts <- strsplit(filename_trio_normalized, ":")[[1]]
                if (length(sel_parts) == 3 && length(filename_parts) == 3) {
                    # Check if pop1 matches and pop2/pop3 match (order-independent)
                    if (sel_parts[1] == filename_parts[1] &&
                        ((sel_parts[2] == filename_parts[2] && sel_parts[3] == filename_parts[3]) ||
                         (sel_parts[2] == filename_parts[3] && sel_parts[3] == filename_parts[2]))) {
                        file_has_trio <- TRUE
                        break
                    }
                }
            }
        }
        
        # Also check column names for PBE column
        if (!file_has_trio) {
            col_names <- peek_column_names(file)
            pbe_cols <- col_names[grepl("^pbe$", col_names, ignore.case=TRUE)]
            if (length(pbe_cols) > 0) {
                # File has PBE column, include it (we'll filter later by trio if needed)
                file_has_trio <- TRUE
            }
        }
        
        if (file_has_trio) {
            files_to_keep <- c(files_to_keep, file)
        }
    }
    
    all_files <- files_to_keep
    if (opts$verbose) {
        cat("Filtered to", length(all_files), "file(s) containing requested sample trios\n")
    }
    if (length(all_files) == 0) {
        stop("No files found containing the specified sample trios")
    }
}

# Note: selected_pairs parsing moved earlier for file filtering

# Read and combine all PBE files using tidyverse
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
    
    # Convert PBE column to numeric (look for "pbe" column, case-insensitive)
    pbe_stat_cols <- grep("^pbe$", colnames(data), ignore.case=TRUE, value=TRUE)
    if (length(pbe_stat_cols) > 0) {
        if (opts$verbose) {
            cat("  Found PBE column:", paste(pbe_stat_cols, collapse=", "), "\n")
        }
        data <- data %>%
            mutate(across(all_of(pbe_stat_cols), as.numeric))
    } else {
        if (opts$verbose) cat("  Warning: No PBE column found (expected column named 'PBE' or 'pbe')\n")
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
    
    # Extract trio from filename and add as metadata
    filename_trio <- extract_trio_from_filename(filename)
    if (!is.null(filename_trio)) {
        filename_trio_normalized <- normalize_trio_name(filename_trio)
        data <- data %>% mutate(source_trio = filename_trio_normalized, source_file = filename)
    } else {
        data <- data %>% mutate(source_trio = NA_character_, source_file = filename)
    }
    
    if (opts$verbose) cat("  Rows:", nrow(data), "\n")
    if (opts$verbose && !is.null(filename_trio)) {
        cat("  Extracted trio from filename:", filename_trio_normalized, "\n")
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

# Parse PBE column
pbe_cols <- colnames(combined_data)[grepl("^pbe$", colnames(combined_data), ignore.case=TRUE)]
if (length(pbe_cols) == 0) {
    stop("No PBE column found in input files. Expected column named 'PBE' or 'pbe'")
}

if (opts$verbose) {
    cat("Found", length(pbe_cols), "PBE column(s)\n")
    cat("PBE columns:", paste(pbe_cols, collapse=", "), "\n")
}

# Parse sample trios from combined data
# Use source_trio column added during file reading, or infer from PBE column
pbe_trios_raw <- list()

# Get unique trios from source_trio column (if present)
if ("source_trio" %in% colnames(combined_data)) {
    unique_source_trios <- unique(combined_data$source_trio[!is.na(combined_data$source_trio)])
    for (trio_name in unique_source_trios) {
        if (!trio_name %in% names(pbe_trios_raw)) {
            pbe_trios_raw[[trio_name]] <- character(0)
        }
        pbe_trios_raw[[trio_name]] <- c(pbe_trios_raw[[trio_name]], "pbe")
    }
}

# If no source_trio found, use generic approach
if (length(pbe_trios_raw) == 0) {
    # Use a single generic trio name
    pbe_trios_raw[["trio_all"]] <- "pbe"
}

# Second pass: for each normalized trio, select the best column (one with most non-NA data)
pbe_trios <- list()
for (trio_name in names(pbe_trios_raw)) {
    cols_for_trio <- pbe_trios_raw[[trio_name]]
    if (length(cols_for_trio) == 1) {
        # Only one column, use it
        pbe_trios[[trio_name]] <- cols_for_trio[1]
    } else {
        # Multiple columns map to same normalized trio
        # Choose the one with the most non-NA data
        if (opts$verbose) {
            cat("  Multiple columns for trio '", trio_name, "': ", paste(cols_for_trio, collapse=", "), "\n", sep="")
        }
        best_col <- NULL
        best_count <- -1
        for (col in cols_for_trio) {
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
            pbe_trios[[trio_name]] <- best_col
            if (opts$verbose && length(cols_for_trio) > 1) {
                cat("    Selected: ", best_col, " (", best_count, " non-NA values)\n", sep="")
            }
        } else {
            # Fallback: use first column
            pbe_trios[[trio_name]] <- cols_for_trio[1]
            if (opts$verbose) {
                cat("    Warning: No valid column found, using first: ", cols_for_trio[1], "\n", sep="")
            }
        }
    }
}

if (length(pbe_trios) == 0) {
    stop("Could not parse sample trios from PBE files")
}

if (opts$verbose) {
    cat("\nNormalized sample trios:\n")
    for (trio_name in names(pbe_trios)) {
        cat("  ", trio_name, " -> ", pbe_trios[[trio_name]], "\n", sep="")
    }
}

# Filter trios if specified (use normalized trio names, handle symmetric ordering)
if (!is.null(selected_trios)) {
    # selected_trios_normalized was created during file filtering
    # If it doesn't exist (no file filtering happened), create it now
    if (!exists("selected_trios_normalized")) {
        selected_trios_normalized <- sapply(selected_trios, function(t) {
            parts <- strsplit(t, ":")[[1]]
            if (length(parts) == 3) {
                return(normalize_trio_name(t))
            }
            return(t)
        })
    }
    # Also check symmetric variants (pop1:pop2:pop3 vs pop1:pop3:pop2)
    selected_trios_expanded <- selected_trios_normalized
    for (trio in selected_trios_normalized) {
        parts <- strsplit(trio, ":")[[1]]
        if (length(parts) == 3) {
            # Add symmetric variant
            symmetric <- paste(parts[1], parts[3], parts[2], sep=":")
            selected_trios_expanded <- c(selected_trios_expanded, symmetric)
        }
    }
    selected_trios_expanded <- unique(selected_trios_expanded)
    
    # Filter trios that match selected trios (exact or symmetric)
    pbe_trios_filtered <- list()
    for (trio_name in names(pbe_trios)) {
        trio_parts <- strsplit(trio_name, ":")[[1]]
        if (length(trio_parts) == 3) {
            # Check if this trio matches any selected trio (exact or symmetric)
            for (sel_trio in selected_trios_expanded) {
                sel_parts <- strsplit(sel_trio, ":")[[1]]
                if (length(sel_parts) == 3) {
                    if (trio_parts[1] == sel_parts[1] &&
                        ((trio_parts[2] == sel_parts[2] && trio_parts[3] == sel_parts[3]) ||
                         (trio_parts[2] == sel_parts[3] && trio_parts[3] == sel_parts[2]))) {
                        pbe_trios_filtered[[trio_name]] <- pbe_trios[[trio_name]]
                        break
                    }
                }
            }
        } else {
            # Generic trio name, include if any selected trio matches pattern
            if (trio_name %in% selected_trios_expanded) {
                pbe_trios_filtered[[trio_name]] <- pbe_trios[[trio_name]]
            }
        }
    }
    pbe_trios <- pbe_trios_filtered
    if (length(pbe_trios) == 0) {
        stop("None of the specified sample trios were found in the data")
    }
}

if (opts$verbose) {
    cat("Processing", length(pbe_trios), "sample trio(s)\n")
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

# Process each sample trio
stats_results <- list()

for (trio_name in names(pbe_trios)) {
    pbe_col <- pbe_trios[[trio_name]]
    if (opts$verbose) {
        cat("\nProcessing sample trio:", trio_name, "\n")
    }
    
    # Extract PBE values for this trio using tidyverse
    # Filter by source_trio if available, otherwise use all data
    if ("source_trio" %in% colnames(combined_data)) {
        trio_data_raw <- combined_data %>%
            filter(is.na(.data$source_trio) | .data$source_trio == trio_name)
    } else {
        trio_data_raw <- combined_data
    }
    
    # Check if PBE column exists
    if (!"pbe" %in% colnames(trio_data_raw)) {
        cat("Warning: PBE column not found in combined_data for trio '", trio_name, "', skipping\n", sep="")
        if (opts$verbose) {
            cat("  Available columns:", paste(head(colnames(combined_data), 20), collapse=", "), if(length(colnames(combined_data)) > 20) "..." else "", "\n")
            cat("  Total rows in combined_data:", nrow(combined_data), "\n")
        }
        next
    }
    
    # Check how many non-NA values exist
    if (opts$verbose) {
        non_na_count <- sum(!is.na(trio_data_raw$pbe))
        cat("  PBE column has ", non_na_count, " non-NA values for trio '", trio_name, "'\n", sep="")
        if (non_na_count == 0) {
            cat("  ERROR: PBE column exists but has no non-NA values!\n")
            next
        }
    }
    
    # Select columns that exist (cum_pos might not exist if add_cumulative_positions failed)
    cols_to_select <- c("chr", "pos", "window_size", "pbe")
    if ("cum_pos" %in% colnames(trio_data_raw)) {
        cols_to_select <- c("chr", "pos", "cum_pos", "window_size", "pbe")
    }
    
    trio_data <- trio_data_raw %>%
        select(all_of(cols_to_select)) %>%
        filter(!is.na(.data$pbe))
    
    # Add cum_pos if it's missing (shouldn't happen, but safety check)
    if (!"cum_pos" %in% colnames(trio_data)) {
        if (opts$verbose) {
            cat("  Warning: cum_pos missing, recalculating...\n")
        }
        trio_data <- add_cumulative_positions(trio_data, chr_lengths)
    }
    
    # Sort by cum_pos to ensure line connects properly (critical when combining multiple files)
    trio_data <- trio_data %>%
        arrange(.data$cum_pos) %>%
        mutate(sample_trio = trio_name)
    
    if (opts$verbose) {
        cat("  Data sorted by cum_pos for proper line connection\n")
        # Check for any duplicate positions that might cause issues
        dup_positions <- sum(duplicated(trio_data$cum_pos))
        if (dup_positions > 0) {
            cat("  Warning: Found", dup_positions, "duplicate cum_pos values (may cause line artifacts)\n")
        }
    }
    
    if (nrow(trio_data) == 0) {
        cat("Warning: No valid PBE data for trio '", trio_name, "', skipping\n", sep="")
        if (opts$verbose) {
            cat("  Total rows in combined_data:", nrow(combined_data), "\n")
            cat("  Non-NA values in column '", pbe_col, "': ", sum(!is.na(combined_data[[pbe_col]])), "\n", sep="")
        }
        next
    }
    
    if (opts$verbose) {
        cat("  Extracted", nrow(trio_data), "rows with valid PBE data\n")
        cat("  PBE range:", min(trio_data$pbe, na.rm=TRUE), "to", max(trio_data$pbe, na.rm=TRUE), "\n")
    }
    
    # Store original data for statistics calculation
    trio_data_original <- trio_data
    
    # Calculate global median and standard deviation for asinh transformation (if needed)
    # Also store log shift if needed for axis label inverse transformation
    global_median <- NULL
    scale_factor <- NULL
    log_shift <- NULL
    if (opts$transform == "asinh") {
        global_median <- median(trio_data_original$pbe, na.rm=TRUE)
        global_sd <- sd(trio_data_original$pbe, na.rm=TRUE)
        
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
            trio_data <- trio_data %>%
                mutate(pbe = asinh((.data$pbe - global_median) / scale_factor))
        } else if (opts$transform == "log") {
            # Check for zeros or negatives
            min_val <- min(trio_data_original$pbe, na.rm=TRUE)
            if (min_val <= 0) {
                # Use log1p to handle zeros/negatives: log1p(x) = log(1 + x)
                # Shift by minimum value to ensure all positive
                log_shift <- abs(min_val) + 1
                trio_data <- trio_data %>%
                    mutate(pbe = log1p(.data$pbe + log_shift))
                if (opts$verbose) cat("Applied log1p transformation with shift:", log_shift, "\n")
            } else {
                log_shift <- 0  # No shift needed
                trio_data <- trio_data %>%
                    mutate(pbe = log(.data$pbe))
                if (opts$verbose) cat("Applied log transformation\n")
            }
        }
    }
    
    # Calculate statistics from ORIGINAL data (before transformation)
    # This ensures reference lines are calculated correctly
    if (opts$`panel-by` == "none" || opts$`panel-by` == "window") {
        # Genome-wide and per-chromosome, optionally by window
        if (opts$`panel-by` == "window" && "window_size" %in% colnames(trio_data_original) && !all(is.na(trio_data_original$window_size))) {
            stats_summary <- trio_data_original %>%
                group_by(.data$sample_trio, .data$window_size, .data$chr) %>%
                summarise(
                    median = median(.data$pbe, na.rm=TRUE),
                    mean = mean(.data$pbe, na.rm=TRUE),
                    sd = sd(.data$pbe, na.rm=TRUE),
                    q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    trio_data_original %>%
                        group_by(.data$sample_trio, .data$window_size) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$pbe, na.rm=TRUE),
                            mean = mean(.data$pbe, na.rm=TRUE),
                            sd = sd(.data$pbe, na.rm=TRUE),
                            q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        } else {
            stats_summary <- trio_data_original %>%
                group_by(.data$sample_trio, .data$chr) %>%
                summarise(
                    median = median(.data$pbe, na.rm=TRUE),
                    mean = mean(.data$pbe, na.rm=TRUE),
                    sd = sd(.data$pbe, na.rm=TRUE),
                    q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    trio_data_original %>%
                        group_by(.data$sample_trio) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$pbe, na.rm=TRUE),
                            mean = mean(.data$pbe, na.rm=TRUE),
                            sd = sd(.data$pbe, na.rm=TRUE),
                            q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
            if (!"window_size" %in% colnames(stats_summary)) {
                stats_summary$window_size <- NA
            }
        }
    } else {
        # Similar logic for other paneling options
        grouping_vars <- c("sample_trio")
        if ("window_size" %in% colnames(trio_data_original) && !all(is.na(trio_data_original$window_size))) {
            # If paneling includes window dimension (both or trio with multiple windows), include it
            if (opts$`panel-by` %in% c("both", "trio")) {
            grouping_vars <- c("sample_trio", "window_size")
            }
        }
        
        # Build group_by expression dynamically
        if (length(grouping_vars) == 1) {
            stats_summary <- trio_data_original %>%
                group_by(.data[[grouping_vars[1]]], .data$chr) %>%
                summarise(
                    median = median(.data$pbe, na.rm=TRUE),
                    mean = mean(.data$pbe, na.rm=TRUE),
                    sd = sd(.data$pbe, na.rm=TRUE),
                    q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    trio_data_original %>%
                        group_by(.data[[grouping_vars[1]]]) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$pbe, na.rm=TRUE),
                            mean = mean(.data$pbe, na.rm=TRUE),
                            sd = sd(.data$pbe, na.rm=TRUE),
                            q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        } else {
            stats_summary <- trio_data_original %>%
                group_by(.data[[grouping_vars[1]]], .data[[grouping_vars[2]]], .data$chr) %>%
                summarise(
                    median = median(.data$pbe, na.rm=TRUE),
                    mean = mean(.data$pbe, na.rm=TRUE),
                    sd = sd(.data$pbe, na.rm=TRUE),
                    q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                    q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                    q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                    q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                    q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                    q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                    .groups="drop"
                ) %>%
                bind_rows(
                    trio_data_original %>%
                        group_by(.data[[grouping_vars[1]]], .data[[grouping_vars[2]]]) %>%
                        summarise(
                            chr = "genome",
                            median = median(.data$pbe, na.rm=TRUE),
                            mean = mean(.data$pbe, na.rm=TRUE),
                            sd = sd(.data$pbe, na.rm=TRUE),
                            q5 = quantile(.data$pbe, 0.05, na.rm=TRUE),
                            q1 = quantile(.data$pbe, 0.01, na.rm=TRUE),
                            q0.2 = quantile(.data$pbe, 0.002, na.rm=TRUE),
                            q95 = quantile(.data$pbe, 0.95, na.rm=TRUE),
                            q99 = quantile(.data$pbe, 0.99, na.rm=TRUE),
                            q99.8 = quantile(.data$pbe, 0.998, na.rm=TRUE),
                            .groups="drop"
                        )
                )
        }
        
        if (!"window_size" %in% colnames(stats_summary)) {
            stats_summary$window_size <- NA
        }
    }
    
    stats_results[[trio_name]] <- stats_summary
    
    # Prepare reference lines data for legend, panel-aware (trio/window_size if present)
    # First, determine which panel columns will actually be used for faceting
    # This must match the faceting logic below
    has_window_size <- "window_size" %in% colnames(trio_data) && !all(is.na(trio_data$window_size))
    if (has_window_size) {
        unique_windows <- unique(trio_data$window_size[!is.na(trio_data$window_size)])
        multiple_windows <- length(unique_windows) > 1
    } else {
        multiple_windows <- FALSE
    }
    unique_trios <- unique(trio_data$sample_trio)
    multiple_trios <- length(unique_trios) > 1
    
    # Determine which panel columns will be used for faceting based on panel-by option
    facet_panel_cols <- character(0)
    if (opts$`panel-by` == "window" && has_window_size) {
        facet_panel_cols <- c("window_size")
        if (multiple_trios) {
            facet_panel_cols <- c(facet_panel_cols, "sample_trio")
        }
    } else if (opts$`panel-by` == "trio" && multiple_trios) {
        facet_panel_cols <- c("sample_trio")
        if (has_window_size && multiple_windows) {
            facet_panel_cols <- c(facet_panel_cols, "window_size")
        }
    } else if (opts$`panel-by` == "both" && has_window_size) {
        facet_panel_cols <- c("sample_trio", "window_size")
    }
    # If panel-by == "none", no panel columns
    
    # Only include panel columns that will actually be used for faceting
    panel_cols <- intersect(facet_panel_cols, colnames(stats_summary))
    
    # Apply transformation to reference lines if needed
    if (opts$transform == "asinh") {
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
    } else if (opts$transform == "log") {
        # For log transformation, apply the same transformation as to the data
        # Use original data to determine shift
        min_val <- min(trio_data_original$pbe, na.rm=TRUE)
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
        cat("    multiple_trios:", multiple_trios, "\n")
        cat("    panel-by option:", opts$`panel-by`, "\n")
        cat("    facet_panel_cols:", paste(facet_panel_cols, collapse=", "), "\n")
        cat("    panel_cols (after intersect):", paste(panel_cols, collapse=", "), "\n")
    }
    
    # Get unique panel combinations from trio_data to ensure ref_lines matches
    if (length(panel_cols) > 0) {
        trio_panels <- trio_data %>%
            select(all_of(panel_cols)) %>%
            distinct()
        if (opts$verbose) {
            cat("    Unique panel combinations in trio_data:", nrow(trio_panels), "\n")
            if (nrow(trio_panels) > 0) {
                print(trio_panels)
            }
        }
    } else {
        trio_panels <- data.frame()
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
            # Ensure we have stats for all panel combinations in trio_data
            if (nrow(trio_panels) > 0 && length(panel_cols) > 0) {
                if (opts$verbose) {
                    cat("    Merging chr_stats with trio_panels by:", paste(panel_cols, collapse=", "), "\n")
                    cat("    chr_stats before merge:", nrow(chr_stats), "rows\n")
                }
                # Merge with trio_panels to ensure all combinations are present
                chr_stats <- trio_panels %>%
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
            # Ensure we have stats for all panel combinations in trio_data
            if (nrow(trio_panels) > 0 && length(panel_cols) > 0) {
                if (opts$verbose) {
                    cat("    Merging genome_stats with trio_panels by:", paste(panel_cols, collapse=", "), "\n")
                    cat("    genome_stats before merge:", nrow(genome_stats), "rows\n")
                }
                # Merge with trio_panels to ensure all combinations are present
                genome_stats <- trio_panels %>%
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
        cat("  Preparing plot with", nrow(trio_data), "data points\n")
        cat("  cum_pos range:", if("cum_pos" %in% colnames(trio_data)) paste(range(trio_data$cum_pos, na.rm=TRUE), collapse=" to ") else "MISSING", "\n")
        cat("  pbe range:", paste(range(trio_data$pbe, na.rm=TRUE), collapse=" to "), "\n")
        cat("  Sample of first 5 rows:\n")
        print(head(trio_data %>% select(chr, pos, cum_pos, pbe, sample_trio), 5))
    }
    
    if (!"cum_pos" %in% colnames(trio_data)) {
        cat("Error: cum_pos column missing for trio '", trio_name, "'. Cannot create plot.\n", sep="")
        next
    }
    
    if (nrow(trio_data) == 0) {
        cat("Error: No data points to plot for trio '", trio_name, "'\n", sep="")
        next
    }
    
    # Check if all PBE values are the same (which would make the line invisible)
    pbe_unique <- length(unique(trio_data$pbe[!is.na(trio_data$pbe)]))
    if (pbe_unique == 1) {
        cat("Warning: All PBE values are identical (", unique(trio_data$pbe[!is.na(trio_data$pbe)])[1], ") for trio '", trio_name, "'\n", sep="")
    }
    
    # Create plot
    # Ensure data is sorted by cum_pos (critical for line plotting)
    trio_data <- trio_data %>% arrange(.data$cum_pos)
    
    # Create plot with explicit grouping to ensure single continuous line
    # When combining data from multiple files, we need to ensure the line connects properly
    p <- ggplot(trio_data, aes(x=cum_pos, y=pbe)) +
        # Add chromosome stripes
        geom_rect(data=chr_stripes, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill=fill),
                  inherit.aes=FALSE, alpha=0.3) +
        scale_fill_identity() +
        # Add data points - use group=1 to ensure single continuous line across all data
        # This is especially important when combining data from multiple files
        geom_line(alpha=0.7, linewidth=0.5, group=1, na.rm=TRUE) +
        labs(
            x="Genomic Position (cumulative)",
            y=if (opts$transform == "asinh") {
                if (is.null(opts$`asinh-scale`)) {
                    "asinh((PBE - median) / SD)"
                } else {
                    paste0("asinh((PBE - ", round(global_median, 4), ") / ", round(scale_factor, 4), ")")
                }
            } else if (opts$transform == "log") {
                "log(PBE)"
            } else {
                "PBE (Population Branch Excess)"
            },
            title=paste("PBE:", trio_name)
        ) +
        theme_bw() +
        theme(
            panel.grid.minor=element_blank(),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1)
        )
    
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
    
    if (opts$transform == "asinh") {
        # Get range of original (untransformed) values
        orig_range <- range(trio_data_original$pbe, na.rm=TRUE)
        
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
        orig_range <- range(trio_data_original$pbe, na.rm=TRUE)
        
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
    # Check if we have multiple window sizes and multiple trios
    has_window_size <- "window_size" %in% colnames(trio_data) && !all(is.na(trio_data$window_size))
    if (has_window_size) {
        unique_windows <- unique(trio_data$window_size[!is.na(trio_data$window_size)])
        multiple_windows <- length(unique_windows) > 1
    } else {
        multiple_windows <- FALSE
    }
    unique_trios <- unique(trio_data$sample_trio)
    multiple_trios <- length(unique_trios) > 1
    
    if (opts$`panel-by` == "window" && has_window_size) {
        # Panel by window size vertically
        if (multiple_trios) {
            # If multiple trios, add horizontal paneling by trio
            p <- p + facet_grid(window_size ~ sample_trio, scales="fixed")
        } else {
            # Single trio, just panel by window size
            p <- p + facet_wrap(~ window_size, scales="fixed", ncol=1)
        }
    } else if (opts$`panel-by` == "trio" && multiple_trios) {
        # Panel by trio vertically
        if (has_window_size && multiple_windows) {
            # If multiple window sizes, add horizontal paneling by window size
            p <- p + facet_grid(sample_trio ~ window_size, scales="fixed")
        } else {
            # Single window size or no window size, just panel by trio
            p <- p + facet_wrap(~ sample_trio, scales="fixed", ncol=1)
        }
    } else if (opts$`panel-by` == "both" && has_window_size) {
        # Explicitly panel by both: trios as rows, window sizes as columns
        p <- p + facet_grid(sample_trio ~ window_size, scales="fixed")
    }
    # If panel-by == "none", no paneling is added
    
    # Add chromosome boundary labels (ordered by length, longest first)
    # Get chromosome labels in the correct order (longest first) - only those in data
    chr_labels <- chr_lengths$chr[chr_lengths$chr %in% unique(trio_data$chr)]
    
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
    trio_safe <- gsub(":", "_", trio_name)
    trio_safe <- gsub("[^A-Za-z0-9_]", "_", trio_safe)
    # Truncate trio name if too long to avoid filename length issues
    if (nchar(trio_safe) > 80) {
        trio_safe <- substr(trio_safe, 1, 80)
    }
    
    # Build filename suffixes with actual values separated by dashes
    # Window size suffix - always include
    if (!is.null(selected_window_size)) {
        window_suffix <- paste0("_w", selected_window_size)
    } else if (has_window_size) {
        # List all window sizes separated by dashes
        unique_windows <- sort(unique(trio_data$window_size[!is.na(trio_data$window_size)]))
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
        unique_chrs <- unique(trio_data$chr[!is.na(trio_data$chr)])
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
        png_file <- file.path(opts$`output-dir`, paste0(prefix, "pbe_", trio_safe, window_suffix, chr_suffix, transform_suffix, ".png"))
        ggsave(png_file, p, width=opts$width, height=opts$height, dpi=300)
        cat("Saved:", png_file, "\n")
    }
    
    if ("pdf" %in% format_parts) {
        pdf_file <- file.path(opts$`output-dir`, paste0(prefix, "pbe_", trio_safe, window_suffix, chr_suffix, transform_suffix, ".pdf"))
        ggsave(pdf_file, p, width=opts$width, height=opts$height)
        cat("Saved:", pdf_file, "\n")
    }
    
    if ("svg" %in% format_parts) {
        svg_file <- file.path(opts$`output-dir`, paste0(prefix, "pbe_", trio_safe, window_suffix, chr_suffix, transform_suffix, ".svg"))
        ggsave(svg_file, p, width=opts$width, height=opts$height, device="svg")
        cat("Saved:", svg_file, "\n")
    }
}

# Combine and save statistics
if (length(stats_results) > 0) {
    all_stats <- bind_rows(stats_results) %>%
        select(sample_trio, window_size, chr, median, mean, sd,
               q5, q1, q0.2, q95, q99, q99.8) %>%
        arrange(sample_trio, window_size, chr)
    
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
    
    unique_trios_stats <- if ("sample_trio" %in% colnames(all_stats)) {
        sort(unique(all_stats$sample_trio[!is.na(all_stats$sample_trio)]))
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
    
    # Trio suffix (for PBE, we have trios instead of pairs)
    if (length(unique_trios_stats) > 1) {
        # Clean trio names
        trio_names_clean <- gsub(":", "_", unique_trios_stats)
        trio_names_clean <- gsub("[^A-Za-z0-9_]", "_", trio_names_clean)
        
        # If more than 3, just show first one and count
        if (length(trio_names_clean) > 3) {
            trio_list <- paste0(trio_names_clean[1], "-and", length(trio_names_clean)-1, "more")
            trio_suffix <- paste0("_trios", trio_list)
        } else {
            trio_list <- paste(trio_names_clean, collapse="-")
            trio_suffix <- paste0("_trios", trio_list)
        }
    } else if (length(unique_trios_stats) == 1) {
        trio_safe <- gsub(":", "_", unique_trios_stats[1])
        trio_safe <- gsub("[^A-Za-z0-9_]", "_", trio_safe)
        # Truncate single trio name if too long
        if (nchar(trio_safe) > 50) {
            trio_safe <- substr(trio_safe, 1, 50)
        }
        trio_suffix <- paste0("_trio", trio_safe)
    } else {
        trio_suffix <- ""
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
    
    stats_file <- file.path(opts$`output-dir`, paste0(prefix, "pbe_statistics", window_suffix, chr_suffix, trio_suffix, transform_suffix_stats, ".csv"))
    write_csv(all_stats, stats_file)
    cat("\nSaved statistics to:", stats_file, "\n")
} else {
    cat("\nWarning: No statistics were calculated\n")
}

cat("\nDone!\n")
