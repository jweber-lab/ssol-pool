#!/bin/bash

###############################################################################
# plot_pbe.sh
# 
# Bash wrapper for plot_pbe.R to create ggplot2 plots of PBE statistics
# with chromosome stripes and summary statistics.
#
# Author: Based on plot_fst.sh
# Usage: See README.md or run with --help
###############################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Default parameters
INPUT_DIR=""
OUTPUT_DIR=""
REFERENCE_GENOME=""
PANEL_BY="both"
PLOT_FORMAT="png"
WIDTH=12
HEIGHT=8
DPI=""
FILE_PREFIX=""
SAMPLE_TRIOS=""
CHROMOSOME=""
WINDOW_SIZE=""
TOP_N_CHROMOSOMES=""
MIN_CHROMOSOME_LENGTH=""
TRANSFORM="none"
ASINH_SCALE=""
INPUT_FORMAT="auto"
Y_VALUE="value"
PLOT_STYLE="line"
OVERLAY_TRIOS=false
VERBOSE=false
RSCRIPT=""
DRY_RUN=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  --input-dir DIR             Directory containing PBE TSV/CSV or HDF5 files (from calculate_pbe.sh or collate)

Optional:
  --input-format FORMAT       Input format: tsv, hdf5, or auto (default: auto)
  --y-value VALUE            Y-axis: value, rank, or quantile (default: value; rank/quantile require HDF5)
  --plot-style STYLE          Plot style: line or line_points (default: line)
  --overlay-trios             Plot all trios in one panel per window (color = trio)
  --output-dir DIR            Output directory for plots (default: ./)
  --reference-genome FILE     Reference genome FASTA file (optional, for chromosome lengths)
  --panel-by OPTION          Paneling option: window, trio, both, or none (default: both) 
  --plot-format FORMAT       Plot format: png, pdf, svg, both, or all (default: png)
  --width N                  Plot width in inches (default: 12)
  --height N                 Plot height in inches (default: 8)
  --dpi N                    DPI for PNG (default: 300)
  --file-prefix PREFIX       Prefix for output files (default: no prefix)
  --sample-trios LIST        Comma-separated list of sample trios to plot (format: pop1:pop2:pop3, default: all trios)
  --chromosome CHR           Single chromosome/scaffold to plot (default: all chromosomes)
  --window-size NUMBER       Single window size to plot (default: all window sizes)
  --top-n-chromosomes N      Plot only the N longest chromosomes/scaffolds (optional)
  --min-chromosome-length N  Plot only chromosomes/scaffolds longer than N bp (optional)
  --transform TRANSFORM      Y-axis transformation: none, asinh (centered on median), or log (default: none)
  --asinh-scale NUMBER       Scale factor for asinh transformation (default: use global standard deviation)
  --verbose                  Enable verbose output for debugging
  --rscript PATH             Path to plot_pbe.R (default: same directory as this script)
  --dry-run                  Preview commands without executing (dry-run mode)
  -h, --help                 Show this help message

Examples:
  # Basic usage:
  $0 \\
    --input-dir /path/to/pbe/output \\
    --output-dir /path/to/plots

  # With reference genome and custom paneling:
  $0 \\
    --input-dir /path/to/pbe/output \\
    --output-dir /path/to/plots \\
    --reference-genome /path/to/reference.fa \\
    --panel-by window \\
    --plot-format all

  # Plot specific sample trios:
  $0 \\
    --input-dir /path/to/pbe/output \\
    --output-dir /path/to/plots \\
    --sample-trios "Echo:Kjer:Cheney,Echo:Myv:Kjer"

  # Plot a specific chromosome:
  $0 \\
    --input-dir /path/to/pbe/output \\
    --output-dir /path/to/plots \\
    --chromosome ptg000624l

Output files:
  - {output_dir}/pbe_{sample_trio}.png/.pdf/.svg: Plots for each sample trio (format depends on --plot-format)
  - {output_dir}/pbe_{sample_trio}_{chromosome}.png/.pdf/.svg: If --chromosome specified
  - {output_dir}/pbe_statistics.csv: Summary statistics (median, mean, q95, q99, q99.8)
  - {output_dir}/pbe_statistics_{chromosome}.csv: If --chromosome specified

Note: PBE is symmetric with respect to pop2 and pop3. If pop1:pop2:pop3 is requested,
      the script will also check for pop1:pop3:pop2 and use whichever is found.

EOF
}

# Function to log messages
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_dry_run() {
    echo -e "${YELLOW}[DRY-RUN]${NC} $1"
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        return 1
    fi
    return 0
}

# Function to check if file exists
check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        exit 1
    fi
}

# Function to check if directory exists
check_dir() {
    if [[ ! -d "$1" ]]; then
        log_error "Directory not found: $1"
        exit 1
    fi
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --reference-genome)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        --panel-by)
            PANEL_BY="$2"
            shift 2
            ;;
        --plot-format)
            PLOT_FORMAT="$2"
            shift 2
            ;;
        --width)
            WIDTH="$2"
            shift 2
            ;;
        --height)
            HEIGHT="$2"
            shift 2
            ;;
        --dpi)
            DPI="$2"
            shift 2
            ;;
        --file-prefix)
            FILE_PREFIX="$2"
            shift 2
            ;;
        --sample-trios)
            SAMPLE_TRIOS="$2"
            shift 2
            ;;
        --chromosome)
            CHROMOSOME="$2"
            shift 2
            ;;
        --window-size)
            WINDOW_SIZE="$2"
            shift 2
            ;;
        --top-n-chromosomes)
            TOP_N_CHROMOSOMES="$2"
            shift 2
            ;;
        --min-chromosome-length)
            MIN_CHROMOSOME_LENGTH="$2"
            shift 2
            ;;
        --transform)
            TRANSFORM="$2"
            shift 2
            ;;
        --asinh-scale)
            ASINH_SCALE="$2"
            shift 2
            ;;
        --input-format)
            INPUT_FORMAT="$2"
            shift 2
            ;;
        --y-value)
            Y_VALUE="$2"
            shift 2
            ;;
        --plot-style)
            PLOT_STYLE="$2"
            shift 2
            ;;
        --overlay-trios)
            OVERLAY_TRIOS=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --rscript)
            RSCRIPT="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" ]]; then
    log_error "Missing required argument: --input-dir"
    usage
    exit 1
fi

# Validate panel-by option
if [[ "$PANEL_BY" != "window" ]] && [[ "$PANEL_BY" != "trio" ]] && \
   [[ "$PANEL_BY" != "both" ]] && [[ "$PANEL_BY" != "none" ]]; then
    log_error "Invalid --panel-by option: $PANEL_BY"
    log_error "Valid options: window, trio, both, none"
    exit 1
fi

# Validate plot format (allow comma-separated or single values)
IFS=',' read -ra FORMAT_ARRAY <<< "$PLOT_FORMAT"
valid_formats=("png" "pdf" "svg" "both" "all")
for fmt in "${FORMAT_ARRAY[@]}"; do
    fmt=$(echo "$fmt" | xargs)  # trim whitespace
    if [[ ! " ${valid_formats[@]} " =~ " ${fmt} " ]]; then
        log_error "Invalid --plot-format option: $fmt"
        log_error "Valid options: png, pdf, svg, both, all (comma-separated)"
        exit 1
    fi
done

# Set default output directory
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="./"
fi

# Find R script
if [[ -z "$RSCRIPT" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    RSCRIPT="${SCRIPT_DIR}/plot_pbe.R"
fi

# Check if R script exists
if [[ ! -f "$RSCRIPT" ]]; then
    log_error "R script not found: $RSCRIPT"
    exit 1
fi

# Check if R is available
if ! check_command "Rscript"; then
    log_error "Rscript not found. Please install R."
    exit 1
fi

# Check input directory
if [[ "$DRY_RUN" == false ]]; then
    check_dir "$INPUT_DIR"
fi

# Check for PBE TSV/CSV files
if [[ "$DRY_RUN" == false ]]; then
    tsv_count=$(find "$INPUT_DIR" -name "*pbe*.tsv" -type f | wc -l)
    csv_count=$(find "$INPUT_DIR" -name "*pbe*.csv" -type f | wc -l)
    total_count=$((tsv_count + csv_count))
    if [[ $total_count -eq 0 ]]; then
        log_warn "No PBE TSV/CSV files found in: $INPUT_DIR"
        log_warn "Expected files matching pattern: *pbe*.tsv or *pbe*.csv"
    else
        log "Found $total_count PBE file(s) (tsv: $tsv_count, csv: $csv_count)"
    fi
fi

# Check reference genome if provided
if [[ -n "$REFERENCE_GENOME" ]] && [[ "$DRY_RUN" == false ]]; then
    if [[ ! -f "$REFERENCE_GENOME" ]]; then
        log_warn "Reference genome file not found: $REFERENCE_GENOME"
        log_warn "Will infer chromosome lengths from data"
    fi
fi

# Create output directory
if [[ "$DRY_RUN" == false ]]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Build R command
R_CMD=(
    Rscript
    "$RSCRIPT"
    --input-dir "$INPUT_DIR"
    --output-dir "$OUTPUT_DIR"
    --panel-by "$PANEL_BY"
    --plot-format "$PLOT_FORMAT"
    --width "$WIDTH"
    --height "$HEIGHT"
)

if [[ -n "$DPI" ]]; then
    R_CMD+=(--dpi "$DPI")
fi

if [[ -n "$REFERENCE_GENOME" ]]; then
    R_CMD+=(--reference-genome "$REFERENCE_GENOME")
fi

if [[ -n "$FILE_PREFIX" ]]; then
    R_CMD+=(--file-prefix "$FILE_PREFIX")
fi

if [[ -n "$SAMPLE_TRIOS" ]]; then
    R_CMD+=(--sample-trios "$SAMPLE_TRIOS")
fi

if [[ -n "$CHROMOSOME" ]]; then
    R_CMD+=(--chromosome "$CHROMOSOME")
fi

if [[ -n "$WINDOW_SIZE" ]]; then
    R_CMD+=(--window-size "$WINDOW_SIZE")
fi

if [[ -n "$TOP_N_CHROMOSOMES" ]]; then
    R_CMD+=(--top-n-chromosomes "$TOP_N_CHROMOSOMES")
fi

if [[ -n "$MIN_CHROMOSOME_LENGTH" ]]; then
    R_CMD+=(--min-chromosome-length "$MIN_CHROMOSOME_LENGTH")
fi

if [[ -n "$TRANSFORM" ]]; then
    R_CMD+=(--transform "$TRANSFORM")
fi

if [[ -n "$ASINH_SCALE" ]]; then
    R_CMD+=(--asinh-scale "$ASINH_SCALE")
fi

if [[ -n "$INPUT_FORMAT" ]]; then
    R_CMD+=(--input-format "$INPUT_FORMAT")
fi

if [[ -n "$Y_VALUE" ]]; then
    R_CMD+=(--y-value "$Y_VALUE")
fi

if [[ -n "$PLOT_STYLE" ]]; then
    R_CMD+=(--plot-style "$PLOT_STYLE")
fi

if [[ "$OVERLAY_TRIOS" == true ]]; then
    R_CMD+=(--overlay-trios)
fi

if [[ "$VERBOSE" == true ]]; then
    R_CMD+=(--verbose)
fi

# Run R script
if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN MODE: Commands will be previewed but not executed"
    log_dry_run ""
    log_dry_run "Would execute:"
    log_dry_run "  ${R_CMD[*]}"
    log_dry_run ""
    log_dry_run "Input directory: $INPUT_DIR"
    log_dry_run "Output directory: $OUTPUT_DIR"
    log_dry_run "Panel by: $PANEL_BY"
    log_dry_run "Plot format: $PLOT_FORMAT"
    if [[ -n "$SAMPLE_TRIOS" ]]; then
        log_dry_run "Sample trios: $SAMPLE_TRIOS"
    fi
    if [[ -n "$CHROMOSOME" ]]; then
        log_dry_run "Chromosome: $CHROMOSOME"
    fi
    if [[ -n "$WINDOW_SIZE" ]]; then
        log_dry_run "Window size: $WINDOW_SIZE"
    fi
    if [[ -n "$TOP_N_CHROMOSOMES" ]]; then
        log_dry_run "Top N chromosomes: $TOP_N_CHROMOSOMES"
    fi
    if [[ -n "$MIN_CHROMOSOME_LENGTH" ]]; then
        log_dry_run "Min chromosome length: $MIN_CHROMOSOME_LENGTH"
    fi
    if [[ -n "$TRANSFORM" ]]; then
        log_dry_run "Transform: $TRANSFORM"
    fi
    if [[ -n "$ASINH_SCALE" ]]; then
        log_dry_run "Asinh scale: $ASINH_SCALE"
    fi
else
    log "Running plot_pbe.R..."
    log "  Input directory: $INPUT_DIR"
    log "  Output directory: $OUTPUT_DIR"
    log "  Panel by: $PANEL_BY"
    log "  Plot format: $PLOT_FORMAT"
    if [[ -n "$SAMPLE_TRIOS" ]]; then
        log "  Sample trios: $SAMPLE_TRIOS"
    fi
    if [[ -n "$CHROMOSOME" ]]; then
        log "  Chromosome: $CHROMOSOME"
    fi
    if [[ -n "$WINDOW_SIZE" ]]; then
        log "  Window size: $WINDOW_SIZE"
    fi
    if [[ -n "$TOP_N_CHROMOSOMES" ]]; then
        log "  Top N chromosomes: $TOP_N_CHROMOSOMES"
    fi
    if [[ -n "$MIN_CHROMOSOME_LENGTH" ]]; then
        log "  Min chromosome length: $MIN_CHROMOSOME_LENGTH"
    fi
    if [[ -n "$TRANSFORM" ]]; then
        log "  Transform: $TRANSFORM"
    fi
    if [[ -n "$ASINH_SCALE" ]]; then
        log "  Asinh scale: $ASINH_SCALE"
    fi
    
    "${R_CMD[@]}" || {
        log_error "R script failed"
        exit 1
    }
    
    log "Plotting complete!"
    log "Output directory: $OUTPUT_DIR"
fi
