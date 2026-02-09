#!/bin/bash

###############################################################################
# calculate_fst.sh
# 
# Calculate per-locus FST (Fixation Index) between sample pairs from poolseq
# data using grenedalf.
#
# FST measures genetic differentiation between populations:
# - FST ranges from 0 (no differentiation) to 1 (complete differentiation)
# - Uses pool-seq corrected estimators to account for sampling bias
# - Computes FST at each position (per-locus) for pairwise comparisons
#
# Author: Based on workflow by JW
# Usage: See README.md or run with --help
###############################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Default parameters
BAM_FILES=()
SAMPLE_INFO_CSV=""
GRENEDALF=""
OUTPUT_DIR=""
FST_METHOD="unbiased-hudson"
WINDOW_TYPE="single"  # single (per-locus) or interval (windowed)
WINDOW_AVERAGE_POLICY="valid-loci"
WINDOW_SIZES=()       # Array for multi-scale analysis (can be specified multiple times)
STEP_SIZES=()          # Array for multi-scale analysis (optional, auto-calculated as window/2 if not provided)
MIN_FREQUENCY=0.01
MIN_COUNT_FILTER=0
MIN_COVERAGE=10
MAX_COVERAGE=1000
MIN_COUNT=2
MIN_TOTAL_READ_DEPTH=10
FILTER_TOTAL_ONLY_BIALLELIC_SNPS=false
FILTER_MASK_TOTAL_FASTA=""
FILTER_MASK_TOTAL_BED=""
POOL_SIZE=50
THREADS=1
DRY_RUN=false
REFERENCE_GENOME=""
FILE_PREFIX=""  # Optional prefix for output files (default: no prefix)
COMPARAND=""
COMPARAND_LIST=""
PARALLEL_WINDOWS=false  # Run multiple window/step combinations in parallel
PARALLEL_MAX_JOBS=0     # Maximum concurrent window jobs (0 = number of CPU cores)
OUTPUT_SEPARATOR="comma"  # Output separator for grenedalf (comma -> .csv, tab -> .tsv)
OVERWRITE_OUTPUT=false  # Allow grenedalf to overwrite existing output files

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Global log file variable (set after argument parsing)
SCRIPT_NAME=$(basename "$0")
LOG_FILE=""
LOG_REDIRECT_ACTIVE=false  # Flag to indicate if stdout/stderr redirection is active

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Input options (choose one):
  -b, --bam FILE              BAM file(s) - can be specified multiple times (minimum 2 required)
  -i, --sample-info FILE      CSV file with sample information (see sample-info.csv)
                              Columns: sample_name, read1_file, read2_file, pool_size, bam_file
                              Note: read1_file and read2_file columns are ignored

Tool options:
  -g, --grenedalf PATH        Path to grenedalf executable (default: check PATH)

Required:
  -o, --output-dir DIR        Output directory for results

FST-specific options:
  --method METHOD             FST estimator method (default: unbiased-hudson)
                              Options: unbiased-nei, unbiased-hudson, kofler, karlsson
  --window-type TYPE          Window type: single (per-locus, default) or interval (windowed)
  --window-size N             Window size in bp for windowed mode (can be specified multiple times for multi-scale)
                               Default: 10000 (single scale). Multiple values enable multi-scale mode.
  --step-size N               Step size in bp for windowed mode (can be specified multiple times for multi-scale)
                               Default: 5000 (single scale). If fewer step sizes than window sizes,
                               remaining step sizes are auto-calculated as window/2.
  --window-average-policy     Window averaging policy (default: valid-loci)
                              Options: window-length, available-loci, valid-loci, valid-snps, sum, provided-loci
  --min-frequency FLOAT       Minimum allele frequency cutoff (default: 0.01 for 1%)
  --min-count-filter INT      Minimum allele count cutoff, alternative to frequency (default: 0 = disabled)

Filtering options:
  --min-coverage N            Minimum coverage per site per sample (default: 10)
  --max-coverage N            Maximum coverage per site per sample (default: 1000)
  --min-count N               Minimum allele count per sample (default: 2)
  --min-total-read-depth N   Minimum total read depth across all samples (default: 10)
  --filter-total-only-biallelic-snps  Filter to only biallelic SNPs (default: false)
  --filter-mask-total-fasta FILE  FASTA mask file for masking regions (vcftools format, optional)
  --filter-mask-total-bed FILE    BED file for masking regions (optional)

Comparison options:
  --comparand SAMPLE          Compute FST between this sample and all others (optional)
  --comparand-list FILE       File with sample pairs to compute FST for (one pair per line, optional)

Other options:
  --reference-genome FILE     Reference genome FASTA file (recommended for chromosome ordering)
  --file-prefix PREFIX        Optional prefix for output files (default: no prefix)
  --separator-char SEP        Output separator: comma, tab, space, semicolon (default: comma)
  --overwrite-output          Allow overwriting existing output files (default: false)
  --pool-size N               Default pool size if not in CSV (default: 50)
  -t, --threads N             Number of threads per grenedalf call (default: 1)
  --parallel-windows          Run multiple window/step combinations in parallel (default: false)
  --parallel-max-jobs N       Maximum concurrent window jobs when using --parallel-windows (default: number of CPU cores)
  --dry-run                   Preview commands without executing (dry-run mode)
  -h, --help                  Show this help message

Examples:
  # Using sample info CSV file (all pairwise comparisons):
  $0 \\
    --sample-info sample-info.csv \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output \\
    --reference-genome /path/to/reference.fa

  # Using BAM files directly (per-locus FST):
  $0 \\
    --bam /path/to/sample1.bam \\
    --bam /path/to/sample2.bam \\
    --bam /path/to/sample3.bam \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output

  # Multi-scale windowed FST:
  $0 \\
    --sample-info sample-info.csv \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output \\
    --window-type interval \\
    --window-size 1000 \\
    --window-size 5000 \\
    --window-size 10000 \\
    --step-size 500 \\
    --step-size 2500 \\
    --step-size 5000 \\
    --pool-size 200 \\
    --reference-genome /path/to/reference.fa

  # Compute FST between one sample and all others:
  $0 \\
    --sample-info sample-info.csv \\
    --comparand Echo_Pool_S1 \\
    --output-dir /path/to/output \\
    --reference-genome /path/to/reference.fa

Output files:
  - {output_dir}/fst_single.{csv|tsv}: Per-locus FST values for each sample pair (when --window-type single)
    (Grenedalf writes .csv; when --separator-char tab, we rename to .tsv)
  - {output_dir}/fst_w{WINDOW}_s{STEP}.{csv|tsv}: Windowed FST values for each sample pair
    (e.g., fst_w10000_s5000.tsv for window=10000, step=5000)
  - Columns: chromosome, position, sample1_sample2_fst, etc.
  - Logs: {output_dir}/log/{prefix_}calculate_fst_YYYYmmdd_HHMMSS.log
  
Note: Sample names are taken from the sample_name column in the CSV file (if provided),
      or derived from BAM file basenames if using --bam directly.

Note: FST requires at least 2 samples. With 3+ samples, all pairwise comparisons are computed by default.

EOF
}

# Function to log messages (to both console and log file)
# Note: If LOG_REDIRECT_ACTIVE is true, stdout/stderr are redirected via tee, so we don't write directly to avoid duplicates
log() {
    local message="${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
    echo -e "$message"
    # Also write directly to log file (without colors) if LOG_FILE is set and redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[$(date +'%Y-%m-%d %H:%M:%S')] $1"
        echo "$clean_message" >> "$LOG_FILE"
    fi
}

log_error() {
    local message="${RED}[ERROR]${NC} $1"
    echo -e "$message" >&2
    # Also write directly to log file if redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[ERROR] $1"
        echo "$clean_message" >> "$LOG_FILE"
    fi
}

log_warn() {
    local message="${YELLOW}[WARN]${NC} $1"
    echo -e "$message"
    # Also write directly to log file if redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[WARN] $1"
        echo "$clean_message" >> "$LOG_FILE"
    fi
}

log_dry_run() {
    local message="${YELLOW}[DRY-RUN]${NC} $1"
    echo -e "$message"
    # Also write directly to log file if redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[DRY-RUN] $1"
        echo "$clean_message" >> "$LOG_FILE"
    fi
}

# Function to execute command or print in dry-run mode
dry_run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: $*"
        return 0
    else
        "$@"
    fi
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        return 1
    fi
    return 0
}

# Function to get number of CPU cores
get_cpu_cores() {
    if command -v nproc &> /dev/null; then
        nproc
    elif [[ -f /proc/cpuinfo ]]; then
        grep -c ^processor /proc/cpuinfo
    elif [[ "$(uname)" == "Darwin" ]]; then
        sysctl -n hw.ncpu
    else
        echo "8"  # Default fallback
    fi
}

# Function to check if file exists
check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        exit 1
    fi
}

# Function to check if directory exists, create if not
check_dir() {
    if [[ ! -d "$1" ]]; then
        if [[ "$DRY_RUN" == true ]]; then
            log_dry_run "Would create directory: $1"
        else
            log "Creating directory: $1"
            mkdir -p "$1"
        fi
    fi
}

# Function to resolve path relative to CSV file location
resolve_path_relative_to_csv() {
    local csv_file="$1"
    local file_path="$2"
    
    # If path is already absolute, return as-is
    if [[ "$file_path" =~ ^/ ]]; then
        echo "$file_path"
        return 0
    fi
    
    # Get CSV directory
    local csv_dir
    csv_dir=$(cd "$(dirname "$csv_file")" && pwd)
    
    # Resolve path relative to CSV directory
    local resolved_path
    resolved_path=$(cd "$csv_dir" && cd "$(dirname "$file_path")" 2>/dev/null && pwd)/$(basename "$file_path") 2>/dev/null || echo "${csv_dir}/${file_path}"
    
    # Normalize path (remove .. and .)
    echo "$resolved_path" | sed 's|/\./|/|g; s|/\.\./|/|g; s|^\./||'
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--bam)
            BAM_FILES+=("$2")
            shift 2
            ;;
        -i|--sample-info)
            SAMPLE_INFO_CSV="$2"
            shift 2
            ;;
        -g|--grenedalf)
            GRENEDALF="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --method)
            FST_METHOD="$2"
            shift 2
            ;;
        --window-type)
            WINDOW_TYPE="$2"
            shift 2
            ;;
        --window-size)
            WINDOW_SIZES+=("$2")
            shift 2
            ;;
        --step-size)
            STEP_SIZES+=("$2")
            shift 2
            ;;
        --window-average-policy)
            WINDOW_AVERAGE_POLICY="$2"
            shift 2
            ;;
        --min-frequency)
            MIN_FREQUENCY="$2"
            shift 2
            ;;
        --min-count-filter)
            MIN_COUNT_FILTER="$2"
            shift 2
            ;;
        --min-coverage)
            MIN_COVERAGE="$2"
            shift 2
            ;;
        --max-coverage)
            MAX_COVERAGE="$2"
            shift 2
            ;;
        --min-count)
            MIN_COUNT="$2"
            shift 2
            ;;
        --min-total-read-depth)
            MIN_TOTAL_READ_DEPTH="$2"
            shift 2
            ;;
        --filter-total-only-biallelic-snps)
            FILTER_TOTAL_ONLY_BIALLELIC_SNPS=true
            shift
            ;;
        --filter-mask-total-fasta)
            FILTER_MASK_TOTAL_FASTA="$2"
            shift 2
            ;;
        --filter-mask-total-bed)
            FILTER_MASK_TOTAL_BED="$2"
            shift 2
            ;;
        --comparand)
            COMPARAND="$2"
            shift 2
            ;;
        --comparand-list)
            COMPARAND_LIST="$2"
            shift 2
            ;;
        --reference-genome)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        --file-prefix)
            FILE_PREFIX="$2"
            shift 2
            ;;
        --separator-char)
            OUTPUT_SEPARATOR="$2"
            shift 2
            ;;
        --overwrite-output)
            OVERWRITE_OUTPUT=true
            shift
            ;;
        --pool-size)
            POOL_SIZE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --parallel-windows)
            PARALLEL_WINDOWS=true
            shift
            ;;
        --parallel-max-jobs)
            PARALLEL_MAX_JOBS="$2"
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
if [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing required argument: --output-dir"
    usage
    exit 1
fi

# Initialize log file (create output/log directory)
check_dir "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
check_dir "$LOG_DIR"
TIMESTAMP=$(date +'%Y%m%d_%H%M%S')
if [[ -n "$FILE_PREFIX" ]]; then
    LOG_FILE="${LOG_DIR}/${FILE_PREFIX}_${SCRIPT_NAME%.*}_${TIMESTAMP}.log"
else
    LOG_FILE="${LOG_DIR}/${SCRIPT_NAME%.*}_${TIMESTAMP}.log"
fi

if [[ "$DRY_RUN" == false ]]; then
    touch "$LOG_FILE"
    {
        echo "=========================================="
        echo "$SCRIPT_NAME"
        echo "Started: $(date +'%Y-%m-%d %H:%M:%S')"
        echo "Log file: $LOG_FILE"
        echo "=========================================="
    } >> "$LOG_FILE"
    log "Log file: $LOG_FILE"

    # Redirect all stdout and stderr to log file (in addition to console)
    # Use file descriptors to preserve original stdout/stderr for console output
    exec 3>&1 4>&2
    # Tee stdout to both console (fd 3) and log file
    exec 1> >(tee -a "$LOG_FILE" >&3)
    # Redirect stderr to stdout (so it also gets teed to console and log file)
    exec 2>&1
    # Set flag so logging functions know redirection is active (avoid duplicate writes)
    LOG_REDIRECT_ACTIVE=true
else
    log_dry_run "Would create log file: $LOG_FILE"
fi

# Check for grenedalf
if [[ -n "$GRENEDALF" ]]; then
    if [[ ! -f "$GRENEDALF" ]] && ! check_command "$GRENEDALF"; then
        log_error "Grenedalf executable not found: $GRENEDALF"
        exit 1
    fi
elif check_command "grenedalf"; then
    GRENEDALF="grenedalf"
    log "Found grenedalf in PATH"
else
    log_error "Grenedalf not found. Please install grenedalf or specify with --grenedalf"
    exit 1
fi

# Validate FST method
case "$FST_METHOD" in
    unbiased-nei|unbiased-hudson|kofler|karlsson)
        ;;
    *)
        log_error "Invalid FST method: $FST_METHOD"
        log_error "Valid methods: unbiased-nei, unbiased-hudson, kofler, karlsson"
        exit 1
        ;;
esac

# Validate output separator
case "$OUTPUT_SEPARATOR" in
    comma|tab|space|semicolon)
        ;;
    *)
        log_error "Invalid separator: $OUTPUT_SEPARATOR"
        log_error "Valid separators: comma, tab, space, semicolon"
        exit 1
        ;;
esac

# Determine output extension based on separator
OUTPUT_EXTENSION="csv"
if [[ "$OUTPUT_SEPARATOR" == "tab" ]]; then
    OUTPUT_EXTENSION="tsv"
fi

# Process sample info CSV if provided
# Initialize arrays (needed even when using --bam)
declare -a SAMPLE_NAMES
declare -a SAMPLE_BAMS
declare -a SAMPLE_POOL_SIZES

if [[ -n "$SAMPLE_INFO_CSV" ]]; then
    if [[ ! -f "$SAMPLE_INFO_CSV" ]]; then
        log_error "Sample info CSV file not found: $SAMPLE_INFO_CSV"
        exit 1
    fi
    
    log "Reading sample information from: $SAMPLE_INFO_CSV"
    
    # Read CSV file (skip header row and comment lines)
    line_num=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        line_num=$((line_num + 1))
        
        # Skip empty lines and comment lines
        [[ -z "$line" ]] && continue
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        
        # Parse CSV line (simple parsing, assumes no quoted commas in values)
        IFS=',' read -r -a fields <<< "$line"
        
        # Trim whitespace from fields
        for i in "${!fields[@]}"; do
            fields[$i]=$(echo "${fields[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        done
        
        # Skip header row (first column is "sample_name" or "Sample" or similar)
        first_col="${fields[0]:-}"
        case "$first_col" in
            [Ss]ample_name|[Ss]ample|[Bb]am_file) continue ;;
        esac
        
        # Expected columns: sample_name, read1_file, read2_file, pool_size, bam_file
        # Note: read1_file and read2_file are ignored (used by process_poolseq.sh)
        if [[ ${#fields[@]} -lt 5 ]]; then
            log_warn "Skipping line $line_num: insufficient columns (expected at least 5: sample_name, read1_file, read2_file, pool_size, bam_file)"
            continue
        fi
        
        sample_name="${fields[0]}"
        read1_file="${fields[1]}"  # Ignored, but present in CSV
        read2_file="${fields[2]}"  # Ignored, but present in CSV
        pool_size="${fields[3]}"
        bam_file="${fields[4]}"
        
        # Resolve BAM file path relative to CSV location
        bam_file=$(resolve_path_relative_to_csv "$SAMPLE_INFO_CSV" "$bam_file")
        
        # Validate required fields
        if [[ -z "$sample_name" ]] || [[ -z "$bam_file" ]] || [[ -z "$pool_size" ]]; then
            log_warn "Skipping line $line_num: missing required fields"
            continue
        fi
        
        # Check if BAM file exists (unless dry-run)
        if [[ "$DRY_RUN" == false ]] && [[ ! -f "$bam_file" ]]; then
            log_warn "BAM file not found for sample $sample_name: $bam_file"
            continue
        fi
        
        SAMPLE_NAMES+=("$sample_name")
        SAMPLE_BAMS+=("$bam_file")
        SAMPLE_POOL_SIZES+=("$pool_size")
        
        log "  Sample: $sample_name, BAM: $bam_file, Pool size: $pool_size"
    done < "$SAMPLE_INFO_CSV"
    
    if [[ ${#SAMPLE_NAMES[@]} -eq 0 ]]; then
        log_error "No valid samples found in CSV file"
        exit 1
    fi
elif [[ ${#BAM_FILES[@]} -gt 0 ]]; then
    # Use BAM files provided on command line
    for bam in "${BAM_FILES[@]}"; do
        if [[ "$DRY_RUN" == false ]]; then
            check_file "$bam"
        fi
        # SAMPLE_NAMES will be empty, will be derived from BAM filename later
        SAMPLE_NAMES+=("")
        SAMPLE_BAMS+=("$bam")
        SAMPLE_POOL_SIZES+=("$POOL_SIZE")
    done
else
    log_error "Must provide either BAM file(s) with --bam or sample info CSV with --sample-info"
    usage
    exit 1
fi

# Validate minimum 2 samples for FST
if [[ ${#SAMPLE_BAMS[@]} -lt 2 ]]; then
    log_error "FST requires at least 2 samples, but only ${#SAMPLE_BAMS[@]} provided"
    exit 1
fi

# Validate window type
if [[ "$WINDOW_TYPE" != "single" ]] && [[ "$WINDOW_TYPE" != "interval" ]]; then
    log_error "Invalid window type: $WINDOW_TYPE (must be 'single' or 'interval')"
    exit 1
fi

# Validate mask files if provided
if [[ -n "$FILTER_MASK_TOTAL_FASTA" ]] && [[ "$DRY_RUN" == false ]]; then
    if [[ ! -f "$FILTER_MASK_TOTAL_FASTA" ]]; then
        log_error "FASTA mask file not found: $FILTER_MASK_TOTAL_FASTA"
        exit 1
    fi
fi

if [[ -n "$FILTER_MASK_TOTAL_BED" ]] && [[ "$DRY_RUN" == false ]]; then
    if [[ ! -f "$FILTER_MASK_TOTAL_BED" ]]; then
        log_error "BED mask file not found: $FILTER_MASK_TOTAL_BED"
        exit 1
    fi
fi

# Warn if both mask files are provided (grenedalf allows only one)
if [[ -n "$FILTER_MASK_TOTAL_FASTA" ]] && [[ -n "$FILTER_MASK_TOTAL_BED" ]]; then
    log_warn "Both FASTA and BED mask files specified. Grenedalf will use only one (FASTA takes precedence if both are valid)"
fi

# Warn if window/step size parameters are provided with single mode
if [[ "$WINDOW_TYPE" == "single" ]]; then
    if [[ ${#WINDOW_SIZES[@]} -gt 0 ]] || [[ ${#STEP_SIZES[@]} -gt 0 ]]; then
        log_warn "Window/step size parameters are ignored when using --window-type single (per-locus FST)"
        if [[ ${#WINDOW_SIZES[@]} -gt 0 ]]; then
            log_warn "  --window-size will be ignored (${#WINDOW_SIZES[@]} value(s) provided)"
        fi
        if [[ ${#STEP_SIZES[@]} -gt 0 ]]; then
            log_warn "  --step-size will be ignored (${#STEP_SIZES[@]} value(s) provided)"
        fi
        log_warn "  Use --window-type interval to enable windowed FST analysis"
    fi
fi

# Parse window sizes and step sizes for multi-scale analysis (only for windowed mode)
declare -a WINDOW_SIZE_ARRAY
declare -a STEP_SIZE_ARRAY
USE_MULTI_SCALE=false

if [[ "$WINDOW_TYPE" == "interval" ]]; then
    # Check if multiple window sizes were provided
    if [[ ${#WINDOW_SIZES[@]} -gt 1 ]]; then
        # Multi-scale mode
        USE_MULTI_SCALE=true
        WINDOW_SIZE_ARRAY=("${WINDOW_SIZES[@]}")
        
        # Initialize STEP_SIZE_ARRAY (empty if no step sizes provided)
        STEP_SIZE_ARRAY=()
        
        # Check if step sizes were provided
        if [[ ${#STEP_SIZES[@]} -gt 0 ]]; then
            STEP_SIZE_ARRAY=("${STEP_SIZES[@]}")
            if [[ ${#STEP_SIZE_ARRAY[@]} -ne ${#WINDOW_SIZE_ARRAY[@]} ]]; then
                log_warn "Number of step sizes (${#STEP_SIZE_ARRAY[@]}) does not match number of window sizes (${#WINDOW_SIZE_ARRAY[@]})"
                log_warn "Will auto-calculate missing step sizes as window/2"
            fi
        fi
        
        # Auto-calculate step sizes for any missing ones
        for i in "${!WINDOW_SIZE_ARRAY[@]}"; do
            if [[ ${#STEP_SIZE_ARRAY[@]} -le $i ]] || [[ -z "${STEP_SIZE_ARRAY[$i]:-}" ]]; then
                # Auto-calculate as window/2
                window_val="${WINDOW_SIZE_ARRAY[$i]}"
                step_val=$((window_val / 2))
                STEP_SIZE_ARRAY[$i]=$step_val
            fi
        done
        
        log "Multi-scale windowed FST analysis: ${#WINDOW_SIZE_ARRAY[@]} window/step combinations"
        for i in "${!WINDOW_SIZE_ARRAY[@]}"; do
            log "  Window ${WINDOW_SIZE_ARRAY[$i]} bp, Step ${STEP_SIZE_ARRAY[$i]} bp"
        done
    elif [[ ${#WINDOW_SIZES[@]} -eq 1 ]]; then
        # Single window size provided
        WINDOW_SIZE_ARRAY=("${WINDOW_SIZES[@]}")
        if [[ ${#STEP_SIZES[@]} -eq 1 ]]; then
            STEP_SIZE_ARRAY=("${STEP_SIZES[@]}")
        elif [[ ${#STEP_SIZES[@]} -eq 0 ]]; then
            # Auto-calculate step as window/2
            window_val="${WINDOW_SIZE_ARRAY[0]}"
            step_val=$((window_val / 2))
            STEP_SIZE_ARRAY=("$step_val")
        else
            log_warn "Multiple step sizes provided but only one window size. Using first step size."
            STEP_SIZE_ARRAY=("${STEP_SIZES[0]}")
        fi
        log "Single-scale windowed FST analysis: Window ${WINDOW_SIZE_ARRAY[0]} bp, Step ${STEP_SIZE_ARRAY[0]} bp"
    else
        # No window sizes provided, use defaults
        WINDOW_SIZE_ARRAY=(10000)
        STEP_SIZE_ARRAY=(5000)
        log "Single-scale windowed FST analysis: Window 10000 bp, Step 5000 bp (defaults)"
    fi
else
    # Per-locus mode (single)
    log "Per-locus FST analysis (window-type: single)"
fi

# Set max jobs for parallel window processing
if [[ "$PARALLEL_WINDOWS" == true ]]; then
    if [[ "$PARALLEL_MAX_JOBS" -eq 0 ]]; then
        PARALLEL_MAX_JOBS=$(get_cpu_cores)
    fi
    log "Parallel window processing enabled: max $PARALLEL_MAX_JOBS concurrent jobs"
fi

# Start processing
if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN MODE: Commands will be previewed but not executed"
    log_dry_run ""
    log_dry_run "Samples to process: ${#SAMPLE_BAMS[@]}"
    for i in "${!SAMPLE_BAMS[@]}"; do
        if [[ ${#SAMPLE_NAMES[@]} -gt $i ]] && [[ -n "${SAMPLE_NAMES[$i]}" ]]; then
            log_dry_run "  Sample $((i+1)): ${SAMPLE_NAMES[$i]} (${SAMPLE_BAMS[$i]})"
        else
            log_dry_run "  Sample $((i+1)): ${SAMPLE_BAMS[$i]}"
        fi
    done
    log_dry_run ""
    log_dry_run "FST method: $FST_METHOD"
    log_dry_run "Min frequency cutoff: $MIN_FREQUENCY"
    log_dry_run "Min coverage per sample: $MIN_COVERAGE"
    log_dry_run "Max coverage per sample: $MAX_COVERAGE"
    log_dry_run "Min allele count per sample: $MIN_COUNT"
    log_dry_run "Min total read depth: $MIN_TOTAL_READ_DEPTH"
    if [[ "$FILTER_TOTAL_ONLY_BIALLELIC_SNPS" == true ]]; then
        log_dry_run "Filtering to biallelic SNPs only: enabled"
    fi
    if [[ -n "$FILTER_MASK_TOTAL_FASTA" ]]; then
        log_dry_run "FASTA mask file: $FILTER_MASK_TOTAL_FASTA"
    fi
    if [[ -n "$FILTER_MASK_TOTAL_BED" ]]; then
        log_dry_run "BED mask file: $FILTER_MASK_TOTAL_BED"
    fi
    log_dry_run ""
fi

if [[ "$WINDOW_TYPE" == "single" ]]; then
    log "Computing per-locus FST for ${#SAMPLE_BAMS[@]} samples"
else
    log "Computing windowed FST for ${#SAMPLE_BAMS[@]} samples"
fi
log "FST method: $FST_METHOD"
log "Min frequency cutoff: $MIN_FREQUENCY"
log "Min coverage per sample: $MIN_COVERAGE"
log "Max coverage per sample: $MAX_COVERAGE"
log "Min allele count per sample: $MIN_COUNT"
log "Min total read depth: $MIN_TOTAL_READ_DEPTH"
if [[ "$FILTER_TOTAL_ONLY_BIALLELIC_SNPS" == true ]]; then
    log "Filtering to biallelic SNPs only: enabled"
fi
if [[ -n "$FILTER_MASK_TOTAL_FASTA" ]]; then
    log "FASTA mask file: $FILTER_MASK_TOTAL_FASTA"
fi
if [[ -n "$FILTER_MASK_TOTAL_BED" ]]; then
    log "BED mask file: $FILTER_MASK_TOTAL_BED"
fi
log "Output separator: $OUTPUT_SEPARATOR (.$OUTPUT_EXTENSION)"
if [[ "$OVERWRITE_OUTPUT" == true ]]; then
    log_warn "Overwriting enabled: existing outputs may be replaced"
else
    log_warn "Grenedalf aborts if output files already exist (use --overwrite-output to allow overwrite)"
fi

# Normalize grenedalf output extension when using tab separator.
# Grenedalf writes .csv regardless of separator; we rename to .tsv for tab.
normalize_output_extension() {
    if [[ "$OUTPUT_SEPARATOR" != "tab" ]]; then
        return 0
    fi
    shopt -s nullglob
    for csv_file in "$OUTPUT_DIR"/*fst*.csv; do
        tsv_file="${csv_file%.csv}.tsv"
        mv -f "$csv_file" "$tsv_file"
    done
    shopt -u nullglob
}

# Handle sample name renaming (to use CSV sample names instead of BAM file paths)
# Create rename samples file if we have sample names from CSV
RENAME_SAMPLES_FILE=""
NEED_RENAME=false
if [[ ${#SAMPLE_NAMES[@]} -gt 0 ]]; then
    # Check if we have any non-empty sample names (from CSV)
    for sample_name in "${SAMPLE_NAMES[@]}"; do
        if [[ -n "$sample_name" ]]; then
            NEED_RENAME=true
            break
        fi
    done
fi

if [[ "$NEED_RENAME" == true ]]; then
    RENAME_SAMPLES_FILE="${OUTPUT_DIR}/.rename_samples.tmp"
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would create rename samples file: $RENAME_SAMPLES_FILE"
        for i in "${!SAMPLE_BAMS[@]}"; do
            if [[ ${#SAMPLE_NAMES[@]} -gt $i ]] && [[ -n "${SAMPLE_NAMES[$i]}" ]]; then
                # Use sample name from CSV
                bam_basename=$(basename "${SAMPLE_BAMS[$i]}" .bam)
                log_dry_run "  $bam_basename -> ${SAMPLE_NAMES[$i]}"
            fi
        done
    else
        > "$RENAME_SAMPLES_FILE"
        for i in "${!SAMPLE_BAMS[@]}"; do
            if [[ ${#SAMPLE_NAMES[@]} -gt $i ]] && [[ -n "${SAMPLE_NAMES[$i]}" ]]; then
                # Grenedalf derives sample names from BAM file basename (without extension)
                # Map that to our desired sample name from CSV
                bam_basename=$(basename "${SAMPLE_BAMS[$i]}" .bam)
                echo "${bam_basename},${SAMPLE_NAMES[$i]}" >> "$RENAME_SAMPLES_FILE"
            fi
        done
        log "Created rename samples file with ${#SAMPLE_NAMES[@]} mappings"
    fi
fi

# Handle pool sizes
# Check if all pool sizes are the same
ALL_SAME_POOL_SIZE=true
FIRST_POOL_SIZE="${SAMPLE_POOL_SIZES[0]}"
for pool_size in "${SAMPLE_POOL_SIZES[@]}"; do
    if [[ "$pool_size" != "$FIRST_POOL_SIZE" ]]; then
        ALL_SAME_POOL_SIZE=false
        break
    fi
done

POOL_SIZES_ARG=""
if [[ "$ALL_SAME_POOL_SIZE" == true ]]; then
    # All samples have the same pool size, use single value
    POOL_SIZES_ARG="$FIRST_POOL_SIZE"
    log "Using uniform pool size: $FIRST_POOL_SIZE"
else
    # Different pool sizes, create temporary file
    POOL_SIZES_FILE="${OUTPUT_DIR}/.pool_sizes.tmp"
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would create pool sizes file: $POOL_SIZES_FILE"
        for i in "${!SAMPLE_NAMES[@]}"; do
            if [[ -n "${SAMPLE_NAMES[$i]}" ]]; then
                log_dry_run "  ${SAMPLE_NAMES[$i]},${SAMPLE_POOL_SIZES[$i]}"
            else
                # Derive sample name from BAM file
                # Remove .bam extension, then remove _All_seq.dedup pattern (from process_poolseq.sh)
                sample_name=$(basename "${SAMPLE_BAMS[$i]}" .bam)
                sample_name="${sample_name%_All_seq.dedup}"  # Remove _All_seq.dedup suffix if present
                sample_name="${sample_name%.dedup}"  # Remove .dedup if still present
                log_dry_run "  $sample_name,${SAMPLE_POOL_SIZES[$i]}"
            fi
        done
    else
        > "$POOL_SIZES_FILE"
        for i in "${!SAMPLE_NAMES[@]}"; do
            if [[ -n "${SAMPLE_NAMES[$i]}" ]]; then
                # Use sample name from CSV (after potential renaming)
                echo "${SAMPLE_NAMES[$i]},${SAMPLE_POOL_SIZES[$i]}" >> "$POOL_SIZES_FILE"
            else
                # Derive sample name from BAM file
                # Remove .bam extension, then remove _All_seq.dedup pattern (from process_poolseq.sh)
                sample_name=$(basename "${SAMPLE_BAMS[$i]}" .bam)
                sample_name="${sample_name%_All_seq.dedup}"  # Remove _All_seq.dedup suffix if present
                sample_name="${sample_name%.dedup}"  # Remove .dedup if still present
                echo "$sample_name,${SAMPLE_POOL_SIZES[$i]}" >> "$POOL_SIZES_FILE"
            fi
        done
        POOL_SIZES_ARG="$POOL_SIZES_FILE"
        log "Created pool sizes file with ${#SAMPLE_POOL_SIZES[@]} samples"
    fi
fi

# Determine window/step combinations to process
if [[ "$WINDOW_TYPE" == "interval" ]]; then
    # Windowed mode: process each window/step combination
    PIDS=()  # Array to track background job PIDs for parallel execution
    
    for scale_idx in "${!WINDOW_SIZE_ARRAY[@]}"; do
        window_size="${WINDOW_SIZE_ARRAY[$scale_idx]}"
        step_size="${STEP_SIZE_ARRAY[$scale_idx]}"
        
        # Always include window/step suffix for interval mode
        file_suffix="_w${window_size}_s${step_size}"

        # Abort early if expected outputs already exist
        if [[ "$DRY_RUN" == false ]] && [[ "$OVERWRITE_OUTPUT" == false ]]; then
            shopt -s nullglob
            # General warning for any FST outputs
            general_outputs=("$OUTPUT_DIR"/*fst*.csv "$OUTPUT_DIR"/*fst*.tsv)
            if [[ ${#general_outputs[@]} -gt 0 ]]; then
                log_warn "Found existing FST output(s) in $OUTPUT_DIR"
            fi

            # Abort only if this run's expected outputs already exist
            if [[ -n "$FILE_PREFIX" ]]; then
                expected_outputs=("$OUTPUT_DIR"/${FILE_PREFIX}*fst*${file_suffix}.csv "$OUTPUT_DIR"/${FILE_PREFIX}*fst*${file_suffix}.tsv)
            else
                expected_outputs=("$OUTPUT_DIR"/*fst*${file_suffix}.csv "$OUTPUT_DIR"/*fst*${file_suffix}.tsv)
            fi
            shopt -u nullglob
            if [[ ${#expected_outputs[@]} -gt 0 ]]; then
                log_warn "Found existing FST output(s) for window=$window_size, step=$step_size in $OUTPUT_DIR"
                log_warn "Expected outputs: ${expected_outputs[*]}"
                log_warn "Grenedalf will abort unless --overwrite-output is set"
                log_error "Aborting due to existing outputs (use --overwrite-output to overwrite)"
                exit 1
            fi
        fi
        
        log "Processing window=$window_size, step=$step_size"
        if [[ "$USE_MULTI_SCALE" == true ]]; then
            log "  Output suffix: $file_suffix"
        fi
        
        # Wait if we've reached max jobs (for parallel mode)
        if [[ "$PARALLEL_WINDOWS" == true ]]; then
            while [[ ${#PIDS[@]} -ge "$PARALLEL_MAX_JOBS" ]]; do
                # Check for completed jobs
                NEW_PIDS=()
                for pid in "${PIDS[@]}"; do
                    if kill -0 "$pid" 2>/dev/null; then
                        NEW_PIDS+=("$pid")
                    fi
                done
                PIDS=("${NEW_PIDS[@]}")
                sleep 1
            done
        fi
        
        # Build grenedalf fst command
        GRENEDALF_CMD=(
            "$GRENEDALF" fst
            --method "$FST_METHOD"
            --window-type interval
            --window-interval-width "$window_size"
            --window-interval-stride "$step_size"
            --window-average-policy "$WINDOW_AVERAGE_POLICY"
            --pool-sizes "$POOL_SIZES_ARG"
            --filter-sample-min-read-depth "$MIN_COVERAGE"
            --filter-sample-max-read-depth "$MAX_COVERAGE"
            --filter-sample-min-count "$MIN_COUNT"
            --separator-char "$OUTPUT_SEPARATOR"
            --out-dir "$OUTPUT_DIR"
            --threads "$THREADS"
        )

        # Allow overwriting outputs if requested
        if [[ "$OVERWRITE_OUTPUT" == true ]]; then
            GRENEDALF_CMD+=(--allow-file-overwriting)
        fi
        
        # Add file suffix for interval mode
        # Note: grenedalf automatically uses .tsv extension when --separator-char tab is used
        if [[ -n "$file_suffix" ]]; then
            GRENEDALF_CMD+=(--file-suffix "$file_suffix")
        fi
        
        # Add sample renaming if we have sample names from CSV
        if [[ -n "$RENAME_SAMPLES_FILE" ]] && [[ -f "$RENAME_SAMPLES_FILE" ]]; then
            GRENEDALF_CMD+=(--rename-samples-list "$RENAME_SAMPLES_FILE")
        fi
        
        # Add frequency or count filter
        if [[ "$MIN_COUNT_FILTER" -gt 0 ]]; then
            GRENEDALF_CMD+=(--filter-total-snp-min-count "$MIN_COUNT_FILTER")
        else
            GRENEDALF_CMD+=(--filter-total-snp-min-frequency "$MIN_FREQUENCY")
        fi
        
        # Add total read depth filter
        if [[ "$MIN_TOTAL_READ_DEPTH" -gt 0 ]]; then
            GRENEDALF_CMD+=(--filter-total-min-read-depth "$MIN_TOTAL_READ_DEPTH")
        fi
        
        # Add biallelic SNP filter
        if [[ "$FILTER_TOTAL_ONLY_BIALLELIC_SNPS" == true ]]; then
            GRENEDALF_CMD+=(--filter-total-only-biallelic-snps)
        fi
        
        # Add mask files if provided
        if [[ -n "$FILTER_MASK_TOTAL_FASTA" ]] && [[ -f "$FILTER_MASK_TOTAL_FASTA" ]]; then
            GRENEDALF_CMD+=(--filter-mask-total-fasta "$FILTER_MASK_TOTAL_FASTA")
        fi
        
        if [[ -n "$FILTER_MASK_TOTAL_BED" ]] && [[ -f "$FILTER_MASK_TOTAL_BED" ]]; then
            GRENEDALF_CMD+=(--filter-mask-total-bed "$FILTER_MASK_TOTAL_BED")
        fi
        
        # Add BAM files
        for bam in "${SAMPLE_BAMS[@]}"; do
            GRENEDALF_CMD+=(--sam-path "$bam")
        done
        
        # Add reference genome if provided
        if [[ -n "$REFERENCE_GENOME" ]]; then
            GRENEDALF_CMD+=(--reference-genome-fasta "$REFERENCE_GENOME")
        fi
        
        # Add file prefix if provided
        if [[ -n "$FILE_PREFIX" ]]; then
            GRENEDALF_CMD+=(--file-prefix "$FILE_PREFIX")
        fi
        
        # Add comparand options if provided
        if [[ -n "$COMPARAND" ]]; then
            GRENEDALF_CMD+=(--comparand "$COMPARAND")
        fi
        
        if [[ -n "$COMPARAND_LIST" ]]; then
            if [[ -n "$COMPARAND" ]]; then
                log_error "Cannot use both --comparand and --comparand-list"
                exit 1
            fi
            GRENEDALF_CMD+=(--comparand-list "$COMPARAND_LIST")
        fi
        
        # Run grenedalf fst (in parallel or sequential)
        if [[ "$PARALLEL_WINDOWS" == true ]] && [[ "$DRY_RUN" == false ]]; then
            log "  Launching grenedalf fst for window=$window_size, step=$step_size (background job)..."
            "${GRENEDALF_CMD[@]}" > "${OUTPUT_DIR}/.fst_w${window_size}_s${step_size}.log" 2>&1 &
            PID=$!
            PIDS+=("$PID")
            log "  Background job PID: $PID"
        else
            log "  Running grenedalf fst for window=$window_size, step=$step_size..."
            dry_run_cmd "${GRENEDALF_CMD[@]}" || {
                if [[ "$DRY_RUN" == false ]]; then
                    log_error "Grenedalf fst command failed for window=$window_size, step=$step_size"
                    # Clean up temp files if created
                    if [[ "$ALL_SAME_POOL_SIZE" == false ]] && [[ -f "$POOL_SIZES_FILE" ]]; then
                        rm -f "$POOL_SIZES_FILE"
                    fi
                    if [[ -n "$RENAME_SAMPLES_FILE" ]] && [[ -f "$RENAME_SAMPLES_FILE" ]]; then
                        rm -f "$RENAME_SAMPLES_FILE"
                    fi
                    exit 1
                fi
            }
            if [[ "$DRY_RUN" == false ]]; then
                normalize_output_extension
            fi
            log "  FST calculation complete for window=$window_size, step=$step_size"
        fi
    done
    
    # Wait for all background jobs to complete (if parallel mode)
    if [[ "$PARALLEL_WINDOWS" == true ]] && [[ "$DRY_RUN" == false ]]; then
        log "Waiting for all ${#PIDS[@]} background jobs to complete..."
        for pid in "${PIDS[@]}"; do
            wait "$pid" || {
                log_error "Background job $pid failed"
                # Clean up temp files if created
                if [[ "$ALL_SAME_POOL_SIZE" == false ]] && [[ -f "$POOL_SIZES_FILE" ]]; then
                    rm -f "$POOL_SIZES_FILE"
                fi
                if [[ -n "$RENAME_SAMPLES_FILE" ]] && [[ -f "$RENAME_SAMPLES_FILE" ]]; then
                    rm -f "$RENAME_SAMPLES_FILE"
                fi
                exit 1
            }
        done
        normalize_output_extension
        log "All window/step combinations completed"
        
        # Clean up log files if successful
        for scale_idx in "${!WINDOW_SIZE_ARRAY[@]}"; do
            window_size="${WINDOW_SIZE_ARRAY[$scale_idx]}"
            step_size="${STEP_SIZE_ARRAY[$scale_idx]}"
            log_file="${OUTPUT_DIR}/.fst_w${window_size}_s${step_size}.log"
            if [[ -f "$log_file" ]]; then
                rm -f "$log_file"
            fi
        done
    fi
else
    # Per-locus mode (single): single run
    # Build grenedalf fst command
    if [[ "$DRY_RUN" == false ]] && [[ "$OVERWRITE_OUTPUT" == false ]]; then
        shopt -s nullglob
        # General warning for any FST outputs
        general_outputs=("$OUTPUT_DIR"/*fst*.csv "$OUTPUT_DIR"/*fst*.tsv)
        if [[ ${#general_outputs[@]} -gt 0 ]]; then
            log_warn "Found existing FST output(s) in $OUTPUT_DIR"
        fi

        # Abort only if this run's expected outputs already exist
        if [[ -n "$FILE_PREFIX" ]]; then
            expected_outputs=("$OUTPUT_DIR"/${FILE_PREFIX}*fst*_single.csv "$OUTPUT_DIR"/${FILE_PREFIX}*fst*_single.tsv)
        else
            expected_outputs=("$OUTPUT_DIR"/*fst*_single.csv "$OUTPUT_DIR"/*fst*_single.tsv)
        fi
        shopt -u nullglob
        if [[ ${#expected_outputs[@]} -gt 0 ]]; then
            log_warn "Found existing per-locus FST output(s) in $OUTPUT_DIR"
            log_warn "Expected outputs: ${expected_outputs[*]}"
            log_warn "Grenedalf will abort unless --overwrite-output is set"
            log_error "Aborting due to existing outputs (use --overwrite-output to overwrite)"
            exit 1
        fi
    fi
    GRENEDALF_CMD=(
        "$GRENEDALF" fst
        --method "$FST_METHOD"
        --window-type single
        --window-average-policy "$WINDOW_AVERAGE_POLICY"
        --pool-sizes "$POOL_SIZES_ARG"
        --filter-sample-min-read-depth "$MIN_COVERAGE"
        --filter-sample-max-read-depth "$MAX_COVERAGE"
        --filter-sample-min-count "$MIN_COUNT"
        --separator-char "$OUTPUT_SEPARATOR"
        --out-dir "$OUTPUT_DIR"
        --threads "$THREADS"
    )

    # Allow overwriting outputs if requested
    if [[ "$OVERWRITE_OUTPUT" == true ]]; then
        GRENEDALF_CMD+=(--allow-file-overwriting)
    fi
    
    # Add frequency or count filter
    if [[ "$MIN_COUNT_FILTER" -gt 0 ]]; then
        GRENEDALF_CMD+=(--filter-total-snp-min-count "$MIN_COUNT_FILTER")
        log "Using min count filter: $MIN_COUNT_FILTER"
    else
        GRENEDALF_CMD+=(--filter-total-snp-min-frequency "$MIN_FREQUENCY")
        log "Using min frequency filter: $MIN_FREQUENCY"
    fi
    
    # Add total read depth filter
    if [[ "$MIN_TOTAL_READ_DEPTH" -gt 0 ]]; then
        GRENEDALF_CMD+=(--filter-total-min-read-depth "$MIN_TOTAL_READ_DEPTH")
        log "Using min total read depth filter: $MIN_TOTAL_READ_DEPTH"
    fi
    
    # Add biallelic SNP filter
    if [[ "$FILTER_TOTAL_ONLY_BIALLELIC_SNPS" == true ]]; then
        GRENEDALF_CMD+=(--filter-total-only-biallelic-snps)
        log "Filtering to biallelic SNPs only"
    fi
    
    # Add mask files if provided
    if [[ -n "$FILTER_MASK_TOTAL_FASTA" ]] && [[ -f "$FILTER_MASK_TOTAL_FASTA" ]]; then
        GRENEDALF_CMD+=(--filter-mask-total-fasta "$FILTER_MASK_TOTAL_FASTA")
        log "Using FASTA mask file: $FILTER_MASK_TOTAL_FASTA"
    fi
    
    if [[ -n "$FILTER_MASK_TOTAL_BED" ]] && [[ -f "$FILTER_MASK_TOTAL_BED" ]]; then
        GRENEDALF_CMD+=(--filter-mask-total-bed "$FILTER_MASK_TOTAL_BED")
        log "Using BED mask file: $FILTER_MASK_TOTAL_BED"
    fi
    
    # Add BAM files
    for bam in "${SAMPLE_BAMS[@]}"; do
        GRENEDALF_CMD+=(--sam-path "$bam")
    done
    
    # Add reference genome if provided
    if [[ -n "$REFERENCE_GENOME" ]]; then
        GRENEDALF_CMD+=(--reference-genome-fasta "$REFERENCE_GENOME")
    fi
    
    # Add file prefix if provided
    if [[ -n "$FILE_PREFIX" ]]; then
        GRENEDALF_CMD+=(--file-prefix "$FILE_PREFIX")
    fi
    
    # Add file suffix for single SNP mode to distinguish from windowed mode
    GRENEDALF_CMD+=(--file-suffix "_single")
    
    # Add sample renaming if we have sample names from CSV
    if [[ -n "$RENAME_SAMPLES_FILE" ]] && [[ -f "$RENAME_SAMPLES_FILE" ]]; then
        GRENEDALF_CMD+=(--rename-samples-list "$RENAME_SAMPLES_FILE")
    fi
    
    # Add comparand options if provided
    if [[ -n "$COMPARAND" ]]; then
        GRENEDALF_CMD+=(--comparand "$COMPARAND")
        log "Computing FST between $COMPARAND and all other samples"
    fi
    
    if [[ -n "$COMPARAND_LIST" ]]; then
        if [[ -n "$COMPARAND" ]]; then
            log_error "Cannot use both --comparand and --comparand-list"
            exit 1
        fi
        GRENEDALF_CMD+=(--comparand-list "$COMPARAND_LIST")
        log "Using comparand list from: $COMPARAND_LIST"
    fi
    
    # Run grenedalf fst
    log "Running grenedalf fst..."
    dry_run_cmd "${GRENEDALF_CMD[@]}" || {
        if [[ "$DRY_RUN" == false ]]; then
            log_error "Grenedalf fst command failed"
            # Clean up temp files if created
            if [[ "$ALL_SAME_POOL_SIZE" == false ]] && [[ -f "$POOL_SIZES_FILE" ]]; then
                rm -f "$POOL_SIZES_FILE"
            fi
            if [[ -n "$RENAME_SAMPLES_FILE" ]] && [[ -f "$RENAME_SAMPLES_FILE" ]]; then
                rm -f "$RENAME_SAMPLES_FILE"
            fi
            exit 1
        fi
    }
    if [[ "$DRY_RUN" == false ]]; then
        normalize_output_extension
    fi
fi

# Clean up temporary files if created
if [[ "$DRY_RUN" == false ]]; then
    if [[ "$ALL_SAME_POOL_SIZE" == false ]] && [[ -f "$POOL_SIZES_FILE" ]]; then
        rm -f "$POOL_SIZES_FILE"
        log "Cleaned up temporary pool sizes file"
    fi
    if [[ -n "$RENAME_SAMPLES_FILE" ]] && [[ -f "$RENAME_SAMPLES_FILE" ]]; then
        rm -f "$RENAME_SAMPLES_FILE"
        log "Cleaned up temporary rename samples file"
    fi
fi

if [[ "$DRY_RUN" == false ]]; then
    # Check for output file
    output_file=$(find "$OUTPUT_DIR" -name "*fst*.${OUTPUT_EXTENSION}" | head -1)
    if [[ -n "$output_file" ]] && [[ -f "$output_file" ]]; then
        log "FST calculation complete"
        log "  Output: $output_file"
        
        # Count sample pairs
        if [[ ${#SAMPLE_BAMS[@]} -eq 2 ]]; then
            log "  Computed FST for 1 sample pair"
        else
            num_pairs=$(( ${#SAMPLE_BAMS[@]} * (${#SAMPLE_BAMS[@]} - 1) / 2 ))
            log "  Computed FST for $num_pairs sample pairs"
        fi
    else
        log_warn "Could not find grenedalf output file in $OUTPUT_DIR"
        log_warn "Expected file with pattern *fst*.${OUTPUT_EXTENSION}"
    fi
fi

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run ""
    log_dry_run "DRY-RUN complete. All checks passed."
    log_dry_run "Run without --dry-run to execute the analysis."
else
    log "Analysis complete!"
    log "Output directory: $OUTPUT_DIR"
fi
