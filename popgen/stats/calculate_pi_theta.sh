#!/bin/bash

###############################################################################
# calculate_pi_theta.sh
# 
# Calculate nucleotide diversity (π), Watterson's theta (θ), and Tajima's D
# from poolseq data using grenedalf.
#
# These statistics measure genetic diversity within populations:
# - π (pi): Average number of pairwise nucleotide differences per site
# - θ (theta): Estimated from the number of segregating sites
# - Tajima's D: Test statistic comparing π and θ to detect selection/demography
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
WINDOW_SIZES=()  # Array for multi-scale analysis (can be specified multiple times)
STEP_SIZES=()    # Array for multi-scale analysis (optional, auto-calculated as window/2 if not provided)
MIN_COVERAGE=10
MAX_COVERAGE=1000
MIN_COUNT=2
POOL_SIZE=50
THREADS=1
DRY_RUN=false
WINDOW_AVERAGE_POLICY="valid-loci"
TAJIMA_D_POLICY="empirical-min-read-depth"
REFERENCE_GENOME=""
FILE_PREFIX=""  # Optional prefix for output files (default: no prefix)
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
  -b, --bam FILE              BAM file(s) - can be specified multiple times
  -i, --sample-info FILE      CSV file with sample information (see sample-info.csv)
                              Columns: sample_name, read1_file, read2_file, pool_size, bam_file
                              Note: read1_file and read2_file columns are ignored (for process_poolseq.sh)

Tool options:
  -g, --grenedalf PATH        Path to grenedalf executable (default: check PATH)

Required:
  -o, --output-dir DIR        Output directory for results

Optional options:
  -w, --window-size N         Window size in bp (can be specified multiple times for multi-scale analysis)
                               Default: 10000 (single scale). Multiple values enable multi-scale mode.
  --step-size N               Step size in bp (can be specified multiple times for multi-scale analysis)
                               Default: 5000 (single scale). If fewer step sizes than window sizes,
                               remaining step sizes are auto-calculated as window/2.
  --min-coverage N            Minimum coverage per site (default: 10)
  --max-coverage N            Maximum coverage per site (default: 1000)
  --min-count N               Minimum allele count (default: 2)
  --pool-size N               Pool size (haploid chromosomes) (default: 50)
  --window-average-policy     Window averaging policy: window-length, available-loci, valid-loci, valid-snps, sum (default: valid-loci)
  --tajima-d-policy           Tajima's D denominator policy: empirical-min-read-depth, provided-min-read-depth, popoolation-bugs, pool-size, uncorrected (default: empirical-min-read-depth)
  --reference-genome FILE     Reference genome FASTA file (optional, for better chromosome ordering)
  --file-prefix PREFIX        Optional prefix for output files (default: no prefix)
  --separator-char SEP        Output separator: comma, tab, space, semicolon (default: comma)
  --overwrite-output          Allow overwriting existing output files (default: false)
  -t, --threads N             Number of threads per grenedalf call (default: 1)
  --parallel-windows          Run multiple window/step combinations in parallel (default: false)
  --parallel-max-jobs N       Maximum concurrent window jobs when using --parallel-windows (default: number of CPU cores)
  --dry-run                   Preview commands without executing (dry-run mode)
  -h, --help                  Show this help message

Examples:
  # Single scale analysis (backward compatible):
  $0 \\
    --bam /path/to/sample.bam \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output \\
    --window-size 10000 \\
    --pool-size 50

  # Multi-scale analysis:
  $0 \\
    --sample-info sample_info.csv \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output \\
    --window-size 1000 \\
    --window-size 5000 \\
    --window-size 10000 \\
    --step-size 500 \\
    --step-size 2500 \\
    --step-size 5000

  # Multi-scale with auto-calculated step sizes (default: window/2):
  $0 \\
    --sample-info sample_info.csv \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output \\
    --window-size 1000 \\
    --window-size 5000 \\
    --window-size 10000

Output files:
  - One file per sample and per window/step: diversity_w{WINDOW}_s{STEP}_sample{SAMPLE}.{csv|tsv}
    (e.g., diversity_w10000_s5000_sampleCheney.tsv)
  - Optional prefix: {prefix}_diversity_w{WINDOW}_s{STEP}_sample{SAMPLE}.{csv|tsv}
  - Single-scale: one file per sample. Multi-scale: one file per sample per (window, step) combination.
  - Grenedalf writes .csv; when --separator-char tab, we rename to .tsv.
  - Logs: {output_dir}/log/{prefix_}calculate_pi_theta_YYYYmmdd_HHMMSS.log

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

# Normalize grenedalf output extension when using tab separator.
# Grenedalf writes .csv regardless of separator; we rename to .tsv for tab.
normalize_output_extension() {
    local target_dir="$1"
    if [[ "$OUTPUT_SEPARATOR" != "tab" ]]; then
        return 0
    fi
    shopt -s nullglob
    for csv_file in "$target_dir"/*diversity*.csv; do
        tsv_file="${csv_file%.csv}.tsv"
        mv -f "$csv_file" "$tsv_file"
    done
    shopt -u nullglob
}

# Add sample, window_size, and step_size columns to grenedalf diversity output file
# Args: output_file sample_name window_size step_size separator
add_sample_column() {
    local output_file="$1"
    local sample_name="$2"
    local window_size="$3"
    local step_size="$4"
    local separator="$5"
    
    if [[ ! -f "$output_file" ]]; then
        return 1
    fi
    
    # Determine separator character
    local sep_char=","
    if [[ "$separator" == "tab" ]]; then
        sep_char=$'\t'
    elif [[ "$separator" == "semicolon" ]]; then
        sep_char=";"
    elif [[ "$separator" == "space" ]]; then
        sep_char=" "
    fi
    
    # Create temporary file
    local temp_file="${output_file}.tmp"
    
    # Use awk to add sample, window_size, and step_size columns as the first columns
    awk -v sample="$sample_name" -v window="$window_size" -v step="$step_size" -v OFS="$sep_char" -v FS="$sep_char" '
        NR == 1 {
            # Header: add "sample", "window_size", "step_size" as first columns
            print "sample" OFS "window_size" OFS "step_size" OFS $0
        }
        NR > 1 {
            # Data rows: add sample name, window_size, step_size as first columns
            print sample OFS window OFS step OFS $0
        }
    ' "$output_file" > "$temp_file" && mv "$temp_file" "$output_file"
    
    return $?
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

# Function to parse CSV file (simple CSV parser, handles quoted fields)
parse_csv_line() {
    local line="$1"
    local IFS=','
    local -a fields
    local in_quotes=false
    local field=""
    
    # Simple CSV parsing - handles quoted fields with commas
    while IFS= read -r -d '' char || [[ -n "$char" ]]; do
        if [[ "$char" == '"' ]]; then
            in_quotes=$((1 - in_quotes))
        elif [[ "$char" == ',' ]] && [[ $in_quotes -eq 0 ]]; then
            fields+=("$field")
            field=""
        else
            field="${field}${char}"
        fi
    done < <(printf '%s\0' "$line")
    fields+=("$field")
    
    echo "${fields[@]}"
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
        -w|--window-size)
            WINDOW_SIZES+=("$2")
            shift 2
            ;;
        --step-size)
            STEP_SIZES+=("$2")
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
        --pool-size)
            POOL_SIZE="$2"
            shift 2
            ;;
        --window-average-policy)
            WINDOW_AVERAGE_POLICY="$2"
            shift 2
            ;;
        --tajima-d-policy)
            TAJIMA_D_POLICY="$2"
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

# Parse window sizes and step sizes for multi-scale analysis
declare -a WINDOW_SIZE_ARRAY
declare -a STEP_SIZE_ARRAY
USE_MULTI_SCALE=false

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
    
    log "Multi-scale analysis mode: ${#WINDOW_SIZE_ARRAY[@]} window/step combinations"
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
    log "Single-scale analysis mode: Window ${WINDOW_SIZE_ARRAY[0]} bp, Step ${STEP_SIZE_ARRAY[0]} bp"
else
    # No window sizes provided, use defaults
    WINDOW_SIZE_ARRAY=(10000)
    STEP_SIZE_ARRAY=(5000)
    log "Single-scale analysis mode: Window 10000 bp, Step 5000 bp (defaults)"
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

# Process sample info CSV if provided
# Initialize arrays (needed even when using --bam)
declare -a SAMPLE_NAMES
declare -a SAMPLE_BAMS
declare -a SAMPLE_POOL_SIZES
declare -a SAMPLE_MIN_COV
declare -a SAMPLE_MAX_COV

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
        SAMPLE_MIN_COV+=("$MIN_COVERAGE")
        SAMPLE_MAX_COV+=("$MAX_COVERAGE")
        
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
        SAMPLE_MIN_COV+=("$MIN_COVERAGE")
        SAMPLE_MAX_COV+=("$MAX_COVERAGE")
    done
else
    log_error "Must provide either BAM file(s) with --bam or sample info CSV with --sample-info"
    usage
    exit 1
fi

# Output separator
log "Output separator: $OUTPUT_SEPARATOR (.$OUTPUT_EXTENSION)"
if [[ "$OVERWRITE_OUTPUT" == true ]]; then
    log_warn "Overwriting enabled: existing outputs may be replaced"
else
    log_warn "Grenedalf aborts if output files already exist (use --overwrite-output to allow overwrite)"
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
        log_dry_run "  Sample $((i+1)): ${SAMPLE_BAMS[$i]}"
    done
    log_dry_run ""
fi

# Process each sample with each window/step combination
PIDS=()  # Array to track background job PIDs for parallel execution

for i in "${!SAMPLE_BAMS[@]}"; do
    bam_file="${SAMPLE_BAMS[$i]}"
    pool_size="${SAMPLE_POOL_SIZES[$i]}"
    min_cov="${SAMPLE_MIN_COV[$i]}"
    max_cov="${SAMPLE_MAX_COV[$i]}"
    
    # Determine sample name for output
    if [[ ${#SAMPLE_NAMES[@]} -gt $i ]] && [[ -n "${SAMPLE_NAMES[$i]}" ]]; then
        sample_name="${SAMPLE_NAMES[$i]}"
    else
        # Extract sample name from BAM file path
        # Remove .bam extension, then remove _All_seq.dedup pattern (from process_poolseq.sh)
        sample_name=$(basename "$bam_file" .bam)
        sample_name="${sample_name%_All_seq.dedup}"  # Remove _All_seq.dedup suffix if present
        sample_name="${sample_name%.dedup}"  # Remove .dedup if still present
    fi
    
    # Process each window/step combination for this sample
    for scale_idx in "${!WINDOW_SIZE_ARRAY[@]}"; do
        window_size="${WINDOW_SIZE_ARRAY[$scale_idx]}"
        step_size="${STEP_SIZE_ARRAY[$scale_idx]}"
        
        # Build file suffix: window, step, and sample (always include for clarity)
        # Format: _w{window}_s{step}_sample{sample}
        file_suffix="_w${window_size}_s${step_size}_sample${sample_name}"
        
        # Build expected output filename (grenedalf creates: {prefix}_diversity{suffix}.{ext})
        # Format: diversity_w{window}_s{step}_sample{sample}.{ext} or {prefix}_diversity_w{window}_s{step}_sample{sample}.{ext}
        if [[ -n "$FILE_PREFIX" ]]; then
            expected_filename="${FILE_PREFIX}_diversity${file_suffix}"
        else
            expected_filename="diversity${file_suffix}"
        fi

        # Abort early if expected outputs already exist
        if [[ "$DRY_RUN" == false ]] && [[ "$OVERWRITE_OUTPUT" == false ]]; then
            shopt -s nullglob
            # General warning for any diversity outputs
            general_outputs=("$OUTPUT_DIR"/*diversity*.csv "$OUTPUT_DIR"/*diversity*.tsv)
            if [[ ${#general_outputs[@]} -gt 0 ]]; then
                log_warn "Found existing diversity output(s) in $OUTPUT_DIR"
            fi

            # Abort only if this run's expected output files actually exist
            shopt -u nullglob
            if [[ -f "$OUTPUT_DIR/${expected_filename}.csv" ]] || [[ -f "$OUTPUT_DIR/${expected_filename}.tsv" ]]; then
                log_warn "Found existing diversity output(s) for sample=$sample_name, window=$window_size, step=$step_size in $OUTPUT_DIR"
                log_warn "Expected outputs: $OUTPUT_DIR/${expected_filename}.csv (or .tsv)"
                log_warn "Grenedalf will abort unless --overwrite-output is set"
                log_error "Aborting due to existing outputs (use --overwrite-output to overwrite)"
                exit 1
            fi
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
        
        log "Processing sample: $sample_name (window=$window_size, step=$step_size)"
        log "  BAM: $bam_file"
        log "  Pool size: $pool_size"
        log "  Window size: $window_size bp"
        log "  Step size: $step_size bp"
        log "  Coverage range: $min_cov - $max_cov"
        log "  Output file: ${expected_filename}.${OUTPUT_EXTENSION}"
        
        # Build grenedalf command
        GRENEDALF_CMD=(
            "$GRENEDALF" diversity
            --sam-path "$bam_file"
            --window-type interval
            --window-interval-width "$window_size"
            --window-interval-stride "$step_size"
            --window-average-policy "$WINDOW_AVERAGE_POLICY"
            --pool-sizes "$pool_size"
            --filter-sample-min-read-depth "$min_cov"
            --filter-sample-max-read-depth "$max_cov"
            --filter-sample-min-count "$MIN_COUNT"
            --tajima-d-denominator-policy "$TAJIMA_D_POLICY"
            --separator-char "$OUTPUT_SEPARATOR"
            --out-dir "$OUTPUT_DIR"
            --threads "$THREADS"
        )

        # Allow overwriting outputs if requested
        if [[ "$OVERWRITE_OUTPUT" == true ]]; then
            GRENEDALF_CMD+=(--allow-file-overwriting)
        fi
        
        # Add file suffix for window/step/sample (always include)
        # Grenedalf writes .csv; we rename to .tsv when --separator-char tab is used
        GRENEDALF_CMD+=(--file-suffix "$file_suffix")
        
        # Add reference genome if provided
        if [[ -n "$REFERENCE_GENOME" ]]; then
            GRENEDALF_CMD+=(--reference-genome-fasta "$REFERENCE_GENOME")
        fi
        
        # Add file prefix if provided (goes before "diversity" base name)
        if [[ -n "$FILE_PREFIX" ]]; then
            GRENEDALF_CMD+=(--file-prefix "$FILE_PREFIX")
        fi
        
        # Run grenedalf (in parallel or sequential)
        if [[ "$PARALLEL_WINDOWS" == true ]] && [[ "$DRY_RUN" == false ]]; then
            log "  Launching grenedalf diversity for window=$window_size, step=$step_size (background job)..."
            "${GRENEDALF_CMD[@]}" > "${OUTPUT_DIR}/.diversity_${sample_name}_w${window_size}_s${step_size}.log" 2>&1 &
            PID=$!
            PIDS+=("$PID")
            log "  Background job PID: $PID"
        else
            log "  Running grenedalf diversity for window=$window_size, step=$step_size..."
            dry_run_cmd "${GRENEDALF_CMD[@]}" || {
                if [[ "$DRY_RUN" == false ]]; then
                    log_error "Grenedalf diversity command failed for sample: $sample_name (window=$window_size, step=$step_size)"
                    exit 1
                fi
            }
            
            if [[ "$DRY_RUN" == false ]]; then
                # Normalize output extension (.csv -> .tsv if using tab separator)
                normalize_output_extension "$OUTPUT_DIR"
                
                # Find output file (grenedalf creates: {prefix}_diversity{suffix}.{ext})
                output_file="${OUTPUT_DIR}/${expected_filename}.${OUTPUT_EXTENSION}"
                
                if [[ ! -f "$output_file" ]]; then
                    # Try alternative pattern (grenedalf might use different extension or naming)
                    output_pattern="${expected_filename}*.${OUTPUT_EXTENSION}"
                    output_file=$(find "$OUTPUT_DIR" -maxdepth 1 -name "$output_pattern" | head -1)
                fi
                
                if [[ -n "$output_file" ]] && [[ -f "$output_file" ]]; then
                    log "  Diversity calculation complete for window=$window_size, step=$step_size"
                    log "    Output: $output_file"
                    
                    # Add sample, window_size, and step_size columns to output file if not already present
                    if ! head -1 "$output_file" | grep -q "^sample"; then
                        log "  Adding sample, window_size, and step_size columns to output..."
                        add_sample_column "$output_file" "$sample_name" "$window_size" "$step_size" "$OUTPUT_SEPARATOR" || {
                            log_warn "Failed to add metadata columns to $output_file"
                        }
                    else
                        log "  Metadata columns already present in output"
                    fi
                else
                    log_warn "Could not find grenedalf output file in $OUTPUT_DIR"
                    log_warn "Expected file with pattern $output_pattern"
                fi
            fi
        fi
    done
done

# Wait for all background jobs to complete (if parallel mode)
if [[ "$PARALLEL_WINDOWS" == true ]] && [[ "$DRY_RUN" == false ]]; then
    log "Waiting for all ${#PIDS[@]} background jobs to complete..."
    for pid in "${PIDS[@]}"; do
        wait "$pid" || {
            log_error "Background job $pid failed"
            exit 1
        }
    done
    log "All window/step combinations completed"
    
    # Normalize output extensions and add sample columns
    normalize_output_extension "$OUTPUT_DIR"
    
    # Clean up log files and verify outputs
    for i in "${!SAMPLE_BAMS[@]}"; do
        if [[ ${#SAMPLE_NAMES[@]} -gt $i ]] && [[ -n "${SAMPLE_NAMES[$i]}" ]]; then
            sample_name="${SAMPLE_NAMES[$i]}"
        else
            sample_name=$(basename "${SAMPLE_BAMS[$i]}" .bam)
            sample_name="${sample_name%_All_seq.dedup}"
            sample_name="${sample_name%.dedup}"
        fi
        
        for scale_idx in "${!WINDOW_SIZE_ARRAY[@]}"; do
            window_size="${WINDOW_SIZE_ARRAY[$scale_idx]}"
            step_size="${STEP_SIZE_ARRAY[$scale_idx]}"
            
            # Build expected filename to match what was used during execution
            file_suffix="_w${window_size}_s${step_size}_sample${sample_name}"
            if [[ -n "$FILE_PREFIX" ]]; then
                expected_filename="${FILE_PREFIX}_diversity${file_suffix}"
            else
                expected_filename="diversity${file_suffix}"
            fi
            
            log_file="${OUTPUT_DIR}/.diversity_${sample_name}_w${window_size}_s${step_size}.log"
            if [[ -f "$log_file" ]]; then
                rm -f "$log_file"
            fi
            
            # Find output file (grenedalf creates: {prefix}_diversity{suffix}.{ext})
            output_file="${OUTPUT_DIR}/${expected_filename}.${OUTPUT_EXTENSION}"
            
            if [[ ! -f "$output_file" ]]; then
                # Try alternative pattern (grenedalf might use different extension or naming)
                output_pattern="${expected_filename}*.${OUTPUT_EXTENSION}"
                output_file=$(find "$OUTPUT_DIR" -maxdepth 1 -name "$output_pattern" | head -1)
            fi
            
            if [[ -n "$output_file" ]] && [[ -f "$output_file" ]]; then
                log "  Verified output for $sample_name (window=$window_size, step=$step_size): $output_file"
                
                # Add sample, window_size, and step_size columns if not already present
                if ! head -1 "$output_file" | grep -q "^sample"; then
                    log "  Adding sample, window_size, and step_size columns to output..."
                    add_sample_column "$output_file" "$sample_name" "$window_size" "$step_size" "$OUTPUT_SEPARATOR" || {
                        log_warn "Failed to add metadata columns to $output_file"
                    }
                fi
            else
                log_warn "Could not find output file for $sample_name (window=$window_size, step=$step_size)"
            fi
        done
    done
fi

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run ""
    log_dry_run "DRY-RUN complete. All checks passed."
    log_dry_run "Run without --dry-run to execute the analysis."
else
    log "Analysis complete!"
    log "Output directory: $OUTPUT_DIR"
    shopt -s nullglob
    output_files=("$OUTPUT_DIR"/*diversity*.${OUTPUT_EXTENSION})
    shopt -u nullglob
    if [[ ${#output_files[@]} -eq 0 ]]; then
        log "  (no diversity output files found)"
    else
        if [[ "$USE_MULTI_SCALE" == true ]]; then
            log "Output files (${#output_files[@]} total; multi-scale: one per sample per window/step):"
            for scale_idx in "${!WINDOW_SIZE_ARRAY[@]}"; do
                w="${WINDOW_SIZE_ARRAY[$scale_idx]}"
                s="${STEP_SIZE_ARRAY[$scale_idx]}"
                scale_files=()
                for f in "${output_files[@]}"; do
                    [[ -f "$f" ]] && [[ "$(basename "$f")" == *"diversity_w${w}_s${s}_"* ]] && scale_files+=("$f")
                done
                if [[ ${#scale_files[@]} -gt 0 ]]; then
                    log "  w=${w} s=${s}:"
                    for f in "${scale_files[@]}"; do
                        log "    $(basename "$f")"
                    done
                fi
            done
        else
            log "Output files:"
            for output_file in "${output_files[@]}"; do
                log "  $(basename "$output_file")"
            done
        fi
    fi
fi
