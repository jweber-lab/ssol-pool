#!/bin/bash

###############################################################################
# plot_fst_cathedral.sh
# 
# Generate FST cathedral plots using grenedalf fst-cathedral and cathedral-plot.
#
# Cathedral plots visualize FST across different window sizes, showing both
# broad and fine-scale differentiation patterns between populations.
#
# Workflow:
# 1. Calls grenedalf fst-cathedral to compute per-pixel FST matrices (CSV/JSON)
# 2. Calls grenedalf cathedral-plot to generate BMP and SVG plots
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
FST_METHOD="unbiased-nei"
CATHEDRAL_WIDTH=1500
CATHEDRAL_HEIGHT=500
MIN_FREQUENCY=0.01
MIN_COUNT_FILTER=0
MIN_COVERAGE=10
MAX_COVERAGE=1000
MIN_COUNT=2
POOL_SIZE=50
THREADS=1
DRY_RUN=false
REFERENCE_GENOME=""
FILE_PREFIX=""  # Optional prefix for output files (default: no prefix)
COMPARAND=""
COMPARAND_LIST=""
COLOR_LIST="inferno"
MIN_VALUE=""
MAX_VALUE=""
CLIP_UNDER=false
CLIP_OVER=false
COLOR_NORMALIZATION="linear"
VERBOSE=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

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
  --method METHOD             FST estimator method (default: unbiased-nei)
                              Options: unbiased-nei, unbiased-hudson
  --min-frequency FLOAT       Minimum allele frequency cutoff (default: 0.01 for 1%)
  --min-count-filter INT      Minimum allele count cutoff, alternative to frequency (default: 0 = disabled)

Filtering options:
  --min-coverage N            Minimum coverage per site per sample (default: 10)
  --max-coverage N            Maximum coverage per site per sample (default: 1000)
  --min-count N               Minimum allele count per sample (default: 2)

Comparison options:
  --comparand SAMPLE          Compute FST between this sample and all others (optional)
  --comparand-list FILE       File with sample pairs to compute FST for (one pair per line, optional)

Plotting options:
  --cathedral-width N         Plot width in pixels (default: 1500)
  --cathedral-height N        Plot height in pixels (default: 500)
  --color-list LIST           Color palette for plots (default: inferno)
                              Can be color name, file with colors, or comma-separated list
  --min-value FLOAT           Minimum value for color scale (default: auto)
  --max-value FLOAT           Maximum value for color scale (default: auto)
  --clip-under                Clip values below min to min value
  --clip-over                 Clip values above max to max value
  --color-normalization TYPE  Color normalization: linear or logarithmic (default: linear)

Other options:
  --reference-genome FILE     Reference genome FASTA file (recommended for chromosome ordering)
  --file-prefix PREFIX        Optional prefix for output files (default: no prefix)
  --pool-size N               Default pool size if not in CSV (default: 50)
  -t, --threads N             Number of threads (default: 1)
  --verbose                   Enable verbose/debug output
  --dry-run                   Preview commands without executing (dry-run mode)
  -h, --help                  Show this help message

Examples:
  # Using sample info CSV file (all pairwise comparisons):
  $0 \\
    --sample-info sample-info.csv \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output \\
    --reference-genome /path/to/reference.fa

  # Using BAM files directly:
  $0 \\
    --bam /path/to/sample1.bam \\
    --bam /path/to/sample2.bam \\
    --grenedalf /usr/bin/grenedalf \\
    --output-dir /path/to/output \\
    --reference-genome /path/to/reference.fa

  # Compute FST between one sample and all others:
  $0 \\
    --sample-info sample-info.csv \\
    --comparand Echo_Pool_S1 \\
    --output-dir /path/to/output \\
    --reference-genome /path/to/reference.fa \\
    --cathedral-width 2000 \\
    --cathedral-height 600

Output files:
  - {output_dir}/*.csv: Per-pixel FST value matrices (from fst-cathedral)
  - {output_dir}/*.json: Plot metadata (from fst-cathedral)
  - {output_dir}/*.bmp: Bitmap plots (from cathedral-plot)
  - {output_dir}/*.svg: Vector plots with axes and legend (from cathedral-plot)
  - Files are named per chromosome and sample pair

Note: FST requires at least 2 samples. Cathedral plots are generated per chromosome and sample pair.

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

log_debug() {
    if [[ "$VERBOSE" == true ]]; then
        echo -e "${YELLOW}[DEBUG]${NC} $1" >&2
    fi
}

# Function to execute command or print in dry-run mode
dry_run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: $*"
        return 0
    else
        if [[ "$VERBOSE" == true ]]; then
            log_debug "Executing: $*"
        fi
        "$@"
    fi
}

# Function to diagnose NaN values in JSON files
# Reports on NaN frequency, location, and context to help understand their cause
diagnose_json_nan() {
    local json_file="$1"
    
    if [[ ! -f "$json_file" ]]; then
        log_error "JSON file not found: $json_file"
        return 1
    fi
    
    # Check for various NaN patterns
    local nan_count=0
    local nan_patterns=()
    
    # Count different NaN patterns
    if grep -q '"nan"' "$json_file" 2>/dev/null; then
        local count=$(grep -o '"nan"' "$json_file" | wc -l)
        nan_count=$((nan_count + count))
        nan_patterns+=("string \"nan\": $count")
    fi
    
    if grep -q ': nan' "$json_file" 2>/dev/null; then
        local count=$(grep -o ': nan' "$json_file" | wc -l)
        nan_count=$((nan_count + count))
        nan_patterns+=("unquoted nan: $count")
    fi
    
    if grep -q ', nan' "$json_file" 2>/dev/null; then
        local count=$(grep -o ', nan' "$json_file" | wc -l)
        nan_count=$((nan_count + count))
        nan_patterns+=("comma-separated nan: $count")
    fi
    
    if [[ $nan_count -gt 0 ]]; then
        log_warn "Found $nan_count NaN value(s) in JSON file: $(basename "$json_file")"
        
        if [[ "$VERBOSE" == true ]]; then
            log_debug "NaN patterns found: ${nan_patterns[*]}"
            
            # Show context around NaN values (first 5 occurrences)
            log_debug "Context around NaN values:"
            grep -n -B 2 -A 2 -i 'nan' "$json_file" | head -n 30 | while IFS= read -r line; do
                log_debug "  $line"
            done
            
            # Check JSON structure to understand where NaN appears
            log_debug "Checking JSON structure..."
            if command -v python3 &> /dev/null; then
                # Try to parse JSON and see where it fails
                if ! python3 -m json.tool "$json_file" > /dev/null 2>&1; then
                    log_debug "JSON file is invalid (likely due to NaN values)"
                    # Try to find the line number where parsing fails
                    python3 -m json.tool "$json_file" 2>&1 | head -n 5 | while IFS= read -r err_line; do
                        log_debug "  JSON error: $err_line"
                    done
                else
                    log_debug "JSON file structure is valid (NaN may be in data values)"
                fi
            fi
            
            # Check file size and line count
            local file_size=$(wc -c < "$json_file")
            local line_count=$(wc -l < "$json_file")
            log_debug "JSON file stats: ${file_size} bytes, ${line_count} lines"
            
            # Show first few lines to understand structure
            log_debug "First 10 lines of JSON file:"
            head -n 10 "$json_file" | while IFS= read -r line; do
                log_debug "  $line"
            done
        fi
        
        # Explain possible causes
        log_warn "NaN values in FST cathedral plots typically indicate:"
        log_warn "  1. Insufficient data (no variants or coverage at those positions)"
        log_warn "  2. Division by zero in FST calculation (identical allele frequencies)"
        log_warn "  3. Invalid allele frequencies (outside 0-1 range)"
        log_warn "  4. Missing data that couldn't be handled"
        
        return 1
    else
        log_debug "No NaN values found in JSON file"
        return 0
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
        --comparand)
            COMPARAND="$2"
            shift 2
            ;;
        --comparand-list)
            COMPARAND_LIST="$2"
            shift 2
            ;;
        --cathedral-width)
            CATHEDRAL_WIDTH="$2"
            shift 2
            ;;
        --cathedral-height)
            CATHEDRAL_HEIGHT="$2"
            shift 2
            ;;
        --color-list)
            COLOR_LIST="$2"
            shift 2
            ;;
        --min-value)
            MIN_VALUE="$2"
            shift 2
            ;;
        --max-value)
            MAX_VALUE="$2"
            shift 2
            ;;
        --clip-under)
            CLIP_UNDER=true
            shift
            ;;
        --clip-over)
            CLIP_OVER=true
            shift
            ;;
        --color-normalization)
            COLOR_NORMALIZATION="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --reference-genome)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        --file-prefix)
            FILE_PREFIX="$2"
            shift 2
            ;;
        --pool-size)
            POOL_SIZE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
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
    unbiased-nei|unbiased-hudson)
        ;;
    *)
        log_error "Invalid FST method: $FST_METHOD"
        log_error "Valid methods: unbiased-nei, unbiased-hudson"
        exit 1
        ;;
esac

# Validate color normalization
case "$COLOR_NORMALIZATION" in
    linear|logarithmic)
        ;;
    *)
        log_error "Invalid color normalization: $COLOR_NORMALIZATION"
        log_error "Valid options: linear, logarithmic"
        exit 1
        ;;
esac

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
    
    # Read CSV file (skip header line)
    line_num=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        line_num=$((line_num + 1))
        
        # Skip header line
        if [[ $line_num -eq 1 ]]; then
            continue
        fi
        
        # Skip empty lines and comment lines
        [[ -z "$line" ]] && continue
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        
        # Parse CSV line (simple parsing, assumes no quoted commas in values)
        IFS=',' read -r -a fields <<< "$line"
        
        # Trim whitespace from fields
        for i in "${!fields[@]}"; do
            fields[$i]=$(echo "${fields[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        done
        
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

# Create output directory
check_dir "$OUTPUT_DIR"

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
    log_dry_run "Cathedral plot dimensions: ${CATHEDRAL_WIDTH}x${CATHEDRAL_HEIGHT} pixels"
    log_dry_run ""
fi

log "Computing FST cathedral plots for ${#SAMPLE_BAMS[@]} samples"
log "FST method: $FST_METHOD"
log "Cathedral plot dimensions: ${CATHEDRAL_WIDTH}x${CATHEDRAL_HEIGHT} pixels"

# Step 1: Run grenedalf fst-cathedral
log "Step 1/2: Computing FST cathedral plot data with grenedalf fst-cathedral..."

if [[ "$VERBOSE" == true ]]; then
    log_debug "FST cathedral computation parameters:"
    log_debug "  Method: $FST_METHOD"
    log_debug "  Pool sizes: $POOL_SIZES_ARG"
    log_debug "  Min coverage: $MIN_COVERAGE"
    log_debug "  Max coverage: $MAX_COVERAGE"
    log_debug "  Min count: $MIN_COUNT"
    if [[ "$MIN_COUNT_FILTER" -gt 0 ]]; then
        log_debug "  Min count filter: $MIN_COUNT_FILTER"
    else
        log_debug "  Min frequency filter: $MIN_FREQUENCY"
    fi
    log_debug "  Cathedral dimensions: ${CATHEDRAL_WIDTH}x${CATHEDRAL_HEIGHT}"
    log_debug "  Threads: $THREADS"
    log_debug "  Output directory: $OUTPUT_DIR"
    if [[ -n "$REFERENCE_GENOME" ]]; then
        log_debug "  Reference genome: $REFERENCE_GENOME"
    fi
    log_debug "  Number of samples: ${#SAMPLE_BAMS[@]}"
fi

# Build grenedalf fst-cathedral command
GRENEDALF_CATHEDRAL_CMD=(
    "$GRENEDALF" fst-cathedral
    --method "$FST_METHOD"
    --pool-sizes "$POOL_SIZES_ARG"
    --filter-sample-min-read-depth "$MIN_COVERAGE"
    --filter-sample-max-read-depth "$MAX_COVERAGE"
    --filter-sample-min-count "$MIN_COUNT"
    --cathedral-width "$CATHEDRAL_WIDTH"
    --cathedral-height "$CATHEDRAL_HEIGHT"
    --out-dir "$OUTPUT_DIR"
    --threads "$THREADS"
    --allow-file-overwriting
)

# Add frequency or count filter
if [[ "$MIN_COUNT_FILTER" -gt 0 ]]; then
    GRENEDALF_CATHEDRAL_CMD+=(--filter-total-snp-min-count "$MIN_COUNT_FILTER")
else
    GRENEDALF_CATHEDRAL_CMD+=(--filter-total-snp-min-frequency "$MIN_FREQUENCY")
fi

# Add BAM files
for bam in "${SAMPLE_BAMS[@]}"; do
    GRENEDALF_CATHEDRAL_CMD+=(--sam-path "$bam")
done

# Add reference genome if provided
if [[ -n "$REFERENCE_GENOME" ]]; then
    GRENEDALF_CATHEDRAL_CMD+=(--reference-genome-fasta "$REFERENCE_GENOME")
fi

# Add file prefix if provided
if [[ -n "$FILE_PREFIX" ]]; then
    GRENEDALF_CATHEDRAL_CMD+=(--file-prefix "$FILE_PREFIX")
fi

# Add comparand options if provided
if [[ -n "$COMPARAND" ]]; then
    GRENEDALF_CATHEDRAL_CMD+=(--comparand "$COMPARAND")
    log "Computing FST between $COMPARAND and all other samples"
fi

if [[ -n "$COMPARAND_LIST" ]]; then
    if [[ -n "$COMPARAND" ]]; then
        log_error "Cannot use both --comparand and --comparand-list"
        exit 1
    fi
    GRENEDALF_CATHEDRAL_CMD+=(--comparand-list "$COMPARAND_LIST")
    log "Using comparand list from: $COMPARAND_LIST"
fi

# Run grenedalf fst-cathedral
log "Running grenedalf fst-cathedral..."

if [[ "$DRY_RUN" == false ]]; then
    if [[ "$VERBOSE" == true ]]; then
        log_debug "Full command: ${GRENEDALF_CATHEDRAL_CMD[*]}"
        log_debug "Starting fst-cathedral computation..."
    fi
    
    if ! "${GRENEDALF_CATHEDRAL_CMD[@]}" 2>&1 | tee /tmp/grenedalf_fst_cathedral_stderr.log; then
        log_error "Grenedalf fst-cathedral command failed"
        if [[ -f /tmp/grenedalf_fst_cathedral_stderr.log ]]; then
            log_error "Error output:"
            cat /tmp/grenedalf_fst_cathedral_stderr.log >&2
            rm -f /tmp/grenedalf_fst_cathedral_stderr.log
        fi
        # Clean up temp file if created
        if [[ "$ALL_SAME_POOL_SIZE" == false ]] && [[ -f "$POOL_SIZES_FILE" ]]; then
            rm -f "$POOL_SIZES_FILE"
        fi
        exit 1
    fi
    rm -f /tmp/grenedalf_fst_cathedral_stderr.log
    
    if [[ "$VERBOSE" == true ]]; then
        log_debug "fst-cathedral computation complete"
        # List generated files
        local csv_count=$(find "$OUTPUT_DIR" -name "*.csv" -type f | wc -l)
        local json_count=$(find "$OUTPUT_DIR" -name "*.json" -type f | wc -l)
        log_debug "Generated files: $csv_count CSV, $json_count JSON"
    fi
else
    dry_run_cmd "${GRENEDALF_CATHEDRAL_CMD[@]}"
fi

log "FST cathedral data computation complete"

# Step 2: Generate plots with grenedalf cathedral-plot
log "Step 2/2: Generating cathedral plots with grenedalf cathedral-plot..."

# Find all JSON files created by fst-cathedral
if [[ "$DRY_RUN" == false ]]; then
    JSON_FILES=($(find "$OUTPUT_DIR" -name "*.json" -type f | sort))
    
    if [[ ${#JSON_FILES[@]} -eq 0 ]]; then
        log_warn "No JSON files found in $OUTPUT_DIR"
        log_warn "Expected JSON files from grenedalf fst-cathedral"
    else
        log "Found ${#JSON_FILES[@]} JSON file(s) to process"
        
        for json_file in "${JSON_FILES[@]}"; do
            log "Processing: $(basename "$json_file")"
            
            # Diagnose NaN values in JSON file
            if ! diagnose_json_nan "$json_file"; then
                log_error "JSON file contains NaN values: $json_file"
                log_error "This will cause grenedalf cathedral-plot to fail"
                log_error "NaN values typically indicate data quality issues:"
                log_error "  - Insufficient coverage/variants at some positions"
                log_error "  - Identical allele frequencies (division by zero)"
                log_error "  - Missing or invalid data"
                log_error ""
                log_error "Consider:"
                log_error "  - Checking coverage filters (--min-coverage, --max-coverage)"
                log_error "  - Adjusting frequency filters (--min-frequency)"
                log_error "  - Reviewing input BAM files for data quality"
                log_error "  - Checking if specific chromosomes/regions have issues"
                log_error ""
                log_error "Skipping this file to prevent crash"
                continue
            fi
            
            # Show JSON file info in verbose mode
            if [[ "$VERBOSE" == true ]]; then
                log_debug "JSON file size: $(wc -c < "$json_file") bytes"
                log_debug "JSON file line count: $(wc -l < "$json_file") lines"
                
                # Validate JSON structure
                if command -v python3 &> /dev/null; then
                    if python3 -m json.tool "$json_file" > /dev/null 2>&1; then
                        log_debug "JSON structure is valid"
                    else
                        log_warn "JSON structure validation failed:"
                        python3 -m json.tool "$json_file" 2>&1 | head -n 5 | while IFS= read -r err_line; do
                            log_debug "  $err_line"
                        done
                    fi
                fi
            fi
            
            # Build grenedalf cathedral-plot command
            GRENEDALF_PLOT_CMD=(
                "$GRENEDALF" cathedral-plot
                --json-path "$json_file"
                --color-list "$COLOR_LIST"
                --color-normalization "$COLOR_NORMALIZATION"
                --out-dir "$OUTPUT_DIR"
                --allow-file-overwriting
            )
            
            # Add min/max value options if provided
            if [[ -n "$MIN_VALUE" ]]; then
                GRENEDALF_PLOT_CMD+=(--min-value "$MIN_VALUE")
            fi
            
            if [[ -n "$MAX_VALUE" ]]; then
                GRENEDALF_PLOT_CMD+=(--max-value "$MAX_VALUE")
            fi
            
            # Add clipping options
            if [[ "$CLIP_UNDER" == true ]]; then
                GRENEDALF_PLOT_CMD+=(--clip-under)
            fi
            
            if [[ "$CLIP_OVER" == true ]]; then
                GRENEDALF_PLOT_CMD+=(--clip-over)
            fi
            
            # Run grenedalf cathedral-plot
            if [[ "$DRY_RUN" == false ]]; then
                if [[ "$VERBOSE" == true ]]; then
                    log_debug "Running command: ${GRENEDALF_PLOT_CMD[*]}"
                fi
                
                # Capture stderr for better error reporting
                if ! "${GRENEDALF_PLOT_CMD[@]}" 2>&1 | tee /tmp/grenedalf_cathedral_plot_stderr.log; then
                    log_error "Grenedalf cathedral-plot command failed for: $json_file"
                    if [[ -f /tmp/grenedalf_cathedral_plot_stderr.log ]]; then
                        log_error "Error output:"
                        cat /tmp/grenedalf_cathedral_plot_stderr.log >&2
                        rm -f /tmp/grenedalf_cathedral_plot_stderr.log
                    fi
                    # Clean up temp file if created
                    if [[ "$ALL_SAME_POOL_SIZE" == false ]] && [[ -f "$POOL_SIZES_FILE" ]]; then
                        rm -f "$POOL_SIZES_FILE"
                    fi
                    exit 1
                fi
                rm -f /tmp/grenedalf_cathedral_plot_stderr.log
            else
                dry_run_cmd "${GRENEDALF_PLOT_CMD[@]}"
            fi
        done
    fi
else
    log_dry_run "Would find JSON files in $OUTPUT_DIR and process them with grenedalf cathedral-plot"
    log_dry_run "  Color list: $COLOR_LIST"
    log_dry_run "  Color normalization: $COLOR_NORMALIZATION"
    if [[ -n "$MIN_VALUE" ]]; then
        log_dry_run "  Min value: $MIN_VALUE"
    fi
    if [[ -n "$MAX_VALUE" ]]; then
        log_dry_run "  Max value: $MAX_VALUE"
    fi
    if [[ "$CLIP_UNDER" == true ]]; then
        log_dry_run "  Clip under: enabled"
    fi
    if [[ "$CLIP_OVER" == true ]]; then
        log_dry_run "  Clip over: enabled"
    fi
fi

# Clean up temporary pool sizes file if created
if [[ "$DRY_RUN" == false ]] && [[ "$ALL_SAME_POOL_SIZE" == false ]] && [[ -f "$POOL_SIZES_FILE" ]]; then
    rm -f "$POOL_SIZES_FILE"
    log "Cleaned up temporary pool sizes file"
fi

if [[ "$DRY_RUN" == false ]]; then
    # Check for output files
    CSV_COUNT=$(find "$OUTPUT_DIR" -name "*.csv" -type f | wc -l)
    BMP_COUNT=$(find "$OUTPUT_DIR" -name "*.bmp" -type f | wc -l)
    SVG_COUNT=$(find "$OUTPUT_DIR" -name "*.svg" -type f | wc -l)
    
    log "Cathedral plot generation complete"
    log "  CSV files (matrices): $CSV_COUNT"
    log "  JSON files (metadata): ${#JSON_FILES[@]}"
    log "  BMP files (plots): $BMP_COUNT"
    log "  SVG files (plots): $SVG_COUNT"
    
    if [[ $BMP_COUNT -eq 0 ]] && [[ $SVG_COUNT -eq 0 ]]; then
        log_warn "No plot files found in $OUTPUT_DIR"
        log_warn "Expected BMP and/or SVG files from grenedalf cathedral-plot"
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
