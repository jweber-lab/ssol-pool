#!/bin/bash

###############################################################################
# identify_outliers.sh
#
# Bash wrapper for identify_outliers.R. Identifies outlier windows from
# diversity (pi, theta, tajima_d), FST, and PBE using quantile thresholds.
# Main options: --hdf5-dir or --diversity-dir/--fst-dir, --output-dir,
# --high-quantile, --low-quantile, --statistics (pi, theta, tajima_d, pbe, fst).
# For output file formats and column naming, see the header of identify_outliers.R.
# Usage: See below or --help
###############################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Default parameters
FST_DIR=""
DIVERSITY_DIR=""
HDF5_DIR=""
OUTPUT_DIR=""
MERGE_DISTANCE="auto"
MIN_DEPTH=""
MAX_DEPTH=""
MIN_MAPPING_QUALITY=""
MAX_MAPPING_QUALITY=""
MIN_SNPS=""
MAX_SNPS=""
MIN_REGION_SNPS=""
MAX_REGION_SNPS=""
REFERENCE_GENOME=""
HIGH_QUANTILE=""
LOW_QUANTILE=""
SEED_HIGH_QUANTILE=""
SEED_LOW_QUANTILE=""
EXPAND_HIGH_QUANTILE=""
EXPAND_LOW_QUANTILE=""
TOP_N_CHROMOSOMES=""
MIN_CHROMOSOME_LENGTH=""
CHROMOSOME=""
REGION_START=""
REGION_END=""
WINDOW_SIZE=""
STATISTICS="pi,theta,tajima_d"
SAMPLE_PAIRS=""
TOP_N_EXTREME=""
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

Required (at least one):
  --fst-dir DIR              Directory containing FST TSV files (from calculate_fst.sh)
  --diversity-dir DIR        Directory containing diversity TSV files (from calculate_pi_theta.sh)
  --hdf5-dir DIR             Directory containing collated HDF5 (diversity_w*.h5, fst_w*.h5)

Required:
  --output-dir DIR           Output directory for results
  Either single-threshold or seed-then-expand quantiles (see below).

Single-threshold mode (require both):
  --high-quantile N          High quantile threshold (e.g., 0.99 for top 1%)
  --low-quantile N           Low quantile threshold (e.g., 0.01 for bottom 1%)

Seed-then-expand mode (require all four; overrides single-threshold when set):
  --seed-high-quantile N     Strict high quantile for seed windows
  --seed-low-quantile N      Strict low quantile for seed windows
  --expand-high-quantile N   Soft high quantile for expanding regions
  --expand-low-quantile N    Soft low quantile for expanding regions

Optional:
  --reference-genome FILE    Reference genome FASTA file (optional, for chromosome lengths)
  --top-n-chromosomes N      Process only the N longest chromosomes/scaffolds (optional)
  --min-chromosome-length N  Process only chromosomes/scaffolds longer than N bp (optional)
  --chromosome CHR           Restrict to a single chromosome/scaffold (optional)
  --region-start N           Restrict to interval start on --chromosome (requires --chromosome)
  --region-end N             Restrict to interval end on --chromosome (requires --chromosome)
  --window-size NUMBER       Filter to specific window size (default: all window sizes)
  --statistics LIST          Comma-separated list: pi, theta, tajima_d, pbe, fst (default: pi,theta,tajima_d)
  --sample-pairs LIST        Comma-separated list of FST pairs to process (default: all)
  --top-n-extreme N          Return only the N most extreme values in each quantile (optional)
  --merge-distance N         Merge nearby windows into regions: bp or 'auto' = max(2*window, 2000) [default: auto]
  --min-depth N              Minimum mean coverage for a window to be a seed or in expansion (optional)
  --max-depth N              Maximum mean coverage for a window to be a seed or in expansion (optional)
  --min-mapping-quality N    Minimum mean mapping quality for a window (optional)
  --max-mapping-quality N    Maximum mean mapping quality for a window (optional)
  --min-snps N               Minimum SNPs per window for seed/expansion (optional)
  --max-snps N               Maximum SNPs per window for seed/expansion (optional)
  --min-region-snps N        Minimum total SNPs in a region to output (optional)
  --max-region-snps N        Maximum total SNPs in a region to output (optional)
  --verbose                  Enable verbose output for debugging
  --rscript PATH             Path to identify_outliers.R (default: same directory as this script)
  --dry-run                  Preview commands without executing (dry-run mode)
  -h, --help                 Show this help message

Examples:
  # Identify outliers from both FST and diversity statistics
  $0 \\
    --fst-dir ../stats/fst_output \\
    --diversity-dir ../stats/diversity_output \\
    --output-dir ./outlier_results \\
    --high-quantile 0.99 \\
    --low-quantile 0.01 \\
    --top-n-chromosomes 10 \\
    --reference-genome ../worm_q10.medaka.purged.fa

  # Process specific statistics and window size
  $0 \\
    --diversity-dir ../stats/diversity_output \\
    --output-dir ./outlier_results \\
    --statistics pi,theta \\
    --window-size 10000 \\
    --min-chromosome-length 500000 \\
    --high-quantile 0.95 \\
    --low-quantile 0.05

Output files:
  - {output_dir}/outlier_windows*.csv: Wide table (chr, start, end, samplename.stat columns, outlier_stat)
  - {output_dir}/outlier_regions*.csv: Merged/expanded regions when --merge-distance or seed-expand used

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
        --fst-dir)
            FST_DIR="$2"
            shift 2
            ;;
        --diversity-dir)
            DIVERSITY_DIR="$2"
            shift 2
            ;;
        --hdf5-dir)
            HDF5_DIR="$2"
            shift 2
            ;;
        --merge-distance)
            MERGE_DISTANCE="$2"
            shift 2
            ;;
        --min-depth)
            MIN_DEPTH="$2"
            shift 2
            ;;
        --max-depth)
            MAX_DEPTH="$2"
            shift 2
            ;;
        --min-mapping-quality)
            MIN_MAPPING_QUALITY="$2"
            shift 2
            ;;
        --max-mapping-quality)
            MAX_MAPPING_QUALITY="$2"
            shift 2
            ;;
        --min-snps)
            MIN_SNPS="$2"
            shift 2
            ;;
        --max-snps)
            MAX_SNPS="$2"
            shift 2
            ;;
        --min-region-snps)
            MIN_REGION_SNPS="$2"
            shift 2
            ;;
        --max-region-snps)
            MAX_REGION_SNPS="$2"
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
        --high-quantile)
            HIGH_QUANTILE="$2"
            shift 2
            ;;
        --low-quantile)
            LOW_QUANTILE="$2"
            shift 2
            ;;
        --seed-high-quantile)
            SEED_HIGH_QUANTILE="$2"
            shift 2
            ;;
        --seed-low-quantile)
            SEED_LOW_QUANTILE="$2"
            shift 2
            ;;
        --expand-high-quantile)
            EXPAND_HIGH_QUANTILE="$2"
            shift 2
            ;;
        --expand-low-quantile)
            EXPAND_LOW_QUANTILE="$2"
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
        --chromosome)
            CHROMOSOME="$2"
            shift 2
            ;;
        --region-start)
            REGION_START="$2"
            shift 2
            ;;
        --region-end)
            REGION_END="$2"
            shift 2
            ;;
        --window-size)
            WINDOW_SIZE="$2"
            shift 2
            ;;
        --statistics)
            STATISTICS="$2"
            shift 2
            ;;
        --sample-pairs)
            SAMPLE_PAIRS="$2"
            shift 2
            ;;
        --top-n-extreme)
            TOP_N_EXTREME="$2"
            shift 2
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
if [[ -z "$FST_DIR" ]] && [[ -z "$DIVERSITY_DIR" ]] && [[ -z "$HDF5_DIR" ]]; then
    log_error "At least one of --fst-dir, --diversity-dir, or --hdf5-dir must be specified"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing required argument: --output-dir"
    usage
    exit 1
fi

# Require EITHER (--high-quantile and --low-quantile) OR (all four seed/expand quantiles)
HAVE_HIGH_LOW=false
HAVE_SEED_EXPAND=false
[[ -n "$HIGH_QUANTILE" && -n "$LOW_QUANTILE" ]] && HAVE_HIGH_LOW=true
[[ -n "$SEED_HIGH_QUANTILE" && -n "$SEED_LOW_QUANTILE" && -n "$EXPAND_HIGH_QUANTILE" && -n "$EXPAND_LOW_QUANTILE" ]] && HAVE_SEED_EXPAND=true

if [[ "$HAVE_HIGH_LOW" == false && "$HAVE_SEED_EXPAND" == false ]]; then
    log_error "Either (--high-quantile and --low-quantile) or (--seed-high-quantile, --seed-low-quantile, --expand-high-quantile, --expand-low-quantile) are required"
    usage
    exit 1
fi

# Validate quantile values when provided (must be between 0 and 1)
if [[ -n "$HIGH_QUANTILE" ]]; then
    if ! awk "BEGIN {exit !($HIGH_QUANTILE >= 0 && $HIGH_QUANTILE <= 1)}"; then
        log_error "Invalid --high-quantile value: $HIGH_QUANTILE (must be between 0 and 1)"
        exit 1
    fi
fi
if [[ -n "$LOW_QUANTILE" ]]; then
    if ! awk "BEGIN {exit !($LOW_QUANTILE >= 0 && $LOW_QUANTILE <= 1)}"; then
        log_error "Invalid --low-quantile value: $LOW_QUANTILE (must be between 0 and 1)"
        exit 1
    fi
fi

if [[ ( -n "$REGION_START" || -n "$REGION_END" ) && -z "$CHROMOSOME" ]]; then
    log_error "--region-start and --region-end require --chromosome"
    exit 1
fi

# Find R script
if [[ -z "$RSCRIPT" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    RSCRIPT="${SCRIPT_DIR}/identify_outliers.R"
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

# Check input directories
if [[ "$DRY_RUN" == false ]]; then
    if [[ -n "$FST_DIR" ]]; then
        check_dir "$FST_DIR"
        tsv_count=$(find "$FST_DIR" -name "*fst*.tsv" -o -name "*fst*.csv" | wc -l)
        if [[ $tsv_count -eq 0 ]]; then
            log_warn "No FST TSV/CSV files found in: $FST_DIR"
            log_warn "Expected files matching pattern: *fst*.tsv or *fst*.csv"
        else
            log "Found $tsv_count FST file(s)"
        fi
    fi
    
    if [[ -n "$DIVERSITY_DIR" ]]; then
        check_dir "$DIVERSITY_DIR"
        tsv_count=$(find "$DIVERSITY_DIR" -name "*diversity*.tsv" -o -name "*diversity*.csv" 2>/dev/null | wc -l)
        if [[ $tsv_count -eq 0 ]]; then
            log_warn "No diversity TSV/CSV files found in: $DIVERSITY_DIR"
            log_warn "Expected files matching pattern: *diversity*.tsv or *diversity*.csv"
        else
            log "Found $tsv_count diversity file(s)"
        fi
    fi
    if [[ -n "$HDF5_DIR" ]]; then
        check_dir "$HDF5_DIR"
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
    --output-dir "$OUTPUT_DIR"
    --statistics "$STATISTICS"
)
if [[ -n "$HIGH_QUANTILE" ]]; then
    R_CMD+=(--high-quantile "$HIGH_QUANTILE")
fi
if [[ -n "$LOW_QUANTILE" ]]; then
    R_CMD+=(--low-quantile "$LOW_QUANTILE")
fi

if [[ -n "$SEED_HIGH_QUANTILE" ]]; then
    R_CMD+=(--seed-high-quantile "$SEED_HIGH_QUANTILE")
fi
if [[ -n "$SEED_LOW_QUANTILE" ]]; then
    R_CMD+=(--seed-low-quantile "$SEED_LOW_QUANTILE")
fi
if [[ -n "$EXPAND_HIGH_QUANTILE" ]]; then
    R_CMD+=(--expand-high-quantile "$EXPAND_HIGH_QUANTILE")
fi
if [[ -n "$EXPAND_LOW_QUANTILE" ]]; then
    R_CMD+=(--expand-low-quantile "$EXPAND_LOW_QUANTILE")
fi

if [[ -n "$FST_DIR" ]]; then
    R_CMD+=(--fst-dir "$FST_DIR")
fi

if [[ -n "$DIVERSITY_DIR" ]]; then
    R_CMD+=(--diversity-dir "$DIVERSITY_DIR")
fi

if [[ -n "$HDF5_DIR" ]]; then
    R_CMD+=(--hdf5-dir "$HDF5_DIR")
fi

if [[ -n "$MERGE_DISTANCE" ]]; then
    R_CMD+=(--merge-distance "$MERGE_DISTANCE")
    [[ -n "$MIN_DEPTH" ]] && R_CMD+=(--min-depth "$MIN_DEPTH")
    [[ -n "$MAX_DEPTH" ]] && R_CMD+=(--max-depth "$MAX_DEPTH")
    [[ -n "$MIN_MAPPING_QUALITY" ]] && R_CMD+=(--min-mapping-quality "$MIN_MAPPING_QUALITY")
    [[ -n "$MAX_MAPPING_QUALITY" ]] && R_CMD+=(--max-mapping-quality "$MAX_MAPPING_QUALITY")
    [[ -n "$MIN_SNPS" ]] && R_CMD+=(--min-snps "$MIN_SNPS")
    [[ -n "$MAX_SNPS" ]] && R_CMD+=(--max-snps "$MAX_SNPS")
    [[ -n "$MIN_REGION_SNPS" ]] && R_CMD+=(--min-region-snps "$MIN_REGION_SNPS")
    [[ -n "$MAX_REGION_SNPS" ]] && R_CMD+=(--max-region-snps "$MAX_REGION_SNPS")
fi

if [[ -n "$REFERENCE_GENOME" ]]; then
    R_CMD+=(--reference-genome "$REFERENCE_GENOME")
fi

if [[ -n "$TOP_N_CHROMOSOMES" ]]; then
    R_CMD+=(--top-n-chromosomes "$TOP_N_CHROMOSOMES")
fi

if [[ -n "$MIN_CHROMOSOME_LENGTH" ]]; then
    R_CMD+=(--min-chromosome-length "$MIN_CHROMOSOME_LENGTH")
fi
if [[ -n "$CHROMOSOME" ]]; then
    R_CMD+=(--chromosome "$CHROMOSOME")
fi
if [[ -n "$REGION_START" ]]; then
    R_CMD+=(--region-start "$REGION_START")
fi
if [[ -n "$REGION_END" ]]; then
    R_CMD+=(--region-end "$REGION_END")
fi

if [[ -n "$WINDOW_SIZE" ]]; then
    R_CMD+=(--window-size "$WINDOW_SIZE")
fi

if [[ -n "$SAMPLE_PAIRS" ]]; then
    R_CMD+=(--sample-pairs "$SAMPLE_PAIRS")
fi

if [[ -n "$TOP_N_EXTREME" ]]; then
    R_CMD+=(--top-n-extreme "$TOP_N_EXTREME")
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
    if [[ -n "$FST_DIR" ]]; then
        log_dry_run "FST directory: $FST_DIR"
    fi
    if [[ -n "$DIVERSITY_DIR" ]]; then
        log_dry_run "Diversity directory: $DIVERSITY_DIR"
    fi
    if [[ -n "$HDF5_DIR" ]]; then
        log_dry_run "HDF5 directory: $HDF5_DIR"
    fi
    log_dry_run "Output directory: $OUTPUT_DIR"
    log_dry_run "High quantile: $HIGH_QUANTILE"
    log_dry_run "Low quantile: $LOW_QUANTILE"
    log_dry_run "Statistics: $STATISTICS"
    if [[ -n "$TOP_N_CHROMOSOMES" ]]; then
        log_dry_run "Top N chromosomes: $TOP_N_CHROMOSOMES"
    fi
    if [[ -n "$MIN_CHROMOSOME_LENGTH" ]]; then
        log_dry_run "Min chromosome length: $MIN_CHROMOSOME_LENGTH"
    fi
    if [[ -n "$WINDOW_SIZE" ]]; then
        log_dry_run "Window size: $WINDOW_SIZE"
    fi
    if [[ -n "$SAMPLE_PAIRS" ]]; then
        log_dry_run "Sample pairs: $SAMPLE_PAIRS"
    fi
    if [[ -n "$TOP_N_EXTREME" ]]; then
        log_dry_run "Top N extreme: $TOP_N_EXTREME"
    fi
    if [[ "$VERBOSE" == true ]]; then
        log_dry_run "Verbose: enabled"
    fi
else
    LOG_DIR="${OUTPUT_DIR}/log"
    mkdir -p "$LOG_DIR"
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    LOG_FILE="${LOG_DIR}/identify_outliers_${TIMESTAMP}.log"
    {
        echo "Command: ${R_CMD[*]}"
        echo "Started: $(date -Iseconds 2>/dev/null || date)"
    } >> "$LOG_FILE"
    log "Log file: $LOG_FILE"
    log "Running identify_outliers.R..."
    if [[ -n "$FST_DIR" ]]; then
        log "  FST directory: $FST_DIR"
    fi
    if [[ -n "$DIVERSITY_DIR" ]]; then
        log "  Diversity directory: $DIVERSITY_DIR"
    fi
    if [[ -n "$HDF5_DIR" ]]; then
        log "  HDF5 directory: $HDF5_DIR"
    fi
    log "  Output directory: $OUTPUT_DIR"
    log "  High quantile: $HIGH_QUANTILE"
    log "  Low quantile: $LOW_QUANTILE"
    log "  Statistics: $STATISTICS"
    if [[ -n "$TOP_N_CHROMOSOMES" ]]; then
        log "  Top N chromosomes: $TOP_N_CHROMOSOMES"
    fi
    if [[ -n "$MIN_CHROMOSOME_LENGTH" ]]; then
        log "  Min chromosome length: $MIN_CHROMOSOME_LENGTH"
    fi
    if [[ -n "$WINDOW_SIZE" ]]; then
        log "  Window size: $WINDOW_SIZE"
    fi
    if [[ -n "$SAMPLE_PAIRS" ]]; then
        log "  Sample pairs: $SAMPLE_PAIRS"
    fi
    if [[ -n "$TOP_N_EXTREME" ]]; then
        log "  Top N extreme: $TOP_N_EXTREME"
    fi
    if [[ "$VERBOSE" == true ]]; then
        log "  Verbose: enabled"
    fi

    "${R_CMD[@]}" 2>&1 | tee -a "$LOG_FILE" || {
        log_error "R script failed"
        exit 1
    }

    for f in "$OUTPUT_DIR"/outlier_*.csv; do
        if [[ -f "$f" ]]; then
            bn=$(basename "$f")
            { echo ""; log "First rows of $bn:"; head -n 4 "$f"; } | tee -a "$LOG_FILE"
        fi
    done

    log "Outlier identification complete!"
    log "Output directory: $OUTPUT_DIR"
fi
