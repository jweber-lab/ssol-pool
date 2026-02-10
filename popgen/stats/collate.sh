#!/bin/bash

###############################################################################
# collate.sh
#
# Wrapper for collate.R: discovers diversity, FST, and PBE TSV/CSV (from
# calculate_pi_theta.sh, calculate_fst.sh, and PBE outputs) and writes
# combined HDF5 plus *_summary.tsv to OUTPUT_DIR. Optional: --seq-qual-dir
# for mean_coverage/mean_mapping_quality, --variant-tsv-dir for n_snps or
# variants.h5. Missing/empty FST or PBE dirs are skipped.
# For full output file list and variable naming (HDF5 groups, TSV columns),
# see the header of collate.R.
# Usage: See below or --help
###############################################################################

set -euo pipefail

DIVERSITY_DIR=""
FST_DIR=""
PBE_DIR=""
SEQ_QUAL_DIR=""
VARIANT_TSV_DIR=""
OUTPUT_DIR=""
RSCRIPT=""
DRY_RUN=false
SINGLE_POSITION_MERGED=false
NO_SUMMARY=false
DROP_ALL_NA=false
VERBOSE=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Input (at least one required):
  --diversity-dir DIR    Directory containing diversity TSV/CSV (e.g. output of calculate_pi_theta.sh)
  --fst-dir DIR          Directory containing FST TSV/CSV (output of calculate_fst.sh)
  --pbe-dir DIR          Directory containing PBE TSV/CSV (output of calculate_pbe)

Optional:
  --seq-qual-dir DIR     Directory with seq_qual_metrics TSV (from seq_qual_metrics.sh)
  --variant-tsv-dir DIR  Directory with variant/sites TSV (e.g. from bcftools query export of BCF)
  --output-dir DIR       Output directory for HDF5 files [default: ./collated]
  --rscript PATH         Path to Rscript [default: Rscript]
  --single-position-merged  Merge per-locus FST, PBE, and variant TSV into single_position.h5
  --no-summary           Do not write companion *_summary.tsv files
  --drop-all-na          Drop rows/sites where every statistic value (FST or PBE) is NA before ranking and collation
  --verbose              Verbose diagnostics (rank/quantile NA checks, distinct-value checks)
  --dry-run              Preview only
  -h, --help             Show this help

Output naming: diversity_w{W}_s{S}.h5, fst_w{W}_s{S}.h5, fst_single.h5, pbe_w{W}_s{S}.h5, pbe_single.h5,
  optional variants.h5. With --single-position-merged: single_position.h5 with /fst, /pbe_trio_*, /variants.
  Companion *_summary.tsv (unless --no-summary): mean, median, variance, skew, q01, q05, q95, q99, n.

Example:
  $0 --diversity-dir pi_theta_out --fst-dir fst_out --pbe-dir pbe_out -o collated
  $0 --diversity-dir stats_out --fst-dir stats_out --pbe-dir stats_out --seq-qual-dir stats_out --output-dir collated
  $0 --fst-dir fst_out --pbe-dir pbe_out --variant-tsv-dir variants_tsv --single-position-merged -o collated
EOF
}

log() { echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
log_warn() { echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S') WARN]${NC} $1"; }
log_dry_run() { echo -e "${YELLOW}[DRY-RUN]${NC} $1"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --diversity-dir) DIVERSITY_DIR="$2"; shift 2 ;;
        --fst-dir) FST_DIR="$2"; shift 2 ;;
        --pbe-dir) PBE_DIR="$2"; shift 2 ;;
        --seq-qual-dir) SEQ_QUAL_DIR="$2"; shift 2 ;;
        --variant-tsv-dir) VARIANT_TSV_DIR="$2"; shift 2 ;;
        --output-dir|-o) OUTPUT_DIR="$2"; shift 2 ;;
        --rscript) RSCRIPT="$2"; shift 2 ;;
        --single-position-merged) SINGLE_POSITION_MERGED=true; shift ;;
        --no-summary) NO_SUMMARY=true; shift ;;
        --drop-all-na) DROP_ALL_NA=true; shift ;;
        --verbose) VERBOSE=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$DIVERSITY_DIR" && -z "$FST_DIR" && -z "$PBE_DIR" ]]; then
    log_error "At least one of --diversity-dir, --fst-dir, --pbe-dir is required"
    usage
    exit 1
fi

: "${OUTPUT_DIR:=./collated}"
: "${RSCRIPT:=Rscript}"

SCRIPT_DIR=$(dirname "$0")
COLLATE_R="${SCRIPT_DIR}/collate.R"
if [[ ! -f "$COLLATE_R" ]]; then
    log_error "collate.R not found: $COLLATE_R"
    exit 1
fi

ARGS=()
[[ -n "$DIVERSITY_DIR" ]] && ARGS+=(--diversity-dir "$DIVERSITY_DIR")
[[ -n "$FST_DIR" ]] && ARGS+=(--fst-dir "$FST_DIR")
[[ -n "$PBE_DIR" ]] && ARGS+=(--pbe-dir "$PBE_DIR")
[[ -n "$SEQ_QUAL_DIR" ]] && ARGS+=(--seq-qual-dir "$SEQ_QUAL_DIR")
[[ -n "$VARIANT_TSV_DIR" ]] && ARGS+=(--variant-tsv-dir "$VARIANT_TSV_DIR")
ARGS+=(--output-dir "$OUTPUT_DIR")
[[ "$SINGLE_POSITION_MERGED" == true ]] && ARGS+=(--single-position-merged)
[[ "$NO_SUMMARY" == true ]] && ARGS+=(--no-summary)
[[ "$DROP_ALL_NA" == true ]] && ARGS+=(--drop-all-na)
[[ "$VERBOSE" == true ]] && ARGS+=(--verbose)

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "Would run: $RSCRIPT $COLLATE_R ${ARGS[*]}"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/collate_${TIMESTAMP}.log"
{
    echo "Command: $RSCRIPT $COLLATE_R ${ARGS[*]}"
    echo "Started: $(date -Iseconds 2>/dev/null || date)"
} >> "$LOG_FILE"
log "Log file: $LOG_FILE"
log "Running collate.R -> $OUTPUT_DIR"
"$RSCRIPT" "$COLLATE_R" "${ARGS[@]}" 2>&1 | tee -a "$LOG_FILE"
R_EXIT=${PIPESTATUS[0]}
if [[ $R_EXIT -ne 0 ]]; then
    log_error "collate.R exited with code $R_EXIT"
    exit "$R_EXIT"
fi
# Log first rows of each *_summary.tsv
for f in "$OUTPUT_DIR"/*_summary.tsv; do
    if [[ -f "$f" ]]; then
        bn=$(basename "$f")
        { echo ""; log "First rows of $bn:"; head -n 5 "$f"; } | tee -a "$LOG_FILE"
    fi
done
log "Done."
