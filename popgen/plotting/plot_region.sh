#!/bin/bash

###############################################################################
# plot_region.sh
#
# Bash wrapper for plot_region.R to create a single stacked figure for a
# genomic region (coverage, π, FST, PBE) with shared x-axis and color key.
#
# Author: ssol-pool
# Usage: See README.md or run with --help
###############################################################################

set -euo pipefail

CHROMOSOME=""
REGION=""
DIVERSITY_DIR=""
FST_DIR=""
PBE_DIR=""
SEQ_QUAL_DIR=""
HDF5_DIR=""
WINDOW_SIZE=""
STEP_SIZE=""
REFERENCE_GENOME=""
OUTPUT_DIR="./"
FILE_PREFIX=""
Y_VALUE="value"
WIDTH=12
HEIGHT=10
DPI=300
PLOT_FORMAT="png"
RSCRIPT=""
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Region (at least one required):
  --chromosome CHR           Plot full chromosome CHR
  --region CHR:START-END     Plot region (e.g. chr1:1000000-2000000)

Input (at least one required; use TSV dirs and/or --hdf5-dir):
  --diversity-dir DIR        Directory with diversity TSV or HDF5
  --fst-dir DIR              Directory with FST TSV or HDF5
  --pbe-dir DIR              Directory with PBE TSV or HDF5
  --seq-qual-dir DIR         Directory with seq_qual TSV (mean_coverage, mean_mapping_quality)
  --hdf5-dir DIR             Collate output dir: diversity_*.h5, fst_*.h5, pbe_*.h5 (aligned windows)

Optional:
  --window-size N            Window size to use (default: first available)
  --step-size N              Step size (optional)
  --reference-genome FILE    Reference genome FASTA (for chromosome lengths)
  --output-dir DIR           Output directory (default: ./)
  --file-prefix PREFIX        Prefix for output filename
  --y-value VALUE            Y-axis for diversity/FST/PBE: value, rank, or quantile (default: value)
  --width N                  Figure width in inches (default: 12)
  --height N                 Figure height in inches (default: 10)
  --dpi N                    DPI for PNG (default: 300)
  --plot-format FORMAT        png, pdf, svg, both, or all (default: png)
  --rscript PATH             Path to plot_region.R (default: same directory as this script)
  --dry-run                  Preview commands without executing
  -h, --help                 Show this help

Output:
  region_plot_CHR_START_END.png (or .pdf/.svg) in --output-dir.
  Panel order: Coverage (if seq_qual/diversity HDF5), π, FST, PBE.

Example:
  $0 --chromosome chr1 --hdf5-dir ./collate_out --reference-genome ref.fa --output-dir ./plots
  $0 --region chr2:5000000-6000000 --diversity-dir ./div --fst-dir ./fst --seq-qual-dir ./sq
EOF
}

log() { echo -e "${GREEN}[plot_region]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[plot_region]${NC} $1"; }
log_error() { echo -e "${RED}[plot_region]${NC} $1"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chromosome)
            CHROMOSOME="$2"
            shift 2
            ;;
        --region)
            REGION="$2"
            shift 2
            ;;
        --diversity-dir)
            DIVERSITY_DIR="$2"
            shift 2
            ;;
        --fst-dir)
            FST_DIR="$2"
            shift 2
            ;;
        --pbe-dir)
            PBE_DIR="$2"
            shift 2
            ;;
        --seq-qual-dir)
            SEQ_QUAL_DIR="$2"
            shift 2
            ;;
        --hdf5-dir)
            HDF5_DIR="$2"
            shift 2
            ;;
        --window-size)
            WINDOW_SIZE="$2"
            shift 2
            ;;
        --step-size)
            STEP_SIZE="$2"
            shift 2
            ;;
        --reference-genome)
            REFERENCE_GENOME="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --file-prefix)
            FILE_PREFIX="$2"
            shift 2
            ;;
        --y-value)
            Y_VALUE="$2"
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
        --plot-format)
            PLOT_FORMAT="$2"
            shift 2
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

if [[ -z "$CHROMOSOME" && -z "$REGION" ]]; then
    log_error "At least one of --chromosome or --region is required"
    usage
    exit 1
fi

if [[ -z "$DIVERSITY_DIR" && -z "$FST_DIR" && -z "$PBE_DIR" && -z "$HDF5_DIR" ]]; then
    log_error "At least one input is required: --diversity-dir, --fst-dir, --pbe-dir, or --hdf5-dir"
    usage
    exit 1
fi

if [[ -z "$RSCRIPT" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    RSCRIPT="${SCRIPT_DIR}/plot_region.R"
fi

if [[ ! -f "$RSCRIPT" ]]; then
    log_error "R script not found: $RSCRIPT"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

R_CMD=(Rscript "$RSCRIPT")
[[ -n "$CHROMOSOME" ]] && R_CMD+=(--chromosome "$CHROMOSOME")
[[ -n "$REGION" ]] && R_CMD+=(--region "$REGION")
[[ -n "$DIVERSITY_DIR" ]] && R_CMD+=(--diversity-dir "$DIVERSITY_DIR")
[[ -n "$FST_DIR" ]] && R_CMD+=(--fst-dir "$FST_DIR")
[[ -n "$PBE_DIR" ]] && R_CMD+=(--pbe-dir "$PBE_DIR")
[[ -n "$SEQ_QUAL_DIR" ]] && R_CMD+=(--seq-qual-dir "$SEQ_QUAL_DIR")
[[ -n "$HDF5_DIR" ]] && R_CMD+=(--hdf5-dir "$HDF5_DIR")
[[ -n "$WINDOW_SIZE" ]] && R_CMD+=(--window-size "$WINDOW_SIZE")
[[ -n "$STEP_SIZE" ]] && R_CMD+=(--step-size "$STEP_SIZE")
[[ -n "$REFERENCE_GENOME" ]] && R_CMD+=(--reference-genome "$REFERENCE_GENOME")
R_CMD+=(--output-dir "$OUTPUT_DIR" --width "$WIDTH" --height "$HEIGHT" --dpi "$DPI" --plot-format "$PLOT_FORMAT")
[[ -n "$FILE_PREFIX" ]] && R_CMD+=(--file-prefix "$FILE_PREFIX")
[[ -n "$Y_VALUE" ]] && R_CMD+=(--y-value "$Y_VALUE")

if [[ "$DRY_RUN" == true ]]; then
    log "DRY-RUN: ${R_CMD[*]}"
    exit 0
fi

"${R_CMD[@]}"
log "Done."
