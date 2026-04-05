#!/bin/bash

###############################################################################
# plot_hist.sh
#
# Bash wrapper for plot_hist.R to plot genome-wide (or subset) distributions
# for a chosen statistic. Default: histogram only; add panels with --panels.
#
# Mirrors plot_region.sh inputs so you can point it at the same directories:
#   - TSV dirs: --diversity-dir / --fst-dir / --pbe-dir / --seq-qual-dir
#   - Collated HDF5 dir: --hdf5-dir
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
OUTPUT_DIR="./"
FILE_PREFIX=""
STAT="pi"
Y_VALUE="value"
TRANSFORM="none"
BINS=60
PANELS="histogram"
HIST_MODE="count"
OVERLAY=false
MAX_GROUPS=30
Q_LINES="0.5,0.95,0.99"
WIDTH=12
HEIGHT=8
DPI=300
PLOT_FORMAT="png"
RSCRIPT=""
VERBOSE=false
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
  cat << EOF
Usage: $0 [OPTIONS]

Subset (optional):
  --chromosome CHR         Subset to chromosome/scaffold CHR
  --region CHR:START-END   Subset to region (e.g. ptg000624l:1-5000000)

Input (at least one required; use TSV dirs and/or --hdf5-dir):
  --diversity-dir DIR      Directory with diversity TSV or HDF5
  --fst-dir DIR            Directory with FST TSV or HDF5
  --pbe-dir DIR            Directory with PBE TSV or HDF5
  --seq-qual-dir DIR       Directory with seq_qual TSV (mean_coverage, mean_mapping_quality)
  --hdf5-dir DIR           Collate output dir: diversity_*.h5, fst_*.h5, pbe_*.h5

Optional:
  --window-size N          Window size to use (default: first available)
  --step-size N            Step size (optional)
  --output-dir DIR         Output directory (default: ./)
  --file-prefix PREFIX     Prefix for output filename
  --stat STAT              Statistic: coverage, mapq, n_snps, pi, theta, tajima_d, fst, pbe (default: pi)
  --y-value VALUE          For diversity/FST/PBE: value, rank, or quantile (default: value)
  --transform TRANSFORM    Transform on stat axis: none, log, or asinh (default: none)
  --bins N                 Histogram bins (default: 60)
  --panels LIST            Comma-separated figure panels: histogram, density, ecdf (default: histogram).
                            \"density\" is a kernel-density curve panel, not the same as --hist-mode.
  --hist-mode MODE         Histogram bar y-axis: count or density (ggplot after_stat(density) for bars only; default: count)
  --overlay                Overlay groups in one panel (color/fill = group); otherwise facet by group
  --max-groups N           Refuse to plot >N groups (default: 30)
  --q-lines LIST           Quantile lines to mark (comma-separated; default: 0.5,0.95,0.99)
  --width N                Figure width in inches (default: 12)
  --height N               Figure height in inches (default: 8)
  --dpi N                  DPI for PNG (default: 300)
  --plot-format FORMAT     png, pdf, svg, both, or all (default: png)
  --rscript PATH           Path to plot_hist.R (default: same directory as this script)
  --verbose                Verbose logging
  --dry-run                Preview command without executing
  -h, --help               Show this help

Examples:
  # Genome-wide pi histogram from collate output:
  $0 --hdf5-dir ./collated --stat pi --output-dir ./plots --overlay

  # Histogram + kernel density + ECDF (stacked figure):
  $0 --hdf5-dir ./collated --stat pi --panels histogram,density,ecdf -o ./plots

  # Bar heights normalized (still histogram panel only; not the density curve panel):
  $0 --hdf5-dir ./collated --stat coverage --hist-mode density -o ./hist/

  # Coverage: log x-axis transform and kernel density panel:
  $0 --seq-qual-dir ./stats_out --stat coverage --transform log --panels histogram,density -o ./plots

  # Subset to chromosome and use ranks:
  $0 --hdf5-dir ./collated --chromosome ptg000624l --stat fst --y-value rank --overlay -o ./plots
EOF
}

log() { echo -e "${GREEN}[plot_hist]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[plot_hist]${NC} $1"; }
log_error() { echo -e "${RED}[plot_hist]${NC} $1"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chromosome) CHROMOSOME="$2"; shift 2 ;;
    --region) REGION="$2"; shift 2 ;;
    --diversity-dir) DIVERSITY_DIR="$2"; shift 2 ;;
    --fst-dir) FST_DIR="$2"; shift 2 ;;
    --pbe-dir) PBE_DIR="$2"; shift 2 ;;
    --seq-qual-dir) SEQ_QUAL_DIR="$2"; shift 2 ;;
    --hdf5-dir) HDF5_DIR="$2"; shift 2 ;;
    --window-size) WINDOW_SIZE="$2"; shift 2 ;;
    --step-size) STEP_SIZE="$2"; shift 2 ;;
    --output-dir|-o) OUTPUT_DIR="$2"; shift 2 ;;
    --file-prefix) FILE_PREFIX="$2"; shift 2 ;;
    --stat) STAT="$2"; shift 2 ;;
    --y-value) Y_VALUE="$2"; shift 2 ;;
    --transform) TRANSFORM="$2"; shift 2 ;;
    --bins) BINS="$2"; shift 2 ;;
    --panels) PANELS="$2"; shift 2 ;;
    --hist-mode) HIST_MODE="$2"; shift 2 ;;
    --overlay) OVERLAY=true; shift ;;
    --max-groups) MAX_GROUPS="$2"; shift 2 ;;
    --q-lines) Q_LINES="$2"; shift 2 ;;
    --width) WIDTH="$2"; shift 2 ;;
    --height) HEIGHT="$2"; shift 2 ;;
    --dpi) DPI="$2"; shift 2 ;;
    --plot-format) PLOT_FORMAT="$2"; shift 2 ;;
    --rscript) RSCRIPT="$2"; shift 2 ;;
    --verbose) VERBOSE=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) log_error "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$DIVERSITY_DIR" && -z "$FST_DIR" && -z "$PBE_DIR" && -z "$HDF5_DIR" && -z "$SEQ_QUAL_DIR" ]]; then
  log_error "At least one input is required: --diversity-dir, --fst-dir, --pbe-dir, --seq-qual-dir, or --hdf5-dir"
  usage
  exit 1
fi

if [[ -z "$RSCRIPT" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  RSCRIPT="${SCRIPT_DIR}/plot_hist.R"
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

R_CMD+=(--output-dir "$OUTPUT_DIR" --width "$WIDTH" --height "$HEIGHT" --dpi "$DPI" --plot-format "$PLOT_FORMAT")
[[ -n "$FILE_PREFIX" ]] && R_CMD+=(--file-prefix "$FILE_PREFIX")
[[ -n "$STAT" ]] && R_CMD+=(--stat "$STAT")
[[ -n "$Y_VALUE" ]] && R_CMD+=(--y-value "$Y_VALUE")
[[ -n "$TRANSFORM" ]] && R_CMD+=(--transform "$TRANSFORM")
R_CMD+=(--bins "$BINS" --panels "$PANELS" --hist-mode "$HIST_MODE" --max-groups "$MAX_GROUPS" --q-lines "$Q_LINES")
[[ "$OVERLAY" == true ]] && R_CMD+=(--overlay)
[[ "$VERBOSE" == true ]] && R_CMD+=(--verbose)

if [[ "$DRY_RUN" == true ]]; then
  log "DRY-RUN: ${R_CMD[*]}"
  exit 0
fi

log "Running plot_hist.R -> $OUTPUT_DIR"
"${R_CMD[@]}"
log "Done."

