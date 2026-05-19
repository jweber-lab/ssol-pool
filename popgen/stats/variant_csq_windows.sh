#!/bin/bash

###############################################################################
# variant_csq_windows.sh
#
# Wrapper for variant_csq_windows.R: aggregate annotated sites into
# variant_features_wW_sS.tsv using FST/diversity template windows.
###############################################################################

set -euo pipefail

SITES_TSV=""
TEMPLATE_DIR=""
TEMPLATE_TSV=""
OUTPUT_DIR="."
RSCRIPT="Rscript"
DRY_RUN=false

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"; }
log_error() { echo "[ERROR] $1" >&2; }

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  --sites-tsv FILE          Annotated TSV from variant_csq.R

One of:
  --template-dir DIR        Directory with *fst*w*_s*.tsv (or diversity window TSVs)
  --template-tsv FILE       Single template TSV

Optional:
  --output-dir DIR          Where to write variant_features_wW_sS.tsv [default: .]
  --rscript PATH            [default: Rscript]
  --dry-run
  -h, --help

Example:
  $0 --sites-tsv variants_csq.tsv --template-dir fst_out --output-dir vf_out
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sites-tsv) SITES_TSV="$2"; shift 2 ;;
        --template-dir) TEMPLATE_DIR="$2"; shift 2 ;;
        --template-tsv) TEMPLATE_TSV="$2"; shift 2 ;;
        --output-dir|-o) OUTPUT_DIR="$2"; shift 2 ;;
        --rscript) RSCRIPT="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
WIN_R="${SCRIPT_DIR}/variant_csq_windows.R"
if [[ ! -f "$WIN_R" ]]; then
    log_error "variant_csq_windows.R not found: $WIN_R"
    exit 1
fi

if [[ -z "$SITES_TSV" || ! -f "$SITES_TSV" ]]; then
    log_error "Missing or invalid --sites-tsv"
    usage
    exit 1
fi
if [[ -z "$TEMPLATE_DIR" && -z "$TEMPLATE_TSV" ]]; then
    log_error "Provide --template-dir or --template-tsv"
    usage
    exit 1
fi

ARGS=(--sites-tsv "$SITES_TSV" --output-dir "$OUTPUT_DIR")
[[ -n "$TEMPLATE_DIR" ]] && ARGS+=(--template-dir "$TEMPLATE_DIR")
[[ -n "$TEMPLATE_TSV" ]] && ARGS+=(--template-tsv "$TEMPLATE_TSV")

if [[ "$DRY_RUN" == true ]]; then
    echo "Would run: $RSCRIPT $WIN_R ${ARGS[*]}"
    exit 0
fi

log "Running variant_csq_windows.R"
"$RSCRIPT" "$WIN_R" "${ARGS[@]}"
log "Done."
