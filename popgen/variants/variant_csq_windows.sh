#!/bin/bash

###############################################################################
# variant_csq_windows.sh
#
# Wrapper for variant_csq_windows.R: aggregate annotated sites into
# variant_features_wW_sS.tsv on grenedalf's interval window grid (built from the reference
# .fai), or on windows listed in FST/diversity template files.
###############################################################################

set -euo pipefail

SITES_TSVS=()
REFERENCE=""
WINDOW_SIZES=()
STEP_SIZES=()
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
  --sites-tsv FILE          Annotated TSV from variant_csq.R; repeat to combine SNP and indel outputs

Windows (recommended: reference grid, identical to diversity / FST / PBE / seq_qual windows):
  --reference FILE          Reference FASTA (with .fai) or .fai; windows from chromosome lengths
  -w, --window-size N       Window size (repeatable; paired with --step-size)
  --step-size N             Step size (repeatable; default: window/2)

Or / additionally:
  --template-dir DIR        Directory with *fst*w*_s*.tsv (or diversity window TSVs); scales come
                            from file names. Without --reference, only the windows in these files
                            are used (chromosomes grenedalf skipped are missing).
  --template-tsv FILE       Single template TSV

Optional:
  --output-dir DIR          Where to write variant_features_wW_sS.tsv [default: .]
  --rscript PATH            [default: Rscript]
  --dry-run
  -h, --help

Example:
  $0 --sites-tsv csq_snps/variants_csq.tsv --sites-tsv csq_indels/variants_csq.tsv \\
     --reference ref.unmasked.fa -w 1000 --step-size 500 -w 10000 --step-size 5000 --output-dir vf_out
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sites-tsv) SITES_TSVS+=("$2"); shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        -w|--window-size) WINDOW_SIZES+=("$2"); shift 2 ;;
        --step-size) STEP_SIZES+=("$2"); shift 2 ;;
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

if [[ ${#SITES_TSVS[@]} -eq 0 ]]; then
    log_error "Missing --sites-tsv"
    usage
    exit 1
fi
for f in "${SITES_TSVS[@]}"; do
    if [[ ! -f "$f" ]]; then
        log_error "Sites TSV not found: $f"
        exit 1
    fi
    if [[ "$f" == *,* ]]; then
        log_error "Sites TSV path contains a comma (not supported): $f"
        exit 1
    fi
done
SITES_TSV=$(IFS=,; echo "${SITES_TSVS[*]}")
if [[ -z "$TEMPLATE_DIR" && -z "$TEMPLATE_TSV" && ( -z "$REFERENCE" || ${#WINDOW_SIZES[@]} -eq 0 ) ]]; then
    log_error "Provide --reference with --window-size (recommended), or --template-dir / --template-tsv"
    usage
    exit 1
fi
if [[ ${#STEP_SIZES[@]} -gt 0 && ${#STEP_SIZES[@]} -ne ${#WINDOW_SIZES[@]} ]]; then
    log_error "--step-size must be given once per --window-size (or omitted for window/2)"
    exit 1
fi

ARGS=(--sites-tsv "$SITES_TSV" --output-dir "$OUTPUT_DIR")
[[ -n "$REFERENCE" ]] && ARGS+=(--reference "$REFERENCE")
[[ ${#WINDOW_SIZES[@]} -gt 0 ]] && ARGS+=(--window-size "$(IFS=,; echo "${WINDOW_SIZES[*]}")")
[[ ${#STEP_SIZES[@]} -gt 0 ]] && ARGS+=(--step-size "$(IFS=,; echo "${STEP_SIZES[*]}")")
[[ -n "$TEMPLATE_DIR" ]] && ARGS+=(--template-dir "$TEMPLATE_DIR")
[[ -n "$TEMPLATE_TSV" ]] && ARGS+=(--template-tsv "$TEMPLATE_TSV")

if [[ "$DRY_RUN" == true ]]; then
    echo "Would run: $RSCRIPT $WIN_R ${ARGS[*]}"
    exit 0
fi

log "Running variant_csq_windows.R"
"$RSCRIPT" "$WIN_R" "${ARGS[@]}"
log "Done."
