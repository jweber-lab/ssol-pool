#!/bin/bash

###############################################################################
# variant_csq_windows.sh
#
# Wrapper for variant_csq_windows.R: aggregate annotated sites into
# variant_features_wW_sS.tsv on the same windows as masked diversity / FST (via
# --template-dir), or on grenedalf's interval grid from the reference .fai.
###############################################################################

set -euo pipefail

SITES_TSVS=()
REFERENCE=""
WINDOW_SIZES=()
STEP_SIZES=()
TEMPLATE_DIR=""
TEMPLATE_TSV=""
FILTER_MASK_TOTAL_BED=""
FILTER_MASK_TOTAL_FASTA=""
OUTPUT_DIR="."
RSCRIPT="Rscript"
DRY_RUN=false
NO_NAN_WINDOWS=true

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"; }
log_error() { echo "[ERROR] $1" >&2; }

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  --sites-tsv FILE          Annotated TSV from variant_csq.R; repeat to combine SNP and indel outputs

Windows (use a stats template so rows match masked diversity / FST for collate):
  --reference FILE          Reference FASTA (with .fai) or .fai; builds the interval grid when a
                            scale has no template, and checks template windows against that grid
  -w, --window-size N       Window size (repeatable; paired with --step-size)
  --step-size N             Step size (repeatable; default: window/2)
  --template-dir DIR        Directory of windowed diversity / FST / PBE tables (searched
                            recursively). Their (chr, start, end) are used even with --reference,
                            so omitted chromosomes and --no-nan-windows gaps match the stats.
  --template-tsv FILE       Single template TSV

Mask (same files as calculate_pi_theta.sh / calculate_fst.sh; sites inside are not counted):
  --filter-mask-total-bed FILE    Repeat / combined mask BED
  --filter-mask-total-fasta FILE  vcftools-style mask FASTA (use one of BED or FASTA)

  --no-nan-windows            When building from --reference, omit windows with n_sites=0
                              (default: on; same idea as grenedalf --no-nan-windows).
                              Template windows are kept as listed in the stats tables.
  --keep-nan-windows          Write empty reference-grid windows (turns off --no-nan-windows)

Optional:
  --output-dir DIR          Where to write variant_features_wW_sS.tsv [default: .]
  --rscript PATH            [default: Rscript]
  --dry-run
  -h, --help

Example:
  $0 --sites-tsv csq_snps/variants_csq.tsv --sites-tsv csq_indels/variants_csq.tsv \\
     --reference ref.unmasked.fa -w 1000 --step-size 500 --template-dir pi_out \\
     --filter-mask-total-bed mask/mask.bed --output-dir vf_out
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
        --filter-mask-total-bed) FILTER_MASK_TOTAL_BED="$2"; shift 2 ;;
        --filter-mask-total-fasta) FILTER_MASK_TOTAL_FASTA="$2"; shift 2 ;;
        --no-nan-windows|--non-nan-windows) NO_NAN_WINDOWS=true; shift ;;
        --keep-nan-windows) NO_NAN_WINDOWS=false; shift ;;
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
if [[ -n "$FILTER_MASK_TOTAL_BED" && -n "$FILTER_MASK_TOTAL_FASTA" ]]; then
    log_error "Use only one of --filter-mask-total-bed and --filter-mask-total-fasta"
    exit 1
fi
if [[ -n "$FILTER_MASK_TOTAL_BED" && ! -f "$FILTER_MASK_TOTAL_BED" ]]; then
    log_error "Mask BED not found: $FILTER_MASK_TOTAL_BED"
    exit 1
fi
if [[ -n "$FILTER_MASK_TOTAL_FASTA" && ! -f "$FILTER_MASK_TOTAL_FASTA" ]]; then
    log_error "Mask FASTA not found: $FILTER_MASK_TOTAL_FASTA"
    exit 1
fi
if [[ -n "$TEMPLATE_DIR" && ! -d "$TEMPLATE_DIR" ]]; then
    log_error "Template directory not found: $TEMPLATE_DIR"
    exit 1
fi
if [[ -n "$TEMPLATE_TSV" && ! -f "$TEMPLATE_TSV" ]]; then
    log_error "Template TSV not found: $TEMPLATE_TSV"
    exit 1
fi

ARGS=(--sites-tsv "$SITES_TSV" --output-dir "$OUTPUT_DIR")
[[ -n "$REFERENCE" ]] && ARGS+=(--reference "$REFERENCE")
[[ ${#WINDOW_SIZES[@]} -gt 0 ]] && ARGS+=(--window-size "$(IFS=,; echo "${WINDOW_SIZES[*]}")")
[[ ${#STEP_SIZES[@]} -gt 0 ]] && ARGS+=(--step-size "$(IFS=,; echo "${STEP_SIZES[*]}")")
[[ -n "$TEMPLATE_DIR" ]] && ARGS+=(--template-dir "$TEMPLATE_DIR")
[[ -n "$TEMPLATE_TSV" ]] && ARGS+=(--template-tsv "$TEMPLATE_TSV")
[[ -n "$FILTER_MASK_TOTAL_BED" ]] && ARGS+=(--filter-mask-total-bed "$FILTER_MASK_TOTAL_BED")
[[ -n "$FILTER_MASK_TOTAL_FASTA" ]] && ARGS+=(--filter-mask-total-fasta "$FILTER_MASK_TOTAL_FASTA")
if [[ "$NO_NAN_WINDOWS" == true ]]; then
    ARGS+=(--no-nan-windows)
else
    ARGS+=(--keep-nan-windows)
fi

if [[ "$DRY_RUN" == true ]]; then
    echo "Would run: $RSCRIPT $WIN_R ${ARGS[*]}"
    exit 0
fi

log "Running variant_csq_windows.R"
"$RSCRIPT" "$WIN_R" "${ARGS[@]}"
log "Done."
