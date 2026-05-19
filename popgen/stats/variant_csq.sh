#!/bin/bash

###############################################################################
# variant_csq.sh
#
# Run bcftools csq on a BCF/VCF, export a query TSV, then variant_csq.R for
# BCSQ parsing and derived region_category / coding_effect columns.
#
# Sibling: variant_csq.R
###############################################################################

set -euo pipefail

BCF_FILE=""
REFERENCE=""
GFF=""
OUTPUT_BCF=""
OUTPUT_TSV=""
OUTPUT_LONG=""
MAP_TSV=""
BCFTOOLS="bcftools"
RSCRIPT="Rscript"
THREADS=1
CSQ_EXTRA=()   # extra args to bcftools csq (e.g. --nc)
DRY_RUN=false
OUTPUT_DIR=""

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log() { echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  -b, --bcf FILE              Input BCF/VCF (indexed recommended)
  -f, --reference FILE        Reference FASTA for codon translation (faidx indexed).
                              Use the **unmasked** assembly that matches your GFF coordinates
                              (not the hard-masked genome used only for alignment).
  -g, --gff FILE              GFF3 for bcftools csq (see bcftools howto)
  -o, --output-bcf FILE       Output BCF path with BCSQ annotations
  OR --output-dir DIR         Output directory; writes calls.csq.bcf and annotated TSV

Outputs:
  --output-tsv FILE           Final annotated site TSV (default: <dir>/variants_csq.tsv with --output-dir)
  --output-long-tsv FILE      Optional long-format BCSQ rows (one row per sub-record)

Optional:
  --map-tsv FILE              Optional mapping table for variant_csq.R (see variant_csq.R --help)
  --bcftools PATH             bcftools binary [default: bcftools]
  --rscript PATH              Rscript [default: Rscript]
  -t, --threads N             Threads for bcftools csq [default: 1]
  --csq-arg ARG               Pass-through to bcftools csq (repeatable)
  --dry-run
  -h, --help

Example:
  $0 --bcf calls.bcf --reference ref.fa --gff genes.gff3 --output-dir csq_out
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bcf) BCF_FILE="$2"; shift 2 ;;
        -f|--reference) REFERENCE="$2"; shift 2 ;;
        -g|--gff) GFF="$2"; shift 2 ;;
        -o|--output-bcf) OUTPUT_BCF="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --output-tsv) OUTPUT_TSV="$2"; shift 2 ;;
        --output-long-tsv) OUTPUT_LONG="$2"; shift 2 ;;
        --map-tsv) MAP_TSV="$2"; shift 2 ;;
        --bcftools) BCFTOOLS="$2"; shift 2 ;;
        --rscript) RSCRIPT="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --csq-arg) CSQ_EXTRA+=("$2"); shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

: "${OUTPUT_DIR:=}"

if [[ -z "$BCF_FILE" || -z "$REFERENCE" || -z "$GFF" ]]; then
    log_error "Missing --bcf, --reference, or --gff"
    usage
    exit 1
fi

if [[ ! -f "$BCF_FILE" ]]; then
    log_error "BCF not found: $BCF_FILE"
    exit 1
fi
if [[ ! -f "$REFERENCE" ]]; then
    log_error "Reference not found: $REFERENCE"
    exit 1
fi
if [[ ! -f "$GFF" ]]; then
    log_error "GFF not found: $GFF"
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
VARIANT_CSQ_R="${SCRIPT_DIR}/variant_csq.R"
if [[ ! -f "$VARIANT_CSQ_R" ]]; then
    log_error "variant_csq.R not found: $VARIANT_CSQ_R"
    exit 1
fi

if [[ -n "$OUTPUT_DIR" && -z "$OUTPUT_BCF" ]]; then
    OUTPUT_BCF="${OUTPUT_DIR}/calls.csq.bcf"
fi
if [[ -z "$OUTPUT_BCF" ]]; then
    log_error "Specify -o/--output-bcf or --output-dir"
    usage
    exit 1
fi

if [[ -n "$OUTPUT_DIR" && -z "$OUTPUT_TSV" ]]; then
    OUTPUT_TSV="${OUTPUT_DIR}/variants_csq.tsv"
elif [[ -z "$OUTPUT_TSV" ]]; then
    base="${OUTPUT_BCF%.bcf}"
    base="${base%.vcf.gz}"
    base="${base%.vcf}"
    OUTPUT_TSV="${base}_annotated.tsv"
fi

mkdir -p "$(dirname "$OUTPUT_BCF")"
mkdir -p "$(dirname "$OUTPUT_TSV")"

QUERY_TMP="$(dirname "$OUTPUT_TSV")/.query_csq_$$.tsv"

if [[ "$DRY_RUN" == true ]]; then
    echo "Would run: $BCFTOOLS csq -f $REFERENCE -g $GFF ${CSQ_EXTRA[*]} --threads $THREADS -Ob -o $OUTPUT_BCF $BCF_FILE"
    echo "Would index: $OUTPUT_BCF"
    echo "Would query -> $QUERY_TMP"
    echo "Would run: $RSCRIPT $VARIANT_CSQ_R --input-tsv $QUERY_TMP --output-tsv $OUTPUT_TSV ..."
    exit 0
fi

command -v "$BCFTOOLS" &>/dev/null || { log_error "bcftools not found: $BCFTOOLS"; exit 1; }

log "Running bcftools csq -> $OUTPUT_BCF"
"$BCFTOOLS" csq -f "$REFERENCE" -g "$GFF" "${CSQ_EXTRA[@]}" --threads "$THREADS" -Ob -o "$OUTPUT_BCF" "$BCF_FILE"
"$BCFTOOLS" index -f "$OUTPUT_BCF"

log "bcftools query -> $QUERY_TMP"
# BCSQ may be absent on some rows; QUAL/FILTER/DP for context
{ printf 'chr\tpos\tref\talt\tqual\tfilter\tbcsq\tdp\n'
  "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/BCSQ\t%INFO/DP\n' "$OUTPUT_BCF"
} > "$QUERY_TMP"

R_ARGS=(--input-tsv "$QUERY_TMP" --output-tsv "$OUTPUT_TSV")
[[ -n "${OUTPUT_LONG:-}" ]] && R_ARGS+=(--output-long-tsv "$OUTPUT_LONG")
[[ -n "${MAP_TSV:-}" ]] && R_ARGS+=(--map-tsv "$MAP_TSV")

log "Running variant_csq.R -> $OUTPUT_TSV"
"$RSCRIPT" "$VARIANT_CSQ_R" "${R_ARGS[@]}"

rm -f "$QUERY_TMP"
log "Done. CSQ BCF: $OUTPUT_BCF"
log "Annotated TSV: $OUTPUT_TSV"
