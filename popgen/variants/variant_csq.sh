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
  -b, --bcf FILE              Input BCF/VCF (indexed recommended). Without FORMAT/GT (e.g. from
                              snp_call.sh), csq runs with -s - -l (per-site, samples ignored)
  -f, --reference FILE        Reference FASTA for codon translation (faidx indexed).
                              Use the **unmasked** assembly that matches your GFF coordinates
                              (not the hard-masked genome used only for alignment).
  -g, --gff FILE              GFF3 for bcftools csq (see bcftools howto). Sites on chromosomes
                              without GFF features (gene-free scaffolds) are kept without BCSQ
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

# csq skips records without genotypes; pool-seq sites (snp_call.sh) have AD but no GT, so
# ignore samples and annotate each record on its own (allele phase across sites is unknown).
if command -v "$BCFTOOLS" &>/dev/null && ! "$BCFTOOLS" view -h "$BCF_FILE" | grep -q '^##FORMAT=<ID=GT,'; then
    log "No FORMAT/GT in $BCF_FILE: annotating sites without genotypes (bcftools csq -s - -l)"
    CSQ_EXTRA=(-s - -l "${CSQ_EXTRA[@]}")
fi

if [[ "$DRY_RUN" == true ]]; then
    echo "Would run: $BCFTOOLS csq -f $REFERENCE -g $GFF ${CSQ_EXTRA[*]} --threads $THREADS -Ob -o $OUTPUT_BCF $BCF_FILE"
    echo "Would index: $OUTPUT_BCF"
    echo "Would query -> $QUERY_TMP"
    echo "Would run: $RSCRIPT $VARIANT_CSQ_R --input-tsv $QUERY_TMP --output-tsv $OUTPUT_TSV ..."
    exit 0
fi

command -v "$BCFTOOLS" &>/dev/null || { log_error "bcftools not found: $BCFTOOLS"; exit 1; }

# csq aborts on records whose chromosome has no GFF features (e.g. gene-free scaffolds). --force
# would also hide REF/fasta mismatches, so instead csq runs only on chromosomes present in the GFF
# and its BCSQ is transferred back onto the full input; other sites are kept without BCSQ.
CONTIG_TMP="$(dirname "$OUTPUT_BCF")/.csq_contigs_$$"
mkdir -p "$CONTIG_TMP"
trap 'rm -rf "$CONTIG_TMP"' EXIT
awk -F'\t' '!/^#/ && NF >= 9 { print $1 }' "$GFF" | sort -u > "${CONTIG_TMP}/gff"
CSQ_INPUT="$BCF_FILE"
if ! "$BCFTOOLS" index -s "$CSQ_INPUT" > "${CONTIG_TMP}/vcf.stats" 2>/dev/null; then
    log "No index for $BCF_FILE: writing an indexed temporary BCF"
    CSQ_INPUT="${CONTIG_TMP}/input.bcf"
    "$BCFTOOLS" view -Ob --threads "$THREADS" -o "$CSQ_INPUT" "$BCF_FILE"
    "$BCFTOOLS" index -f "$CSQ_INPUT"
    "$BCFTOOLS" index -s "$CSQ_INPUT" > "${CONTIG_TMP}/vcf.stats"
fi
cut -f1 "${CONTIG_TMP}/vcf.stats" | sort -u > "${CONTIG_TMP}/vcf"
n_vcf=$(wc -l < "${CONTIG_TMP}/vcf")
n_missing=$(comm -23 "${CONTIG_TMP}/vcf" "${CONTIG_TMP}/gff" | wc -l)
n_shared=$(comm -12 "${CONTIG_TMP}/vcf" "${CONTIG_TMP}/gff" | wc -l)

if [[ "$n_shared" -eq 0 ]]; then
    log_error "No chromosome names shared between $BCF_FILE and $GFF (e.g. VCF: $(head -1 "${CONTIG_TMP}/vcf"), GFF: $(head -1 "${CONTIG_TMP}/gff")); check naming or use --csq-arg --unify-chr-names"
    exit 1
elif [[ "$n_missing" -eq 0 ]]; then
    log "Running bcftools csq -> $OUTPUT_BCF"
    "$BCFTOOLS" csq -f "$REFERENCE" -g "$GFF" "${CSQ_EXTRA[@]}" --threads "$THREADS" -Ob -o "$OUTPUT_BCF" "$CSQ_INPUT"
else
    log "${n_missing} of ${n_vcf} chromosome(s) in $BCF_FILE have no GFF features (e.g. $(comm -23 "${CONTIG_TMP}/vcf" "${CONTIG_TMP}/gff" | head -3 | paste -sd, -)); their sites are kept without BCSQ"
    comm -12 "${CONTIG_TMP}/vcf" "${CONTIG_TMP}/gff" | awk -v OFS='\t' '{ print $1, 1, 2147483647 }' > "${CONTIG_TMP}/targets.tsv"
    log "Running bcftools csq on ${n_shared} chromosome(s) with GFF features"
    "$BCFTOOLS" view -R "${CONTIG_TMP}/targets.tsv" -Ou "$CSQ_INPUT" \
        | "$BCFTOOLS" csq -f "$REFERENCE" -g "$GFF" "${CSQ_EXTRA[@]}" --threads "$THREADS" -Ob -o "${CONTIG_TMP}/genic.bcf" -
    "$BCFTOOLS" index -f "${CONTIG_TMP}/genic.bcf"
    log "Transferring BCSQ onto all sites -> $OUTPUT_BCF"
    "$BCFTOOLS" annotate -a "${CONTIG_TMP}/genic.bcf" -c INFO/BCSQ --threads "$THREADS" -Ob -o "$OUTPUT_BCF" "$CSQ_INPUT"
fi
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
