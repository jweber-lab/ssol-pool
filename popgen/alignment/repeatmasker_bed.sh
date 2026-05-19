#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# repeatmasker_bed.sh
#
# Convert RepeatMasker .out (whitespace-delimited, not .align) to 3-column BED
# for grenedalf --filter-mask-total-bed and variant_call --exclude-regions-bed.
###############################################################################

MERGE=1                  # 1 = merge overlapping intervals (requires bedtools)
SORT=1                   # 1 = sort BED (always sorted before merge; see below)
MIN_LEN=1                # minimum interval length to keep
STRIP_PAREN=1            # remove parentheses from coordinates if present
CHR_PREFIX=""            # set to "chr" to force prefix, "" to leave unchanged
FILTER_CLASS=""          # awk regex on column 11 (e.g. "LINE" matches LINE/L1)

usage() {
    cat <<EOF
Usage: $0 -i repeatmasker.out -o repeats.bed [options]

Convert RepeatMasker .out to BED (0-based, half-open). Merged output is suitable
for --filter-mask-total-bed and --exclude-regions-bed in popgen/stats/.

Options:
  -i FILE               RepeatMasker .out input (required)
  -o FILE               Output BED path (required; parent directory created if needed)
  --no-merge            Do not merge intervals (bedtools not required)
  --no-sort             Do not sort when --no-merge is set (ignored if merging)
  --min-len INT         Minimum interval length in bp (default: 1)
  --chr-prefix STR      Add prefix to chromosome names if missing (e.g. "chr")
  --filter-class STR    Keep repeats whose class/family (col 11) matches this regex
  -h, --help            Show this help

Notes:
  - Merge uses bedtools merge and requires sorted input; sorting is always applied
    before merge even if --no-sort is set (a warning is printed).
  - --filter-class uses awk regex match (e.g. "^LINE/" for LINE only).
EOF
}

IN=""
OUT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) IN="$2"; shift 2 ;;
        -o) OUT="$2"; shift 2 ;;
        --no-merge) MERGE=0; shift ;;
        --no-sort) SORT=0; shift ;;
        --min-len) MIN_LEN="$2"; shift 2 ;;
        --chr-prefix) CHR_PREFIX="$2"; shift 2 ;;
        --filter-class) FILTER_CLASS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$IN" || -z "$OUT" ]]; then
    echo "ERROR: -i and -o are required" >&2
    usage
    exit 1
fi

if [[ ! -f "$IN" ]]; then
    echo "ERROR: Input file not found: $IN" >&2
    exit 1
fi

if [[ "$MIN_LEN" -lt 1 ]] 2>/dev/null; then
    echo "ERROR: --min-len must be a positive integer" >&2
    exit 1
fi

OUT_DIR=$(dirname "$OUT")
if [[ "$OUT_DIR" != "." ]] && [[ -n "$OUT_DIR" ]]; then
    mkdir -p "$OUT_DIR"
fi

echo "[repeatmasker_bed] input=$IN output=$OUT merge=$MERGE sort=$SORT min_len=$MIN_LEN chr_prefix=${CHR_PREFIX:-<none>} filter_class=${FILTER_CLASS:-<none>}"

TMP=$(mktemp)
trap 'rm -f "$TMP" "${TMP}.sorted"' EXIT

# RepeatMasker .out columns (typical):
#  5 = sequence (chrom), 6 = start (1-based), 7 = end (inclusive), 11 = class/family
awk -v min_len="$MIN_LEN" \
    -v strip_paren="$STRIP_PAREN" \
    -v chr_prefix="$CHR_PREFIX" \
    -v filter_class="$FILTER_CLASS" '
BEGIN { OFS="\t" }

NR <= 3 { next }

{
    chrom = $5
    start = $6
    end   = $7
    class = $11

    if (strip_paren) {
        gsub(/[()]/, "", start)
        gsub(/[()]/, "", end)
    }

    if (start !~ /^[0-9]+$/ || end !~ /^[0-9]+$/) next

    start = start + 0
    end   = end + 0

    if (end < start) next

    bed_start = start - 1
    bed_end   = end

    if ((bed_end - bed_start) < min_len) next

    if (filter_class != "" && class !~ filter_class) next

    if (chr_prefix != "") {
        if (chrom !~ "^" chr_prefix) {
            chrom = chr_prefix chrom
        }
    }

    print chrom, bed_start, bed_end
}
' "$IN" > "$TMP"

do_sort=false
if [[ "$MERGE" -eq 1 ]]; then
    if [[ "$SORT" -eq 0 ]]; then
        echo "[repeatmasker_bed] WARNING: --no-sort ignored; bedtools merge requires sorted BED" >&2
    fi
    do_sort=true
elif [[ "$SORT" -eq 1 ]]; then
    do_sort=true
fi

if [[ "$do_sort" == true ]]; then
    sort -k1,1 -k2,2n "$TMP" > "${TMP}.sorted"
    mv "${TMP}.sorted" "$TMP"
fi

if [[ "$MERGE" -eq 1 ]]; then
    if ! command -v bedtools >/dev/null 2>&1; then
        echo "ERROR: bedtools required for merging (install bedtools or use --no-merge)" >&2
        exit 1
    fi
    bedtools merge -i "$TMP" > "$OUT"
else
    cp "$TMP" "$OUT"
fi

echo "Output written to: $OUT"
echo "Intervals: $(wc -l < "$OUT" | tr -d ' ')"
