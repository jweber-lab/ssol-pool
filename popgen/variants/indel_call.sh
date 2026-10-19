#!/bin/bash

###############################################################################
# indel_call.sh
#
# Indel site discovery with freebayes --pooled-continuous (no ploidy model;
# alleles are called from read counts), filtered by per-pool read counts the
# same way snp_call.sh filters SNPs, so the two outputs line up.
# Pipeline: freebayes (indels only, parallel over genome chunks) ->
# bcftools norm (split multiallelic, left-align) -> count filters ->
# <prefix>.tsv and <prefix>.bcf (FORMAT/AD, FORMAT/DP; no GT).
#
# freebayes QUAL is not informative for low-frequency pooled alleles (real 2-10%
# indels score ~0), so sites are filtered on counts and frequencies, not QUAL.
#
# Each BAM must carry exactly one @RG SM (process_poolseq.sh sets it); freebayes
# assigns reads to pools by SM. With --sample-info, output columns use sample_name.
#
# Author: ssol-poolseq
# Usage: See below or run with --help
###############################################################################

set -euo pipefail

BAM_FILES=()
SAMPLE_INFO_CSV=""
REFERENCE_FASTA=""
OUTPUT_DIR=""
PREFIX="indel_sites"
FREEBAYES=""
BCFTOOLS=""
THREADS=1
CHUNK_SIZE=1000000
MIN_COVERAGE=10
MAX_COVERAGE=1000
MIN_COUNT=2
MIN_FREQUENCY=0.01
MIN_COUNT_FILTER=0
MIN_TOTAL_READ_DEPTH=10
MIN_MAPQ=20
MIN_BASEQ=0
FILTER_MASK_TOTAL_BED=""
KEEP_RAW=false
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Indel site discovery with freebayes --pooled-continuous:
BAM(s) -> freebayes (indels) -> bcftools norm -> count filters -> <prefix>.tsv, <prefix>.bcf

Input (choose one):
  -i, --sample-info FILE      CSV with columns: sample_name, read1_file, read2_file, pool_size, bam_file
                              (same file as calculate_fst.sh / snp_call.sh)
  -b, --bam FILE              BAM file(s) - can be specified multiple times

Required:
  -f, --reference FILE        Reference FASTA (indexed: .fai); same unmasked assembly as alignment
  -o, --output-dir DIR        Output directory

Filters (defaults as snp_call.sh / calculate_fst.sh):
  --min-coverage N            Minimum read depth per pool; shallower pools are set missing (default: 10)
  --max-coverage N            Maximum read depth per pool; deeper pools are set missing (default: 1000).
                              freebayes also skips positions above N x number of pools (collapsed repeats)
  --min-count N               Minimum reads for an allele in a pool; lower counts are set to 0 (default: 2)
  --min-frequency F           Minimum total minor allele frequency (default: 0.01)
  --min-count-filter N        Minimum total minor allele count instead of --min-frequency (default: 0 = off)
  --min-total-read-depth N    Minimum read depth summed over non-missing pools (default: 10)
  --filter-mask-total-bed FILE    Drop indels overlapping these intervals (e.g. RepeatMasker *.out.bed)
  --min-mapq N                Minimum read mapping quality (default: 20)
  --min-baseq N               Minimum base quality for freebayes (default: 0, freebayes default)

freebayes:
  -t, --threads N             Genome chunks run in parallel (default: 1)
  --chunk-size N              Chunk size in bp (default: 1000000)

Optional:
  --prefix NAME               Output file prefix (default: indel_sites; keep "sites", "calls" or
                              "variant" in it so collate.sh --variant-tsv-dir finds the TSV)
  --keep-raw                  Keep the unfiltered freebayes VCF as <prefix>.raw.vcf.gz
  --freebayes PATH            Path to freebayes (default: freebayes in PATH)
  --bcftools PATH             Path to bcftools (default: bcftools in PATH)
  --dry-run                   Preview commands without executing
  -h, --help                  Show this help

Outputs (in --output-dir):
  <prefix>.tsv    chr, pos, ref, alt, dp, alt_freq, indel_len, then <sample>_AD (ref,alt) and
                  <sample>_alt_freq (same layout as snp_call.sh plus indel_len). alt_freq = ALT reads /
                  read depth, so split alleles of a multiallelic site are each relative to all reads
  <prefix>.bcf    Indels with INFO/DP, INFO/AF, FORMAT/AD, FORMAT/DP per pool (no GT); input for
                  variant_csq.sh

Example:
  $0 --sample-info ../../sample_info.csv --reference ref.unmasked.fa \\
     --filter-mask-total-bed genome.out.bed --output-dir indels_out --threads 8
EOF
}

log() { echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_dry_run() { echo -e "${YELLOW}[DRY-RUN]${NC} $1"; }

check_command() {
    command -v "$1" &>/dev/null || { log_error "$1 not found"; exit 1; }
}

resolve_path_relative_to_csv() {
    local csv_file="$1"
    local file_path="$2"
    if [[ "$file_path" == /* ]]; then
        echo "$file_path"
    else
        echo "$(cd "$(dirname "$csv_file")" && pwd)/${file_path}"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--sample-info) SAMPLE_INFO_CSV="$2"; shift 2 ;;
        -b|--bam) BAM_FILES+=("$2"); shift 2 ;;
        -f|--reference) REFERENCE_FASTA="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --freebayes) FREEBAYES="$2"; shift 2 ;;
        --bcftools) BCFTOOLS="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --chunk-size) CHUNK_SIZE="$2"; shift 2 ;;
        --min-coverage) MIN_COVERAGE="$2"; shift 2 ;;
        --max-coverage) MAX_COVERAGE="$2"; shift 2 ;;
        --min-count) MIN_COUNT="$2"; shift 2 ;;
        --min-frequency) MIN_FREQUENCY="$2"; shift 2 ;;
        --min-count-filter) MIN_COUNT_FILTER="$2"; shift 2 ;;
        --min-total-read-depth) MIN_TOTAL_READ_DEPTH="$2"; shift 2 ;;
        --filter-mask-total-bed) FILTER_MASK_TOTAL_BED="$2"; shift 2 ;;
        --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
        --min-baseq) MIN_BASEQ="$2"; shift 2 ;;
        --keep-raw) KEEP_RAW=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -n "$SAMPLE_INFO_CSV" && ${#BAM_FILES[@]} -gt 0 ]] || [[ -z "$SAMPLE_INFO_CSV" && ${#BAM_FILES[@]} -eq 0 ]]; then
    log_error "Provide exactly one input: --sample-info or --bam"
    usage
    exit 1
fi
if [[ -z "$REFERENCE_FASTA" || -z "$OUTPUT_DIR" ]]; then
    log_error "Missing --reference or --output-dir"
    usage
    exit 1
fi
for f in "$REFERENCE_FASTA" "$FILTER_MASK_TOTAL_BED"; do
    if [[ -n "$f" && ! -f "$f" ]]; then
        log_error "File not found: $f"
        exit 1
    fi
done
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ && "$CHUNK_SIZE" =~ ^[1-9][0-9]*$ ]]; then
    log_error "--threads and --chunk-size must be positive integers (got: $THREADS, $CHUNK_SIZE)"
    exit 1
fi

# Resolve samples from CSV (sample_name in column 1, bam_file in column 5)
SAMPLE_NAMES=()
if [[ -n "$SAMPLE_INFO_CSV" ]]; then
    if [[ ! -f "$SAMPLE_INFO_CSV" ]]; then
        log_error "Sample info file not found: $SAMPLE_INFO_CSV"
        exit 1
    fi
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        IFS=',' read -r -a fields <<< "$(tr -d '"' <<< "$line")"
        case "${fields[0]:-}" in
            [Ss]ample_name|[Ss]ample|[Bb]am_file) continue ;;
        esac
        if [[ ${#fields[@]} -lt 5 || -z "${fields[0]}" || -z "${fields[4]}" ]]; then
            log_warn "Skipping CSV line without sample_name and bam_file: $line"
            continue
        fi
        bam=$(resolve_path_relative_to_csv "$SAMPLE_INFO_CSV" "${fields[4]}")
        if [[ "$DRY_RUN" == false && ! -f "$bam" ]]; then
            log_warn "BAM not found for sample ${fields[0]}: $bam"
            continue
        fi
        SAMPLE_NAMES+=("${fields[0]}")
        BAM_FILES+=("$bam")
    done < "$SAMPLE_INFO_CSV"
    if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
        log_error "No valid samples found in $SAMPLE_INFO_CSV"
        exit 1
    fi
fi

: "${FREEBAYES:=freebayes}"
: "${BCFTOOLS:=bcftools}"
check_command "$FREEBAYES"
check_command "$BCFTOOLS"
check_command samtools

# freebayes groups reads by @RG SM: each BAM needs exactly one SM, distinct across BAMs.
# RENAME maps SM -> output name (CSV sample_name, else SM).
declare -A RENAME=()
declare -A SEEN_SM=()
if [[ "$DRY_RUN" == false ]]; then
    for i in "${!BAM_FILES[@]}"; do
        bam="${BAM_FILES[$i]}"
        [[ -f "$bam" ]] || { log_error "BAM not found: $bam"; exit 1; }
        if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" && ! -f "${bam}.csi" ]]; then
            log "BAM index not found for $bam; running samtools index"
            if ! samtools index -@ "$THREADS" "$bam"; then
                log_error "samtools index failed for $bam (BAM must be coordinate-sorted and its directory writable)"
                exit 1
            fi
        elif [[ -f "${bam}.bai" && "${bam}.bai" -ot "$bam" ]] || [[ -f "${bam}.csi" && "${bam}.csi" -ot "$bam" ]]; then
            log_warn "BAM index for $bam is older than the BAM; re-run samtools index if the BAM was rewritten"
        fi
        mapfile -t sms < <(samtools view -H "$bam" | awk -F'\t' '/^@RG/ {for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) print substr($i, 4)}' | sort -u)
        if [[ ${#sms[@]} -ne 1 ]]; then
            log_error "$bam has ${#sms[@]} distinct @RG SM values (need exactly 1). Set one with: samtools addreplacerg -r 'ID:<name>' -r 'SM:<name>'"
            exit 1
        fi
        sm="${sms[0]}"
        if [[ -n "${SEEN_SM[$sm]:-}" ]]; then
            log_error "@RG SM '$sm' appears in both ${SEEN_SM[$sm]} and $bam; freebayes would merge them into one pool"
            exit 1
        fi
        SEEN_SM[$sm]="$bam"
        RENAME[$sm]="${SAMPLE_NAMES[$i]:-$sm}"
    done
fi

if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
    log "Indexing reference with samtools faidx..."
    [[ "$DRY_RUN" == true ]] && log_dry_run "Would execute: samtools faidx $REFERENCE_FASTA" || samtools faidx "$REFERENCE_FASTA"
fi

OUT_TSV="${OUTPUT_DIR}/${PREFIX}.tsv"
OUT_BCF="${OUTPUT_DIR}/${PREFIX}.bcf"
OUT_RAW="${OUTPUT_DIR}/${PREFIX}.raw.vcf.gz"
WORK_DIR="${OUTPUT_DIR}/.${PREFIX}_work"

SKIP_COVERAGE=$((MAX_COVERAGE * ${#BAM_FILES[@]}))
FB_CMD=(
    "$FREEBAYES" -f "$REFERENCE_FASTA"
    --pooled-continuous
    -I -X
    --min-alternate-fraction "$MIN_FREQUENCY"
    --min-alternate-count "$MIN_COUNT"
    --min-mapping-quality "$MIN_MAPQ"
    --min-base-quality "$MIN_BASEQ"
    --skip-coverage "$SKIP_COVERAGE"
)
# 0-based half-open chunks (freebayes --region uses BED-style coordinates)
mapfile -t REGIONS < <(awk -F'\t' -v size="$CHUNK_SIZE" '{for (s = 0; s < $2; s += size) printf "%s:%d-%d\n", $1, s, (s + size < $2 ? s + size : $2)}' "${REFERENCE_FASTA}.fai")

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "Would run ${#REGIONS[@]} freebayes chunks, $THREADS at a time, e.g.:"
    log_dry_run "  ${FB_CMD[*]} --region ${REGIONS[0]:-<chr:start-end>} ${BAM_FILES[*]}"
    log_dry_run "Would run: bcftools norm -f $REFERENCE_FASTA -m -any | bcftools view -v indels -> count filters -> $OUT_TSV, $OUT_BCF"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
trap 'kill $(jobs -rp) 2>/dev/null || true; rm -rf "$WORK_DIR"' EXIT

log "Running freebayes --pooled-continuous (indels only) on ${#BAM_FILES[@]} pools: ${#REGIONS[@]} chunks, $THREADS in parallel"
log "  ${FB_CMD[*]} --region <chunk> <bams>"

run_chunk() {
    local i="$1"
    "${FB_CMD[@]}" --region "${REGIONS[$i]}" "${BAM_FILES[@]}" > "${WORK_DIR}/chunk_${i}.vcf" 2> "${WORK_DIR}/chunk_${i}.log" \
        || { echo "${REGIONS[$i]}" > "${WORK_DIR}/chunk_${i}.failed"; }
}

for i in "${!REGIONS[@]}"; do
    while [[ $(jobs -rp | wc -l) -ge $THREADS ]]; do
        wait -n || true
    done
    run_chunk "$i" &
done
wait

if compgen -G "${WORK_DIR}/chunk_*.failed" > /dev/null; then
    for f in "${WORK_DIR}"/chunk_*.failed; do
        log_error "freebayes failed on $(cat "$f"):"
        tail -5 "${f%.failed}.log" >&2
    done
    exit 1
fi

# Chunks are in reference order: header from the first, records from all.
RAW_VCF="${WORK_DIR}/raw.vcf"
{
    grep '^#' "${WORK_DIR}/chunk_0.vcf" | sed '/^##FORMAT=<ID=GQ,/s/Type=Float/Type=Integer/'
    for i in "${!REGIONS[@]}"; do
        grep -v '^#' "${WORK_DIR}/chunk_${i}.vcf" || true
    done
} > "$RAW_VCF"
if [[ "$KEEP_RAW" == true ]]; then
    "$BCFTOOLS" view -Oz -o "$OUT_RAW" "$RAW_VCF"
    "$BCFTOOLS" index -f "$OUT_RAW"
fi

# Split multiallelic records (per-allele AO), left-align/trim, drop chunk-boundary duplicates,
# keep indels (optionally outside the mask), then filter on counts.
NORM_VCF="${WORK_DIR}/norm.bcf"
"$BCFTOOLS" norm -f "$REFERENCE_FASTA" -m -any -Ou "$RAW_VCF" 2> "${WORK_DIR}/norm.log" \
    | "$BCFTOOLS" norm -d exact -Ou 2>> "${WORK_DIR}/norm.log" \
    | "$BCFTOOLS" view -v indels -Ob -o "$NORM_VCF"
"$BCFTOOLS" index -f "$NORM_VCF"
if [[ -n "$FILTER_MASK_TOTAL_BED" ]]; then
    "$BCFTOOLS" view -T "^${FILTER_MASK_TOTAL_BED}" --targets-overlap record -Ob -o "${WORK_DIR}/masked.bcf" "$NORM_VCF"
    NORM_VCF="${WORK_DIR}/masked.bcf"
fi

SAMPLES_FILE="${WORK_DIR}/samples.txt"
"$BCFTOOLS" query -l "$NORM_VCF" | while read -r sm; do echo "${RENAME[$sm]:-$sm}"; done > "$SAMPLES_FILE"

VCF_TMP="${WORK_DIR}/sites.vcf"
"$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%RO\t%AO]\n' "$NORM_VCF" \
| awk -F'\t' -v OFS='\t' -v tsv_out="$OUT_TSV" -v vcf_out="$VCF_TMP" -v samples_file="$SAMPLES_FILE" \
    -v fai="${REFERENCE_FASTA}.fai" -v min_cov="$MIN_COVERAGE" -v max_cov="$MAX_COVERAGE" \
    -v min_count="$MIN_COUNT" -v min_freq="$MIN_FREQUENCY" -v min_count_filter="$MIN_COUNT_FILTER" \
    -v min_total="$MIN_TOTAL_READ_DEPTH" '
BEGIN {
    n = 0
    while ((getline s < samples_file) > 0) name[++n] = s
    h = "chr\tpos\tref\talt\tdp\talt_freq\tindel_len"
    for (s = 1; s <= n; s++) h = h "\t" name[s] "_AD\t" name[s] "_alt_freq"
    print h > tsv_out
    print "##fileformat=VCFv4.2" > vcf_out
    print "##source=indel_call.sh (freebayes --pooled-continuous read counts; genotypes not called)" > vcf_out
    while ((getline line < fai) > 0) {
        split(line, c, "\t")
        print "##contig=<ID=" c[1] ",length=" c[2] ">" > vcf_out
    }
    print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth summed over pools passing filters\">" > vcf_out
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT reads / read depth over pools passing filters\">" > vcf_out
    print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read counts for REF and ALT (freebayes RO, AO)\">" > vcf_out
    print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth (freebayes DP)\">" > vcf_out
    h = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (s = 1; s <= n; s++) h = h "\t" name[s]
    print h > vcf_out
}
{
    tr = 0; ta = 0; td = 0; any = 0
    for (s = 1; s <= n; s++) {
        dp = $(3 * s + 2); ro = $(3 * s + 3); ao = $(3 * s + 4)
        miss[s] = (dp == "." || ro == "." || ao == "." || dp + 0 < min_cov || dp + 0 > max_cov)
        if (miss[s]) continue
        ro += 0; ao += 0; dp += 0
        d[s] = (dp > ro + ao ? dp : ro + ao)
        if (ro < min_count) ro = 0
        if (ao < min_count) ao = 0
        r[s] = ro; a[s] = ao
        tr += ro; ta += ao; td += d[s]; any = 1
    }
    if (!any || td < min_total || tr == 0 || ta == 0) next
    af = ta / td
    if (min_count_filter > 0) {
        if ((ta < td - ta ? ta : td - ta) < min_count_filter) next
    } else if ((af < 1 - af ? af : 1 - af) < min_freq) next
    tsv = $1 OFS $2 OFS $3 OFS $4 OFS td OFS sprintf("%.6g", af) OFS (length($4) - length($3))
    vcf = $1 OFS $2 OFS "." OFS $3 OFS $4 OFS "." OFS "PASS" OFS "DP=" td ";AF=" sprintf("%.6g", af) OFS "AD:DP"
    for (s = 1; s <= n; s++) {
        if (miss[s]) { tsv = tsv OFS "." OFS "NA"; vcf = vcf OFS ".:."; continue }
        tsv = tsv OFS r[s] "," a[s] OFS (d[s] > 0 ? sprintf("%.6g", a[s] / d[s]) : "NA")
        vcf = vcf OFS r[s] "," a[s] ":" d[s]
    }
    print tsv > tsv_out
    print vcf > vcf_out
    kept++
}
END { print kept + 0 > "/dev/stderr" }' 2> "${WORK_DIR}/kept.txt"

log "Writing BCF: $OUT_BCF"
"$BCFTOOLS" view -Ob -o "$OUT_BCF" "$VCF_TMP"
"$BCFTOOLS" index -f "$OUT_BCF"

log "Done. $(cat "${WORK_DIR}/kept.txt") indel sites."
log "  TSV: $OUT_TSV (chr, pos, ref, alt, dp, alt_freq, indel_len, <sample>_AD, <sample>_alt_freq)"
log "  BCF: $OUT_BCF (FORMAT/AD, FORMAT/DP; no GT)"
if [[ "$KEEP_RAW" == true ]]; then
    log "  Raw freebayes VCF: $OUT_RAW"
fi
log "Next step: variant_csq.sh --bcf $OUT_BCF --reference <ref> --gff <gff3> --output-dir <dir>"
