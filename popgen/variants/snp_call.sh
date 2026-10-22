#!/bin/bash

###############################################################################
# snp_call.sh
#
# SNP site discovery from grenedalf allele counts (no genotype calling).
# Pipeline: grenedalf sync (same filters as calculate_fst.sh) -> SNP sites as
# sync, TSV (per-pool ref,alt counts and alt frequency), and BCF (FORMAT/AD).
#
# Unlike genotype calling (bcftools call supports ploidy 1 or 2 only), this keeps every
# position that passes grenedalf's read-count filters, including low-frequency
# alleles, so the site list matches the SNPs used for FST/diversity.
#
# Sites are polymorphic across the pools (grenedalf --filter-total-only-snps):
# positions where every pool is fixed for the same non-reference base are not
# included. ALT is the most frequent non-reference base in total; any third
# allele is ignored in the TSV/BCF (still visible in the sync output).
# GT is not written; the BCF carries FORMAT/AD and FORMAT/DP for each pool.
#
# grenedalf writes every covered position (filtered ones as missing), so its
# output is streamed through a FIFO and only passing positions are kept.
#
# Author: ssol-poolseq
# Usage: See below or run with --help
###############################################################################

set -euo pipefail

BAM_FILES=()
SYNC_FILE=""
SAMPLE_INFO_CSV=""
REFERENCE_FASTA=""
OUTPUT_DIR=""
PREFIX="snp_sites"
GRENEDALF=""
BCFTOOLS=""
THREADS=1
MIN_COVERAGE=10
MAX_COVERAGE=1000
MIN_COUNT=2
MIN_FREQUENCY=0.01
MIN_COUNT_FILTER=0
MIN_TOTAL_READ_DEPTH=10
MIN_MAPQ=0
MIN_BASEQ=0
BIALLELIC_ONLY=false
FILTER_MASK_TOTAL_BED=""
FILTER_MASK_TOTAL_FASTA=""
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

SNP site discovery from grenedalf allele counts:
BAM(s) or sync -> grenedalf sync (filtered) -> <prefix>.sync, <prefix>.tsv, <prefix>.bcf

Input (choose one):
  -i, --sample-info FILE      CSV with columns: sample_name, read1_file, read2_file, pool_size, bam_file
                              (same file as calculate_fst.sh; sample_name is used for output columns)
  -b, --bam FILE              BAM file(s) - can be specified multiple times (names from BAM basename)
  -s, --sync FILE             Existing multi-sample sync file (filters are applied again)

Required:
  -f, --reference FILE        Reference FASTA (indexed: .fai); same unmasked assembly as alignment
  -o, --output-dir DIR        Output directory

Filters (defaults match calculate_fst.sh, so sites match the SNPs used for FST):
  --min-coverage N            Minimum read depth per sample (default: 10)
  --max-coverage N            Maximum read depth per sample (default: 1000)
  --min-count N               Minimum count for a base to count as an allele, per sample (default: 2)
  --min-frequency F           Minimum total minor allele frequency (default: 0.01)
  --min-count-filter N        Minimum total alt count instead of --min-frequency (default: 0 = use frequency)
  --min-total-read-depth N    Minimum read depth summed across samples (default: 10)
  --biallelic-only            Keep only positions with exactly two alleles in total
  --filter-mask-total-bed FILE    Mask intervals in BED (e.g. RepeatMasker *.out.bed)
  --filter-mask-total-fasta FILE  Mask as vcftools-style FASTA (use one of BED or FASTA)
  --min-mapq N                Minimum read mapping quality, BAM input only (default: 0, as calculate_fst.sh)
  --min-baseq N               Minimum base quality, BAM input only (default: 0, as calculate_fst.sh)

Optional:
  --prefix NAME               Output file prefix (default: snp_sites; keep "sites", "calls" or
                              "variant" in it so collate.sh --variant-tsv-dir finds the TSV)
  -g, --grenedalf PATH        Path to grenedalf (default: grenedalf in PATH)
  --bcftools PATH             Path to bcftools (default: bcftools in PATH)
  -t, --threads N             Threads for grenedalf (default: 1)
  --dry-run                   Preview commands without executing
  -h, --help                  Show this help

Outputs (in --output-dir):
  <prefix>.sync   Passing positions only (grenedalf sync; failing samples marked .:.:.:.:.:.)
  <prefix>.tsv    chr, pos, ref, alt, dp, alt_freq, then <sample>_AD (ref,alt) and <sample>_alt_freq
  <prefix>.bcf    Sites with FORMAT/AD and FORMAT/DP per pool (no GT); input for variant_csq.sh

Example:
  $0 --sample-info ../../sample_info.csv --reference ref.unmasked.fa \\
     --filter-mask-total-bed genome.out.bed --output-dir snps_out --threads 8
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
        -s|--sync) SYNC_FILE="$2"; shift 2 ;;
        -f|--reference) REFERENCE_FASTA="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        -g|--grenedalf) GRENEDALF="$2"; shift 2 ;;
        --bcftools) BCFTOOLS="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --min-coverage) MIN_COVERAGE="$2"; shift 2 ;;
        --max-coverage) MAX_COVERAGE="$2"; shift 2 ;;
        --min-count) MIN_COUNT="$2"; shift 2 ;;
        --min-frequency) MIN_FREQUENCY="$2"; shift 2 ;;
        --min-count-filter) MIN_COUNT_FILTER="$2"; shift 2 ;;
        --min-total-read-depth) MIN_TOTAL_READ_DEPTH="$2"; shift 2 ;;
        --biallelic-only) BIALLELIC_ONLY=true; shift ;;
        --filter-mask-total-bed) FILTER_MASK_TOTAL_BED="$2"; shift 2 ;;
        --filter-mask-total-fasta) FILTER_MASK_TOTAL_FASTA="$2"; shift 2 ;;
        --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
        --min-baseq) MIN_BASEQ="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

# Validate arguments
n_inputs=0
[[ -n "$SAMPLE_INFO_CSV" ]] && n_inputs=$((n_inputs + 1))
[[ ${#BAM_FILES[@]} -gt 0 ]] && n_inputs=$((n_inputs + 1))
[[ -n "$SYNC_FILE" ]] && n_inputs=$((n_inputs + 1))
if [[ $n_inputs -ne 1 ]]; then
    log_error "Provide exactly one input: --sample-info, --bam, or --sync"
    usage
    exit 1
fi
if [[ -z "$REFERENCE_FASTA" || -z "$OUTPUT_DIR" ]]; then
    log_error "Missing --reference or --output-dir"
    usage
    exit 1
fi
if [[ -n "$FILTER_MASK_TOTAL_BED" && -n "$FILTER_MASK_TOTAL_FASTA" ]]; then
    log_error "Use only one of --filter-mask-total-bed and --filter-mask-total-fasta"
    exit 1
fi
for f in "$REFERENCE_FASTA" "$FILTER_MASK_TOTAL_BED" "$FILTER_MASK_TOTAL_FASTA" "$SYNC_FILE"; do
    if [[ -n "$f" && ! -f "$f" ]]; then
        log_error "File not found: $f"
        exit 1
    fi
done

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
if [[ "$DRY_RUN" == false ]]; then
    for bam in "${BAM_FILES[@]}"; do
        [[ -f "$bam" ]] || { log_error "BAM not found: $bam"; exit 1; }
    done
fi

: "${GRENEDALF:=grenedalf}"
: "${BCFTOOLS:=bcftools}"
check_command "$GRENEDALF"
check_command "$BCFTOOLS"
check_command samtools

if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
    log "Indexing reference with samtools faidx..."
    [[ "$DRY_RUN" == true ]] && log_dry_run "Would execute: samtools faidx $REFERENCE_FASTA" || samtools faidx "$REFERENCE_FASTA"
fi

mkdir -p "$OUTPUT_DIR"
OUT_SYNC="${OUTPUT_DIR}/${PREFIX}.sync"
OUT_TSV="${OUTPUT_DIR}/${PREFIX}.tsv"
OUT_BCF="${OUTPUT_DIR}/${PREFIX}.bcf"
WORK_DIR="${OUTPUT_DIR}/.${PREFIX}_work"

# grenedalf 0.6.3 segfaults once the data reach a chromosome that has no interval in a BED
# mask, so every such chromosome gets a 1-bp interval on its last base.
if [[ -n "$FILTER_MASK_TOTAL_BED" && "$DRY_RUN" == false ]]; then
    padded_mask="${OUTPUT_DIR}/.${PREFIX}_mask_all_chromosomes.bed"
    awk -F'\t' -v OFS='\t' '
        NR == FNR { if ($0 !~ /^(#|track|browser)/ && NF >= 3) { seen[$1] = 1; print $1, $2, $3 }; next }
        !($1 in seen) && $2 > 0 { print $1, $2 - 1, $2 }' "$FILTER_MASK_TOTAL_BED" "${REFERENCE_FASTA}.fai" > "$padded_mask"
    n_orig=$(awk '$0 !~ /^(#|track|browser)/ && NF >= 3 { n++ } END { print n+0 }' "$FILTER_MASK_TOTAL_BED")
    n_padded=$(( $(wc -l < "$padded_mask") - n_orig ))
    if (( n_padded > 0 )); then
        log "Mask has no interval on $n_padded chromosome(s); masking their last base so grenedalf accepts the mask ($padded_mask)"
        FILTER_MASK_TOTAL_BED="$padded_mask"
    fi
fi

# Build grenedalf sync command
GRENEDALF_CMD=("$GRENEDALF" sync)
if [[ -n "$SYNC_FILE" ]]; then
    GRENEDALF_CMD+=(--sync-path "$SYNC_FILE")
else
    for bam in "${BAM_FILES[@]}"; do
        GRENEDALF_CMD+=(--sam-path "$bam")
    done
    GRENEDALF_CMD+=(--sam-min-map-qual "$MIN_MAPQ" --sam-min-base-qual "$MIN_BASEQ")
fi
GRENEDALF_CMD+=(
    --reference-genome-fasta "$REFERENCE_FASTA"
    --filter-sample-min-read-depth "$MIN_COVERAGE"
    --filter-sample-max-read-depth "$MAX_COVERAGE"
    --filter-sample-min-count "$MIN_COUNT"
    --filter-total-min-read-depth "$MIN_TOTAL_READ_DEPTH"
)
if [[ "$BIALLELIC_ONLY" == true ]]; then
    GRENEDALF_CMD+=(--filter-total-only-biallelic-snps)
else
    GRENEDALF_CMD+=(--filter-total-only-snps)
fi
if [[ "$MIN_COUNT_FILTER" -gt 0 ]]; then
    GRENEDALF_CMD+=(--filter-total-snp-min-count "$MIN_COUNT_FILTER")
else
    GRENEDALF_CMD+=(--filter-total-snp-min-frequency "$MIN_FREQUENCY")
fi
[[ -n "$FILTER_MASK_TOTAL_BED" ]] && GRENEDALF_CMD+=(--filter-mask-total-bed "$FILTER_MASK_TOTAL_BED")
[[ -n "$FILTER_MASK_TOTAL_FASTA" ]] && GRENEDALF_CMD+=(--filter-mask-total-fasta "$FILTER_MASK_TOTAL_FASTA")
if [[ ${#SAMPLE_NAMES[@]} -gt 0 ]]; then
    GRENEDALF_CMD+=(--rename-samples-list "${WORK_DIR}/rename_samples.csv")
fi
GRENEDALF_CMD+=(--out-dir "$WORK_DIR" --allow-file-overwriting --threads "$THREADS")

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "Would execute: ${GRENEDALF_CMD[*]}"
    log_dry_run "  (output streamed via FIFO ${WORK_DIR}/sync.sync; passing positions -> $OUT_SYNC, $OUT_TSV)"
    log_dry_run "Would execute: $BCFTOOLS view -Ob -o $OUT_BCF <sites VCF> && $BCFTOOLS index -f $OUT_BCF"
    exit 0
fi

rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT

if [[ ${#SAMPLE_NAMES[@]} -gt 0 ]]; then
    for i in "${!BAM_FILES[@]}"; do
        echo "$(basename "${BAM_FILES[$i]}" .bam),${SAMPLE_NAMES[$i]}"
    done > "${WORK_DIR}/rename_samples.csv"
fi

FIFO="${WORK_DIR}/sync.sync"
VCF_TMP="${WORK_DIR}/sites.vcf"
mkfifo "$FIFO"

# Convert passing sync positions to sync, TSV, and VCF. Sync base order: A:T:C:G:N:D.
awk -F'\t' -v OFS='\t' -v sync_out="$OUT_SYNC" -v tsv_out="$OUT_TSV" -v vcf_out="$VCF_TMP" \
    -v fai="${REFERENCE_FASTA}.fai" '
function write_headers(   s, h) {
    h = "chr\tpos\tref\talt\tdp\talt_freq"
    for (s = 1; s <= n; s++) h = h "\t" name[s] "_AD\t" name[s] "_alt_freq"
    print h > tsv_out
    print "##fileformat=VCFv4.2" > vcf_out
    print "##source=snp_call.sh (grenedalf sync read counts; genotypes not called)" > vcf_out
    while ((getline line < fai) > 0) {
        split(line, c, "\t")
        print "##contig=<ID=" c[1] ",length=" c[2] ">" > vcf_out
    }
    print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth (A,C,G,T) summed over pools passing filters\">" > vcf_out
    print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT/(REF+ALT) read-count frequency over pools passing filters\">" > vcf_out
    print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read counts for REF and ALT\">" > vcf_out
    print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth (A,C,G,T)\">" > vcf_out
    h = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    for (s = 1; s <= n; s++) h = h "\t" name[s]
    print h > vcf_out
    headers_done = 1
}
BEGIN { base[1] = "A"; base[2] = "T"; base[3] = "C"; base[4] = "G"; idx["A"] = 1; idx["T"] = 2; idx["C"] = 3; idx["G"] = 4 }
/^#/ {
    if (!headers_done) {
        n = NF - 3
        for (s = 1; s <= n; s++) name[s] = $(s + 3)
        print > sync_out
        write_headers()
    }
    next
}
{
    if (!headers_done) {
        n = NF - 3
        for (s = 1; s <= n; s++) name[s] = "S" s
        write_headers()
    }
    r = idx[toupper($3)]
    if (!r) { bad_ref++; next }
    any = 0
    for (b = 1; b <= 4; b++) tot[b] = 0
    for (s = 1; s <= n; s++) {
        miss[s] = ($(s + 3) ~ /^\./)
        if (miss[s]) continue
        any = 1
        split($(s + 3), c, ":")
        for (b = 1; b <= 4; b++) { cnt[s, b] = c[b] + 0; tot[b] += cnt[s, b] }
    }
    if (!any) next
    a = 0
    for (b = 1; b <= 4; b++) if (b != r && (a == 0 || tot[b] > tot[a])) a = b
    if (tot[a] == 0) { no_alt++; next }
    print > sync_out
    dp = tot[1] + tot[2] + tot[3] + tot[4]
    af = tot[a] / (tot[r] + tot[a])
    tsv = $1 OFS $2 OFS base[r] OFS base[a] OFS dp OFS sprintf("%.6g", af)
    vcf = $1 OFS $2 OFS "." OFS base[r] OFS base[a] OFS "." OFS "PASS" OFS "DP=" dp ";AF=" sprintf("%.6g", af) OFS "AD:DP"
    for (s = 1; s <= n; s++) {
        if (miss[s]) { tsv = tsv OFS "." OFS "NA"; vcf = vcf OFS ".:."; continue }
        sr = cnt[s, r]; sa = cnt[s, a]
        sdp = cnt[s, 1] + cnt[s, 2] + cnt[s, 3] + cnt[s, 4]
        tsv = tsv OFS sr "," sa OFS ((sr + sa) > 0 ? sprintf("%.6g", sa / (sr + sa)) : "NA")
        vcf = vcf OFS sr "," sa ":" sdp
    }
    print tsv > tsv_out
    print vcf > vcf_out
    kept++
}
END {
    if (!headers_done) exit 1
    printf "%d\t%d\t%d\n", kept, bad_ref, no_alt > "/dev/stderr"
}' < "$FIFO" 2> "${WORK_DIR}/convert.stats" &
CONVERT_PID=$!

log "Running grenedalf sync (filters as calculate_fst.sh) -> passing positions only"
log "  ${GRENEDALF_CMD[*]}"
set +e
"${GRENEDALF_CMD[@]}" > "${WORK_DIR}/grenedalf.log" 2>&1
GRENEDALF_STATUS=$?
# If grenedalf failed before opening the FIFO, the converter is still blocked on open();
# opening read-write never blocks and closing it sends EOF.
if [[ $GRENEDALF_STATUS -ne 0 ]]; then
    exec 3<>"$FIFO"
    exec 3>&-
fi
wait "$CONVERT_PID"
CONVERT_STATUS=$?
set -e

if [[ $GRENEDALF_STATUS -ne 0 ]]; then
    log_error "grenedalf sync failed (exit $GRENEDALF_STATUS). Log:"
    tail -20 "${WORK_DIR}/grenedalf.log" >&2
    exit 1
fi
if [[ $CONVERT_STATUS -ne 0 ]]; then
    log_error "Converting grenedalf output failed (exit $CONVERT_STATUS)"
    cat "${WORK_DIR}/convert.stats" >&2
    exit 1
fi

read -r N_KEPT N_BAD_REF N_NO_ALT < "${WORK_DIR}/convert.stats"
[[ "${N_BAD_REF:-0}" -gt 0 ]] && log_warn "Skipped ${N_BAD_REF} positions with a non-ACGT reference base"
[[ "${N_NO_ALT:-0}" -gt 0 ]] && log_warn "Skipped ${N_NO_ALT} positions with no non-reference base"

log "Writing BCF: $OUT_BCF"
"$BCFTOOLS" view -Ob -o "$OUT_BCF" "$VCF_TMP"
"$BCFTOOLS" index -f "$OUT_BCF"

log "Done. ${N_KEPT:-0} sites."
log "  Sync: $OUT_SYNC"
log "  TSV:  $OUT_TSV (chr, pos, ref, alt, dp, alt_freq, <sample>_AD, <sample>_alt_freq)"
log "  BCF:  $OUT_BCF (FORMAT/AD, FORMAT/DP; no GT)"
log "Next step: variant_csq.sh --bcf $OUT_BCF --reference <ref> --gff <gff3> --output-dir <dir>"
