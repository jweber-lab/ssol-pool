#!/bin/bash

###############################################################################
# variant_call.sh
#
# Generate BCF (variant calls) from BAM(s) using bcftools mpileup and call.
# Pipeline: bcftools mpileup -> bcftools call -> BCF.
# Suitable for pool-seq when multiple BAMs are provided (one sample per BAM).
# Uses bcftools call -m (multiallelic caller) for rare variants and >2 alleles.
#
# Note: Does not use any .mpileup file from process_poolseq.sh. That file is
# samtools mpileup (text) for popoolation2 sync; this script uses bcftools
# mpileup (reads BAM directly, outputs BCF) and produces different output.
#
# Typical next step (post-call filtering):
#   bcftools view -i 'QUAL>=20 && DP>10' calls.bcf -Ob -o calls_filtered.bcf
# Pool-seq: often use --ploidy 1; for very high coverage increase --max-depth.
#
# Optional: --also-tsv writes a TSV (chr, pos, ref, alt, qual, dp, per-sample AD) for collate and AF.
# Convert-only: --bcf-to-tsv FILE converts an existing BCF to TSV without running call.
#
# Author: ssol-poolseq
# Usage: See below or run with --help
###############################################################################

set -euo pipefail

BAM_FILES=()
SAMPLE_INFO_CSV=""
REFERENCE_FASTA=""
OUTPUT_BCF=""
OUTPUT_TSV=""
OUTPUT_DIR=""
BCF_TO_TSV_INPUT=""   # When set, only convert BCF to TSV (no mpileup/call)
BCFTOOLS=""
THREADS=1
MIN_MAPQ=20
MIN_BASEQ=20
MAX_DEPTH=1000
PLOIDY=""
SKIP_INDELS=false
INDELS_ONLY=false
ALSO_TSV=false
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Variant calling (default): BAM(s) -> mpileup -> call -> BCF

Input (choose one):
  -b, --bam FILE              BAM file(s) - can be specified multiple times
  -i, --sample-info FILE      CSV with columns: sample_name, bam_file

Required (unless using --bcf-to-tsv):
  -f, --reference FILE        Reference genome FASTA (must be indexed: .fai)
  -o, --output FILE           Output BCF path (e.g. calls.bcf)
  OR
  --output-dir DIR            Output directory; writes calls.bcf inside DIR

Convert-only (no BAM/reference needed):
  --bcf-to-tsv FILE           Convert existing BCF/VCF to TSV (chr, pos, ref, alt, qual, dp, per-sample AD)
  -o, --output FILE           Output TSV path (required with --bcf-to-tsv)
  OR --output-dir DIR         Writes variants.tsv inside DIR

Optional:
  -g, --bcftools PATH         Path to bcftools (default: bcftools in PATH)
  -t, --threads N             Threads for mpileup (default: 1)
  --min-mapq N                Minimum mapping quality (default: 20)
  --min-baseq N               Minimum base quality (default: 20)
  --max-depth N               Max reads per input per position (default: 1000)
  --ploidy N                  Ploidy for call (default: 2; use 1 for haploid pools)
  --skip-indels               Output SNPs only (bcftools call -V indels)
  --indels-only               Output indels only (bcftools call -V snps)
  --also-tsv                  Also write TSV (for collate) alongside BCF
  --output-tsv FILE           TSV path when using --also-tsv (default: <bcf_base>.tsv)
  --dry-run                   Preview commands without executing
  -h, --help                  Show this help

Example:
  $0 --reference ref.fa --bam s1.bam --bam s2.bam -o calls.bcf
  $0 --reference ref.fa --sample-info sample_info.csv --output-dir out --ploidy 1 --also-tsv
  $0 --bcf-to-tsv calls.bcf -o calls.tsv
  $0 --bcf-to-tsv calls.bcf --output-dir out
  Typical next step: bcftools view -i 'QUAL>=20 && DP>10' calls.bcf -Ob -o calls_filtered.bcf
EOF
}

log() { echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_dry_run() { echo -e "${YELLOW}[DRY-RUN]${NC} $1"; }

dry_run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: $*"
        return 0
    fi
    "$@"
}

check_command() {
    command -v "$1" &>/dev/null || { log_error "$1 not found"; exit 1; }
}

# Write BCF to TSV with chr, pos, ref, alt, qual, dp, and per-sample AD (ref,alt for AF).
# AD columns are named {sample}_AD; value is "ref,alt" so AF = alt/(ref+alt).
write_bcf_to_tsv() {
    local bcf="$1"
    local tsv="$2"
    local bcftools="${3:-bcftools}"
    local samples
    readarray -t samples < <("$bcftools" query -l "$bcf" 2>/dev/null || true)
    local header=$'chr\tpos\tref\talt\tqual\tdp'
    local s
    for s in "${samples[@]}"; do
        [[ -z "$s" ]] && continue
        header="${header}"$'\t'"${s}_AD"
    done
    { printf '%s\n' "$header"; "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP[\t%AD]\n' "$bcf"; } > "$tsv"
}

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bam) BAM_FILES+=("$2"); shift 2 ;;
        -i|--sample-info) SAMPLE_INFO_CSV="$2"; shift 2 ;;
        -f|--reference) REFERENCE_FASTA="$2"; shift 2 ;;
        -o|--output) OUTPUT_BCF="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --output-tsv) OUTPUT_TSV="$2"; shift 2 ;;
        --bcf-to-tsv) BCF_TO_TSV_INPUT="$2"; shift 2 ;;
        -g|--bcftools) BCFTOOLS="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
        --min-baseq) MIN_BASEQ="$2"; shift 2 ;;
        --max-depth) MAX_DEPTH="$2"; shift 2 ;;
        --ploidy) PLOIDY="$2"; shift 2 ;;
        --skip-indels) SKIP_INDELS=true; shift ;;
        --indels-only) INDELS_ONLY=true; shift ;;
        --also-tsv) ALSO_TSV=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ "$SKIP_INDELS" == true && "$INDELS_ONLY" == true ]]; then
    log_error "Cannot use both --skip-indels and --indels-only"
    usage
    exit 1
fi

# Resolve output paths: convert-only uses TSV path; call mode uses BCF path (and optionally TSV)
if [[ -n "$BCF_TO_TSV_INPUT" ]]; then
    # Convert-only mode: need output path for TSV
    if [[ -n "$OUTPUT_DIR" && -z "$OUTPUT_BCF" ]]; then
        OUTPUT_TSV="${OUTPUT_DIR}/variants.tsv"
    elif [[ -n "$OUTPUT_BCF" ]]; then
        OUTPUT_TSV="$OUTPUT_BCF"
        OUTPUT_BCF=""
    fi
    if [[ -z "$OUTPUT_TSV" ]]; then
        log_error "With --bcf-to-tsv, specify -o FILE or --output-dir DIR for TSV output"
        usage
        exit 1
    fi
else
    # Call mode: need reference and BCF output
    if [[ -z "$REFERENCE_FASTA" ]]; then
        log_error "Missing --reference"
        usage
        exit 1
    fi
    if [[ -n "$OUTPUT_DIR" && -z "$OUTPUT_BCF" ]]; then
        OUTPUT_BCF="${OUTPUT_DIR}/calls.bcf"
    fi
    if [[ -z "$OUTPUT_BCF" ]]; then
        log_error "Missing --output or --output-dir"
        usage
        exit 1
    fi
    if [[ "$ALSO_TSV" == true && -z "$OUTPUT_TSV" ]]; then
        OUTPUT_TSV="${OUTPUT_BCF%.bcf}"
        OUTPUT_TSV="${OUTPUT_TSV%.vcf}.tsv"
    fi
fi

# Convert-only: just BCF -> TSV and exit
if [[ -n "$BCF_TO_TSV_INPUT" ]]; then
    if [[ ! -f "$BCF_TO_TSV_INPUT" ]]; then
        log_error "BCF file not found: $BCF_TO_TSV_INPUT"
        exit 1
    fi
    : "${BCFTOOLS:=bcftools}"
    check_command "$BCFTOOLS"
    mkdir -p "$(dirname "$OUTPUT_TSV")"
    log "Converting BCF to TSV: $BCF_TO_TSV_INPUT -> $OUTPUT_TSV"
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "write_bcf_to_tsv (chr, pos, ref, alt, qual, dp, per-sample _AD) -> $OUTPUT_TSV"
    else
        write_bcf_to_tsv "$BCF_TO_TSV_INPUT" "$OUTPUT_TSV" "$BCFTOOLS"
        log "Done. Output: $OUTPUT_TSV (columns: chr, pos, ref, alt, qual, dp, <sample>_AD; AD is ref,alt for AF)"
    fi
    exit 0
fi

# Resolve BAMs from CSV if provided
if [[ -n "$SAMPLE_INFO_CSV" ]]; then
    if [[ ! -f "$SAMPLE_INFO_CSV" ]]; then
        log_error "Sample info file not found: $SAMPLE_INFO_CSV"
        exit 1
    fi
    csv_dir=$(dirname "$SAMPLE_INFO_CSV")
    while IFS=',' read -r col1 col2 col3 col4 col5 rest; do
        # Skip header row (first column is sample_name, Sample, bam_file, etc.)
        case "$col1" in
            [Ss]ample_name|[Ss]ample|[Bb]am_file|"# Sample") continue ;;
        esac
        bam="${col5}"
        [[ -z "$bam" ]] && bam="${col2}"
        [[ "$bam" != /* ]] && bam="${csv_dir}/${bam}"
        if [[ -n "$bam" && -f "$bam" ]]; then
            BAM_FILES+=("$bam")
        fi
    done < <(grep -v '^#' "$SAMPLE_INFO_CSV" | tr -d '"')
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    log_error "No BAM files specified. Use --bam or --sample-info."
    usage
    exit 1
fi

if [[ -z "$BCFTOOLS" ]]; then
    if command -v bcftools &>/dev/null; then
        BCFTOOLS=bcftools
    else
        log_error "bcftools not found. Install or set --bcftools."
        exit 1
    fi
fi

check_command "$BCFTOOLS"
check_command samtools

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    log_error "Reference not found: $REFERENCE_FASTA"
    exit 1
fi

if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
    log "Indexing reference with samtools faidx..."
    dry_run_cmd samtools faidx "$REFERENCE_FASTA" || true
fi

mkdir -p "$(dirname "$OUTPUT_BCF")"

# Build mpileup: bcftools mpileup -f ref -Q min_baseq -q min_mapq -d max_depth -Ou bams...
# Note: bcftools uses --threads (samtools uses -@)
MPILEUP_CMD=(
    "$BCFTOOLS" mpileup
    -f "$REFERENCE_FASTA"
    -Q "$MIN_BASEQ"
    -q "$MIN_MAPQ"
    -d "$MAX_DEPTH"
    -Ou
    --threads "$THREADS"
    "${BAM_FILES[@]}"
)

# Build call: bcftools call -m -v -Ou [--ploidy 1] [-V indels|snps]
# -m: multiallelic caller (needed for pool-seq rare variants and >2 alleles)
CALL_CMD=(
    "$BCFTOOLS" call
    -m
    -v
    -Ou
)
[[ -n "$PLOIDY" ]] && CALL_CMD+=(--ploidy "$PLOIDY")
[[ "$SKIP_INDELS" == true ]] && CALL_CMD+=(-V indels)
[[ "$INDELS_ONLY" == true ]] && CALL_CMD+=(-V snps)

log "Running mpileup | call -> $OUTPUT_BCF"
if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "${MPILEUP_CMD[@]} | ${CALL_CMD[@]} | $BCFTOOLS view -Ob -o $OUTPUT_BCF"
else
    "${MPILEUP_CMD[@]}" | "${CALL_CMD[@]}" | "$BCFTOOLS" view -Ob -o "$OUTPUT_BCF"
fi

log "Indexing BCF..."
dry_run_cmd "$BCFTOOLS" index -f "$OUTPUT_BCF" || true

if [[ "$ALSO_TSV" == true ]] && [[ -n "$OUTPUT_TSV" ]]; then
    log "Writing TSV for collate (with per-sample AD for AF): $OUTPUT_TSV"
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "write_bcf_to_tsv (chr, pos, ref, alt, qual, dp, per-sample _AD) -> $OUTPUT_TSV"
    else
        write_bcf_to_tsv "$OUTPUT_BCF" "$OUTPUT_TSV" "$BCFTOOLS"
        log "TSV written: $OUTPUT_TSV (columns: chr, pos, ref, alt, qual, dp, <sample>_AD; AD is ref,alt for AF)"
    fi
fi

log "Done. Output: $OUTPUT_BCF"
log "Typical next step (post-call filter): bcftools view -i 'QUAL>=20 && DP>10' $OUTPUT_BCF -Ob -o <filtered.bcf>"
