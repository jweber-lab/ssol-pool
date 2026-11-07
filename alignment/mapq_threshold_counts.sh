#!/bin/bash

###############################################################################
# mapq_threshold_counts.sh
#
# Report number of mapped reads at multiple MAPQ thresholds per BAM (histogram
# then cumulative-sum; one pass, no per-read threshold loop). Use to choose
# MAPQ cutoffs that balance sequence quality and mapped depth.
#
# Author: ssol-poolseq
# Usage: See below or run with --help
###############################################################################

set -euo pipefail

BAM_FILES=()
SAMPLE_NAMES=()
SAMPLE_INFO_CSV=""
OUTPUT_DIR=""
MAPQ_THRESHOLDS="0,5,10,15,20,25,30,40,50"
PRIMARY_ONLY=false
REFERENCE_FASTA=""
READ_LENGTH=150
THREADS=1
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Input (choose one):
  -b, --bam FILE              BAM file(s) - can be specified multiple times
  -i, --sample-info FILE      CSV with columns: sample_name, bam_file (or bam_file path)

Required:
  -o, --output-dir DIR        Output directory for TSV files

Optional:
  --mapq-thresholds LIST      Comma-separated MAPQ thresholds (default: 0,5,10,15,20,25,30,40,50)
  --primary-only              Count only primary alignments (exclude secondary/supplementary)
  --reference FILE            Reference FASTA (with .fai) for optional estimated_mean_depth
  --read-length N             Approximate read length for depth estimate (default: 150)
  -t, --threads N             Threads for samtools view (default: 1)
  --dry-run                   Preview commands without executing
  -h, --help                  Show this help

Output:
  Per sample: {output_dir}/{sample}/mapq_threshold_counts.tsv
  Combined:  {output_dir}/mapq_threshold_counts_all.tsv
  Columns:   mapq_threshold, n_mapped_reads, fraction_retained [, estimated_mean_depth]

Example:
  $0 --sample-info sample_info.csv -o mapq_stats
  $0 --bam s1.bam --bam s2.bam -o mapq_stats --mapq-thresholds "0,10,20,30"
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

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        -b|--bam) BAM_FILES+=("$2"); shift 2 ;;
        -i|--sample-info) SAMPLE_INFO_CSV="$2"; shift 2 ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --mapq-thresholds) MAPQ_THRESHOLDS="$2"; shift 2 ;;
        --primary-only) PRIMARY_ONLY=true; shift ;;
        --reference) REFERENCE_FASTA="$2"; shift 2 ;;
        --read-length) READ_LENGTH="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing --output-dir"
    usage
    exit 1
fi

# Resolve BAMs and sample names from CSV if provided
if [[ -n "$SAMPLE_INFO_CSV" ]]; then
    if [[ ! -f "$SAMPLE_INFO_CSV" ]]; then
        log_error "Sample info file not found: $SAMPLE_INFO_CSV"
        exit 1
    fi
    csv_dir=$(dirname "$SAMPLE_INFO_CSV")
    while IFS=',' read -r col1 col2 col3 col4 col5 rest; do
        case "$col1" in
            [Ss]ample_name|[Ss]ample|[Bb]am_file|"# Sample") continue ;;
        esac
        name="${col1}"
        bam="${col5}"
        [[ -z "$bam" ]] && bam="${col2}"
        [[ "$bam" != /* ]] && bam="${csv_dir}/${bam}"
        if [[ -n "$bam" && -f "$bam" ]]; then
            SAMPLE_NAMES+=("$name")
            BAM_FILES+=("$bam")
        fi
    done < <(grep -v '^#' "$SAMPLE_INFO_CSV" | tr -d '"')
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    log_error "No BAM files specified. Use --bam or --sample-info."
    usage
    exit 1
fi

check_command samtools
mkdir -p "$OUTPUT_DIR"

# Reference length for optional depth estimate (sum of column 2 in .fai)
REF_LENGTH=0
if [[ -n "$REFERENCE_FASTA" ]] && [[ -f "${REFERENCE_FASTA}.fai" ]]; then
    REF_LENGTH=$(awk '{s+=$2} END {print s+0}' "${REFERENCE_FASTA}.fai")
fi

# Build samtools view flags: -F 4 = exclude unmapped; -F 2304 = exclude secondary+supplementary if --primary-only
VIEW_FLAGS="-F 4"
[[ "$PRIMARY_ONLY" == true ]] && VIEW_FLAGS="-F 4 -F 2304"

# Combined output (append per sample)
COMBINED_TSV="${OUTPUT_DIR}/mapq_threshold_counts_all.tsv"
COMBINED_HEADER_WRITTEN=false

for i in "${!BAM_FILES[@]}"; do
    bam="${BAM_FILES[$i]}"
    if [[ ${#SAMPLE_NAMES[@]} -gt i && -n "${SAMPLE_NAMES[$i]:-}" ]]; then
        sample_name="${SAMPLE_NAMES[$i]}"
    else
        sample_name=$(basename "$bam" .bam)
        sample_name="${sample_name%_All_seq.dedup}"
        sample_name="${sample_name%.dedup}"
    fi
    sample_out="${OUTPUT_DIR}/${sample_name}"
    mkdir -p "$sample_out"
    out_tsv="${sample_out}/mapq_threshold_counts.tsv"

    log "Processing sample: $sample_name"
    log "  BAM: $bam"

    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would run: samtools view $VIEW_FLAGS -@ $THREADS $bam | awk (histogram + cumulative counts) -> $out_tsv"
        continue
    fi

    # One pass: build MAPQ histogram (0..60), then in END compute per-threshold counts from cumulative sum.
    # Awk receives thresholds as comma-separated; outputs mapq_threshold, n_mapped_reads, fraction_retained [, estimated_mean_depth]
    awk -v thresh_list="$MAPQ_THRESHOLDS" -v ref_len="$REF_LENGTH" -v read_len="$READ_LENGTH" '
        BEGIN {
            OFS = "\t"
            n = split(thresh_list, T, ",")
            for (i = 1; i <= n; i++) thresh[i] = T[i] + 0
            max_bin = 60
        }
        {
            mq = $5 + 0
            if (mq < 0) mq = 0
            if (mq > max_bin) mq = max_bin
            hist[mq]++
        }
        END {
            total = 0
            for (j = 0; j <= max_bin; j++) total += hist[j] + 0
            for (i = 1; i <= n; i++) {
                t = thresh[i]
                cnt = 0
                for (j = t; j <= max_bin; j++) cnt += hist[j] + 0
                frac = (total > 0) ? (cnt / total) : 0
                depth = (ref_len > 0 && read_len > 0) ? (cnt * read_len / ref_len) : ""
                if (depth != "")
                    print t, cnt, frac, depth
                else
                    print t, cnt, frac
            }
        }
    ' <(samtools view $VIEW_FLAGS -@ "$THREADS" "$bam") > "${out_tsv}.tmp"

    # Write per-sample TSV with header (include estimated_mean_depth column only if we have ref)
    if [[ $REF_LENGTH -gt 0 ]]; then
        echo -e "mapq_threshold\tn_mapped_reads\tfraction_retained\testimated_mean_depth" > "$out_tsv"
    else
        echo -e "mapq_threshold\tn_mapped_reads\tfraction_retained" > "$out_tsv"
    fi
    cat "${out_tsv}.tmp" >> "$out_tsv"
    rm -f "${out_tsv}.tmp"

    # Append to combined TSV (with sample column)
    if [[ "$COMBINED_HEADER_WRITTEN" == false ]]; then
        if [[ $REF_LENGTH -gt 0 ]]; then
            echo -e "sample\tmapq_threshold\tn_mapped_reads\tfraction_retained\testimated_mean_depth" > "$COMBINED_TSV"
        else
            echo -e "sample\tmapq_threshold\tn_mapped_reads\tfraction_retained" > "$COMBINED_TSV"
        fi
        COMBINED_HEADER_WRITTEN=true
    fi
    awk -v sn="$sample_name" 'BEGIN{OFS="\t"} {print sn, $0}' "$out_tsv" | tail -n +2 >> "$COMBINED_TSV"

    log "  Wrote $out_tsv"
done

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "Would write combined TSV: $COMBINED_TSV"
    log_dry_run "DRY-RUN complete."
else
    log "Done. Per-sample TSVs under $OUTPUT_DIR/<sample>/; combined: $COMBINED_TSV"
fi
