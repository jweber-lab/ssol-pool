#!/bin/bash

###############################################################################
# seq_qual_metrics.sh
#
# Compute per-window averages for coverage and mapping quality using samtools.
# Output TSV compatible with collate: chr, start, end, sample, mean_coverage,
# mean_mapping_quality. Coordinates are 1-based inclusive (genome convention).
#
# Uses:
#   - samtools depth: per-base depth, aggregated to windowed mean coverage
#   - samtools view: aligned reads with MAPQ, aggregated to windowed mean MAPQ
#
# Author: ssol-poolseq
# Usage: See below or run with --help
###############################################################################

set -euo pipefail

BAM_FILES=()
SAMPLE_INFO_CSV=""
OUTPUT_DIR=""
REFERENCE_GENOME=""
WINDOW_SIZES=()
STEP_SIZES=()
THREADS=1
DRY_RUN=false
SAMPLE_NAMES=()

# Log file (set after output-dir is known)
SCRIPT_NAME=$(basename "$0")
LOG_FILE=""
LOG_REDIRECT_ACTIVE=false

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
  --reference-genome FILE     Reference genome FASTA (optional; use .fai for chromosome order in output)
  -w, --window-size N         Window size in bp (can be specified multiple times; default: 10000)
  --step-size N               Step size in bp (can be specified multiple times; default: 5000 or window/2)
  -t, --threads N             Threads for samtools (default: 1)
  --dry-run                   Preview commands without executing
  -h, --help                  Show this help

Output (per sample, per window/step combination):
  {output_dir}/{sample}/seq_qual_metrics_w{W}_s{S}.tsv  (one file per --window-size)
  Columns: chr, start, end, sample, mean_coverage, mean_mapping_quality
  (start/end are 1-based inclusive)
  Logs: {output_dir}/log/seq_qual_metrics_YYYYmmdd_HHMMSS.log

Example:
  $0 --bam sample1.bam --bam sample2.bam -o out --window-size 10000 --step-size 5000
  $0 --sample-info sample_info.csv -o out --window-size 10000 --reference-genome /path/to/ref.fa
  $0 --bam s.bam -o out --window-size 1000 --window-size 10000 --step-size 500 --step-size 5000
EOF
}

log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" >> "$LOG_FILE"
    fi
}
log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[ERROR] $1" >> "$LOG_FILE"
    fi
}
log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[WARN] $1" >> "$LOG_FILE"
    fi
}
log_dry_run() {
    echo -e "${YELLOW}[DRY-RUN]${NC} $1"
    if [[ -n "$LOG_FILE" && -n "$1" && "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[DRY-RUN] $1" >> "$LOG_FILE"
    fi
}

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
        --reference-genome) REFERENCE_GENOME="$2"; shift 2 ;;
        -w|--window-size) WINDOW_SIZES+=("$2"); shift 2 ;;
        --step-size) STEP_SIZES+=("$2"); shift 2 ;;
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

# Initialize log file
mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +'%Y%m%d_%H%M%S')
LOG_FILE="${LOG_DIR}/${SCRIPT_NAME%.*}_${TIMESTAMP}.log"
if [[ "$DRY_RUN" == false ]]; then
    touch "$LOG_FILE"
    {
        echo "=========================================="
        echo "$SCRIPT_NAME"
        echo "Started: $(date +'%Y-%m-%d %H:%M:%S')"
        echo "Log file: $LOG_FILE"
        echo "=========================================="
    } >> "$LOG_FILE"
    log "Log file: $LOG_FILE"
    exec 3>&1 4>&2
    exec 1> >(tee -a "$LOG_FILE" >&3)
    exec 2>&1
    LOG_REDIRECT_ACTIVE=true
fi

# Resolve BAMs and sample names from CSV if provided
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
        name="${col1}"
        bam="${col5}"
        [[ -z "$bam" ]] && bam="${col2}"
        [[ "$bam" != /* ]] && bam="${csv_dir}/${bam}"
        if [[ -f "$bam" ]]; then
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

# Build window/step arrays (default: one scale 10000/5000; else step = window/2 when not specified)
declare -a WINDOW_SIZE_ARRAY
declare -a STEP_SIZE_ARRAY
if [[ ${#WINDOW_SIZES[@]} -eq 0 ]]; then
    WINDOW_SIZE_ARRAY=(10000)
    STEP_SIZE_ARRAY=(5000)
    log "Using default window=10000, step=5000 (step = window/2)"
elif [[ ${#WINDOW_SIZES[@]} -gt 0 ]]; then
    WINDOW_SIZE_ARRAY=("${WINDOW_SIZES[@]}")
    STEP_SIZE_ARRAY=()
    for i in "${!WINDOW_SIZE_ARRAY[@]}"; do
        if [[ ${#STEP_SIZES[@]} -gt i && -n "${STEP_SIZES[$i]:-}" ]]; then
            STEP_SIZE_ARRAY[$i]="${STEP_SIZES[$i]}"
            log "Window/step scale $((i+1)): window=${WINDOW_SIZE_ARRAY[$i]}, step=${STEP_SIZES[$i]} (from --step-size)"
        else
            STEP_SIZE_ARRAY[$i]=$((${WINDOW_SIZE_ARRAY[$i]} / 2))
            log "Window/step scale $((i+1)): window=${WINDOW_SIZE_ARRAY[$i]}, step=${STEP_SIZE_ARRAY[$i]} (default: step = window/2)"
        fi
    done
fi
if [[ ${#WINDOW_SIZE_ARRAY[@]} -gt 1 ]]; then
    log "Multi-scale mode: ${#WINDOW_SIZE_ARRAY[@]} window/step combinations"
fi

# Start processing
if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN MODE: Commands will be previewed but not executed"
    log_dry_run ""
    log_dry_run "Samples to process: ${#BAM_FILES[@]}"
    for i in "${!BAM_FILES[@]}"; do
        log_dry_run "  Sample $((i+1)): ${BAM_FILES[$i]}"
    done
    log_dry_run "Window/step combinations: ${#WINDOW_SIZE_ARRAY[@]}"
    log_dry_run ""
fi

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

    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
        log "Sample: $sample_name — indexing BAM..."
        dry_run_cmd samtools index -@ "$THREADS" "$bam" || true
    fi

    for scale_idx in "${!WINDOW_SIZE_ARRAY[@]}"; do
        WINDOW_SIZE="${WINDOW_SIZE_ARRAY[$scale_idx]}"
        STEP_SIZE="${STEP_SIZE_ARRAY[$scale_idx]}"
        out_tsv="${sample_out}/seq_qual_metrics_w${WINDOW_SIZE}_s${STEP_SIZE}.tsv"

        # Coordinates: 1-based inclusive (genome convention). Window start = int((pos-1)/step)*step + 1; end = start + window_size - 1.
        log "Processing sample: $sample_name (window=$WINDOW_SIZE, step=$STEP_SIZE)"
        log "  BAM: $bam"
        log "  Window size: $WINDOW_SIZE bp"
        log "  Step size: $STEP_SIZE bp"

        # 1) Coverage: samtools depth (1-based pos) -> windowed mean. Windows 1-based inclusive.
        log "  Computing per-window mean coverage..."
        depth_tsv=$(mktemp -t depth.XXXXXX.tsv)
        if [[ "$DRY_RUN" != true ]]; then
            # Window start 1-based: int((pos-1)/s)*s + 1; end inclusive: start + w - 1
            samtools depth -aa -@ "$THREADS" "$bam" | awk -v w="$WINDOW_SIZE" -v s="$STEP_SIZE" '
        BEGIN { OFS="\t" }
        {
            chr=$1; pos=$2+0; depth=$3+0
            win = int((pos-1)/s)*s + 1
            key = chr "\t" win
            sum[key] += depth
            cnt[key] += 1
        }
        END {
            for (k in sum) {
                split(k, a, "\t")
                print a[1], a[2], a[2]+w-1, sum[k]/cnt[k]
            }
        }
        ' | sort -k1,1 -k2,2n > "$depth_tsv"
    else
        log_dry_run "samtools depth -aa $bam | awk ... | sort > depth_tsv"
    fi

    # 2) MAPQ: samtools view -> chr, start, end, mapq -> windowed mean
    log "  Computing per-window mean mapping quality..."
    mapq_tsv=$(mktemp -t mapq.XXXXXX.tsv)
    if [[ "$DRY_RUN" != true ]]; then
        # BAM pos is 1-based; window start = int((pos-1)/s)*s+1; overlap when win <= end && (win+w-1) >= pos
        samtools view -F 4 -@ "$THREADS" "$bam" | awk -v w="$WINDOW_SIZE" -v s="$STEP_SIZE" '
        function ref_len(cig) {
            r = 0
            n = length(cig)
            num = ""
            for (i = 1; i <= n; i++) {
                c = substr(cig, i, 1)
                if (c ~ /[0-9]/) num = num c
                else if (c ~ /[MDN=X]/) { r += num+0; num = "" }
                else num = ""
            }
            return r
        }
        {
            chr = $3
            pos = $4+0
            cig = $6
            mapq = $5+0
            len = ref_len(cig)
            end = pos + len - 1
            win = int((pos-1)/s)*s + 1
            while (win <= end) {
                if ((win + w - 1) >= pos) {
                    key = chr "\t" win
                    sum[key] += mapq
                    cnt[key] += 1
                }
                win += s
            }
        }
        END {
            for (k in sum) {
                split(k, a, "\t")
                print a[1], a[2], a[2]+w-1, sum[k]/cnt[k]
            }
        }
        ' | sort -k1,1 -k2,2n > "$mapq_tsv"
    else
        log_dry_run "samtools view -F 4 $bam | awk ... | sort > mapq_tsv"
    fi

    # 3) Merge coverage and MAPQ by (chr, start); add sample column. Output: chr, start, end (1-based inclusive), sample, mean_coverage, mean_mapping_quality.
    if [[ "$DRY_RUN" != true ]]; then
        (
            printf '%s\n' 'chr	start	end	sample	mean_coverage	mean_mapping_quality'
            awk -v sn="$sample_name" '
                BEGIN { OFS="\t" }
                NR==FNR {
                    depth[$1"\t"$2]=$4
                    end[$1"\t"$2]=$3
                    next
                }
                {
                    key=$1"\t"$2
                    if (key in depth)
                        print $1, $2, end[key], sn, depth[key], $4
                }
            ' "$depth_tsv" "$mapq_tsv" | sort -k1,1 -k2,2n
        ) > "$out_tsv"
        rm -f "$depth_tsv" "$mapq_tsv"
        # Optional: sort by reference chromosome order (from .fai)
        # Prepend order to header too so cut -f2- does not drop chr from the header
        if [[ -n "$REFERENCE_GENOME" ]] && [[ -f "${REFERENCE_GENOME}.fai" ]]; then
            log "  Sorting by reference chromosome order"
            awk 'NR==FNR{ord[$1]=NR; next} FNR==1{print 0, $0; next} {print (ord[$1]?ord[$1]:999999), $0}' \
                "${REFERENCE_GENOME}.fai" "$out_tsv" | sort -k1,1n -k3,3n -k4,4n | cut -f2- > "${out_tsv}.tmp" && mv "${out_tsv}.tmp" "$out_tsv"
        fi
        log "  Wrote $out_tsv"
    else
        log_dry_run "Merge depth and MAPQ -> $out_tsv"
    fi
    done
done

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run ""
    log_dry_run "DRY-RUN complete. All checks passed."
    log_dry_run "Run without --dry-run to execute the analysis."
else
    log "Done. Output under $OUTPUT_DIR"
fi
