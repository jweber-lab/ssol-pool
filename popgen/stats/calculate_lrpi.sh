#!/bin/bash

###############################################################################
# calculate_lrpi.sh
#
# Log ratio of nucleotide diversity between two populations:
#   LRPI = ln(pi_A / pi_B)
# from per-sample diversity tables written by calculate_pi_theta.sh.
#
# This is a directed two-population statistic. Swapping A and B negates it.
# pi must be finite and greater than 0 in both populations; otherwise the
# window is NA. Windows present in only one sample are dropped.
#
# Computation is done in R (calculate_lrpi.R). Output is TSV, CSV, or HDF5.
#
# Usage: See README.md or run with --help
###############################################################################

set -euo pipefail

DIVERSITY_DIR=""
DIVERSITY_A=""
DIVERSITY_B=""
POP1_NAME=""
POP2_NAME=""
COMPARAND=""
ALL_PAIRS=false
LRPI_COMPARISONS_FILE=""
SAMPLE_INFO_CSV=""
COMPUTE_DUP_COMPS=false
OUTPUT_DIR=""
FILE_PREFIX=""
OUTPUT_FORMAT="tsv"
DROP_ALL_NA=false
DRY_RUN=false
VERBOSE=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m'

SCRIPT_NAME=$(basename "$0")
LOG_FILE=""
LOG_REDIRECT_ACTIVE=false

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Input (choose one):
  -d, --diversity-dir DIR     Directory of calculate_pi_theta.sh tables
                              (diversity_w{W}_s{S}_sample{NAME}.tsv or .csv).
                              Each window/step is joined separately.
  --diversity-a FILE          One diversity table for population 1 (with --diversity-b)
  --diversity-b FILE          One diversity table for population 2

Comparison (choose one). LRPI = ln(pi_pop1 / pi_pop2). Order matters:
  --pop1-name NAME            Numerator population
  --pop2-name NAME            Denominator population
  --lrpi-comparisons FILE     CSV with columns pop1, pop2 (one directed pair per row).
                              The same ordered pair is skipped unless --compute-dup-comps.
                              The swapped pair is a different comparison.
  --comparand NAME            ln(pi_NAME / pi_other) for every other sample in the directory
  --all-pairs                 Every unordered pair once. Names are sorted, so pop1 < pop2
                              and the value is ln(pi_pop1 / pi_pop2). Swap the names to flip the sign.

Required:
  -o, --output-dir DIR        Output directory

Optional:
  --sample-info FILE          Warn if a comparison name is missing from the sample_name column
  --compute-dup-comps         Keep duplicate ordered pairs from the comparisons file
  --file-prefix PREFIX        Prefix for output files
  --output-format FMT         tsv, csv, or hdf5 (default: tsv)
  --drop-all-na               Drop windows where LRPI is NA
  --dry-run                   Preview only
  --verbose                   Log each pair and, in R, join counts and non-positive pi
  -h, --help                  Show this help

Output (one file per diversity scale and directed pair):
  {output_dir}/{prefix_}{label}_{pop1}_{pop2}_lrpi.{tsv|csv}
  label is w{W}_s{S} from the diversity filename.
  Columns: chr, pos, start, end, pop1, pop2, pi_pop1, pi_pop2, lrpi
  hdf5: {output_dir}/{prefix_}{label}_{pop1}_{pop2}_lrpi.h5 (group /lrpi)
  Logs: {output_dir}/log/{prefix_}calculate_lrpi_YYYYmmdd_HHMMSS.log

Examples:
  $0 --diversity-dir pi_out --pop1-name Echo --pop2-name Kjer -o lrpi_out

  $0 --diversity-dir pi_out --comparand Echo -o lrpi_out

  $0 --diversity-dir pi_out --all-pairs -o lrpi_out

  $0 --diversity-dir pi_out --lrpi-comparisons pairs.csv --sample-info ../../sample_info.csv -o lrpi_out

EOF
}

log() {
    local message="${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
    echo -e "$message"
    if [[ -n "$LOG_FILE" ]] && [[ -f "$LOG_FILE" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" >> "$LOG_FILE"
    fi
}

log_error() {
    local message="${RED}[ERROR]${NC} $1"
    echo -e "$message" >&2
    if [[ -n "$LOG_FILE" ]] && [[ -f "$LOG_FILE" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[ERROR] $1" >> "$LOG_FILE"
    fi
}

log_warn() {
    local message="${YELLOW}[WARN]${NC} $1"
    echo -e "$message"
    if [[ -n "$LOG_FILE" ]] && [[ -f "$LOG_FILE" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[WARN] $1" >> "$LOG_FILE"
    fi
}

log_dry_run() {
    local message="${YELLOW}[DRY-RUN]${NC} $1"
    echo -e "$message"
    if [[ -n "$LOG_FILE" ]] && [[ -f "$LOG_FILE" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        echo "[DRY-RUN] $1" >> "$LOG_FILE"
    fi
}

check_command() {
    command -v "$1" &>/dev/null
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        exit 1
    fi
}

check_dir() {
    if [[ ! -d "$1" ]]; then
        if [[ "$DRY_RUN" == true ]]; then
            log_dry_run "Would create directory: $1"
        else
            log "Creating directory: $1"
            mkdir -p "$1"
        fi
    fi
}

sanitize_pop_name() {
    echo "$1" | sed 's/[^a-zA-Z0-9_-]/_/g'
}

# Print "window step sample" for a diversity basename, or return 1.
parse_diversity_base() {
    local base="$1"
    base="${base%.csv}"
    base="${base%.tsv}"
    if [[ "$base" =~ diversity_w([0-9]+)_s([0-9]+)_sample(.+)$ ]]; then
        echo "${BASH_REMATCH[1]} ${BASH_REMATCH[2]} ${BASH_REMATCH[3]}"
        return 0
    fi
    return 1
}

# Fill DIV_PATHS, DIV_SAMPLES, DIV_LABELS from a directory.
discover_diversity_files() {
    local dir="$1"
    DIV_PATHS=()
    DIV_SAMPLES=()
    DIV_LABELS=()
    local f parsed w s sample label key
    declare -A SEEN_SCALE_SAMPLE
    while IFS= read -r -d '' f; do
        [[ -f "$f" ]] || continue
        parsed=$(parse_diversity_base "$(basename "$f")") || {
            log_warn "Skipping diversity file with an unrecognized name: $f"
            continue
        }
        read -r w s sample <<< "$parsed"
        label="w${w}_s${s}"
        key="${label}	${sample}"
        if [[ -n "${SEEN_SCALE_SAMPLE[$key]:-}" ]]; then
            log_error "Two diversity files for ${sample} at ${label}: ${SEEN_SCALE_SAMPLE[$key]} and $f"
            exit 1
        fi
        SEEN_SCALE_SAMPLE[$key]="$f"
        DIV_PATHS+=("$f")
        DIV_SAMPLES+=("$sample")
        DIV_LABELS+=("$label")
    done < <(find "$dir" -maxdepth 1 \( -name "*diversity*.tsv" -o -name "*diversity*.csv" \) -print0 2>/dev/null)
}

# Comparisons CSV: columns pop1, pop2. Directed. Sets PAIR_POP1, PAIR_POP2.
parse_lrpi_comparisons() {
    local csv="$1"
    PAIR_POP1=()
    PAIR_POP2=()
    declare -A SEEN_KEYS
    local line_num=0 pop1_col=-1 pop2_col=-1
    while IFS= read -r line || [[ -n "$line" ]]; do
        line_num=$((line_num + 1))
        [[ -z "$line" ]] && continue
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        IFS=',' read -r -a F <<< "$line"
        local i
        for i in "${!F[@]}"; do
            F[$i]=$(echo "${F[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        done
        if [[ $line_num -eq 1 ]]; then
            for i in "${!F[@]}"; do
                local h
                h=$(echo "${F[$i]}" | tr '[:upper:]' '[:lower:]')
                case "$h" in
                    pop1) pop1_col=$i ;;
                    pop2) pop2_col=$i ;;
                esac
            done
            if [[ $pop1_col -lt 0 ]] || [[ $pop2_col -lt 0 ]]; then
                log_error "lrpi comparisons CSV must have header columns pop1, pop2 (found: $line)"
                return 1
            fi
            continue
        fi
        [[ ${#F[@]} -le $pop1_col ]] || [[ ${#F[@]} -le $pop2_col ]] && continue
        local p1="${F[$pop1_col]}" p2="${F[$pop2_col]}"
        [[ -z "$p1" ]] || [[ -z "$p2" ]] && continue
        if [[ "$p1" == "$p2" ]]; then
            log_error "Comparison pop1 and pop2 are the same ($p1) on line $line_num"
            return 1
        fi
        local key="${p1}	${p2}"
        if [[ -n "${SEEN_KEYS[$key]:-}" ]] && [[ "$COMPUTE_DUP_COMPS" != true ]]; then
            log_warn "Skipping duplicate ordered pair pop1=$p1 pop2=$p2 (line $line_num). Use --compute-dup-comps to keep it."
            continue
        fi
        SEEN_KEYS[$key]=1
        PAIR_POP1+=("$p1")
        PAIR_POP2+=("$p2")
    done < "$csv"
    if [[ ${#PAIR_POP1[@]} -eq 0 ]]; then
        log_error "No valid pairs in $csv"
        return 1
    fi
    return 0
}

# Unique labels, single-style sort not needed; sort labels for stable output.
unique_labels() {
    printf '%s\n' "${DIV_LABELS[@]}" | LC_ALL=C sort -u
}

samples_at_label() {
    local label="$1" i
    for i in "${!DIV_LABELS[@]}"; do
        [[ "${DIV_LABELS[$i]}" == "$label" ]] && printf '%s\n' "${DIV_SAMPLES[$i]}"
    done | LC_ALL=C sort -u
}

file_for() {
    local sample="$1" label="$2" i
    for i in "${!DIV_PATHS[@]}"; do
        if [[ "${DIV_SAMPLES[$i]}" == "$sample" && "${DIV_LABELS[$i]}" == "$label" ]]; then
            echo "${DIV_PATHS[$i]}"
            return 0
        fi
    done
    return 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--diversity-dir) DIVERSITY_DIR="$2"; shift 2 ;;
        --diversity-a) DIVERSITY_A="$2"; shift 2 ;;
        --diversity-b) DIVERSITY_B="$2"; shift 2 ;;
        --pop1-name) POP1_NAME="$2"; shift 2 ;;
        --pop2-name) POP2_NAME="$2"; shift 2 ;;
        --comparand) COMPARAND="$2"; shift 2 ;;
        --all-pairs) ALL_PAIRS=true; shift ;;
        --lrpi-comparisons) LRPI_COMPARISONS_FILE="$2"; shift 2 ;;
        --sample-info) SAMPLE_INFO_CSV="$2"; shift 2 ;;
        --compute-dup-comps) COMPUTE_DUP_COMPS=true; shift ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --file-prefix) FILE_PREFIX="$2"; shift 2 ;;
        --output-format) OUTPUT_FORMAT="$2"; shift 2 ;;
        --drop-all-na) DROP_ALL_NA=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        --verbose) VERBOSE=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing required argument: --output-dir"
    usage
    exit 1
fi

mode_n=0
[[ -n "$LRPI_COMPARISONS_FILE" ]] && mode_n=$((mode_n + 1))
[[ -n "$COMPARAND" ]] && mode_n=$((mode_n + 1))
[[ "$ALL_PAIRS" == true ]] && mode_n=$((mode_n + 1))
[[ -n "$POP1_NAME" || -n "$POP2_NAME" ]] && mode_n=$((mode_n + 1))
if [[ "$mode_n" -eq 0 ]]; then
    log_error "Choose a comparison: --pop1-name and --pop2-name, --lrpi-comparisons, --comparand, or --all-pairs"
    usage
    exit 1
fi
if [[ "$mode_n" -gt 1 ]]; then
    log_error "Use only one of --pop1-name/--pop2-name, --lrpi-comparisons, --comparand, or --all-pairs"
    exit 1
fi
if [[ -n "$POP1_NAME" || -n "$POP2_NAME" ]]; then
    if [[ -z "$POP1_NAME" || -z "$POP2_NAME" ]]; then
        log_error "Both --pop1-name and --pop2-name are required"
        exit 1
    fi
    if [[ "$POP1_NAME" == "$POP2_NAME" ]]; then
        log_error "pop1 and pop2 must be different populations"
        exit 1
    fi
fi

input_n=0
[[ -n "$DIVERSITY_DIR" ]] && input_n=$((input_n + 1))
[[ -n "$DIVERSITY_A" || -n "$DIVERSITY_B" ]] && input_n=$((input_n + 1))
if [[ "$input_n" -ne 1 ]]; then
    log_error "Provide --diversity-dir, or both --diversity-a and --diversity-b"
    exit 1
fi
if [[ -n "$DIVERSITY_A" || -n "$DIVERSITY_B" ]]; then
    if [[ -z "$DIVERSITY_A" || -z "$DIVERSITY_B" ]]; then
        log_error "Both --diversity-a and --diversity-b are required"
        exit 1
    fi
    if [[ "$ALL_PAIRS" == true || -n "$COMPARAND" || -n "$LRPI_COMPARISONS_FILE" ]]; then
        log_error "--diversity-a/--diversity-b is one pair; use --pop1-name and --pop2-name"
        exit 1
    fi
fi

case "$OUTPUT_FORMAT" in
    tsv|csv|hdf5) ;;
    *)
        log_error "Invalid --output-format: $OUTPUT_FORMAT (use tsv, csv, or hdf5)"
        exit 1
        ;;
esac

check_dir "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
check_dir "$LOG_DIR"
TIMESTAMP=$(date +'%Y%m%d_%H%M%S')
if [[ -n "$FILE_PREFIX" ]]; then
    LOG_FILE="${LOG_DIR}/${FILE_PREFIX}_${SCRIPT_NAME%.*}_${TIMESTAMP}.log"
else
    LOG_FILE="${LOG_DIR}/${SCRIPT_NAME%.*}_${TIMESTAMP}.log"
fi

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
else
    log_dry_run "Would create log file: $LOG_FILE"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/calculate_lrpi.R"
if [[ "$DRY_RUN" == false ]] && [[ ! -f "$R_SCRIPT" ]]; then
    log_error "R script not found: $R_SCRIPT"
    exit 1
fi
check_command "Rscript" || { log_error "Rscript not found. Please install R and Rscript."; exit 1; }

declare -a DIV_PATHS=() DIV_SAMPLES=() DIV_LABELS=()
declare -a PAIR_POP1=() PAIR_POP2=()

if [[ -n "$DIVERSITY_DIR" ]]; then
    if [[ "$DRY_RUN" == false ]]; then
        [[ -d "$DIVERSITY_DIR" ]] || { log_error "Diversity directory not found: $DIVERSITY_DIR"; exit 1; }
        discover_diversity_files "$DIVERSITY_DIR"
        [[ ${#DIV_PATHS[@]} -gt 0 ]] || { log_error "No diversity_w*_s*_sample* tables in $DIVERSITY_DIR"; exit 1; }
        log "Found ${#DIV_PATHS[@]} diversity file(s) in $DIVERSITY_DIR"
    else
        if [[ -d "$DIVERSITY_DIR" ]]; then
            discover_diversity_files "$DIVERSITY_DIR"
        fi
        if [[ ${#DIV_PATHS[@]} -eq 0 ]]; then
            DIV_PATHS=("${DIVERSITY_DIR}/diversity_w1000_s500_sampleEcho.tsv" "${DIVERSITY_DIR}/diversity_w1000_s500_sampleKjer.tsv")
            DIV_SAMPLES=("Echo" "Kjer")
            DIV_LABELS=("w1000_s500" "w1000_s500")
            log_dry_run "Would discover diversity files in $DIVERSITY_DIR"
        fi
    fi
fi

if [[ -n "$LRPI_COMPARISONS_FILE" ]]; then
    if [[ "$DRY_RUN" == false ]]; then
        check_file "$LRPI_COMPARISONS_FILE"
        parse_lrpi_comparisons "$LRPI_COMPARISONS_FILE" || exit 1
    elif [[ -f "$LRPI_COMPARISONS_FILE" ]]; then
        parse_lrpi_comparisons "$LRPI_COMPARISONS_FILE" || exit 1
    else
        PAIR_POP1=("Echo")
        PAIR_POP2=("Kjer")
        log_dry_run "Would load pairs from $LRPI_COMPARISONS_FILE"
    fi
    log "Loaded ${#PAIR_POP1[@]} directed pair(s)"
elif [[ -n "$POP1_NAME" ]]; then
    PAIR_POP1=("$POP1_NAME")
    PAIR_POP2=("$POP2_NAME")
fi

warn_unknown_samples() {
    [[ -n "$SAMPLE_INFO_CSV" && -f "$SAMPLE_INFO_CSV" && "$DRY_RUN" == false ]] || return 0
    declare -A VALID_NAMES
    local saw_header=false line sn
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        if [[ "$saw_header" == false ]]; then
            saw_header=true
            continue
        fi
        IFS=',' read -r sn _ <<< "$line"
        sn=$(echo "$sn" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        [[ -n "$sn" ]] && VALID_NAMES[$sn]=1
    done < "$SAMPLE_INFO_CSV"
    local i p
    for i in "${!PAIR_POP1[@]}"; do
        for p in "${PAIR_POP1[$i]}" "${PAIR_POP2[$i]}"; do
            [[ -z "${VALID_NAMES[$p]:-}" ]] && log_warn "'$p' is not in the sample_name column of $SAMPLE_INFO_CSV"
        done
    done
    if [[ -n "$COMPARAND" && -z "${VALID_NAMES[$COMPARAND]:-}" ]]; then
        log_warn "'$COMPARAND' is not in the sample_name column of $SAMPLE_INFO_CSV"
    fi
}

# Pairs that do not depend on the scale are checked once. comparand and all-pairs
# are expanded inside the scale loop.
if [[ ${#PAIR_POP1[@]} -gt 0 ]]; then
    warn_unknown_samples
fi

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN MODE: commands will be previewed but not executed"
    log_dry_run "Output dir: $OUTPUT_DIR, format: $OUTPUT_FORMAT"
    if [[ -n "$DIVERSITY_A" ]]; then
        log_dry_run "Would run ${POP1_NAME} / ${POP2_NAME} on $DIVERSITY_A and $DIVERSITY_B"
    else
        log_dry_run "Diversity files: ${#DIV_PATHS[@]}"
        if [[ "$ALL_PAIRS" == true ]]; then
            log_dry_run "Would run every unordered pair at each scale (pop1 < pop2)"
        elif [[ -n "$COMPARAND" ]]; then
            log_dry_run "Would run ln(pi_${COMPARAND} / pi_other) at each scale"
        else
            local_i=0
            for local_i in "${!PAIR_POP1[@]}"; do
                log_dry_run "  Pair: ${PAIR_POP1[$local_i]} / ${PAIR_POP2[$local_i]}"
            done
        fi
    fi
    exit 0
fi

run_one() {
    local file_a="$1" file_b="$2" p1="$3" p2="$4" label="$5"
    local s1 s2 out_prefix
    s1=$(sanitize_pop_name "$p1")
    s2=$(sanitize_pop_name "$p2")
    out_prefix="${FILE_PREFIX}${FILE_PREFIX:+_}${label}_${s1}_${s2}"
    local R_ARGS=(
        --diversity-a "$file_a"
        --diversity-b "$file_b"
        --pop1-name "$p1"
        --pop2-name "$p2"
        --output-dir "$OUTPUT_DIR"
        --output-format "$OUTPUT_FORMAT"
        --file-prefix "$out_prefix"
    )
    [[ "$DROP_ALL_NA" == true ]] && R_ARGS+=(--drop-all-na)
    [[ "$VERBOSE" == true ]] && R_ARGS+=(--verbose)
    log "Running LRPI for $p1 / $p2 at $label (prefix: $out_prefix)"
    [[ "$VERBOSE" == true ]] && log "  Rscript $R_SCRIPT ${R_ARGS[*]}"
    if ! Rscript "$R_SCRIPT" "${R_ARGS[@]}"; then
        log_error "R script failed for $p1 / $p2 at $label"
        exit 1
    fi
}

if [[ -n "$DIVERSITY_A" ]]; then
    check_file "$DIVERSITY_A"
    check_file "$DIVERSITY_B"
    parsed=$(parse_diversity_base "$(basename "$DIVERSITY_A")") || parsed=""
    if [[ -n "$parsed" ]]; then
        read -r w s _ <<< "$parsed"
        label="w${w}_s${s}"
    else
        label="diversity"
    fi
    run_one "$DIVERSITY_A" "$DIVERSITY_B" "$POP1_NAME" "$POP2_NAME" "$label"
    log "Done. Results in $OUTPUT_DIR"
    exit 0
fi

while IFS= read -r label; do
    [[ -n "$label" ]] || continue
    mapfile -t scale_samples < <(samples_at_label "$label")
    declare -a SCALE_P1 SCALE_P2
    SCALE_P1=()
    SCALE_P2=()
    if [[ "$ALL_PAIRS" == true ]]; then
        if [[ ${#scale_samples[@]} -lt 2 ]]; then
            log_warn "Scale $label has fewer than two samples; skipping"
            continue
        fi
        for (( i=0; i<${#scale_samples[@]}-1; i++ )); do
            for (( j=i+1; j<${#scale_samples[@]}; j++ )); do
                SCALE_P1+=("${scale_samples[$i]}")
                SCALE_P2+=("${scale_samples[$j]}")
            done
        done
    elif [[ -n "$COMPARAND" ]]; then
        found=false
        for s in "${scale_samples[@]}"; do
            if [[ "$s" == "$COMPARAND" ]]; then
                found=true
            fi
        done
        if [[ "$found" == false ]]; then
            log_error "Comparand $COMPARAND has no diversity file at $label"
            exit 1
        fi
        for s in "${scale_samples[@]}"; do
            [[ "$s" == "$COMPARAND" ]] && continue
            SCALE_P1+=("$COMPARAND")
            SCALE_P2+=("$s")
        done
    else
        for i in "${!PAIR_POP1[@]}"; do
            fa=$(file_for "${PAIR_POP1[$i]}" "$label" || true)
            fb=$(file_for "${PAIR_POP2[$i]}" "$label" || true)
            if [[ -z "$fa" || -z "$fb" ]]; then
                log_error "Missing diversity table at $label for ${PAIR_POP1[$i]} / ${PAIR_POP2[$i]}"
                exit 1
            fi
            SCALE_P1+=("${PAIR_POP1[$i]}")
            SCALE_P2+=("${PAIR_POP2[$i]}")
        done
    fi
    for i in "${!SCALE_P1[@]}"; do
        fa=$(file_for "${SCALE_P1[$i]}" "$label")
        fb=$(file_for "${SCALE_P2[$i]}" "$label")
        run_one "$fa" "$fb" "${SCALE_P1[$i]}" "${SCALE_P2[$i]}" "$label"
    done
    unset SCALE_P1 SCALE_P2
done < <(unique_labels)

log "Done. Results in $OUTPUT_DIR"
