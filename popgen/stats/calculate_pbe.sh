#!/bin/bash

###############################################################################
# calculate_pbe.sh
#
# Calculate Population Branch Statistic (PBS) and Population Branch Excess
# (PBE) from per-locus FST data calculated by calculate_fst.sh.
#
# PBS measures population-specific branch length:
# PBS1 = (T12 + T13 - T23) / 2
# where Tij = -ln(1 - Fst_ij) is the branch length between populations i and j
#
# PBE uses median-based scaling to isolate population-specific local adaptation:
# PBE1 = PBS1 - Expectation(PBS1) = PBS1 - (T23 * median(PBS1))/median(T23)
# where medians are whole-genome or chromosome-specific (see --median-scope).
#
#
# Computation is done in R (calculate_pbe.R) in a vectorized way.
# Output can be written as TSV, CSV, or HDF5.
#
# Usage: See README.md or run with --help
###############################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Default parameters
FST_FILE=""
FST_DIR=""
POP1_NAME=""
POP2_NAME=""
POP3_NAME=""
PBE_COMPARISONS_FILE=""  # CSV with pop1, pop2, pop3 columns (alternative to single --pop1/2/3-name)
SAMPLE_INFO_CSV=""       # Optional: validate comparison names against this (e.g. sample_info.csv)
COMPUTE_DUP_COMPS=false  # If set, compute PBE for duplicate/redundant comparisons (same pop1 and pop2/pop3 pair)
OUTPUT_DIR=""
FILE_PREFIX=""
OUTPUT_FORMAT="tsv"
MEDIAN_SCOPE="genome"   # genome or chromosome
DRY_RUN=false

# Colors, used in logging messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

# Global log file variable (set after argument parsing)
SCRIPT_NAME=$(basename "$0")
LOG_FILE=""
LOG_REDIRECT_ACTIVE=false  # Flag to indicate if stdout/stderr redirection is active

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Input options (choose one):
  -f, --fst-file FILE         Single FST TSV/CSV file from grenedalf/calculate_fst.sh
  -d, --fst-dir DIR           Directory containing FST output: discover all *fst*.tsv and *fst*.csv,
                              run PBS/PBE for each, and encode window/step in output names (e.g. pbe_w1000_s500.tsv)

Population comparison (choose one):
  --pop1-name NAME            Population 1 (target/focal) — single comparison
  --pop2-name NAME            Population 2 (outgroup/reference)
  --pop3-name NAME            Population 3 (outgroup/reference)
  --pbe-comparisons FILE      CSV with columns pop1, pop2, pop3 (one comparison per row).
                              Names should match samples (e.g. in sample_info.csv).
                              Duplicate/redundant rows (same pop1 and same pop2/pop3 pair,
                              order of pop2/pop3 does not matter) are skipped with a warning
                              unless --compute-dup-comps is set.

Required:
  -o, --output-dir DIR        Output directory for results

Optional:
  --sample-info FILE          Optional CSV (e.g. sample_info.csv) to validate comparison
                              names; warns if pop1/pop2/pop3 not in sample_name column
  --compute-dup-comps         Compute PBE for duplicate/redundant comparisons (default: skip)
  --file-prefix PREFIX        Prefix for output files (default: none)
  --output-format FMT         Output format: tsv, csv, or hdf5 (default: tsv)
  --median-scope SCOPE        Medians for PBE: genome (default) or chromosome
  --dry-run                   Preview only, do not run
  -h, --help                  Show this help

Output files (format depends on --output-format):
  - With --fst-file: {output_dir}/{prefix_}pbs.{tsv|csv}, {prefix_}pbe.{tsv|csv}
    (prefix derived from filename, e.g. w1000_s500 or single)
  - With --fst-dir: one PBS/PBE set per discovered FST file, e.g. {output_dir}/{prefix_}w1000_s500_pbs.{tsv|csv}
  - With --pbe-comparisons: prefix includes label and pop names, e.g. {output_dir}/{prefix_}w1000_s500_Echo_Kjer_Cheney_pbe.{tsv|csv}
  - hdf5:    {output_dir}/{prefix_}pbe.h5 (single file with PBS and PBE)
  - Logs: {output_dir}/log/{prefix_}calculate_pbe_YYYYmmdd_HHMMSS.log

Examples:

  # Single FST file (output prefix derived from filename, e.g. w1000_s500 or single):
  $0 --fst-file fst/fst_w1000_s500.tsv --pop1-name Echo --pop2-name Kjer --pop3-name Cheney -o pbe_results

  # FST directory: discover all *fst*.tsv/csv and run PBS/PBE for each (outputs e.g. w1000_s500_pbe.tsv):
  $0 --fst-dir fst/ --pop1-name Echo --pop2-name Kjer --pop3-name Cheney -o pbe_results --output-format csv

  # Multiple comparisons from CSV (pop1, pop2, pop3 columns; names match sample_info):
  $0 --fst-dir fst/ --pbe-comparisons pbe_comparisons.csv -o pbe_results --sample-info sample_info.csv

  # Outputting to HDF5 instead of tsv/csv:
  $0 --fst-file fst.tsv --pop1-name A --pop2-name B --pop3-name C -o out --output-format hdf5

Note: Requires FST data for all comparisons between the three populations
(e.g., pop1_pop2, pop1_pop3, pop2_pop3 pairs).

EOF
}

# Function to log messages (to both console and log file)
# Note: If LOG_REDIRECT_ACTIVE is true, stdout/stderr are redirected via tee, so we don't write directly to avoid duplicates
log() {
    local message="${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
    echo -e "$message"
    # Also write directly to log file (without colors) if LOG_FILE is set and redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[$(date +'%Y-%m-%d %H:%M:%S')] $1"
        echo "$clean_message" >> "$LOG_FILE"
    fi
}

log_error() {
    local message="${RED}[ERROR]${NC} $1"
    echo -e "$message" >&2
    # Also write directly to log file if redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[ERROR] $1"
        echo "$clean_message" >> "$LOG_FILE"
    fi
}

log_warn() {
    local message="${YELLOW}[WARN]${NC} $1"
    echo -e "$message"
    # Also write directly to log file if redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[WARN] $1"
        echo "$clean_message" >> "$LOG_FILE"
    fi
}

log_dry_run() {
    local message="${YELLOW}[DRY-RUN]${NC} $1"
    echo -e "$message"
    # Also write directly to log file if redirection is not active
    if [[ -n "$LOG_FILE" ]] && [[ -n "$1" ]] && [[ "$LOG_REDIRECT_ACTIVE" == false ]]; then
        local clean_message="[DRY-RUN] $1"
        echo "$clean_message" >> "$LOG_FILE"
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
    if ! command -v "$1" &>/dev/null; then
        return 1
    fi
    return 0
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

# Derive output label from FST filename: w{W}_s{S}, single, or sanitized basename
derive_fst_label() {
    local base="$1"
    base="${base%.csv}"
    base="${base%.tsv}"
    if [[ "$base" =~ fst_w([0-9]+)_s([0-9]+)$ ]]; then
        echo "w${BASH_REMATCH[1]}_s${BASH_REMATCH[2]}"
    elif [[ "$base" =~ _single$ ]] || [[ "$base" == *"fst_single" ]]; then
        echo "single"
    else
        # Fallback: strip leading path and fst, sanitize
        echo "${base}" | sed 's/^.*fst//; s/^[-_]//; s/[^a-zA-Z0-9_-]/_/g; s/^$/run/'
    fi
}

# Discover all FST files in dir; sets FST_FILES and FST_LABELS (caller must declare)
# Order: single first, then windowed by label (w1000_s500, w5000_s2500, ...)
discover_fst_files() {
    local dir="$1"
    FST_FILES=()
    FST_LABELS=()
    local f label
    while IFS= read -r -d '' f; do
        [[ -f "$f" ]] || continue
        label=$(derive_fst_label "$(basename "$f")")
        FST_FILES+=("$f")
        FST_LABELS+=("$label")
    done < <(find "$dir" -maxdepth 1 \( -name "*fst*.tsv" -o -name "*fst*.csv" \) -print0 2>/dev/null)
    # Sort: single first, then by label
    local i j tmp li lj
    for (( i=0; i < ${#FST_FILES[@]}; i++ )); do
        for (( j=i+1; j < ${#FST_FILES[@]}; j++ )); do
            li="${FST_LABELS[$i]}"; lj="${FST_LABELS[$j]}"
            swap=false
            if [[ "$lj" == "single" ]] && [[ "$li" != "single" ]]; then swap=true; fi
            if [[ "$swap" == false ]] && [[ "$li" != "single" ]] && [[ "$lj" != "single" ]] && [[ "$li" > "$lj" ]]; then swap=true; fi
            if [[ "$swap" == true ]]; then
                tmp="${FST_FILES[$i]}"; FST_FILES[$i]="${FST_FILES[$j]}"; FST_FILES[$j]="$tmp"
                tmp="${FST_LABELS[$i]}"; FST_LABELS[$i]="${FST_LABELS[$j]}"; FST_LABELS[$j]="$tmp"
            fi
        done
    done
}

# Parse pbe_comparisons CSV (columns pop1, pop2, pop3). Sets COMPARISON_POP1, COMPARISON_POP2, COMPARISON_POP3.
# Skips duplicate/redundant rows (same pop1 and same {pop2,pop3} pair) unless COMPUTE_DUP_COMPS is true.
# Returns 0 on success, 1 on error.
parse_pbe_comparisons() {
    local csv="$1"
    COMPARISON_POP1=()
    COMPARISON_POP2=()
    COMPARISON_POP3=()
    declare -A SEEN_KEYS
    local line_num=0
    local pop1_col=-1 pop2_col=-1 pop3_col=-1

    while IFS= read -r line || [[ -n "$line" ]]; do
        line_num=$((line_num + 1))
        [[ -z "$line" ]] && continue
        [[ "$line" =~ ^[[:space:]]*# ]] && continue

        # Parse CSV line (simple: no quoted commas)
        IFS=',' read -r -a F <<< "$line"
        for i in "${!F[@]}"; do
            F[$i]=$(echo "${F[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        done

        if [[ $line_num -eq 1 ]]; then
            # Header: find column indices (case-insensitive)
            for i in "${!F[@]}"; do
                local h="${F[$i]}"
                h=$(echo "$h" | tr '[:upper:]' '[:lower:]')
                case "$h" in
                    pop1) pop1_col=$i ;;
                    pop2) pop2_col=$i ;;
                    pop3) pop3_col=$i ;;
                esac
            done
            if [[ $pop1_col -lt 0 ]] || [[ $pop2_col -lt 0 ]] || [[ $pop3_col -lt 0 ]]; then
                log_error "pbe_comparisons CSV must have header columns: pop1, pop2, pop3 (found in: $line)"
                return 1
            fi
            continue
        fi

        [[ ${#F[@]} -le $pop1_col ]] || [[ ${#F[@]} -le $pop2_col ]] || [[ ${#F[@]} -le $pop3_col ]] && continue
        local p1="${F[$pop1_col]}" p2="${F[$pop2_col]}" p3="${F[$pop3_col]}"
        [[ -z "$p1" ]] || [[ -z "$p2" ]] || [[ -z "$p3" ]] && continue

        # Canonical key: pop1 + min(pop2,pop3) + max(pop2,pop3) — PBE is symmetric in pop2/pop3
        local k2 k3
        if [[ "$p2" < "$p3" ]] || [[ "$p2" == "$p3" ]]; then
            k2="$p2"; k3="$p3"
        else
            k2="$p3"; k3="$p2"
        fi
        local key="${p1}	${k2}	${k3}"

        if [[ -n "${SEEN_KEYS[$key]:-}" ]] && [[ "$COMPUTE_DUP_COMPS" != true ]]; then
            log_warn "Skipping duplicate/redundant comparison (pop2/pop3 order does not matter): pop1=$p1 pop2=$p2 pop3=$p3 (line $line_num). Use --compute-dup-comps to compute anyway."
            continue
        fi
        SEEN_KEYS[$key]=1
        COMPARISON_POP1+=("$p1")
        COMPARISON_POP2+=("$p2")
        COMPARISON_POP3+=("$p3")
    done < "$csv"

    if [[ ${#COMPARISON_POP1[@]} -eq 0 ]]; then
        log_error "No valid comparisons found in $csv (or all skipped as duplicates)"
        return 1
    fi
    return 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--fst-file)   FST_FILE="$2"; shift 2 ;;
        -d|--fst-dir)    FST_DIR="$2";  shift 2 ;;
        --pop1-name)     POP1_NAME="$2"; shift 2 ;;
        --pop2-name)     POP2_NAME="$2"; shift 2 ;;
        --pop3-name)     POP3_NAME="$2"; shift 2 ;;
        --pbe-comparisons) PBE_COMPARISONS_FILE="$2"; shift 2 ;;
        --sample-info)   SAMPLE_INFO_CSV="$2"; shift 2 ;;
        --compute-dup-comps) COMPUTE_DUP_COMPS=true; shift ;;
        -o|--output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --file-prefix)   FILE_PREFIX="$2"; shift 2 ;;
        --output-format) OUTPUT_FORMAT="$2"; shift 2 ;;
        --median-scope)  MEDIAN_SCOPE="$2"; shift 2 ;;
        --dry-run)       DRY_RUN=true; shift ;;
        -h|--help)       usage; exit 0 ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$OUTPUT_DIR" ]]; then
    log_error "Missing required argument: --output-dir"
    usage
    exit 1
fi

# Initialize log file (create output/log directory)
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

    # Redirect all stdout and stderr to log file (in addition to console)
    # Use file descriptors to preserve original stdout/stderr for console output
    exec 3>&1 4>&2
    # Tee stdout to both console (fd 3) and log file
    exec 1> >(tee -a "$LOG_FILE" >&3)
    # Redirect stderr to stdout (so it also gets teed to console and log file)
    exec 2>&1
    # Set flag so logging functions know redirection is active (avoid duplicate writes)
    LOG_REDIRECT_ACTIVE=true
else
    log_dry_run "Would create log file: $LOG_FILE"
fi

# Comparison triples: either from --pbe-comparisons file or single --pop1/2/3-name
declare -a COMPARISON_POP1 COMPARISON_POP2 COMPARISON_POP3
if [[ -n "$PBE_COMPARISONS_FILE" ]]; then
    if [[ -n "$POP1_NAME" ]] || [[ -n "$POP2_NAME" ]] || [[ -n "$POP3_NAME" ]]; then
        log_warn "Ignoring --pop1-name/--pop2-name/--pop3-name when --pbe-comparisons is set"
    fi
    if [[ "$DRY_RUN" == false ]]; then
        check_file "$PBE_COMPARISONS_FILE"
        parse_pbe_comparisons "$PBE_COMPARISONS_FILE" || exit 1
        log "Loaded ${#COMPARISON_POP1[@]} comparison(s) from $PBE_COMPARISONS_FILE"
    else
        if [[ -f "$PBE_COMPARISONS_FILE" ]]; then
            parse_pbe_comparisons "$PBE_COMPARISONS_FILE" 2>/dev/null || {
                COMPARISON_POP1=("Echo"); COMPARISON_POP2=("Kjer"); COMPARISON_POP3=("Cheney")
            }
        else
            COMPARISON_POP1=("Echo"); COMPARISON_POP2=("Kjer"); COMPARISON_POP3=("Cheney")
            log_dry_run "Would load comparisons from $PBE_COMPARISONS_FILE"
        fi
    fi
    # Optional: validate names against sample_info
    if [[ -n "$SAMPLE_INFO_CSV" ]] && [[ -f "$SAMPLE_INFO_CSV" ]] && [[ "$DRY_RUN" == false ]]; then
        declare -A VALID_NAMES
        local saw_header=false
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
        for i in "${!COMPARISON_POP1[@]}"; do
            local p1="${COMPARISON_POP1[$i]}" p2="${COMPARISON_POP2[$i]}" p3="${COMPARISON_POP3[$i]}"
            [[ -z "${VALID_NAMES[$p1]:-}" ]] && log_warn "pop1 '$p1' not in sample_info sample_name column"
            [[ -z "${VALID_NAMES[$p2]:-}" ]] && log_warn "pop2 '$p2' not in sample_info sample_name column"
            [[ -z "${VALID_NAMES[$p3]:-}" ]] && log_warn "pop3 '$p3' not in sample_info sample_name column"
        done
    fi
else
    if [[ -z "$POP1_NAME" ]] || [[ -z "$POP2_NAME" ]] || [[ -z "$POP3_NAME" ]]; then
        log_error "Missing population names: use --pop1-name, --pop2-name, --pop3-name or --pbe-comparisons FILE"
        usage
        exit 1
    fi
    COMPARISON_POP1=("$POP1_NAME")
    COMPARISON_POP2=("$POP2_NAME")
    COMPARISON_POP3=("$POP3_NAME")
fi

# Resolve FST file(s) and labels
declare -a FST_FILES FST_LABELS
if [[ -n "$FST_FILE" ]]; then
    [[ "$DRY_RUN" == false ]] && check_file "$FST_FILE"
    FST_FILES=("$FST_FILE")
    FST_LABELS=("$(derive_fst_label "$(basename "$FST_FILE")")")
    log "Using FST file: $FST_FILE (output label: ${FST_LABELS[0]})"
elif [[ -n "$FST_DIR" ]]; then
    if [[ "$DRY_RUN" == false ]]; then
        [[ -d "$FST_DIR" ]] || { log_error "FST directory not found: $FST_DIR"; exit 1; }
        discover_fst_files "$FST_DIR"
        [[ ${#FST_FILES[@]} -gt 0 ]] || { log_error "No *fst*.tsv or *fst*.csv in $FST_DIR"; exit 1; }
        log "Found ${#FST_FILES[@]} FST file(s) in $FST_DIR"
        for i in "${!FST_FILES[@]}"; do
            log "  ${FST_FILES[$i]} -> ${FST_LABELS[$i]}"
        done
    else
        discover_fst_files "$FST_DIR" 2>/dev/null || true
        if [[ ${#FST_FILES[@]} -eq 0 ]]; then
            FST_FILES=("${FST_DIR}/fst_w1000_s500.csv")
            FST_LABELS=("w1000_s500")
        fi
        log_dry_run "Would discover FST files in $FST_DIR and run for each"
    fi
else
    log_error "Provide either --fst-file or --fst-dir"
    usage
    exit 1
fi

# Validate output format and median-scope
case "$OUTPUT_FORMAT" in
    tsv|csv|hdf5) ;;
    *)
        log_error "Invalid --output-format: $OUTPUT_FORMAT (use tsv, csv, or hdf5)"
        exit 1
        ;;
esac
case "${MEDIAN_SCOPE:-genome}" in
    genome|chromosome) ;;
    *)
        log_error "Invalid --median-scope: ${MEDIAN_SCOPE:-genome} (use genome or chromosome)"
        exit 1
        ;;
esac
MEDIAN_SCOPE="${MEDIAN_SCOPE:-genome}"

# Locate R script (same directory as this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/calculate_pbe.R"
if [[ "$DRY_RUN" == false ]] && [[ ! -f "$R_SCRIPT" ]]; then
    log_error "R script not found: $R_SCRIPT"
    exit 1
fi

# Require R
check_command "Rscript" || { log_error "Rscript not found. Please install R and Rscript."; exit 1; }

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN MODE: Commands will be previewed but not executed"
    log_dry_run "Input: ${#FST_FILES[@]} FST file(s), ${#COMPARISON_POP1[@]} comparison(s)"
    for i in "${!FST_FILES[@]}"; do
        log_dry_run "  ${FST_FILES[$i]} -> ${FST_LABELS[$i]}"
    done
    for i in "${!COMPARISON_POP1[@]}"; do
        log_dry_run "  Comparison: ${COMPARISON_POP1[$i]} / ${COMPARISON_POP2[$i]} / ${COMPARISON_POP3[$i]}"
    done
    log_dry_run "Output dir: $OUTPUT_DIR, format: $OUTPUT_FORMAT"
    log_dry_run "Would run Rscript $R_SCRIPT for each (FST file x comparison) with --file-prefix <label>_<pop1>_<pop2>_<pop3>"
    exit 0
fi

# Sanitize population name for use in filenames
sanitize_pop_name() {
    echo "$1" | sed 's/[^a-zA-Z0-9_-]/_/g'
}

# Run R for each (FST file x comparison)
for fi in "${!FST_FILES[@]}"; do
    fst_path="${FST_FILES[$fi]}"
    label="${FST_LABELS[$fi]}"
    for ci in "${!COMPARISON_POP1[@]}"; do
        p1="${COMPARISON_POP1[$ci]}"
        p2="${COMPARISON_POP2[$ci]}"
        p3="${COMPARISON_POP3[$ci]}"
        s1=$(sanitize_pop_name "$p1")
        s2=$(sanitize_pop_name "$p2")
        s3=$(sanitize_pop_name "$p3")
        out_prefix="${FILE_PREFIX}${FILE_PREFIX:+_}${label}_${s1}_${s2}_${s3}"
        R_ARGS=(
            --fst-file "$fst_path"
            --pop1-name "$p1"
            --pop2-name "$p2"
            --pop3-name "$p3"
            --output-dir "$OUTPUT_DIR"
            --output-format "$OUTPUT_FORMAT"
            --file-prefix "$out_prefix"
            --median-scope "$MEDIAN_SCOPE"
        )
        log "Running PBS/PBE for $fst_path | $p1 / $p2 / $p3 (prefix: $out_prefix)"
        if ! Rscript "$R_SCRIPT" "${R_ARGS[@]}"; then
            log_error "R script failed for $fst_path ($p1 / $p2 / $p3)"
            exit 1
        fi
    done
done

log "Done. Results in $OUTPUT_DIR"
