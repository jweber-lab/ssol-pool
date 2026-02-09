#!/bin/bash

###############################################################################
# process_poolseq.sh
# 
# Pipeline for processing Illumina 2x150 paired-end poolseq data from
# Schistocephalus solidus genomic DNA.
#
# This script processes pooled sequencing data through adapter trimming,
# read merging (for overlapping reads), mapping, deduplication, and variant
# calling preparation. BAM @RG SM (sample) is set from --sample-name (or
# sample_name in sample_info.csv when using --sample-info) so downstream
# tools (e.g. variant_call.sh, bcftools) see matching sample names.
#
# Author: Based on workflow by JW
# Usage: See README.md or run with --help
###############################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Default parameters
SAMPLE_NAME=""
R1_INPUT=""
R2_INPUT=""
REFERENCE=""
OUTPUT_DIR=""
WORK_DIR=""
ADAPTERS_FA=""
BBTOOLS_DIR=""
POPOOLATION2_JAR=""
GRENEDALF=""
USE_GRENEDALF_SYNC=false
THREADS=8
MIN_OVERLAP=20
MIN_OVERLAP0=15
TRIMQ=20
MPILEUP_MIN_QUAL=20
SKIP_SYNC=false
DRY_RUN=false
PARALLEL=false
PARALLEL_MAX_JOBS=0
SAMPLE_INFO_CSV=""
USE_TMUX=""
USE_SCREEN=""
MULTIPLEXER=""
FILE_PREFIX=""  # Optional prefix for output files (default: no prefix)
LOG_DIR_PARENT=""  # Parent log directory (passed from parallel mode to child jobs)
BBTOOLS_MEMORY=""  # Maximum Java heap size for BBtools (e.g., "8g" or "16000m"), empty = auto

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required options:
  -n, --sample-name NAME      Sample name (e.g., Echo_Pool_S1)
  -1, --read1 FILE            Path to R1 fastq.gz file
  -2, --read2 FILE            Path to R2 fastq.gz file
  -r, --reference FILE        Path to reference genome (FASTA, can be .gz)
  -o, --output-dir DIR        Output directory for final results (default: ./)

Required for full pipeline:
  -a, --adapters FILE         Path to adapters.fa file (BBtools resource)
                              Optional if auto-detected relative to BBtools location
  -b, --bbtools-dir DIR       Directory containing BBtools scripts (bbduk.sh, bbmerge.sh)
                              Optional if bbduk.sh and bbmerge.sh are in PATH

Optional:
  -w, --work-dir DIR          Working directory for intermediate files (default: output-dir/work)
  -t, --threads N             Number of threads (default: 8)
  -p, --popoolation2 JAR      Path to popoolation2 mpileup2sync.jar (optional, for sync file generation)
  -g, --grenedalf PATH        Path to grenedalf executable (optional, default: check PATH)
  --use-grenedalf-sync        Use grenedalf for sync file generation instead of popoolation2
  --min-overlap N             Minimum overlap for bbmerge (default: 20)
  --min-overlap0 N            Minimum overlap0 for bbmerge (default: 15)
  --trimq N                   Quality threshold for trimming in merge step (default: 20)
  --mpileup-min-qual N        Minimum quality for mpileup (default: 20)
  --file-prefix PREFIX        Optional prefix for output files (default: no prefix)
  --bbtools-memory SIZE       Maximum Java heap size for BBtools (e.g., "8g" or "16000m")
                              Default: auto-detect available memory and allocate 80% divided by parallel jobs
  --skip-sync                 Skip sync file generation (default: false)
  --dry-run                   Preview commands without executing (dry-run mode)
  --parallel                  Enable parallel processing for multiple samples
  --sample-info FILE          CSV file with sample information (see sample-info.csv)
                              Columns: sample_name, read1_file, read2_file, pool_size
                              BAM @RG SM (sample) is set from sample_name so downstream
                              tools (e.g. variant_call.sh, bcftools) see matching names.
  --parallel-max-jobs N       Maximum concurrent jobs (default: number of CPU cores)
  --use-tmux                  Force use of tmux for parallel jobs
  --use-screen                Force use of screen for parallel jobs
  -h, --help                  Show this help message

Example:
  $0 \\
    --sample-name Echo_Pool_S1 \\
    --read1 /path/to/Echo_Pool_S1_L001_R1_001.fastq.gz \\
    --read2 /path/to/Echo_Pool_S1_L001_R2_001.fastq.gz \\
    --reference /path/to/GCA_017591395.1_ASM1759139v1_genomic.fna.gz \\
    --output-dir /path/to/output \\
    --adapters /path/to/bbmap/resources/adapters.fa \\
    --bbtools-dir /path/to/bbtools \\
    --threads 16

EOF
}

# Global log file variable (set after argument parsing)
LOG_FILE=""
LOG_REDIRECT_ACTIVE=false  # Flag to indicate if stdout/stderr redirection is active

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

# Function to execute command or print in dry-run mode
dry_run_cmd() {
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: $*"
        return 0
    else
        "$@"
    fi
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "$1 is not installed or not in PATH"
        exit 1
    fi
}

# Function to check if file exists
check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        exit 1
    fi
}

# Function to check if directory exists, create if not
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

# Function to resolve path relative to CSV file location
resolve_path_relative_to_csv() {
    local csv_file="$1"
    local file_path="$2"
    
    # If path is already absolute, return as-is
    if [[ "$file_path" =~ ^/ ]]; then
        echo "$file_path"
        return 0
    fi
    
    # Get CSV directory
    local csv_dir
    csv_dir=$(cd "$(dirname "$csv_file")" && pwd)
    
    # Resolve path relative to CSV directory
    local resolved_path
    resolved_path=$(cd "$csv_dir" && cd "$(dirname "$file_path")" 2>/dev/null && pwd)/$(basename "$file_path") 2>/dev/null || echo "${csv_dir}/${file_path}"
    
    # Normalize path (remove .. and .)
    echo "$resolved_path" | sed 's|/\./|/|g; s|/\.\./|/|g; s|^\./||'
}

# Function to get relative path from CSV location
get_relative_path_from_csv() {
    local csv_file="$1"
    local target_file="$2"
    
    # Get absolute paths
    local csv_dir
    csv_dir=$(cd "$(dirname "$csv_file")" && pwd)
    local target_abs
    target_abs=$(cd "$(dirname "$target_file")" && pwd)/$(basename "$target_file")
    
    # If target is in the same directory as CSV, return just the filename
    if [[ "$(dirname "$target_abs")" == "$csv_dir" ]]; then
        echo "$(basename "$target_abs")"
        return 0
    fi
    
    # Calculate relative path using available tools
    local rel_path
    if command -v python3 &> /dev/null; then
        rel_path=$(python3 -c "import os; print(os.path.relpath('$target_abs', '$csv_dir'))" 2>/dev/null)
    elif command -v perl &> /dev/null; then
        rel_path=$(perl -MFile::Spec -e "print File::Spec->abs2rel('$target_abs', '$csv_dir')" 2>/dev/null)
    elif command -v realpath &> /dev/null; then
        rel_path=$(realpath --relative-to="$csv_dir" "$target_abs" 2>/dev/null)
    fi
    
    # If all methods failed, try manual calculation
    if [[ -z "$rel_path" ]] || [[ "$rel_path" == "$target_abs" ]]; then
        # Manual relative path calculation
        local csv_parts=($(echo "$csv_dir" | tr '/' ' '))
        local target_parts=($(echo "$target_abs" | tr '/' ' '))
        local common_len=0
        local i=0
        
        # Find common prefix
        while [[ $i -lt ${#csv_parts[@]} ]] && [[ $i -lt ${#target_parts[@]} ]] && [[ "${csv_parts[$i]}" == "${target_parts[$i]}" ]]; do
            common_len=$((i + 1))
            i=$((i + 1))
        done
        
        # Build relative path
        local up_levels=$((${#csv_parts[@]} - common_len))
        rel_path=""
        for ((j=0; j<up_levels; j++)); do
            rel_path="${rel_path}../"
        done
        for ((j=common_len; j<${#target_parts[@]}; j++)); do
            if [[ -n "$rel_path" ]] && [[ "${rel_path: -1}" != "/" ]]; then
                rel_path="${rel_path}/"
            fi
            rel_path="${rel_path}${target_parts[$j]}"
        done
    fi
    
    echo "$rel_path"
}

# Function to update CSV file with BAM file path
# Supports file locking for parallel/safe updates
update_csv_with_bam() {
    local csv_file="$1"
    local sample_name="$2"
    local bam_file="$3"
    local use_lock="${4:-false}"  # Optional: use file locking (for parallel mode)
    
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would update CSV file $csv_file with bam_file=$bam_file for sample $sample_name"
        return 0
    fi
    
    if [[ ! -f "$csv_file" ]]; then
        log_warn "CSV file not found: $csv_file, skipping update"
        return 1
    fi
    
    # Internal function to do the actual update (called with or without lock)
    _do_update() {
        # Convert BAM path to relative path (relative to CSV file location)
        local rel_bam_file
        rel_bam_file=$(get_relative_path_from_csv "$csv_file" "$bam_file")
        
        # Create temporary file
        local tmp_file="${csv_file}.tmp"
        
        # Read CSV and update
        local sample_found=false
        local bam_file_col_index=-1  # -1 means bam_file column doesn't exist yet
        local BAM_COL_POSITION=4     # Insert bam_file as 5th column (0-indexed: 4)
        
        while IFS= read -r line || [[ -n "$line" ]]; do
            # Skip comment lines
            if [[ "$line" =~ ^# ]]; then
                echo "$line" >> "$tmp_file"
                continue
            fi
            
            # Check if this is the header line
            if [[ "$line" =~ ^[[:space:]]*sample_name ]]; then
                # Parse header to find bam_file column index
                IFS=',' read -r -a header_fields <<< "$line"
                
                # Trim whitespace and find bam_file column
                local i
                for i in "${!header_fields[@]}"; do
                    # Trim whitespace (leading and trailing)
                    header_fields[$i]=$(echo "${header_fields[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                    if [[ "${header_fields[$i]}" == "bam_file" ]]; then
                        bam_file_col_index=$i
                        break
                    fi
                done
                
                # If bam_file column doesn't exist, insert it at position 4 (5th column)
                if [[ $bam_file_col_index -eq -1 ]]; then
                    # Split header: first 4 columns, then rest
                    local header_prefix=""
                    local header_suffix=""
                    local j
                    for j in "${!header_fields[@]}"; do
                        if [[ $j -lt $BAM_COL_POSITION ]]; then
                            [[ -n "$header_prefix" ]] && header_prefix+=","
                            header_prefix+="${header_fields[$j]}"
                        elif [[ $j -ge $BAM_COL_POSITION ]]; then
                            [[ -n "$header_suffix" ]] && header_suffix+=","
                            header_suffix+="${header_fields[$j]}"
                        fi
                    done
                    # Construct new header with bam_file inserted
                    if [[ -n "$header_suffix" ]]; then
                        echo "${header_prefix},bam_file,${header_suffix}" >> "$tmp_file"
                    else
                        echo "${header_prefix},bam_file" >> "$tmp_file"
                    fi
                    bam_file_col_index=$BAM_COL_POSITION
                else
                    echo "$line" >> "$tmp_file"
                fi
                continue
            fi
            
            # Parse CSV line
            IFS=',' read -r -a fields <<< "$line"
            
            # Trim whitespace from all fields
            local i
            for i in "${!fields[@]}"; do
                fields[$i]=$(echo "${fields[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
            done
            
            # Check if this is the sample we're looking for
            if [[ "${fields[0]}" == "$sample_name" ]]; then
                sample_found=true
                
                # If bam_file column exists in header, update it; otherwise insert at position 4
                if [[ $bam_file_col_index -ge 0 ]]; then
                    # Ensure fields array is large enough
                    while [[ ${#fields[@]} -le $bam_file_col_index ]]; do
                        fields+=("")
                    done
                    # Update bam_file column
                    fields[$bam_file_col_index]="$rel_bam_file"
                else
                    # Insert bam_file at position 4 (5th column)
                    # Build new array: first 4 fields, bam_file, then remaining fields
                    local new_fields=("${fields[@]:0:$BAM_COL_POSITION}")
                    new_fields+=("$rel_bam_file")
                    if [[ ${#fields[@]} -gt $BAM_COL_POSITION ]]; then
                        new_fields+=("${fields[@]:$BAM_COL_POSITION}")
                    fi
                    fields=("${new_fields[@]}")
                    bam_file_col_index=$BAM_COL_POSITION
                fi
                
                # Reconstruct entire line with all fields preserved
                IFS=','
                echo "${fields[*]}" >> "$tmp_file"
                unset IFS
            else
                # Keep other lines as-is, but ensure consistent column count if header was updated
                if [[ $bam_file_col_index -ge 0 ]] && [[ ${#fields[@]} -le $bam_file_col_index ]]; then
                    # Header has bam_file column but this row doesn't - extend array
                    while [[ ${#fields[@]} -le $bam_file_col_index ]]; do
                        fields+=("")
                    done
                    IFS=','
                    echo "${fields[*]}" >> "$tmp_file"
                    unset IFS
                else
                    # Preserve line as-is
                    echo "$line" >> "$tmp_file"
                fi
            fi
        done < "$csv_file"
        
        # Replace original file with updated one
        mv "$tmp_file" "$csv_file"
        
        if [[ "$sample_found" == true ]]; then
            log "Updated CSV file with BAM path for sample: $sample_name"
            return 0
        else
            log_warn "Sample $sample_name not found in CSV file"
            return 1
        fi
    }
    
    # Execute update with or without file locking
    if [[ "$use_lock" == true ]]; then
        # Use file locking to prevent race conditions in parallel mode
        # Create lock file based on CSV file path
        local lock_file="${csv_file}.lock"
        
        # Try to acquire lock (wait up to 60 seconds)
        if command -v flock &> /dev/null; then
            # Use flock for file locking
            (
                exec 200>"$lock_file"
                flock -w 60 200 || {
                    log_warn "Could not acquire lock for CSV update (timeout), retrying without lock..."
                    _do_update
                    exit $?
                }
                _do_update
                exit $?
            )
            return $?
        else
            # Fallback: use mkdir as atomic lock (works on most filesystems)
            local lock_dir="${csv_file}.lockdir"
            local max_attempts=60
            local attempt=0
            
            while [[ $attempt -lt $max_attempts ]]; do
                if mkdir "$lock_dir" 2>/dev/null; then
                    # Got the lock
                    trap "rmdir '$lock_dir' 2>/dev/null" EXIT
                    _do_update
                    local result=$?
                    rmdir "$lock_dir" 2>/dev/null
                    trap - EXIT
                    return $result
                else
                    # Lock held, wait and retry
                    sleep 1
                    attempt=$((attempt + 1))
                fi
            done
            
            log_warn "Could not acquire lock for CSV update (timeout after ${max_attempts}s), proceeding without lock..."
            _do_update
            return $?
        fi
    else
        # No locking needed (single-sample mode)
        _do_update
        return $?
    fi
}

# Function to detect available multiplexer
detect_multiplexer() {
    if [[ -n "$USE_TMUX" ]] && [[ "$USE_TMUX" == true ]]; then
        if command -v tmux &> /dev/null; then
            MULTIPLEXER="tmux"
            return 0
        else
            log_error "tmux requested but not found"
            return 1
        fi
    elif [[ -n "$USE_SCREEN" ]] && [[ "$USE_SCREEN" == true ]]; then
        if command -v screen &> /dev/null; then
            MULTIPLEXER="screen"
            return 0
        else
            log_error "screen requested but not found"
            return 1
        fi
    else
        # Auto-detect: prefer tmux over screen
        if command -v tmux &> /dev/null; then
            MULTIPLEXER="tmux"
            return 0
        elif command -v screen &> /dev/null; then
            MULTIPLEXER="screen"
            return 0
        else
            MULTIPLEXER="background"
            log_warn "No multiplexer found (tmux/screen), using background jobs"
            return 0
        fi
    fi
}

# Function to get number of CPU cores
get_cpu_cores() {
    if command -v nproc &> /dev/null; then
        nproc
    elif [[ -f /proc/cpuinfo ]]; then
        grep -c ^processor /proc/cpuinfo
    elif [[ "$(uname)" == "Darwin" ]]; then
        sysctl -n hw.ncpu
    else
        echo "8"  # Default fallback
    fi
}

# Function to get available system memory in bytes
get_available_memory_bytes() {
    if [[ "$(uname)" == "Linux" ]]; then
        # Linux: use MemAvailable from /proc/meminfo (more accurate than free)
        if [[ -f /proc/meminfo ]]; then
            # Get MemAvailable in KB, convert to bytes
            mem_kb=$(grep "^MemAvailable:" /proc/meminfo | awk '{print $2}')
            if [[ -n "$mem_kb" ]]; then
                echo $((mem_kb * 1024))
                return 0
            fi
        fi
        # Fallback: use free command
        if command -v free &> /dev/null; then
            mem_kb=$(free -b | grep "^Mem:" | awk '{print $7}')  # Available memory in bytes
            if [[ -n "$mem_kb" ]]; then
                echo "$mem_kb"
                return 0
            fi
        fi
    elif [[ "$(uname)" == "Darwin" ]]; then
        # macOS: use vm_stat
        if command -v vm_stat &> /dev/null; then
            # Get page size
            page_size=$(vm_stat | head -1 | awk '{print $8}' | sed 's/\.//')
            # Get free pages
            free_pages=$(vm_stat | grep "Pages free" | awk '{print $3}' | sed 's/\.//')
            # Calculate free memory in bytes
            echo $((free_pages * page_size))
            return 0
        fi
    fi
    # Fallback: return 0 (will trigger default)
    echo "0"
}

# Function to calculate BBtools memory per job
# Takes total available memory and number of parallel jobs
# Returns memory string in format suitable for -Xmx (e.g., "8g" or "16000m")
calculate_bbtools_memory_per_job() {
    local total_mem_bytes=$1
    local num_jobs=$2
    
    if [[ -z "$total_mem_bytes" ]] || [[ "$total_mem_bytes" -eq 0 ]]; then
        # If we can't detect memory, return empty (use BBtools default)
        echo ""
        return 0
    fi
    
    if [[ -z "$num_jobs" ]] || [[ "$num_jobs" -eq 0 ]]; then
        num_jobs=1
    fi
    
    # Calculate 80% of available memory, divided by number of jobs
    # Convert to MB for easier calculation
    local total_mem_mb=$((total_mem_bytes / 1024 / 1024))
    local mem_per_job_mb=$((total_mem_mb * 80 / 100 / num_jobs))
    
    # Ensure minimum of 1GB per job
    if [[ $mem_per_job_mb -lt 1024 ]]; then
        mem_per_job_mb=1024
    fi
    
    # Convert to GB if >= 1GB, otherwise use MB
    if [[ $mem_per_job_mb -ge 1024 ]]; then
        local mem_gb=$((mem_per_job_mb / 1024))
        echo "${mem_gb}g"
    else
        echo "${mem_per_job_mb}m"
    fi
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--sample-name)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -1|--read1)
            R1_INPUT="$2"
            shift 2
            ;;
        -2|--read2)
            R2_INPUT="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -a|--adapters)
            ADAPTERS_FA="$2"
            shift 2
            ;;
        -b|--bbtools-dir)
            BBTOOLS_DIR="$2"
            shift 2
            ;;
        -w|--work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -p|--popoolation2)
            POPOOLATION2_JAR="$2"
            shift 2
            ;;
        -g|--grenedalf)
            GRENEDALF="$2"
            shift 2
            ;;
        --use-grenedalf-sync)
            USE_GRENEDALF_SYNC=true
            shift
            ;;
        --min-overlap)
            MIN_OVERLAP="$2"
            shift 2
            ;;
        --min-overlap0)
            MIN_OVERLAP0="$2"
            shift 2
            ;;
        --trimq)
            TRIMQ="$2"
            shift 2
            ;;
        --mpileup-min-qual)
            MPILEUP_MIN_QUAL="$2"
            shift 2
            ;;
        --skip-sync)
            SKIP_SYNC=true
            shift
            ;;
        --file-prefix)
            FILE_PREFIX="$2"
            shift 2
            ;;
        --bbtools-memory)
            BBTOOLS_MEMORY="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --parallel)
            PARALLEL=true
            shift
            ;;
        --sample-info)
            SAMPLE_INFO_CSV="$2"
            shift 2
            ;;
        --parallel-max-jobs)
            PARALLEL_MAX_JOBS="$2"
            shift 2
            ;;
        --use-tmux)
            USE_TMUX=true
            shift
            ;;
        --use-screen)
            USE_SCREEN=true
            shift
            ;;
        --log-dir)
            LOG_DIR_PARENT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Check if reference is indexed (do this BEFORE parallel mode to avoid race conditions)
# This check happens early so indexing is done once before splitting into parallel jobs
if [[ -n "$REFERENCE" ]]; then
    if [[ "$DRY_RUN" == false ]]; then
        # Validate reference file exists
        if [[ ! -f "$REFERENCE" ]]; then
            log_error "Reference genome file not found: $REFERENCE"
            exit 1
        fi
        
        # Check if already indexed
        if [[ ! -f "${REFERENCE}.bwt" ]] && [[ ! -f "${REFERENCE%.gz}.bwt" ]]; then
            log "Reference genome not indexed. Creating BWA index..."
            log "This will be done once before parallel jobs start (if using parallel mode)"
            bwa index "$REFERENCE"
            log "Reference genome indexing complete"
        else
            log "Reference genome already indexed"
        fi
    else
        log_dry_run "Would check if reference is indexed: $REFERENCE"
        if [[ ! -f "${REFERENCE}.bwt" ]] && [[ ! -f "${REFERENCE%.gz}.bwt" ]]; then
            log_dry_run "Would create BWA index for: $REFERENCE"
        fi
    fi
fi

# Handle parallel processing mode
if [[ "$PARALLEL" == true ]]; then
    if [[ -z "$SAMPLE_INFO_CSV" ]]; then
        log_error "Parallel mode requires --sample-info FILE"
        usage
        exit 1
    fi
    
    if [[ ! -f "$SAMPLE_INFO_CSV" ]]; then
        log_error "Sample info CSV file not found: $SAMPLE_INFO_CSV"
        exit 1
    fi
    
    # Validate required common parameters for parallel mode
    if [[ -z "$REFERENCE" ]]; then
        log_error "Parallel mode requires --reference FILE (common to all samples)"
        usage
        exit 1
    fi
    
    # Initialize log file for parallel launcher
    LOG_DIR="./log"
    check_dir "$LOG_DIR"
    
    TIMESTAMP=$(date +'%Y%m%d_%H%M%S')
    if [[ -n "$FILE_PREFIX" ]]; then
        LOG_FILE="${LOG_DIR}/${FILE_PREFIX}_parallel_${TIMESTAMP}.log"
    else
        LOG_FILE="${LOG_DIR}/parallel_${TIMESTAMP}.log"
    fi
    
    if [[ "$DRY_RUN" == false ]]; then
        touch "$LOG_FILE"
        # Write initial header to log file
        {
            echo "=========================================="
            echo "process_poolseq.sh - Parallel Launcher"
            echo "Started: $(date +'%Y-%m-%d %H:%M:%S')"
            echo "Log file: $LOG_FILE"
            echo "=========================================="
        } >> "$LOG_FILE"
        log "Parallel launcher log file: $LOG_FILE"
    else
        log_dry_run "Would create log file: $LOG_FILE"
    fi
    
    # Detect multiplexer
    detect_multiplexer || exit 1
    log "Using $MULTIPLEXER for parallel processing"
    
    # Set max jobs if not specified
    if [[ "$PARALLEL_MAX_JOBS" -eq 0 ]]; then
        PARALLEL_MAX_JOBS=$(get_cpu_cores)
    fi
    log "Maximum concurrent jobs: $PARALLEL_MAX_JOBS"
    
    # Detect available memory and calculate BBtools memory per job (unless manually specified)
    if [[ -z "$BBTOOLS_MEMORY" ]]; then
        if [[ "$DRY_RUN" == false ]]; then
            AVAILABLE_MEM_BYTES=$(get_available_memory_bytes)
            if [[ -n "$AVAILABLE_MEM_BYTES" ]] && [[ "$AVAILABLE_MEM_BYTES" -gt 0 ]]; then
                AVAILABLE_MEM_GB=$((AVAILABLE_MEM_BYTES / 1024 / 1024 / 1024))
                log "Detected available system memory: ${AVAILABLE_MEM_GB} GB"
                
                # Calculate memory per job (80% of total, divided by number of jobs)
                BBTOOLS_MEMORY=$(calculate_bbtools_memory_per_job "$AVAILABLE_MEM_BYTES" "$PARALLEL_MAX_JOBS")
                if [[ -n "$BBTOOLS_MEMORY" ]]; then
                    log "BBtools memory per job (80% of available / $PARALLEL_MAX_JOBS jobs): $BBTOOLS_MEMORY"
                else
                    log_warn "Could not calculate BBtools memory, using BBtools default (may cause memory errors)"
                fi
            else
                log_warn "Could not detect available system memory, using BBtools default (may cause memory errors)"
                BBTOOLS_MEMORY=""
            fi
        else
            log_dry_run "Would detect available memory and calculate BBtools memory per job"
            BBTOOLS_MEMORY=""
        fi
    else
        log "Using manually specified BBtools memory: $BBTOOLS_MEMORY"
    fi
    
    # Read sample list and launch jobs
    JOB_COUNT=0
    PIDS=()
    
    line_num=0
    header_skipped=false
    while IFS= read -r line || [[ -n "$line" ]]; do
        line_num=$((line_num + 1))
        
        # Skip empty lines and comments
        [[ -z "$line" ]] && continue
        [[ "$line" =~ ^# ]] && continue
        
        # Parse CSV line (simple parsing, assumes no quoted commas in values)
        IFS=',' read -r -a fields <<< "$line"
        
        # Trim whitespace from fields
        for i in "${!fields[@]}"; do
            fields[$i]=$(echo "${fields[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        done
        
        # CSV format: skip header line (check if first field is "sample_name" case-insensitively)
        if [[ "$header_skipped" == false ]]; then
            first_field_lower=$(echo "${fields[0]}" | tr '[:upper:]' '[:lower:]')
            if [[ "$first_field_lower" == "sample_name" ]]; then
                header_skipped=true
                continue
            fi
        fi
        
        # Expected columns: sample_name, read1_file, read2_file, pool_size
        if [[ ${#fields[@]} -lt 4 ]]; then
            log_warn "Skipping line $line_num: insufficient columns (expected at least 4: sample_name, read1_file, read2_file, pool_size)"
            continue
        fi
        
        sample_name="${fields[0]}"
        r1_file="${fields[1]}"
        r2_file="${fields[2]}"
        pool_size="${fields[3]}"
        
        # Resolve file paths relative to CSV location
        r1_file=$(resolve_path_relative_to_csv "$SAMPLE_INFO_CSV" "$r1_file")
        r2_file=$(resolve_path_relative_to_csv "$SAMPLE_INFO_CSV" "$r2_file")
        
        # Use OUTPUT_DIR parameter (with default ./) and create per-sample subdirectory
        if [[ -z "$OUTPUT_DIR" ]]; then
            output_dir="./${sample_name}"
        else
            output_dir="${OUTPUT_DIR}/${sample_name}"
        fi
        
        # Validate required fields
        if [[ -z "$sample_name" ]] || [[ -z "$r1_file" ]] || [[ -z "$r2_file" ]] || [[ -z "$pool_size" ]]; then
            log_warn "Skipping line $line_num: missing required fields"
            continue
        fi
        
        # Wait if we've reached max jobs
        while [[ ${#PIDS[@]} -ge "$PARALLEL_MAX_JOBS" ]]; do
            # Check for completed jobs
            NEW_PIDS=()
            for pid in "${PIDS[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    NEW_PIDS+=("$pid")
                fi
            done
            PIDS=("${NEW_PIDS[@]}")
            sleep 1
        done
        
        # Build command for this sample
        SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"
        CMD="$SCRIPT_DIR/$SCRIPT_NAME"
        CMD="$CMD --sample-name \"$sample_name\""
        CMD="$CMD --read1 \"$r1_file\""
        CMD="$CMD --read2 \"$r2_file\""
        CMD="$CMD --reference \"$REFERENCE\""
        CMD="$CMD --output-dir \"$output_dir\""
        
        # Add common arguments
        [[ -n "$ADAPTERS_FA" ]] && CMD="$CMD --adapters \"$ADAPTERS_FA\""
        [[ -n "$BBTOOLS_DIR" ]] && CMD="$CMD --bbtools-dir \"$BBTOOLS_DIR\""
        [[ -n "$WORK_DIR" ]] && CMD="$CMD --work-dir \"$WORK_DIR\""
        [[ -n "$THREADS" ]] && CMD="$CMD --threads \"$THREADS\""
        [[ -n "$POPOOLATION2_JAR" ]] && CMD="$CMD --popoolation2 \"$POPOOLATION2_JAR\""
        [[ -n "$GRENEDALF" ]] && CMD="$CMD --grenedalf \"$GRENEDALF\""
        [[ -n "$FILE_PREFIX" ]] && CMD="$CMD --file-prefix \"$FILE_PREFIX\""
        # Pass BBtools memory setting to child jobs (either auto-calculated or manually specified)
        [[ -n "$BBTOOLS_MEMORY" ]] && CMD="$CMD --bbtools-memory \"$BBTOOLS_MEMORY\""
        [[ "$USE_GRENEDALF_SYNC" == true ]] && CMD="$CMD --use-grenedalf-sync"
        [[ "$SKIP_SYNC" == true ]] && CMD="$CMD --skip-sync"
        # Pass CSV file path to each job so it can update it when BAM is created
        [[ -n "$SAMPLE_INFO_CSV" ]] && CMD="$CMD --sample-info \"$SAMPLE_INFO_CSV\""
        # Pass log directory to child jobs so all logs go to the same place
        # Convert to absolute path so it works regardless of where child job runs from
        if [[ "$DRY_RUN" == false ]]; then
            LOG_DIR_ABS=$(cd "$LOG_DIR" && pwd)
        else
            # In dry-run, use the path as-is (it will be created)
            LOG_DIR_ABS="$LOG_DIR"
        fi
        CMD="$CMD --log-dir \"$LOG_DIR_ABS\""
        [[ "$DRY_RUN" == true ]] && CMD="$CMD --dry-run"
        
        # Detect and initialize conda environment if available
        # Check for conda in common locations or PATH
        CONDA_INIT=""
        
        # First, try to detect environment from PATH (most reliable)
        ENV_NAME=""
        if [[ "$PATH" =~ /opt/conda/envs/([^/:]+) ]]; then
            ENV_NAME="${BASH_REMATCH[1]}"
        elif [[ "$PATH" =~ /([^/:]+)/envs/([^/:]+) ]]; then
            CONDA_BASE="${BASH_REMATCH[1]}"
            ENV_NAME="${BASH_REMATCH[2]}"
        elif [[ -n "$CONDA_DEFAULT_ENV" ]]; then
            ENV_NAME="$CONDA_DEFAULT_ENV"
        fi
        
        # Determine conda base directory
        CONDA_BASE=""
        if [[ -n "$CONDA_PREFIX" ]]; then
            # Extract base from prefix (e.g., /opt/conda/envs/pop-gen -> /opt/conda)
            if [[ "$CONDA_PREFIX" =~ ^(.+)/envs/[^/]+$ ]]; then
                CONDA_BASE="${BASH_REMATCH[1]}"
            fi
        fi
        
        # Try common conda locations
        if [[ -z "$CONDA_BASE" ]]; then
            if [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
                CONDA_BASE="/opt/conda"
            elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
                CONDA_BASE="$HOME/anaconda3"
            elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
                CONDA_BASE="$HOME/miniconda3"
            elif command -v conda &> /dev/null; then
                # Try to find conda base from conda command
                CONDA_BASE=$(conda info --base 2>/dev/null || echo "")
            fi
        fi
        
        # Build conda initialization command
        if [[ -n "$CONDA_BASE" ]] && [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
            CONDA_INIT="source \"$CONDA_BASE/etc/profile.d/conda.sh\""
            if [[ -n "$ENV_NAME" ]]; then
                CONDA_INIT="${CONDA_INIT} && conda activate \"$ENV_NAME\""
            fi
            CONDA_INIT="${CONDA_INIT} && "
        fi
        
        # Launch job based on multiplexer
        SESSION_NAME="poolseq_${sample_name}"
        SESSION_NAME="${SESSION_NAME//[^a-zA-Z0-9_]/_}"  # Sanitize session name
        
        log "Launching job for sample: $sample_name (session: $SESSION_NAME)"
        if [[ -n "$CONDA_INIT" ]]; then
            log "  Conda environment will be initialized for this job"
        fi
        
        if [[ "$MULTIPLEXER" == "tmux" ]]; then
            if [[ "$DRY_RUN" == true ]]; then
                log_dry_run "Would create tmux session: $SESSION_NAME"
                log_dry_run "Command: $CMD"
            else
                # Wrap command with conda initialization if available
                FULL_CMD="${CONDA_INIT}$CMD"
                tmux new-session -d -s "$SESSION_NAME" "bash -c '$FULL_CMD 2>&1; echo \"Job completed with exit code: \$?\"'"
                log "Job launched in tmux session: $SESSION_NAME"
                log "View with: tmux attach -t $SESSION_NAME"
                log "Log file will be created by child job in: ${output_dir}/log/"
            fi
        elif [[ "$MULTIPLEXER" == "screen" ]]; then
            if [[ "$DRY_RUN" == true ]]; then
                log_dry_run "Would create screen session: $SESSION_NAME"
                log_dry_run "Command: $CMD"
            else
                # Wrap command with conda initialization if available
                FULL_CMD="${CONDA_INIT}$CMD"
                screen -dmS "$SESSION_NAME" bash -c "$FULL_CMD 2>&1; echo \"Job completed with exit code: \$?\""
                log "Job launched in screen session: $SESSION_NAME"
                log "View with: screen -r $SESSION_NAME"
                log "Log file will be created by child job in: ${output_dir}/log/"
            fi
        else
            # Background job
            if [[ "$DRY_RUN" == true ]]; then
                log_dry_run "Would launch background job"
                log_dry_run "Command: $CMD"
            else
                # Wrap command with conda initialization if available
                FULL_CMD="${CONDA_INIT}$CMD"
                bash -c "$FULL_CMD" > /dev/null 2>&1 &
                PID=$!
                PIDS+=("$PID")
                log "Job launched in background (PID: $PID)"
                log "Log file will be created by child job in: ${output_dir}/log/"
            fi
        fi
        
        JOB_COUNT=$((JOB_COUNT + 1))
    done < "$SAMPLE_INFO_CSV"
    
    # Wait for all background jobs if using background mode
    if [[ "$MULTIPLEXER" == "background" ]] && [[ "$DRY_RUN" == false ]]; then
        log "Waiting for all $JOB_COUNT jobs to complete..."
        for pid in "${PIDS[@]}"; do
            wait "$pid"
        done
        log "All jobs completed"
        
        # Post-processing: Verify CSV updates and optionally update any missing BAM paths
        if [[ -n "$SAMPLE_INFO_CSV" ]] && [[ -f "$SAMPLE_INFO_CSV" ]]; then
            log "Verifying CSV file updates..."
            # Read CSV and check which samples have BAM files
            missing_bams=0
            updated_bams=0
            line_num=0
            header_skipped=false
            while IFS= read -r line || [[ -n "$line" ]]; do
                line_num=$((line_num + 1))
                [[ -z "$line" ]] && continue
                [[ "$line" =~ ^# ]] && continue
                
                IFS=',' read -r -a fields <<< "$line"
                for i in "${!fields[@]}"; do
                    fields[$i]=$(echo "${fields[$i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
                done
                
                # Skip header line (check if first field is "sample_name" case-insensitively)
                if [[ "$header_skipped" == false ]]; then
                    first_field_lower=$(echo "${fields[0]}" | tr '[:upper:]' '[:lower:]')
                    if [[ "$first_field_lower" == "sample_name" ]]; then
                        header_skipped=true
                        continue
                    fi
                fi
                
                sample_name="${fields[0]}"
                bam_file=""
                if [[ ${#fields[@]} -gt 4 ]]; then
                    bam_file="${fields[4]}"
                fi
                
                # Determine expected output directory for this sample
                sample_output_dir=""
                if [[ -z "$OUTPUT_DIR" ]]; then
                    sample_output_dir="./${sample_name}"
                else
                    sample_output_dir="${OUTPUT_DIR}/${sample_name}"
                fi
                
                # Check if BAM file path is missing or file doesn't exist
                if [[ -z "$bam_file" ]] || [[ "$bam_file" == "" ]]; then
                    # Try to find BAM file in expected output location
                    expected_bam=""
                    if [[ -n "$FILE_PREFIX" ]]; then
                        expected_bam="${sample_output_dir}/${FILE_PREFIX}_${sample_name}_All_seq.dedup.bam"
                    else
                        expected_bam="${sample_output_dir}/${sample_name}_All_seq.dedup.bam"
                    fi
                    
                    if [[ -f "$expected_bam" ]]; then
                        log "Found missing BAM file for $sample_name, updating CSV..."
                        update_csv_with_bam "$SAMPLE_INFO_CSV" "$sample_name" "$expected_bam" true
                        updated_bams=$((updated_bams + 1))
                    else
                        missing_bams=$((missing_bams + 1))
                    fi
                else
                    # BAM path is in CSV, verify file exists
                    resolved_bam=$(resolve_path_relative_to_csv "$SAMPLE_INFO_CSV" "$bam_file")
                    if [[ -f "$resolved_bam" ]]; then
                        updated_bams=$((updated_bams + 1))
                    else
                        # Path in CSV but file doesn't exist - try to find it
                        expected_bam=""
                        if [[ -n "$FILE_PREFIX" ]]; then
                            expected_bam="${sample_output_dir}/${FILE_PREFIX}_${sample_name}_All_seq.dedup.bam"
                        else
                            expected_bam="${sample_output_dir}/${sample_name}_All_seq.dedup.bam"
                        fi
                        
                        if [[ -f "$expected_bam" ]]; then
                            log "BAM path in CSV is incorrect for $sample_name, updating..."
                            update_csv_with_bam "$SAMPLE_INFO_CSV" "$sample_name" "$expected_bam" true
                            updated_bams=$((updated_bams + 1))
                        else
                            missing_bams=$((missing_bams + 1))
                        fi
                    fi
                fi
            done < "$SAMPLE_INFO_CSV"
            
            if [[ $updated_bams -gt 0 ]]; then
                log "CSV verification complete: $updated_bams sample(s) have BAM files recorded"
            fi
            if [[ $missing_bams -gt 0 ]]; then
                log_warn "CSV verification: $missing_bams sample(s) are missing BAM file paths (jobs may still be running)"
            fi
        fi
    else
        log "Launched $JOB_COUNT jobs"
        if [[ "$MULTIPLEXER" == "tmux" ]]; then
            log "List sessions with: tmux ls"
            log "Attach to session with: tmux attach -t <session_name>"
            if [[ -n "$SAMPLE_INFO_CSV" ]]; then
                log "Note: CSV file ($SAMPLE_INFO_CSV) will be updated by each job as it completes"
                log "Each job uses file locking to safely update the CSV concurrently"
            fi
        elif [[ "$MULTIPLEXER" == "screen" ]]; then
            log "List sessions with: screen -ls"
            log "Attach to session with: screen -r <session_name>"
            if [[ -n "$SAMPLE_INFO_CSV" ]]; then
                log "Note: CSV file ($SAMPLE_INFO_CSV) will be updated by each job as it completes"
                log "Each job uses file locking to safely update the CSV concurrently"
            fi
        fi
    fi
    
    exit 0
fi

# Validate required arguments for single sample mode
# Set default output directory if not provided
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="./"
fi

# Initialize log file for single-sample mode
# (Parallel mode initializes its log file earlier, before launching jobs)
# Use parent log directory if provided (from parallel mode), otherwise create in output dir
if [[ -n "$LOG_DIR_PARENT" ]]; then
    LOG_DIR="$LOG_DIR_PARENT"
    check_dir "$LOG_DIR"
else
    LOG_DIR="${OUTPUT_DIR}/log"
    check_dir "$LOG_DIR"
fi

# Generate timestamped log filename
TIMESTAMP=$(date +'%Y%m%d_%H%M%S')
if [[ -n "$FILE_PREFIX" ]]; then
    LOG_FILE="${LOG_DIR}/${FILE_PREFIX}_${SAMPLE_NAME}_${TIMESTAMP}.log"
else
    LOG_FILE="${LOG_DIR}/${SAMPLE_NAME}_${TIMESTAMP}.log"
fi

# Create log file
if [[ "$DRY_RUN" == false ]]; then
    touch "$LOG_FILE"
    # Write initial header to log file
    {
        echo "=========================================="
        echo "process_poolseq.sh - Single Sample Mode"
        echo "Sample: $SAMPLE_NAME"
        echo "Started: $(date +'%Y-%m-%d %H:%M:%S')"
        echo "Log file: $LOG_FILE"
        echo "=========================================="
    } >> "$LOG_FILE"
    log "Log file: $LOG_FILE"
    
    # Redirect all stdout and stderr to log file (in addition to console)
    # This captures all command output (from bwa, samtools, etc.) in the log file
    # Use file descriptors to preserve original stdout/stderr for console output
    exec 3>&1 4>&2
    # Tee stdout to both console (fd 3) and log file
    exec 1> >(tee -a "$LOG_FILE" >&3)
    # Redirect stderr to stdout (so it also gets teed to console and log file)
    exec 2>&1
    # Set flag so logging functions know redirection is active (avoid duplicate writes)
    LOG_REDIRECT_ACTIVE=true
fi

# Validate required arguments
if [[ -z "$SAMPLE_NAME" ]] || [[ -z "$R1_INPUT" ]] || [[ -z "$R2_INPUT" ]] || \
   [[ -z "$REFERENCE" ]]; then
    log_error "Missing required arguments"
    usage
    exit 1
fi

# Validate threads against available CPU cores (single-sample mode)
AVAILABLE_CPUS=$(get_cpu_cores)
if [[ "$THREADS" -gt "$AVAILABLE_CPUS" ]]; then
    log_warn "Requested threads ($THREADS) exceeds available CPU cores ($AVAILABLE_CPUS)"
    log_warn "Reducing to $AVAILABLE_CPUS to prevent over-subscription"
    THREADS=$AVAILABLE_CPUS
fi
log "Available CPU cores: $AVAILABLE_CPUS"
log "Using threads: $THREADS"

# Detect available memory and set BBtools memory for single-sample mode (if not manually specified)
if [[ -z "$BBTOOLS_MEMORY" ]] && [[ "$DRY_RUN" == false ]]; then
    AVAILABLE_MEM_BYTES=$(get_available_memory_bytes)
    if [[ -n "$AVAILABLE_MEM_BYTES" ]] && [[ "$AVAILABLE_MEM_BYTES" -gt 0 ]]; then
        AVAILABLE_MEM_GB=$((AVAILABLE_MEM_BYTES / 1024 / 1024 / 1024))
        log "Detected available system memory: ${AVAILABLE_MEM_GB} GB"
        
        # Use 80% of available memory for single-sample mode
        BBTOOLS_MEMORY=$(calculate_bbtools_memory_per_job "$AVAILABLE_MEM_BYTES" "1")
        if [[ -n "$BBTOOLS_MEMORY" ]]; then
            log "BBtools memory (80% of available): $BBTOOLS_MEMORY"
        fi
    fi
fi

# Set default work directory if not provided
if [[ -z "$WORK_DIR" ]]; then
    WORK_DIR="${OUTPUT_DIR}/work"
fi

# Check required commands (skip in dry-run mode for faster preview)
if [[ "$DRY_RUN" == false ]]; then
    log "Checking required software..."
    check_command "bwa"
    check_command "samtools"
    check_command "java"
else
    log_dry_run "Would check required software: bwa, samtools, java"
fi

# Check grenedalf if specified or if using grenedalf sync
if [[ -n "$GRENEDALF" ]]; then
    if [[ ! -f "$GRENEDALF" ]] && ! command -v "$GRENEDALF" &> /dev/null; then
        log_error "Grenedalf executable not found: $GRENEDALF"
        exit 1
    fi
elif [[ "$USE_GRENEDALF_SYNC" == true ]]; then
    # Try to find grenedalf in PATH
    if command -v "grenedalf" &> /dev/null; then
        GRENEDALF="grenedalf"
        log "Found grenedalf in PATH: $GRENEDALF"
    else
        log_error "Grenedalf not found in PATH and --grenedalf not specified"
        exit 1
    fi
fi

# Check required files
# Note: REFERENCE was already validated and indexed above (before parallel mode)
check_file "$R1_INPUT"
check_file "$R2_INPUT"
# REFERENCE already checked above, but validate again for single-sample mode safety
if [[ -n "$REFERENCE" ]] && [[ "$DRY_RUN" == false ]]; then
    if [[ ! -f "$REFERENCE" ]]; then
        log_error "Reference genome file not found: $REFERENCE"
        exit 1
    fi
fi

# Check BBtools and detect adapters.fa location
if [[ -n "$ADAPTERS_FA" ]]; then
    check_file "$ADAPTERS_FA"
fi

# Try to find BBtools scripts
if [[ -z "$BBTOOLS_DIR" ]]; then
    # Check if bbduk.sh and bbmerge.sh are in PATH
    if command -v bbduk.sh &> /dev/null && command -v bbmerge.sh &> /dev/null; then
        BBDUK="bbduk.sh"
        BBMERGE="bbmerge.sh"
        log "Found BBtools scripts in PATH: $BBDUK, $BBMERGE"
        
        # Try to auto-detect adapters.fa location relative to bbduk.sh
        if [[ -z "$ADAPTERS_FA" ]]; then
            # Get full path to bbduk.sh
            BBDUK_PATH=$(command -v bbduk.sh)
            BBDUK_DIR=$(dirname "$BBDUK_PATH")
            
            # Common locations for adapters.fa relative to bbduk.sh:
            # 1. ../opt/bbmap-*/resources/adapters.fa (conda style)
            # 2. ../resources/adapters.fa (direct relative)
            # 3. ../../resources/adapters.fa (alternative)
            
            # Try conda-style path first (most common)
            if [[ -d "$BBDUK_DIR/../opt" ]]; then
                # Look for bbmap-* directories
                for bbmap_dir in "$BBDUK_DIR"/../opt/bbmap-*/resources/adapters.fa; do
                    if [[ -f "$bbmap_dir" ]]; then
                        ADAPTERS_FA="$bbmap_dir"
                        log "Auto-detected adapters.fa: $ADAPTERS_FA"
                        break
                    fi
                done
            fi
            
            # Try direct relative path if not found
            if [[ -z "$ADAPTERS_FA" ]] && [[ -f "$BBDUK_DIR/../resources/adapters.fa" ]]; then
                ADAPTERS_FA="$BBDUK_DIR/../resources/adapters.fa"
                log "Auto-detected adapters.fa: $ADAPTERS_FA"
            fi
            
            # Try alternative relative path
            if [[ -z "$ADAPTERS_FA" ]] && [[ -f "$BBDUK_DIR/../../resources/adapters.fa" ]]; then
                ADAPTERS_FA="$BBDUK_DIR/../../resources/adapters.fa"
                log "Auto-detected adapters.fa: $ADAPTERS_FA"
            fi
            
            if [[ -z "$ADAPTERS_FA" ]]; then
                log_warn "Could not auto-detect adapters.fa file"
                log_warn "Searched relative to: $BBDUK_PATH"
                log_warn "Please specify --adapters FILE or ensure adapters.fa is in a standard location"
            fi
        fi
    else
        if [[ -n "$ADAPTERS_FA" ]]; then
            log_error "BBtools directory must be provided with --bbtools-dir when using adapters"
            log_error "Alternatively, ensure bbduk.sh and bbmerge.sh are in your PATH"
            exit 1
        fi
    fi
else
    # Use provided directory
    BBDUK="${BBTOOLS_DIR}/bbduk.sh"
    BBMERGE="${BBTOOLS_DIR}/bbmerge.sh"
    if [[ "$DRY_RUN" == false ]]; then
        check_file "$BBDUK"
        check_file "$BBMERGE"
    else
        log_dry_run "Would check BBtools scripts: $BBDUK, $BBMERGE"
    fi
    
    # Try to auto-detect adapters.fa in BBtools directory if not provided
    if [[ -z "$ADAPTERS_FA" ]]; then
        # Common locations within BBtools directory:
        # 1. resources/adapters.fa (most common)
        # 2. opt/bbmap-*/resources/adapters.fa (conda style)
        if [[ -f "${BBTOOLS_DIR}/resources/adapters.fa" ]]; then
            ADAPTERS_FA="${BBTOOLS_DIR}/resources/adapters.fa"
            log "Auto-detected adapters.fa: $ADAPTERS_FA"
        elif [[ -d "${BBTOOLS_DIR}/opt" ]]; then
            # Look for bbmap-* directories
            for bbmap_dir in "${BBTOOLS_DIR}"/opt/bbmap-*/resources/adapters.fa; do
                if [[ -f "$bbmap_dir" ]]; then
                    ADAPTERS_FA="$bbmap_dir"
                    log "Auto-detected adapters.fa: $ADAPTERS_FA"
                    break
                fi
            done
        fi
    fi
fi

# Validate adapters file if we have one (either provided or auto-detected)
if [[ -n "$ADAPTERS_FA" ]]; then
    if [[ "$DRY_RUN" == false ]]; then
        check_file "$ADAPTERS_FA"
    else
        log_dry_run "Would check adapters file: $ADAPTERS_FA"
    fi
fi

# Reference indexing was already done above (before parallel mode split)
# This check is redundant but kept for safety in single-sample mode
if [[ "$DRY_RUN" == false ]] && [[ -n "$REFERENCE" ]]; then
    if [[ ! -f "${REFERENCE}.bwt" ]] && [[ ! -f "${REFERENCE%.gz}.bwt" ]]; then
        log_warn "Reference genome not indexed (should have been done earlier)"
        log "Creating BWA index now..."
        bwa index "$REFERENCE"
    fi
fi

# Create directories
check_dir "$OUTPUT_DIR"
check_dir "$WORK_DIR"
CLEANED_DIR="${WORK_DIR}/cleaned"
MERGED_DIR="${WORK_DIR}/merged"
MAPPED_DIR="${WORK_DIR}/mapped"
DEDUP_DIR="${WORK_DIR}/dedup"
check_dir "$CLEANED_DIR"
check_dir "$MERGED_DIR"
check_dir "$MAPPED_DIR"
check_dir "$DEDUP_DIR"

# Set output file names (with prefix if provided)
if [[ -n "$FILE_PREFIX" ]]; then
    FINAL_BAM="${OUTPUT_DIR}/${FILE_PREFIX}_${SAMPLE_NAME}_All_seq.dedup.bam"
    MPILEUP="${OUTPUT_DIR}/${FILE_PREFIX}_${SAMPLE_NAME}_All_seq.dedup.mpileup"
    SYNC_FILE="${OUTPUT_DIR}/${FILE_PREFIX}_${SAMPLE_NAME}.sync"
else
    FINAL_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_All_seq.dedup.bam"
    MPILEUP="${OUTPUT_DIR}/${SAMPLE_NAME}_All_seq.dedup.mpileup"
    SYNC_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.sync"
fi

# Intermediate files (no prefix needed, in work directories)
R1_CLEAN="${CLEANED_DIR}/${SAMPLE_NAME}_R1_001.clean.fastq.gz"
R2_CLEAN="${CLEANED_DIR}/${SAMPLE_NAME}_R2_001.clean.fastq.gz"
R1_QUAL="${CLEANED_DIR}/${SAMPLE_NAME}_R1_001.merge_qual_correct.fq"
R2_QUAL="${CLEANED_DIR}/${SAMPLE_NAME}_R2_001.merge_qual_correct.fq"
MERGED_FQ="${MERGED_DIR}/${SAMPLE_NAME}.merge_qual_correct.fq"
UNMERGED_R1="${MERGED_DIR}/${SAMPLE_NAME}_R1_001.qual_correct.fq"
UNMERGED_R2="${MERGED_DIR}/${SAMPLE_NAME}_R2_001.qual_correct.fq"
MERGED_SAM="${MAPPED_DIR}/${SAMPLE_NAME}_merged.sam"
UNMERGED_SAM="${MAPPED_DIR}/${SAMPLE_NAME}_unmerged.sam"
MERGED_BAM="${DEDUP_DIR}/${SAMPLE_NAME}_merged.dedup.bam"
UNMERGED_BAM="${DEDUP_DIR}/${SAMPLE_NAME}_unmerged.dedup.bam"

# Start processing
if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "DRY-RUN MODE: Commands will be previewed but not executed"
    log_dry_run "Sample: $SAMPLE_NAME"
    log_dry_run "R1: $R1_INPUT"
    log_dry_run "R2: $R2_INPUT"
    log_dry_run "Reference: $REFERENCE"
    log_dry_run "Output directory: $OUTPUT_DIR"
    log_dry_run "Working directory: $WORK_DIR"
    log_dry_run "Threads: $THREADS"
    log_dry_run ""
    log_dry_run "Output files that would be created:"
    log_dry_run "  - $FINAL_BAM"
    log_dry_run "  - $MPILEUP"
    if [[ "$SKIP_SYNC" == false ]]; then
        log_dry_run "  - $SYNC_FILE"
    fi
    log_dry_run ""
fi
log "Starting poolseq processing pipeline for sample: $SAMPLE_NAME"
log "R1: $R1_INPUT"
log "R2: $R2_INPUT"
log "Reference: $REFERENCE"
log "Output directory: $OUTPUT_DIR"
log "Working directory: $WORK_DIR"
log "Threads: $THREADS"

# Step 1: Trim adapters
if [[ -n "$ADAPTERS_FA" ]]; then
    log "Step 1/7: Trimming adapters..."
    if [[ ! -f "$R1_CLEAN" ]] || [[ ! -f "$R2_CLEAN" ]]; then
        # Build BBduk command with optional memory limit
        BBDUK_CMD=("$BBDUK")
        if [[ -n "$BBTOOLS_MEMORY" ]]; then
            BBDUK_CMD+=(-Xmx"$BBTOOLS_MEMORY")
        fi
        BBDUK_CMD+=(
            in1="$R1_INPUT"
            in2="$R2_INPUT"
            out1="$R1_CLEAN"
            out2="$R2_CLEAN"
            ref="$ADAPTERS_FA"
            ktrim=r
            ftm=5
            k=23
            mink=11
            hdist=2
            tbo
            tpe
            threads="$THREADS"
        )
        dry_run_cmd "${BBDUK_CMD[@]}"
        log "Adapter trimming complete"
    else
        log "Adapter trimming already done, skipping..."
    fi
    R1_INPUT="$R1_CLEAN"
    R2_INPUT="$R2_CLEAN"
else
    log "Step 1/7: Skipping adapter trimming (no adapters file provided)"
fi

# Step 2: Quality correction and overlap assessment
if [[ -n "$ADAPTERS_FA" ]]; then
    log "Step 2/7: Quality correction and overlap assessment..."
    if [[ ! -f "$R1_QUAL" ]] || [[ ! -f "$R2_QUAL" ]]; then
        # Build BBmerge command with optional memory limit
        BBMERGE_CMD=("$BBMERGE")
        if [[ -n "$BBTOOLS_MEMORY" ]]; then
            BBMERGE_CMD+=(-Xmx"$BBTOOLS_MEMORY")
        fi
        BBMERGE_CMD+=(
            in1="$R1_INPUT"
            in2="$R2_INPUT"
            minoverlap="$MIN_OVERLAP"
            minoverlap0="$MIN_OVERLAP0"
            out1="$R1_QUAL"
            out2="$R2_QUAL"
            ecco
            mix
            threads="$THREADS"
        )
        dry_run_cmd "${BBMERGE_CMD[@]}"
        log "Quality correction complete"
    else
        log "Quality correction already done, skipping..."
    fi
    R1_INPUT="$R1_QUAL"
    R2_INPUT="$R2_QUAL"
else
    log "Step 2/7: Skipping quality correction (no BBtools provided)"
    R1_QUAL="$R1_INPUT"
    R2_QUAL="$R2_INPUT"
fi

# Step 3: Merge overlapping reads
if [[ -n "$ADAPTERS_FA" ]]; then
    log "Step 3/7: Merging overlapping reads..."
    if [[ ! -f "$MERGED_FQ" ]] || [[ ! -f "$UNMERGED_R1" ]] || [[ ! -f "$UNMERGED_R2" ]]; then
        # Build BBmerge command with optional memory limit
        BBMERGE_CMD=("$BBMERGE")
        if [[ -n "$BBTOOLS_MEMORY" ]]; then
            BBMERGE_CMD+=(-Xmx"$BBTOOLS_MEMORY")
        fi
        BBMERGE_CMD+=(
            trimq="$TRIMQ"
            in1="$R1_INPUT"
            in2="$R2_INPUT"
            minoverlap="$MIN_OVERLAP"
            minoverlap0="$MIN_OVERLAP0"
            out="$MERGED_FQ"
            outu1="$UNMERGED_R1"
            outu2="$UNMERGED_R2"
            threads="$THREADS"
        )
        dry_run_cmd "${BBMERGE_CMD[@]}"
        log "Read merging complete"
    else
        log "Read merging already done, skipping..."
    fi
else
    log "Step 3/7: Skipping read merging (no BBtools provided)"
    log_warn "Assuming reads are already processed or merging not needed"
    MERGED_FQ=""
    UNMERGED_R1="$R1_INPUT"
    UNMERGED_R2="$R2_INPUT"
fi

# Step 4: Map reads with BWA mem
log "Step 4/7: Mapping reads with BWA mem..."

# Map merged reads if they exist
if [[ -n "$MERGED_FQ" ]] && [[ -f "$MERGED_FQ" ]]; then
    if [[ ! -f "$MERGED_SAM" ]]; then
        log "Mapping merged reads..."
        if [[ "$DRY_RUN" == true ]]; then
            log_dry_run "Would execute: bwa mem -M -t $THREADS $REFERENCE $MERGED_FQ > $MERGED_SAM"
        else
            bwa mem -M -t "$THREADS" "$REFERENCE" "$MERGED_FQ" > "$MERGED_SAM"
        fi
        log "Merged reads mapping complete"
    else
        log "Merged reads already mapped, skipping..."
    fi
fi

# Map unmerged paired-end reads
if [[ ! -f "$UNMERGED_SAM" ]]; then
    log "Mapping unmerged paired-end reads..."
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: bwa mem -M -t $THREADS $REFERENCE $UNMERGED_R1 $UNMERGED_R2 > $UNMERGED_SAM"
    else
        bwa mem -M -t "$THREADS" "$REFERENCE" "$UNMERGED_R1" "$UNMERGED_R2" > "$UNMERGED_SAM"
    fi
    log "Unmerged reads mapping complete"
else
    log "Unmerged reads already mapped, skipping..."
fi

# Step 5: Process SAM files (collate, fixmate, sort, markdup)
log "Step 5/7: Processing SAM files (deduplication)..."

# Process merged SAM if it exists
if [[ -f "$MERGED_SAM" ]] && [[ ! -f "$MERGED_BAM" ]]; then
    log "Processing merged SAM..."
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: samtools collate | samtools fixmate | samtools sort | samtools markdup"
        log_dry_run "  Input: $MERGED_SAM"
        log_dry_run "  Output: $MERGED_BAM"
    else
        samtools collate -@ "$THREADS" -O -u "$MERGED_SAM" | \
            samtools fixmate -@ "$THREADS" -m -u - - | \
            samtools sort -@ "$THREADS" -u - | \
            samtools markdup -@ "$THREADS" -r - "$MERGED_BAM"
    fi
    log "Merged BAM processing complete"
fi

# Process unmerged SAM
if [[ ! -f "$UNMERGED_BAM" ]]; then
    log "Processing unmerged SAM..."
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: samtools collate | samtools fixmate | samtools sort | samtools markdup"
        log_dry_run "  Input: $UNMERGED_SAM"
        log_dry_run "  Output: $UNMERGED_BAM"
    else
        samtools collate -@ "$THREADS" -O -u "$UNMERGED_SAM" | \
            samtools fixmate -@ "$THREADS" -m -u - - | \
            samtools sort -@ "$THREADS" -u - | \
            samtools markdup -@ "$THREADS" -r - "$UNMERGED_BAM"
    fi
    log "Unmerged BAM processing complete"
fi

# Step 6: Merge BAM files
# Read group @RG SM (sample) is set from SAMPLE_NAME so that BAM sample names match
# sample_info.csv (parallel mode) or --sample-name (single mode). Downstream tools
# (e.g. variant_call.sh, bcftools) then see the same sample names.
log "Step 6/7: Merging BAM files..."
if [[ ! -f "$FINAL_BAM" ]]; then
    # Create read group header: ID and SM = SAMPLE_NAME for downstream sample naming
    RG_HEADER="@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tLB:${SAMPLE_NAME}\tPL:ILLUMINA"
    RG_FILE="${WORK_DIR}/${SAMPLE_NAME}_rg.txt"
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would create read group file: $RG_FILE"
    else
        echo -e "$RG_HEADER" > "$RG_FILE"
    fi
    
    # Merge BAMs (use -h to set header; then addreplacerg so all reads get RG:Z:SAMPLE_NAME)
    if [[ -f "$MERGED_BAM" ]]; then
        if [[ "$DRY_RUN" == true ]]; then
            log_dry_run "Would execute: samtools merge -@ $THREADS -h $RG_FILE $FINAL_BAM $MERGED_BAM $UNMERGED_BAM"
            log_dry_run "Would execute: samtools addreplacerg -r '<RG_HEADER>' -o ${FINAL_BAM}.tmp $FINAL_BAM (assign SM to all reads)"
        else
            samtools merge -@ "$THREADS" -h "$RG_FILE" "$FINAL_BAM" "$MERGED_BAM" "$UNMERGED_BAM"
            samtools addreplacerg -r "$RG_HEADER" -o "${FINAL_BAM}.tmp" "$FINAL_BAM"
            mv "${FINAL_BAM}.tmp" "$FINAL_BAM"
        fi
    else
        # If no merged BAM, just use unmerged (and set read group)
        if [[ "$DRY_RUN" == true ]]; then
            log_dry_run "Would copy: $UNMERGED_BAM to $FINAL_BAM"
            log_dry_run "Would execute: samtools addreplacerg -r '<RG_HEADER>' -o ${FINAL_BAM}.tmp $FINAL_BAM"
            log_dry_run "Would move: ${FINAL_BAM}.tmp to $FINAL_BAM"
        else
            cp "$UNMERGED_BAM" "$FINAL_BAM"
            samtools addreplacerg -r "$RG_HEADER" -o "${FINAL_BAM}.tmp" "$FINAL_BAM"
            mv "${FINAL_BAM}.tmp" "$FINAL_BAM"
        fi
    fi
    
    # Index final BAM
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: samtools index -@ $THREADS $FINAL_BAM"
    else
        samtools index -@ "$THREADS" "$FINAL_BAM"
    fi
    log "BAM merging complete"
    # For MAPQ threshold vs depth assessment, see mapq_threshold_counts.sh (e.g. run on this BAM with -i sample_info -o dir).
    # Update CSV file with BAM path if CSV file is provided
    # Use file locking if we might be running in parallel (when SAMPLE_INFO_CSV is provided,
    # it could be a parallel job updating the same CSV file)
    if [[ -n "$SAMPLE_INFO_CSV" ]]; then
        # Always use locking when SAMPLE_INFO_CSV is provided to be safe
        # (harmless overhead in single-sample mode, critical for parallel mode)
        update_csv_with_bam "$SAMPLE_INFO_CSV" "$SAMPLE_NAME" "$FINAL_BAM" true
    fi
else
    log "BAM files already merged, skipping..."
fi

# Step 7: Create mpileup
log "Step 7/7: Creating mpileup..."
if [[ ! -f "$MPILEUP" ]]; then
    if [[ "$DRY_RUN" == true ]]; then
        log_dry_run "Would execute: samtools mpileup -B -q $MPILEUP_MIN_QUAL $FINAL_BAM > $MPILEUP"
    else
        samtools mpileup -B -q "$MPILEUP_MIN_QUAL" "$FINAL_BAM" > "$MPILEUP"
    fi
    log "Mpileup creation complete"
else
    log "Mpileup already exists, skipping..."
fi

# Optional: Create sync file
if [[ "$SKIP_SYNC" == false ]]; then
    if [[ ! -f "$SYNC_FILE" ]]; then
        # Use grenedalf if requested and available
        if [[ "$USE_GRENEDALF_SYNC" == true ]] && [[ -n "$GRENEDALF" ]]; then
            log "Creating sync file with grenedalf..."
            # Grenedalf sync command (exact options to be verified)
            # Assuming: grenedalf sync --sam-path <BAM> --output <SYNC_FILE>
            if command -v "$GRENEDALF" &> /dev/null || [[ -f "$GRENEDALF" ]]; then
                dry_run_cmd "$GRENEDALF" sync \
                    --sam-path "$FINAL_BAM" \
                    --output "$SYNC_FILE" \
                    --threads "$THREADS" \
                    --min-base-qual "$MPILEUP_MIN_QUAL" || {
                    if [[ "$DRY_RUN" == false ]]; then
                        log_error "Grenedalf sync failed"
                        exit 1
                    fi
                }
                log "Sync file creation complete with grenedalf"
            else
                log_error "Grenedalf executable not found: $GRENEDALF"
                exit 1
            fi
        # Fall back to popoolation2 if available
        elif [[ -n "$POPOOLATION2_JAR" ]]; then
            log "Creating sync file with popoolation2..."
            if [[ ! -f "$POPOOLATION2_JAR" ]]; then
                log_warn "Popoolation2 JAR not found, skipping sync file creation"
            else
                dry_run_cmd java -jar "$POPOOLATION2_JAR" \
                    --input "$MPILEUP" \
                    --output "$SYNC_FILE" \
                    --fastq-type sanger \
                    --min-qual "$MPILEUP_MIN_QUAL" \
                    --threads "$THREADS"
                log "Sync file creation complete with popoolation2"
            fi
        else
            log "Skipping sync file creation (no tool specified: use --use-grenedalf-sync or --popoolation2)"
        fi
    else
        log "Sync file already exists, skipping..."
    fi
fi

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run ""
    log_dry_run "DRY-RUN complete. All checks passed."
    log_dry_run "Run without --dry-run to execute the pipeline."
else
    log "Pipeline complete!"
    log "Final outputs:"
    log "  BAM file: $FINAL_BAM"
    log "  Mpileup: $MPILEUP"
    if [[ -f "$SYNC_FILE" ]]; then
        log "  Sync file: $SYNC_FILE"
    fi
    log "Intermediate files are in: $WORK_DIR"
fi

