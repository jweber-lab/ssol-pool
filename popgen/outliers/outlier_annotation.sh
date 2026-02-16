#!/bin/bash

###############################################################################
# outlier_annotation.sh
#
# Wrapper for outlier_annotation.R: annotate outlier regions using BLAST
# against multiple DBs (from config file). Per-region BLAST, then map hits
# to GFF for genes/GO. Output: outlier_regions_genes.csv.
#
# BLAST config: YAML or JSON with list of { db_path, name, gff_path }.
###############################################################################

set -euo pipefail

REGIONS=""
REFERENCE=""
BLAST_CONFIG=""
OUTPUT_DIR="."
RSCRIPT="Rscript"
SAMTOOLS="samtools"
BLAST_CMD="blastn"
THREADS=1
PARALLEL_DBS=1
VERBOSE="false"

# Log with timestamp (to stderr so tee still captures R stdout)
log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] $1" >&2
}

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  --regions FILE         Regions CSV/TSV (from identify_outliers: outlier_regions*.csv).
                          Extra columns (e.g. outlier_stat, n_windows, region_mean_*) are preserved in output.
  --reference FILE       Reference genome FASTA (must have .fai index)
  --blast-config FILE    BLAST config YAML or JSON (list of db_path, name, gff_path)

Optional:
  --output-dir DIR       Output directory [default: .]
  --rscript PATH         Path to Rscript [default: Rscript]
  --samtools PATH        Path to samtools [default: samtools]
  --blast-cmd NAME       BLAST command: blastn, blastp, etc. [default: blastn]
  --threads N            Threads per BLAST run [default: 1]
  --parallel-dbs N       Run N BLAST DBs in parallel (Unix/macOS; 1=sequential) [default: 1]
  --verbose              Enable verbose and debug output (logged)
  -h, --help             Show this help

Example BLAST config (YAML, blast_config.yml):
  databases:
    - db_path: /path/to/nr
      name: nr
      gff_path: /path/to/annot.gff
    - db_path: /path/to/uniprot
      name: uniprot
      gff_path: /path/to/genes.gff
      annotation_tsv: /path/to/gene_table.tsv   # optional: join by gene_id
  Output CSV includes gene_id, product, go_terms, strand, BLAST length/coords, etc.
  Missing GFF or annotation files produce warnings only.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --regions) REGIONS="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --blast-config) BLAST_CONFIG="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --rscript) RSCRIPT="$2"; shift 2 ;;
        --samtools) SAMTOOLS="$2"; shift 2 ;;
        --blast-cmd) BLAST_CMD="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --parallel-dbs) PARALLEL_DBS="$2"; shift 2 ;;
        --verbose) VERBOSE="true"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$REGIONS" || -z "$REFERENCE" || -z "$BLAST_CONFIG" ]]; then
    echo "ERROR: --regions, --reference, and --blast-config are required" >&2
    usage
    exit 1
fi

SCRIPT_DIR=$(dirname "$0")
R_FILE="${SCRIPT_DIR}/outlier_annotation.R"
if [[ ! -f "$R_FILE" ]]; then
    echo "ERROR: outlier_annotation.R not found: $R_FILE" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/outlier_annotation_${TIMESTAMP}.log"
{
    echo "Command: $RSCRIPT $R_FILE --regions $REGIONS --reference $REFERENCE --blast-config $BLAST_CONFIG --output-dir $OUTPUT_DIR --samtools $SAMTOOLS --blast-cmd $BLAST_CMD --threads $THREADS --parallel-dbs $PARALLEL_DBS${VERBOSE:+ --verbose}"
    echo "Started: $(date -Iseconds 2>/dev/null || date)"
} >> "$LOG_FILE"
log "Log file: $LOG_FILE"
log "Running outlier_annotation.R..."
log "  Regions: $REGIONS"
log "  Reference: $REFERENCE"
log "  BLAST config: $BLAST_CONFIG"
log "  Output dir: $OUTPUT_DIR"
log "  BLAST command: $BLAST_CMD  threads: $THREADS  parallel-dbs: $PARALLEL_DBS"
[[ "$VERBOSE" == true ]] && log "  Verbose: enabled"

R_EXTRA=()
[[ "$VERBOSE" == true ]] && R_EXTRA+=(--verbose)

"$RSCRIPT" "$R_FILE" \
    --regions "$REGIONS" \
    --reference "$REFERENCE" \
    --blast-config "$BLAST_CONFIG" \
    --output-dir "$OUTPUT_DIR" \
    --samtools "$SAMTOOLS" \
    --blast-cmd "$BLAST_CMD" \
    --threads "$THREADS" \
    --parallel-dbs "$PARALLEL_DBS" \
    "${R_EXTRA[@]}" 2>&1 | tee -a "$LOG_FILE" || {
    log "R script failed (see $LOG_FILE)"
    exit 1
}
