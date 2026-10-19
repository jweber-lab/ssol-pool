#!/bin/bash

###############################################################################
# outlier_annotation.sh
#
# Wrapper for outlier_annotation.R: annotate outlier regions from
#   - consequence calls (variant_csq.sh variants_csq.tsv, or collate variants.h5 /
#     single_position.h5) and the assembly GFF: sites, consequences and genes per region.
#     Output: outlier_csq_sites_<suffix>.csv, outlier_csq_genes_<suffix>.csv,
#     outlier_csq_summary_<suffix>.csv, outlier_csq_background_<suffix>.csv.
#   - BLAST against multiple DBs (from config file), mapping hits to each DB's GFF.
#     Output: outlier_genes_<suffix>.csv, outlier_genes_summary_<suffix>.csv.
#
# BLAST config: YAML or JSON with list of { db_path, name, gff_path }.
###############################################################################

set -euo pipefail

REGIONS=""
REFERENCE=""
BLAST_CONFIG=""
VARIANT_CSQ=()
GFF=""
OUTPUT_DIR="."
RSCRIPT="Rscript"
SAMTOOLS="samtools"
BLAST_CMD="blastn"
THREADS=1
PARALLEL_DBS=1
BLAST_TASK=""
EVALUE=""
MIN_QUERY_COVERAGE=""
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
  At least one of --variant-csq, --gff, --blast-config.

Consequence calls (same assembly as the regions; no BLAST):
  --variant-csq PATH     variants_csq.tsv file or directory from variant_csq.sh, or collate
                          variants.h5 / single_position.h5 (repeatable: e.g. SNP and indel csq dirs)
  --gff FILE             Assembly GFF3 given to variant_csq.sh: lists genes overlapping each region
                          and adds product/GO to the genes named by csq

BLAST:
  --blast-config FILE    BLAST config YAML or JSON (list of db_path, name, gff_path)
  --reference FILE       Reference genome FASTA (must have .fai index); required with --blast-config

Optional:
  --output-dir DIR       Output directory [default: .]
  --rscript PATH         Path to Rscript [default: Rscript]
  --samtools PATH        Path to samtools [default: samtools]
  --blast-cmd NAME       BLAST command: blastn, blastp, etc. [default: blastn]
  --threads N            Threads per BLAST run [default: 1]
  --parallel-dbs N       Run N BLAST DBs in parallel (Unix/macOS; 1=sequential) [default: 1]
  --blast-task TASK      BLAST -task (e.g. blastn for sensitive search; default megablast can miss divergent hits)
  --evalue EVAL          BLAST -evalue threshold (e.g. 10)
  --min-query-coverage F Only keep genes from subject sequences that align at least fraction F (0–1) of the query region [default: no filter]
  --verbose              Enable verbose and debug output (logged)
  -h, --help             Show this help

Consequence outputs (in --output-dir, suffix from the regions file name):
  outlier_csq_sites_*.csv       One row per (region, site): region columns, then site columns
                                (variant_type, coding_effect, genes, per-pool AD/alt_freq, FST if from HDF5)
  outlier_csq_genes_*.csv       One row per (region, gene): n_sites, n_protein_altering, coding_effects,
                                overlaps_region, GFF gene_name/product/GO/coordinates
  outlier_csq_summary_*.csv     One row per region: site counts by type, n_effect_*, n_region_*,
                                genes_csq, genes_protein_altering, genes_gff_overlap, products_gff
  outlier_csq_background_*.csv  Genome-wide counts by consequence (all sites and non-outlier sites),
                                for plotting/plot_outlier_csq.sh

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

Example:
  $0 --regions outliers/outlier_regions_w1000_....csv \\
     --variant-csq csq_snps --variant-csq csq_indels --gff annotation.gff3 \\
     --blast-config blast_config.yml --reference ref.fa --output-dir outlier_annotation
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --regions) REGIONS="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --blast-config) BLAST_CONFIG="$2"; shift 2 ;;
        --variant-csq) VARIANT_CSQ+=("$2"); shift 2 ;;
        --gff) GFF="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --rscript) RSCRIPT="$2"; shift 2 ;;
        --samtools) SAMTOOLS="$2"; shift 2 ;;
        --blast-cmd) BLAST_CMD="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --parallel-dbs) PARALLEL_DBS="$2"; shift 2 ;;
        --blast-task) BLAST_TASK="$2"; shift 2 ;;
        --evalue) EVALUE="$2"; shift 2 ;;
        --min-query-coverage) MIN_QUERY_COVERAGE="$2"; shift 2 ;;
        --verbose) VERBOSE="true"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$REGIONS" ]]; then
    echo "ERROR: --regions is required" >&2
    usage
    exit 1
fi
if [[ -z "$BLAST_CONFIG" && -z "$GFF" && ${#VARIANT_CSQ[@]} -eq 0 ]]; then
    echo "ERROR: give at least one of --variant-csq, --gff, --blast-config" >&2
    usage
    exit 1
fi
if [[ -n "$BLAST_CONFIG" && -z "$REFERENCE" ]]; then
    echo "ERROR: --reference is required with --blast-config" >&2
    exit 1
fi
for p in "${VARIANT_CSQ[@]}"; do
    if [[ "$p" == *,* ]]; then
        echo "ERROR: --variant-csq paths cannot contain commas: $p" >&2
        exit 1
    fi
done

SCRIPT_DIR=$(dirname "$0")
R_FILE="${SCRIPT_DIR}/outlier_annotation.R"
if [[ ! -f "$R_FILE" ]]; then
    echo "ERROR: outlier_annotation.R not found: $R_FILE" >&2
    exit 1
fi

R_ARGS=(--regions "$REGIONS" --output-dir "$OUTPUT_DIR")
if [[ ${#VARIANT_CSQ[@]} -gt 0 ]]; then
    R_ARGS+=(--variant-csq "$(IFS=,; echo "${VARIANT_CSQ[*]}")")
fi
[[ -n "$GFF" ]] && R_ARGS+=(--gff "$GFF")
if [[ -n "$BLAST_CONFIG" ]]; then
    R_ARGS+=(--blast-config "$BLAST_CONFIG" --reference "$REFERENCE" --samtools "$SAMTOOLS"
             --blast-cmd "$BLAST_CMD" --threads "$THREADS" --parallel-dbs "$PARALLEL_DBS")
    [[ -n "$BLAST_TASK" ]] && R_ARGS+=(--blast-task "$BLAST_TASK")
    [[ -n "$EVALUE" ]] && R_ARGS+=(--evalue "$EVALUE")
    [[ -n "$MIN_QUERY_COVERAGE" ]] && R_ARGS+=(--min-query-coverage "$MIN_QUERY_COVERAGE")
fi
[[ "$VERBOSE" == true ]] && R_ARGS+=(--verbose)

mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/outlier_annotation_${TIMESTAMP}.log"
{
    echo "Command: $RSCRIPT $R_FILE ${R_ARGS[*]}"
    echo "Started: $(date -Iseconds 2>/dev/null || date)"
} >> "$LOG_FILE"
log "Log file: $LOG_FILE"
log "Running outlier_annotation.R..."
log "  Regions: $REGIONS"
[[ ${#VARIANT_CSQ[@]} -gt 0 ]] && log "  Consequence tables: ${VARIANT_CSQ[*]}"
[[ -n "$GFF" ]] && log "  Assembly GFF: $GFF"
if [[ -n "$BLAST_CONFIG" ]]; then
    log "  Reference: $REFERENCE"
    log "  BLAST config: $BLAST_CONFIG"
    log "  BLAST command: $BLAST_CMD  threads: $THREADS  parallel-dbs: $PARALLEL_DBS"
    [[ -n "$BLAST_TASK" ]] && log "  blast-task: $BLAST_TASK"
    [[ -n "$EVALUE" ]] && log "  evalue: $EVALUE"
    [[ -n "$MIN_QUERY_COVERAGE" ]] && log "  min-query-coverage: $MIN_QUERY_COVERAGE"
fi
log "  Output dir: $OUTPUT_DIR"
[[ "$VERBOSE" == true ]] && log "  Verbose: enabled"

"$RSCRIPT" "$R_FILE" "${R_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE" || {
    log "R script failed (see $LOG_FILE)"
    exit 1
}
