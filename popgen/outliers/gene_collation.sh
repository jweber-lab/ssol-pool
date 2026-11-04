#!/bin/bash

###############################################################################
# gene_collation.sh
#
# Wrapper for gene_collation.R. Reads outlier_annotation.sh CSVs and the
# assembly GFF, and writes gene-centric tables (one gene per row, plus
# gene×region and gene×site).
#
# Usage: see --help
###############################################################################

set -euo pipefail

GFF=""
ANNOTATION_DIR=""
SUFFIX=""
SITES=""
GENES=""
SUMMARY=""
BLAST_GENES=""
INCLUDE_NON_OUTLIER=false
REFERENCE=""
BLAST_CONFIG=""
SWISSPROT_DB=""
SAMTOOLS="samtools"
BLAST_CMD=""
BLASTX_CMD=""
BLAST_TASK=""
THREADS=""
PARALLEL_DBS=""
EVALUE=""
MAX_TARGET_SEQS=""
MAX_HSPS=""
OUTPUT_DIR="."
RSCRIPT="Rscript"
VERBOSE=false
DRY_RUN=false

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log() { echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" >&2; }
log_error() { echo -e "${RED}[ERROR]${NC} $1" >&2; }
log_dry_run() { echo -e "${YELLOW}[DRY-RUN]${NC} $1"; }

usage() {
    cat << EOF
Usage: $0 --gff FILE --annotation-dir DIR [options]

Required:
  --gff FILE               Assembly GFF3 (the one given to variant_csq.sh)
  --annotation-dir DIR     Directory of outlier_annotation.sh output
                           (or pass --sites / --genes / --summary instead)

Optional:
  --suffix STR             Pick one annotation run when the directory has several
  --sites FILE             outlier_csq_sites_*.csv
  --genes FILE             outlier_csq_genes_*.csv
  --summary FILE           outlier_csq_summary_*.csv
  --blast-genes FILE       outlier_genes_*.csv (BLAST homologs, joined by region)
  --include-non-outlier    Also list GFF genes that overlap no outlier region
                           (and BLAST those genes too, when --blast-config is set)
  --reference FILE         Assembly FASTA, indexed. Required with --blast-config or --swissprot-db
  --blast-config FILE      YAML/JSON databases (db_path, name, gff_path, optional annotation_tsv).
                           Each outlier gene is searched; the best HSP is named from that GFF
  --swissprot-db PREFIX    Swiss-Prot protein BLAST database. Spliced CDS are searched with blastx
  --samtools PATH          samtools [default: samtools]
  --blast-cmd CMD          BLAST program for --blast-config [default: blastn]
  --blastx-cmd CMD         blastx program for --swissprot-db [default: blastx]
  --blast-task TASK        Nucleotide BLAST -task [default: blastn]. Pass megablast to restore the old default
  --threads N              Threads per BLAST process [default: 4]
  --parallel-dbs N         Databases searched at once [default: 1]
  --evalue N               BLAST -evalue [default: 1e-5]
  --max-target-seqs N      BLAST -max_target_seqs [default: 1]
  --max-hsps N             BLAST -max_hsps [default: 1]
  --output-dir DIR         Output directory [default: .]
  --rscript PATH           Rscript [default: Rscript]
  --verbose                Log which annotation files and how many GFF genes were read
  --dry-run
  -h, --help

Outputs (suffix copied from the annotation files):
  outlier_gene_summary_*.csv   One row per GFF gene: coordinates, n_regions,
                               n_windows, qseqids, then n_sites, coding_effects, outlier_stat.
                               Region BLAST homologs stay at the end, then one gene/product/GO/evalue
                               column set per database, homolog_genes/products/go, and with
                               --swissprot-db the top blastx hit (gene, product, organism, accession, evalue)
  outlier_gene_regions_*.csv   One row per (gene, outlier region)
  outlier_gene_sites_*.csv     One row per (gene, site); sites assigned by
                               position in the gene and/or csq gene name

Example:
  $0 --gff annotation.gff3 --annotation-dir outlier_annotation --output-dir gene_collation
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gff) GFF="$2"; shift 2 ;;
        --annotation-dir) ANNOTATION_DIR="$2"; shift 2 ;;
        --suffix) SUFFIX="$2"; shift 2 ;;
        --sites) SITES="$2"; shift 2 ;;
        --genes) GENES="$2"; shift 2 ;;
        --summary) SUMMARY="$2"; shift 2 ;;
        --blast-genes) BLAST_GENES="$2"; shift 2 ;;
        --include-non-outlier) INCLUDE_NON_OUTLIER=true; shift ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --blast-config) BLAST_CONFIG="$2"; shift 2 ;;
        --swissprot-db) SWISSPROT_DB="$2"; shift 2 ;;
        --samtools) SAMTOOLS="$2"; shift 2 ;;
        --blast-cmd) BLAST_CMD="$2"; shift 2 ;;
        --blastx-cmd) BLASTX_CMD="$2"; shift 2 ;;
        --blast-task) BLAST_TASK="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --parallel-dbs) PARALLEL_DBS="$2"; shift 2 ;;
        --evalue) EVALUE="$2"; shift 2 ;;
        --max-target-seqs) MAX_TARGET_SEQS="$2"; shift 2 ;;
        --max-hsps) MAX_HSPS="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --rscript) RSCRIPT="$2"; shift 2 ;;
        --verbose) VERBOSE=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$GFF" ]]; then
    log_error "--gff is required"
    usage
    exit 1
fi
if [[ -z "$ANNOTATION_DIR" && -z "$SITES" && -z "$GENES" && -z "$SUMMARY" ]]; then
    log_error "Give --annotation-dir or at least one of --sites, --genes, --summary"
    usage
    exit 1
fi
if [[ ! -f "$GFF" ]]; then
    log_error "GFF not found: $GFF"
    exit 1
fi
if [[ -n "$ANNOTATION_DIR" && ! -d "$ANNOTATION_DIR" ]]; then
    log_error "Annotation directory not found: $ANNOTATION_DIR"
    exit 1
fi
if [[ ( -n "$BLAST_CONFIG" || -n "$SWISSPROT_DB" ) && -z "$REFERENCE" ]]; then
    log_error "--reference is required with --blast-config or --swissprot-db"
    exit 1
fi
if [[ -n "$SWISSPROT_DB" && ! -e "${SWISSPROT_DB}.phr" && ! -e "${SWISSPROT_DB}.pin" && ! -e "${SWISSPROT_DB}.psq" && ! -e "${SWISSPROT_DB}.pdb" ]]; then
    log_error "Swiss-Prot BLAST database not found at prefix: $SWISSPROT_DB"
    exit 1
fi
if [[ -n "$BLAST_CONFIG" && ! -f "$BLAST_CONFIG" ]]; then
    log_error "BLAST config not found: $BLAST_CONFIG"
    exit 1
fi
if [[ -n "$REFERENCE" && ! -f "$REFERENCE" ]]; then
    log_error "Reference not found: $REFERENCE"
    exit 1
fi
if [[ -n "$THREADS" && ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    log_error "--threads must be an integer >= 1"
    exit 1
fi
if [[ -n "$PARALLEL_DBS" && ! "$PARALLEL_DBS" =~ ^[1-9][0-9]*$ ]]; then
    log_error "--parallel-dbs must be an integer >= 1"
    exit 1
fi
if [[ -n "$MAX_TARGET_SEQS" && ! "$MAX_TARGET_SEQS" =~ ^[1-9][0-9]*$ ]]; then
    log_error "--max-target-seqs must be an integer >= 1"
    exit 1
fi
if [[ -n "$MAX_HSPS" && ! "$MAX_HSPS" =~ ^[1-9][0-9]*$ ]]; then
    log_error "--max-hsps must be an integer >= 1"
    exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
R_FILE="${SCRIPT_DIR}/gene_collation.R"
if [[ ! -f "$R_FILE" ]]; then
    log_error "gene_collation.R not found: $R_FILE"
    exit 1
fi
if ! command -v "$RSCRIPT" >/dev/null 2>&1; then
    log_error "Rscript not found: $RSCRIPT"
    exit 1
fi

R_ARGS=(--gff "$GFF" --output-dir "$OUTPUT_DIR")
[[ -n "$ANNOTATION_DIR" ]] && R_ARGS+=(--annotation-dir "$ANNOTATION_DIR")
[[ -n "$SUFFIX" ]] && R_ARGS+=(--suffix "$SUFFIX")
[[ -n "$SITES" ]] && R_ARGS+=(--sites "$SITES")
[[ -n "$GENES" ]] && R_ARGS+=(--genes "$GENES")
[[ -n "$SUMMARY" ]] && R_ARGS+=(--summary "$SUMMARY")
[[ -n "$BLAST_GENES" ]] && R_ARGS+=(--blast-genes "$BLAST_GENES")
[[ "$INCLUDE_NON_OUTLIER" == true ]] && R_ARGS+=(--include-non-outlier)
[[ -n "$REFERENCE" ]] && R_ARGS+=(--reference "$REFERENCE")
[[ -n "$BLAST_CONFIG" ]] && R_ARGS+=(--blast-config "$BLAST_CONFIG")
[[ -n "$SWISSPROT_DB" ]] && R_ARGS+=(--swissprot-db "$SWISSPROT_DB")
[[ -n "$BLAST_CONFIG" || -n "$SWISSPROT_DB" ]] && R_ARGS+=(--samtools "$SAMTOOLS")
[[ -n "$BLAST_CMD" ]] && R_ARGS+=(--blast-cmd "$BLAST_CMD")
[[ -n "$BLASTX_CMD" ]] && R_ARGS+=(--blastx-cmd "$BLASTX_CMD")
[[ -n "$BLAST_TASK" ]] && R_ARGS+=(--blast-task "$BLAST_TASK")
[[ -n "$THREADS" ]] && R_ARGS+=(--threads "$THREADS")
[[ -n "$PARALLEL_DBS" ]] && R_ARGS+=(--parallel-dbs "$PARALLEL_DBS")
[[ -n "$EVALUE" ]] && R_ARGS+=(--evalue "$EVALUE")
[[ -n "$MAX_TARGET_SEQS" ]] && R_ARGS+=(--max-target-seqs "$MAX_TARGET_SEQS")
[[ -n "$MAX_HSPS" ]] && R_ARGS+=(--max-hsps "$MAX_HSPS")
[[ "$VERBOSE" == true ]] && R_ARGS+=(--verbose)

if [[ "$DRY_RUN" == true ]]; then
    log_dry_run "Would run: $RSCRIPT $R_FILE ${R_ARGS[*]}"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/log"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/gene_collation_${TIMESTAMP}.log"
{
    echo "Command: $RSCRIPT $R_FILE ${R_ARGS[*]}"
    echo "Started: $(date -Iseconds 2>/dev/null || date)"
} >> "$LOG_FILE"
log "Log file: $LOG_FILE"
log "Running gene_collation.R..."
log "  GFF: $GFF"
[[ -n "$ANNOTATION_DIR" ]] && log "  Annotation dir: $ANNOTATION_DIR"
log "  Output dir: $OUTPUT_DIR"

"$RSCRIPT" "$R_FILE" "${R_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE" || {
    log_error "R script failed (see $LOG_FILE)"
    exit 1
}

for f in "$OUTPUT_DIR"/outlier_gene_*.csv; do
    if [[ -f "$f" ]]; then
        bn=$(basename "$f")
        { echo ""; log "First rows of $bn:"; head -n 4 "$f"; } | tee -a "$LOG_FILE"
    fi
done
log "Gene collation complete."
