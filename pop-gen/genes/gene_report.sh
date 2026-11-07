#!/bin/bash

###############################################################################
# gene_report.sh
#
# Add genes from one or more gene tables into one output directory. Genes
# already present are left in place. Two summaries in that directory are
# rewritten from the single-site and windowed HDF5 statistics. New gene
# folders hold sequence, windows, sites, and proteome blastx hits.
#
# Usage: see --help
###############################################################################

set -euo pipefail

GENE_FILES=()
HDF5_DIR=""
REGIONS=""
SITES=""
OUTPUT_DIR="gene_report"
REFERENCE=""
GFF=""
PROTEOME_CONFIG=""
SAMTOOLS="samtools"
BLASTX_CMD=""
THREADS=""
EVALUE=""
MAX_TARGET_SEQS=""
MAX_HSPS=""
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
Usage: $0 --genes FILE [--genes FILE ...] [options]

Required:
  --genes FILE             Gene table with gene_id, chr, gene_start, gene_end.
                           Repeat to combine lists. A later run with the same
                           --output-dir adds only genes that are not already there.

Optional:
  --hdf5-dir DIR           Collated HDF5 directory. Scores every diversity, fst,
                           pbe, qibe, lrpi, and single_position file, at every
                           window size and at single sites
  --regions FILE           outlier_gene_regions_*.csv, copied into new gene folders
                           Default: sibling of each --genes file
  --sites FILE             outlier_gene_sites_*.csv, copied into new gene folders
                           Default: sibling of each --genes file
  --output-dir DIR         Output directory [default: gene_report]
  --reference FILE         Assembly FASTA, indexed. Writes gene.fa for new genes
  --gff FILE               Assembly GFF3. With --reference, writes cds.fa
  --proteome-config FILE   YAML/JSON protein databases (db_path, name, annotation_tsv, id_map).
                           Requires --reference and --gff. blastx keeps every returned HSP
  --samtools PATH          samtools [default: samtools]
  --blastx-cmd CMD         blastx program [default: blastx]
  --threads N              Threads per blastx [default: 4]
  --evalue N               blastx -evalue [default: 1e-5]
  --max-target-seqs N      blastx -max_target_seqs [default: 25]
  --max-hsps N             blastx -max_hsps [default: 5]
  --rscript PATH           Rscript [default: Rscript]
  --verbose
  --dry-run
  -h, --help

Outputs, in --output-dir:
  gene_summary_brief.csv   One row per gene already in the directory: id, chr,
                           strand, coordinates, n_snps, n_indels, and four lists
                           of stats passing the 0.99 and 0.999 high and low tails
  gene_summary_full.csv    The same genes, plus per-stat pass counts, maximum
                           overlapping value and quantile, and consequence counts
  <gene_id>/               Created only for genes that were not already present:
                           gene_stats.csv, windows.csv, sites.csv, gene.fa,
                           cds.fa, and blast_<name>.csv
                           gene_stats.csv columns: statistic, average, max, max_quantile.
                           Sequence-quality measures (coverage, mapping quality, breadth)
                           are included. The statistic is <scale>|<name>.

A stat key is <scale>|<name>, for example w1000_s500|Echo.pi or
single|Cheney:Echo.fst. High 0.99 means quantile >= 0.99. High 0.999 means
quantile >= 0.999. Low 0.99 means quantile <= 0.01. Low 0.999 means quantile
<= 0.001. SNP and indel counts come from variants.h5 or single_position.h5
when that file is present.

Example:
  $0 --genes outlier_gene_summary_fst.csv --genes outlier_gene_summary_qibe.csv \\
     --hdf5-dir ../stats/hdf5s \\
     --reference assembly.fa --gff braker.gff3 \\
     --proteome-config proteomes.yml --output-dir gene_report
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genes) GENE_FILES+=("$2"); shift 2 ;;
        --hdf5-dir) HDF5_DIR="$2"; shift 2 ;;
        --regions) REGIONS="$2"; shift 2 ;;
        --sites) SITES="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --gff) GFF="$2"; shift 2 ;;
        --proteome-config) PROTEOME_CONFIG="$2"; shift 2 ;;
        --samtools) SAMTOOLS="$2"; shift 2 ;;
        --blastx-cmd) BLASTX_CMD="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --evalue) EVALUE="$2"; shift 2 ;;
        --max-target-seqs) MAX_TARGET_SEQS="$2"; shift 2 ;;
        --max-hsps) MAX_HSPS="$2"; shift 2 ;;
        --rscript) RSCRIPT="$2"; shift 2 ;;
        --verbose) VERBOSE=true; shift ;;
        --dry-run) DRY_RUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ ${#GENE_FILES[@]} -eq 0 ]]; then
    log_error "--genes is required"
    usage
    exit 1
fi
for g in "${GENE_FILES[@]}"; do
    if [[ ! -f "$g" ]]; then
        log_error "Gene table not found: $g"
        exit 1
    fi
done
if [[ -n "$HDF5_DIR" && ! -d "$HDF5_DIR" ]]; then
    log_error "HDF5 directory not found: $HDF5_DIR"
    exit 1
fi
if [[ -n "$PROTEOME_CONFIG" && ( -z "$REFERENCE" || -z "$GFF" ) ]]; then
    log_error "--reference and --gff are required with --proteome-config"
    exit 1
fi
if [[ -n "$PROTEOME_CONFIG" && ! -f "$PROTEOME_CONFIG" ]]; then
    log_error "Proteome config not found: $PROTEOME_CONFIG"
    exit 1
fi
if [[ -n "$REFERENCE" && ! -f "$REFERENCE" ]]; then
    log_error "Reference not found: $REFERENCE"
    exit 1
fi
if [[ -n "$GFF" && ! -f "$GFF" ]]; then
    log_error "GFF not found: $GFF"
    exit 1
fi
if [[ -n "$REGIONS" && ! -f "$REGIONS" ]]; then
    log_error "Regions table not found: $REGIONS"
    exit 1
fi
if [[ -n "$SITES" && ! -f "$SITES" ]]; then
    log_error "Sites table not found: $SITES"
    exit 1
fi
if [[ -n "$THREADS" && ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    log_error "--threads must be an integer >= 1"
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
R_FILE="${SCRIPT_DIR}/gene_report.R"
if [[ ! -f "$R_FILE" ]]; then
    log_error "gene_report.R not found: $R_FILE"
    exit 1
fi
if ! command -v "$RSCRIPT" >/dev/null 2>&1; then
    log_error "Rscript not found: $RSCRIPT"
    exit 1
fi

R_ARGS=(--output-dir "$OUTPUT_DIR")
for g in "${GENE_FILES[@]}"; do
    R_ARGS+=(--genes "$g")
done
[[ -n "$HDF5_DIR" ]] && R_ARGS+=(--hdf5-dir "$HDF5_DIR")
[[ -n "$REGIONS" ]] && R_ARGS+=(--regions "$REGIONS")
[[ -n "$SITES" ]] && R_ARGS+=(--sites "$SITES")
[[ -n "$REFERENCE" ]] && R_ARGS+=(--reference "$REFERENCE")
[[ -n "$GFF" ]] && R_ARGS+=(--gff "$GFF")
[[ -n "$PROTEOME_CONFIG" ]] && R_ARGS+=(--proteome-config "$PROTEOME_CONFIG" --samtools "$SAMTOOLS")
[[ -n "$REFERENCE" && -z "$PROTEOME_CONFIG" ]] && R_ARGS+=(--samtools "$SAMTOOLS")
[[ -n "$BLASTX_CMD" ]] && R_ARGS+=(--blastx-cmd "$BLASTX_CMD")
[[ -n "$THREADS" ]] && R_ARGS+=(--threads "$THREADS")
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
LOG_FILE="${LOG_DIR}/gene_report_${TIMESTAMP}.log"
{
    echo "Command: $RSCRIPT $R_FILE ${R_ARGS[*]}"
    echo "Started: $(date -Iseconds 2>/dev/null || date)"
} >> "$LOG_FILE"
log "Log file: $LOG_FILE"
log "Running gene_report.R..."

"$RSCRIPT" "$R_FILE" "${R_ARGS[@]}" 2>&1 | tee -a "$LOG_FILE" || {
    log_error "R script failed (see $LOG_FILE)"
    exit 1
}
log "Gene report complete."
