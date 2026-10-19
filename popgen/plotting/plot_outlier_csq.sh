#!/bin/bash

###############################################################################
# plot_outlier_csq.sh
#
# Wrapper for plot_outlier_csq.R: consequence composition and enrichment of sites in
# outlier regions per outlier statistic, against the genome-wide background.
# Inputs: outlier_csq_sites_*.csv and outlier_csq_background_*.csv from
# ../outliers/outlier_annotation.sh --variant-csq. All options are passed to the R script.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_FILE="${SCRIPT_DIR}/plot_outlier_csq.R"
RSCRIPT="${RSCRIPT:-Rscript}"

usage() {
    cat << EOF
Usage: $0 (--annotation-dir DIR | --sites FILE [--background FILE]) [OPTIONS]

  --annotation-dir DIR   outlier_annotation.sh output dir (one outlier_csq_sites_*.csv)
  --sites FILE           outlier_csq_sites_*.csv
  --background FILE      outlier_csq_background_*.csv [default: same suffix as --sites]
  --category CAT         coding_effect or region_category [default: coding_effect]
  --variant-type MODE    split (SNP/indel panels), all, snp, or indel [default: split]
  --drop-none            Drop sites with no consequence record
  --min-sites N          Minimum sites per statistic for enrichment [default: 5]
  --alpha Q              BH q threshold for filled points [default: 0.05]
  -o, --output-dir DIR   Output directory [default: ./]
  --file-prefix PREFIX   Prefix for output files
  --width N --height N   Inches [default: 10 x 6]
  --dpi N                PNG DPI [default: 300]
  --plot-format FORMAT   png, pdf, svg, both, or all [default: png]
  --dry-run              Print the command only
  -h, --help             Show this help

Outputs: outlier_csq_composition_<cat>.png, outlier_csq_enrichment_<cat>.png/.csv

Example:
  $0 --annotation-dir ../outliers/outlier_annotation -o csq_plots
  $0 --annotation-dir ../outliers/outlier_annotation --category region_category --variant-type all -o csq_plots
EOF
}

ARGS=()
DRY_RUN=false
for a in "$@"; do
    case "$a" in
        -h|--help) usage; exit 0 ;;
        --dry-run) DRY_RUN=true ;;
        *) ARGS+=("$a") ;;
    esac
done
[[ ${#ARGS[@]} -eq 0 ]] && { usage; exit 1; }
[[ -f "$R_FILE" ]] || { echo "ERROR: $R_FILE not found" >&2; exit 1; }

CMD=("$RSCRIPT" "$R_FILE" "${ARGS[@]}")
if [[ "$DRY_RUN" == true ]]; then
    echo "DRY-RUN: ${CMD[*]}"
    exit 0
fi
"${CMD[@]}"
