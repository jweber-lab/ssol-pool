#!/bin/bash

###############################################################################
# plot_fst_csq.sh
#
# Wrapper for plot_fst_csq.R: per-site FST by consequence class from collate output
# (single_position.h5, or variants.h5 + fst_single.h5; collate needs --fst-dir with
# fst_single.tsv and --variant-tsv-dir/--variant-csq). All options are passed to the R script.
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_FILE="${SCRIPT_DIR}/plot_fst_csq.R"
RSCRIPT="${RSCRIPT:-Rscript}"

usage() {
    cat << EOF
Usage: $0 --hdf5-dir DIR [OPTIONS]

  --hdf5-dir DIR            Collate output directory
  --category CAT            coding_effect or region_category [default: coding_effect]
  --y-value VALUE           value or quantile [default: value]
  --sample-pairs LIST       Comma-separated pairs, e.g. Echo:Kjer,Echo:Myv [default: all]
  --chromosome CHR          Restrict to one chromosome
  --regions FILE            identify_outliers outlier_regions_*.csv: split sites inside/outside regions
  --region-stat STAT        Regions whose outlier_stat includes STAT, or 'any' [default: fst]
  --reference-category CAT  Category tested against (Wilcoxon) [default: synonymous]
  --min-sites N             Minimum sites to draw/test a class [default: 3]
  --drop-none               Drop sites with no consequence record
  --manhattan               Also plot per-site FST along the genome, coloured by consequence
  -o, --output-dir DIR      Output directory [default: ./]
  --file-prefix PREFIX      Prefix for output files
  --width N --height N      Inches [default: 12 x 7]
  --dpi N                   PNG DPI [default: 300]
  --plot-format FORMAT      png, pdf, svg, both, or all [default: png]
  --dry-run                 Print the command only
  -h, --help                Show this help

Outputs: fst_csq_<cat>.png/.csv, fst_csq_manhattan_<cat>.png (with --manhattan)

Example:
  $0 --hdf5-dir ../stats/collated --regions ../outliers/out/outlier_regions_w1000_....csv --manhattan -o fst_csq
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
