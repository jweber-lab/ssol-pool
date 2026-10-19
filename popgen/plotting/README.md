# Plotting Scripts for Pool-Seq Analysis

This directory contains scripts for visualizing population genetics statistics from pool-seq data analysis.

## Overview

Plotting workflows available:

1. **FST Cathedral Plots** (`plot_fst_cathedral.sh`): Multi-scale FST cathedral plots using grenedalf
2. **Diversity Statistics** (`plot_diversity.R` / `plot_diversity.sh`): π, θ, and Tajima's D
3. **Region Stack** (`plot_region.R` / `plot_region.sh`): One stacked figure per region (coverage, π, FST, PBE)
4. **FST Statistics** (`plot_fst.R` / `plot_fst.sh`): FST plots with paneling and overlay options
5. **PBE Statistics** (`plot_pbe.R` / `plot_pbe.sh`): PBE plots with paneling and overlay options
6. **Distributions** (`plot_hist.R` / `plot_hist.sh`): Histogram / density / ECDF of one statistic (coverage, mapq, n_snps, pi, theta, tajima_d, fst, pbe)
7. **Outlier consequences** (`plot_outlier_csq.R` / `plot_outlier_csq.sh`): Consequence composition and enrichment of sites in outlier regions, per outlier statistic
8. **Per-site FST vs consequence** (`plot_fst_csq.R` / `plot_fst_csq.sh`): Per-site FST by consequence class, optionally inside vs outside outlier regions, plus a genome-wide view

**Collate HDF5 input:** the window plots (2–6) read `diversity_w*.h5`, `fst_w*.h5`, `pbe_w*.h5`. Per-site `fst_single.h5` / `pbe_single.h5` are skipped when windowed files exist; per-site FST is plotted by `plot_fst_csq` (from `single_position.h5` or `variants.h5` + `fst_single.h5`). String columns from collate (`genes`, `genes_protein_altering`, `ref`, `alt`, `coding_effect`, `<sample>_AD`) are read as text and column names such as `Echo:Kjer.fst` keep their form.

## Scripts

### 1. plot_fst_cathedral.sh

Generates FST cathedral plots using `grenedalf fst-cathedral` and `grenedalf cathedral-plot`. Cathedral plots visualize FST across different window sizes, showing both broad and fine-scale differentiation patterns.

#### Usage

```bash
./plot_fst_cathedral.sh \
  --sample-info ../sample-info.csv \
  --output-dir ./cathedral_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --cathedral-width 2000 \
  --cathedral-height 600
```

#### Key Options

- `--sample-info FILE`: CSV file with sample information (columns: sample_name, read1_file, read2_file, pool_size, bam_file)
- `--bam FILE`: BAM files (can specify multiple times, minimum 2 required)
- `--output-dir DIR`: Output directory for results (required)
- `--reference-genome FILE`: Reference genome FASTA file (recommended)
- `--comparand SAMPLE`: Compute FST between this sample and all others
- `--cathedral-width N`: Plot width in pixels (default: 1500)
- `--cathedral-height N`: Plot height in pixels (default: 500)
- `--method METHOD`: FST method: unbiased-nei or unbiased-hudson (default: unbiased-nei)
- `--file-prefix PREFIX`: Optional prefix for output files
- `--dry-run`: Preview commands without executing

#### Output Files

- `*.csv`: Per-pixel FST value matrices (from `grenedalf fst-cathedral`)
- `*.json`: Plot metadata (from `grenedalf fst-cathedral`)
- `*.bmp`: Bitmap plots (from `grenedalf cathedral-plot`)
- `*.svg`: Vector plots with axes and legend (from `grenedalf cathedral-plot`)
- Files are named per chromosome and sample pair

#### Example

```bash
# Generate cathedral plots for all pairwise comparisons
./plot_fst_cathedral.sh \
  --sample-info ../sample-info.csv \
  --output-dir ./cathedral_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --pool-size 200 \
  --threads 8

# Generate cathedral plot for specific sample pair
./plot_fst_cathedral.sh \
  --sample-info ../sample-info.csv \
  --output-dir ./cathedral_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --comparand Echo \
  --cathedral-width 2000 \
  --cathedral-height 600
```

---

### 2. plot_diversity.R / plot_diversity.sh

Creates ggplot2 plots of diversity statistics (π, θ, Tajima's D) with chromosome stripes, paneling options, and summary statistics.

#### Usage

**Using the wrapper script (recommended):**
```bash
./plot_diversity.sh \
  --input-dir ../stats/diversity_output \
  --output-dir ./diversity_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --panel-by both \
  --plot-format both
```

**Direct R script usage:**
```bash
Rscript plot_diversity.R \
  --input-dir ../stats/diversity_output \
  --output-dir ./diversity_plots \
  --panel-by both \
  --statistics pi,theta,tajima_d
```

#### Key Options

- `--input-dir DIR`: Directory containing diversity TSV or HDF5 files (from `calculate_pi_theta.sh` or collate) (required)
- `--input-format FORMAT`: `tsv`, `hdf5`, or `auto` (default: auto)
- `--y-value VALUE`: Y-axis: `value`, `rank`, or `quantile` (default: value; rank/quantile require HDF5)
- `--plot-style STYLE`: `line` or `line_points` (default: line)
- `--overlay-samples`: Plot all samples in one panel per window (color = sample)
- `--output-dir DIR`: Output directory for plots (default: ./)
- `--reference-genome FILE`: Reference genome FASTA file (optional, for chromosome lengths)
- `--panel-by OPTION`: Paneling option: `window`, `sample`, `both`, or `none` (default: `both`)
- `--statistics LIST`: Comma-separated list: `pi`, `theta`, `tajima_d` (default: all)
- `--plot-format FORMAT`: Plot format: `png`, `pdf`, or `both` (default: `png`)
- `--width N`, `--height N`: Plot dimensions in inches (default: 12, 8)
- `--dpi N`: DPI for PNG (default: 300)
- `--file-prefix PREFIX`: Prefix for output files (default: no prefix)
- `--dry-run`: Preview commands without executing (wrapper only)

#### Output Files

- `diversity_{statistic}.png` (or `.pdf`): Plots for each statistic (π, θ, Tajima's D)
- `diversity_statistics.csv`: Summary statistics with columns:
  - `statistic`: pi, theta, or tajima_d
  - `sample`: Sample name
  - `window_size`: Window size in bp (or NA if not applicable)
  - `chromosome`: Chromosome name or "genome" for genome-wide stats
  - `median`: Median value
  - `mean`: Mean value
  - `q95`: 95th percentile value

#### Features

- **Chromosome Stripes**: Alternating white/grey stripes for each chromosome
- **Paneling**: Options to panel by window size, sample name, or both
- **Reference Lines**: Horizontal dashed lines showing:
  - Blue: Median (genome-wide)
  - Red: Mean (genome-wide)
  - Orange: 95th percentile (genome-wide)
- **Statistics**: Calculated both genome-wide and per-chromosome

#### Example

```bash
# Plot all statistics with paneling by sample and window size
./plot_diversity.sh \
  --input-dir ../stats/diversity_output \
  --output-dir ./diversity_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --panel-by both \
  --plot-format both

# Plot only π and θ, paneled by window size
./plot_diversity.sh \
  --input-dir ../stats/diversity_output \
  --output-dir ./diversity_plots \
  --statistics pi,theta \
  --panel-by window \
  --plot-format png
```

---

### 3. plot_region.R / plot_region.sh

Creates a **single stacked figure** for a genomic region (or full chromosome) with panels for coverage, π, FST, and PBE, shared x-axis, and a shared color key for sample/pair/trio.

#### Usage

```bash
./plot_region.sh \
  --chromosome chr1 \
  --hdf5-dir ../stats/collate_out \
  --reference-genome ../ref.fa \
  --output-dir ./region_plots
```

Or a specific region and TSV inputs:

```bash
./plot_region.sh \
  --region chr2:5000000-6000000 \
  --diversity-dir ../stats/diversity \
  --fst-dir ../stats/fst \
  --seq-qual-dir ../alignment/seq_qual \
  --output-dir ./region_plots
```

#### Key Options

- **Region**: `--chromosome CHR` (full chromosome) and/or `--region CHR:START-END`
- **Inputs** (at least one required): `--diversity-dir`, `--fst-dir`, `--pbe-dir`, `--seq-qual-dir` (TSV), and/or `--hdf5-dir` (collate output with `diversity_*.h5`, `fst_*.h5`, `pbe_*.h5`)
- `--window-size N`: Use only this window size (optional)
- `--reference-genome FILE`: For chromosome lengths
- `--y-value VALUE`: `value`, `rank`, or `quantile` for diversity/FST/PBE (default: value)
- `--width N`, `--height N`: Figure size in inches (default: 12, 10)
- `--dpi N`: DPI for PNG (default: 300)
- `--plot-format FORMAT`: png, pdf, svg, both, or all

#### Output

- One figure per run: `region_plot_CHR_START_END.png` (or `.pdf`/`.svg`) in `--output-dir`
- Panel order: Coverage (if seq_qual or diversity HDF5 with mean_coverage), π, FST, PBE. Only panels with data are included.
- Panel labels (a), (b), (c) and shared legend

#### Dependencies

- R packages: `ggplot2`, `dplyr`, `tidyr`, `readr`, `optparse`, **`patchwork`** (for combining panels). Optional: `hdf5r` when using `--hdf5-dir`.

---

### 4. plot_fst.R / plot_fst.sh

Creates ggplot2 plots of FST statistics with chromosome stripes, paneling options, and summary statistics for each sample pair.

#### Usage

**Using the wrapper script (recommended):**
```bash
./plot_fst.sh \
  --input-dir ../stats/fst_output \
  --output-dir ./fst_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --panel-by both \
  --plot-format both
```

**Direct R script usage:**
```bash
Rscript plot_fst.R \
  --input-dir ../stats/fst_output \
  --output-dir ./fst_plots \
  --panel-by window \
  --sample-pairs "Echo:Myv,Echo:Cheney"
```

#### Key Options

- `--input-dir DIR`: Directory containing FST TSV or HDF5 files (from `calculate_fst.sh` or collate) (required)
- `--input-format FORMAT`: `tsv`, `hdf5`, or `auto` (default: auto)
- `--y-value VALUE`: Y-axis: `value`, `rank`, or `quantile` (default: value)
- `--plot-style STYLE`: `line` or `line_points` (default: line)
- `--overlay-pairs`: Plot all pairs in one panel per window (color = pair)
- `--output-dir DIR`: Output directory for plots (default: ./)
- `--reference-genome FILE`: Reference genome FASTA file (optional, for chromosome lengths)
- `--panel-by OPTION`: Paneling option: `window`, `pair`, `both`, or `none` (default: `both`)
- `--plot-format FORMAT`: Plot format: `png`, `pdf`, or `both` (default: `png`)
- `--width N`, `--height N`: Plot dimensions in inches (default: 12, 8)
- `--dpi N`: DPI for PNG (default: 300)
- `--file-prefix PREFIX`: Prefix for output files (default: no prefix)
- `--sample-pairs LIST`: Comma-separated list of sample pairs to plot (default: all pairs)
- `--dry-run`: Preview commands without executing (wrapper only)

#### Output Files

- `fst_{sample_pair}.png` (or `.pdf`): Plots for each sample pair (e.g., `fst_Echo_Myv.png`)
- `fst_statistics.csv`: Summary statistics with columns:
  - `sample_pair`: Sample pair name (format: `sample1:sample2`)
  - `window_size`: Window size in bp (or NA if not applicable)
  - `chromosome`: Chromosome name or "genome" for genome-wide stats
  - `median`: Median FST value
  - `mean`: Mean FST value
  - `q95`: 95th percentile FST value

#### Features

- **Chromosome Stripes**: Alternating white/grey stripes for each chromosome
- **Paneling**: Options to panel by window size, sample pair, or both
- **Reference Lines**: Horizontal dashed lines showing:
  - Blue: Median (genome-wide)
  - Red: Mean (genome-wide)
  - Orange: 95th percentile (genome-wide)
- **Statistics**: Calculated both genome-wide and per-chromosome for each sample pair

#### Example

```bash
# Plot all sample pairs with paneling by window size and pair
./plot_fst.sh \
  --input-dir ../stats/fst_output \
  --output-dir ./fst_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --panel-by both \
  --plot-format both

# Plot specific sample pairs, paneled by window size
./plot_fst.sh \
  --input-dir ../stats/fst_output \
  --output-dir ./fst_plots \
  --panel-by window \
  --sample-pairs "Echo:Myv,Echo:Cheney,Myv:Cheney" \
  --plot-format png
```

---

### 6. plot_hist.R / plot_hist.sh

Genome-wide (or `--chromosome` / `--region` subset) distribution of one statistic from the same inputs as `plot_region` (TSV dirs or `--hdf5-dir`). Panels: `--panels histogram,density,ecdf`; `--overlay` colours groups in one panel; `--compare A,B` adds a bin-proportion difference panel. Run `./plot_hist.sh --help` for all options.

```bash
./plot_hist.sh --hdf5-dir ../stats/collated --stat fst --panels histogram,ecdf --overlay -o ./hist_plots
```

---

### 7. plot_outlier_csq.R / plot_outlier_csq.sh

Which kinds of variants sit in the outlier regions of each statistic, compared with the rest of the genome. Reads `outlier_csq_sites_*.csv` and `outlier_csq_background_*.csv` written by [`outlier_annotation.sh --variant-csq`](../outliers/README_identify_outliers_workflow.md#outlier-annotation-outlier_annotationsh).

- A region's `outlier_stat` can list several statistics (`fst:pbe:pi` when regions are merged across statistics); its sites count toward each. Each site is counted once per statistic.
- **Composition** (`outlier_csq_composition_<cat>.png`): stacked share of sites per consequence class for each statistic, next to *all sites* and *non-outlier sites*, with SNPs and indels in separate panels (`--variant-type split`; or `all`, `snp`, `indel`).
- **Enrichment** (`outlier_csq_enrichment_<cat>.png` / `.csv`): per statistic and class, log2 odds ratio of the class among that statistic's outlier sites vs non-outlier sites (Haldane-corrected, 95% CI), Fisher exact p and BH q (per variant type). Filled points have q < `--alpha`.
- `--category coding_effect` (default: nonsense, frameshift, missense, inframe_indel, other_coding, synonymous, not_coding, none) or `region_category` (cds, splice, UTRs, intron, ...). `none` = no csq record (intergenic); `--drop-none` removes it.

```bash
./plot_outlier_csq.sh --annotation-dir ../outliers/outlier_annotation -o csq_plots
./plot_outlier_csq.sh --annotation-dir ../outliers/outlier_annotation --category region_category --variant-type all -o csq_plots
```

---

### 8. plot_fst_csq.R / plot_fst_csq.sh

Per-site FST (SNPs) by consequence class, one panel per pair, from collate output: `single_position.h5` (`--single-position-merged`) or `variants.h5` + `fst_single.h5`. Collate needs `--fst-dir` with `fst_single.tsv` and `--variant-tsv-dir` / `--variant-csq` so sites carry `coding_effect` / `region_category`.

- **Violin + box** (`fst_csq_<cat>.png`): FST (or `--y-value quantile`) per class, site counts under each. With `--regions outlier_regions_*.csv`, each class is split into sites inside vs outside outlier regions (regions whose `outlier_stat` includes `--region-stat`, default `fst`; `any` for all).
- **Table** (`fst_csq_<cat>.csv`): n, mean, median, q95 per pair x class (x inside/outside); Wilcoxon p and BH q vs `--reference-category` (default `synonymous`) and, with `--regions`, inside vs outside.
- **Genome view** (`--manhattan`, `fst_csq_manhattan_<cat>.png`): per-site FST along the genome; sites with a coding consequence coloured by class, the rest grey; outlier regions shaded.

```bash
./plot_fst_csq.sh --hdf5-dir ../stats/collated \
  --regions ../outliers/out/outlier_regions_w1000_..._merge_only.csv --manhattan -o fst_csq
```

---

## Input File Formats

### Diversity TSV Files (for plot_diversity.R)

Expected format from `calculate_pi_theta.sh`:
- Tab-separated values
- Columns: `chromosome` (or `chr`), `position` (or `pos`, `start`), `pi`, `theta`, `tajima_d` (or `tajima-d`)
- Files named: `diversity.tsv` or `diversity_w{WINDOW}_s{STEP}.tsv`
- Sample names can be inferred from directory structure or filename

### FST TSV Files (for plot_fst.R)

Expected format from `calculate_fst.sh`:
- Tab-separated values
- Columns: `chromosome` (or `chr`), `position` (or `pos`, `start`), and FST columns
- FST columns named: `sample1:sample2.fst` (colon-separated) or `sample1_sample2_fst` (underscore-separated)
- Files named: `fst.tsv` or `fst_w{WINDOW}_s{STEP}.tsv`

---

## Dependencies

### Required Software

- **R** (version 3.6+): For R plotting scripts
- **grenedalf**: For cathedral plot generation
- **Bash** (version 4+): For wrapper scripts

### R Packages

The R scripts require the following packages (install with `install.packages()` if needed):

- `ggplot2`: For plotting
- `dplyr`: For data manipulation
- `readr`: For reading TSV files
- `optparse`: For command-line argument parsing
- `tidyr`: For data reshaping (plot_fst.R, plot_region.R, plot_outlier_csq.R, plot_fst_csq.R)
- `tibble`: For HDF5 reading (plot_common.R)
- `patchwork`: For combining panels in plot_region.R and plot_hist.R

Install all required packages:
```r
install.packages(c("ggplot2", "dplyr", "readr", "optparse", "tidyr", "tibble", "patchwork"))
```

When using HDF5 input (e.g. collate output or `--input-format hdf5` / `--hdf5-dir`), install `hdf5r`:
```r
install.packages("hdf5r")
```

### Optional

- **Reference genome FASTA**: For accurate chromosome length determination
  - If provided, the scripts will look for `.fai` (FASTA index) or `.dict` (SAM dictionary) files
  - If not provided, chromosome lengths are inferred from the data

---

## Workflow Examples

### Complete Analysis Pipeline

1. **Process samples** (generate BAM files):
```bash
./process_poolseq.sh \
  --sample-info sample-info.csv \
  --reference ../worm_q10.medaka.purged.fa \
  --output-dir ./interm \
  --parallel
```

2. **Calculate diversity statistics**:
```bash
./stats/calculate_pi_theta.sh \
  --sample-info sample-info.csv \
  --output-dir ./stats/diversity \
  --window-size 1000 \
  --window-size 10000 \
  --window-size 100000 \
  --reference-genome ../worm_q10.medaka.purged.fa
```

3. **Calculate FST**:
```bash
./stats/calculate_fst.sh \
  --sample-info sample-info.csv \
  --output-dir ./stats/fst \
  --window-type interval \
  --window-size 1000 \
  --window-size 10000 \
  --window-size 100000 \
  --reference-genome ../worm_q10.medaka.purged.fa
```

4. **Generate plots**:
```bash
# Diversity plots
./plotting/plot_diversity.sh \
  --input-dir ./stats/diversity \
  --output-dir ./plotting/diversity_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --panel-by both \
  --plot-format both

# FST plots
./plotting/plot_fst.sh \
  --input-dir ./stats/fst \
  --output-dir ./plotting/fst_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --panel-by both \
  --plot-format both

# FST cathedral plots
./plotting/plot_fst_cathedral.sh \
  --sample-info sample-info.csv \
  --output-dir ./plotting/cathedral_plots \
  --reference-genome ../worm_q10.medaka.purged.fa \
  --cathedral-width 2000 \
  --cathedral-height 600
```

---

## Troubleshooting

### Common Issues

1. **No input files found**
   - Check that input directory path is correct
   - Verify that TSV files match expected naming patterns (`*diversity*.tsv` or `*fst*.tsv`)
   - Ensure files are readable

2. **R package errors**
   - Install missing packages: `install.packages("package_name")`
   - Check R version: `R --version` (should be 3.6+)

3. **Chromosome stripe issues**
   - Provide `--reference-genome` for accurate chromosome lengths
   - Or ensure data contains chromosome and position columns

4. **FST column parsing errors**
   - Check that FST columns end with `.fst`
   - Verify column naming format (colon-separated: `sample1:sample2.fst` or underscore-separated: `sample1_sample2_fst`)

5. **Missing statistics in output**
   - Check that input files contain the expected columns
   - Verify that data has valid (non-NA) values for the statistics

---

## Figure standards

All plotting scripts (plot_diversity, plot_fst, plot_pbe, plot_region, plot_hist, plot_outlier_csq, plot_fst_csq) use:

- **Theme**: `theme_bw(base_size = 11, base_family = "sans")` for consistent typography
- **DPI**: PNG output uses 300 DPI by default (override with `--dpi N` in wrappers)
- **Panel labels**: `plot_region` adds (a), (b), (c) to stacked panels
- **Color**: A colorblind-friendly qualitative palette is used for sample/pair/trio in overlay and plot_region; consequence classes have fixed colours and severity order (`CSQ_EFFECT_LEVELS`, `CSQ_EFFECT_COLORS` in `plot_common.R`)
- **Dimensions**: `--width` and `--height` are in inches; vector (PDF/SVG) preferred for submission

---

## Notes

- All scripts support `--dry-run` mode for previewing commands (wrapper scripts)
- Output file names include prefixes if `--file-prefix` is specified
- Chromosome boundaries can be determined from reference genome (`.fai` or `.dict` files) or inferred from data
- Statistics are calculated both genome-wide and per-chromosome
- Reference lines (median, mean, 95th percentile) are shown for genome-wide statistics only
- Paneling options allow flexible visualization of multi-scale or multi-sample data

---

## Author

Based on workflow by JW
