# identify_outliers workflow

Step-by-step description of the **identify_outliers.sh** and **identify_outliers.R** pipeline, and recent changes.

---

## Shell script (identify_outliers.sh)

1. **Parse arguments**  
   Options (--hdf5-dir, --output-dir, --high-quantile, --low-quantile, --statistics, --window-size, --chromosome, etc.) are parsed. At least one of --fst-dir, --diversity-dir, or --hdf5-dir is required; --output-dir is required; either (--high-quantile and --low-quantile) or all four seed/expand quantiles are required.

2. **Validate**  
   Quantiles in [0,1], directories exist, Rscript and the R script path are found. Optional checks on FST/diversity dirs for presence of files.

3. **Build R command**  
   An `Rscript identify_outliers.R ...` command is built with all passed-through flags (output-dir, statistics, quantiles, merge-distance, merge-across-samples, window-size, chromosome, verbose, etc.).

4. **Run R script**  
   The command is run; stdout/stderr are teed into a timestamped log file under `{output-dir}/log/`.

5. **Post-run**  
   First few lines of each `outlier_*.csv` in the output dir are printed (and appended to the log).

---

## R script (identify_outliers.R) — main steps

### 1. Parse options and validate

- Options (paths, quantiles, statistics list, merge settings, filters) are read.
- Either single-threshold mode (--high-quantile, --low-quantile) or **seed-then-expand** (--seed-high/low, --expand-high/low) is enforced.
- Statistics list (e.g. pi, theta, tajima_d, pbe, fst) is parsed.

### 2. Read input data

- **HDF5 (--hdf5-dir):**  
  Diversity (`diversity_w*_s*.h5`), FST (`fst_w*.h5`), PBE (`pbe_*.h5`) are read. Optional **--window-size** restricts to that window size. Columns include chr, start, end, sample (or pair/trio), stat values, *_quantile, and where present: mean_coverage, mean_mapping_quality, n_snps.
- **TSV fallback:**  
  If --hdf5-dir not set, --diversity-dir and/or --fst-dir are used to read TSV/CSV with the same logical structure.
- PBE and FST can get mean_coverage etc. from diversity (by chr, start, end) when available.

### 3. Chromosome scope

- **Chromosome lengths** from data or from --reference-genome (FASTA index or dict).
- **Filters:** --top-n-chromosomes, --min-chromosome-length, --chromosome (single chr), and optionally --region-start / --region-end on that chr.
- Diversity, FST, and PBE are restricted to the selected chromosomes and, if set, to windows overlapping [region_start, region_end].

### 4. Per-statistic outlier detection

For each requested statistic:

- **Diversity (pi, theta, tajima_d):**  
  For each sample and each window size, windows are flagged as high or low outliers by quantile (or by seed-then-expand). Seed-then-expand: seed windows at strict quantile (e.g. 0.999/0.001), then expand with softer quantile (e.g. 0.99/0.01), then merge by gap (merge_distance). Optional depth/mapping-quality/SNP filters apply to seeds and expansion.
- **FST:**  
  Same logic per sample_pair.
- **PBE:**  
  Same logic per trio.

Each result is standardized (same columns: chr, start, end, value, quantile, sample or sample_pair, window_size, etc.) and stored in **outlier_results** (and, in expand mode, expanded regions in **expanded_region_results**).

### 5. Combine outliers into one wide table

- **all_outliers** = bind_rows of all outlier_results (long form: one row per window × statistic × sample/pair).
- A **value_col** is defined per row (e.g. `Cheney.pi`, `Echo:Myv.fst`).
- **wide_long** = (chr, start, end, value_col, value) from all_outliers, distinct on (chr, start, end, value_col).
- **pivot_wider** turns value_col into columns → one row per (chr, start, end) with columns like Cheney.pi, Echo.pi, …, and their _quantile columns. **window_meta** (window_size, mean_coverage, etc.) is joined from the first row per window.
- Initially only the sample(s) that passed the threshold have non-NA in their stat/quantile columns.

### 6. Fill step (recent)

- **Goal:** For every outlier window, fill in the **stat value and quantile for every sample** from the unfiltered data, not only for the sample that passed the threshold.
- **Keys:** When the wide table has **window_size**, the fill uses keys **(chr, start, end, window_size)** so the same window definition is used as for outlier detection.
- **Source data for fill:**  
  Diversity, FST, and PBE are filtered to the **same window_size(s)** as in the wide table and deduped (one row per (chr, start, end, sample) or per pair/trio) so the fill does not pull in other window sizes or duplicate rows.
- **build_windows_by_value_col** now returns tibbles that include **window_size** (when present in the source), so the fill pipeline can join on it.
- **full_long** = long table of (chr, start, end, window_size when present, value_col, value, quantile) for outlier windows only (semi_join with keys).
- **pivot_wider** uses **id_cols = join_by** (chr, start, end, and window_size when present) so there is exactly one row per window in full_values and full_quants → **one-to-one** join, no many-to-many.
- The wide table is then updated by left_join with this filled table so all sample columns get their value and quantile from the unfiltered data.

### 7. Build regions (if merge or seed-expand)

- If seed-expand was used, **expanded_region_results** are combined and optionally merged across samples (merge_distance).
- Regions get **n_windows**, **region_mean_*** (simple mean over overlapping windows), **region_max_***, **region_max_*_quantile** via **add_region_stat_summaries**.
- Optional region-level filters (e.g. min-region-snps) are applied.
- **outlier_regions_*.csv** is written (and optionally a merge-only regions file when not in expand mode).

### 8. Write outputs

- **outlier_windows_*.csv:** Wide table (chr, start, end, window_size, …, outlier_stat, outlier_direction, then all samplename.stat and samplename.stat_quantile columns).
- **outlier_regions_*.csv:** One row per merged/expanded region with region summaries.
- Filename suffix encodes window size, chromosome(s), quantiles, merge/expand, and optional region.

---

## Recent changes (summary)

1. **Fill from unfiltered data**  
   Stat and quantile columns are filled for **all** samples from the raw diversity/FST/PBE data so you see values even for samples that did not pass the threshold for that window.

2. **window_size in fill path**  
   To avoid many-to-many joins and duplicate rows:
   - Fill uses **(chr, start, end, window_size)** as the key when the wide table has window_size.
   - Source data for the fill are restricted to the **same window_size(s)** as the outlier table and deduped.
   - **build_windows_by_value_col** now passes through **window_size** in its returned tibbles.
   - **pivot_wider** in the fill step uses **id_cols = join_by** so there is one row per (chr, start, end, window_size).

3. **Verbose diagnostics for fill**  
   With --verbose, the script prints [fill] lines: wide_outliers row count and unique keys, diversity_data row count and unique (chr, start, end, sample), window_size in use, full_long and full_values/full_quants row counts and unique keys, and a note when duplicate keys would cause a many-to-many join.

4. **Output naming**  
   - **identify_outliers:** Removed duplicate scaffold in region filename when using --chromosome (whole chr): no longer append _chr_<chr> when chr_suffix already encodes the chromosome.
   - **outlier_annotation** (separate script): Output files renamed from outlier_regions_genes_* to **outlier_genes_*** and **outlier_genes_summary_***.

5. **Summary CSV (outlier_annotation)**  
   - List columns use semicolon separators; scalar columns first, long list columns (gene_ids, products) last; unique products only; no CDS/gene/mRNA duplicate rows in the summary (one row per region × db, unique gene_id and product).

The in-script comments and this README are the main documentation for the workflow and these changes.
