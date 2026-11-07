# Collation

`collate.sh` joins the population-genetic tables from [pop-gen/stats/](../stats/) into HDF5. It reads diversity, FST, PBE, QIBE, and LRPI TSV/CSV, ranks each statistic, and can attach sequence-quality windows from [alignment/seq_qual_metrics.sh](../../alignment/seq_qual_metrics.sh) and variant tables from [pop-gen/variants/](../variants/).

`collate.sh` is the wrapper. `collate.R` does the work. Missing or empty FST, PBE, QIBE, and LRPI directories are skipped. At least one of the statistic directories is required.

Windowed inputs share the grenedalf interval grid documented in [Window boundaries](../stats/README.md#window-boundaries): start = 1 + k·step, and the last window is the first whose end reaches the chromosome end.

## Requirements

- **R** with **dplyr**, **readr**, **optparse**, **tidyr**, **hdf5r**, and **tibble**

## collate.sh / collate.R

Collate writes HDF5 with a `/windows` group for windowed statistics and a `/sites` group for single sites. It computes an integer rank (1..N) and a 0–1 quantile per statistic: within (sample, window size) for diversity, within the window for FST, within each pair file for LRPI, within each trio file for PBE, and within each quartet file for QIBE.

- **Input**: `--diversity-dir`, `--fst-dir`, `--pbe-dir`, `--qibe-dir`, `--lrpi-dir` (at least one); optional `--seq-qual-dir`, `--variant-tsv-dir`, `--variant-csq`, **`--variant-feature-dir`** (directory containing `variant_features_wW_sS.tsv` from `variant_csq_windows`). `--qibe-dir` is the `calculate_qibe.sh` output; collate reads `*_qibe.tsv` (not `*_qibs.tsv`). `--lrpi-dir` is the `calculate_lrpi.sh` output; collate reads `*_lrpi.tsv`.
- **Scale matching**: diversity / FST / PBE / QIBE / LRPI are grouped by exact `w{W}_s{S}`. Seq-qual and variant-feature tables prefer that exact filename; if it is missing, collate uses a finer same-window table whose step divides the analysis step and keeps rows with `(start - 1) %% S == 0` (and `end - start + 1 == W`). So `seq_qual_metrics_w1000_s250.tsv` can fill `diversity_w1000_s500`. A different window width, or a step that is not a multiple of the source step, is not used.
- **Output**: In `--output-dir` (default `./collated`): `diversity_w{W}_s{S}.h5`, `fst_w{W}_s{S}.h5`, `pbe_w{W}_s{S}.h5`, `qibe_w{W}_s{S}.h5`, `lrpi_w{W}_s{S}.h5` (group `/windows`) and, unless `--single-position-merged`, `fst_single.h5` / `pbe_single.h5` / `qibe_single.h5` / `lrpi_single.h5` (group `/sites`). Columns are chr, start, end, sample or pairs (pop1, pop2 for LRPI, plus `pi_pop1` and `pi_pop2`; pop1, pop2, pop3 for PBE; pop1–pop4 for QIBE), stat values, and `*_rank`, `*_quantile`. LRPI is ranked within each pair file.
- **Window consequences** (`--variant-feature-dir` only): matching-scale `variant_features_*.tsv` columns are **left-joined** on `(chr, start, end)` into **windowed** diversity, FST, PBE, QIBE, and LRPI, adding `n_sites`, `n_sites_snp`, `n_sites_indel`, `n_csq_*`, `n_region_*`, `n_effect_*`, and `genes` / `genes_protein_altering`. They are not added to `fst_single.h5`, `pbe_single.h5`, `qibe_single.h5`, `lrpi_single.h5`, or `single_position.h5`. Build the feature tables with the same windows as the stats (`variant_csq_windows.sh --template-dir`, and `--no-nan-windows` on by default in diversity and FST); a left join still runs if a window is missing, and that row's counts are NA. Collate's window check names the mismatch.
- **Site consequences** (`variants.h5`, group `/sites`): `--variant-tsv-dir` and `--variant-csq` do not touch window HDF5. `--variant-tsv-dir` is a directory with `snp_sites.tsv` and/or `indel_sites.tsv` ([`snp_call.sh`](../variants/README.md#snp_callsh) / [`indel_call.sh`](../variants/README.md#indel_callsh)). Rows are stacked, one per site and ALT allele, with `variant_type` (`snp` or `indel`), `ref`, `alt`, `dp`, `alt_freq`, `indel_len` (indels) and per-pool `<sample>_AD` (`ref,alt`) and `<sample>_alt_freq`. `--variant-csq` (repeatable in `collate.sh`; comma-separated in `collate.R`) is a site TSV from [`variant_csq.sh`](../variants/README.md#variant_csqsh--variant_csqr) or a directory of them. A file path is used as given. A directory contributes TSVs whose names contain `variant`, `sites`, `calls`, or `csq` (so `variants.csq.tsv` and `variants_csq.tsv` both match); `*_long.tsv` is ignored. Tables with a `bcsq` or `csq_consequences` column are annotation: `bcsq`, `csq_consequences`, `region_category`, `coding_effect`, `genes`, and `genes_protein_altering` are joined onto the sites by `chr`, `pos`, `ref`, `alt` (by `chr`, `pos` only if `ref`/`alt` are absent). The same discovery applies inside `--variant-tsv-dir`. With no site TSV, the csq table is the site table.
- **`--single-position-merged`**: window HDF5 files are unchanged (and still receive `--variant-feature-dir`). `fst_single.h5`, `pbe_single.h5`, `qibe_single.h5`, `lrpi_single.h5`, and `variants.h5` are replaced by one `single_position.h5` (`/sites`): the site table, including csq columns, with per-locus FST, PBE, QIBE, and LRPI as wide columns (`Pop1:Pop2.fst`, `Pop1:Pop2:Pop3.pbe`, `Pop1:Pop2:Pop3:Pop4.qibe`, `Pop1:Pop2.lrpi`, with `_rank` / `_quantile`) on SNP rows only. An indel at the same position as a SNP does not inherit the SNP's values. Per-locus FST positions that are not in the site table are dropped when a site table is present. With no site or csq table, the rows are the per-locus FST (or PBE, QIBE, or LRPI) positions and have no consequence columns.
- **Masking and windows**: diversity and FST windows carry **`n_masked`** (masked bp, from grenedalf `total.masked`). Every windowed input is checked against the shared grid; see [Window boundaries](../stats/README.md#window-boundaries).
- **Options**: `--drop-all-na` drops rows or sites where every statistic value (every pair's FST, or PBE, QIBE, or LRPI) is NA before ranking. `--no-summary` skips the companion `*_summary.tsv` files. `--output-dir` defaults to `./collated`.
- **Sample names** keep their case from the input (FST pairs such as `Echo:Kjer`, `<sample>_AD` columns); standard column names (`chr`, `pos`, `ref`, ...) are lowercased.

```bash
./collate.sh --diversity-dir pi_out --fst-dir fst_out --pbe-dir pbe_out \
  --qibe-dir qibe_out --lrpi-dir lrpi_out --seq-qual-dir seq_qual_mask \
  --variant-tsv-dir sites --variant-csq csq_snps --variant-csq csq_indels \
  --variant-feature-dir variant_features -o collated
```

The statistic directories are outputs of `calculate_pi_theta.sh`, `calculate_fst.sh`, `calculate_pbe.sh`, `calculate_qibe.sh`, and `calculate_lrpi.sh` in [pop-gen/stats/](../stats/). An end-to-end run, from alignment through those scripts and into collate, is in [pop-gen/stats/README.md](../stats/README.md#end-to-end-example).
