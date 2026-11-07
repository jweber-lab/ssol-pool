# ssol-pool

This repository deals with **analysis of pool-sequenced samples** listed in **sample_info.csv**, run on the **stickleback lab server** of **Jesse Weber's lab**.

## Single source of truth for `sample_info.csv`

The sample manifest used across the repo has one canonical location:

- **Path**: `sample_info.csv` at the **repository root** (i.e. `sample_info.csv` relative to the project root).

Scripts and documentation should reference this path. Any change to it requires updating all call sites and this README.

The other key input for downstream analyses is:

- **PBE comparisons**: `pop-gen/stats/pbe_comparisons.csv` — defines population trios for Population Branch Excess (PBE) calculations.

See [SAMPLES.md](SAMPLES.md) for a brief description of the populations (Vancouver, Alaska, Norway, Iceland).

---



## Popgen pipelines (overview)

The pipeline flow is: **alignment** → **stats** (diversity, FST, PBE, QIBE, LRPI) and **variants** (calling, consequence annotation) → **collation** into HDF5 → **plotting**, with **outliers** and **genes** as a side branch. Conventions: Bash wrapper scripts (`.sh`) call R scripts (`.R`) in the same directory; inputs are driven by `sample_info.csv` and, where relevant, `pbe_comparisons.csv`.

### Directory layout


| Directory                              | Description                           | Details                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| -------------------------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [alignment/](alignment/)               | Alignment and pool-seq processing     | Adapter trimming (BBtools), merging, BWA mapping, deduplication, read groups, MAPQ filtering (unfiltered BAM kept), per-window alignment quality and the combined repeat + alignment-quality mask. Main scripts: `process_poolseq.sh`, `mapq_threshold_counts.sh`, `seq_qual_metrics.sh`, `build_mask.sh`, `repeatmasker_bed.sh`. Full details: [alignment/README.md](alignment/README.md).                                                                 |
| [pop-gen/stats/](pop-gen/stats/)       | Diversity, FST, PBE, QIBE, LRPI       | π, θ, Tajima's D (grenedalf/popoolation2), FST, PBS/PBE, QIBS/QIBE, and ln(π_A/π_B). Scripts: `calculate_pi_theta.sh`, `calculate_fst.sh`, `calculate_pbe.sh`, `calculate_qibe.sh`, `calculate_lrpi.sh`. Full details: [pop-gen/stats/README.md](pop-gen/stats/README.md). |
| [pop-gen/collation/](pop-gen/collation/) | HDF5 collation                      | Joins diversity, FST, PBE, QIBE, and LRPI tables, with optional sequence-quality and variant tables, into windowed and single-site HDF5. Scripts: `collate.sh` / `collate.R`. Full details: [pop-gen/collation/README.md](pop-gen/collation/README.md). |
| [pop-gen/variants/](pop-gen/variants/) | Variant calling and annotation        | SNP calling from grenedalf read counts (`snp_call`, same filters as FST), indel calling with freebayes `--pooled-continuous` (`indel_call`), **bcftools csq** annotation (`variant_csq`), and windowed variant-feature counts (`variant_csq_windows`) for collate. Scripts: `snp_call.sh`, `indel_call.sh`, `variant_csq.sh` / `variant_csq.R`, `variant_csq_windows.sh` / `variant_csq_windows.R`. Full details: [pop-gen/variants/README.md](pop-gen/variants/README.md). |
| [pop-gen/plotting/](pop-gen/plotting/) | R/ggplot2 plots                       | Diversity, FST, and PBE visualizations. Scripts: `plot_fst_cathedral.sh`, `plot_diversity.sh` / `plot_diversity.R`, `plot_fst.sh` / `plot_fst.R`, `plot_pbe.sh` / `plot_pbe.R`. Full details: [pop-gen/plotting/README.md](pop-gen/plotting/README.md).                                                                                                                                                                                                                     |
| [pop-gen/outliers/](pop-gen/outliers/) | Outlier identification and annotation | Identify outlier windows from diversity, FST, PBE, QIBE, and LRPI, annotate those intervals with csq/GFF and optional BLAST, then collate by assembly gene. Scripts: `identify_outliers.sh` / `identify_outliers.R`, `outlier_annotation.sh` / `outlier_annotation.R`, `gene_collation.sh` / `gene_collation.R`. Full details: [pop-gen/outliers/README.md](pop-gen/outliers/README.md). |
| [pop-gen/genes/](pop-gen/genes/)       | Per-gene reports                      | Combine gene lists into one directory, score every recorded HDF5 statistic over each gene, and keep sequence and proteome blastx results for genes added in that run. Script: `gene_report.sh` / `gene_report.R`. Full details: [pop-gen/genes/README.md](pop-gen/genes/README.md). |


---



## Conventions

- **Bash + R**: Each `.sh` wrapper invokes an R script in the same directory (e.g. `calculate_pbe.sh` calls `calculate_pbe.R`).
- **Inputs**: Use `sample_info.csv` at the repository root and `pop-gen/stats/pbe_comparisons.csv` for PBE; do not move or rename these files without updating this README and all script and documentation references.
- `help/`: Gitignored directory for local tool documentation (e.g. grenedalf wiki) so coders and agents can reference it offline; not versioned.

