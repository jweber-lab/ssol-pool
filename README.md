# ssol-pool

This repository deals with **analysis of pool-sequenced samples** listed in **sample_info.csv**, run on the **stickleback lab server** of **Jesse Weber's lab**. Pipelines are organized under `popgen/`. This README is the **single source of truth** for canonical paths (see below).

## Single source of truth for `sample_info.csv`

The sample manifest used across the repo has one canonical location:

- **Path**: `sample_info.csv` at the **repository root** (i.e. `sample_info.csv` relative to the project root).

Scripts and documentation should reference this path. Any change to it requires updating all call sites and this README.

The other key input for downstream analyses is:

- **PBE comparisons**: `popgen/stats/pbe_comparisons.csv` — defines population trios for Population Branch Excess (PBE) calculations.

See [SAMPLES.md](SAMPLES.md) for a brief description of the populations (Vancouver, Alaska, Norway, Iceland).

---

## Popgen pipelines (overview)

The pipeline flow is: **alignment** → **stats** (diversity, FST, PBE, collation) and **variants** (calling, consequence annotation) → **plotting**, with **outliers** (identification and annotation) as a side branch. Conventions: Bash wrapper scripts (`.sh`) call R scripts (`.R`) in the same directory; inputs are driven by `sample_info.csv` and, where relevant, `pbe_comparisons.csv`.

### Directory layout

| Directory | Description | Details |
|-----------|-------------|---------|
| [popgen/alignment/](popgen/alignment/) | Alignment and pool-seq processing | Adapter trimming (BBtools), merging, BWA mapping, deduplication, BAM merging, mpileup/sync. Main scripts: `process_poolseq.sh`, `mapq_threshold_counts.sh`, `seq_qual_metrics.sh`. Full details: [popgen/alignment/README.md](popgen/alignment/README.md). |
| [popgen/stats/](popgen/stats/) | Diversity, FST, PBE, collation | π, θ, Tajima's D (grenedalf/popoolation2), FST, PBS/PBE, and collation of TSV/CSV into HDF5. Scripts: `calculate_pi_theta.sh`, `calculate_fst.sh`, `calculate_pbe.sh`, `collate.sh` / `collate.R`. Full details: [popgen/stats/README.md](popgen/stats/README.md). |
| [popgen/variants/](popgen/variants/) | Variant calling and annotation | SNP calling from grenedalf read counts (`snp_call`, same filters as FST), indel calling with freebayes `--pooled-continuous` (`indel_call`), **bcftools csq** annotation (`variant_csq`), and windowed variant-feature counts (`variant_csq_windows`) for collate. Scripts: `snp_call.sh`, `indel_call.sh`, `variant_csq.sh` / `variant_csq.R`, `variant_csq_windows.sh` / `variant_csq_windows.R`. Full details: [popgen/variants/README.md](popgen/variants/README.md). |
| [popgen/plotting/](popgen/plotting/) | R/ggplot2 plots | Diversity, FST, and PBE visualizations. Scripts: `plot_fst_cathedral.sh`, `plot_diversity.sh` / `plot_diversity.R`, `plot_fst.sh` / `plot_fst.R`, `plot_pbe.sh` / `plot_pbe.R`. Full details: [popgen/plotting/README.md](popgen/plotting/README.md). |
| [popgen/outliers/](popgen/outliers/) | Outlier identification and annotation | Identify outlier windows and annotate with BLAST. Scripts: `identify_outliers.sh` / `identify_outliers.R`, `outlier_annotation.sh` / `outlier_annotation.R`; `blast_config.example.yml` for BLAST configuration. |

---

## Conventions

- **Bash + R**: Each `.sh` wrapper invokes an R script in the same directory (e.g. `calculate_pbe.sh` calls `calculate_pbe.R`).
- **Inputs**: Use `sample_info.csv` at the repository root and `popgen/stats/pbe_comparisons.csv` for PBE; do not move or rename these files without updating this README and all script and documentation references.
- **`help/`**: Gitignored directory for local tool documentation (e.g. grenedalf wiki) so coders and agents can reference it offline; not versioned.
