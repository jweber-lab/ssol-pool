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

The pipeline flow is: **alignment** → **stats** (diversity, FST, PBE, variant calling, collation) → **plotting**, with **outliers** (identification and annotation) as a side branch. Conventions: Bash wrapper scripts (`.sh`) call R scripts (`.R`) in the same directory; inputs are driven by `sample_info.csv` and, where relevant, `pbe_comparisons.csv`.

### Directory layout

| Directory | Description | Details |
|-----------|-------------|---------|
| [popgen/alignment/](popgen/alignment/) | Alignment and pool-seq processing | Adapter trimming (BBtools), merging, BWA mapping, deduplication, BAM merging, mpileup/sync. Main scripts: `process_poolseq.sh`, `mapq_threshold_counts.sh`, `seq_qual_metrics.sh`. Full details: [popgen/alignment/README.md](popgen/alignment/README.md). |
| [popgen/stats/](popgen/stats/) | Diversity, FST, PBE, variant calling, collation | π, θ, Tajima's D (grenedalf/popoolation2), FST, PBS/PBE, bcftools variant calling, and collation of TSV/CSV into HDF5. Scripts: `calculate_pi_theta.sh`, `calculate_fst.sh`, `calculate_pbe.sh`, `variant_call.sh`, `collate.sh` / `collate.R`. Full details: [popgen/stats/README.md](popgen/stats/README.md). |
| [popgen/plotting/](popgen/plotting/) | R/ggplot2 plots | Diversity, FST, and PBE visualizations. Scripts: `plot_fst_cathedral.sh`, `plot_diversity.sh` / `plot_diversity.R`, `plot_fst.sh` / `plot_fst.R`, `plot_pbe.sh` / `plot_pbe.R`. Full details: [popgen/plotting/README.md](popgen/plotting/README.md). |
| [popgen/outliers/](popgen/outliers/) | Outlier identification and annotation | Identify outlier windows and annotate with BLAST. Scripts: `identify_outliers.sh` / `identify_outliers.R`, `outlier_annotation.sh` / `outlier_annotation.R`; `blast_config.example.yml` for BLAST configuration. |

---

## Conventions

- **Bash + R**: Each `.sh` wrapper invokes an R script in the same directory (e.g. `calculate_pbe.sh` calls `calculate_pbe.R`).
- **Inputs**: Use `sample_info.csv` at the repository root and `popgen/stats/pbe_comparisons.csv` for PBE; do not move or rename these files without updating this README and all script and documentation references.
- **`help/`**: Gitignored directory for local tool documentation (e.g. grenedalf wiki) so coders and agents can reference it offline; not versioned.
