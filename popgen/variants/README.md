# Variant Calling and Annotation for Poolseq Data

This directory contains scripts for finding variant sites in the processed pool-seq BAMs, annotating their functional consequences, and summarizing those annotations per window so they can be joined with the population genetics statistics in [popgen/stats/](../stats/).

1. **Site discovery**:
   - **`snp_call.sh`** (SNPs): BAM(s) → SNP sites from **grenedalf** read counts, using the same filters as `calculate_fst.sh`, so the sites are the SNPs behind FST/diversity. Writes sync, TSV, and BCF (per-pool allele depths, no genotypes).
   - **`indel_call.sh`** (indels): BAM(s) → indel sites from **freebayes** `--pooled-continuous`, filtered on read counts with the same thresholds as `snp_call.sh`. Writes TSV and BCF in the same layout.

   Neither uses genotype calling, which misses low-frequency alleles in pools (see [Why not bcftools call](#why-not-bcftools-call)).
2. **`variant_csq.sh` / `variant_csq.R`**: **bcftools csq** consequence annotation of any of these BCFs → site-level TSV.
3. **`variant_csq_windows.sh` / `variant_csq_windows.R`**: per-window counts of annotated sites, on the same window tiling as grenedalf FST/diversity output.

## Requirements

- **grenedalf** - Required for `snp_call.sh` (`conda install bioconda::grenedalf`)
- **freebayes** (≥ 1.3) - Required for `indel_call.sh` (`conda install bioconda::freebayes`)
- **bcftools** - Required for all scripts here (csq for `variant_csq.sh`)
- **samtools** - Required for `snp_call.sh` and `indel_call.sh` (reference indexing, BAM headers)
- **R** with packages **dplyr**, **readr**, **optparse**, **tidyr**, **tibble** - Required for `variant_csq.R` and `variant_csq_windows.R`

## Reference assemblies and repeat masking

These scripts follow the repository-wide repeat-masking policy described in [popgen/stats/README.md](../stats/README.md#repeat-masking-and-reference-assemblies):

- **Variants** (`snp_call.sh` / `indel_call.sh`): use the **same unmasked** reference as alignment; drop masked sites with the combined `mask.bed` from [`build_mask.sh`](../alignment/README.md#combined-mask-build_masksh) (repeats plus poorly mapping windows; a plain RepeatMasker BED from [`repeatmasker_bed.sh`](../alignment/repeatmasker_bed.sh) also works), via **`--filter-mask-total-bed`**. Use the same file as diversity/FST so sites and window statistics exclude the same positions.
- **Functional annotation** (`variant_csq.sh`): **`--reference` must be the unmasked assembly that matches your GFF** (not hard-masked). Run csq on a BCF that still has the same `chr`/`pos` as that assembly.

### Example

```bash
# Pooled SNP sites (unmasked ref; same mask and filters as diversity/FST)
./snp_call.sh --sample-info ../../sample_info.csv --reference ref.unmasked.fa \
  --filter-mask-total-bed mask/mask.bed --output-dir snps_out --threads 8

# Pooled indel sites (same filters and mask)
./indel_call.sh --sample-info ../../sample_info.csv --reference ref.unmasked.fa \
  --filter-mask-total-bed mask/mask.bed --output-dir indels_out --threads 8

# CSQ (unmasked ref + GFF; BCFs from the steps above)
./variant_csq.sh --bcf snps_out/snp_sites.bcf --reference ref.unmasked.fa \
  --gff annotation.gff3 --output-dir csq_snps
./variant_csq.sh --bcf indels_out/indel_sites.bcf --reference ref.unmasked.fa \
  --gff annotation.gff3 --output-dir csq_indels

# Windowed variant features from both tables, on the same grid as FST/diversity (for collate)
./variant_csq_windows.sh --sites-tsv csq_snps/variants_csq.tsv \
  --sites-tsv csq_indels/variants_csq.tsv --reference ref.unmasked.fa \
  -w 1000 --step-size 500 --output-dir variant_features
```

### Into collate and outliers

- **Single-site HDF5**: [`collate.sh`](../stats/README.md#collatesh--collater) `--variant-tsv-dir` (a directory holding `snp_sites.tsv` and `indel_sites.tsv`) with `--variant-csq csq_snps --variant-csq csq_indels` writes `variants.h5`, or `single_position.h5` with `--single-position-merged`. That file has one row per site, with `variant_type`, the per-pool `_AD` / `_alt_freq` columns, and the csq columns, including `genes` and `genes_protein_altering`. Single-position FST and PBE are added to SNP rows only.
- **Window HDF5**: `--variant-feature-dir variant_features` adds the `n_sites*` and `n_csq_*` / `n_region_*` / `n_effect_*` counts and the `genes` / `genes_protein_altering` lists to `diversity_w*.h5`, `fst_w*.h5` and `pbe_w*.h5`.
- **Outliers**: `identify_outliers.R` ranks windows only on the diversity, FST and PBE quantiles, so variant counts and gene lists are not carried into `outlier_windows_*.csv` / `outlier_regions_*.csv`. [`outlier_annotation.sh`](../outliers/README_identify_outliers_workflow.md#outlier-annotation-outlier_annotationsh) `--variant-csq csq_snps --variant-csq csq_indels --gff annotation.gff3` adds the sites, consequences and genes inside each outlier region. [`gene_collation.sh`](../outliers/README_identify_outliers_workflow.md#gene-collation-gene_collationsh) then lists those results by GFF gene.
- **Plots**: [`plot_outlier_csq.sh`](../plotting/README.md#7-plot_outlier_csqr--plot_outlier_csqsh) shows the consequence composition and enrichment of each outlier statistic's sites. [`plot_fst_csq.sh`](../plotting/README.md#8-plot_fst_csqr--plot_fst_csqsh) shows per-site FST by consequence class from the collated HDF5.

**Traceability:** `region_category` / `coding_effect` in csq output are plotting summaries; **`bcsq`** and **`csq_consequences`** retain full bcftools csq detail.

## Scripts

### snp_call.sh

Finds pooled SNP sites from **grenedalf** read counts, without genotype calling. grenedalf `sync` reads the BAMs (or an existing sync) and applies the same filters and defaults as `calculate_fst.sh`, so the output sites are the SNPs that FST and diversity are computed on. grenedalf writes every covered position, with filtered ones marked missing, so its output is streamed through a FIFO and only passing positions are written to disk.

- **Input** (one of): `--sample-info` (same CSV as `calculate_fst.sh`; `sample_name` becomes the column/sample name), `--bam` (repeatable; names from BAM basenames), or `--sync` (existing multi-sample sync; filters are applied again). Plus `--reference` (unmasked FASTA, `.fai` indexed) and `--output-dir`.
- **Filters** (defaults as `calculate_fst.sh`): `--min-coverage` 10 and `--max-coverage` 1000 (read depth per sample), `--min-count` 2 (per-sample count for a base to count as an allele), `--min-frequency` 0.01 (total minor allele frequency; or `--min-count-filter N` for a total alt count instead), `--min-total-read-depth` 10, `--biallelic-only`, `--filter-mask-total-bed` / `--filter-mask-total-fasta`. `--min-mapq` / `--min-baseq` default to 0, as in `calculate_fst.sh`; filter the BAMs at alignment (`--bam-min-mapq 20`) or raise these here.
- **Output** (`--prefix`, default `snp_sites`; keep `sites`, `calls` or `variant` in a custom prefix so `collate.sh --variant-tsv-dir` finds the TSV):
  - `snp_sites.sync`: passing positions only (grenedalf sync; a pool failing the per-sample filters at a passing site is `.:.:.:.:.:.`).
  - `snp_sites.tsv`: `chr`, `pos`, `ref`, `alt`, `dp`, `alt_freq` (total), then `<sample>_AD` (`ref,alt`) and `<sample>_alt_freq` per pool. Usable as a `collate.sh --variant-tsv-dir` input.
  - `snp_sites.bcf`: same sites with INFO/DP, INFO/AF and per-pool FORMAT/AD, FORMAT/DP; no GT. Input for `variant_csq.sh`, which detects the missing genotypes and annotates each site on its own.

**What counts as a site:** a position that is polymorphic across the pools after filtering. Positions where every pool is fixed for the same non-reference base are not sites (they are not variable within your samples, and FST excludes them too). ALT is the most frequent non-reference base in total; a third allele, if present, is dropped from the TSV/BCF but kept in the sync. Indels are not reported; use [`indel_call.sh`](#indel_callsh).

**Example:**
```bash
./snp_call.sh --sample-info ../../sample_info.csv --reference ref.unmasked.fa \
  --filter-mask-total-bed genome.out.bed --output-dir snps_out --threads 8
```

### indel_call.sh

Finds pooled indel sites with **freebayes** `--pooled-continuous`, which estimates allele frequencies from read counts instead of fitting a genotype at a fixed ploidy. SNP and MNP observations are discarded (`-I -X`), so freebayes calls only indels, and complex alleles are decomposed by `bcftools norm`. freebayes QUAL is close to 0 for real indels at 2–10% in pools, so sites are filtered on per-pool read counts with the same thresholds and defaults as `snp_call.sh`, not on QUAL.

- **Input** (one of): `--sample-info` (same CSV as `snp_call.sh`; `sample_name` becomes the column/sample name) or `--bam` (repeatable; names from the BAM `@RG SM`). Plus `--reference` (unmasked FASTA, `.fai` indexed) and `--output-dir`. Each BAM must have exactly one `@RG SM`, distinct across BAMs (`process_poolseq.sh` sets this), because freebayes assigns reads to pools by SM. BAMs must be coordinate-sorted and indexed for the per-chunk `--region` calls; a missing `.bai`/`.csi` is created with `samtools index` (next to the BAM, so its directory must be writable). `snp_call.sh` streams whole BAMs through grenedalf and does not need an index.
- **Filters** (defaults as `snp_call.sh`): `--min-coverage` 10 and `--max-coverage` 1000 (a pool outside this depth is set missing at that site), `--min-count` 2 (reads per pool for an allele to count; also freebayes `--min-alternate-count`), `--min-frequency` 0.01 (total minor allele frequency; also freebayes `--min-alternate-fraction`; or `--min-count-filter N`), `--min-total-read-depth` 10, `--filter-mask-total-bed`. freebayes also skips positions deeper than `--max-coverage` × number of pools (collapsed repeats). `--min-mapq` defaults to 20 and `--min-baseq` to 0.
- **Parallelism**: `--threads N` runs N freebayes processes on `--chunk-size` (default 1 Mb) genome chunks. Indels spanning a chunk boundary are called once (duplicates are removed after normalization).
- **Output** (`--prefix`, default `indel_sites`; keep `sites`, `calls` or `variant` in a custom prefix for collate):
  - `indel_sites.tsv`: `chr`, `pos`, `ref`, `alt`, `dp`, `alt_freq` (total), `indel_len` (positive = insertion), then `<sample>_AD` (`ref,alt`) and `<sample>_alt_freq` per pool.
  - `indel_sites.bcf`: same sites with INFO/DP, INFO/AF and per-pool FORMAT/AD, FORMAT/DP; no GT. Input for `variant_csq.sh`.
  - `--keep-raw` also keeps the unfiltered freebayes VCF (`indel_sites.raw.vcf.gz`).

Multiallelic sites (e.g. several lengths of a microsatellite) are split into one row per ALT allele. `<sample>_AD` holds the reads carrying the reference allele and that row's ALT allele, and `alt_freq` is ALT reads divided by read depth (all reads at the site), so the frequencies of the alleles at one site do not depend on each other. If `snp_sites.tsv` and `indel_sites.tsv` are in the same `collate.sh --variant-tsv-dir`, collate counts both together.

**Simulation** (two 200-chromosome pools at ~300x, see [Why not bcftools call](#why-not-bcftools-call)): recovered every planted indel from 2% frequency up, with frequencies within a few percent of the truth. That includes all three alleles at a CA-repeat site, a 10% deletion present in only one pool, and the site straddling a chunk boundary with 4 threads. bcftools mpileup and call found 2 of the 9 sites. A SNP control site was correctly excluded.

**Example:**
```bash
./indel_call.sh --sample-info ../../sample_info.csv --reference ref.unmasked.fa \
  --filter-mask-total-bed genome.out.bed --output-dir indels_out --threads 8
```

### Running SNP and indel calling together

`snp_call.sh` and `indel_call.sh` are independent, so run them at the same time (e.g. in two tmux windows) and split the cores between them. `snp_call.sh --threads` sets grenedalf's worker threads; `indel_call.sh --threads` runs that many single-threaded freebayes processes on genome chunks and scales close to linearly. freebayes is much slower per base than grenedalf, so give `indel_call.sh` the larger share (e.g. 4 and 8 of 12 cores). `variant_csq.sh --threads` only adds bcftools compression threads; csq itself is single-threaded and fast.

### Why not bcftools call

Genotype callers (`bcftools call`, GATK HaplotypeCaller) assign each sample a genotype under a fixed ploidy. `bcftools call` supports only ploidy 1 or 2, so it cannot model a pool of *n* chromosomes, and pooled allele frequencies are lost:

- **Diploid model:** only 0%, 50% and 100% alt frequencies are available, so an allele at a few percent in every pool is explained as sequencing error. In a simulation with two pools at ~400x, it called SNPs at 25% frequency and above and missed every site at 1–10%; raising the prior (`bcftools call -P`, up to 0.99) did not change this. For indels (two 200-chromosome pools at ~300x) it called 2 of 9 planted sites, missing everything at 2–10%, including a 10% deletion present in only one pool.
- **Haploid model:** worse; mixed sites must be all-ref or all-alt, and it called only the 90% and 100% sites.
- **GATK HaplotypeCaller at ploidy = pool size:** drops alleles at multiallelic sites (the genotype count explodes), runs about 17× slower, and produced many false positives.

On the same simulations, `snp_call.sh` recovered every polymorphic SNP from 2% up (the 1% site fell just under the shared 0.01 minimum-frequency filter, as for FST), and `indel_call.sh` recovered every indel from 2% up.

### variant_csq.sh / variant_csq.R

Runs **bcftools csq** on a BCF/VCF using the reference FASTA and a **GFF3** suitable for csq (see [bcftools consequence calling](https://samtools.github.io/bcftools/howtos/csq-calling.html)). Produces an annotated BCF and a **site-level TSV** with:

- **`bcsq`**: raw `INFO/BCSQ` from csq (full detail).
- **`csq_consequences`**: semicolon-separated list of verbatim consequence tokens parsed from BCSQ.
- **`region_category`**, **`coding_effect`**: derived columns for plotting, mapped from csq tokens via a built-in table (extend or override with `--map-tsv` in `variant_csq.R`). For variants with multiple transcript annotations, the single-row categories use the **highest-severity** consequence for these two columns only; the raw `bcsq` string is unchanged.
- **`genes`**: semicolon-separated gene names from all BCSQ records at the site; **`genes_protein_altering`**: only the genes whose record changes the protein (missense, frameshift, stop gained/lost, start lost, splice, inframe indels).
- **`variant_type`**: `snp` or `indel`, from the REF/ALT lengths.
- Optional **`--output-long-tsv`**: one row per BCSQ sub-record (`variant_csq.sh` passes through to R). Pointer records (`@pos`, which refer to a compound change annotated at another site) are skipped, and multi-term consequences such as `missense&inframe_altering` are split before mapping.

**Reference:** use the **unmasked** FASTA that matches the GFF (see [Reference assemblies and repeat masking](#reference-assemblies-and-repeat-masking) above). The BCF may be repeat-filtered; coordinates must still match that assembly.

**Scaffolds without genes:** bcftools csq aborts on sites whose chromosome has no GFF features. `variant_csq.sh` runs csq only on chromosomes present in the GFF and transfers BCSQ back onto the full BCF, so sites on gene-free scaffolds are kept (empty `bcsq`) while the REF/fasta consistency check stays active (`--force` is not used). If no chromosome names are shared between the BCF and GFF, the script stops (naming mismatch; see `--csq-arg --unify-chr-names`).

**Input without genotypes:** csq normally builds haplotypes from `GT` and skips records without one. If the input has no FORMAT/GT (as from `snp_call.sh` and `indel_call.sh`), `variant_csq.sh` adds `-s - -l`: samples are ignored and each site is annotated on its own, since which pooled alleles share a chromosome is unknown. Two SNPs in one codon are then annotated separately rather than as a combined change.

**Example:**
```bash
./variant_csq.sh --bcf snps_out/snp_sites.bcf --reference ref.unmasked.fa --gff genes.gff3 --output-dir csq_out
# Annotated sites: csq_out/variants_csq.tsv ; CSQ BCF: csq_out/calls.csq.bcf
```

### variant_csq_windows.sh / variant_csq_windows.R

Aggregates annotated sites into **per-window count tables** on grenedalf's interval grid ([Window boundaries](../stats/README.md#window-boundaries)). Preferred: **`--reference`** (FASTA or `.fai`, the same one given to diversity/FST) with **`-w` / `--step-size`** (repeatable, paired; step defaults to window/2). The grid is built from chromosome lengths, so every window is written, including those without sites. Alternatively, **`--template-dir`** (or **`--template-tsv`**) takes existing windowed FST or diversity TSVs whose names include `w{W}_s{S}` (e.g. `fst_w1000_s500.tsv`): with `--reference` the grid is built for those scales; without it the template's own windows are reused, which inherits any gaps in that table. Sites are assigned to windows by interval search per chromosome, so memory stays flat on large scaffolds. Output files: **`variant_features_w{W}_s{S}.tsv`** in `--output-dir`, with **`n_sites`**, **`n_sites_snp`** / **`n_sites_indel`**, dynamic **`n_csq_*`** (verbatim consequence tokens), **`n_region_*`**, and **`n_effect_*`** columns, plus **`genes`** / **`genes_protein_altering`** (sorted, semicolon-separated union over the window's sites), suitable for joining in collate. `--sites-tsv` can be repeated to combine the SNP and indel csq tables.

**Example:**
```bash
./variant_csq_windows.sh --sites-tsv csq_snps/variants_csq.tsv --sites-tsv csq_indels/variants_csq.tsv \
  --reference ref.unmasked.fa -w 1000 --step-size 500 -w 10000 --step-size 5000 --output-dir variant_features
```

The output directory can be passed to `collate.sh --variant-feature-dir` in [popgen/stats/](../stats/README.md#collatesh--collater), which left-joins these columns onto the diversity, FST and PBE windows.
