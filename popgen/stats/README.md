# Population Genetics Analyses for Poolseq Data

This directory contains scripts for calculating population genetics statistics from processed poolseq data. These analyses use **grenedalf** as the primary tool (with direct BAM file input) and **popoolation2** as an optional fallback (for sync file compatibility).

## Overview

The analysis scripts calculate:

1. **π (Pi)**: Nucleotide diversity - the average number of pairwise nucleotide differences per site
2. **θ (Theta)**: Watterson's theta - estimated from the number of segregating sites
3. **Tajima's D**: Test statistic comparing π and θ to detect selection and demographic events
4. **PBS**: Population Branch Statistic — population-specific branch length from three pairwise FSTs
5. **PBE**: Population Branch Excess — PBS minus the median-scaled expectation from the sister-pair branch (`T23`)

## Requirements

### Software Dependencies

**Primary Tool (Recommended):**
- **grenedalf** - Primary tool for population genetics analyses
  - Available via Bioconda: `conda install bioconda::grenedalf`
  - GitHub: https://github.com/lczech/grenedalf
  - Wiki: https://github.com/lczech/grenedalf/wiki
  - Advantages: 200x faster, direct BAM input, lower memory usage

**Fallback Tool (Optional):**
- **popoolation2** - Fallback for backward compatibility
  - Available at: https://sourceforge.net/projects/popoolation2/
  - Required scripts: `pi-sliding.pl`, `theta-sliding.pl`, `fst-sliding.pl`
- **Java** - Required for popoolation2
- **Perl** - Required for running popoolation2 scripts

**Common Requirements:**
- **samtools** - Required for BAM file processing (if using grenedalf)
- **R** with packages **dplyr**, **readr**, **optparse**, **tidyr**, **hdf5r**, **tibble** - Required for `calculate_pbe.R` and `collate.R` (requirements for the variant scripts are in [popgen/variants/README.md](../variants/README.md))

### Input Files

**Primary Method (grenedalf):**
- **BAM files** - Direct input from the main processing pipeline (`process_poolseq.sh`)
  - One BAM file per population for diversity and FST
  - PBS/PBE are computed from FST tables (`calculate_pbe.sh`), not from BAMs
  - BAM files should be sorted and indexed

**Fallback Method (popoolation2):**
- **Sync files** - Generated from the main processing pipeline (`process_poolseq.sh`)
  - Sync files contain allele frequency information in the format:
    ```
    chr    pos    ref    pop1    pop2    pop3    ...
    ```
  - Where each population column contains allele counts in the format: `A:T:G:C:N:del`

## Repeat masking and reference assemblies

After **RepeatModeler / RepeatMasker**, you may have soft-masked and hard-masked FASTA files. This pipeline follows a single coherent policy:

1. **Alignment** ([`process_poolseq.sh`](../alignment/process_poolseq.sh)): use the **unmasked** reference. Optionally apply **`--bam-min-mapq 20 --bam-primary-only`** so the BAMs passed to stats exclude low-MAPQ and non-primary alignments.
2. **Diversity (π, θ, Tajima's D)** (`calculate_pi_theta.sh`): BAM input; exclude masked positions with **`--filter-mask-total-bed`** (same BED as below). Use the combined `mask.bed` from [`build_mask.sh`](../alignment/README.md#combined-mask-build_masksh) (repeats plus poorly mapping windows), or a plain RepeatMasker BED.
3. **FST** (`calculate_fst.sh`): same **`--filter-mask-total-bed`** (or `--filter-mask-total-fasta`) as diversity. Per-locus (`--window-type single`) and windowed FST keep **biallelic SNPs** and drop all-NA windows/sites by default (`--filter-total-only-biallelic-snps`, `--no-nan-windows`). Turn those off with `--no-filter-total-only-biallelic-snps` / `--keep-nan-windows`. **`--no-extra-columns`** drops grenedalf's per-pair filter-count columns.
4. **Variants** ([`snp_call.sh`](../variants/README.md#snp_callsh) for SNPs, [`indel_call.sh`](../variants/README.md#indel_callsh) for indels): **same unmasked** reference as alignment; same **`--filter-mask-total-bed`** as diversity/FST.
5. **Functional annotation** ([`variant_csq.sh`](../variants/README.md#variant_csqsh--variant_csqr)): **`--reference` must be the unmasked assembly that matches your GFF** (not hard-masked). Run csq on a BCF that still has the same `chr`/`pos` as that assembly.

**PBE** (`calculate_pbe.sh`) uses FST outputs, so repeat handling follows step 3.

Masked positions count as missing in grenedalf, so they are removed from the valid-loci denominator instead of being treated as invariant. Collate keeps grenedalf's `total.masked` as **`n_masked`** (masked bp per window) in the diversity and FST HDF5, so windows that are mostly masked can be filtered downstream.

Grenedalf expects the BED mask in **reference (`.fai`) order**. `calculate_pi_theta.sh`, `calculate_fst.sh` and `snp_call.sh` rewrite `--filter-mask-total-bed` into that order and add a 1 bp interval on chromosomes that have none (grenedalf 0.6.3 otherwise segfaults). Tiny leftover scaffolds often have no repeat or seq-qual interval, so they used to be appended after the last masked chromosome and abort with `Invalid sorting order` (for example `scaffold_171` then `scaffold_170`).

### Window boundaries

Every windowed table must use the same `(chr, start, end)` so collate's exact joins line up. All scripts use grenedalf's interval grid: start = 1 + k·step, with the last window being the first whose end reaches the chromosome end (the end is not clipped).

- **Diversity and FST require `--reference-genome`** for interval windows (FST also requires it whenever a mask is given). Without it, grenedalf ends each sample's windows at its last covered position, so samples and pairs get different grids, and with a mask grenedalf 0.6.3 crashes partway through. With a reference, all samples and pairs share one grid. Chromosomes with no data at all are omitted. Fully masked or empty windows are dropped by default (`--no-nan-windows` in `calculate_pi_theta.sh` and `calculate_fst.sh`); `--keep-nan-windows` writes them as `nan`.
- **PBE** keeps FST's `start` / `end`. Older PBE files with only `pos` get `end = start + window - 1` in collate.
- **`seq_qual_metrics.sh`** builds the same grid from the `.fai` and fails if BAM lengths differ from it.
- **`variant_csq_windows.sh`** builds the same grid from `--reference` and `-w` / `--step-size`, and by default omits windows with no sites (`--no-nan-windows`). Pass **`--template-dir`** (diversity or FST output) so the rows are exactly the masked-stat windows. **`--filter-mask-total-bed`** is the same mask as diversity/FST; sites inside it are not counted.
- **Collate** checks each table against the grid and against the first table at the same scale (on shared chromosomes), and prints `Window check ... matches ...` or a `WARNING window check` naming the off-grid or unmatched windows. For `--seq-qual-dir` and `--variant-feature-dir`, an exact `w{W}_s{S}` table is used when present; otherwise collate keeps the analysis-grid windows from a finer same-window table when the analysis step is a multiple of the source step (e.g. `seq_qual_metrics_w1000_s250.tsv` supplies the 500 bp starts for diversity `w1000_s500`).

Use the same reference FASTA (or its `.fai`) and the same window width for every script. Analysis step sizes should match, or be an integer multiple of, the seq_qual / variant-feature step.

### End-to-end example

Assume `ref.unmasked.fa`, `mask/mask.bed` from [`build_mask.sh`](../alignment/README.md#combined-mask-build_masksh), and BAMs from alignment.

```bash
# 1–2. Alignment and mask (in popgen/alignment/)
./process_poolseq.sh --reference ref.unmasked.fa --bam-min-mapq 20 --bam-primary-only ...
./seq_qual_metrics.sh --sample-info ../../sample_info.csv --reference-genome ref.unmasked.fa \
  -w 1000 --step-size 250 -o seq_qual_mask
./build_mask.sh --reference ref.unmasked.fa --repeatmasker-out ref.fa.out \
  --seq-qual-dir seq_qual_mask -o mask

# 3. Diversity + FST (same mask and reference, so windows and masked bp match)
./calculate_pi_theta.sh --sample-info ../../sample_info.csv --reference-genome ref.unmasked.fa \
  --filter-mask-total-bed mask/mask.bed --window-size 1000 --step-size 500 -o pi_out ...
./calculate_fst.sh --sample-info ../../sample_info.csv --reference-genome ref.unmasked.fa \
  --window-type interval --window-size 1000 --step-size 500 \
  --filter-mask-total-bed mask/mask.bed -o fst_out ...
./calculate_pbe.sh --fst-dir fst_out --pop1-name Echo --pop2-name Kjer --pop3-name Cheney -o pbe_out

# 4. Pooled SNP and indel sites (unmasked ref; same repeat BED and filters as FST) - scripts in popgen/variants/
../variants/snp_call.sh --sample-info ../../sample_info.csv --reference ref.unmasked.fa \
  --filter-mask-total-bed mask/mask.bed --output-dir snps_out
../variants/indel_call.sh --sample-info ../../sample_info.csv --reference ref.unmasked.fa \
  --filter-mask-total-bed mask/mask.bed --output-dir indels_out

# 5. CSQ (unmasked ref + GFF; BCFs from step 4)
../variants/variant_csq.sh --bcf snps_out/snp_sites.bcf --reference ref.unmasked.fa \
  --gff annotation.gff3 --output-dir csq_snps
../variants/variant_csq.sh --bcf indels_out/indel_sites.bcf --reference ref.unmasked.fa \
  --gff annotation.gff3 --output-dir csq_indels

# 6. Windowed variant features on the same grid, then collate
../variants/variant_csq_windows.sh --sites-tsv csq_snps/variants_csq.tsv \
  --sites-tsv csq_indels/variants_csq.tsv --reference ref.unmasked.fa \
  -w 1000 --step-size 500 --template-dir pi_out \
  --filter-mask-total-bed mask/mask.bed --output-dir variant_features
mkdir -p sites && cp snps_out/snp_sites.tsv indels_out/indel_sites.tsv sites/
./collate.sh --diversity-dir pi_out --fst-dir fst_out --pbe-dir pbe_out --seq-qual-dir seq_qual_mask \
  --variant-tsv-dir sites --variant-csq csq_snps --variant-csq csq_indels \
  --variant-feature-dir variant_features -o collated
```

**Traceability:** `region_category` / `coding_effect` in csq output are plotting summaries; **`bcsq`** and **`csq_consequences`** retain full bcftools csq detail. See [popgen/alignment/README.md](../alignment/README.md) for BAM-filter details.

## Scripts

**Note:** Per-window coverage and mapping quality are computed by `seq_qual_metrics.sh` in [popgen/alignment/](popgen/alignment/). Its output TSV can be passed to collate via `--seq-qual-dir`. A mask-scale table (e.g. `-w 1000 --step-size 250`) is enough when the analysis window matches and the analysis step is a multiple of that step.

**Variant sites and annotation** (`snp_call.sh`, `indel_call.sh`, `variant_csq.sh`, `variant_csq_windows.sh`) live in [popgen/variants/](../variants/README.md); their site tables and consequence annotations are collated with `--variant-tsv-dir` / `--variant-csq`, and their per-window feature tables with `--variant-feature-dir`.

### collate.sh / collate.R

Collates diversity, FST, and PBE TSV/CSV into **HDF5** with `/windows` group(s). Computes integer rank (1..N) and 0–1 quantile per stat within (sample, window_size) for diversity, and within window for FST/PBE.

- **Input**: `--diversity-dir`, `--fst-dir`, `--pbe-dir` (at least one); optional `--seq-qual-dir`, `--variant-tsv-dir`, `--variant-csq`, **`--variant-feature-dir`** (directory containing `variant_features_wW_sS.tsv` from `variant_csq_windows`).
- **Scale matching**: diversity / FST / PBE are grouped by exact `w{W}_s{S}`. Seq-qual and variant-feature tables prefer that exact filename; if it is missing, collate uses a finer same-window table whose step divides the analysis step and keeps rows with `(start - 1) %% S == 0` (and `end - start + 1 == W`). So `seq_qual_metrics_w1000_s250.tsv` can fill `diversity_w1000_s500`. A different window width, or a step that is not a multiple of the source step, is not used.
- **Output**: In `--output-dir`: `diversity_w{W}_s{S}.h5`, `fst_w{W}_s{S}.h5`, `pbe_w{W}_s{S}.h5` (and `fst_single.h5` / `pbe_single.h5` for per-locus). Each HDF5 has group `/windows` (or `/sites`) with chr, start, end, sample/pairs (pop1, pop2, pop3 for PBE), stat values, and `*_rank`, `*_quantile`. With `--variant-feature-dir`, matching-scale `variant_features_*.tsv` columns are **left-joined** on `(chr, start, end)` so window HDF5 also carries `n_sites`, `n_sites_snp`, `n_sites_indel`, `n_csq_*`, `n_region_*`, `n_effect_*` counts and `genes` / `genes_protein_altering` (semicolon-separated gene names with a site in the window) alongside population-genetics stats.
- **Variant sites** (`variants.h5`, group `/sites`): `--variant-tsv-dir` points at a directory with `snp_sites.tsv` and/or `indel_sites.tsv` ([`snp_call.sh`](../variants/README.md#snp_callsh) / [`indel_call.sh`](../variants/README.md#indel_callsh)). Their rows are stacked, one per site and ALT allele, with `variant_type` (`snp` or `indel`), `ref`, `alt`, `dp`, `alt_freq`, `indel_len` (indels) and per-pool `<sample>_AD` (`ref,alt`) and `<sample>_alt_freq`. `--variant-csq` (repeatable in `collate.sh`; comma-separated in `collate.R`) takes `variants_csq.tsv` files or their directories from [`variant_csq.sh`](../variants/README.md#variant_csqsh--variant_csqr); their `bcsq`, `csq_consequences`, `region_category`, `coding_effect`, `genes` and `genes_protein_altering` columns are joined onto the sites by `chr`, `pos`, `ref`, `alt`. A `variants_csq.tsv` inside `--variant-tsv-dir` is treated the same way, and with no site TSVs the csq tables become the sites. `*_long.tsv` files are ignored.
- **`--single-position-merged`**: replaces `fst_single.h5`, `pbe_single.h5` and `variants.h5` with one `single_position.h5` (`/sites`): the variant table above with single-position FST and PBE added as wide columns (`Pop1:Pop2.fst`, `Pop1:Pop2:Pop3.pbe`, with `_rank` / `_quantile`). FST and PBE are joined only onto SNP rows, so an indel at the same position as a SNP does not inherit the SNP's values; positions with FST but no called site are not included. Without variant TSVs, the rows are the FST (or PBE) positions.
- **Masking and windows**: diversity and FST windows carry **`n_masked`** (masked bp, from grenedalf `total.masked`). Every windowed input is checked against the shared grid; see [Window boundaries](#window-boundaries).
- **Options**: `--drop-all-na` — drop rows/sites where every statistic value (e.g. every pair's FST, or PBE) is NA/NaN before ranking and collation.
- **Sample names** keep their case from the input (FST pairs such as `Echo:Kjer`, `<sample>_AD` columns); standard column names (`chr`, `pos`, `ref`, ...) are lowercased.
- **Run**: `./collate.sh --diversity-dir pi_out --fst-dir fst_out --pbe-dir pbe_out -o collated`

### calculate_pi_theta.sh

Calculates nucleotide diversity (π), Watterson's theta (θ), and Tajima's D from BAM files (grenedalf) or sync files (popoolation2 fallback).

**Repeat masking:** pass **`--filter-mask-total-bed`** (RepeatMasker BED) so grenedalf excludes repeat intervals from windowed diversity—the same option as `calculate_fst.sh`. BAMs should come from alignment to the unmasked reference (see [Repeat masking](#repeat-masking-and-reference-assemblies)).

#### Usage

**Primary method with grenedalf (recommended):**
```bash
./calculate_pi_theta.sh \
  --bam /path/to/sample.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir /path/to/output \
  --window-size 10000 \
  --pool-size 50
```

**Fallback method with popoolation2:**
```bash
./calculate_pi_theta.sh \
  --sync /path/to/populations.sync \
  --popoolation2 /path/to/popoolation2 \
  --output-dir /path/to/output \
  --window-size 1000 \
  --step-size 1000
```

#### Options

**Input (choose one):**
- `-b, --bam FILE`: BAM file(s) - can be specified multiple times for multiple populations (primary method)
- `-s, --sync FILE`: Sync file (fallback method)

**Tool selection:**
- `-g, --grenedalf PATH`: Path to grenedalf executable (default: check PATH)
- `-p, --popoolation2 DIR`: Directory containing popoolation2 scripts (for fallback)
- `--use-popoolation2`: Force use of popoolation2 instead of grenedalf

**Required:**
- `-o, --output-dir DIR`: Output directory for results

**Optional:**
- `-w, --window-size N`: Window size in bp (default: 1000)
- `--step-size N`: Step size in bp (default: 1000)
- `--min-coverage N`: Minimum coverage per site (default: 4)
- `--max-coverage N`: Maximum coverage per site (default: 200)
- `--min-count N`: Minimum allele count (default: 2, popoolation2 only)
- `--pool-size N`: Pool size (haploid chromosomes) for Tajima's D calculation (default: 50)
- `--reference-genome FILE`: Reference FASTA (required with grenedalf); fixes the window grid across samples (see [Window boundaries](#window-boundaries))
- `--filter-mask-total-bed FILE`: Exclude intervals in BED from windowed diversity (`mask.bed` from `build_mask.sh`, or a RepeatMasker `*.out.bed`)
- `--filter-mask-total-fasta FILE`: Alternative grenedalf mask FASTA (use one of BED or FASTA)
- `--no-nan-windows`: Omit windows where every value is NA (default: on). `--keep-nan-windows` writes those rows.
- `-t, --threads N`: Number of threads (default: 1)
- `--dry-run`: Preview commands without executing (dry-run mode)

#### Output Files

- `pi.txt`: Nucleotide diversity (π) per window
  - Columns: chr, pos, covered_sites, pi
- `theta.txt`: Watterson's theta (θ) per window
  - Columns: chr, pos, covered_sites, theta
- `tajima_d.txt`: Tajima's D per window
  - Columns: chr, pos, covered_sites, pi, theta, Tajima_D
- `pi_snps.txt`: SNPs used for π calculation
- `theta_snps.txt`: SNPs used for θ calculation

#### Examples

**Using grenedalf with BAM file (recommended):**
```bash
./calculate_pi_theta.sh \
  --bam ../output/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir ./pi_theta_results \
  --window-size 10000 \
  --step-size 5000 \
  --pool-size 50 \
  --threads 8
```

**Using popoolation2 with sync file (fallback):**
```bash
./calculate_pi_theta.sh \
  --sync ../output/Echo_Kjer.sync \
  --popoolation2 /opt/popoolation2 \
  --output-dir ./pi_theta_results \
  --window-size 10000 \
  --step-size 5000 \
  --min-coverage 10 \
  --max-coverage 100
```

### calculate_pbe.sh / calculate_pbe.R

Computes Population Branch Statistic (PBS) and Population Branch Excess (PBE) from **FST tables** produced by `calculate_fst.sh` (grenedalf). It does not read BAMs or sync files. Names in `--pop1-name` / `--pop2-name` / `--pop3-name` must match the FST pair columns (grenedalf’s `PopA:PopB.fst`).

**PBS = (T12 + T13 − T23) / 2**, with **Tij = −ln(1 − FSTij)** (FST ≥ 1 is capped). Pop1 is the focal branch; pop2 and pop3 are the sister pair.

**PBE = PBS − T23 × median(PBS) / median(T23)**. Medians are genome-wide (`--median-scope genome`, default) or per chromosome. This is an excess relative to the sister-pair branch, not a Z-score.

Windows are whatever the FST file used: `start` / `end` are copied when present so collate can join on the shared grid. Run FST with the same `--reference-genome`, `--window-size`, and `--step-size` as diversity. Per-locus `fst_single` is biallelic SNPs with all-NA sites omitted by default in `calculate_fst.sh` (add `--no-extra-columns` there if you also want to drop per-pair count columns).

#### Usage

```bash
# One trio, every *fst*.tsv/csv in the FST directory (single and each w{W}_s{S})
./calculate_pbe.sh --fst-dir fst_out \
  --pop1-name Echo --pop2-name Kjer --pop3-name Cheney -o pbe_out

# One FST file
./calculate_pbe.sh --fst-file fst_out/fst_w1000_s500.tsv \
  --pop1-name Echo --pop2-name Kjer --pop3-name Cheney -o pbe_out

# Several trios (CSV header: pop1,pop2,pop3). Duplicate pop1 + {pop2,pop3} pairs are skipped
# unless --compute-dup-comps. --sample-info warns if a name is missing from sample_name.
./calculate_pbe.sh --fst-dir fst_out --pbe-comparisons pbe_comparisons.csv \
  --sample-info ../../sample_info.csv -o pbe_out
```

#### Options

**Input (choose one):**
- `-f, --fst-file FILE`: one FST TSV/CSV from `calculate_fst.sh`
- `-d, --fst-dir DIR`: directory of `*fst*.tsv` / `*fst*.csv` (top level)

**Trio (choose one):**
- `--pop1-name` / `--pop2-name` / `--pop3-name`: single comparison
- `--pbe-comparisons FILE`: CSV with columns `pop1`, `pop2`, `pop3`

**Required:** `-o, --output-dir DIR`

**Optional:** `--sample-info FILE`, `--compute-dup-comps`, `--file-prefix PREFIX`, `--output-format tsv|csv|hdf5` (default tsv), `--median-scope genome|chromosome`, `--drop-all-na` (drop rows with any of the three FSTs missing; skip the file if none remain), `--dry-run`

#### Output files

For each FST file × trio (label is `w{W}_s{S}` or `single` from the FST name):

- `{label}_{pop1}_{pop2}_{pop3}_pbs.tsv` — chr, pos, start, end (if present), pop1, pop2, pop3, three FST columns, T12, T13, T23, PBS
- `{label}_{pop1}_{pop2}_{pop3}_pbe.tsv` — same plus PBE. **Collate reads these** (`--pbe-dir`)
- `--output-format hdf5`: `{label}_{pop1}_{pop2}_{pop3}_pbs_pbe.h5` (not used by collate)
- Logs: `{output_dir}/log/calculate_pbe_YYYYmmdd_HHMMSS.log`

The FST table must contain all three pairs. Population names keep their case.

## Workflow Example

### Complete Analysis Pipeline

1. **Process individual samples** (using main pipeline):
```bash
for sample in Echo_Pool_S1 Kjer_Pool_S4 Outgroup_Pool_S2; do
  ../process_poolseq.sh \
    --sample-name "$sample" \
    --read1 "/path/to/${sample}_R1.fastq.gz" \
    --read2 "/path/to/${sample}_R2.fastq.gz" \
    --reference /path/to/reference.fna.gz \
    --output-dir "/path/to/output/${sample}" \
    --adapters /path/to/adapters.fa \
    --bbtools-dir /path/to/bbtools \
    --grenedalf /usr/bin/grenedalf \
    --use-grenedalf-sync \
    --threads 16
done
```

2. **Calculate π and θ using grenedalf (recommended)**:
```bash
./calculate_pi_theta.sh \
  --bam /path/to/output/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir ./diversity_analysis \
  --window-size 10000 \
  --step-size 5000 \
  --pool-size 50 \
  --threads 8
```

3. **Calculate FST, then PBS/PBE** (same reference, mask, and window/step as diversity):
```bash
./calculate_fst.sh \
  --sample-info ../../sample_info.csv \
  --reference-genome /path/to/ref.unmasked.fa \
  --window-type interval --window-size 10000 --step-size 5000 \
  --filter-mask-total-bed mask/mask.bed \
  --output-dir ./fst_out \
  --threads 8
./calculate_pbe.sh \
  --fst-dir ./fst_out \
  --pop1-name Echo --pop2-name Kjer --pop3-name Outgroup \
  --output-dir ./pbe_out
```

### Alternative Workflow with Popoolation2 (Fallback)

If grenedalf is not available, you can use popoolation2 with sync files:

1. **Process samples and create sync files**:
```bash
for sample in Echo_Pool_S1 Kjer_Pool_S4 Outgroup_Pool_S2; do
  ../process_poolseq.sh \
    --sample-name "$sample" \
    --read1 "/path/to/${sample}_R1.fastq.gz" \
    --read2 "/path/to/${sample}_R2.fastq.gz" \
    --reference /path/to/reference.fna.gz \
    --output-dir "/path/to/output/${sample}" \
    --adapters /path/to/adapters.fa \
    --bbtools-dir /path/to/bbtools \
    --popoolation2 /opt/popoolation2/mpileup2sync.jar
done
```

2. **Create combined sync file**:
```bash
# Create combined mpileup
samtools mpileup -B \
  /path/to/output/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  /path/to/output/Kjer_Pool_S4/Kjer_Pool_S4_All_seq.dedup.bam \
  /path/to/output/Outgroup_Pool_S2/Outgroup_Pool_S2_All_seq.dedup.bam \
  -o combined.mpileup

# Convert to sync
java -jar /opt/popoolation2/mpileup2sync.jar \
  --input combined.mpileup \
  --output combined.sync \
  --fastq-type sanger \
  --min-qual 20 \
  --threads 8
```

3. **Calculate π and θ with popoolation2**:
```bash
./calculate_pi_theta.sh \
  --sync combined.sync \
  --popoolation2 /opt/popoolation2 \
  --output-dir ./diversity_analysis \
  --window-size 10000 \
  --step-size 5000
```

4. **PBS/PBE** in this repo always come from `calculate_fst.sh` tables (grenedalf), not from popoolation2. After a popoolation2 diversity run, still compute FST with `calculate_fst.sh` on the BAMs, then `calculate_pbe.sh --fst-dir fst_out ...`.

## Dry-Run Mode

Both analysis scripts support `--dry-run` mode for testing configurations before running full analyses. In dry-run mode, scripts will:
- Check input file existence
- Show all commands that would be executed
- Display output files that would be created
- Skip actual command execution

**Example:**
```bash
./calculate_pi_theta.sh \
  --bam /path/to/sample.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir /path/to/output \
  --window-size 10000 \
  --dry-run
```

This is particularly useful for:
- Verifying BAM file paths and formats
- Testing grenedalf/popoolation2 configurations
- Understanding what output files will be generated
- Debugging parameter settings

## Interpreting Results

### π (Nucleotide Diversity)

- **High π**: High genetic diversity, may indicate balancing selection or large effective population size
- **Low π**: Low genetic diversity, may indicate selective sweeps, bottlenecks, or small effective population size
- Compare π values between populations to identify regions with differential diversity

### θ (Watterson's Theta)

- Estimates diversity from the number of segregating sites
- Should be similar to π under neutral evolution
- Discrepancies between π and θ can indicate selection or demographic events

### Tajima's D

- **D = 0**: Consistent with neutral evolution and constant population size
- **D > 0**: Excess of intermediate-frequency alleles
  - May indicate balancing selection, population structure, or population contraction
- **D < 0**: Excess of rare alleles
  - May indicate selective sweeps, population expansion, or purifying selection
- **Significance**: |D| > 2 is often considered significant, but interpretation depends on demographic history

### PBS (Population Branch Statistic)

- **Positive PBS**: Longer branch in pop1 than expected from the pop2–pop3 split; population-specific divergence
- **High PBS**: Candidate regions for pop1-specific selection
- **Negative PBS**: Shared variation or gene flow involving pop1
- PBS scales with background FST, so raw PBS is harder to compare across the genome than PBE

### PBE (Population Branch Excess)

- **Not a Z-score.** PBE = PBS − T23 × median(PBS) / median(T23), with medians genome-wide or per chromosome
- Positive PBE: more pop1 branch length than that background relationship predicts
- Negative PBE: less pop1-specific divergence than expected given T23
- Downstream outlier calling (`identify_outliers.sh`) uses the collated PBE quantiles, not a fixed ±2 cutoff

## Notes

- Window sizes should be chosen based on your genome size and expected signal
  - Smaller windows (1-5 kb): Higher resolution, more noise
  - Larger windows (10-50 kb): Lower resolution, less noise
- Step size controls overlap between windows
  - Step size = window size: No overlap
  - Step size < window size: Overlapping windows (smoother results)
- Coverage filters are important for poolseq data
  - Too low: Poor quality estimates
  - Too high: May indicate repetitive regions or mapping artifacts
- PBS requires three populations with a clear evolutionary relationship
  - Pop1: Target population (e.g., adapted)
  - Pop2: Comparison population (e.g., control)
  - Pop3: Outgroup (e.g., ancestral)

## Troubleshooting

### Sync File Format Issues

If you get errors about sync file format, verify:
1. File is tab-delimited
2. First three columns are: chr, pos, ref
3. Population columns contain allele counts in format: `A:T:G:C:N:del`

### Memory Issues

For large genomes:
- Use larger window sizes to reduce memory usage
- Process chromosomes separately
- Increase Java heap size: `export JAVA_OPTS="-Xmx8g"`

### PBS/PBE errors

- `Could not find FST columns for all three pairs`: names must match the FST header (`Echo:Kjer.fst`). Use the `sample_name` values from `sample_info.csv`
- `--fst-dir` sees no files: `calculate_fst.sh` writes `fst_w{W}_s{S}.{csv|tsv}` or `fst_single.{csv|tsv}` in that directory
- Collate missing PBE: pass `--pbe-dir` to the directory of `*_pbe.tsv` (not the HDF5 from `--output-format hdf5`)

## Citation

If you use these scripts in your research, please cite:
- This repository
- Popoolation2: Kofler et al. (2011) Bioinformatics 27:3436-3437
- Original PBS method: Yi et al. (2010) Science 328:1120-1124

