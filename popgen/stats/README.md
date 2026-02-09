# Population Genetics Analyses for Poolseq Data

This directory contains scripts for calculating population genetics statistics from processed poolseq data. These analyses use **grenedalf** as the primary tool (with direct BAM file input) and **popoolation2** as an optional fallback (for sync file compatibility).

## Overview

The analysis scripts calculate:

1. **π (Pi)**: Nucleotide diversity - the average number of pairwise nucleotide differences per site
2. **θ (Theta)**: Watterson's theta - estimated from the number of segregating sites
3. **Tajima's D**: Test statistic comparing π and θ to detect selection and demographic events
4. **PBS**: Population Branch Statistic - measures population-specific branch length, useful for detecting selection
5. **PBE**: Population Branch Excess - Z-score normalization of PBS for identifying outliers

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
- **awk** - Required for PBS/PBE calculations
- **samtools** - Required for BAM file processing (if using grenedalf)
- **bcftools** - Required for variant_call.sh (VCF/BCF generation)
- **R** with packages **dplyr**, **readr**, **optparse**, **tidyr**, **hdf5r** - Required for collate.R

### Input Files

**Primary Method (grenedalf):**
- **BAM files** - Direct input from the main processing pipeline (`process_poolseq.sh`)
  - One BAM file per population for diversity calculations
  - Three BAM files (one per population) for PBS calculations
  - BAM files should be sorted and indexed

**Fallback Method (popoolation2):**
- **Sync files** - Generated from the main processing pipeline (`process_poolseq.sh`)
  - Sync files contain allele frequency information in the format:
    ```
    chr    pos    ref    pop1    pop2    pop3    ...
    ```
  - Where each population column contains allele counts in the format: `A:T:G:C:N:del`

## Scripts

### seq_qual_metrics.sh

Computes per-window averages for **coverage** and **mapping quality** using samtools. Output TSV is compatible with collate (chr, start, end, sample, mean_coverage, mean_mapping_quality).

- **Input**: BAM file(s) via `--bam` or `--sample-info` CSV (columns: sample_name, bam_file).
- **Output**: `{output_dir}/{sample}/seq_qual_metrics_w{W}_s{S}.tsv`.
- **Options**: `--reference-genome` (optional; use .fai for chromosome order in output), `--window-size`, `--step-size`, `--threads`.

### variant_call.sh

Generates **BCF** (variant calls) from BAM(s) using **bcftools** mpileup and call. Pipeline: bcftools mpileup → bcftools call → BCF. Uses **bcftools call -m** (multiallelic caller) for pool-seq rare variants and sites with more than two alleles.

- **Input**: BAM file(s), reference FASTA (indexed with samtools faidx).
- **Output**: BCF file (e.g. `calls.bcf`), indexed.
- **Options**: `--output` or `--output-dir`, `--ploidy` (e.g. 1 for haploid pools), `--min-mapq` (default 20), `--min-baseq` (default 20), `--max-depth` (default 1000; increase for very high-coverage pool-seq), `--skip-indels` (SNPs only), `--indels-only` (indels only).

**Typical next step (post-call filtering):**
```bash
bcftools view -i 'QUAL>=20 && DP>10' calls.bcf -Ob -o calls_filtered.bcf
bcftools index -f calls_filtered.bcf
```

**Pool-seq examples:**
- Haploid pools: use `--ploidy 1`.
- Stricter quality: `bcftools view -i 'QUAL>=30 && DP>20' calls.bcf -Ob -o calls_strict.bcf`
- Exclude very low-frequency alleles (e.g. AF &lt; 0.01): `bcftools view -i 'QUAL>=20 && DP>10 && INFO/AF>0.01' calls.bcf -Ob -o calls_af01.bcf` (requires AF in INFO; add with bcftools +fill-tags or similar if needed).
- For very high coverage, increase `--max-depth` (e.g. 2000) or set it based on average depth from `seq_qual_metrics.sh` output.

### collate.sh / collate.R

Collates diversity, FST, and PBE TSV/CSV into **HDF5** with `/windows` group(s). Computes integer rank (1..N) and 0–1 quantile per stat within (sample, window_size) for diversity, and within window for FST/PBE.

- **Input**: `--diversity-dir`, `--fst-dir`, `--pbe-dir` (at least one); optional `--seq-qual-dir`.
- **Output**: In `--output-dir`: `diversity_w{N}.h5`, `fst_w{N}.h5`, `pbe.h5` (one group per trio). Each HDF5 has group `/windows` (or `/windows_trio_*` for PBE) with chr, start, end, sample/pairs, stat values, and `*_rank`, `*_quantile`.
- **Run**: `./collate.sh --diversity-dir pi_out --fst-dir fst_out --pbe-dir pbe_out -o collated`

### calculate_pi_theta.sh

Calculates nucleotide diversity (π), Watterson's theta (θ), and Tajima's D from BAM files (grenedalf) or sync files (popoolation2 fallback).

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

### calculate_pbs.sh

Calculates Population Branch Statistic (PBS) and Population Branch Excess (PBE) from three populations using BAM files (grenedalf) or sync files (popoolation2 fallback). PBS measures population-specific branch length and is useful for detecting signals of selection.

PBS is calculated as: **PBS = (T12 + T13 - T23) / 2**

Where Tij = -ln(1 - Fst_ij) is the branch length between populations i and j.

PBE is calculated as: **PBE = (PBS - mean(PBS)) / sd(PBS)**

PBE is a Z-score normalization of PBS, making it easier to identify outliers across the genome.

#### Usage

**Primary method with grenedalf (recommended):**
```bash
./calculate_pbs.sh \
  --bam-pop1 /path/to/pop1.bam \
  --bam-pop2 /path/to/pop2.bam \
  --bam-pop3 /path/to/pop3.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir /path/to/output \
  --pop1-name Echo \
  --pop2-name Kjer \
  --pop3-name Outgroup \
  --window-size 10000
```

**Fallback method with popoolation2:**
```bash
./calculate_pbs.sh \
  --sync /path/to/populations.sync \
  --popoolation2 /path/to/popoolation2 \
  --output-dir /path/to/output \
  --pop1-index 1 \
  --pop2-index 2 \
  --pop3-index 3 \
  --pop1-name Echo \
  --pop2-name Kjer \
  --pop3-name Outgroup
```

#### Options

**Input (choose one):**
- `--bam-pop1 FILE`: BAM file for population 1 (primary method)
- `--bam-pop2 FILE`: BAM file for population 2 (primary method)
- `--bam-pop3 FILE`: BAM file for population 3 (primary method)
- `-s, --sync FILE`: Sync file (fallback method)
- `--pop1-index N`: Column index (1-based) for population 1 in sync file (fallback method)
- `--pop2-index N`: Column index (1-based) for population 2 in sync file (fallback method)
- `--pop3-index N`: Column index (1-based) for population 3 in sync file (fallback method)

**Tool selection:**
- `-g, --grenedalf PATH`: Path to grenedalf executable (default: check PATH)
- `-p, --popoolation2 DIR`: Directory containing popoolation2 scripts (for fallback)
- `--use-popoolation2`: Force use of popoolation2 instead of grenedalf

**Required:**
- `-o, --output-dir DIR`: Output directory for results

**Optional:**
- `--pop1-name NAME`: Name for population 1 (default: Pop1)
- `--pop2-name NAME`: Name for population 2 (default: Pop2)
- `--pop3-name NAME`: Name for population 3 (default: Pop3)
- `-w, --window-size N`: Window size in bp (default: 1000)
- `--step-size N`: Step size in bp (default: 1000)
- `--min-coverage N`: Minimum coverage per site (default: 4)
- `--max-coverage N`: Maximum coverage per site (default: 200)
- `--min-count N`: Minimum allele count (default: 2, popoolation2 only)
- `--pool-size N`: Pool size (haploid chromosomes) (default: 50)
- `-t, --threads N`: Number of threads (default: 1)
- `--dry-run`: Preview commands without executing (dry-run mode)

#### Understanding Population Indices

Sync files have the format:
```
chr    pos    ref    pop1    pop2    pop3    ...
```

The first three columns are: chromosome, position, reference allele. Population data starts at column 4. Therefore:
- Column 4 = Population index 1
- Column 5 = Population index 2
- Column 6 = Population index 3

However, the script uses 1-based indexing for populations, so:
- `--pop1-index 1` refers to column 4 (first population)
- `--pop2-index 2` refers to column 5 (second population)
- `--pop3-index 3` refers to column 6 (third population)

#### Output Files

- `pbs.txt`: PBS values per window
  - Columns: chr, pos, Fst_pop1_pop2, Fst_pop1_pop3, Fst_pop2_pop3, T12, T13, T23, PBS
- `pbe.txt`: Population Branch Excess (PBE) per window
  - Columns: chr, pos, Fst_pop1_pop2, Fst_pop1_pop3, Fst_pop2_pop3, T12, T13, T23, PBS, PBE
- `fst_pop1_pop2.txt`: Fst between population 1 and 2
- `fst_pop1_pop3.txt`: Fst between population 1 and 3
- `fst_pop2_pop3.txt`: Fst between population 2 and 3
- `fst_*_snps.txt`: SNP-level Fst files for each pair

#### Examples

**Using grenedalf with BAM files (recommended):**
```bash
./calculate_pbs.sh \
  --bam-pop1 ../output/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  --bam-pop2 ../output/Kjer_Pool_S4/Kjer_Pool_S4_All_seq.dedup.bam \
  --bam-pop3 ../output/Outgroup_Pool_S2/Outgroup_Pool_S2_All_seq.dedup.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir ./pbs_results \
  --pop1-name Echo \
  --pop2-name Kjer \
  --pop3-name Outgroup \
  --window-size 10000 \
  --step-size 5000 \
  --threads 8
```

**Using popoolation2 with sync file (fallback):**
```bash
./calculate_pbs.sh \
  --sync ../output/Echo_Kjer_Outgroup.sync \
  --popoolation2 /opt/popoolation2 \
  --output-dir ./pbs_results \
  --pop1-index 1 \
  --pop2-index 2 \
  --pop3-index 3 \
  --pop1-name Echo \
  --pop2-name Kjer \
  --pop3-name Outgroup \
  --window-size 10000 \
  --step-size 5000
```

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

3. **Calculate PBS using grenedalf (recommended)**:
```bash
./calculate_pbs.sh \
  --bam-pop1 /path/to/output/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  --bam-pop2 /path/to/output/Kjer_Pool_S4/Kjer_Pool_S4_All_seq.dedup.bam \
  --bam-pop3 /path/to/output/Outgroup_Pool_S2/Outgroup_Pool_S2_All_seq.dedup.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir ./pbs_analysis \
  --pop1-name Echo \
  --pop2-name Kjer \
  --pop3-name Outgroup \
  --window-size 10000 \
  --step-size 5000 \
  --threads 8
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

4. **Calculate PBS with popoolation2**:
```bash
./calculate_pbs.sh \
  --sync combined.sync \
  --popoolation2 /opt/popoolation2 \
  --output-dir ./pbs_analysis \
  --pop1-index 1 \
  --pop2-index 2 \
  --pop3-index 3 \
  --pop1-name Echo \
  --pop2-name Kjer \
  --pop3-name Outgroup \
  --window-size 10000 \
  --step-size 5000
```

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

- **Positive PBS**: Indicates population-specific divergence, potentially due to selection
- **High PBS values**: Regions with strong population-specific selection
- **Negative PBS**: Indicates shared ancestry or gene flow
- Typically, PBS > 2-3 standard deviations above the mean are considered outliers

### PBE (Population Branch Excess)

- **PBE is a Z-score**: Values are normalized to have mean=0 and sd=1
- **PBE > 2 or 3**: Strong outliers, likely under selection in the target population
- **PBE < -2 or -3**: Regions with unusually low divergence
- **Advantage over PBS**: Easier to compare across different genomic regions and identify consistent outliers

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

### PBS Calculation Errors

- Ensure all three Fst files are generated successfully
- Check that population indices match your sync file structure
- Verify windows are aligned across all three Fst files

## Citation

If you use these scripts in your research, please cite:
- This repository
- Popoolation2: Kofler et al. (2011) Bioinformatics 27:3436-3437
- Original PBS method: Yi et al. (2010) Science 328:1120-1124

