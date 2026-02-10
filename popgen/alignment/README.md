# Poolseq Processing Pipeline for Schistocephalus solidus

This repository contains a reproducible pipeline for processing Illumina 2x150 paired-end poolseq data from genomic DNA of pooled and homogenized *Schistocephalus solidus* individuals.

## Overview

This pipeline was initially developed by Jesse Weber (JW) to identify Fst outliers between Kjer and Echo populations, which were then used to generate primers for genotyping F1/F2 hybrid worms. The pipeline processes pooled sequencing data through the following steps:

1. **Adapter trimming** (BBtools bbduk)
2. **Quality correction and overlap assessment** (BBtools bbmerge)
3. **Read merging** for overlapping paired-end reads
4. **Mapping** with BWA mem
5. **Deduplication** (samtools collate, fixmate, sort, markdup)
6. **BAM merging** of merged and unmerged reads
7. **Variant calling preparation** (mpileup and optional sync file generation)

## Citation

If you use this pipeline in your research, please cite this repository and the original workflow.

## Requirements

### Software Dependencies

**Required:**
- **BWA** (Burrows-Wheeler Aligner) - for read mapping
- **samtools** - for BAM file processing
- **BBtools** - for adapter trimming and read merging
  - Download from: https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/

**For Population Genetics Analyses (Recommended):**
- **grenedalf** - Primary tool for population genetics statistics (200x faster than popoolation2)
  - Install via Bioconda: `conda install bioconda::grenedalf`
  - GitHub: https://github.com/lczech/grenedalf
  - Wiki: https://github.com/lczech/grenedalf/wiki
  - Advantages: Direct BAM input, faster processing, lower memory usage

**Optional (Fallback):**
- **Java** - required for popoolation2 (optional, for backward compatibility)
- **popoolation2** (optional) - for sync file generation and legacy analyses
  - Available at: https://sourceforge.net/projects/popoolation2/

### Reference Genome

First pass mapping was done to an temporary genome file:
- `worm_q10.medaka.purged.fa`

The script will automatically index the reference genome if not already indexed.

## Installation

1. Clone this repository:
```bash
git clone <repository-url>
cd ssol-poolseq
```

2. Ensure all required software is installed and in your PATH:
```bash
# Check BWA
bwa

# Check samtools
samtools --version

# Check BBtools (if using)
bbduk.sh --version

# Check grenedalf (recommended for analyses)
grenedalf --version

# Or install grenedalf via conda:
conda install bioconda::grenedalf
```

3. Make the script executable (if not already):
```bash
chmod +x process_poolseq.sh
```

## Usage

### Basic Usage

**With grenedalf for sync file generation (recommended):**
```bash
./process_poolseq.sh \
  --sample-name Echo_Pool_S1 \
  --read1 /path/to/Echo_Pool_S1_L001_R1_001.fastq.gz \
  --read2 /path/to/Echo_Pool_S1_L001_R2_001.fastq.gz \
  --reference /path/to/GCA_017591395.1_ASM1759139v1_genomic.fna.gz \
  --output-dir /path/to/output \
  --adapters /path/to/bbmap/resources/adapters.fa \
  --bbtools-dir /path/to/bbtools \
  --grenedalf /usr/bin/grenedalf \
  --use-grenedalf-sync \
  --threads 16
```

**With popoolation2 for sync file generation (fallback):**
```bash
./process_poolseq.sh \
  --sample-name Echo_Pool_S1 \
  --read1 /path/to/Echo_Pool_S1_L001_R1_001.fastq.gz \
  --read2 /path/to/Echo_Pool_S1_L001_R2_001.fastq.gz \
  --reference /path/to/GCA_017591395.1_ASM1759139v1_genomic.fna.gz \
  --output-dir /path/to/output \
  --adapters /path/to/bbmap/resources/adapters.fa \
  --bbtools-dir /path/to/bbtools \
  --popoolation2 /path/to/popoolation2/mpileup2sync.jar \
  --threads 16
```

### Full Example with All Options

```bash
./process_poolseq.sh \
  --sample-name Echo_Pool_S1 \
  --read1 /path/to/Echo_Pool_S1_L001_R1_001.fastq.gz \
  --read2 /path/to/Echo_Pool_S1_L001_R2_001.fastq.gz \
  --reference /path/to/GCA_017591395.1_ASM1759139v1_genomic.fna.gz \
  --output-dir /path/to/output \
  --work-dir /path/to/work \
  --adapters /path/to/bbmap/resources/adapters.fa \
  --bbtools-dir /path/to/bbtools \
  --popoolation2 /path/to/popoolation2/mpileup2sync.jar \
  --threads 16 \
  --min-overlap 20 \
  --min-overlap0 15 \
  --trimq 20 \
  --mpileup-min-qual 20
```

### Command-Line Options

#### Required Options

- `-n, --sample-name NAME`: Sample name (e.g., Echo_Pool_S1)
- `-1, --read1 FILE`: Path to R1 fastq.gz file
- `-2, --read2 FILE`: Path to R2 fastq.gz file
- `-r, --reference FILE`: Path to reference genome (FASTA, can be .gz)
- `-o, --output-dir DIR`: Output directory for final results

#### Required for Full Pipeline

- `-a, --adapters FILE`: Path to adapters.fa file (BBtools resource)
- `-b, --bbtools-dir DIR`: Directory containing BBtools scripts (bbduk.sh, bbmerge.sh)

#### Optional Options

- `-w, --work-dir DIR`: Working directory for intermediate files (default: `output-dir/work`)
- `-t, --threads N`: Number of threads (default: 8)
- `-p, --popoolation2 JAR`: Path to popoolation2 mpileup2sync.jar (optional, for sync file generation)
- `--min-overlap N`: Minimum overlap for bbmerge (default: 20)
- `--min-overlap0 N`: Minimum overlap0 for bbmerge (default: 15)
- `--trimq N`: Quality threshold for trimming in merge step (default: 20)
- `--mpileup-min-qual N`: Minimum quality for mpileup (default: 20)
- `--skip-sync`: Skip sync file generation (default: false)
- `--dry-run`: Preview commands without executing (useful for testing)
- `--parallel`: Enable parallel processing for multiple samples
- `--sample-list FILE`: File with sample configurations for parallel processing
- `--parallel-max-jobs N`: Maximum concurrent jobs (default: number of CPU cores)
- `--use-tmux`: Force use of tmux for parallel jobs
- `--use-screen`: Force use of screen for parallel jobs
- `-h, --help`: Show help message

### Output Files

The pipeline generates the following output files in the specified output directory:

- `{SAMPLE_NAME}_All_seq.dedup.bam`: Final merged and deduplicated BAM file
- `{SAMPLE_NAME}_All_seq.dedup.bam.bai`: BAM index file
- `{SAMPLE_NAME}_All_seq.dedup.mpileup`: Mpileup file for variant calling
- `{SAMPLE_NAME}.sync`: Sync file (if popoolation2 is provided and `--skip-sync` is not used)

Intermediate files are stored in the work directory (default: `output-dir/work/`):
- `cleaned/`: Adapter-trimmed and quality-corrected reads
- `merged/`: Merged reads and unmerged read pairs
- `mapped/`: SAM files from BWA mapping
- `dedup/`: Deduplicated BAM files

## Pipeline Steps

### Step 1: Adapter Trimming
Removes Illumina adapters using BBtools `bbduk.sh` with the following parameters:
- `ktrim=r`: Trim adapters from the right (3' end)
- `ftm=5`: Force trim mode
- `k=23`, `mink=11`: Kmer size for matching
- `hdist=2`: Hamming distance
- `tbo`, `tpe`: Trim based on pair overlap

### Step 2: Quality Correction
Performs quality correction and overlap assessment using BBtools `bbmerge.sh` with `ecco` and `mix` options. This step corrects genotype calls for overlapping reads to account for sequencing errors.

**Example output statistics:**
- Echo: 29.864% of pairs joined, 70.136% no solution
- Kjer: 40.863% of pairs joined, 59.137% no solution
- Average insert size: ~224 bp
- Insert range: 20-286 bp

### Step 3: Read Merging
Merges overlapping paired-end reads using `bbmerge.sh` with quality trimming (`trimq=20`). This is critical for poolseq data because overlapping reads cannot be counted as independent observations.

### Step 4: Mapping
Maps both merged and unmerged reads to the reference genome using BWA mem with the `-M` flag (mark shorter split hits as secondary).

### Step 5: Deduplication
Processes SAM files through:
1. `samtools collate`: Collates reads by name
2. `samtools fixmate`: Fixes mate information
3. `samtools sort`: Sorts reads
4. `samtools markdup`: Marks duplicates

### Step 6: BAM Merging
Merges the deduplicated BAM files from merged and unmerged reads into a single BAM file with proper read group headers.

### Step 7: Mpileup Generation
Creates mpileup files for downstream variant calling and population genetics analyses.

### Quality threshold assessment (MAPQ)

To choose a MAPQ cutoff that balances sequence quality and mapped depth, run the companion script **`mapq_threshold_counts.sh`** on your final BAMs. It reports the number of mapped reads at multiple MAPQ thresholds (one pass per BAM, histogram then cumulative-sum).

**Run on existing BAMs (e.g. after process_poolseq):**
```bash
./mapq_threshold_counts.sh \
  -i /path/to/sample_info.csv \
  -o mapq_stats \
  --mapq-thresholds "0,5,10,15,20,25,30,40,50"
```

Optional: `--reference ref.fa` (with `ref.fa.fai`) adds an `estimated_mean_depth` column; `--primary-only` counts only primary alignments.

**Output:** Per-sample TSV under `{output_dir}/{sample}/mapq_threshold_counts.tsv` and a combined `mapq_threshold_counts_all.tsv` with columns: `mapq_threshold`, `n_mapped_reads`, `fraction_retained` (and optionally `estimated_mean_depth`). Plot these vs threshold to pick a cutoff; the pipeline mpileup (Step 7) and downstream variant calling typically use **MAPQ ≥ 20** by default.

### Per-window coverage and mapping quality (seq_qual_metrics.sh)

**`seq_qual_metrics.sh`** computes per-window averages for **coverage** and **mapping quality** using samtools. Output TSV is compatible with collate (in popgen/stats/) for merging into diversity HDF5 via `--seq-qual-dir`.

- **Input**: BAM file(s) via `--bam` or `--sample-info` CSV (columns: sample_name, bam_file).
- **Output**: `{output_dir}/{sample}/seq_qual_metrics_w{W}_s{S}.tsv` (chr, start, end, sample, mean_coverage, mean_mapping_quality).
- **Options**: `--reference-genome` (optional; use .fai for chromosome order), `--window-size`, `--step-size` (default: window/2), `--threads`. Logs: `{output_dir}/log/seq_qual_metrics_YYYYmmdd_HHMMSS.log`.

Use the same window/step as your diversity runs (e.g. `--window-size 1000 --step-size 500`) so collate can join qual data to diversity windows.

### Optional: Sync File Generation
Sync files can be generated using either grenedalf (recommended, faster) or popoolation2 (fallback, for backward compatibility). If `--use-grenedalf-sync` is specified and grenedalf is available, grenedalf will be used. Otherwise, if popoolation2 JAR is provided, it will be used. Sync files are optional if using grenedalf directly with BAM files for downstream analyses.

## Processing Multiple Samples

### Parallel Processing (Recommended)

The script supports parallel processing using tmux or screen sessions. This allows you to process multiple samples simultaneously while monitoring each job independently.

**1. Create a sample list file** (tab-separated format):
```bash
# sample_list.txt
# Format: sample_name<TAB>r1_file<TAB>r2_file<TAB>output_dir
Echo_Pool_S1	/path/to/Echo_Pool_S1_R1.fastq.gz	/path/to/Echo_Pool_S1_R2.fastq.gz	/path/to/output/Echo_Pool_S1
Kjer_Pool_S4	/path/to/Kjer_Pool_S4_R1.fastq.gz	/path/to/Kjer_Pool_S4_R2.fastq.gz	/path/to/output/Kjer_Pool_S4
Outgroup_Pool_S2	/path/to/Outgroup_Pool_S2_R1.fastq.gz	/path/to/Outgroup_Pool_S2_R2.fastq.gz	/path/to/output/Outgroup_Pool_S2
```

**2. Run with parallel processing:**
```bash
./process_poolseq.sh \
  --parallel \
  --sample-list sample_list.txt \
  --reference /path/to/reference.fna.gz \
  --adapters /path/to/adapters.fa \
  --bbtools-dir /path/to/bbtools \
  --grenedalf /usr/bin/grenedalf \
  --use-grenedalf-sync \
  --threads 16 \
  --parallel-max-jobs 4
```

**Options:**
- `--parallel`: Enable parallel processing mode
- `--sample-list FILE`: Path to tab-separated sample list file
- `--parallel-max-jobs N`: Maximum concurrent jobs (default: number of CPU cores)
- `--use-tmux`: Force use of tmux (auto-detected if available)
- `--use-screen`: Force use of screen (auto-detected if available)

**Monitoring parallel jobs:**
- **tmux**: `tmux ls` to list sessions, `tmux attach -t poolseq_<sample_name>` to attach
- **screen**: `screen -ls` to list sessions, `screen -r poolseq_<sample_name>` to attach
- Each job logs to: `{output_dir}/{sample_name}.log`

### Sequential Processing (Alternative)

To process multiple samples sequentially, you can create a simple loop:

```bash
for sample in Echo_Pool_S1 Kjer_Pool_S4; do
  ./process_poolseq.sh \
    --sample-name "$sample" \
    --read1 "/path/to/${sample}_L001_R1_001.fastq.gz" \
    --read2 "/path/to/${sample}_L001_R2_001.fastq.gz" \
    --reference /path/to/reference.fna.gz \
    --output-dir "/path/to/output/${sample}" \
    --adapters /path/to/adapters.fa \
    --bbtools-dir /path/to/bbtools \
    --threads 16
done
```

## Creating Combined Mpileup for Multiple Populations

After processing individual samples, you can create a combined mpileup for population comparisons:

```bash
samtools mpileup -B \
  /path/to/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  /path/to/Kjer_Pool_S4/Kjer_Pool_S4_All_seq.dedup.bam \
  -o Echo_Kjer_combined.mpileup
```

## Downstream Analyses

After processing your samples, you can perform population genetics analyses using the scripts in the `analyses/` directory. These scripts use **grenedalf** as the primary tool (with direct BAM input) and **popoolation2** as an optional fallback.

### Available Analysis Scripts

The `analyses/` directory contains scripts for calculating:

1. **π (Pi) and θ (Theta)**: Nucleotide diversity and Watterson's theta
   - Script: `analyses/calculate_pi_theta.sh`
   - Calculates genetic diversity statistics within populations
   - Supports both grenedalf (BAM input) and popoolation2 (sync file input)

2. **PBS (Population Branch Statistic)**: Population-specific branch length
   - Script: `analyses/calculate_pbs.sh`
   - Useful for detecting signals of selection
   - Supports both grenedalf (BAM input) and popoolation2 (sync file input)

### Quick Start with Grenedalf (Recommended)

**Calculate diversity statistics directly from BAM files:**
```bash
cd analyses
./calculate_pi_theta.sh \
  --bam ../output/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  --grenedalf /usr/bin/grenedalf \
  --output-dir ./diversity_results \
  --window-size 10000 \
  --pool-size 50 \
  --threads 8
```

**Calculate PBS directly from BAM files:**
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
  --threads 8
```

### Alternative: Using Popoolation2 with Sync Files (Fallback)

If grenedalf is not available, you can use popoolation2 with sync files:

1. **Create combined sync file** (if analyzing multiple populations):
```bash
# Create combined mpileup
samtools mpileup -B \
  /path/to/Echo_Pool_S1/Echo_Pool_S1_All_seq.dedup.bam \
  /path/to/Kjer_Pool_S4/Kjer_Pool_S4_All_seq.dedup.bam \
  /path/to/Outgroup_Pool_S2/Outgroup_Pool_S2_All_seq.dedup.bam \
  -o combined.mpileup

# Convert to sync
java -jar /path/to/popoolation2/mpileup2sync.jar \
  --input combined.mpileup \
  --output combined.sync \
  --fastq-type sanger \
  --min-qual 20 \
  --threads 8
```

2. **Calculate diversity statistics with popoolation2**:
```bash
cd analyses
./calculate_pi_theta.sh \
  --sync ../combined.sync \
  --popoolation2 /path/to/popoolation2 \
  --output-dir ./diversity_results \
  --window-size 10000
```

3. **Calculate PBS with popoolation2**:
```bash
./calculate_pbs.sh \
  --sync ../combined.sync \
  --popoolation2 /path/to/popoolation2 \
  --output-dir ./pbs_results \
  --pop1-index 1 \
  --pop2-index 2 \
  --pop3-index 3 \
  --pop1-name Echo \
  --pop2-name Kjer \
  --pop3-name Outgroup
```

For detailed documentation, see [`analyses/README.md`](analyses/README.md).

## Dry-Run Mode

All scripts support `--dry-run` mode, which previews commands and file operations without executing them. This is useful for:
- Testing script configurations before running
- Verifying input files and paths
- Understanding what the pipeline will do
- Debugging configuration issues

**Example:**
```bash
./process_poolseq.sh \
  --sample-name Echo_Pool_S1 \
  --read1 /path/to/Echo_Pool_S1_R1.fastq.gz \
  --read2 /path/to/Echo_Pool_S1_R2.fastq.gz \
  --reference /path/to/reference.fna.gz \
  --output-dir /path/to/output \
  --adapters /path/to/adapters.fa \
  --bbtools-dir /path/to/bbtools \
  --dry-run
```

Dry-run mode will:
- Check input file existence
- Show all commands that would be executed
- Display output files that would be created
- Skip actual command execution
- Return exit code 0 if all checks pass

## Notes

- The pipeline is designed for Illumina 2x150 paired-end reads with insert sizes <300bp
- Overlapping reads are merged to avoid double-counting in pooled samples
- The script includes checkpointing: if intermediate files exist, those steps are skipped
- All intermediate files are preserved in the work directory for troubleshooting
- Parallel processing uses tmux or screen for job management (auto-detected)
- Each parallel job creates its own log file: `{output_dir}/{sample_name}.log`

## Troubleshooting

### Reference Genome Not Indexed
The script will automatically index the reference genome if needed. This may take some time for large genomes.

### Memory Issues
If you encounter memory issues:
- Reduce the number of threads
- Process samples sequentially rather than in parallel
- Ensure sufficient disk space for intermediate files

### BBtools Not Found
Make sure the BBtools directory path is correct and contains `bbduk.sh` and `bbmerge.sh`.

## License

[Add your license information here]

## Contact

[Add contact information here]

## Acknowledgments

Original workflow developed by JW for identifying Fst outliers between Kjer and Echo populations of *Schistocephalus solidus*.
