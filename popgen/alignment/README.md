# Poolseq Processing Pipeline for Schistocephalus solidus

This repository contains a reproducible pipeline for processing Illumina 2x150 paired-end poolseq data from genomic DNA of pooled and homogenized *Schistocephalus solidus* individuals.

## Overview

This pipeline was initially developed by Jesse Weber (JW) to identify Fst outliers between Kjer and Echo populations, which were then used to generate primers for genotyping F1/F2 hybrid worms. The pipeline processes pooled sequencing data through the following steps:

1. **Adapter trimming** (BBtools bbduk)
2. **Quality correction and overlap assessment** (BBtools bbmerge)
3. **Read merging** for overlapping paired-end reads
4. **Mapping** with BWA mem
5. **Deduplication** (samtools collate, fixmate, sort, markdup)
6. **BAM merging** of merged and unmerged reads, read-group tagging, optional MAPQ / primary-only filtering
7. **Optional** mpileup and sync file

Companion scripts in this folder then measure per-window alignment quality (`seq_qual_metrics.sh`), build the repeat + alignment-quality mask (`build_mask.sh`, `repeatmasker_bed.sh`) and help choose a MAPQ cutoff (`mapq_threshold_counts.sh`).

## Citation

If you use this pipeline in your research, please cite this repository and the original workflow.

## Requirements

### Software Dependencies

**Required:**
- **BWA** (Burrows-Wheeler Aligner) - for read mapping
- **samtools** - for BAM file processing
- **BBtools** (and Java) - for adapter trimming and read merging
  - Download from: https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/
- **tmux** or **screen** - for `--parallel` (falls back to background jobs)

All of these (except tmux/screen) are in the `popgen` conda environment.

**For Population Genetics Analyses (Recommended):**
- **grenedalf** - Primary tool for population genetics statistics (200x faster than popoolation2)
  - Install via Bioconda: `conda install bioconda::grenedalf`
  - GitHub: https://github.com/lczech/grenedalf
  - Wiki: https://github.com/lczech/grenedalf/wiki
  - Advantages: Direct BAM input, faster processing, lower memory usage

**Optional:**
- **popoolation2** - only for `--popoolation2` sync file generation
  - Available at: https://sourceforge.net/projects/popoolation2/

### Reference genome and repeat masking

**Recommended practice for this pipeline:** align reads to the **unmasked** assembly, then exclude repeats when computing population-genetics statistics—not by hard-masking the reference before BWA.

| Assembly | Role in this repo |
|----------|-------------------|
| **Unmasked** (or soft-masked; BWA treats both the same) | `--reference` in `process_poolseq.sh`; same FASTA for `snp_call.sh` / `indel_call.sh` |
| **Combined mask BED** (`mask.bed` from [`build_mask.sh`](#combined-mask-build_masksh): repeats + low-quality alignment windows) | `--filter-mask-total-bed` in `calculate_pi_theta.sh` / `calculate_fst.sh` / `snp_call.sh` / `indel_call.sh` |

**`.out` → BED:** convert RepeatMasker output with [`repeatmasker_bed.sh`](repeatmasker_bed.sh) (merges overlaps with bedtools by default):

```bash
./repeatmasker_bed.sh -i genome.fa.out -o genome.out.bed
# If BAM/VCF use "chr" prefix but .out does not:
./repeatmasker_bed.sh -i genome.fa.out -o genome.out.bed --chr-prefix chr
```

| **Unmasked + GFF** | `--reference` in `variant_csq.sh` (codon translation; must match GFF coordinates) |

**Why not hard-mask before alignment?** Masking repeats as `N` before mapping can push reads to wrong loci (paralogs), adding false signal. It is better to map honestly, then drop low-confidence alignments and repeat intervals at the stats stage.

**Soft vs hard masking.** A soft-masked FASTA is only a carrier for the repeat annotation: BWA, samtools, grenedalf and bcftools ignore case, so it behaves as unmasked. Hard masking (`N`) is never needed: every downstream step takes the mask as a BED, and grenedalf counts masked bp as missing, so they leave the valid-loci denominator rather than counting as invariant sites. `variant_csq.sh` must have the unmasked sequence to translate codons.

Recommended order:

1. Align to the unmasked (or soft-masked) reference (`process_poolseq.sh`).
2. Convert repeats to BED (`repeatmasker_bed.sh`, or pass `.out` / soft-masked FASTA straight to `build_mask.sh`).
3. Run `seq_qual_metrics.sh` at a finely overlapping scale for the mask (e.g. `-w 1000 --step-size 250`).
4. Build one mask: `build_mask.sh` → `mask.bed`.
5. Pass the same `mask.bed` and the same `--reference-genome` to diversity, FST, SNP and indel calling.
6. Run `seq_qual_metrics.sh` and `variant_csq_windows.sh` at the analysis scales with the same reference, then collate (see [window boundaries](../stats/README.md#window-boundaries)).

**Post-alignment BAM filtering** (`process_poolseq.sh`, after merge):

- `--bam-min-mapq 20` — keep reads with MAPQ ≥ 20 in the final BAM.
- `--bam-primary-only` — remove secondary and supplementary alignments (`samtools view -F 2304`).
- The unfiltered BAM is kept as `*.dedup.unfiltered.bam` for `seq_qual_metrics.sh` / `build_mask.sh`.

Example (recommended flags on top of your usual run):

```bash
./process_poolseq.sh \
  --reference /path/to/ref.unmasked.fa \
  --bam-min-mapq 20 \
  --bam-primary-only \
  ... other options ...
```

First pass mapping in this project used:

- `worm_q10.medaka.purged.fa`

The script will automatically index the reference genome with BWA and `samtools faidx` if not already indexed.

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

### Recommended run (all samples)

Start the launcher inside tmux (it stays running until every sample has a slot), with the `popgen` conda environment active so BBtools, BWA and samtools are found and each job re-activates it:

```bash
conda activate popgen
./process_poolseq.sh --parallel \
  --sample-info ../../sample_info.csv \
  --reference /path/to/ref.unmasked.fa \
  --output-dir /path/to/output \
  --bam-min-mapq 20 --bam-primary-only \
  --threads 6 --parallel-max-jobs 2
```

- `--sample-info` CSV columns: `sample_name, read1_file, read2_file, pool_size` (relative paths are resolved from the CSV's directory). Each finished sample writes its final BAM path into a `bam_file` column (5th column), which `seq_qual_metrics.sh`, `calculate_pi_theta.sh`, `calculate_fst.sh`, `snp_call.sh` and `indel_call.sh` read.
- Each sample runs in its own tmux (or screen) session `poolseq_<sample>` with `--threads` cores. At most `--parallel-max-jobs` run at once (default: CPU cores / threads); the launcher waits for a free slot before starting the next sample. BBtools memory defaults to 80% of available RAM divided by the number of concurrent jobs (`--bbtools-memory` to override).
- Logs: `./log/<sample>_<time>.log` per sample, `./log/parallel_<time>.log` for the launcher, and `./log/poolseq_<sample>.exit` with each job's exit code (0 = success). Monitor with `tmux ls` / `tmux attach -t poolseq_<sample>`.
- BBtools and `adapters.fa` are found from the conda environment when `--bbtools-dir` / `--adapters` are not given. Without an adapters file, steps 1-3 (trimming, error correction, merging) are skipped.

### One sample

```bash
./process_poolseq.sh \
  --sample-name Echo_Pool_S1 \
  --read1 /path/to/Echo_Pool_S1_L001_R1_001.fastq.gz \
  --read2 /path/to/Echo_Pool_S1_L001_R2_001.fastq.gz \
  --reference /path/to/ref.unmasked.fa \
  --output-dir /path/to/output/Echo_Pool_S1 \
  --bam-min-mapq 20 --bam-primary-only \
  --threads 16
```

### Re-running and resuming

Re-running the same command resumes an interrupted sample. Steps 1-5 record completion in `{work_dir}/.done/` and are re-run if their marker is missing (an interrupted step is redone, never skipped with a truncated file); once one step re-runs, all later steps re-run too. If the merged BAM already exists, steps 1-5 are skipped entirely; delete it to redo the alignment. BAMs, the mpileup and the sync file are written to temporary names and renamed when complete.

### Command-Line Options

#### Required

- `-n, --sample-name NAME`: Sample name; also the BAM read group (`@RG ID/SM`). Not needed with `--parallel` (taken from the CSV).
- `-1, --read1 FILE` / `-2, --read2 FILE`: R1 / R2 fastq.gz (single-sample mode)
- `-r, --reference FILE`: Reference FASTA, unmasked or soft-masked, never hard-masked (see [Reference genome and repeat masking](#reference-genome-and-repeat-masking)). Use an uncompressed FASTA: the same file is given to the downstream scripts, which need its `.fai` (created here if missing). BWA and `.fai` indexes are built once before parallel jobs start.
- `-o, --output-dir DIR`: Output directory (per-sample subdirectories in parallel mode)

#### Read processing

- `-a, --adapters FILE` / `-b, --bbtools-dir DIR`: adapters.fa and BBtools directory (default: found in PATH / conda env)
- `--min-overlap N` (20), `--min-overlap0 N` (15): bbmerge overlap settings
- `--trimq N` (20): quality trimming in the merge step
- `--bbtools-memory SIZE`: Java heap for BBtools (e.g. `8g`)

#### Final BAM

- `--bam-min-mapq N`: keep reads with MAPQ >= N in the final BAM (default 0 = off; 20 recommended, since grenedalf diversity/FST have no MAPQ filter of their own)
- `--bam-primary-only`: drop secondary and supplementary alignments
- With either filter, the unfiltered BAM is kept as `*.dedup.unfiltered.bam`. `seq_qual_metrics.sh` uses it automatically, because the mask's MAPQ-0 signal is exactly what the filter removes (a repeat collapsed onto one locus shows as zero coverage in the filtered BAM, indistinguishable from a deletion).

#### Optional outputs

- `--mpileup`: also write a per-sample mpileup (off by default; no script in this repo reads it)
- `--mpileup-min-qual N` (20): minimum base quality for the mpileup and sync file
- `--use-grenedalf-sync` (with `-g, --grenedalf PATH` or grenedalf in PATH): write a per-sample sync file from the final BAM
- `-p, --popoolation2 JAR`: write the sync file with popoolation2 mpileup2sync instead (implies `--mpileup`)
- `--skip-sync`: never write a sync file

#### Running

- `-t, --threads N` (8): threads per sample
- `-w, --work-dir DIR`: intermediate files (default `{output_dir}/work`)
- `--parallel`, `--sample-info FILE`, `--parallel-max-jobs N`, `--use-tmux`, `--use-screen`: see [Recommended run](#recommended-run-all-samples)
- `--file-prefix PREFIX`: prefix for output file names
- `--dry-run`: print the commands without running them

### Output Files

In the output directory (per sample):

| File | Contents |
|------|----------|
| `{sample}_All_seq.dedup.bam` (+ `.bai`) | Final BAM: deduplicated, read group `SM:{sample}`, filtered if `--bam-min-mapq` / `--bam-primary-only`. Recorded in the CSV `bam_file` column. |
| `{sample}_All_seq.dedup.unfiltered.bam` (+ `.bai`) | Only with a filter: the same BAM before filtering, for `seq_qual_metrics.sh` / `build_mask.sh` |
| `{sample}_All_seq.dedup.mpileup` | Only with `--mpileup` / `--popoolation2` |
| `{sample}.sync` | Only with `--use-grenedalf-sync` / `--popoolation2`; column named `{sample}` |

Intermediate files in the work directory (safe to delete once the final BAM exists): `cleaned/` (trimmed and error-corrected reads), `merged/` (merged reads and unmerged pairs), `mapped/` (uncompressed SAM from BWA, the largest files), `dedup/` (deduplicated BAMs), `.done/` (step markers).

## Pipeline Steps

### Step 1: Adapter Trimming
Removes Illumina adapters using BBtools `bbduk.sh` with the following parameters:
- `ktrim=r`: Trim adapters from the right (3' end)
- `ftm=5`: Force trim mode
- `k=23`, `mink=11`: Kmer size for matching
- `hdist=2`: Hamming distance
- `tbo`, `tpe`: Trim based on pair overlap

### Step 2: Quality Correction
Performs quality correction and overlap assessment using BBtools `bbmerge.sh` with `ecco` and `mix` options. This step corrects base calls in the overlapping part of read pairs to account for sequencing errors.

**Example output statistics:**
- Echo: 29.864% of pairs joined, 70.136% no solution
- Kjer: 40.863% of pairs joined, 59.137% no solution
- Average insert size: ~224 bp
- Insert range: 20-286 bp

### Step 3: Read Merging
Merges overlapping paired-end reads using `bbmerge.sh` with quality trimming (`trimq=20`). This is critical for poolseq data because overlapping reads cannot be counted as independent observations.

### Step 4: Mapping
Maps merged reads (single-end) and unmerged pairs separately with BWA mem `-M` (shorter split hits marked secondary, removed by `--bam-primary-only`).

### Step 5: Deduplication
Each SAM goes through `samtools collate` → `fixmate -m` → `sort` → `markdup -r` (duplicates removed).

### Step 6: BAM Merging, read groups and filtering
The two deduplicated BAMs are merged and every read is tagged with one read group (`ID`, `SM` = sample name; `samtools addreplacerg -w`), which freebayes (`indel_call.sh`) needs to assign reads to pools. A re-run re-tags an existing BAM whose reads lack `RG:Z`. With `--bam-min-mapq` / `--bam-primary-only`, the filtered final BAM is then derived from the unfiltered one (and rebuilt whenever the unfiltered BAM is newer). The final BAM path is written to the sample CSV.

### Step 7: Mpileup and sync (optional)
Written only on request (`--mpileup`, `--use-grenedalf-sync`, `--popoolation2`); the downstream scripts in this repo read the BAMs directly.

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

**Output:** Per-sample TSV under `{output_dir}/{sample}/mapq_threshold_counts.tsv` and a combined `mapq_threshold_counts_all.tsv` with columns: `mapq_threshold`, `n_mapped_reads`, `fraction_retained` (and optionally `estimated_mean_depth`). Plot these vs threshold to pick a cutoff for `--bam-min-mapq` (20 is the usual choice, and `indel_call.sh` defaults to it). Run it on the `*.dedup.unfiltered.bam` files (via `--bam`), since a filtered BAM has nothing below its cutoff.

### Per-window coverage and mapping quality (seq_qual_metrics.sh)

**`seq_qual_metrics.sh`** computes per-window alignment-quality metrics with samtools + awk. Output TSV is compatible with collate (in popgen/stats/) for merging into diversity HDF5 via `--seq-qual-dir`.

- **Input**: BAM file(s) via `--bam` or `--sample-info` CSV (columns: sample_name, bam_file). When `X.unfiltered.bam` exists beside a listed `X.bam` (kept by `process_poolseq.sh` with a MAPQ / primary-only filter), it is used instead, so `frac_mapq0` and `mean_coverage` count the reads the filter removed; `--no-unfiltered` turns this off. Missing BAM indexes are created.
- **Output**: `{output_dir}/{sample}/seq_qual_metrics_w{W}_s{S}.tsv`, one per window/step pair. Reads with flags 0x704 (unmapped, secondary, QC-fail, duplicate) are ignored; supplementary alignments add coverage but are not counted as reads.

| Column | Meaning |
|--------|---------|
| `mean_coverage` | Aligned bases (CIGAR M/=/X) per window bp, all MAPQ — comparable to grenedalf's default (MAPQ ≥ 0) |
| `mean_mapping_quality` | Base-weighted mean MAPQ |
| `mean_coverage_hq` | Same as `mean_coverage`, restricted to MAPQ ≥ `--min-mapq` (default 20) |
| `frac_mapq0` | Fraction of aligned bases with MAPQ 0 (multi-mapping / repeats / paralogs) |
| `breadth` | Fraction of window bp with MAPQ ≥ `--min-mapq` depth ≥ `--min-depth` (default 1) |
| `rel_coverage` | `mean_coverage` / sample median window coverage (≈2 suggests collapsed duplication, ≈0 deletion/absence) |
| `n_reads` | Primary reads starting in the window |
| `frac_proper_pair` | Proper-pair fraction among paired primary reads starting in the window |
| `frac_mate_other_chr` | Fraction of those whose mate maps to another chromosome (SV / misassembly signal) |

- **Windows** follow grenedalf's interval grid (start = 1 + k·step, last window is the first reaching the chromosome end; coverage at the end is averaged over the bp inside the chromosome), so collate's exact (chr, start, end) join matches. Windows overlap correctly when step < window.
- **Parallelism**: the genome is cut into ~4×`--threads` chunks (≥ 1 Mb; small contigs packed together). `--threads` chunks per sample run concurrently, each a `samtools view` + `samtools depth` pass into bins of gcd(all window/step sizes) bp. All scales are then built from one sorted bin table, so additional `-w/--step-size` pairs cost almost nothing. `--max-jobs` runs samples in parallel; total processes ≈ 2 × threads × max-jobs.
- **Options**: `--reference-genome` (recommended: the .fai sets chromosome order, and BAM `@SQ` lengths are checked against it so a BAM aligned to another assembly fails instead of producing shifted windows), `--window-size`/`--step-size` (repeatable, paired; default step = window/2), `--min-mapq`, `--min-depth`, `--threads`, `--max-jobs`, `--dry-run`. Logs: `{output_dir}/log/seq_qual_metrics_YYYYmmdd_HHMMSS.log`. A failure in any chunk or sample makes the script exit non-zero.

Use the same window/step as your diversity runs (e.g. `--window-size 1000 --step-size 500`) so collate can join qual data to diversity windows.

### Combined mask (build_mask.sh)

**`build_mask.sh`** merges repeat annotations and poorly mapping windows into one BED for `--filter-mask-total-bed`.

```bash
./seq_qual_metrics.sh --sample-info ../../sample_info.csv --reference-genome ref.unmasked.fa \
  -w 1000 --step-size 250 -o seq_qual_mask -t 8 -j 4
./build_mask.sh --reference ref.unmasked.fa --repeatmasker-out ref.fa.out \
  --seq-qual-dir seq_qual_mask --output-dir mask
```

- **Repeat sources** (repeatable, any combination): `--repeat-bed`, `--repeatmasker-out` (converted with `repeatmasker_bed.sh`), `--softmask-fasta` (lowercase runs).
- **Alignment-quality source**: `--seq-qual-dir` with `seq_qual_metrics.sh` output. By default the scale with the smallest step is used. A window is flagged in a sample when `frac_mapq0 > 0.25`, `rel_coverage > 2` (collapsed duplication), `rel_coverage < 0.2` (deletion or unmappable), or `breadth < 0.5`; `--min-mean-mapq` is off by default. A window flagged in any pool is flagged for all pools (`--min-sample-fraction 0`), which is conservative and suits a few pools. With many pools, raise it (e.g. 0.5) so one pool's dropout does not remove the region everywhere; `seq_qual_flagged_windows.tsv` shows how many pools flagged each window.
- **From windows to bases**: the metrics are window averages, so a problem straddling a boundary between non-overlapping windows is split between them and can stay under the threshold in both (a 500 bp all-MAPQ-0 stretch across the 1000 bp boundary gives 0.25 in each window and is missed). Overlapping windows catch it, because one window contains it whole. Masking every flagged window in full would then over-mask, so a base is masked only when at least `--min-window-vote` (0.5) of the windows covering it are flagged. With a step of w/4, the mask typically covers the problem plus about one step on each side. With a step of w/2, each base is in only two windows, so the vote equals masking any flagged window; with w == s it misses straddling problems (both cases log a warning). `--min-window-vote 0.01` masks every flagged window in full. Smaller windows give a finer mask but noisier metrics.
- **Outputs** in `--output-dir`:

| File | Contents |
|------|----------|
| `mask.bed` | Merged mask, in reference (.fai) order, clipped to chromosome ends |
| `mask_sources.bed` | Unmerged intervals with a 4th column `repeat` or `seq_qual` |
| `seq_qual_flagged_windows.tsv` | Masked windows: `n_flagged` / `n_samples` and the reasons |
| `mask_summary.tsv` | Per chromosome and genome: length, repeat bp, seq_qual bp, masked bp, fraction masked |

Chromosomes absent from the reference are dropped with a note. `--merge-distance` and `--min-length` control merging and tiny intervals.

## Processing Multiple Samples Sequentially

Without `--parallel`, loop over samples (passing `--sample-info` still records each BAM in the CSV):

```bash
for sample in Echo_Pool_S1 Kjer_Pool_S4; do
  ./process_poolseq.sh \
    --sample-name "$sample" \
    --read1 "/path/to/${sample}_L001_R1_001.fastq.gz" \
    --read2 "/path/to/${sample}_L001_R2_001.fastq.gz" \
    --reference /path/to/ref.unmasked.fa \
    --output-dir "/path/to/output/${sample}" \
    --sample-info ../../sample_info.csv \
    --bam-min-mapq 20 --bam-primary-only \
    --threads 12
done
```

## Downstream Analyses

After alignment, continue with the [recommended order](#reference-genome-and-repeat-masking): `seq_qual_metrics.sh` and `build_mask.sh` here, then diversity, FST and PBE in [popgen/stats/](../stats/README.md) and SNP / indel calling in [popgen/variants/](../variants/README.md). All of them take the same `--sample-info` CSV (its `bam_file` column) and the same reference FASTA; grenedalf reads the BAMs directly, so no mpileup or sync file is needed.

## Dry-Run Mode

`--dry-run` checks inputs and prints every command without running it:

```bash
./process_poolseq.sh --parallel --sample-info ../../sample_info.csv \
  --reference /path/to/ref.unmasked.fa --output-dir /path/to/output \
  --bam-min-mapq 20 --bam-primary-only --dry-run
```

## Notes

- Designed for Illumina 2x150 paired-end reads with inserts < 300 bp; overlapping pairs are merged so they are not counted twice in pooled allele counts.
- Disk: the uncompressed SAMs in `work/mapped/` are the largest files; with a filter, each sample keeps two final BAMs (filtered and unfiltered).

## Troubleshooting

### Reference Genome Not Indexed
The script builds the BWA index and `.fai` if needed (once, before parallel jobs start). For a gzipped reference, `samtools faidx` only works on bgzip-compressed files; decompress the FASTA instead.

### Memory Issues
- Lower `--parallel-max-jobs` (BBtools memory per job is 80% of available RAM divided by it) or set `--bbtools-memory`
- Reduce `--threads`
- Check disk space for intermediate files

### BBtools Not Found
Activate the conda environment that has BBtools, or give `--bbtools-dir` with `bbduk.sh` and `bbmerge.sh`.

### A parallel job failed
Check `./log/poolseq_<sample>.exit` and `./log/<sample>_<time>.log`. Re-running the same launcher command resumes failed samples where they stopped and skips finished ones (a `poolseq_<sample>` tmux session that is still open is not relaunched).

## License

[Add your license information here]

## Contact

[Add contact information here]

## Acknowledgments

Original workflow developed by JW for identifying Fst outliers between Kjer and Echo populations of *Schistocephalus solidus*.
