# Outliers

Three scripts, in order:

1. **`identify_outliers.sh`** flags windows from collate HDF5 (or from diversity/FST TSVs) and merges nearby windows into regions.
2. **`outlier_annotation.sh`** attaches consequence calls and assembly genes to those regions. BLAST against other assemblies is optional.
3. **`gene_collation.sh`** turns the annotation tables around so each row is an assembly gene.

Coordinates are 1-based inclusive, the same grid as collate (`start = 1 + k·step`, end not clipped).

---

## 1. identify_outliers.sh

Reads `diversity_w*_s*.h5`, `fst_w*_s*.h5`, and `pbe_w*_s*.h5` from `--hdf5-dir` (collate output). Without HDF5, `--diversity-dir` and/or `--fst-dir` read the TSVs instead. PBE is HDF5-only.

Quantiles are taken from the collate `*_quantile` columns when those columns are present (one distribution per sample, pair, or trio at that window and step). Otherwise they are computed the same way. A window is an outlier when its quantile is at or above `--high-quantile` or at or below `--low-quantile`. Diversity stats are per sample; FST is per pair; PBE is per trio. After flagging, every sample's value and quantile is filled in for that window, not only the sample that crossed the threshold.

`--window-size` keeps one window width. If that width was collated at more than one step, pass `--step-size` too; otherwise the overlapping grids would be pooled and the script stops.

**Defaults that are already on** when you only pass the two quantiles:

- `--merge-distance auto` merges windows into regions. The gap is `max(2 × window, 2000)` bp, so 2000 bp at a 1000 bp window.
- `--merge-across-samples true` builds one set of regions across samples, pairs, and trios. `false` merges within each sample or pair.
- Set `--merge-distance 0` to write windows only.

Seed-then-expand (instead of the single threshold) needs all four of `--seed-high-quantile`, `--seed-low-quantile`, `--expand-high-quantile`, `--expand-low-quantile`. Seeds are the strict tail; neighboring windows in the soft tail are added, then merged with the same gap.

Optional filters apply to seeds and to expansion, and to the single-threshold windows: `--min-depth`, `--max-depth`, `--min-mapping-quality`, `--max-mapping-quality`, `--min-snps`, `--max-snps`, `--min-total-sites`, `--max-total-sites`. Region filters are `--min-region-snps` and `--max-region-snps`. Depth and mapping quality come from the seq-qual columns collate joined onto diversity. They are not applied unless you set them.

Other limits: `--chromosome`, `--region-start` / `--region-end`, `--top-n-chromosomes`, `--min-chromosome-length`, `--sample-pairs` (FST pairs; default all), `--top-n-extreme N` (keep the N most extreme windows in each tail).

### First pass

For the collate in `../stats/hdf5s` at 1000 bp windows and 500 bp steps:

```bash
./identify_outliers.sh \
  --hdf5-dir ../stats/hdf5s \
  --output-dir direct-threshold \
  --high-quantile 0.99 \
  --low-quantile 0.01 \
  --window-size 1000 \
  --step-size 500 \
  --statistics pi,theta,tajima_d,fst,pbe \
  --verbose
```

`--step-size 500` is the analysis grid. It is required only when `hdf5s` also contains another step at window 1000; with a single step the script proceeds without it. The rest of a first pass is already the default: regions are written, merged across samples, with no depth or mapping-quality filter. `--statistics` must list `fst` and `pbe` or those files are not read. If no `pbe_*.h5` matches the window and step, PBE is skipped rather than failing the run. `single_position.h5` is not read here.

Low quantiles on FST and PBE flag the bottom tail (unusually low differentiation or PBE), not only the high tail.

### Outputs

Written under `--output-dir`. The filename encodes window, step, chromosome scope, quantiles, and `_across_samples` or `_within_samples`, plus `_merge_only` or `_seed_expand`.

- `outlier_windows_*.csv`: one row per outlier window. `outlier_stat` and `outlier_direction` (`high` / `low`, several joined with `:`) say what crossed the threshold. Columns are `sample.stat` and `sample.stat_quantile` (diversity), `sample1:sample2.fst`, `sample1:sample2:sample3.pbe`. Coverage, mapping quality, `n_snps`, `n_total`, and `n_masked` are included when collate had them.
- `outlier_regions_*.csv`: one row per merged region. `region_start` / `region_end` span the windows. `n_windows`, `region_mean_*`, `region_max_*`, and `region_max_*_quantile` summarize the windows inside the region. `outlier_annotation.sh` reads this file or the windows file.

Logs: `{output-dir}/log/identify_outliers_*.log`. `--verbose` prints the scale, whether the collate quantile column was used, and the fill counts (row counts and whether a join would duplicate windows).

---

## 2. outlier_annotation.sh

Annotates one interval table from `identify_outliers.sh`. Give `--variant-csq`, `--gff`, or both. BLAST is separate and optional.

`outlier_regions_*.csv` uses `region_start` / `region_end` (the merged span). `outlier_windows_*.csv` uses `start` / `end`; each window is annotated on its own and `n_windows` is set to 1. A site inside two overlapping windows is listed once per window. The output suffix drops either `outlier_regions_` or `outlier_windows_` from the input filename.

```bash
./outlier_annotation.sh \
  --regions direct-threshold/outlier_windows_w1000_s500_....csv \
  --variant-csq ../stats/hdf5s/single_position.h5 \
  --gff /path/to/annotation.gff3 \
  --output-dir direct-threshold/annotation
```

The same command with `outlier_regions_w1000_s500_....csv` annotates the merged spans instead of each window.

**Consequence calls** (`--variant-csq`, repeatable) are sites on this assembly, so they are not BLASTed. Each path is a `variants_csq.tsv` from `variant_csq.sh`, a directory of those TSVs, or a collate HDF5. `*_long.tsv` is ignored. With `--single-position-merged`, collate writes `single_position.h5` (group `/sites`) and does not write `variants.h5`; pass that file, or the csq TSV directories. HDF5 adds per-pool `*_AD` / `*_alt_freq` and per-site FST/PBE when those columns were collated. TSV input has the csq columns only.

**Assembly GFF** (`--gff`) is the GFF given to `variant_csq.sh`. It lists every gene whose interval overlaps the region, including genes with no variant, and adds `gene_name`, `product`, and `go_terms` to the genes named by csq.

Outputs, suffix taken from the interval filename after `outlier_regions_` or `outlier_windows_`:

| File | One row per | Contents |
| --- | --- | --- |
| `outlier_csq_sites_*.csv` | region × site | Lead columns `qseqid` (`chr:start-end`), `region_chr`, `region_start`, `region_end`, `outlier_stat`, `n_windows`, then the site (`chr`, `pos`, `ref`, `alt`, `variant_type`, `region_category`, `coding_effect`, `genes`, `genes_protein_altering`, `csq_consequences`, `bcsq`). HDF5 also adds the frequency and FST/PBE columns. |
| `outlier_csq_genes_*.csv` | region × gene | Same region lead columns, then `gene`, `n_sites`, `n_protein_altering`, `coding_effects` (for example `missense:3;synonymous:2`), `overlaps_region`, and GFF `gene_id`, `gene_name`, `product`, `go_terms`, `strand`, `gene_start`, `gene_end`. |
| `outlier_csq_summary_*.csv` | region | Site counts (`n_sites`, `n_sites_snp`, `n_sites_indel`, `n_protein_altering`, `n_effect_*`, `n_region_*`), the region columns, then `genes_csq`, `genes_protein_altering`, `genes_gff_overlap`, `products_gff`. |
| `outlier_csq_background_*.csv` | class | Genome-wide counts by (`set`, `variant_type`, `region_category`, `coding_effect`) for `set` = `all_sites` and `non_outlier` (sites in no region). `plot_outlier_csq.sh` uses this for composition and enrichment. |

Region columns from identify (`outlier_stat`, `n_windows`, `region_mean_*`, …) are copied through. Their names are lowercased.

**BLAST** (`--blast-config` plus `--reference`) extracts each region from the reference FASTA (`.fai` required) and searches each database in the config. Hits are mapped through that database's GFF, not the assembly GFF. Use this for homologs in other annotated genomes. Optional `--blast-task blastn` (the default megablast misses divergent hits), `--evalue`, `--min-query-coverage` (fraction of the region the subject must cover), `--threads`, `--parallel-dbs`.

```yaml
databases:
  - db_path: /path/to/nt
    name: nt
    gff_path: /path/to/other.gff3
    annotation_tsv: /path/to/gene_table.tsv   # optional; joined on gene_id
```

A missing GFF or annotation table warns and does not stop the run. BLAST outputs are `outlier_genes_*.csv` (one row per hit) and `outlier_genes_summary_*.csv` (one row per region × database).

---

## 3. gene_collation.sh

Annotation is one row per region or per site. Gene collation uses the assembly GFF as the gene list and attaches the regions, sites, and BLAST homologs that overlap each gene.

```bash
./gene_collation.sh \
  --gff /path/to/annotation.gff3 \
  --annotation-dir direct-threshold/annotation \
  --output-dir direct-threshold/genes
```

`--gff` is the same file passed to `variant_csq.sh` and to `outlier_annotation.sh`. Gene rows are the `gene` features, or `mRNA` / `CDS` if the file has no `gene` rows. `gene_id` is the `gene_id` attribute, then `Name`, then `ID`. The name used to match csq calls (`csq_key`) is `Name`, or `ID` with a leading `gene:` removed. `product` and `go_terms` come from `product` / `description` / `Note` and from `Ontology_term` or `Dbxref` GO ids.

`--annotation-dir` is the `outlier_annotation.sh` output. The script picks `outlier_csq_sites_*.csv`, `outlier_csq_summary_*.csv`, and when present `outlier_csq_genes_*.csv` and BLAST `outlier_genes_*.csv` (not `outlier_genes_summary_*.csv`). `--suffix` selects one run when the directory holds several. Alternatively pass `--sites`, `--genes`, `--summary`, and `--blast-genes` yourself.

Sites are assigned to a gene when the position falls inside the gene (`assigned_by=position`) and/or the csq `genes` field names it (`assigned_by=csq`, or `csq+position` when both hit). A csq gene name that is not in this GFF is dropped; `--verbose` lists those names. Intergenic sites stay in the site table from annotation and do not appear here. BLAST homologs are genes from other assemblies; they join by overlapping region (`qseqid`), not by `gene_id`.

`--include-non-outlier` also writes GFF genes that overlap no outlier region (`in_outlier=FALSE`).

Outputs, suffix copied from the annotation files:

- `outlier_gene_summary_*.csv`: one row per gene. GFF identity and coordinates, `n_regions`, `n_sites`, `n_protein_altering`, `coding_effects`, `outlier_stat`, `qseqids`, and BLAST homologs collapsed across those regions (`blast_gene_ids`, `blast_products`, `best_evalue`).
- `outlier_gene_regions_*.csv`: one row per gene × outlier region, with that region's stats and the site counts for the gene in the region.
- `outlier_gene_sites_*.csv`: one row per gene × site. `qseqids` are semicolon-collapsed when the site sits in more than one region. `assigned_by` records whether the link was position, csq, or both.
