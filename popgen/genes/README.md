# Genes

`gene_report.sh` collects genes from one or more `gene_collation.sh` summaries into a single output directory. A gene that already has a folder there is left in place, so later runs, or one run with several `--genes` files, add only the genes that are new. The two summary tables are rewritten each run from every gene in the directory.

```bash
./gene_report.sh \
  --genes ../outliers/single-stat/fst/genes/outlier_gene_summary_w10000.csv \
  --genes ../outliers/single-stat/qibe/genes/outlier_gene_summary_w10000.csv \
  --hdf5-dir ../stats/hdf5s \
  --reference /path/to/ssol-nuclear.fa \
  --gff /path/to/braker.gff3 \
  --proteome-config proteomes.yml \
  --output-dir gene_report
```

`--genes` is a table with `gene_id`, `chr`, `gene_start`, and `gene_end`. Repeat it to combine lists. A sibling `outlier_gene_regions_*.csv` and `outlier_gene_sites_*.csv` are copied into each new gene folder.

Stat columns come from `--hdf5-dir`, not from the outlier run that supplied the gene. Every `diversity_*.h5`, `fst_*.h5`, `pbe_*.h5`, `qibe_*.h5`, `lrpi_*.h5`, and `single_position.h5` is read, including every window size and the single-site files. When `single_position.h5` is present it supplies the single-site statistics and the variant calls, and the separate `*_single.h5` and `variants.h5` files are not read again. The summaries are refreshed if those HDF5 files change.

A window overlaps a gene when it is on the same chromosome and the intervals overlap. A single site overlaps when its position is inside the gene. The stat key is `<scale>|<name>`, for example `w1000_s500|Echo.pi` or `single|Cheney:Echo.fst`.

`gene_summary_brief.csv` has one row per gene: `gene_id`, `chr`, `strand`, `gene_start`, `gene_end`, `n_snps`, `n_indels`, and four lists of stat keys that pass a tail:

- `stats_high_0.99`: quantile at least 0.99
- `stats_high_0.999`: quantile at least 0.999
- `stats_low_0.99`: quantile at most 0.01
- `stats_low_0.999`: quantile at most 0.001

`gene_summary_full.csv` repeats the identity and variant counts, then `n_sites`, `n_protein_altering`, `coding_effects`, `region_categories`, and `csq_consequences`. For each stat key it then has six columns: `n_high_0.99`, `n_high_0.999`, `n_low_0.99`, `n_low_0.999`, `max_value`, and `max_quantile`. The two maxima are independent, so the highest value and the highest quantile can come from different windows. A count is the number of overlapping windows, or of overlapping single sites on the `single` scale. A blank count means that stat had no overlapping window or site. A zero means it overlapped and none passed that tail.

`n_snps` and `n_indels` count variable sites in the gene span from `single_position.h5` or `variants.h5`. Consequence columns are counted from the same file. Without that file, the counts fall back to the gene table and its sites file.

Each new `<gene_id>/` folder contains:

- `gene_id.txt` and `gene_summary_row.csv`, the input row
- `windows.csv` and, when present, `sites.csv`
- `gene.fa`, the genomic span, when `--reference` is set
- `cds.fa`, the spliced CDS of the longest transcript, when `--gff` is also set
- `blast_<name>.csv` for each proteome, every HSP returned by `blastx` with gene, product, and GO filled from that proteome's annotation table

Proteome `blastx` defaults are `--evalue 1e-5`, `--max-target-seqs 25`, and `--max-hsps 5`. A requested target count below 5 is raised for the search. The proteome config is the same YAML or JSON used by `gene_collation.sh`. Sequence and blastx run only for genes added in that invocation.
