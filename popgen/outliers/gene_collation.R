#!/usr/bin/env Rscript

###############################################################################
# gene_collation.R
#
# Reorganize outlier_annotation.sh output by assembly GFF gene, instead of by
# variant or region. The GFF is the gene universe; annotation CSVs supply the
# outlier regions, csq sites, and optional BLAST homologs (joined by region).
#
# Inputs: --gff plus either --annotation-dir (discovers outlier_csq_*.csv and
# optional outlier_genes_*.csv) or explicit --sites / --genes / --summary /
# --blast-genes. --suffix picks one annotation run when several exist.
#
# Outputs (suffix from the annotation files, or "genes"):
#   outlier_gene_summary_<suffix>.csv  One row per GFF gene that overlaps an
#     outlier region (or every gene with --include-non-outlier). GFF fields,
#     n_regions, n_windows, qseqids, outlier_stat, site counts / coding_effects, and BLAST
#     homologs collapsed across overlapping regions (other-assembly gene_ids).
#     With --blast-config, each outlier gene sequence is searched and the best
#     HSP per database adds <name>_gene/product/go/evalue plus homolog_genes,
#     homolog_products, and homolog_go. Nucleotide task defaults to blastn.
#     With --proteome-config, blastx of the spliced CDS searches each protein
#     database and adds <name>_gene/product/go/accession/evalue plus the
#     collapsed proteome_genes, proteome_products, and proteome_go columns.
#   outlier_gene_regions_<suffix>.csv  One row per (gene, region).
#   outlier_gene_sites_<suffix>.csv    One row per (gene, site): sites whose
#     position falls in the gene and/or whose csq `genes` field names it.
#     qseqids are semicolon-collapsed when a site sits in several regions.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
  library(tidyr)
})

coalesce_null <- function(x, y) if (is.null(x)) y else x

extract_attr_first <- function(attr, pattern) {
  out <- rep(NA_character_, length(attr))
  ok <- !is.na(attr)
  if (!any(ok)) return(out)
  m <- regmatches(attr[ok], regexec(pattern, attr[ok], perl = TRUE))
  out[ok] <- vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  out
}

extract_gff_gene_id <- function(attr) {
  out <- gsub(".*gene_id=([^;]+).*", "\\1", attr)
  out <- ifelse(is.na(out) | out == attr, gsub(".*Name=([^;]+).*", "\\1", attr), out)
  out <- ifelse(is.na(out) | out == attr, gsub(".*ID=([^;]+).*", "\\1", attr), out)
  out
}

extract_gff_attr_optional <- function(attr) {
  gene_name <- extract_attr_first(attr, "^(?:.*;)?(?:Name|gene_name)=([^;]+)")
  product <- extract_attr_first(attr, "^(?:.*;)?(?:product|description|Note)=([^;]+)")
  go_raw <- extract_attr_first(attr, "^(?:.*;)?Ontology_term=([^;]+)")
  go_from_dbxref <- extract_attr_first(attr, "^(?:.*;)?Dbxref=([^;]+)")
  go_terms <- vapply(seq_along(attr), function(i) {
    a <- c(go_raw[i], go_from_dbxref[i])
    a <- a[!is.na(a) & nchar(a) > 0]
    if (length(a) == 0) return(NA_character_)
    gos <- unique(unlist(regmatches(a, gregexpr("GO:[0-9]+", a))))
    if (length(gos) == 0) return(NA_character_)
    paste(gos, collapse = ";")
  }, character(1))
  list(gene_name = gene_name, product = product, go_terms = go_terms)
}

read_assembly_genes <- function(gff_path, verbose = FALSE) {
  if (is.null(gff_path) || !file.exists(gff_path)) stop("GFF not found: ", gff_path)
  d <- tryCatch({
    read_tsv(gff_path, comment = "#",
             col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
             show_col_types = FALSE, na = ".")
  }, error = function(e) stop("Could not read GFF ", gff_path, ": ", conditionMessage(e)))
  if (nrow(d) == 0) stop("GFF is empty: ", gff_path)
  d <- d %>% filter(.data$type %in% c("gene", "mRNA", "CDS"))
  if (nrow(d) == 0) stop("No gene/mRNA/CDS features in ", gff_path)
  d$start <- as.numeric(d$start)
  d$end <- as.numeric(d$end)
  d <- d %>% filter(!is.na(.data$start), !is.na(.data$end))
  d$gene_id <- extract_gff_gene_id(d$attr)
  d$gene_id <- ifelse(is.na(d$gene_id) | !nzchar(trimws(d$gene_id)),
                      paste0(d$seqid, ":", d$start, "-", d$end), d$gene_id)
  opt <- extract_gff_attr_optional(d$attr)
  d$gene_name <- opt$gene_name
  d$product <- opt$product
  d$go_terms <- opt$go_terms
  id_attr <- extract_attr_first(d$attr, "^(?:.*;)?ID=([^;]+)")
  d$csq_key <- dplyr::coalesce(extract_attr_first(d$attr, "^(?:.*;)?Name=([^;]+)"),
                               sub("^gene:", "", id_attr))
  g <- if (any(d$type == "gene")) d %>% filter(.data$type == "gene") else d
  first_non_na <- function(x) {
    x <- x[!is.na(x) & !(is.character(x) & !nzchar(x))]
    if (length(x) == 0) return(NA)
    x[[1]]
  }
  out <- g %>%
    group_by(.data$seqid, .data$gene_id) %>%
    summarise(
      gene_start = min(.data$start), gene_end = max(.data$end),
      csq_key = first_non_na(.data$csq_key),
      gene_name = first_non_na(.data$gene_name),
      product = first_non_na(.data$product),
      go_terms = first_non_na(.data$go_terms),
      strand = first_non_na(.data$strand),
      .groups = "drop"
    ) %>%
    mutate(chr = .data$seqid, gene_length = .data$gene_end - .data$gene_start + 1L)
  if (verbose) message("GFF genes: ", nrow(out), " from ", gff_path)
  out
}

read_table_auto <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  first_line <- readLines(path, n = 1)
  if (length(first_line) == 0) return(NULL)
  delim <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
  if (delim == "\t") read_tsv(path, show_col_types = FALSE) else read_csv(path, show_col_types = FALSE)
}

find_col <- function(d, candidates) {
  nms <- names(d)
  low <- tolower(nms)
  for (cnd in candidates) {
    hit <- which(low == tolower(cnd))
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

rename_known <- function(d, mapping) {
  for (i in seq_along(mapping)) {
    src <- find_col(d, mapping[[i]])
    dest <- names(mapping)[i]
    if (!is.na(src) && src != dest && !dest %in% names(d)) names(d)[names(d) == src] <- dest
  }
  d
}

split_semicolon <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(character(0))
  v <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  v[nzchar(v)]
}

join_unique <- function(x) {
  v <- sort(unique(split_semicolon(x)))
  if (length(v) == 0) NA_character_ else paste(v, collapse = ";")
}

effect_counts_string <- function(x) {
  x <- x[!is.na(x) & x != "none" & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  tb <- sort(table(x), decreasing = TRUE)
  paste0(names(tb), ":", as.integer(tb), collapse = ";")
}

safe_min <- function(x) {
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

safe_sum <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else sum(x)
}

fmt_bp <- function(x) format(x, scientific = FALSE, trim = TRUE)

read_blast_config <- function(path) {
  if (!file.exists(path)) stop("BLAST config file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yml", "yaml")) {
    if (!requireNamespace("yaml", quietly = TRUE))
      stop("Package 'yaml' required for YAML config. Install with: install.packages(\"yaml\")")
    cfg <- yaml::read_yaml(path)
  } else if (ext == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE))
      stop("Package 'jsonlite' required for JSON config. Install with: install.packages(\"jsonlite\")")
    cfg <- jsonlite::read_json(path)
  } else {
    stop("BLAST config must be YAML (.yml/.yaml) or JSON (.json)")
  }
  dbs <- coalesce_null(cfg$databases, cfg)
  if (!is.list(dbs)) dbs <- list(dbs)
  if (length(dbs) > 0 && !is.null(names(dbs)) && !"db_path" %in% names(dbs)) dbs <- list(dbs)
  if (length(dbs) == 0) stop("BLAST config has no databases")
  dbs
}

db_field <- function(db, ...) {
  for (nm in list(...)) {
    if (!is.null(db[[nm]]) && !is.na(db[[nm]]) && nzchar(trimws(as.character(db[[nm]]))))
      return(trimws(as.character(db[[nm]])))
  }
  NULL
}

safe_db_name <- function(name, i) {
  s <- gsub("[^A-Za-z0-9]+", "_", name)
  s <- gsub("^_+|_+$", "", s)
  if (!nzchar(s)) paste0("db", i) else s
}

normalize_sseqid <- function(x) sub("^[^|]+\\|([^|]+)\\|?$", "\\1", x)

# One FASTA record per gene. Headers are g1, g2, ... so BLAST qseqid matches a lookup.
write_gene_fasta <- function(genes, reference, out_fasta, samtools) {
  fai <- paste0(reference, ".fai")
  if (!file.exists(fai)) stop("Reference index not found: ", fai, " (run samtools faidx)")
  reg_vec <- paste0(genes$chr, ":", fmt_bp(genes$gene_start), "-", fmt_bp(genes$gene_end))
  seq <- system2(samtools, c("faidx", reference, reg_vec), stdout = TRUE, stderr = FALSE)
  status <- attr(seq, "status")
  if (!is.null(status) && status != 0L) stop("samtools faidx failed for gene sequences")
  hdr <- which(startsWith(seq, ">"))
  if (length(hdr) != nrow(genes))
    stop("samtools faidx returned ", length(hdr), " sequences for ", nrow(genes), " genes")
  seq[hdr] <- paste0(">", genes$qseqid)
  writeLines(seq, out_fasta)
  out_fasta
}

# NCBI warns that -max_target_seqs below 5 can drop the best hit during the search.
# The summary still keeps one HSP per gene. A requested value above 5 is left as given.
search_target_seqs <- function(max_target_seqs) max(as.integer(max_target_seqs), 5L)

run_gene_blast <- function(query_fasta, db_path, out_tsv, blast_cmd, threads, evalue, max_target_seqs, max_hsps, blast_task,
                           outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue") {
  extra <- c("-evalue", as.character(evalue),
             "-max_target_seqs", as.character(search_target_seqs(max_target_seqs)),
             "-max_hsps", as.character(max_hsps))
  if (!is.null(blast_task) && nzchar(blast_task)) extra <- c("-task", shQuote(blast_task), extra)
  cmd <- c(shQuote(blast_cmd),
           "-query", shQuote(query_fasta),
           "-db", shQuote(db_path),
           extra,
           "-outfmt", shQuote(outfmt),
           "-num_threads", as.character(as.integer(threads)),
           "-out", shQuote(out_tsv))
  script <- tempfile(pattern = "gene_blast_", fileext = ".sh")
  on.exit(unlink(script, force = TRUE), add = TRUE)
  writeLines(paste(cmd, collapse = " "), script)
  result <- system2("sh", script, stdout = TRUE, stderr = TRUE)
  exit <- attr(result, "status")
  if (is.null(exit)) exit <- 0L
  if (exit != 0L) {
    message("BLAST failed (", blast_cmd, " exit ", exit, "): ", paste(result, collapse = " "))
    return(FALSE)
  }
  TRUE
}

fill_blank <- function(x, y) {
  y <- if (is.null(y) || length(y) == 0) NA_character_ else as.character(y[[1]])
  if (is.na(y) || !nzchar(y)) return(x)
  if (is.na(x) || !nzchar(as.character(x))) y else x
}

# Best HSP per query, then the subject GFF gene with the longest overlap.
top_hit_annotations <- function(blast_tsv, qmap, subject_genes, anno) {
  empty <- tibble(gene_id = character(), gene = character(), product = character(), go = character(), evalue = numeric())
  if (!file.exists(blast_tsv) || file.info(blast_tsv)$size == 0) return(empty)
  b <- read_tsv(blast_tsv,
                col_names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue"),
                show_col_types = FALSE)
  if (nrow(b) == 0) return(empty)
  b$evalue_num <- suppressWarnings(as.numeric(b$evalue))
  b$length <- as.numeric(b$length)
  b$hit_start <- pmin(as.numeric(b$sstart), as.numeric(b$send))
  b$hit_end <- pmax(as.numeric(b$sstart), as.numeric(b$send))
  b$sseqid_norm <- normalize_sseqid(as.character(b$sseqid))
  b <- b[order(b$evalue_num, -b$length), , drop = FALSE]
  b <- b[!duplicated(b$qseqid), , drop = FALSE]
  genes_by_seq <- if (is.null(subject_genes) || nrow(subject_genes) == 0) list() else split(subject_genes, subject_genes$chr)
  gene <- product <- go <- rep(NA_character_, nrow(b))
  for (i in seq_len(nrow(b))) {
    g <- genes_by_seq[[b$sseqid_norm[i]]]
    if (is.null(g) || nrow(g) == 0) next
    ov <- pmin(g$gene_end, b$hit_end[i]) - pmax(g$gene_start, b$hit_start[i]) + 1
    if (!any(!is.na(ov) & ov > 0)) next
    pick <- g[which.max(ifelse(is.na(ov) | ov < 1, -Inf, ov)), , drop = FALSE]
    gene[i] <- if (!is.na(pick$gene_name) && nzchar(pick$gene_name)) pick$gene_name else pick$gene_id
    product[i] <- pick$product
    go[i] <- pick$go_terms
    if (!is.null(anno) && nrow(anno) > 0 && "gene_id" %in% names(anno)) {
      hit <- match(pick$gene_id, anno$gene_id)
      if (!is.na(hit)) {
        if ("gene_name" %in% names(anno)) gene[i] <- fill_blank(gene[i], anno$gene_name[hit])
        if ("product" %in% names(anno)) product[i] <- fill_blank(product[i], anno$product[hit])
        go_col <- intersect(c("go_terms", "go"), names(anno))
        if (length(go_col)) go[i] <- fill_blank(go[i], anno[[go_col[1]]][hit])
      }
    }
  }
  asm_id <- qmap$gene_id[match(b$qseqid, qmap$qseqid)]
  tibble(gene_id = asm_id, gene = gene, product = product, go = go, evalue = as.character(b$evalue))
}

blast_outlier_genes <- function(genes, reference, config_path, out_dir, samtools, blast_cmd, threads,
                                evalue, max_target_seqs, max_hsps, blast_task, parallel_dbs, verbose) {
  genes <- genes[!is.na(genes$chr) & !is.na(genes$gene_start) & !is.na(genes$gene_end), , drop = FALSE]
  genes <- genes[!duplicated(genes$gene_id), , drop = FALSE]
  if (nrow(genes) == 0) return(NULL)
  genes$qseqid <- paste0("g", seq_len(nrow(genes)))
  qmap <- genes[, c("qseqid", "gene_id")]
  query_fa <- file.path(out_dir, "genes_query.fa")
  write_gene_fasta(genes, reference, query_fa, samtools)
  if (verbose) message("BLAST query: ", nrow(genes), " genes in ", query_fa)
  targets <- search_target_seqs(max_target_seqs)
  if (targets > as.integer(max_target_seqs))
    message("Searching with -max_target_seqs ", targets,
            " (values below 5 can drop the best hit). One hit per gene is kept.")
  config <- read_blast_config(config_path)
  one_db <- function(i) {
    db <- config[[i]]
    db_path <- db_field(db, "db_path")
    if (is.null(db_path)) {
      message("Skipping BLAST entry ", i, ": db_path missing")
      return(NULL)
    }
    name <- safe_db_name(coalesce_null(db_field(db, "name"), paste0("db", i)), i)
    if (name == "homolog") name <- "homolog_db"
    gff_path <- db_field(db, "gff_path")
    anno_path <- db_field(db, "annotation_tsv", "annotation_table")
    out_tsv <- file.path(out_dir, paste0("blast_genes_", name, ".tsv"))
    message("BLAST ", name, ": ", nrow(genes), " genes against ", db_path,
            " (max_target_seqs=", max_target_seqs, ", max_hsps=", max_hsps, ", evalue=", evalue, ")")
    ok <- run_gene_blast(query_fa, db_path, out_tsv, blast_cmd, threads, evalue, max_target_seqs, max_hsps, blast_task)
    if (!ok) return(NULL)
    subject <- if (!is.null(gff_path) && file.exists(gff_path)) {
      tryCatch(read_assembly_genes(gff_path, verbose = verbose), error = function(e) {
        message("Warning: could not read ", gff_path, ": ", conditionMessage(e))
        NULL
      })
    } else {
      if (!is.null(gff_path)) message("Warning: GFF not found, names left blank: ", gff_path)
      NULL
    }
    anno <- NULL
    if (!is.null(anno_path)) {
      if (file.exists(anno_path)) {
        anno <- tryCatch(read_table_auto(anno_path) %>% rename_all(tolower), error = function(e) {
          message("Warning: could not read ", anno_path, ": ", conditionMessage(e))
          NULL
        })
      } else {
        message("Warning: annotation table not found: ", anno_path)
      }
    }
    hits <- top_hit_annotations(out_tsv, qmap, subject, anno)
    out <- tibble(gene_id = genes$gene_id)
    if (nrow(hits) > 0) out <- out %>% left_join(hits %>% select("gene_id", "gene", "product", "go", "evalue"), by = "gene_id")
    if (!"gene" %in% names(out)) out$gene <- NA_character_
    if (!"product" %in% names(out)) out$product <- NA_character_
    if (!"go" %in% names(out)) out$go <- NA_character_
    if (!"evalue" %in% names(out)) out$evalue <- NA_character_
    renamed <- tibble(gene_id = out$gene_id)
    renamed[[paste0(name, "_gene")]] <- out$gene
    renamed[[paste0(name, "_product")]] <- out$product
    renamed[[paste0(name, "_go")]] <- out$go
    renamed[[paste0(name, "_evalue")]] <- out$evalue
    renamed
  }
  n_parallel <- max(1L, min(as.integer(parallel_dbs), length(config)))
  parts <- if (n_parallel > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_along(config), one_db, mc.cores = n_parallel)
  } else {
    lapply(seq_along(config), one_db)
  }
  parts <- parts[vapply(parts, function(p) is.data.frame(p) && nrow(p) > 0, logical(1))]
  if (length(parts) == 0) return(NULL)
  hom <- parts[[1]]
  if (length(parts) > 1L) {
    for (p in parts[-1]) hom <- hom %>% full_join(p, by = "gene_id")
  }
  gene_cols <- grep("_gene$", names(hom), value = TRUE)
  product_cols <- grep("_product$", names(hom), value = TRUE)
  go_cols <- grep("_go$", names(hom), value = TRUE)
  # The collapsed names themselves must not be picked up by the patterns above.
  row_set <- function(cols) {
    if (length(cols) == 0) return(rep(NA_character_, nrow(hom)))
    apply(hom[, cols, drop = FALSE], 1, join_unique)
  }
  hom$homolog_genes <- row_set(gene_cols)
  hom$homolog_products <- row_set(product_cols)
  hom$homolog_go <- row_set(go_cols)
  hom
}

revcomp_dna <- function(s) {
  if (is.na(s) || !nzchar(s)) return(s)
  chartr("ACGTNacgtn", "TGCANtgcan", paste(rev(strsplit(s, "", fixed = TRUE)[[1]]), collapse = ""))
}

# Spliced CDS for each gene, longest transcript. Headers are g1, g2, ...
write_cds_fasta <- function(genes, gff_path, reference, samtools, out_fasta) {
  d <- read_tsv(gff_path, comment = "#",
                col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
                show_col_types = FALSE, na = ".")
  d <- d %>% filter(.data$type %in% c("gene", "mRNA", "CDS"))
  d$start <- as.numeric(d$start)
  d$end <- as.numeric(d$end)
  d <- d %>% filter(!is.na(.data$start), !is.na(.data$end))
  if (!any(d$type == "CDS")) return(NULL)
  feat_id <- extract_attr_first(d$attr, "^(?:.*;)?ID=([^;]+)")
  feat_name <- extract_attr_first(d$attr, "^(?:.*;)?Name=([^;]+)")
  feat_gene <- extract_gff_gene_id(d$attr)
  feat_parent <- extract_attr_first(d$attr, "^(?:.*;)?Parent=([^;]+)")
  key_to_gene <- character(0)
  add_key <- function(key, gene) {
    ok <- !is.na(key) & nzchar(key) & !is.na(gene) & nzchar(gene)
    if (!any(ok)) return()
    key_to_gene[key[ok]] <<- gene[ok]
  }
  is_gene <- d$type == "gene"
  if (any(is_gene)) {
    add_key(feat_gene[is_gene], feat_gene[is_gene])
    add_key(feat_id[is_gene], feat_gene[is_gene])
    add_key(sub("^gene:", "", feat_id[is_gene]), feat_gene[is_gene])
    add_key(feat_name[is_gene], feat_gene[is_gene])
  }
  is_tx <- d$type == "mRNA"
  tx_gene <- rep(NA_character_, nrow(d))
  if (any(is_tx)) {
    parents <- strsplit(ifelse(is.na(feat_parent), "", feat_parent), ",", fixed = TRUE)
    tx_gene[is_tx] <- vapply(which(is_tx), function(i) {
      hit <- key_to_gene[trimws(parents[[i]])]
      hit <- hit[!is.na(hit) & nzchar(hit)]
      if (length(hit)) hit[[1]] else NA_character_
    }, character(1))
    add_key(feat_id[is_tx], tx_gene[is_tx])
  }
  is_cds <- d$type == "CDS"
  cds <- d[is_cds, , drop = FALSE]
  cds$tx <- feat_parent[is_cds]
  cds$gene_id <- ifelse(!is.na(feat_parent[is_cds]) & feat_parent[is_cds] %in% names(key_to_gene),
                        key_to_gene[feat_parent[is_cds]],
                        ifelse(!is.na(feat_gene[is_cds]) & feat_gene[is_cds] %in% names(key_to_gene),
                               key_to_gene[feat_gene[is_cds]], feat_gene[is_cds]))
  # Parent can list several ids. Use the first one that maps.
  unresolved <- is.na(cds$gene_id) | !nzchar(cds$gene_id)
  if (any(unresolved)) {
    cds$gene_id[unresolved] <- vapply(which(unresolved), function(i) {
      parts <- trimws(strsplit(ifelse(is.na(cds$tx[i]), "", cds$tx[i]), ",", fixed = TRUE)[[1]])
      hit <- key_to_gene[parts]
      hit <- hit[!is.na(hit) & nzchar(hit)]
      if (length(hit)) hit[[1]] else NA_character_
    }, character(1))
  }
  cds <- cds[!is.na(cds$gene_id) & cds$gene_id %in% genes$gene_id, , drop = FALSE]
  if (nrow(cds) == 0) return(NULL)
  cds$tx <- ifelse(is.na(cds$tx) | !nzchar(cds$tx), cds$gene_id, cds$tx)
  cds$width <- cds$end - cds$start + 1
  tx_len <- cds %>% group_by(.data$gene_id, .data$tx) %>% summarise(width = sum(.data$width), .groups = "drop")
  best_tx <- tx_len %>% group_by(.data$gene_id) %>% slice_max(.data$width, n = 1, with_ties = FALSE)
  cds <- cds %>% inner_join(best_tx %>% select("gene_id", "tx"), by = c("gene_id", "tx"))
  cds <- cds[order(match(cds$gene_id, genes$gene_id), cds$start), , drop = FALSE]
  pieces <- lapply(split(cds, cds$gene_id), function(ex) {
    if (any(ex$strand == "-")) ex <- ex[order(-ex$start), , drop = FALSE]
    reg <- paste0(ex$seqid, ":", fmt_bp(ex$start), "-", fmt_bp(ex$end))
    raw <- system2(samtools, c("faidx", reference, reg), stdout = TRUE, stderr = FALSE)
    if (!is.null(attr(raw, "status")) && attr(raw, "status") != 0L) return(NA_character_)
    hdr <- which(startsWith(raw, ">"))
    if (length(hdr) != length(reg)) return(NA_character_)
    seqs <- character(length(hdr))
    for (i in seq_along(hdr)) {
      end <- if (i == length(hdr)) length(raw) else hdr[i + 1] - 1L
      seqs[i] <- paste(raw[(hdr[i] + 1L):end], collapse = "")
    }
    if (any(ex$strand == "-")) seqs <- vapply(seqs, revcomp_dna, character(1))
    paste(seqs, collapse = "")
  })
  keep_ids <- names(pieces)[!is.na(pieces) & nzchar(unlist(pieces, use.names = FALSE))]
  if (length(keep_ids) == 0) return(NULL)
  qmap <- tibble(qseqid = paste0("g", seq_along(keep_ids)), gene_id = keep_ids)
  seqs <- unname(unlist(pieces[keep_ids]))
  lines <- unlist(Map(c, paste0(">", qmap$qseqid), seqs), use.names = FALSE)
  writeLines(lines, out_fasta)
  qmap
}

parse_protein_title <- function(sseqid, stitle) {
  sseqid <- ifelse(is.na(sseqid), "", as.character(sseqid))
  stitle <- ifelse(is.na(stitle), "", as.character(stitle))
  full <- ifelse(nzchar(stitle) & nzchar(sseqid) & !grepl(sseqid, stitle, fixed = TRUE),
                 paste(sseqid, stitle), ifelse(nzchar(stitle), stitle, sseqid))
  gn <- rep(NA_character_, length(full))
  has_sym <- grepl("gene_symbol:", full, fixed = TRUE)
  gn[has_sym] <- sub(".*\\bgene_symbol:(\\S+).*", "\\1", full[has_sym])
  has_gn <- is.na(gn) & grepl("GN=", full, fixed = TRUE)
  gn[has_gn] <- sub(".*\\bGN=(\\S+).*", "\\1", full[has_gn])
  has_gene <- is.na(gn) & grepl("gene:", full, fixed = TRUE)
  gn[has_gene] <- sub(".*\\bgene:(\\S+).*", "\\1", full[has_gene])
  prod <- rep(NA_character_, length(full))
  has_desc <- grepl("description:", full, fixed = TRUE)
  if (any(has_desc)) {
    desc <- sub(".*\\bdescription:", "", full[has_desc])
    prod[has_desc] <- trimws(sub("\\s+\\[Source:.*$", "", desc))
  }
  has_os <- is.na(prod) & grepl("OS=", full, fixed = TRUE)
  if (any(has_os)) {
    tmp <- sub("\\s+OS=.*$", "", full[has_os])
    tmp <- ifelse(grepl("|", tmp, fixed = TRUE), sub("^\\S+\\s+", "", tmp), tmp)
    prod[has_os] <- tmp
  }
  prod[!nzchar(prod)] <- NA_character_
  org <- ifelse(grepl("OS=", full, fixed = TRUE),
                sub(".*\\bOS=(.*?)\\s+(?:OX=|GN=|PE=|SV=).*$", "\\1", full), NA_character_)
  acc <- ifelse(grepl("|", sseqid, fixed = TRUE), sub("^[^|]*\\|([^|]+)\\|.*", "\\1", sseqid), sseqid)
  acc <- sub("\\s.*$", "", acc)
  entry <- sub("^.*\\|", "", sseqid)
  gn <- ifelse(is.na(gn) & nzchar(entry) & grepl("_", entry, fixed = TRUE) & !grepl("\\s", entry),
               sub("_[A-Za-z0-9]+$", "", entry), gn)
  tibble(gene = gn, product = prod, organism = org, accession = acc)
}

present_cols <- function(d, candidates) {
  if (is.null(d) || ncol(d) == 0) return(character(0))
  low <- gsub("[^a-z0-9]+", "_", tolower(names(d)))
  cands <- gsub("[^a-z0-9]+", "_", tolower(candidates))
  names(d)[low %in% cands]
}

first_present <- function(d, candidates) {
  if (is.null(d) || ncol(d) == 0) return(NA_character_)
  low <- gsub("[^a-z0-9]+", "_", tolower(names(d)))
  for (cnd in candidates) {
    hit <- which(low == gsub("[^a-z0-9]+", "_", tolower(cnd)))
    if (length(hit)) return(names(d)[hit[1]])
  }
  NA_character_
}

merge_proteome_rec <- function(old, new) {
  one <- function(a, b) {
    a <- if (is.null(a) || length(a) == 0 || is.na(a[[1]]) || !nzchar(as.character(a[[1]]))) NA_character_ else as.character(a[[1]])
    b <- if (is.null(b) || length(b) == 0 || is.na(b[[1]]) || !nzchar(as.character(b[[1]]))) NA_character_ else as.character(b[[1]])
    if (is.na(a)) b else a
  }
  go <- sort(unique(c(split_semicolon(as.character(old$go)), split_semicolon(as.character(new$go)))))
  list(gene = one(old$gene, new$gene), product = one(old$product, new$product),
       go = if (length(go) == 0) NA_character_ else paste(go, collapse = ";"))
}

put_proteome_key <- function(lookup, key, rec) {
  key <- as.character(key)
  if (is.na(key) || !nzchar(key)) return()
  keys <- unique(c(key, sub("\\.\\d+$", "", key)))
  if (grepl("|", key, fixed = TRUE)) {
    parts <- strsplit(key, "|", fixed = TRUE)[[1]]
    keys <- unique(c(keys, parts, sub("\\.\\d+$", "", parts)))
  }
  for (k in keys) {
    if (is.na(k) || !nzchar(k)) next
    if (exists(k, lookup, inherits = FALSE)) {
      assign(k, merge_proteome_rec(get(k, lookup, inherits = FALSE), rec), envir = lookup)
    } else {
      assign(k, rec, envir = lookup)
    }
  }
}

# Annotation rows win over the FASTA title. id_map connects a BLAST sseqid to that row.
load_proteome_lookup <- function(anno_path, map_path, label) {
  lookup <- new.env(parent = emptyenv())
  anno <- NULL
  if (!is.null(anno_path)) {
    if (!file.exists(anno_path)) {
      message("Warning: proteome annotation table not found: ", anno_path)
    } else {
      anno <- read_table_auto(anno_path)
      if (!is.null(anno)) names(anno) <- tolower(names(anno))
    }
  }
  if (!is.null(anno) && nrow(anno) > 0) {
    id_cols <- present_cols(anno, c("gene_id", "protein_id", "accession", "transcript_id", "seq_id",
                                    "entry", "gene_stable_id", "protein_stable_id", "id"))
    gene_col <- first_present(anno, c("gene_name", "gene_symbol", "symbol", "gene_names", "name"))
    product_col <- first_present(anno, c("product", "protein_name", "protein_names", "description", "function"))
    go_col <- first_present(anno, c("go_terms", "go_term_accession", "go_term_accessions", "go", "go_ids", "go_id", "gene_ontology", "gene_ontology_ids"))
    if (length(id_cols) == 0) {
      message("Warning: ", label, " annotation has no gene/protein id column: ", paste(names(anno), collapse = ", "))
    } else {
      message(label, " annotation: id=", paste(id_cols, collapse = "/"),
              ", gene=", ifelse(is.na(gene_col), "(none)", gene_col),
              ", product=", ifelse(is.na(product_col), "(none)", product_col),
              ", go=", ifelse(is.na(go_col), "(none)", go_col))
      gene <- if (is.na(gene_col)) rep(NA_character_, nrow(anno)) else as.character(anno[[gene_col]])
      product <- if (is.na(product_col)) rep(NA_character_, nrow(anno)) else as.character(anno[[product_col]])
      go <- if (is.na(go_col)) rep(NA_character_, nrow(anno)) else gsub("|", ";", as.character(anno[[go_col]]), fixed = TRUE)
      for (i in seq_len(nrow(anno))) {
        rec <- list(gene = gene[[i]], product = product[[i]], go = go[[i]])
        for (col in id_cols) put_proteome_key(lookup, anno[[col]][[i]], rec)
      }
    }
  }
  if (!is.null(map_path)) {
    if (!file.exists(map_path)) {
      message("Warning: proteome id map not found: ", map_path)
    } else {
      map <- read_table_auto(map_path)
      if (!is.null(map) && nrow(map) > 0) {
        names(map) <- tolower(names(map))
        from_cols <- present_cols(map, c("seq_id", "sseqid", "protein_id", "transcript_id", "accession",
                                         "fasta_id", "entry", "protein_stable_id", "id"))
        to_col <- first_present(map, c("gene_id", "annotation_id", "gene_stable_id", "gene_name", "gene_symbol"))
        from_cols <- setdiff(from_cols, to_col)
        if (length(from_cols) == 0) from_cols <- setdiff(names(map), to_col)[1]
        if (is.na(to_col) || length(from_cols) == 0 || is.na(from_cols[[1]])) {
          message("Warning: ", label, " id map needs a sequence id and a gene id: ", paste(names(map), collapse = ", "))
        } else {
          message(label, " id map: ", paste(from_cols, collapse = "/"), " -> ", to_col)
          for (i in seq_len(nrow(map))) {
            target <- as.character(map[[to_col]][[i]])
            rec <- if (!is.na(target) && nzchar(target) && exists(target, lookup, inherits = FALSE)) {
              get(target, lookup, inherits = FALSE)
            } else {
              list(gene = target, product = NA_character_, go = NA_character_)
            }
            for (col in from_cols) put_proteome_key(lookup, map[[col]][[i]], rec)
          }
        }
      }
    }
  }
  lookup
}

lookup_proteome <- function(lookup, sseqid, accession) {
  keys <- unique(c(sseqid, accession, sub("\\.\\d+$", "", c(sseqid, accession))))
  parts <- unlist(lapply(keys, function(x) if (grepl("|", x, fixed = TRUE)) strsplit(x, "|", fixed = TRUE)[[1]] else character()))
  keys <- unique(c(keys, parts))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  for (k in keys) {
    if (exists(k, lookup, inherits = FALSE)) return(get(k, lookup, inherits = FALSE))
  }
  NULL
}

fill_from_lookup <- function(title_value, rec_value) {
  if (is.null(rec_value) || length(rec_value) == 0 || is.na(rec_value) || !nzchar(as.character(rec_value)))
    return(title_value)
  as.character(rec_value)
}

blast_proteomes <- function(genes, gff_path, reference, config_path, out_dir, samtools, blastx_cmd, threads,
                            evalue, max_target_seqs, max_hsps, parallel_dbs, verbose) {
  genes <- genes[!is.na(genes$gene_id), , drop = FALSE]
  genes <- genes[!duplicated(genes$gene_id), , drop = FALSE]
  if (nrow(genes) == 0) return(NULL)
  config <- read_blast_config(config_path)
  query_fa <- file.path(out_dir, "genes_cds_query.fa")
  qmap <- write_cds_fasta(genes, gff_path, reference, samtools, query_fa)
  if (is.null(qmap) || nrow(qmap) == 0)
    message("Proteome blastx: no CDS features for the outlier genes")
  else if (verbose)
    message("Proteome query: ", nrow(qmap), " spliced CDS in ", query_fa,
            " (", nrow(genes) - nrow(qmap), " genes had no CDS)")
  targets <- search_target_seqs(max_target_seqs)
  if (targets > as.integer(max_target_seqs))
    message("Proteome search uses -max_target_seqs ", targets, ". One hit per proteome is kept.")
  one_db <- function(i) {
    db <- config[[i]]
    db_path <- db_field(db, "db_path")
    name <- safe_db_name(coalesce_null(db_field(db, "name"), paste0("proteome", i)), i)
    if (name == "proteome") name <- "proteome_db"
    empty <- tibble(gene_id = genes$gene_id)
    empty[[paste0(name, "_gene")]] <- NA_character_
    empty[[paste0(name, "_product")]] <- NA_character_
    empty[[paste0(name, "_go")]] <- NA_character_
    empty[[paste0(name, "_accession")]] <- NA_character_
    empty[[paste0(name, "_evalue")]] <- NA_character_
    if (is.null(db_path)) {
      message("Skipping proteome ", name, ": db_path missing")
      return(empty)
    }
    db_files <- paste0(db_path, c(".phr", ".pin", ".psq", ".pdb", ".pal"))
    if (!any(file.exists(db_files))) {
      message("Skipping proteome ", name, ": BLAST database not found at ", db_path)
      return(empty)
    }
    if (is.null(qmap) || nrow(qmap) == 0) return(empty)
    out_tsv <- file.path(out_dir, paste0("blast_genes_", name, ".tsv"))
    message("blastx ", nrow(qmap), " CDS against ", name, " (", db_path, ", evalue=", evalue, ")")
    ok <- run_gene_blast(query_fa, db_path, out_tsv, blastx_cmd, threads, evalue, max_target_seqs, max_hsps,
                         blast_task = NULL,
                         outfmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle")
    if (!ok || !file.exists(out_tsv) || file.info(out_tsv)$size == 0) return(empty)
    b <- read_tsv(out_tsv,
                  col_names = c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle"),
                  show_col_types = FALSE)
    if (nrow(b) == 0) return(empty)
    b$evalue_num <- suppressWarnings(as.numeric(b$evalue))
    b$bitscore <- suppressWarnings(as.numeric(b$bitscore))
    b <- b[order(b$evalue_num, -b$bitscore), , drop = FALSE]
    b <- b[!duplicated(b$qseqid), , drop = FALSE]
    parsed <- parse_protein_title(b$sseqid, b$stitle)
    lookup <- load_proteome_lookup(db_field(db, "annotation_tsv", "annotation_table"),
                                   db_field(db, "id_map", "id_map_tsv"), name)
    gene <- product <- go <- accession <- character(nrow(b))
    for (j in seq_len(nrow(b))) {
      rec <- lookup_proteome(lookup, b$sseqid[[j]], parsed$accession[[j]])
      gene[[j]] <- fill_from_lookup(parsed$gene[[j]], if (is.null(rec)) NA_character_ else rec$gene)
      product[[j]] <- fill_from_lookup(parsed$product[[j]], if (is.null(rec)) NA_character_ else rec$product)
      go[[j]] <- fill_from_lookup(NA_character_, if (is.null(rec)) NA_character_ else rec$go)
      accession[[j]] <- parsed$accession[[j]]
    }
    out <- empty
    hit_row <- match(qmap$gene_id[match(b$qseqid, qmap$qseqid)], out$gene_id)
    ok_hit <- !is.na(hit_row)
    out[[paste0(name, "_gene")]][hit_row[ok_hit]] <- gene[ok_hit]
    out[[paste0(name, "_product")]][hit_row[ok_hit]] <- product[ok_hit]
    out[[paste0(name, "_go")]][hit_row[ok_hit]] <- go[ok_hit]
    out[[paste0(name, "_accession")]][hit_row[ok_hit]] <- accession[ok_hit]
    out[[paste0(name, "_evalue")]][hit_row[ok_hit]] <- as.character(b$evalue)[ok_hit]
    out
  }
  n_parallel <- max(1L, min(as.integer(parallel_dbs), length(config)))
  parts <- if (n_parallel > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(seq_along(config), one_db, mc.cores = n_parallel)
  } else {
    lapply(seq_along(config), one_db)
  }
  parts <- parts[vapply(parts, function(p) is.data.frame(p) && nrow(p) > 0, logical(1))]
  if (length(parts) == 0) return(NULL)
  hom <- parts[[1]]
  if (length(parts) > 1L) {
    for (p in parts[-1]) hom <- hom %>% left_join(p, by = "gene_id")
  }
  gene_cols <- grep("_gene$", names(hom), value = TRUE)
  product_cols <- grep("_product$", names(hom), value = TRUE)
  go_cols <- grep("_go$", names(hom), value = TRUE)
  row_set <- function(cols) {
    if (length(cols) == 0) return(rep(NA_character_, nrow(hom)))
    apply(hom[, cols, drop = FALSE], 1, join_unique)
  }
  hom$proteome_genes <- row_set(gene_cols)
  hom$proteome_products <- row_set(product_cols)
  hom$proteome_go <- row_set(go_cols)
  hom
}

annotation_suffix <- function(path) {
  bn <- tools::file_path_sans_ext(basename(path))
  sub("^outlier_csq_(sites|genes|summary)_", "", bn, perl = TRUE)
}

pick_one <- function(files, kind, suffix, verbose) {
  if (length(files) == 0) return(NULL)
  if (!is.null(suffix) && nzchar(suffix)) {
    files <- files[grepl(suffix, basename(files), fixed = TRUE)]
    if (length(files) == 0) return(NULL)
  }
  if (length(files) > 1L) {
    stop("Several ", kind, " files found. Pass --suffix to choose one:\n  ",
         paste(basename(files), collapse = "\n  "))
  }
  if (verbose) message("Using ", kind, ": ", files)
  files
}

discover_annotation <- function(dir, suffix, verbose) {
  if (is.null(dir) || !dir.exists(dir)) return(list())
  list(
    sites = pick_one(list.files(dir, pattern = "^outlier_csq_sites_.*\\.(csv|tsv)$", full.names = TRUE),
                     "sites", suffix, verbose),
    genes = pick_one(list.files(dir, pattern = "^outlier_csq_genes_.*\\.(csv|tsv)$", full.names = TRUE),
                     "csq genes", suffix, verbose),
    summary = pick_one(list.files(dir, pattern = "^outlier_csq_summary_.*\\.(csv|tsv)$", full.names = TRUE),
                       "summary", suffix, verbose),
    blast = pick_one(
      grep("outlier_genes_summary_", list.files(dir, pattern = "^outlier_genes_.*\\.(csv|tsv)$", full.names = TRUE),
           value = TRUE, invert = TRUE),
      "BLAST genes", suffix, verbose)
  )
}

# Interval overlap: left rows vs right intervals, by chromosome.
overlap_hits <- function(left_chr, left_start, left_end, right) {
  if (length(left_chr) == 0 || nrow(right) == 0)
    return(tibble(left_row = integer(), region_row = integer()))
  hits <- list()
  by_chr <- split(seq_along(left_chr), left_chr)
  for (i in seq_len(nrow(right))) {
    idx <- by_chr[[as.character(right$chr[i])]]
    if (is.null(idx)) next
    keep <- left_end[idx] >= right$region_start[i] & left_start[idx] <= right$region_end[i]
    if (any(keep)) hits[[length(hits) + 1L]] <- tibble(left_row = idx[keep], region_row = i)
  }
  bind_rows(hits)
}

gene_keys_row <- function(gene_id, csq_key, gene_name) {
  v <- unique(c(gene_id, sub("^gene:", "", gene_id), csq_key, gene_name))
  v[!is.na(v) & nzchar(v)]
}

# ---------------------------------------------------------------------------
# gene_report.R sources this file for the functions above. Set
# GENE_COLLATION_FUNCS_ONLY to skip the collation run.
if (!isTRUE(get0("GENE_COLLATION_FUNCS_ONLY", ifnotfound = FALSE))) {
option_list <- list(
  make_option(c("--gff"), type = "character", default = NULL,
              help = "Assembly GFF3 (the one given to variant_csq.sh / outlier_annotation.sh)"),
  make_option(c("--annotation-dir"), type = "character", default = NULL,
              help = "Directory of outlier_annotation.sh output (outlier_csq_*.csv, optional outlier_genes_*.csv)"),
  make_option(c("--suffix"), type = "character", default = NULL,
              help = "Pick annotation files whose names contain this suffix (required when several runs share a directory)"),
  make_option(c("--sites"), type = "character", default = NULL, help = "outlier_csq_sites_*.csv"),
  make_option(c("--genes"), type = "character", default = NULL, help = "outlier_csq_genes_*.csv"),
  make_option(c("--summary"), type = "character", default = NULL, help = "outlier_csq_summary_*.csv"),
  make_option(c("--blast-genes"), type = "character", default = NULL, help = "outlier_genes_*.csv (BLAST hits)"),
  make_option(c("--include-non-outlier"), action = "store_true", default = FALSE,
              help = "Keep GFF genes that overlap no outlier region. With --blast-config, those genes are queried too"),
  make_option(c("--reference"), type = "character", default = NULL,
              help = "Assembly FASTA (indexed). Required with --blast-config or --proteome-config"),
  make_option(c("--blast-config"), type = "character", default = NULL,
              help = "YAML/JSON list of nucleotide databases (db_path, name, gff_path, optional annotation_tsv). Each outlier gene span is searched; the best HSP is named from that genome's GFF"),
  make_option(c("--proteome-config"), type = "character", default = NULL,
              help = "YAML/JSON list of protein databases (db_path, name, optional annotation_tsv, id_map). Spliced CDS are searched with blastx"),
  make_option(c("--blastx-cmd"), type = "character", default = "blastx", help = "blastx program for --proteome-config [default: %default]"),
  make_option(c("--samtools"), type = "character", default = "samtools", help = "samtools for faidx [default: %default]"),
  make_option(c("--blast-cmd"), type = "character", default = "blastn", help = "BLAST program for --blast-config [default: %default]"),
  make_option(c("--blast-task"), type = "character", default = "blastn",
              help = "BLAST -task for nucleotide databases [default: %default]. megablast misses divergent hits; pass megablast to restore it"),
  make_option(c("--threads"), type = "integer", default = 4L, help = "Threads per BLAST process [default: %default]"),
  make_option(c("--parallel-dbs"), type = "integer", default = 1L, help = "How many databases to search at once [default: %default]"),
  make_option(c("--evalue"), type = "double", default = 1e-5, help = "BLAST -evalue [default: %default]"),
  make_option(c("--max-target-seqs"), type = "integer", default = 1L, help = "BLAST -max_target_seqs [default: %default]"),
  make_option(c("--max-hsps"), type = "integer", default = 1L, help = "BLAST -max_hsps [default: %default]"),
  make_option(c("--output-dir"), type = "character", default = ".", help = "Output directory [default: %default]"),
  make_option(c("--verbose"), action = "store_true", default = FALSE)
)
opts <- parse_args(OptionParser(
  usage = "usage: %prog --gff FILE --annotation-dir DIR [options]",
  option_list = option_list
))
if (is.null(opts$gff)) stop("--gff is required")
if (is.null(opts$`annotation-dir`) && is.null(opts$sites) && is.null(opts$genes) && is.null(opts$summary))
  stop("Give --annotation-dir or at least one of --sites, --genes, --summary")
if ((!is.null(opts$`blast-config`) || !is.null(opts$`proteome-config`)) && is.null(opts$reference))
  stop("--reference is required with --blast-config or --proteome-config")
if (!is.null(opts$threads) && opts$threads < 1L) stop("--threads must be >= 1")
if (!is.null(opts$`parallel-dbs`) && opts$`parallel-dbs` < 1L) stop("--parallel-dbs must be >= 1")
if (!is.null(opts$`max-target-seqs`) && opts$`max-target-seqs` < 1L) stop("--max-target-seqs must be >= 1")
if (!is.null(opts$`max-hsps`) && opts$`max-hsps` < 1L) stop("--max-hsps must be >= 1")

found <- discover_annotation(opts$`annotation-dir`, opts$suffix, opts$verbose)
sites_path <- coalesce_null(opts$sites, found$sites)
genes_path <- coalesce_null(opts$genes, found$genes)
summary_path <- coalesce_null(opts$summary, found$summary)
blast_path <- coalesce_null(opts$`blast-genes`, found$blast)
if (is.null(sites_path) && is.null(genes_path) && is.null(summary_path))
  stop("No outlier_csq_sites/genes/summary files found")

suffix <- opts$suffix
if (is.null(suffix) || !nzchar(suffix)) {
  src <- coalesce_null(summary_path, coalesce_null(sites_path, genes_path))
  suffix <- annotation_suffix(src)
}
if (!nzchar(suffix)) suffix <- "genes"

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)
asm <- read_assembly_genes(opts$gff, verbose = opts$verbose)

sites <- read_table_auto(sites_path)
genes_ann <- read_table_auto(genes_path)
summary_ann <- read_table_auto(summary_path)
blast <- read_table_auto(blast_path)

sites <- if (is.null(sites)) sites else rename_known(sites, list(
  qseqid = "qseqid", region_chr = c("region_chr", "chr"),
  region_start = "region_start", region_end = "region_end",
  chr = "chr", pos = "pos"
))
genes_ann <- if (is.null(genes_ann)) genes_ann else rename_known(genes_ann, list(
  qseqid = "qseqid", gene = c("gene", "gene_id", "csq_key"),
  region_chr = c("region_chr", "chr"), region_start = "region_start", region_end = "region_end"
))
summary_ann <- if (is.null(summary_ann)) summary_ann else rename_known(summary_ann, list(
  qseqid = "qseqid", region_chr = c("region_chr", "chr"),
  region_start = "region_start", region_end = "region_end"
))
blast <- if (is.null(blast)) blast else rename_known(blast, list(qseqid = "qseqid"))

# Regions: prefer summary (one row per region, with pass-through stats); else unique from sites/genes.
regions <- NULL
if (!is.null(summary_ann) && nrow(summary_ann) > 0) {
  regions <- summary_ann
} else {
  src <- bind_rows(
    if (!is.null(sites) && nrow(sites) > 0) sites %>% select(any_of(c("qseqid", "region_chr", "region_start", "region_end", "outlier_stat", "n_windows"))) else NULL,
    if (!is.null(genes_ann) && nrow(genes_ann) > 0) genes_ann %>% select(any_of(c("qseqid", "region_chr", "region_start", "region_end", "outlier_stat", "n_windows"))) else NULL
  )
  if (!is.null(src) && nrow(src) > 0) {
    regions <- src %>%
      mutate(qseqid = if ("qseqid" %in% names(.)) .data$qseqid else paste0(.data$region_chr, ":", .data$region_start, "-", .data$region_end)) %>%
      distinct(.data$qseqid, .keep_all = TRUE)
  }
}
if (is.null(regions) || nrow(regions) == 0) stop("No outlier regions in the annotation files")
if (!"region_chr" %in% names(regions) && "chr" %in% names(regions)) regions$region_chr <- regions$chr
if (!all(c("region_chr", "region_start", "region_end") %in% names(regions)))
  stop("Annotation regions need region_chr, region_start, region_end (or qseqid plus those columns)")
regions$region_start <- as.numeric(regions$region_start)
regions$region_end <- as.numeric(regions$region_end)
if (!"qseqid" %in% names(regions)) {
  regions$qseqid <- paste0(regions$region_chr, ":", regions$region_start, "-", regions$region_end)
} else {
  miss <- is.na(regions$qseqid) | !nzchar(as.character(regions$qseqid))
  if (any(miss))
    regions$qseqid[miss] <- paste0(regions$region_chr[miss], ":", regions$region_start[miss], "-", regions$region_end[miss])
}
regions$chr <- regions$region_chr
if (opts$verbose) message("Regions: ", nrow(regions))

# GFF genes overlapping each region.
olap <- overlap_hits(asm$chr, asm$gene_start, asm$gene_end, regions)
gene_regions <- if (nrow(olap) == 0) {
  tibble(gene_id = character(), qseqid = character())
} else {
  bind_cols(
    asm[olap$left_row, c("gene_id", "csq_key", "gene_name", "product", "go_terms", "strand",
                         "chr", "gene_start", "gene_end", "gene_length"), drop = FALSE],
    regions[olap$region_row, , drop = FALSE] %>% select(-any_of(c("chr")))
  ) %>% mutate(overlaps_region = TRUE)
}

# Csq-named genes that missed the GFF overlap (name mismatch): add from outlier_csq_genes.
if (!is.null(genes_ann) && nrow(genes_ann) > 0 && "gene" %in% names(genes_ann)) {
  key_map <- lapply(seq_len(nrow(asm)), function(i) gene_keys_row(asm$gene_id[i], asm$csq_key[i], asm$gene_name[i]))
  names(key_map) <- asm$gene_id
  key_to_id <- unlist(lapply(names(key_map), function(id) setNames(rep(id, length(key_map[[id]])), key_map[[id]])))
  genes_ann$gene_id_mapped <- unname(key_to_id[as.character(genes_ann$gene)])
  extra <- genes_ann %>% filter(is.na(.data$gene_id_mapped))
  matched <- genes_ann %>% filter(!is.na(.data$gene_id_mapped))
  if (nrow(matched) > 0) {
    have <- paste(gene_regions$gene_id, gene_regions$qseqid, sep = "\t")
    add <- matched %>%
      transmute(gene_id = .data$gene_id_mapped, qseqid = .data$qseqid) %>%
      filter(!paste(.data$gene_id, .data$qseqid, sep = "\t") %in% have)
    if (nrow(add) > 0) {
      add_asm <- asm %>% semi_join(add, by = "gene_id")
      add_reg <- regions %>% semi_join(add, by = "qseqid")
      add_rows <- add %>%
        left_join(add_asm, by = "gene_id") %>%
        left_join(add_reg %>% select(-any_of("chr")), by = "qseqid") %>%
        mutate(overlaps_region = TRUE)
      gene_regions <- bind_rows(gene_regions, add_rows)
    }
  }
  if (nrow(extra) > 0 && opts$verbose)
    message("Csq gene names not in GFF (dropped from gene tables): ",
            paste(unique(extra$gene), collapse = ", "))
}

# Sites -> genes by position and/or csq name.
gene_sites <- tibble()
if (!is.null(sites) && nrow(sites) > 0 && "pos" %in% names(sites)) {
  sites$pos <- as.numeric(sites$pos)
  if (!"chr" %in% names(sites) && "region_chr" %in% names(sites)) sites$chr <- sites$region_chr
  site_dummy_end <- sites$pos
  pos_hits <- overlap_hits(sites$chr, sites$pos, site_dummy_end,
                           asm %>% transmute(chr = .data$chr, region_start = .data$gene_start, region_end = .data$gene_end))
  pos_pairs <- if (nrow(pos_hits) == 0) tibble(site_row = integer(), gene_id = character(), assigned_by = character()) else
    tibble(site_row = pos_hits$left_row, gene_id = asm$gene_id[pos_hits$region_row], assigned_by = "position")

  name_pairs <- tibble(site_row = integer(), gene_id = character(), assigned_by = character())
  if ("genes" %in% names(sites)) {
    key_to_id <- unlist(lapply(seq_len(nrow(asm)), function(i) {
      k <- gene_keys_row(asm$gene_id[i], asm$csq_key[i], asm$gene_name[i])
      setNames(rep(asm$gene_id[i], length(k)), k)
    }))
    for (i in seq_len(nrow(sites))) {
      names_i <- split_semicolon(sites$genes[i])
      ids <- unique(unname(key_to_id[names_i]))
      ids <- ids[!is.na(ids)]
      if (length(ids)) name_pairs <- bind_rows(name_pairs, tibble(site_row = i, gene_id = ids, assigned_by = "csq"))
    }
  }
  pairs <- bind_rows(pos_pairs, name_pairs) %>%
    group_by(.data$site_row, .data$gene_id) %>%
    summarise(assigned_by = paste(sort(unique(.data$assigned_by)), collapse = "+"), .groups = "drop")
  if (nrow(pairs) > 0) {
    site_keep <- setdiff(names(sites), c("gene_id", "gene_name", "csq_key", "product", "go_terms",
                                         "strand", "gene_start", "gene_end", "gene_length", "assigned_by"))
    gene_sites <- bind_cols(
      asm[match(pairs$gene_id, asm$gene_id), c("gene_id", "csq_key", "gene_name", "product", "go_terms",
                                               "strand", "chr", "gene_start", "gene_end", "gene_length"), drop = FALSE] %>%
        rename(gene_chr = chr),
      sites[pairs$site_row, site_keep, drop = FALSE],
      tibble(assigned_by = pairs$assigned_by)
    )
    # One row per (gene, site); collapse regions.
    site_key <- intersect(c("chr", "pos", "ref", "alt"), names(gene_sites))
    grp <- c("gene_id", site_key)
    list_cols <- intersect(c("qseqid", "outlier_stat"), names(gene_sites))
    other <- setdiff(names(gene_sites), c(grp, list_cols, "assigned_by",
                                          "region_chr", "region_start", "region_end", "n_windows"))
    gene_sites <- gene_sites %>%
      group_by(across(all_of(grp))) %>%
      summarise(
        qseqids = if ("qseqid" %in% names(gene_sites)) join_unique(.data$qseqid) else NA_character_,
        outlier_stat = if ("outlier_stat" %in% names(gene_sites)) join_unique(.data$outlier_stat) else NA_character_,
        assigned_by = paste(sort(unique(.data$assigned_by)), collapse = "+"),
        across(any_of(other), ~ {
          v <- .[!is.na(.)]
          if (is.character(v)) v <- v[nzchar(v)]
          if (length(v) == 0) NA else v[[1]]
        }),
        .groups = "drop"
      )
  }
}

# Per-(gene, region) site counts from gene_sites qseqids.
if (nrow(gene_sites) > 0 && "qseqids" %in% names(gene_sites) && nrow(gene_regions) > 0) {
  gs_long <- gene_sites %>%
    mutate(qseqid = strsplit(coalesce(.data$qseqids, ""), ";", fixed = TRUE)) %>%
    unnest(cols = "qseqid") %>%
    mutate(qseqid = trimws(.data$qseqid)) %>%
    filter(nzchar(.data$qseqid))
  per <- gs_long %>%
    group_by(.data$gene_id, .data$qseqid) %>%
    summarise(
      n_sites = dplyr::n(),
      n_sites_snp = if ("variant_type" %in% names(gs_long)) sum(.data$variant_type == "snp", na.rm = TRUE) else NA_integer_,
      n_sites_indel = if ("variant_type" %in% names(gs_long)) sum(.data$variant_type == "indel", na.rm = TRUE) else NA_integer_,
      n_protein_altering = if ("genes_protein_altering" %in% names(gs_long))
        sum(!is.na(.data$genes_protein_altering) & nzchar(.data$genes_protein_altering)) else 0L,
      coding_effects = if ("coding_effect" %in% names(gs_long)) effect_counts_string(.data$coding_effect) else NA_character_,
      .groups = "drop"
    )
  gene_regions <- gene_regions %>%
    select(-any_of(c("n_sites", "n_sites_snp", "n_sites_indel", "n_protein_altering", "coding_effects"))) %>%
    left_join(per, by = c("gene_id", "qseqid"))
}
if (nrow(gene_regions) > 0) {
  for (col in c("n_sites", "n_sites_snp", "n_sites_indel", "n_protein_altering")) {
    if (col %in% names(gene_regions)) gene_regions[[col]] <- coalesce(as.integer(gene_regions[[col]]), 0L)
  }
}

# BLAST homologs attach by qseqid (other assembly), not by gene_id.
blast_by_region <- NULL
if (!is.null(blast) && nrow(blast) > 0 && "qseqid" %in% names(blast)) {
  blast_by_region <- blast %>%
    group_by(.data$qseqid) %>%
    summarise(
      n_blast_hits = dplyr::n(),
      blast_dbs = if ("db" %in% names(blast)) join_unique(.data$db) else NA_character_,
      blast_gene_ids = if ("gene_id" %in% names(blast)) join_unique(.data$gene_id) else NA_character_,
      blast_products = if ("product" %in% names(blast)) join_unique(.data$product) else NA_character_,
      best_evalue = if ("evalue" %in% names(blast)) safe_min(as.numeric(.data$evalue)) else NA_real_,
      .groups = "drop"
    )
  gene_regions <- gene_regions %>% left_join(blast_by_region, by = "qseqid")
}

# One row per gene.
if (opts$`include-non-outlier`) {
  gene_ids <- asm$gene_id
} else {
  gene_ids <- unique(c(gene_regions$gene_id, if (nrow(gene_sites) > 0) gene_sites$gene_id else character(0)))
}
gene_summary <- asm %>% filter(.data$gene_id %in% gene_ids)
if (nrow(gene_regions) > 0) {
  reg_agg <- gene_regions %>%
    group_by(.data$gene_id) %>%
    summarise(
      n_regions = dplyr::n_distinct(.data$qseqid),
      qseqids = join_unique(.data$qseqid),
      outlier_stat = if ("outlier_stat" %in% names(gene_regions)) join_unique(.data$outlier_stat) else NA_character_,
      n_windows = if ("n_windows" %in% names(gene_regions)) safe_sum(as.numeric(.data$n_windows)) else NA_real_,
      n_blast_hits = if ("n_blast_hits" %in% names(gene_regions)) safe_sum(.data$n_blast_hits) else NA_real_,
      blast_dbs = if ("blast_dbs" %in% names(gene_regions)) join_unique(.data$blast_dbs) else NA_character_,
      blast_gene_ids = if ("blast_gene_ids" %in% names(gene_regions)) join_unique(.data$blast_gene_ids) else NA_character_,
      blast_products = if ("blast_products" %in% names(gene_regions)) join_unique(.data$blast_products) else NA_character_,
      best_evalue = if ("best_evalue" %in% names(gene_regions)) safe_min(.data$best_evalue) else NA_real_,
      .groups = "drop"
    )
  gene_summary <- gene_summary %>% left_join(reg_agg, by = "gene_id")
}
if (nrow(gene_sites) > 0) {
  site_agg <- gene_sites %>%
    group_by(.data$gene_id) %>%
    summarise(
      n_sites = dplyr::n(),
      n_sites_snp = if ("variant_type" %in% names(gene_sites)) sum(.data$variant_type == "snp", na.rm = TRUE) else NA_integer_,
      n_sites_indel = if ("variant_type" %in% names(gene_sites)) sum(.data$variant_type == "indel", na.rm = TRUE) else NA_integer_,
      n_protein_altering = if ("genes_protein_altering" %in% names(gene_sites))
        sum(!is.na(.data$genes_protein_altering) & nzchar(.data$genes_protein_altering)) else 0L,
      coding_effects = if ("coding_effect" %in% names(gene_sites)) effect_counts_string(.data$coding_effect) else NA_character_,
      .groups = "drop"
    )
  gene_summary <- gene_summary %>% left_join(site_agg, by = "gene_id")
}
if (!"n_regions" %in% names(gene_summary)) gene_summary$n_regions <- 0L
if (!"n_sites" %in% names(gene_summary)) gene_summary$n_sites <- 0L
gene_summary$n_regions <- coalesce(as.integer(gene_summary$n_regions), 0L)
gene_summary$n_sites <- coalesce(as.integer(gene_summary$n_sites), 0L)
gene_summary$in_outlier <- gene_summary$n_regions > 0 | gene_summary$n_sites > 0
if ("n_protein_altering" %in% names(gene_summary))
  gene_summary$n_protein_altering <- coalesce(as.integer(gene_summary$n_protein_altering), 0L)

extra_last <- character(0)
if (!is.null(opts$`blast-config`)) {
  queries <- if (isTRUE(opts$`include-non-outlier`)) gene_summary else gene_summary[gene_summary$in_outlier %in% TRUE, , drop = FALSE]
  hom <- blast_outlier_genes(
    queries, opts$reference, opts$`blast-config`, opts$`output-dir`,
    opts$samtools, opts$`blast-cmd`, opts$threads, opts$evalue,
    opts$`max-target-seqs`, opts$`max-hsps`, opts$`blast-task`, opts$`parallel-dbs`, opts$verbose
  )
  if (!is.null(hom) && nrow(hom) > 0) {
    gene_summary <- gene_summary %>% left_join(hom, by = "gene_id")
    extra_last <- setdiff(names(hom), "gene_id")
  }
}
if (!is.null(opts$`proteome-config`)) {
  queries <- if (isTRUE(opts$`include-non-outlier`)) gene_summary else gene_summary[gene_summary$in_outlier %in% TRUE, , drop = FALSE]
  sw <- blast_proteomes(
    queries, opts$gff, opts$reference, opts$`proteome-config`, opts$`output-dir`,
    opts$samtools, opts$`blastx-cmd`, opts$threads, opts$evalue,
    opts$`max-target-seqs`, opts$`max-hsps`, opts$`parallel-dbs`, opts$verbose
  )
  if (!is.null(sw) && nrow(sw) > 0) {
    gene_summary <- gene_summary %>% left_join(sw, by = "gene_id")
    extra_last <- c(extra_last, setdiff(names(sw), "gene_id"))
  }
}

lead <- c("gene_id", "gene_name", "csq_key", "product", "go_terms", "strand",
          "chr", "gene_start", "gene_end", "gene_length", "in_outlier",
          "n_regions", "n_windows", "qseqids",
          "n_sites", "n_sites_snp", "n_sites_indel", "n_protein_altering",
          "coding_effects", "outlier_stat")
list_last <- c(intersect(c("blast_gene_ids", "blast_products", "blast_dbs"), names(gene_summary)), extra_last)
gene_summary <- gene_summary %>%
  select(-any_of("seqid")) %>%
  select(any_of(lead), setdiff(names(gene_summary), c(lead, list_last, "seqid")), any_of(list_last))
gene_summary <- if ("n_protein_altering" %in% names(gene_summary)) {
  gene_summary %>% arrange(desc(.data$n_protein_altering), desc(.data$n_sites), .data$chr, .data$gene_start)
} else {
  gene_summary %>% arrange(desc(.data$n_sites), .data$chr, .data$gene_start)
}

if (nrow(gene_regions) > 0) {
  gr_lead <- c("gene_id", "gene_name", "csq_key", "chr", "gene_start", "gene_end",
               "qseqid", "region_chr", "region_start", "region_end", "outlier_stat",
               "n_windows", "overlaps_region", "n_sites", "n_protein_altering", "coding_effects")
  gene_regions <- gene_regions %>%
    select(any_of(gr_lead), setdiff(names(gene_regions), gr_lead)) %>%
    arrange(.data$chr, .data$gene_start, .data$region_start)
}

if (nrow(gene_sites) > 0) {
  gs_lead <- c("gene_id", "gene_name", "csq_key", "gene_chr", "gene_start", "gene_end",
               "assigned_by", "chr", "pos", "ref", "alt", "variant_type",
               "region_category", "coding_effect", "genes", "genes_protein_altering",
               "qseqids", "outlier_stat")
  # gene_chr may be named chr from asm join after collapse
  if (!"gene_chr" %in% names(gene_sites) && "chr" %in% names(gene_sites)) {
    # site chr is the variant chromosome; gene interval is gene_start/end on the same chr
  }
  gene_sites <- gene_sites %>%
    select(any_of(gs_lead), setdiff(names(gene_sites), gs_lead)) %>%
    arrange(.data$gene_id, .data$chr, .data$pos)
}

out_summary <- file.path(opts$`output-dir`, paste0("outlier_gene_summary_", suffix, ".csv"))
out_regions <- file.path(opts$`output-dir`, paste0("outlier_gene_regions_", suffix, ".csv"))
out_sites <- file.path(opts$`output-dir`, paste0("outlier_gene_sites_", suffix, ".csv"))
write_csv(gene_summary, out_summary, na = "")
write_csv(gene_regions, out_regions, na = "")
write_csv(gene_sites, out_sites, na = "")
message("Wrote ", out_summary, " (", nrow(gene_summary), " genes)")
message("Wrote ", out_regions, " (", nrow(gene_regions), " gene-region rows)")
message("Wrote ", out_sites, " (", nrow(gene_sites), " gene-site rows)")
message("Done.")
}
