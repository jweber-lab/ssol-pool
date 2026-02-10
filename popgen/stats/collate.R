#!/usr/bin/env Rscript

###############################################################################
# collate.R
#
# Purpose: Collate diversity, FST, PBE, and optional variant TSV into HDF5.
# Uses (window, step) for grouping; combines multiple inputs per scale with
# sample/pop1,pop2/pop1,pop2,pop3 as columns. Missing or empty FST/PBE input
# dirs are skipped without failing.
# Usage: Rscript collate.R --diversity-dir DIR [--fst-dir DIR] [--pbe-dir DIR]
#        --output-dir DIR [--seq-qual-dir DIR] [--variant-tsv-dir DIR] [options]
# Options: --drop-all-na  Drop rows/sites where every statistic value is NA before ranking/collation.
#
# Output files:
#   diversity_w{W}_s{S}.h5  (windowed only; single is not supported)
#     Group /windows: chr, start, end, sample, pi, theta, tajima_d,
#     *_rank, *_quantile, mean_coverage, mean_mapping_quality, n_snps (total.passed), n_total (total.passed + total.invariant)
#   diversity_*_summary.tsv  Companion summary (mean, median, variance, skew, quantiles)
#
#   fst_w{W}_s{S}.h5  (or fst_single.h5)
#     Group /windows or /sites: chr, start, end, pos, pop1, pop2, fst,
#     fst_rank, fst_quantile, n_snps (optional)
#   fst_*_summary.tsv  Companion summary by pop1, pop2
#
#   pbe_w{W}_s{S}.h5  (or pbe_single.h5)
#     Group /windows or /sites: chr, start, end, pos, pop1, pop2, pop3, pbe,
#     pbe_rank, pbe_quantile (optional)
#   pbe_*_summary.tsv  Companion summary by pop1, pop2, pop3
#
#   variants.h5  Group /sites: chr, pos, start, end (from variant TSV)
#   single_position.h5  (when --single-position-merged) One group /sites: variant-centric
#     table with chr, pos, start, end, variant cols, plus FST/PBE in wide form:
#     Sample1:Sample2.fst, Sample1:Sample2.fst_rank, Sample1:Sample2.fst_quantile;
#     Sample1:Sample2:Sample3.pbe, .pbe_rank, .pbe_quantile. single_position_summary.tsv
#     summarizes all numeric columns.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(readr)
  library(tidyr)
})

if (!requireNamespace("hdf5r", quietly = TRUE)) {
  stop("collate.R requires package 'hdf5r'. Install with: install.packages(\"hdf5r\")")
}

# ---- Helpers ----
detect_delim <- function(path) {
  first <- readLines(path, n = 1L)
  if (grepl("\t", first) && !grepl(",", first)) return("\t")
  if (grepl(",", first)) return(",")
  "\t"
}

read_tab <- function(path, delim = NULL, na = c("", "NA", "nan", "NaN", "NAN")) {
  if (is.null(delim)) delim <- detect_delim(path)
  if (delim == "\t") read_tsv(path, na = na, show_col_types = FALSE)
  else read_csv(path, na = na, show_col_types = FALSE)
}

# Rank and quantile within groups (0-1 quantile: (rank - 0.5) / n to avoid 0/1)
# When n==0 (all NA in group), quantile is set to NA to avoid division by zero.
add_rank_quantile <- function(data, group_vars, value_col) {
  if (length(group_vars) == 0) {
    r <- rank(data[[value_col]], ties.method = "average", na.last = "keep")
    n <- sum(!is.na(data[[value_col]]))
    data[[paste0(value_col, "_rank")]] <- r
    data[[paste0(value_col, "_quantile")]] <- if (n > 0) (r - 0.5) / n else NA_real_
    return(data)
  }
  data %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(
      !!paste0(value_col, "_rank") := rank(.data[[value_col]], ties.method = "average", na.last = "keep"),
      !!paste0(value_col, "_quantile") := {
        n_val <- sum(!is.na(.data[[value_col]]))
        rv <- rank(.data[[value_col]], ties.method = "average", na.last = "keep")
        ifelse(n_val > 0, (rv - 0.5) / n_val, NA_real_)
      }
    ) %>%
    ungroup()
}

# Skewness (moment-based, base R)
skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (is.na(s) || s <= 0) return(NA_real_)
  mean((x - m)^3) / (s^3)
}

# Print HDF5 content preview to stdout (path, groups, primary group columns, first 2 rows) for logging
h5_preview <- function(path, group_name = "windows") {
  if (!file.exists(path)) return(invisible(NULL))
  h5 <- tryCatch(hdf5r::H5File$new(path, mode = "r"), error = function(e) NULL)
  if (is.null(h5)) return(invisible(NULL))
  on.exit(h5$close_all(), add = TRUE)
  root_names <- names(h5)
  message("HDF5 preview: ", path)
  message("  Groups: ", paste(root_names, collapse = ", "))
  if (!group_name %in% root_names) return(invisible(NULL))
  grp <- h5[[group_name]]
  grp_names <- names(grp)
  n <- if (length(grp_names) > 0) length(grp[[grp_names[1L]]][]) else 0
  message("  Group /", group_name, " datasets: ", paste(grp_names, collapse = ", "), " (n=", n, ")")
  if (n == 0) return(invisible(NULL))
  preview <- lapply(grp_names, function(nm) {
    v <- grp[[nm]][]
    if (is.character(v)) head(v, 2) else as.numeric(head(v, 2))
  })
  names(preview) <- grp_names
  preview_df <- as.data.frame(preview, stringsAsFactors = FALSE)
  message("  Preview (first 2 rows):")
  print(preview_df)
  invisible(NULL)
}

# Summary companion: mean, median, variance, skew, q01, q05, q95, q99, n per stat (and optional group)
# If summary_path is NULL and return_df is TRUE, returns the summary tibble; otherwise writes to summary_path.
write_summary_companion <- function(data, summary_path, stat_cols = NULL, group_cols = character(0), quantile_probs = c(0.01, 0.05, 0.95, 0.99), return_df = FALSE) {
  num_cols <- names(data)[vapply(data, is.numeric, logical(1L))]
  num_cols <- num_cols[!grepl("_rank$|_quantile$", num_cols)]
  if (length(num_cols) == 0) return(if (return_df) NULL else invisible(NULL))
  if (!is.null(stat_cols)) num_cols <- intersect(num_cols, stat_cols)
  group_cols <- intersect(group_cols, names(data))
  num_cols <- num_cols[!num_cols %in% group_cols]
  if (length(num_cols) == 0) return(if (return_df) NULL else invisible(NULL))
  long <- data %>%
    pivot_longer(all_of(num_cols), names_to = "stat", values_to = "value") %>%
    filter(!is.na(value))
  if (length(group_cols) > 0) {
    long <- long %>% group_by(stat, across(all_of(group_cols)))
  } else {
    long <- long %>% group_by(stat)
  }
  qnames <- paste0("q", sub("\\.", "", sprintf("%.2f", quantile_probs)))
  sum_exprs <- list(
    mean = quote(mean(value, na.rm = TRUE)),
    median = quote(median(value, na.rm = TRUE)),
    variance = quote(var(value, na.rm = TRUE)),
    skew = quote(skewness(value)),
    n = quote(n())
  )
  for (i in seq_along(quantile_probs)) {
    sum_exprs[[qnames[i]]] <- bquote(quantile(value, probs = .(quantile_probs[i]), na.rm = TRUE))
  }
  out <- long %>%
    summarise(!!!sum_exprs, .groups = "drop")
  if (return_df) return(out)
  if (!is.null(summary_path)) {
    readr::write_tsv(out, summary_path)
    message("Wrote ", summary_path)
  }
  invisible(NULL)
}

# Infer window vs single and (window_size, step_size) from table content when filename is ambiguous
infer_scale_from_window_table <- function(path) {
  d <- tryCatch(read_tab(path), error = function(e) NULL)
  if (is.null(d) || nrow(d) < 2) return(list(scale = "single", window_size = NA_real_, step_size = NA_real_))
  d <- d %>% rename_all(tolower)
  start_col <- intersect(c("start", "begin"), names(d))[1L]
  end_col <- intersect(c("end", "stop"), names(d))[1L]
  if (is.na(start_col) || is.na(end_col)) return(list(scale = "single", window_size = NA_real_, step_size = NA_real_))
  start_v <- as.numeric(d[[start_col]])
  end_v <- as.numeric(d[[end_col]])
  span <- end_v - start_v + 1
  span <- span[!is.na(span)]
  if (length(span) == 0 || all(span <= 1)) return(list(scale = "single", window_size = NA_real_, step_size = NA_real_))
  window_size <- as.numeric(names(which.max(table(round(span)))))
  if (length(window_size) == 0 || is.na(window_size)) return(list(scale = "single", window_size = NA_real_, step_size = NA_real_))
  stride <- diff(sort(unique(start_v)))
  stride <- stride[stride > 0]
  step_size <- if (length(stride) > 0) as.numeric(names(which.max(table(round(stride))))) else window_size
  if (is.na(step_size)) step_size <- window_size
  list(scale = paste0("w", window_size, "_s", step_size), window_size = window_size, step_size = step_size)
}

# ---- Diversity ----
# Parse window and step from filename (e.g. diversity_w1000_s500.tsv -> 1000, 500)
parse_window_step <- function(basename_file) {
  m <- regmatches(basename_file, regexpr("w(\\d+)_s(\\d+)", basename_file, ignore.case = TRUE))
  if (length(m) == 0) return(list(window_size = NA_real_, step_size = NA_real_))
  list(
    window_size = as.numeric(sub("w([0-9]+)_s.*", "\\1", m)),
    step_size = as.numeric(sub("w[0-9]+_s([0-9]+).*", "\\1", m))
  )
}

find_diversity_files <- function(diversity_dir) {
  dirs <- list.dirs(diversity_dir, recursive = FALSE)
  out <- list()
  for (d in dirs) {
    sample_name <- basename(d)
    # Normalize subdir-derived sample name: e.g. diversity_w1000_s500_sampleCheney -> Cheney
    if (grepl("^diversity_w[0-9]+_s[0-9]+_sample(.+)$", sample_name, ignore.case = TRUE)) {
      sample_name <- sub("^diversity_w[0-9]+_s[0-9]+_sample(.+)$", "\\1", sample_name, ignore.case = TRUE)
    } else if (grepl("^diversity_w[0-9]+_s[0-9]+_(.+)$", sample_name, ignore.case = TRUE)) {
      sample_name <- sub("^diversity_w[0-9]+_s[0-9]+_(.+)$", "\\1", sample_name, ignore.case = TRUE)
    }
    files <- list.files(d, pattern = ".*diversity.*\\.(tsv|csv)$", full.names = TRUE, ignore.case = TRUE)
    for (f in files) {
      if (grepl("single", basename(f), ignore.case = TRUE)) {
        message("Diversity: ignoring file (diversity is never single): ", f)
        next
      }
      ps <- parse_window_step(basename(f))
      out[[length(out) + 1L]] <- list(path = f, sample = sample_name, window_size = ps$window_size, step_size = ps$step_size)
    }
  }
  files <- list.files(diversity_dir, pattern = ".*diversity.*\\.(tsv|csv)$", full.names = TRUE, ignore.case = TRUE)
  for (f in files) {
    if (grepl("single", basename(f), ignore.case = TRUE)) {
      message("Diversity: ignoring file (diversity is never single): ", f)
      next
    }
    ps <- parse_window_step(basename(f))
    base <- tools::file_path_sans_ext(basename(f))
    sample_name <- base
    if (grepl("^diversity_w[0-9]+_s[0-9]+_sample(.+)$", base, ignore.case = TRUE)) {
      sample_name <- sub("^diversity_w[0-9]+_s[0-9]+_sample(.+)$", "\\1", base, ignore.case = TRUE)
    } else if (grepl("^diversity_w[0-9]+_s[0-9]+_(.+)$", base, ignore.case = TRUE)) {
      sample_name <- sub("^diversity_w[0-9]+_s[0-9]+_(.+)$", "\\1", base, ignore.case = TRUE)
    } else {
      sample_name <- sub("^([A-Za-z0-9_]+)_diversity.*", "\\1", base, ignore.case = TRUE)
      if (sample_name == base) sample_name <- "unknown"
    }
    out[[length(out) + 1L]] <- list(path = f, sample = sample_name, window_size = ps$window_size, step_size = ps$step_size)
  }
  out
}

read_one_diversity <- function(path, sample_name, window_size, step_size) {
  d <- read_tab(path)
  d <- d %>% rename_all(tolower)
  chr_col <- intersect(c("chr", "chromosome", "chrom"), names(d))[1L]
  pos_col <- intersect(c("position", "pos", "start", "end"), names(d))[1L]
  if (is.na(chr_col) || is.na(pos_col)) return(NULL)
  d <- d %>% rename(chr = !!chr_col, pos = !!pos_col)
  d$chr <- as.character(d$chr)
  d$pos <- as.numeric(d$pos)
  pi_col <- grep("\\.theta_pi$", names(d), value = TRUE)[1L]
  th_col <- grep("\\.theta_watterson$", names(d), value = TRUE)[1L]
  td_col <- grep("\\.tajimas?_d$", names(d), value = TRUE)[1L]
  if (!is.na(pi_col)) d$pi <- as.numeric(d[[pi_col]])
  if (!is.na(th_col)) d$theta <- as.numeric(d[[th_col]])
  if (!is.na(td_col)) d$tajima_d <- as.numeric(d[[td_col]])
  passed_col <- intersect(c("total.passed", "total_passed"), names(d))[1L]
  invariant_col <- intersect(c("total.invariant", "total_invariant"), names(d))[1L]
  if (!is.na(passed_col)) {
    d$n_snps <- as.numeric(d[[passed_col]])
  }
  if (!is.na(passed_col) && !is.na(invariant_col)) {
    d$n_total <- as.numeric(d[[passed_col]]) + as.numeric(d[[invariant_col]])
  }
  # Use sample from CSV if present and non-NA; otherwise use name derived from path/filename
  if (!("sample" %in% names(d) && any(!is.na(d$sample)))) {
    d$sample <- sample_name
  } else {
    d$sample <- as.character(d$sample)
  }
  d$window_size <- window_size
  d$step_size <- step_size
  d %>% select(any_of(c("chr", "pos", "sample", "window_size", "step_size", "pi", "theta", "tajima_d", "n_snps", "n_total")))  # nolint: object_usage_linter
}

collate_diversity <- function(diversity_dir, output_dir, seq_qual_dir = NULL, variant_tsv_dir = NULL, write_summary = TRUE) {
  infos <- find_diversity_files(diversity_dir)
  if (length(infos) == 0) return(invisible(NULL))
  scale_key <- vapply(infos, function(x) {
    if (is.na(x$window_size) && is.na(x$step_size)) "single"
    else paste0("w", x$window_size, "_s", x$step_size)
  }, "")
  by_scale <- split(infos, scale_key)
  if ("single" %in% names(by_scale)) {
    message("Diversity: single-scale files are not supported; ignoring ", length(by_scale[["single"]]), " file(s).")
    by_scale[["single"]] <- NULL
  }
  for (key in names(by_scale)) {
    scale_infos <- by_scale[[key]]
    tbls <- lapply(scale_infos, function(x) read_one_diversity(x$path, x$sample, x$window_size, x$step_size))
    tbls <- tbls[!vapply(tbls, is.null, logical(1L))]
    if (length(tbls) == 0) next
    d <- bind_rows(tbls)
    d <- d %>% filter(!is.na(chr), !is.na(pos))
    if (!"start" %in% names(d)) d$start <- d$pos
    if (!"end" %in% names(d)) d$end <- d$pos
    if (!"n_snps" %in% names(d)) {
      message("Diversity: no n_snps column (e.g. total.passed) in input; consider re-running with grenedalf output that includes these counts.")
    }
    for (stat in c("pi", "theta", "tajima_d")) {
      if (stat %in% names(d)) d <- add_rank_quantile(d, c("sample", "window_size", "step_size"), stat)
    }
    if (!is.null(seq_qual_dir) && dir.exists(seq_qual_dir)) {
      sq_files <- list.files(seq_qual_dir, pattern = "seq_qual_metrics.*\\.tsv$", full.names = TRUE, recursive = TRUE)
      sq_list <- lapply(sq_files, function(f) {
        x <- read_tab(f)
        if (!all(c("chr", "start", "end", "sample", "mean_coverage", "mean_mapping_quality") %in% names(x))) return(NULL)
        x
      })
      sq_list <- sq_list[!vapply(sq_list, is.null, logical(1L))]
      if (length(sq_list) > 0) {
        sq <- bind_rows(sq_list)
        d <- d %>% left_join(sq %>% select(all_of(c("chr", "start", "end", "sample", "mean_coverage", "mean_mapping_quality"))),
                             by = c("chr", "start", "end", "sample"))
      }
    }
    out_path <- file.path(output_dir, paste0("diversity_", key, ".h5"))
    write_diversity_h5(d, out_path)
    message("Wrote ", out_path)
    h5_preview(out_path, "windows")
    if (write_summary) write_summary_companion(d, sub("\\.h5$", "_summary.tsv", out_path), group_cols = c("sample", "window_size", "step_size"))
  }
  invisible(NULL)
}

write_diversity_h5 <- function(d, path) {
  h5 <- hdf5r::H5File$new(path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  grp <- h5$create_group("windows")
  grp[["chr"]] <- as.character(d$chr)
  grp[["start"]] <- as.numeric(d$start)
  grp[["end"]] <- as.numeric(d$end)
  grp[["sample"]] <- as.character(d$sample)
  for (col in c("pi", "theta", "tajima_d", "pi_rank", "pi_quantile", "theta_rank", "theta_quantile", "tajima_d_rank", "tajima_d_quantile", "mean_coverage", "mean_mapping_quality", "n_snps", "n_total")) {
    if (col %in% names(d)) grp[[col]] <- as.numeric(d[[col]])
  }
  invisible()
}

# ---- FST ----
# FST SNP count comes from total.passed. Total (passed+invariant) or total.invariant could be
# considered in future for consistency with diversity (see help/grenedalf.wiki/Output.md).
find_fst_files <- function(fst_dir) {
  files <- list.files(fst_dir, pattern = ".*fst.*\\.(tsv|csv)$", full.names = TRUE, ignore.case = TRUE)
  out <- list()
  for (f in files) {
    b <- basename(f)
    is_single <- grepl("_single\\.|_single$", tools::file_path_sans_ext(b), ignore.case = TRUE)
    if (is_single) {
      out[[length(out) + 1L]] <- list(path = f, window_size = NA_real_, step_size = NA_real_, scale = "single")
    } else {
      ps <- parse_window_step(b)
      if (is.na(ps$window_size) || is.na(ps$step_size)) {
        inf <- infer_scale_from_window_table(f)
        out[[length(out) + 1L]] <- list(path = f, window_size = inf$window_size, step_size = inf$step_size, scale = inf$scale)
      } else {
        scale <- paste0("w", ps$window_size, "_s", ps$step_size)
        out[[length(out) + 1L]] <- list(path = f, window_size = ps$window_size, step_size = ps$step_size, scale = scale)
      }
    }
  }
  out
}

process_one_fst <- function(path, window_size = NULL, step_size = NULL, drop_all_na = FALSE) {
  d <- read_tab(path)
  d <- d %>% rename_all(tolower)
  chr_col <- intersect(c("chr", "chromosome", "chrom"), names(d))[1L]
  pos_col <- intersect(c("position", "pos", "start", "end"), names(d))[1L]
  if (is.na(chr_col) || is.na(pos_col)) return(NULL)
  d <- d %>% rename(chr = !!chr_col, pos = !!pos_col)
  d$chr <- as.character(d$chr)
  d$pos <- as.numeric(d$pos)
  if (!"start" %in% names(d)) d$start <- d$pos
  if (!"end" %in% names(d)) d$end <- d$pos
  if (!is.null(window_size)) d$window_size <- as.numeric(window_size)
  if (!is.null(step_size)) d$step_size <- as.numeric(step_size)
  n_snps_col <- intersect(c("total.passed", "total_passed"), names(d))[1L]
  if (!is.na(n_snps_col)) d$n_snps <- as.numeric(d[[n_snps_col]])
  fst_cols <- grep("\\.fst$|_fst$", names(d), value = TRUE, ignore.case = TRUE)
  for (fc in fst_cols) {
    d[[fc]] <- as.numeric(d[[fc]])
    pair_name <- gsub("\\.fst$|_fst$", "", fc, ignore.case = TRUE)
    names(d)[names(d) == fc] <- paste0("fst_", pair_name)
  }
  if (drop_all_na) {
    fst_stat_cols <- grep("^fst_.+", names(d), value = TRUE)
    fst_stat_cols <- fst_stat_cols[!grepl("_rank$|_quantile$", fst_stat_cols)]
    if (length(fst_stat_cols) > 0) {
      d <- d[rowSums(!is.na(d[, fst_stat_cols, drop = FALSE])) > 0, , drop = FALSE]
    }
  }
  d
}

# Reshape wide FST table (one column per pair) to long format with pop1, pop2, fst; rank/quantile added later by (pop1, pop2, window_size, step_size)
fst_wide_to_long <- function(d) {
  fst_cols <- grep("^fst_.+", names(d), value = TRUE)
  fst_cols <- fst_cols[!grepl("_rank$|_quantile$", fst_cols)]
  if (length(fst_cols) == 0) return(d)
  base_cols <- c("chr", "start", "end", "pos", "n_snps", "window_size", "step_size")
  base_cols <- intersect(base_cols, names(d))
  out_list <- list()
  for (fc in fst_cols) {
    pair_name <- sub("^fst_", "", fc)
    parts <- strsplit(pair_name, ":")[[1]]
    if (length(parts) < 2) parts <- strsplit(pair_name, "_")[[1]]
    pop1 <- if (length(parts) >= 1) parts[1] else pair_name
    pop2 <- if (length(parts) >= 2) parts[2] else ""
    row <- d[, base_cols, drop = FALSE]
    row$pop1 <- pop1
    row$pop2 <- pop2
    row$fst <- as.numeric(d[[fc]])
    out_list[[length(out_list) + 1L]] <- row
  }
  bind_rows(out_list)
}

write_fst_h5 <- function(d, out_path, group_name = "windows") {
  h5 <- hdf5r::H5File$new(out_path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  grp <- h5$create_group(group_name)
  for (col in names(d)) grp[[col]] <- if (is.character(d[[col]])) as.character(d[[col]]) else as.numeric(d[[col]])
  h5$close_all()
  on.exit()
  message("Wrote ", out_path)
}

collate_fst <- function(fst_dir, output_dir, single_position_merged = FALSE, write_summary = TRUE, single_position_data = NULL, variant_tsv_dir = NULL, drop_all_na = FALSE) {
  infos <- find_fst_files(fst_dir)
  if (length(infos) == 0) return(invisible(NULL))
  by_scale <- split(infos, vapply(infos, function(x) x$scale, ""))
  for (scale in names(by_scale)) {
    scale_infos <- by_scale[[scale]]
    long_list <- lapply(scale_infos, function(x) {
      d <- process_one_fst(x$path, x$window_size, x$step_size, drop_all_na = drop_all_na)
      if (is.null(d)) return(NULL)
      fst_wide_to_long(d)
    })
    long_list <- long_list[!vapply(long_list, is.null, logical(1L))]
    if (length(long_list) == 0) next
    d_long <- bind_rows(long_list)
    if ("fst" %in% names(d_long)) {
      group_vars <- c("pop1", "pop2", "window_size", "step_size")
      group_vars <- intersect(group_vars, names(d_long))
      d_long <- add_rank_quantile(d_long, group_vars, "fst")
    }
    if (!"n_snps" %in% names(d_long)) {
      message("FST: no n_snps column (e.g. total.passed) in input; consider re-running with grenedalf output that includes these counts.")
    }
    is_single <- scale == "single"
    if (is_single && single_position_merged && !is.null(single_position_data)) {
      single_position_data$fst <- d_long
      next
    }
    out_path <- file.path(output_dir, if (is_single) "fst_single.h5" else paste0("fst_", scale, ".h5"))
    grp_name <- if (is_single) "sites" else "windows"
    write_fst_h5(d_long, out_path, grp_name)
    h5_preview(out_path, grp_name)
    if (write_summary) write_summary_companion(d_long, sub("\\.h5$", "_summary.tsv", out_path), group_cols = c("pop1", "pop2"))
  }
  invisible(NULL)
}

# ---- PBE ----
# PBE filenames: w1000_s500_Echo_Kjer_Cheney_pbe.tsv or single_Echo_Kjer_Cheney_pbe.tsv or pbe_Echo_Kjer_Cheney.tsv
find_pbe_files <- function(pbe_dir) {
  files <- list.files(pbe_dir, pattern = ".*pbe\\.(tsv|csv)$", full.names = TRUE, ignore.case = TRUE)
  out <- list()
  for (i in seq_along(files)) {
    f <- files[i]
    base <- tools::file_path_sans_ext(basename(f))
    base <- sub("_pbe$", "", base, ignore.case = TRUE)
    from_filename_single <- grepl("^single[_\\-]|_single$", base, ignore.case = TRUE)
    ps <- parse_window_step(paste0(basename(f)))
    has_window_in_name <- !is.na(ps$window_size) && !is.na(ps$step_size)
    if (from_filename_single) {
      window_size <- NA_real_
      step_size <- NA_real_
      scale <- "single"
    } else if (has_window_in_name) {
      window_size <- ps$window_size
      step_size <- ps$step_size
      scale <- paste0("w", window_size, "_s", step_size)
    } else {
      inf <- infer_scale_from_window_table(f)
      scale <- inf$scale
      window_size <- inf$window_size
      step_size <- inf$step_size
    }
    is_single <- scale == "single"
    trio_name <- base
    if (grepl("^w[0-9]+_s[0-9]+[_\\-](.+)$", trio_name, ignore.case = TRUE)) {
      trio_name <- sub("^w[0-9]+_s[0-9]+[_\\-](.+)$", "\\1", trio_name, ignore.case = TRUE)
    } else if (grepl("^single[_\\-](.+)$", trio_name, ignore.case = TRUE)) {
      trio_name <- sub("^single[_\\-](.+)$", "\\1", trio_name, ignore.case = TRUE)
    } else if (grepl("^pbe[_\\-](.+)$", trio_name, ignore.case = TRUE)) {
      trio_name <- sub("^pbe[_\\-](.+)$", "\\1", trio_name, ignore.case = TRUE)
    } else if (nchar(trio_name) == 0 || !grepl("[A-Za-z]", trio_name)) {
      trio_name <- paste0("trio_", i)
    }
    out[[length(out) + 1L]] <- list(path = f, window_size = window_size, step_size = step_size, scale = scale, trio_name = trio_name, is_single = is_single)
  }
  out
}

# Parse trio name (e.g. "Echo_Kjer_Cheney" or "Echo-Kjer-Cheney") into pop1, pop2, pop3
parse_trio_pops <- function(trio_name) {
  parts <- strsplit(trio_name, "_|-")[[1]]
  if (length(parts) >= 3) {
    list(pop1 = parts[1], pop2 = parts[2], pop3 = parts[3])
  } else {
    list(pop1 = trio_name, pop2 = NA_character_, pop3 = NA_character_)
  }
}

read_one_pbe <- function(path, trio_name, drop_all_na = FALSE) {
  d <- read_tab(path)
  d <- d %>% rename_all(tolower)
  chr_col <- intersect(c("chr", "chromosome"), names(d))[1L]
  pos_col <- intersect(c("pos", "position", "start"), names(d))[1L]
  if (is.na(chr_col) || is.na(pos_col)) return(NULL)
  d <- d %>% rename(chr = !!chr_col, pos = !!pos_col)
  d$chr <- as.character(d$chr)
  d$pos <- as.numeric(d$pos)
  if (!"start" %in% names(d)) d$start <- d$pos
  if (!"end" %in% names(d)) d$end <- d$pos
  if (drop_all_na && "pbe" %in% names(d)) {
    d <- d %>% filter(!is.na(.data$pbe))
  }
  if ("pbe" %in% names(d)) d <- add_rank_quantile(d, character(0), "pbe")
  pops <- parse_trio_pops(trio_name)
  d$pop1 <- pops$pop1
  d$pop2 <- pops$pop2
  d$pop3 <- pops$pop3
  d$trio_name <- trio_name
  d
}

write_pbe_groups_h5 <- function(tbl_list, out_path, group_prefix = "windows_trio", sites = FALSE) {
  inner_name <- if (sites) "sites_trio" else "windows_trio"
  h5 <- hdf5r::H5File$new(out_path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  for (x in tbl_list) {
    d <- x$data
    trio_safe <- gsub("[^A-Za-z0-9_]", "_", substr(x$trio_name, 1L, 50L))
    grp_name <- paste0(inner_name, "_", trio_safe)
    grp <- h5$create_group(grp_name)
    for (col in names(d)) grp[[col]] <- if (is.character(d[[col]])) as.character(d[[col]]) else as.numeric(d[[col]])
  }
  h5$close_all()
  on.exit()
  message("Wrote ", out_path)
}

# Write PBE in long format: single /windows group with chr, start, end, pop1, pop2, pop3, pbe, pbe_rank, pbe_quantile
write_pbe_long_h5 <- function(tbl_list, out_path, group_name = "windows", sites = FALSE) {
  if (length(tbl_list) == 0) return(invisible(NULL))
  all_d <- bind_rows(lapply(tbl_list, function(x) x$data))
  keep_cols <- c("chr", "start", "end", "pos", "pop1", "pop2", "pop3", "pbe")
  if ("pbe_rank" %in% names(all_d)) keep_cols <- c(keep_cols, "pbe_rank")
  if ("pbe_quantile" %in% names(all_d)) keep_cols <- c(keep_cols, "pbe_quantile")
  keep_cols <- intersect(keep_cols, names(all_d))
  all_d <- all_d[, keep_cols]
  h5 <- hdf5r::H5File$new(out_path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  grp <- h5$create_group(group_name)
  for (col in names(all_d)) grp[[col]] <- if (is.character(all_d[[col]])) as.character(all_d[[col]]) else as.numeric(all_d[[col]])
  h5$close_all()
  on.exit()
  message("Wrote ", out_path)
}

collate_pbe <- function(pbe_dir, output_dir, single_position_merged = FALSE, write_summary = TRUE, single_position_data = NULL, drop_all_na = FALSE) {
  infos <- find_pbe_files(pbe_dir)
  if (length(infos) == 0) return(invisible(NULL))
  by_scale <- split(infos, vapply(infos, function(x) x$scale, ""))
  for (scale in names(by_scale)) {
    scale_infos <- by_scale[[scale]]
    tbl_list <- lapply(scale_infos, function(x) {
      d <- read_one_pbe(x$path, x$trio_name, drop_all_na = drop_all_na)
      if (is.null(d)) return(NULL)
      list(trio_name = x$trio_name, data = d)
    })
    tbl_list <- tbl_list[!vapply(tbl_list, is.null, logical(1L))]
    if (length(tbl_list) == 0) next
    is_single <- scale == "single"
    if (is_single && single_position_merged && !is.null(single_position_data)) {
      single_position_data$pbe <- tbl_list
      next
    }
    out_path <- file.path(output_dir, if (is_single) "pbe_single.h5" else paste0("pbe_", scale, ".h5"))
    grp_name <- if (is_single) "sites" else "windows"
    write_pbe_long_h5(tbl_list, out_path, grp_name, sites = is_single)
    h5_preview(out_path, grp_name)
    if (write_summary && length(tbl_list) > 0) {
      all_d <- bind_rows(lapply(tbl_list, function(x) x$data))
      write_summary_companion(all_d, sub("\\.h5$", "_summary.tsv", out_path), group_cols = c("pop1", "pop2", "pop3"))
    }
  }
  invisible(NULL)
}

# ---- Variant TSV (single-position from BCF export) ----
find_variant_tsv <- function(variant_dir) {
  list.files(variant_dir, pattern = ".*(variant|sites|calls).*\\.(tsv|csv)$", full.names = TRUE, ignore.case = TRUE)
}

# Load variant positions (chr, pos) from variant TSV directory for per-window SNP count fallback
load_variant_positions <- function(variant_tsv_dir) {
  files <- find_variant_tsv(variant_tsv_dir)
  if (length(files) == 0) return(NULL)
  tbls <- lapply(files, function(f) {
    x <- read_tab(f) %>% rename_all(tolower)
    chr_col <- intersect(c("chr", "chromosome", "chrom", "contig"), names(x))[1L]
    pos_col <- intersect(c("pos", "position", "start"), names(x))[1L]
    if (is.na(chr_col) || is.na(pos_col)) return(NULL)
    x %>% rename(chr = !!chr_col, pos = !!pos_col) %>% select(chr, pos)
  })
  tbls <- tbls[!vapply(tbls, is.null, logical(1L))]
  if (length(tbls) == 0) return(NULL)
  d <- bind_rows(tbls)
  d$chr <- as.character(d$chr)
  d$pos <- as.numeric(d$pos)
  d %>% filter(!is.na(chr), !is.na(pos))
}

# WIP: Assessing SNP counts from the variants TSV is too inaccurate; do not use from diversity or FST collation.
# Add n_snps to window table (chr, start, end) by counting variant positions per window. Left here for possible future use.
add_n_snps_from_variants <- function(window_df, variant_tsv_dir) {
  vars <- load_variant_positions(variant_tsv_dir)
  if (is.null(vars) || nrow(vars) == 0) return(window_df)
  windows <- window_df %>% distinct(.data$chr, .data$start, .data$end)
  n_snps_df <- vars %>%
    inner_join(windows, by = "chr") %>%
    filter(.data$pos >= .data$start, .data$pos <= .data$end) %>%
    count(chr, start, end, name = "n_snps")
  window_df %>% left_join(n_snps_df, by = c("chr", "start", "end"))
}

collate_variant_tsv <- function(variant_tsv_dir, output_dir, single_position_merged = FALSE, write_summary = TRUE, single_position_data = NULL) {
  files <- find_variant_tsv(variant_tsv_dir)
  if (length(files) == 0) return(invisible(NULL))
  tbls <- lapply(files, function(f) {
    d <- read_tab(f)
    d %>% rename_all(tolower)
  })
  d <- bind_rows(tbls)
  if (nrow(d) == 0) return(invisible(NULL))
  chr_col <- intersect(c("chr", "chromosome", "chrom", "contig"), names(d))[1L]
  pos_col <- intersect(c("pos", "position", "start"), names(d))[1L]
  if (is.na(chr_col) || is.na(pos_col)) {
    message("Variant TSV: no chr/pos columns found, skipping")
    return(invisible(NULL))
  }
  d <- d %>% rename(chr = !!chr_col, pos = !!pos_col)
  d$chr <- as.character(d$chr)
  d$pos <- as.numeric(d$pos)
  if (!"start" %in% names(d)) d$start <- d$pos
  if (!"end" %in% names(d)) d$end <- d$pos
  depth_cols <- grep("depth|coverage|^dp$", names(d), ignore.case = TRUE, value = TRUE)
  depth_cols <- depth_cols[!grepl("_rank$|_quantile$", depth_cols)]
  quality_cols <- grep("qual|quality|mapping|mapq", names(d), ignore.case = TRUE, value = TRUE)
  quality_cols <- quality_cols[!grepl("_rank$|_quantile$", quality_cols)]
  for (col in c(depth_cols, quality_cols)) {
    if (col %in% names(d) && is.numeric(d[[col]])) d <- add_rank_quantile(d, character(0), col)
  }
  if (single_position_merged && !is.null(single_position_data)) {
    single_position_data$variants <- d
    return(invisible(NULL))
  }
  out_path <- file.path(output_dir, "variants.h5")
  h5 <- hdf5r::H5File$new(out_path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  grp <- h5$create_group("sites")
  for (col in names(d)) grp[[col]] <- if (is.character(d[[col]])) as.character(d[[col]]) else as.numeric(d[[col]])
  h5$close_all()
  on.exit()
  message("Wrote ", out_path)
  h5_preview(out_path, "sites")
  if (write_summary) write_summary_companion(d, file.path(output_dir, "variants_summary.tsv"), group_cols = character(0))
  invisible(NULL)
}

# Build base table for single_position merged: variants, or distinct positions from fst/pbe
single_position_base <- function(single_position_data) {
  if (!is.null(single_position_data$variants) && nrow(single_position_data$variants) > 0) {
    return(single_position_data$variants)
  }
  if (!is.null(single_position_data$fst) && nrow(single_position_data$fst) > 0) {
    return(single_position_data$fst %>% distinct(.data$chr, .data$pos, .data$start, .data$end))
  }
  if (!is.null(single_position_data$pbe) && length(single_position_data$pbe) > 0 && nrow(single_position_data$pbe[[1L]]$data) > 0) {
    return(single_position_data$pbe[[1L]]$data %>% distinct(.data$chr, .data$pos, .data$start, .data$end))
  }
  NULL
}

# Write single_position.h5: one table (group /sites) with variants + FST + PBE in wide format (sample1:sample2.fst, sample1:sample2:sample3.pbe, etc.)
write_single_position_merged <- function(single_position_data, output_dir, write_summary = TRUE) {
  base <- single_position_base(single_position_data)
  if (is.null(base)) return(invisible(NULL))
  if (!"pos" %in% names(base)) base$pos <- base$start
  merged <- base
  if (!is.null(single_position_data$fst) && nrow(single_position_data$fst) > 0) {
    fst <- single_position_data$fst
    if (!"pos" %in% names(fst)) fst$pos <- fst$start
    fst <- fst %>% mutate(pair = paste0(.data$pop1, ":", .data$pop2))
    fst_wide <- fst %>%
      select(.data$chr, .data$pos, .data$pair, .data$fst, .data$fst_rank, .data$fst_quantile) %>%
      pivot_wider(names_from = .data$pair, values_from = c(.data$fst, .data$fst_rank, .data$fst_quantile),
                  names_glue = "{pair}.{.value}")
    merged <- merged %>% left_join(fst_wide, by = c("chr", "pos"))
  }
  if (!is.null(single_position_data$pbe) && length(single_position_data$pbe) > 0) {
    for (x in single_position_data$pbe) {
      pbe_d <- x$data
      if (!"pos" %in% names(pbe_d)) pbe_d$pos <- pbe_d$start
      trio_id <- paste0(pbe_d$pop1[1L], ":", pbe_d$pop2[1L], ":", pbe_d$pop3[1L])
      pbe_wide <- pbe_d %>%
        select(.data$chr, .data$pos, .data$pbe, .data$pbe_rank, .data$pbe_quantile) %>%
        mutate(trio = trio_id) %>%
        pivot_wider(names_from = .data$trio, values_from = c(.data$pbe, .data$pbe_rank, .data$pbe_quantile),
                    names_glue = "{trio}.{.value}")
      merged <- merged %>% left_join(pbe_wide, by = c("chr", "pos"))
    }
  }
  h5_path <- file.path(output_dir, "single_position.h5")
  h5 <- hdf5r::H5File$new(h5_path, mode = "w")
  on.exit(h5$close_all(), add = TRUE)
  grp <- h5$create_group("sites")
  for (col in names(merged)) {
    if (is.character(merged[[col]])) grp[[col]] <- as.character(merged[[col]]) else grp[[col]] <- as.numeric(merged[[col]])
  }
  h5$close_all()
  on.exit()
  message("Wrote ", h5_path)
  h5_preview(h5_path, "sites")
  if (write_summary) {
    s <- write_summary_companion(merged, NULL, group_cols = character(0), return_df = TRUE)
    if (!is.null(s)) {
      readr::write_tsv(s, file.path(output_dir, "single_position_summary.tsv"))
      message("Wrote ", file.path(output_dir, "single_position_summary.tsv"))
    }
  }
  invisible(NULL)
}

# ---- MAIN ----
option_list <- list(
  make_option(c("--diversity-dir"), type = "character", default = NULL, help = "Directory containing diversity TSV/CSV (e.g. output of calculate_pi_theta.sh)"),
  make_option(c("--fst-dir"), type = "character", default = NULL, help = "Directory containing FST TSV/CSV"),
  make_option(c("--pbe-dir"), type = "character", default = NULL, help = "Directory containing PBE TSV/CSV"),
  make_option(c("--seq-qual-dir"), type = "character", default = NULL, help = "Optional: directory with seq_qual_metrics TSV"),
  make_option(c("--variant-tsv-dir"), type = "character", default = NULL, help = "Optional: directory with variant/sites TSV (e.g. from bcftools query export)"),
  make_option(c("--output-dir"), type = "character", default = ".", help = "Output directory for HDF5 files [default: %default]"),
  make_option(c("--single-position-merged"), action = "store_true", default = FALSE, help = "Merge per-locus FST, PBE, and variant TSV into single_position.h5"),
  make_option(c("--no-summary"), action = "store_true", default = FALSE, help = "Do not write companion *_summary.tsv files"),
  make_option(c("--drop-all-na"), action = "store_true", default = FALSE, help = "Drop rows/sites where every statistic value (FST or PBE) is NA before ranking and collation"),
  make_option(c("--verbose"), action = "store_true", default = FALSE, help = "Verbose")
)
parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list)
opts <- parse_args(parser)

if (is.null(opts$`diversity-dir`) && is.null(opts$`fst-dir`) && is.null(opts$`pbe-dir`)) {
  stop("At least one of --diversity-dir, --fst-dir, --pbe-dir is required")
}

dir.create(opts$`output-dir`, showWarnings = FALSE, recursive = TRUE)
write_summary <- !opts$`no-summary`
single_position_merged <- opts$`single-position-merged`
drop_all_na <- opts$`drop-all-na`
single_position_data <- if (single_position_merged) list(fst = NULL, pbe = NULL, variants = NULL) else NULL

if (!is.null(opts$`diversity-dir`) && dir.exists(opts$`diversity-dir`)) {
  collate_diversity(opts$`diversity-dir`, opts$`output-dir`, opts$`seq-qual-dir`, opts$`variant-tsv-dir`, write_summary = write_summary)
}
if (!is.null(opts$`fst-dir`) && dir.exists(opts$`fst-dir`)) {
  collate_fst(opts$`fst-dir`, opts$`output-dir`, single_position_merged = single_position_merged, write_summary = write_summary, single_position_data = single_position_data, variant_tsv_dir = opts$`variant-tsv-dir`, drop_all_na = drop_all_na)
}
if (!is.null(opts$`pbe-dir`) && dir.exists(opts$`pbe-dir`)) {
  collate_pbe(opts$`pbe-dir`, opts$`output-dir`, single_position_merged = single_position_merged, write_summary = write_summary, single_position_data = single_position_data, drop_all_na = drop_all_na)
}
if (!is.null(opts$`variant-tsv-dir`) && dir.exists(opts$`variant-tsv-dir`)) {
  collate_variant_tsv(opts$`variant-tsv-dir`, opts$`output-dir`, single_position_merged = single_position_merged, write_summary = write_summary, single_position_data = single_position_data)
}
if (single_position_merged && !is.null(single_position_data) &&
    (!is.null(single_position_data$fst) || (!is.null(single_position_data$pbe) && length(single_position_data$pbe) > 0) || !is.null(single_position_data$variants))) {
  write_single_position_merged(single_position_data, opts$`output-dir`, write_summary = write_summary)
}

message("Collate done.")
