## ======================================================================
## analysis/revision_tables.R
## Builds the revision tables and figures from the benchmark CSVs.
##
## Reviewer comments: 2, 9, minor 5 (GRN tables with CIs), minor 9 (coverage table), minor 6 (runtime vs edges), 5 and 6 (compute cost)
## Revision notes: Uses re-run results when present, otherwise the submitted results.
## How to run: Rscript analysis/revision_tables.R  -> analysis/out/
## Full write-up: docs/REVISION_REPORT.md
## ======================================================================

## ======================================================================
## Revision tables and figures from benchmark CSVs
## ----------------------------------------------------------------------
## Run from anywhere inside the repository:
##   Rscript analysis/revision_tables.R
##
## Reads whichever result files exist, preferring the re-runs made with
## the revision code (benchmarks/*/results/) and falling back to the
## submitted runs (benchmarks/*/results_submitted/). Writes to
## analysis/out/:
##   grn_summary.csv          mean +/- SE of every metric, every method and
##                            cell (AUROC and AUPR for all cells; SE on
##                            precision)                        [R2 c2, c9]
##   grn_paired.csv           mean paired difference, 95% bootstrap CI,
##                            two-sided Wilcoxon, BH within one family per
##                            metric                            [R2 c9, m5]
##   grn_table2_f1.md         Table 2/3 replacement: F1 mean +/- SE and
##                            ASCEND - competitor mean diff [95% CI]
##   coverage_table.csv       coverage (and precision on resolved pairs)
##                            for ASCEND and CBL                 [R2 m9]
##   runtime_vs_edges.pdf     wall-clock runtime against number of claimed
##                            edges, per method                  [R2 m6]
##   compute_cost.csv / .pdf  ASCEND vs CBL: CI tests, time per test and
##                            conditioning-set size             [R2 c5, c6]
## ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
ROOT <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
source(file.path(ROOT, "R", "stats_utils.R"))
B <- function(...) file.path(ROOT, "benchmarks", ...)

OUT <- file.path(ROOT, "analysis", "out")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
first_existing <- function(...) { f <- Filter(file.exists, c(...)); if (length(f)) f[1] else NA_character_ }
say <- function(...) cat(sprintf(...), "\n")

## ---------------------------------------------------------------------
## 1. GRN benchmark (ASCEND vs GENIE3 / ARACNe / WGCNA, ...)
## ---------------------------------------------------------------------
grn_file <- first_existing(B("grn", "results", "benchmark_v3_main_raw.csv"),
                           B("grn", "results_submitted", "benchmark_v2_main_raw.csv"))
if (!is.na(grn_file)) {
  say("[grn] %s", grn_file)
  g <- fread(grn_file)
  cells <- c("n", "sp", "r2")
  mets  <- intersect(c("aupr", "auroc", "aupr_cont", "auroc_cont", "f1", "precision",
                       "recall", "coverage", "orient_acc", "direction_acc", "dir_f1",
                       "shd", "runtime_s"), names(g))
  summ <- summary_table(g, mets, cells)
  fwrite(summ, file.path(OUT, "grn_summary.csv"))

  pmets <- intersect(c("f1", "precision", "recall", "aupr", "auroc", "aupr_cont", "auroc_cont"), names(g))
  pt <- paired_table(g, pmets, cells, ref = "ASCEND", alternative = "two.sided",
                     family_by = "metric")
  fwrite(pt, file.path(OUT, "grn_paired.csv"))

  # markdown: F1 mean +/- SE per method, and paired diff [CI] vs each competitor
  comps <- setdiff(unique(g$method), "ASCEND")
  f1s <- as.data.table(summ)[, .(n, sp, r2, method, txt = fmt_mean_se(f1_mean, f1_se))]
  f1w <- dcast(f1s, r2 + sp + n ~ method, value.var = "txt")
  pd  <- as.data.table(pt)[metric == "f1", .(n, sp, r2, competitor,
                                             txt = fmt_diff_ci(mean_diff, ci_lo, ci_hi, q_bh))]
  pdw <- dcast(pd, r2 + sp + n ~ competitor, value.var = "txt")
  setnames(pdw, comps, paste0("diff vs ", comps))
  tab <- merge(f1w, pdw, by = c("r2", "sp", "n"))
  md  <- c(paste0("| ", paste(names(tab), collapse = " | "), " |"),
           paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"),
           apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")),
           "",
           sprintf("Differences: mean of per-replicate (ASCEND - competitor) F1 with 95%% bootstrap CI; stars: BH q within the family of all %d F1 comparisons (*** < 0.001, ** < 0.01, * < 0.05), two-sided paired Wilcoxon.",
                   pt$family_size[pt$metric == "f1"][1]))
  writeLines(md, file.path(OUT, "grn_table2_f1.md"))
}

## ---------------------------------------------------------------------
## 2. Coverage table, ASCEND vs CBL
## ---------------------------------------------------------------------
cbl_file <- first_existing(B("cbl", "results", "results_ascend_vs_cbl_v2.csv"),
                           B("cbl", "results_submitted", "results_ascend_vs_cbl.csv"))
if (!is.na(cbl_file)) {
  say("[cbl] %s", cbl_file)
  cb <- fread(cbl_file)[status == "ok"]
  pcol <- if ("dir_precision" %in% names(cb)) "dir_precision" else "precision"
  cov_tab <- cb[, .(runs = .N,
                    coverage_mean = mean(coverage, na.rm = TRUE),
                    coverage_se   = sd(coverage, na.rm = TRUE) / sqrt(.N),
                    precision_mean = mean(get(pcol), na.rm = TRUE)),
                by = .(method, n, d_x, d_z)]
  setorder(cov_tab, n, d_x, d_z, method)
  fwrite(cov_tab, file.path(OUT, "coverage_table.csv"))

  ## computational work: CI tests, time per test, conditioning-set size
  ## (only in re-runs, which record time_per_test)
  if ("time_per_test" %in% names(cb)) {
    ab <- cb[, .(runs = .N, time_sec = mean(time_sec), n_ci_tests = mean(n_ci_tests),
                 time_per_test = mean(time_per_test, na.rm = TRUE),
                 mean_cond_size = mean(mean_cond_size, na.rm = TRUE),
                 max_pair_cond = mean(max_pair_cond, na.rm = TRUE),
                 dir_f1 = mean(dir_f1, na.rm = TRUE)),
             by = .(method, n, d_x, d_z)]
    fwrite(ab, file.path(OUT, "compute_cost.csv"))
    abl <- melt(ab[n == 1024 & d_x == 5], id.vars = c("method", "d_z"),
                measure.vars = c("time_sec", "n_ci_tests", "time_per_test", "mean_cond_size"))
    p <- ggplot(abl, aes(d_z, value, colour = method)) + geom_line() + geom_point() +
      facet_wrap(~variable, scales = "free_y") + scale_y_log10() +
      labs(x = expression(d[z]), y = NULL, colour = NULL) + theme_bw()
    ggsave(file.path(OUT, "compute_cost.pdf"), p, width = 8, height = 5)
  }
}

## ---------------------------------------------------------------------
## 3. Runtime against discovered edge count
## ---------------------------------------------------------------------
cd_file <- first_existing(B("causal_grid", "ascend_benchmark_v3_merged.csv"),
                          B("causal_grid", "results_submitted", "ascend_benchmark_merged.csv"))
if (!is.na(cd_file)) {
  say("[causal] %s", cd_file)
  cd <- fread(cd_file)[status == "ok"]
  # claimed edges = number of positive directed claims
  cd[, edges := if ("dir_claims" %in% names(cd)) dir_claims else tp + fp]
  cd <- cd[is.finite(elapsed_sec) & elapsed_sec > 0]
  p <- ggplot(cd, aes(pmax(edges, 0.5), elapsed_sec, colour = method)) +
    geom_point(alpha = 0.25, size = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Discovered (claimed) ancestral edges", y = "Wall-clock time (s)",
         colour = NULL) + theme_bw()
  ggsave(file.path(OUT, "runtime_vs_edges.pdf"), p, width = 6, height = 4.5)
  slopes <- cd[edges > 0, .(loglog_slope = coef(lm(log(elapsed_sec) ~ log(edges)))[2], runs = .N),
               by = method]
  fwrite(slopes, file.path(OUT, "runtime_vs_edges_slopes.csv"))
}

say("[done] outputs in %s", OUT)
