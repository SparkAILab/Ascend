## ======================================================================
## analysis/scoring_audit.R
## How much the scoring bugs in the submitted benchmark scripts changed
## the reported numbers.
##
## Reviewer comments: 2 (metrics), 1 (fair baselines), 9 (statistics).
## What it does:
##   A. GRN benchmark (paper Tables 1-3). On the submitted simulator and
##      settings (n = 2000, d_z = 20, d_x = 15, all six (sp, R2) cells),
##      scores the same ASCEND output twice: with the submitted functions
##      (*_legacy, which read only i < j) and with the fixed ones. Reports
##      ASCEND F1, the matched edge budget K_match given to the
##      competitors, and "direction accuracy" (legacy) vs orientation
##      accuracy.
##   B. Causal grid (paper Fig. 6). Reads the submitted merged results and
##      counts runs per method and status. GES and PC "errors" come from
##      the cpdag_to_ancestral() crash on undirected CPDAG edges (fixed in
##      benchmarks/causal_grid/hpc_run.R), so their submitted averages only
##      cover the runs whose CPDAG happened to be fully oriented.
## Run: Rscript analysis/scoring_audit.R     (~2 min; SEEDS=10 default)
## Writes analysis/out/scoring_audit_grn.csv, scoring_audit_grid_status.csv
## Full write-up: docs/REVISION_REPORT.md, section 1
## ======================================================================
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
suppressPackageStartupMessages({ library(data.table); library(dplyr) })
OUT <- file.path(ROOT, "analysis", "out"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEEDS <- seq_len(as.integer(Sys.getenv("SEEDS", "10")))

## load the GRN script's functions (everything before its "Run" section);
## competitor packages are not needed for this audit
L <- readLines(file.path(ROOT, "benchmarks", "grn", "run_grn_benchmark.R"))
L <- L[seq_len(grep("^## --- Run ---", L)[1] - 1L)]
L <- sub("PRROC::pr.curve", "pr_curve_ap", L, fixed = TRUE)
pr_curve_ap <- function(scores.class0, scores.class1, curve = FALSE) {     # average precision
  sc <- c(scores.class0, scores.class1); lb <- c(rep(1, length(scores.class0)), rep(0, length(scores.class1)))
  lb <- lb[order(-sc)]; list(auc.integral = sum(cumsum(lb) / seq_along(lb) * lb) / sum(lb))
}
tmp <- tempfile(fileext = ".R"); writeLines(L, tmp)
library <- function(pkg, ...) { p <- as.character(substitute(pkg))
  if (requireNamespace(p, quietly = TRUE)) suppressPackageStartupMessages(base::library(p, character.only = TRUE)) }
sys.source(tmp, envir = globalenv()); rm(library)

rows <- list()
for (sp in c(0.5, 0.7, 0.9)) for (r2 in c(0.5, 0.7)) for (s in SEEDS) {
  sim <- sim_dat(n = 2000, d_z = 20, d_x = 15, r2 = r2, lin_pr = 1, sp = sp,
                 p_cross = 0.2, x_effect = 0.8, z_scale = FALSE, seed = 42 + s)
  gt <- prepare_gt(sim); x <- gt$xlabs
  M <- ascend(sim, verbose = FALSE)
  b_new <- ascend_binary_skeleton(M, x); b_old <- ascend_binary_skeleton_legacy(M, x)
  f_new <- prf_from_binary(b_new, gt$A_skel); f_old <- prf_from_binary(b_old, gt$A_skel)
  if (is.null(rownames(M))) dimnames(M) <- list(x, x)
  rows[[length(rows) + 1]] <- data.frame(
    sp = sp, r2 = r2, seed = s,
    K_true = gt$K, K_match_submitted = sum(b_old[upper.tri(b_old)]), K_match_fixed = sum(b_new[upper.tri(b_new)]),
    f1_submitted = f_old$f1, f1_fixed = f_new$f1,
    precision_submitted = f_old$precision, precision_fixed = f_new$precision,
    recall_submitted = f_old$recall, recall_fixed = f_new$recall,
    direction_acc_submitted = ascend_direction_acc(M, x, gt$A_anc),
    orient_acc_fixed = eval_ancestral(M, gt$A_anc)$orient_acc)
}
A <- do.call(rbind, rows)
write.csv(A, file.path(OUT, "scoring_audit_grn.csv"), row.names = FALSE)
cat("== A. GRN benchmark: ASCEND scored with the submitted vs fixed functions ==\n")
print(aggregate(. ~ sp + r2, A[, setdiff(names(A), "seed")], mean), digits = 3)

g <- fread(file.path(ROOT, "benchmarks", "causal_grid", "results_submitted", "ascend_benchmark_merged.csv"))
B <- dcast(g[, .N, by = .(method, status)], method ~ status, value.var = "N", fill = 0)
B[, error_share := round(error / (error + ok), 3)]
fwrite(B, file.path(OUT, "scoring_audit_grid_status.csv"))
cat("\n== B. Submitted causal grid: runs by method and status ==\n"); print(B)
