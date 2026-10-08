## ======================================================================
## analysis/prescreen_scaling.R
## How ASCEND's cost and accuracy scale with the background size d_z,
## and what the background prescreen (prescreen = 0.30) costs in accuracy.
##
## Reviewer comments: 5 (complexity, empirical side), 6 (number of tests
##   vs cost per test), 4 (a strict prescreen behaves like missing
##   background).
## What it does: fixed simulator, d_x = 10, n in {300, 1024},
##   d_z in {25, 50, 100, 200}. ASCEND (guarded conditioning set) is run
##   with the default prescreen = 0.30 and with the prescreen off
##   (prescreen = 1). Records wall time, CI tests by type, time per test,
##   conditioning-set sizes, and the shared accuracy metrics.
## Run: Rscript analysis/prescreen_scaling.R      (~10-15 min; SEEDS=5 default)
## Writes analysis/out/prescreen_scaling.csv
## Full write-up: docs/REVISION_REPORT.md, sections 4 and 7
## ======================================================================
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
suppressPackageStartupMessages(library(R.utils))
for (f in c("ascend.R", "eval_metrics.R")) source(file.path(ROOT, "R", f))
OUT <- file.path(ROOT, "analysis", "out"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEEDS <- 100L + seq_len(as.integer(Sys.getenv("SEEDS", "5")))

rows <- list()
for (n in c(300, 1024)) for (dz in c(25, 50, 100, 200)) for (s in SEEDS) {
  sim <- sim_dat(n = n, d_z = dz, d_x = 10, r2 = 0.5, sp = 0.7, p_cross = 0.1,
                 x_effect = 0.9, z_scale = TRUE, seed = s)
  T <- truth_from_adj(sim$adj_xx)
  for (pre in c(0.30, 1)) {
    t0 <- Sys.time()
    M <- tryCatch(withTimeout(ascend(sim, prescreen = pre, verbose = FALSE), timeout = 600, onTimeout = "error"),
                  error = function(e) NULL)
    el <- as.numeric(Sys.time() - t0, units = "secs")
    base <- data.frame(n = n, d_z = dz, seed = s, prescreen = pre, time_sec = el)
    if (is.null(M)) { rows[[length(rows) + 1]] <- base; next }
    st <- attr(M, "stats"); e <- eval_ancestral(M, T)
    rows[[length(rows) + 1]] <- cbind(base,
      n_ci_tests = st$n_ci_tests, n_pair_tests = st$n_pair_tests, n_witness_tests = st$n_witness_tests,
      n_mb_tests = st$n_mb_tests, ci_time_sec = st$ci_time_sec, time_per_test = st$ci_time_sec / st$n_ci_tests,
      mean_cond_size = st$mean_cond_size, mean_pair_cond = st$mean_pair_cond, max_pair_cond = st$max_pair_cond,
      as.data.frame(e[c("dir_precision", "dir_recall", "dir_f1", "orient_acc", "coverage", "shd")]))
    cat(sprintf("n=%4d d_z=%3d seed=%d prescreen=%.2f  %.1fs  F1=%.2f\n", n, dz, s, pre, el, e$dir_f1))
  }
}
res <- do.call(dplyr::bind_rows, rows)
write.csv(res, file.path(OUT, "prescreen_scaling.csv"), row.names = FALSE)
options(width = 200)
print(aggregate(cbind(time_sec, n_ci_tests, time_per_test, mean_pair_cond, max_pair_cond, dir_f1, orient_acc) ~
                  n + d_z + prescreen, res, mean, na.action = na.pass), digits = 3)
