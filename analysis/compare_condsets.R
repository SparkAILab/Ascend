## ======================================================================
## Guarded vs full conditioning sets: a quick, self-contained comparison
## ----------------------------------------------------------------------
## Same algorithm, same Fisher-z test, same data. Only the conditioning
## set differs:
##   guarded : nearest ancestors of both endpoints (the ASCEND design)
##   full    : Z U known non-descendants of both endpoints (CBL's set)
## Reports wall-clock time, CI tests (split by kind), time per test,
## conditioning-set size, and the shared accuracy metrics.
##
## Run from the repo root (base R + R.utils only):
##   Rscript analysis/compare_condsets.R                # default grid, ~5-10 min
##   SEEDS=10 Rscript analysis/compare_condsets.R       # more seeds
## Writes analysis/out/compare_condsets.csv (per run) and prints means.
## ======================================================================

suppressPackageStartupMessages(library(R.utils))
source("ascend.R"); source("eval_metrics.R")

SEEDS   <- 100L + seq_len(as.integer(Sys.getenv("SEEDS", "5")))
TIMEOUT <- as.numeric(Sys.getenv("TIMEOUT", "300"))
# (n, d_z) pairs: the low-n / high-d_z end is where the size of the
# conditioning set matters most for test power
GRID <- data.frame(n   = c(300, 300, 300, 1024, 1024, 1024),
                   d_z = c( 50, 150, 250,   50,  100,  200))
D_X <- 10L

one <- function(sim, truth, cs) {
  t0 <- Sys.time()
  M  <- tryCatch(withTimeout(ascend(sim, cond_set = cs, verbose = FALSE),
                             timeout = TIMEOUT, onTimeout = "error"),
                 error = function(e) NULL)
  el <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(M)) return(data.frame(cond_set = cs, status = "timeout/error", time_s = el))
  s <- attr(M, "stats"); e <- eval_ancestral(M, truth)
  data.frame(cond_set = cs, status = "ok", time_s = el,
             ci_tests = s$n_ci_tests, blanket_tests = s$n_mb_tests,
             us_per_test = 1e6 * s$ci_time_sec / max(1, s$n_ci_tests),
             mean_S = s$mean_cond_size, max_S_ij = s$max_pair_cond,
             dir_precision = e$dir_precision, dir_recall = e$dir_recall,
             dir_f1 = e$dir_f1, orient_acc = e$orient_acc, coverage = e$coverage,
             shd = e$shd)
}

rows <- list()
for (g in seq_len(nrow(GRID))) for (seed in SEEDS) {
  n <- GRID$n[g]; dz <- GRID$d_z[g]
  sim <- sim_dat(n = n, d_z = dz, d_x = D_X, r2 = 0.5, sp = 0.7, p_cross = 0.10,
                 x_effect = 0.9, z_scale = TRUE, seed = seed)
  truth <- truth_from_adj(sim$adj_xx)
  for (cs in c("guarded", "full"))
    rows[[length(rows) + 1]] <- cbind(n = n, d_z = dz, seed = seed, one(sim, truth, cs))
  cat(sprintf("n=%d d_z=%d seed=%d done\n", n, dz, seed))
}
res <- do.call(dplyr::bind_rows, rows)
dir.create("analysis/out", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "analysis/out/compare_condsets.csv", row.names = FALSE)

num <- c("time_s", "ci_tests", "blanket_tests", "us_per_test", "mean_S", "max_S_ij",
         "dir_precision", "dir_recall", "dir_f1", "orient_acc", "coverage", "shd")
ok  <- res[res$status == "ok", ]
agg <- aggregate(ok[, num], by = ok[, c("n", "d_z", "cond_set")],
                 FUN = function(v) mean(v, na.rm = TRUE))
options(width = 200)
print(agg[order(agg$n, agg$d_z, agg$cond_set), ], digits = 3, row.names = FALSE)
