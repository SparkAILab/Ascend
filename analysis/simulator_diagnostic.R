## ======================================================================
## analysis/simulator_diagnostic.R
## Why every re-run uses the fixed simulator, sim_dat(z_scale = TRUE).
##
## Reviewer comments: 5 and 6 indirectly (the submitted scaling claims rest
##   on this simulator), and the whole simulation section.
## What it does: for d_z from 10 to 300 (d_x = 20, n = 2000), simulates data
##   with the submitted simulator (z_scale = FALSE) and the fixed one
##   (z_scale = TRUE) from the same seeds and records
##     - the largest background variance, max var(Z)
##     - the median |partial correlation| of each true X -> X edge given
##       all of Z (the signal ASCEND has to detect)
##     - ASCEND's true positives and directed F1.
##   It also reads the submitted 81-cell grid and reports ASCEND's mean
##   true positives by d_z / d_x (the collapse at d_z = 3 d_x).
## Run: Rscript analysis/simulator_diagnostic.R     (~3-5 min; SEEDS=3 default)
## Writes analysis/out/simulator_diagnostic.csv and
##        analysis/out/submitted_grid_tp_by_dz.csv
## Full write-up: docs/REVISION_REPORT.md, section 2
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
SEEDS <- 100L + seq_len(as.integer(Sys.getenv("SEEDS", "3")))

pcor_edges <- function(sim) {                 # |partial cor| of true X->X edges given all Z
  X <- as.matrix(sim$dat[, grep("^x", names(sim$dat))]); Z <- as.matrix(sim$dat[, grep("^z", names(sim$dat))])
  R <- resid(lm(X ~ Z)); C <- cor(R)
  e <- which(sim$adj_xx == 1, arr.ind = TRUE)    # adj_xx[child, parent]
  if (!nrow(e)) return(NA_real_)
  median(abs(C[e]))
}

rows <- list()
for (dz in c(10, 20, 40, 60, 100, 200, 300)) for (zs in c(FALSE, TRUE)) for (s in SEEDS) {
  sim <- sim_dat(n = 2000, d_z = dz, d_x = 20, r2 = 0.5, sp = 2/3, p_cross = 0.2,
                 x_effect = 0.9, z_scale = zs, seed = s)
  zv <- apply(sim$dat[, grep("^z", names(sim$dat))], 2, var)
  M <- tryCatch(withTimeout(ascend(sim, verbose = FALSE), timeout = 300, onTimeout = "error"),
                error = function(e) NULL)
  e <- if (is.null(M)) NULL else eval_ancestral(M, truth_from_adj(sim$adj_xx))
  rows[[length(rows) + 1]] <- data.frame(
    d_z = dz, simulator = if (zs) "fixed (z_scale = TRUE)" else "submitted (z_scale = FALSE)",
    seed = s, max_var_z = max(zv), median_edge_pcor = pcor_edges(sim),
    dir_tp = if (is.null(e)) NA else e$dir_tp, dir_f1 = if (is.null(e)) NA else e$dir_f1,
    n_true_rel = if (is.null(e)) NA else e$n_true_rel)
  cat(sprintf("d_z=%3d %-28s seed=%d max var(Z)=%9.3g pcor=%.3f TP=%s\n", dz,
              rows[[length(rows)]]$simulator, s, max(zv), rows[[length(rows)]]$median_edge_pcor,
              rows[[length(rows)]]$dir_tp))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT, "simulator_diagnostic.csv"), row.names = FALSE)
print(aggregate(cbind(max_var_z, median_edge_pcor, dir_tp, dir_f1) ~ d_z + simulator, res, mean,
                na.action = na.pass), digits = 3)

## the submitted 81-cell grid: ASCEND true positives by d_z / d_x
g <- read.csv(file.path(ROOT, "benchmarks", "causal_grid", "results_submitted", "ascend_benchmark_merged.csv"))
g <- g[g$method == "ASCEND" & g$status == "ok", ]
g$dz_over_dx <- g$d_z / g$d_x
tp <- aggregate(cbind(tp, f1) ~ dz_over_dx + d_x, g, mean, na.action = na.pass)
write.csv(tp, file.path(OUT, "submitted_grid_tp_by_dz.csv"), row.names = FALSE)
print(tp, digits = 3)
