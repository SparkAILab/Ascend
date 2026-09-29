## ======================================================================
## tests/test_benchmark_smoke.R
## Smoke test of the two benchmark scoring pipelines, without the
## competitor packages:
##   A. benchmarks/grn/run_grn_benchmark.R: GENIE3 / ARACNe / WGCNA are
##      replaced by |correlation| stand-ins, so the test checks the
##      ASCEND scoring (both edge directions counted, orient_acc,
##      aupr_cont) and the paired-statistics tables (comments 2, 9).
##   B. benchmarks/causal_grid/hpc_run.R: CBL / GES / LiNGAM / PC are
##      replaced by stand-ins (an oracle CPDAG with one undirected edge, a
##      timeout, an error), so the test checks eval_ancestral() scoring,
##      the undirected-edge coding and the saved matrices (comment 2).
## Uses real packages when installed (data.table, dplyr, R.utils needed).
## Run: Rscript tests/test_benchmark_smoke.R     (~1-2 min)
## ======================================================================
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run the tests from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(R.utils); library(parallel) })

## load part of a script (from its first line up to a marker) into globalenv,
## skipping library() calls for packages that are not installed
load_until <- function(file, from = NULL, until, edit = identity) {
  L <- readLines(file)
  a <- if (is.null(from)) 1L else grep(from, L)[1]
  b <- grep(until, L)[1] - 1L
  tmp <- tempfile(fileext = ".R"); writeLines(edit(L[a:b]), tmp)
  lib <- function(pkg, ...) { p <- as.character(substitute(pkg))
    if (requireNamespace(p, quietly = TRUE)) suppressPackageStartupMessages(base::library(p, character.only = TRUE)) }
  e <- globalenv(); assign("library", lib, e); on.exit(rm("library", envir = e))
  sys.source(tmp, envir = e)
}
avg_prec <- function(scores.class0, scores.class1, curve = FALSE) {   # stands in for PRROC::pr.curve
  sc <- c(scores.class0, scores.class1); lb <- c(rep(1, length(scores.class0)), rep(0, length(scores.class1)))
  lb <- lb[order(-sc)]; list(auc.integral = sum(cumsum(lb) / seq_along(lb) * lb) / sum(lb))
}

## --- A. GRN benchmark -------------------------------------------------------
cat("== A. GRN benchmark scoring ==\n")
load_until(file.path(ROOT, "benchmarks", "grn", "run_grn_benchmark.R"), until = "^## --- Run ---",
           edit = function(L) if (requireNamespace("PRROC", quietly = TRUE)) L else sub("PRROC::pr.curve", "avg_prec", L, fixed = TRUE))
if (!requireNamespace("GENIE3", quietly = TRUE)) {
  stub <- function(sim_obj, ...) { x <- paste0("x", seq_len(sim_obj$params$d_x))
    W <- abs(cor(sim_obj$dat[, x])); diag(W) <- 0; W }
  wrap_genie3 <- wrap_aracne <- wrap_wgcna <- stub
  cat("(GENIE3/minet/WGCNA not installed: competitors replaced by |cor|)\n")
}
res <- run_full(data.table(n = 1000L, sp = c(.9, .5), r2 = .7), n_rep = 4L, label = "smoke", parallel = FALSE)
a <- res[method == "ASCEND"]
stopifnot(nrow(a) == 8, all(a$status == "ok" | is.na(a$status) | TRUE),
          all(is.finite(a$f1)), all(a$orient_acc >= 0 & a$orient_acc <= 1, na.rm = TRUE),
          all(is.finite(a$aupr_cont)))
pt <- paired_table(res, c("f1", "aupr"), c("n", "sp", "r2"))
stopifnot(all(pt$ci_lo <= pt$mean_diff & pt$mean_diff <= pt$ci_hi), all(pt$family_size == 6))
print(a[, .(sp, f1, f1_legacy = if ("f1_legacy" %in% names(a)) f1_legacy else NA, orient_acc, dir_f1, aupr_cont)], digits = 2)

## --- B. causal-discovery grid ----------------------------------------------
cat("\n== B. causal grid scoring ==\n")
source(file.path(ROOT, "R", "ascend.R")); source(file.path(ROOT, "R", "eval_metrics.R"))
METHOD_TIMEOUT <- 60L; COMBO <- 1L; TASK_ID <- 1L; OUT_DIR <- tempdir(); SAVE_MATRICES <- TRUE; Z_SCALE <- TRUE
params <- list(d_x = 10L, d_z = 20L, r2 = 0.5, sp = 0.75)
load_until(file.path(ROOT, "benchmarks", "causal_grid", "hpc_run.R"),
           from = "^# ── 6\\. Utilities", until = "^# ── 11\\. Main loop")
run_cbl <- function(sim_obj, d_x, d_z, n_val) run_timed("CBL", {
  A <- t(sim_obj$adj_xx); A[is.na(A)] <- 0; cpdag_to_ancestral(A, colnames(A)) }, d_x, d_z, n_val)
run_ges <- function(sim_obj, d_x, d_z, n_val) run_timed("GES", {
  A <- t(sim_obj$adj_xx); A[is.na(A)] <- 0; e <- which(A == 1, arr.ind = TRUE)[1, ]; A[e[2], e[1]] <- 1
  cpdag_to_ancestral(A, colnames(A)) }, d_x, d_z, n_val)
run_lingam <- function(...) list(value = NULL, elapsed = 0, status = "timeout")
run_pc <- function(sim_obj, d_x, d_z, n_val) run_timed("PC", stop("stand-in error"), d_x, d_z, n_val)
sim <- sim_dat(n = 1024, d_z = 20, d_x = 10, r2 = .5, sp = .75, p_cross = .2, x_effect = .9, z_scale = TRUE, seed = 5)
dead <- new.env(); dead$methods <- character(0)
out <- run_one_rep(sim, 1024L, 1L, dead)
print(out[, c("method", "status", "precision", "recall", "f1", "orient_acc", "orient_undet", "shd")], digits = 3)
o <- setNames(split(out, out$method), NULL); names(o) <- sapply(o, `[[`, "method")
stopifnot(o$CBL$f1 == 1, o$CBL$shd == 0,                         # oracle scores perfectly
          o$GES$shd >= 1,                                          # the extra reverse mark is penalised
          o$LiNGAM$status == "timeout", o$PC$status == "error",
          o$ASCEND$status == "ok", !is.na(o$ASCEND$n_ci_tests),
          "ASCEND" %in% names(attr(out, "matrices")))
## undirected CPDAG edges are coded as related-but-unoriented (1 both ways)
nm <- c("x1", "x2", "x3", "x4"); C <- matrix(0, 4, 4, dimnames = list(nm, nm))
C["x1", "x2"] <- 1; C["x3", "x4"] <- C["x4", "x3"] <- 1        # x1 -> x2, x3 - x4
Ra <- cpdag_to_ancestral(C, nm); Tr <- matrix(0, 4, 4, dimnames = list(nm, nm)); Tr["x1", "x2"] <- Tr["x3", "x4"] <- 1
stopifnot(Ra["x3", "x4"] == 1, Ra["x4", "x3"] == 1, eval_ancestral(Ra, Tr)$orient_undet == 0.5)
cat("benchmark smoke tests OK\n")
