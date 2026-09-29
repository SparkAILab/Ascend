## ======================================================================
## benchmarks/cbl/run_ascend_vs_cbl.R
## ASCEND vs CBL sweeps over n, d_x and d_z (paper Fig. 5).
##
## Reviewer comments: 6 (CI tests vs per-test cost), 5 (empirical scaling), minor 9 (coverage),
##          2 (shared metrics)
## Revision notes: Records for both methods: CI tests (ASCEND split into pair / witness /
##   Markov-blanket tests; CBL counted per l0 lasso verdict), time inside
##   the tests, time per test, conditioning-set sizes, true max in-degree,
##   and the shared eval_ancestral() metrics. Uses the fixed simulator.
## How to run: Rscript benchmarks/cbl/run_ascend_vs_cbl.R    (resumable)
##   BENCH_QUICK=1 ... (smoke test), BENCH_SEEDS, BENCH_TIMEOUT, BENCH_METHODS, BENCH_OUT, Z_SCALE
##   Writes benchmarks/cbl/results/results_ascend_vs_cbl_v2.csv
## Full write-up: docs/REVISION_REPORT.md
## ======================================================================

## ============================================================================
## run_ascend_vs_cbl.R
##
## Benchmark of ASCEND against CBL across three one-dimensional parameter
## sweeps of the simulator, sharing a common default point.
##
## Sweeps (default point: n = 1024, d_x = 5, d_z = 10):
##   n   in {256, 512, 1024, 2048, 4096}
##   d_x in {5, 10, 15, 20}
##   d_z in {10, 20, 30, 40, 50}
##
## Each (method, configuration) cell is run over 5 seeds with a per-cell
## timeout (default 1 hour). Results are appended to the output CSV after
## every cell, so an interrupted run can be resumed without recomputation.
##
## Usage:
##   Rscript benchmarks/cbl/run_ascend_vs_cbl.R
##
## Environment variables:
##   BENCH_QUICK=1     smoke-test mode: n = 256 only, one seed
##   BENCH_SEEDS=N     override the number of seeds per cell
##   BENCH_TIMEOUT=S   override the per-cell timeout, in seconds
## ============================================================================

suppressPackageStartupMessages({
  library(R.utils)
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

## Locate the repository root (the folder holding R/ascend.R), so the
## script runs from the repo root or from any sub-folder.
ROOT <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
ASCEND_PATH <- file.path(ROOT, "R", "ascend.R")
CBL_PATH    <- file.path(ROOT, "R", "cbl.R")
EVAL_PATH   <- file.path(ROOT, "R", "eval_metrics.R")
RESULTS_CSV <- Sys.getenv("BENCH_OUT",
                          file.path(ROOT, "benchmarks", "cbl", "results", "results_ascend_vs_cbl_v2.csv"))
dir.create(dirname(RESULTS_CSV), recursive = TRUE, showWarnings = FALSE)
METHODS     <- strsplit(Sys.getenv("BENCH_METHODS", "ascend,cbl"), ",")[[1]]

## ----------------------------------------------------------------------
## Sourcing
## ----------------------------------------------------------------------
## Both files guard their demonstration/example-run block with
## `if (sys.nframe() == 0)`, which is only true when the file is run as
## the top-level script (e.g. via Rscript), not when it's source()'d from
## another script. A plain source() is therefore safe here.

cat("[setup] sourcing ASCEND...\n"); flush.console()
source(ASCEND_PATH, local = FALSE, chdir = TRUE)

cat("[setup] sourcing CBL...\n"); flush.console()
source(CBL_PATH, local = FALSE, chdir = TRUE)
source(EVAL_PATH, local = FALSE, chdir = TRUE)

## ----------------------------------------------------------------------
## Independence-test accounting
## ----------------------------------------------------------------------
## To compare the two methods on a common scale, we count the number of
## conditional independence verdicts each one extracts from the data:
##
##   ASCEND: one verdict per call to ci_pval().
##   CBL:    each l0(x, y, ...) call fits one regression yielding a
##           length-ncol(x) selection vector, i.e. one independence
##           verdict per feature in x conditional on the rest.
##
## Both functions are wrapped in place so that call counts can be reset
## and read out per run. n_ci_tests is the apples-to-apples count used
## for cross-method comparison; n_l0_calls is CBL-specific auxiliary
## bookkeeping.

.work_counter <- new.env(parent = emptyenv())

reset_counters <- function() {
  .work_counter$ascend_ci_calls <- 0
  .work_counter$cbl_l0_calls    <- 0
  .work_counter$cbl_ci_tests    <- 0
  .work_counter$cbl_l0_time     <- 0
}
reset_counters()

local({
  orig_ci <- ci_pval
  assign("ci_pval",
         function(...) {
           .work_counter$ascend_ci_calls <- .work_counter$ascend_ci_calls + 1
           orig_ci(...)
         },
         envir = globalenv())
})

local({
  orig_l0 <- l0
  assign("l0",
         function(x, y, f, prms) {
           p <- if (is.null(dim(x))) 1L else ncol(x)
           .work_counter$cbl_l0_calls <- .work_counter$cbl_l0_calls + 1
           .work_counter$cbl_ci_tests <- .work_counter$cbl_ci_tests + p
           t0  <- Sys.time()
           out <- orig_l0(x, y, f, prms)
           .work_counter$cbl_l0_time <- .work_counter$cbl_l0_time +
             as.numeric(Sys.time() - t0, units = "secs")
           out
         },
         envir = globalenv())
})

is_ascend <- function(method) method == "ascend"

get_ci_tests <- function(method) {
  if (is_ascend(method)) .work_counter$ascend_ci_calls else .work_counter$cbl_ci_tests
}

get_aux_count <- function(method) {
  if (is_ascend(method)) NA_real_ else .work_counter$cbl_l0_calls
}

## ----------------------------------------------------------------------
## Empty result row (shared by timeout/error paths)
## ----------------------------------------------------------------------

## Metric columns come from eval_metrics.R so they match every other benchmark.
na_metric_row <- function() {
  m <- eval_ancestral(matrix(NA_real_, 2, 2, dimnames = list(c("a", "b"), c("a", "b"))),
                      matrix(0, 2, 2, dimnames = list(c("a", "b"), c("a", "b"))))
  as.data.table(lapply(m, function(v) NA_real_))
}

empty_result_row <- function(method, n, d_x, d_z, seed, status,
                             time_sec = NA_real_, mem_mb = NA_real_,
                             n_ci_tests = NA_real_, n_l0_calls = NA_real_,
                             error_msg = NA_character_) {
  data.table(
    method = method, n = n, d_x = d_x, d_z = d_z, seed = seed,
    status = status,
    time_sec = time_sec, mem_mb = mem_mb,
    n_ci_tests = n_ci_tests, n_l0_calls = n_l0_calls,
    ci_time_sec = NA_real_, time_per_test = NA_real_,
    n_pair_tests = NA_real_, n_witness_tests = NA_real_, n_mb_tests = NA_real_,
    mean_cond_size = NA_real_, max_cond_size = NA_real_,
    mean_pair_cond = NA_real_, max_pair_cond = NA_real_, n_sweeps = NA_real_,
    true_max_indeg = NA_real_, true_n_edges = NA_real_,
    na_metric_row(),
    error_msg = error_msg
  )
}

## ----------------------------------------------------------------------
## Single (method, configuration, seed) run
## ----------------------------------------------------------------------
## Generates one simulated dataset, runs the chosen method inside a
## withTimeout() wrapper, and returns one row of metrics. Timeouts and
## errors are caught and reported with status "timeout" / "error" rather
## than aborting the sweep.

run_one_cell <- function(method, n, d_x, d_z, seed, timeout_sec) {
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  
  sim_obj <- sim_dat(
    n        = n,
    d_z      = d_z,
    d_x      = d_x,
    r2       = 0.5,
    lin_pr   = 1,
    sp       = 0.3,
    p_cross  = 0.15,
    x_effect = 0.9,
    z_scale  = Sys.getenv("Z_SCALE", "1") == "1",   # unit-variance background
    seed     = seed
  )
  amat_true <- truth_from_adj(sim_obj$adj_xx)
  adj0 <- sim_obj$adj_xx; adj0[is.na(adj0)] <- 0
  true_max_indeg <- max(rowSums(adj0)); true_n_edges <- sum(adj0)
  if (is.null(sim_obj$params$lin_pr)) sim_obj$params$lin_pr <- 1
  
  reset_counters()
  gc(reset = TRUE, full = TRUE)
  t0   <- proc.time()[["elapsed"]]
  mem0 <- sum(gc()[, "used"])
  
  est <- tryCatch(
    withTimeout({
      switch(method,
             ascend = ascend(sim_obj, alpha = 0.05, alpha_mb = 0.05,
                             fdr = TRUE, min_votes = 1, verbose = FALSE),
             cbl    = cbl_fn(sim_obj, gamma = 0.5, maxiter = 100, B = 50),
             stop("Unknown method: ", method)
      )
    }, timeout = timeout_sec, onTimeout = "error"),
    error = function(e) structure(NA, error_msg = conditionMessage(e))
  )
  
  t1        <- proc.time()[["elapsed"]]
  mem1      <- sum(gc()[, "used"])
  elapsed   <- t1 - t0
  mem_delta <- max(0, mem1 - mem0)
  
  if (is.null(dim(est))) {
    msg    <- attr(est, "error_msg") %||% ""
    status <- if (grepl("reached elapsed time limit|TimeoutException", msg)) "timeout" else "error"
    return(empty_result_row(
      method, n, d_x, d_z, seed, status,
      time_sec = elapsed, mem_mb = mem_delta,
      n_ci_tests = get_ci_tests(method), n_l0_calls = get_aux_count(method),
      error_msg = msg %||% NA_character_
    ))
  }
  
  # CBL returns m[descendant, ancestor]; everything here is [ancestor, descendant]
  if (method == "cbl") est <- cbl_as_ancestral(est)
  if (is.null(rownames(est))) dimnames(est) <- dimnames(amat_true)
  st   <- attr(est, "stats")
  n_ci <- get_ci_tests(method)
  ci_t <- if (!is.null(st)) st$ci_time_sec else .work_counter$cbl_l0_time
  getst <- function(k) if (!is.null(st) && !is.null(st[[k]])) st[[k]] else NA_real_
  
  data.table(
    method = method, n = n, d_x = d_x, d_z = d_z, seed = seed,
    status = "ok",
    time_sec = elapsed, mem_mb = mem_delta,
    n_ci_tests = n_ci, n_l0_calls = get_aux_count(method),
    ci_time_sec = ci_t, time_per_test = if (n_ci > 0) ci_t / n_ci else NA_real_,
    n_pair_tests = getst("n_pair_tests"), n_witness_tests = getst("n_witness_tests"),
    n_mb_tests = getst("n_mb_tests"),
    mean_cond_size = getst("mean_cond_size"), max_cond_size = getst("max_cond_size"),
    mean_pair_cond = getst("mean_pair_cond"), max_pair_cond = getst("max_pair_cond"),
    n_sweeps = getst("n_sweeps"),
    true_max_indeg = true_max_indeg, true_n_edges = true_n_edges,
    as.data.table(eval_ancestral(est, amat_true)),
    error_msg = NA_character_
  )
}

## ----------------------------------------------------------------------
## Sweep definition
## ----------------------------------------------------------------------

DEFAULT_DX <- 5L
DEFAULT_DZ <- 10L
DEFAULT_N  <- 1024L

quick <- Sys.getenv("BENCH_QUICK", "0") == "1"

SEEDS       <- as.integer(Sys.getenv("BENCH_SEEDS", if (quick) "1" else "5"))
TIMEOUT_SEC <- as.numeric(Sys.getenv("BENCH_TIMEOUT", if (quick) "120" else "3600"))

n_vals  <- if (quick) 256L else c(256L, 512L, 1024L, 2048L, 4096L)
dx_vals <- if (quick) 5L   else c(5L, 10L, 15L, 20L)
dz_vals <- if (quick) 10L  else c(10L, 20L, 30L, 40L, 50L)

sweep_n  <- data.table(n = n_vals,    d_x = DEFAULT_DX, d_z = DEFAULT_DZ)
sweep_dx <- data.table(n = DEFAULT_N, d_x = dx_vals,    d_z = DEFAULT_DZ)
sweep_dz <- data.table(n = DEFAULT_N, d_x = DEFAULT_DX, d_z = dz_vals)
configs  <- unique(rbindlist(list(sweep_n, sweep_dx, sweep_dz)))

plan <- CJ(method  = METHODS,
           cfg_idx = seq_len(nrow(configs)),
           seed    = 100L + seq_len(SEEDS))
plan <- merge(plan, configs[, .(cfg_idx = .I, n, d_x, d_z)], by = "cfg_idx")
setorder(plan, n, d_x, d_z, method, seed)

cat(sprintf("[plan] %d cells total (%d configs x %d methods x %d seeds)\n",
            nrow(plan), nrow(configs), length(METHODS), SEEDS))
cat(sprintf("[plan] per-cell timeout: %.0f s\n", TIMEOUT_SEC))
flush.console()

## ----------------------------------------------------------------------
## Resume support
## ----------------------------------------------------------------------
## A cell is skipped if it already appears in the CSV with status "ok" or
## "timeout". Errored cells are retried, since errors may be transient
## (e.g. a resource limit hit under concurrent load).

already_done <- function() {
  if (!file.exists(RESULTS_CSV)) return(NULL)
  fread(RESULTS_CSV)
}

done <- already_done()
if (!is.null(done) && nrow(done) > 0) {
  cat(sprintf("[resume] %d rows already in %s; skipping those cells.\n",
              nrow(done), RESULTS_CSV))
}

is_done <- function(method, n, d_x, d_z, seed) {
  if (is.null(done) || nrow(done) == 0) return(FALSE)
  any(done$method == method & done$n == n & done$d_x == d_x &
        done$d_z == d_z & done$seed == seed &
        done$status %in% c("ok", "timeout"))
}

## ----------------------------------------------------------------------
## Main loop
## ----------------------------------------------------------------------

header_written <- file.exists(RESULTS_CSV) && file.size(RESULTS_CSV) > 0L
total_cells    <- nrow(plan)
t_sweep0       <- proc.time()[["elapsed"]]

for (k in seq_len(total_cells)) {
  row <- plan[k]
  
  if (is_done(row$method, row$n, row$d_x, row$d_z, row$seed)) {
    cat(sprintf("[%4d/%4d] skip (cached): %-11s n=%-5d d_x=%-3d d_z=%-3d seed=%d\n",
                k, total_cells, row$method, row$n, row$d_x, row$d_z, row$seed))
    flush.console()
    next
  }
  
  cat(sprintf("[%4d/%4d] run : %-11s n=%-5d d_x=%-3d d_z=%-3d seed=%d ... ",
              k, total_cells, row$method, row$n, row$d_x, row$d_z, row$seed))
  flush.console()
  
  res <- tryCatch(
    run_one_cell(row$method, row$n, row$d_x, row$d_z, row$seed, TIMEOUT_SEC),
    error = function(e) {
      empty_result_row(row$method, row$n, row$d_x, row$d_z, row$seed,
                       status = "error", error_msg = conditionMessage(e))
    }
  )
  
  fwrite(res, RESULTS_CSV, append = header_written)
  header_written <- TRUE
  
  if (res$status == "ok") {
    cat(sprintf("OK  %7.2fs  dirF1=%.3f  SHD=%g  CI=%g  t/test=%.2gs\n",
                res$time_sec, res$dir_f1 %||% NA_real_,
                res$shd %||% NA_real_, res$n_ci_tests %||% NA_real_,
                res$time_per_test %||% NA_real_))
  } else {
    cat(sprintf("%s after %.1fs  (%s)\n",
                toupper(res$status), res$time_sec %||% NA_real_,
                substr(res$error_msg %||% "", 1, 80)))
  }
  flush.console()
}

t_sweep1 <- proc.time()[["elapsed"]]
cat(sprintf("\n[done] total wall time: %.1f min   results -> %s\n",
            (t_sweep1 - t_sweep0) / 60, RESULTS_CSV))