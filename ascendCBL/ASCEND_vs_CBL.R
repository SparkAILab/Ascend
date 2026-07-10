## ============================================================================
## benchmark_ascend_vs_cbl.R
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
##   Rscript benchmark_ascend_vs_cbl.R
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

ASCEND_PATH <- "../ascend.R"
CBL_PATH    <- "../cbl.R"
RESULTS_CSV <- "results_ascend_vs_cbl.csv"

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
           orig_l0(x, y, f, prms)
         },
         envir = globalenv())
})

get_ci_tests <- function(method) {
  if (method == "ascend") .work_counter$ascend_ci_calls else .work_counter$cbl_ci_tests
}

get_aux_count <- function(method) {
  if (method == "ascend") NA_real_ else .work_counter$cbl_l0_calls
}

## ----------------------------------------------------------------------
## Structural Hamming distance and unresolved-pair fraction
## ----------------------------------------------------------------------
## evaluate() (defined in ascend.R) already yields precision,
## recall, F1, coverage, and TP/FP/FN/TN. We add two further metrics:
##
##   shd_full/shd_resolved: SHD over the upper triangle of the estimated
##     ancestral matrix vs. the truth. shd_full penalises unresolved (NA)
##     pairs at true positives; shd_resolved counts only cells where the
##     method made an explicit positive or negative claim.
##   unresolved_frac: fraction of truth-defined pairs left unresolved.

compute_shd <- function(estimated, truth) {
  nms <- rownames(truth)
  est <- estimated[nms, nms, drop = FALSE]
  d   <- nrow(truth)
  
  full <- 0L; resolved <- 0L; n_res <- 0L
  
  for (i in seq_len(d - 1)) {
    for (j in (i + 1):d) {
      t_ij <- truth[i, j]
      if (is.na(t_ij)) next
      
      e_ij <- est[i, j]
      e_ji <- est[j, i]
      
      true_pos <- (t_ij == 1)
      is_na    <- is.na(e_ij) && is.na(e_ji)
      pred_pos <- (!is.na(e_ij) && e_ij %in% c(0.5, 1)) ||
        (!is.na(e_ji) && e_ji %in% c(0.5, 1))
      
      if (is_na) {
        if (true_pos) full <- full + 1L
      } else {
        n_res <- n_res + 1L
        if (true_pos != pred_pos) {
          full     <- full + 1L
          resolved <- resolved + 1L
        }
      }
    }
  }
  
  list(shd_full = full, shd_resolved = resolved, n_resolved = n_res)
}

unresolved_frac <- function(estimated, truth) {
  nms <- rownames(truth)
  est <- estimated[nms, nms, drop = FALSE]
  d   <- nrow(truth)
  
  total <- 0L; na_count <- 0L
  
  for (i in seq_len(d - 1)) {
    for (j in (i + 1):d) {
      if (is.na(truth[i, j])) next
      total <- total + 1L
      if (is.na(est[i, j])) na_count <- na_count + 1L
    }
  }
  
  if (total == 0L) NA_real_ else na_count / total
}

## ----------------------------------------------------------------------
## Empty result row (shared by timeout/error paths)
## ----------------------------------------------------------------------

empty_result_row <- function(method, n, d_x, d_z, seed, status,
                             time_sec = NA_real_, mem_mb = NA_real_,
                             n_ci_tests = NA_real_, n_l0_calls = NA_real_,
                             error_msg = NA_character_) {
  data.table(
    method = method, n = n, d_x = d_x, d_z = d_z, seed = seed,
    status = status,
    time_sec = time_sec, mem_mb = mem_mb,
    n_ci_tests = n_ci_tests, n_l0_calls = n_l0_calls,
    precision = NA_real_, recall = NA_real_, f1 = NA_real_,
    dir_acc = NA_real_, coverage = NA_real_,
    tp = NA_real_, fp = NA_real_, fn = NA_real_, tn = NA_real_,
    eval_unresolved = NA_real_, unresolved_frac = NA_real_,
    shd_full = NA_real_, shd_resolved = NA_real_, n_resolved = NA_real_,
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
    seed     = seed
  )
  amat_true <- true_ancestral(sim_obj$adj_xx)
  
  reset_counters()
  gc(reset = TRUE, full = TRUE)
  t0   <- proc.time()[["elapsed"]]
  mem0 <- sum(gc()[, "used"])
  
  est <- tryCatch(
    withTimeout({
      switch(method,
             ascend = ascend(sim_obj, alpha = 0.05, alpha_mb = 0.05,
                             fdr = TRUE, min_votes = 1),
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
  
  ev  <- evaluate(est, amat_true, verbose = FALSE)
  shd <- compute_shd(est, amat_true)
  uf  <- unresolved_frac(est, amat_true)
  
  data.table(
    method = method, n = n, d_x = d_x, d_z = d_z, seed = seed,
    status = "ok",
    time_sec = elapsed, mem_mb = mem_delta,
    n_ci_tests = get_ci_tests(method), n_l0_calls = get_aux_count(method),
    precision = ev$precision, recall = ev$recall, f1 = ev$f1,
    dir_acc = ev$dir_acc, coverage = ev$coverage,
    tp = ev$tp, fp = ev$fp, fn = ev$fn, tn = ev$tn,
    eval_unresolved = ev$unresolved, unresolved_frac = uf,
    shd_full = shd$shd_full, shd_resolved = shd$shd_resolved,
    n_resolved = shd$n_resolved,
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

plan <- CJ(method  = c("ascend", "cbl"),
           cfg_idx = seq_len(nrow(configs)),
           seed    = 100L + seq_len(SEEDS))
plan <- merge(plan, configs[, .(cfg_idx = .I, n, d_x, d_z)], by = "cfg_idx")
setorder(plan, n, d_x, d_z, method, seed)

cat(sprintf("[plan] %d cells total (%d configs x 2 methods x %d seeds)\n",
            nrow(plan), nrow(configs), SEEDS))
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
    cat(sprintf("[%4d/%4d] skip (cached): %-6s n=%-5d d_x=%-3d d_z=%-3d seed=%d\n",
                k, total_cells, row$method, row$n, row$d_x, row$d_z, row$seed))
    flush.console()
    next
  }
  
  cat(sprintf("[%4d/%4d] run : %-6s n=%-5d d_x=%-3d d_z=%-3d seed=%d ... ",
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
    cat(sprintf("OK  %7.2fs  F1=%.3f  SHD=%g  CI=%g\n",
                res$time_sec, res$f1 %||% NA_real_,
                res$shd_full %||% NA_real_, res$n_ci_tests %||% NA_real_))
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