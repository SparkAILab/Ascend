## Benchmark of ASCEND against CBL across three one-dimensional sweeps
## of the simulator parameters, with a shared default point.
##
## Sweeps (default point n = 1024, d_x = 5, d_z = 10):
##   n   in {256, 512, 1024, 2048, 4096}
##   d_x in {5, 10, 15, 20}
##   d_z in {10, 20, 30, 40, 50}
##
## Each cell is run with 5 seeds and a per-cell timeout (default 1 hour).
## Results are appended to the output CSV after every cell, so an
## interrupted run can be resumed.
##
## Usage:
##   Rscript benchmark.R
##
## Environment variables:
##   BENCH_QUICK=1     run only n = 256 with one seed, for smoke testing
##   BENCH_SEEDS=N     override the seed count
##   BENCH_TIMEOUT=S   override the per-cell timeout in seconds

suppressPackageStartupMessages({
  library(R.utils)
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a


#source("ascend.R")
#source("cbl.R")

ASCEND_PATH <- "ascend.R"
CBL_PATH    <- "cbl.R"
RESULTS_CSV <- "results_12.csv"

# Source a file with its trailing example-run block removed. The trimmed
# copy is written next to the original so that relative source() calls
# inside the file (cbl_fixed.R calls source('shah_ss.R') at the top)
# resolve correctly.
source_without_example <- function(path, marker_regex) {
  txt <- readLines(path, warn = FALSE)
  hit <- grep(marker_regex, txt)
  if (length(hit) == 0) {
    warning(sprintf("Could not find example marker in %s; sourcing whole file.", path))
    source(path, local = FALSE)
    return(invisible())
  }
  cutoff <- hit[1] - 1L
  src_dir <- dirname(normalizePath(path))
  tmp <- tempfile(tmpdir = src_dir, fileext = ".R")
  writeLines(txt[seq_len(cutoff)], tmp)
  on.exit(unlink(tmp), add = TRUE)
  source(tmp, local = FALSE, chdir = TRUE)
  invisible()
}

cat("[setup] sourcing ASCEND...\n"); flush.console()
source_without_example(ASCEND_PATH, "^# 8\\. Example run")

cat("[setup] sourcing CBL...\n"); flush.console()
source_without_example(CBL_PATH, "^sim_obj <- sim_dat\\(")

## --- Counting independence verdicts ---------------------------------------
##
## To compare the two methods on the same scale, we count the number of
## conditional independence verdicts each one extracts from the data:
##
##   ASCEND: one verdict per call to ci_test_pval().
##   CBL:    each l0(x, y, ...) call fits one regression that yields a
##           length-ncol(x) selection vector, i.e. one independence
##           verdict per feature in x conditional on the rest.
##
## We wrap both functions in place to maintain three counters per run.
## n_ci_tests is the apples-to-apples count.

.work_counter <- new.env(parent = emptyenv())
.work_counter$ascend_ci_calls <- 0
.work_counter$cbl_l0_calls    <- 0
.work_counter$cbl_ci_tests    <- 0

reset_counters <- function() {
  .work_counter$ascend_ci_calls <- 0
  .work_counter$cbl_l0_calls    <- 0
  .work_counter$cbl_ci_tests    <- 0
}

local({
  orig_ci <- ci_test_pval
  assign("ci_test_pval",
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
  if (method == "ascend") .work_counter$ascend_ci_calls
  else .work_counter$cbl_ci_tests
}
get_aux_count <- function(method) {
  if (method == "ascend") NA_real_
  else .work_counter$cbl_l0_calls
}

## --- Structural Hamming distance and unresolved-pair fraction -------------
##
## ASCEND's evaluate_ancestral() yields precision, recall, F1, coverage,
## and TP/FP/FN/TN. We add two additional metrics here.
##
## SHD compares the upper triangle of the estimated ancestral matrix to
## the truth. shd_full penalises unresolved (NA) pairs at true positives;
## shd_resolved skips NA pairs and only scores cells where the method
## made a positive or negative claim.

compute_shd <- function(estimated, truth) {
  nms <- rownames(truth)
  est <- estimated[nms, nms, drop = FALSE]
  d   <- nrow(truth)
  full <- 0; resolved <- 0; n_res <- 0
  for (i in seq_len(d - 1)) for (j in (i + 1):d) {
    t_ij <- truth[i, j]
    if (is.na(t_ij)) next
    e_ij <- est[i, j]; e_ji <- est[j, i]
    true_pos <- (t_ij == 1)
    is_na    <- is.na(e_ij) && is.na(e_ji)
    pred_pos <- (!is.na(e_ij) && e_ij %in% c(0.5, 1)) ||
      (!is.na(e_ji) && e_ji %in% c(0.5, 1))
    if (is_na) {
      if (true_pos) full <- full + 1
    } else {
      n_res <- n_res + 1
      if (true_pos != pred_pos) {
        full     <- full + 1
        resolved <- resolved + 1
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
  for (i in seq_len(d - 1)) for (j in (i + 1):d) {
    if (is.na(truth[i, j])) next
    total <- total + 1L
    if (is.na(est[i, j])) na_count <- na_count + 1L
  }
  if (total == 0L) NA_real_ else na_count / total
}

## --- Single (method, configuration, seed) run -----------------------------
##
## Generates one simulated dataset, runs the chosen method inside a
## withTimeout() wrapper, and returns one row of metrics. Timeouts and
## errors are caught and returned with status = "timeout" or "error".

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
  amat_true <- get_true_ancestral_matrix(sim_obj$adj_xx)
  
  reset_counters()
  gc(reset = TRUE, full = TRUE)
  t0 <- proc.time()[["elapsed"]]
  mem0 <- sum(gc()[, "used"])
  
  est <- tryCatch(
    withTimeout({
      if (method == "ascend") {
        ascend_fn(sim_obj, maxiter = 10, alpha = 0.05,
                  alpha_mb_start = 0.20, alpha_mb_floor = 0.05,
                  alpha_decay = 0.70, fdr_correction = TRUE,
                  min_votes = 1)
      } else if (method == "cbl") {
        cbl_fn(sim_obj, gamma = 0.5, maxiter = 100, B = 50)
      } else {
        stop("Unknown method: ", method)
      }
    }, timeout = timeout_sec, onTimeout = "error"),
    error = function(e) {
      structure(NA, error_msg = conditionMessage(e))
    }
  )
  
  t1 <- proc.time()[["elapsed"]]
  mem1 <- sum(gc()[, "used"])
  elapsed <- t1 - t0
  mem_delta <- max(0, mem1 - mem0)
  
  if (is.null(dim(est))) {
    status <- if (grepl("reached elapsed time limit|TimeoutException",
                        attr(est, "error_msg") %||% "")) "timeout" else "error"
    return(data.table(
      method = method, n = n, d_x = d_x, d_z = d_z, seed = seed,
      status = status,
      time_sec = elapsed, mem_mb = mem_delta,
      n_ci_tests = get_ci_tests(method),
      n_l0_calls = get_aux_count(method),
      precision = NA_real_, recall_resolved = NA_real_, recall_overall = NA_real_,
      f1 = NA_real_, accuracy = NA_real_, coverage = NA_real_,
      tp = NA_real_, fp = NA_real_, fn = NA_real_, tn = NA_real_,
      unresolved_frac = NA_real_,
      shd_full = NA_real_, shd_resolved = NA_real_, n_resolved = NA_real_,
      error_msg = attr(est, "error_msg") %||% NA_character_
    ))
  }
  
  ev  <- evaluate_ancestral(est, amat_true, verbose = FALSE)
  shd <- compute_shd(est, amat_true)
  uf  <- unresolved_frac(est, amat_true)
  
  data.table(
    method = method, n = n, d_x = d_x, d_z = d_z, seed = seed,
    status = "ok",
    time_sec = elapsed, mem_mb = mem_delta,
    n_ci_tests = get_ci_tests(method),
    n_l0_calls = get_aux_count(method),
    precision = ev$precision, recall_resolved = ev$recall_resolved,
    recall_overall = ev$recall_overall, f1 = ev$f1, accuracy = ev$accuracy,
    coverage = ev$coverage,
    tp = ev$tp, fp = ev$fp, fn = ev$fn, tn = ev$tn,
    unresolved_frac = uf,
    shd_full = shd$shd_full, shd_resolved = shd$shd_resolved,
    n_resolved = shd$n_resolved,
    error_msg = NA_character_
  )
}

## --- Sweep definition -----------------------------------------------------

DEFAULT_DX <- 5L
DEFAULT_DZ <- 10L
DEFAULT_N  <- 1024L

quick <- Sys.getenv("BENCH_QUICK", "0") == "1"
SEEDS <- as.integer(Sys.getenv("BENCH_SEEDS",
                               if (quick) "1" else "5"))
TIMEOUT_SEC <- as.numeric(Sys.getenv("BENCH_TIMEOUT",
                                     if (quick) "120" else "3600"))

n_vals  <- if (quick) c(256L) else c(256L, 512L, 1024L, 2048L, 4096L)
dx_vals <- if (quick) c(5L)   else c(5L, 10L, 15L, 20L)
dz_vals <- if (quick) c(10L)  else c(10L, 20L, 30L, 40L, 50L)

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

## --- Resume support -------------------------------------------------------
## A cell is skipped if it already appears in the CSV with status "ok"
## or "timeout". Errored cells are retried since they may be transient.

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

## --- Main loop ------------------------------------------------------------

header_written <- file.exists(RESULTS_CSV) && file.size(RESULTS_CSV) > 0L
total_cells <- nrow(plan)
t_sweep0    <- proc.time()[["elapsed"]]

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
      data.table(
        method = row$method, n = row$n, d_x = row$d_x, d_z = row$d_z,
        seed = row$seed, status = "error",
        time_sec = NA_real_, mem_mb = NA_real_,
        n_ci_tests = NA_real_, n_l0_calls = NA_real_,
        precision = NA_real_, recall_resolved = NA_real_, recall_overall = NA_real_,
        f1 = NA_real_, accuracy = NA_real_, coverage = NA_real_,
        tp = NA_real_, fp = NA_real_, fn = NA_real_, tn = NA_real_,
        unresolved_frac = NA_real_,
        shd_full = NA_real_, shd_resolved = NA_real_, n_resolved = NA_real_,
        error_msg = conditionMessage(e)
      )
    }
  )
  
  fwrite(res, RESULTS_CSV, append = header_written)
  header_written <- TRUE
  
  if (res$status == "ok") {
    cat(sprintf("OK  %7.2fs  F1=%.3f  SHD=%g  CI=%g\n",
                res$time_sec, res$f1 %||% NA_real_,
                res$shd_full %||% NA_real_,
                res$n_ci_tests %||% NA_real_))
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