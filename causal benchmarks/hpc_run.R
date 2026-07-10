#!/usr/bin/env Rscript
# =============================================================================
# hpc_run.R — ASCEND benchmark, KCL CREATE HPC  (729-task version)
#
# Usage (SLURM array, task id 1 .. 81*length(N_VEC)):
#   Rscript hpc_run.R <task_id>
#
# TASK LAYOUT
#   There are NCOMBO parameter combinations (see make_grid) and length(N_VEC)
#   sample sizes. A task does ALL replicates for ONE (combo, n) pair. Tasks are
#   laid out so that the NCOMBO combos at a given n are CONTIGUOUS:
#
#       task_id 1..NCOMBO            -> n = N_VEC[1]   (all combos)
#       task_id NCOMBO+1..2*NCOMBO   -> n = N_VEC[2]   (all combos)
#       ...
#
#   This lets submit_all.sh give each n-band its own --time / --mem, so the
#   cheap small-n tasks don't sit under a wall-time sized for n = 131072.
#   Run `Rscript hpc_run.R count` to print "NCOMBO NVALS NTASKS".
#
#   Decode:
#       n_index = ((task_id - 1) %/% NCOMBO) + 1   # which N_VEC entry
#       combo   = ((task_id - 1) %%  NCOMBO) + 1   # which GRID row
#
# Parameter grid (81 combos):
#   r2  : 1/2, 2/3, 3/4      d_x : 20, 40, 80
#   d_z : d_x, 2*d_x, 3*d_x  sp  : 1/2, 2/3, 3/4
#
# Methods: ASCEND, CBL, GES, LiNGAM, PC  (Framework 1: ancestral recovery)
#
# Each method has a hard per-replicate timeout (METHOD_TIMEOUT seconds).
# Once a method times out it is skipped for the rest of this task, and a marker
# file is written so larger-n tasks of the SAME combo skip it too (best effort).
#
# Results saved immediately after each replicate — safe to restart.
#
# Output:
#   results/job_<combo>/params.rds          — parameters for this combo
#   results/job_<combo>/n<n>_rep<r>.rds     — data.frame, one row per method
#   results/job_<combo>/timeout_<METHOD>.rds— min n at which <METHOD> timed out
#   results/job_<combo>/done_n<n>.rds       — written when this (combo,n) finishes
#
# RESULT COLUMNS (per method per replicate)
#   job, task_id, method, n, rep, d_x, d_z, r2, sp, n_over_dz,
#   precision, recall, f1, accuracy, coverage,
#   tp, fp, fn, tn, unres_tp, unres_tn,
#   elapsed_sec,        <- wall-clock time for THIS method on THIS replicate
#   status              <- "ok" | "timeout" | "error" | "skipped_timeout" | "failed_subset"
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(pcalg)
  library(bnlearn)
  library(igraph)
  library(data.table)
  library(matrixStats)
  library(dplyr)
  library(R.utils)
})

# ── 1. Constants ──────────────────────────────────────────────────────────────

METHOD_TIMEOUT <- 3600L          # 1 hour per method per replicate
N_REP          <- 20L
N_VEC          <- 2L^(9:17)      # 512 … 131072.  EDIT THIS LINE to cap (e.g. 2L^(9:16));
# the array size derives from it automatically.
P_CROSS        <- 0.20
X_EFFECT       <- 0.90

# ── 2. Parameter grid (81 combinations) ──────────────────────────────────────

make_grid <- function() {
  grid <- expand.grid(
    r2      = c(1/2, 2/3, 3/4),
    d_x     = c(20L, 40L, 80L),
    dz_mult = c(1L,  2L,  3L),
    sp      = c(1/2, 2/3, 3/4),
    stringsAsFactors = FALSE
  )
  grid$d_z     <- grid$dz_mult * grid$d_x
  grid$dz_mult <- NULL
  grid <- grid[order(grid$d_x, grid$d_z, grid$r2, grid$sp), ]
  rownames(grid) <- seq_len(nrow(grid))
  grid
}
GRID    <- make_grid()
NCOMBO  <- nrow(GRID)
stopifnot(NCOMBO == 81L)
N_TASKS <- NCOMBO * length(N_VEC)

# ── 3. Decode task id from command line ───────────────────────────────────────

args <- commandArgs(trailingOnly=TRUE)
# `Rscript hpc_run.R count` prints "NCOMBO NVALS NTASKS" — used by submit_all.sh
# so the array size tracks the grid automatically when you edit make_grid/N_VEC.
if (length(args) >= 1L && identical(args[1L], "count")) {
  cat(sprintf("%d %d %d\n", NCOMBO, length(N_VEC), N_TASKS)); quit(save="no")
}
if (length(args) < 1L)
  stop(sprintf("Usage: Rscript hpc_run.R <task_id>  [1-%d]", N_TASKS))
TASK_ID <- as.integer(args[1L])
if (is.na(TASK_ID) || TASK_ID < 1L || TASK_ID > N_TASKS)
  stop(sprintf("task_id must be 1-%d, got: %s", N_TASKS, args[1L]))

N_INDEX  <- ((TASK_ID - 1L) %/% NCOMBO) + 1L
COMBO    <- ((TASK_ID - 1L) %%  NCOMBO) + 1L
N_VAL    <- N_VEC[N_INDEX]
params   <- as.list(GRID[COMBO, ])

cat(sprintf("\n[Task %d] combo=%d  n=%d  |  r2=%.4f  d_x=%d  d_z=%d  sp=%.4f\n",
            TASK_ID, COMBO, N_VAL, params$r2, params$d_x, params$d_z, params$sp))

# ── 4. Output directory (one per combo; tasks differ by the n<n>_rep* files) ──

OUT_DIR <- file.path("results", sprintf("job_%03d", COMBO))
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
saveRDS(params, file.path(OUT_DIR, "params.rds"))

# ── 5. Source algorithm files ─────────────────────────────────────────────────
# Source the rewritten ASCEND whether it is named ascend.R or ascend_complete.R.
ascend_src <- Filter(file.exists, c("ascend.R", "ascend_complete.R"))
if (!length(ascend_src))
  stop("Cannot find ascend.R or ascend_complete.R in the working directory")
source(ascend_src[1])              # sim_dat, ascend, true_ancestral, evaluate, ...
cat(sprintf("[Setup] Sourced ASCEND from: %s\n", ascend_src[1]))

source("cbl.R")              # cbl_fn is found at ../cbl.R

for (fn in c("sim_dat", "ascend", "cbl_fn", "minD")) {
  if (!exists(fn)) stop(sprintf("Required function not found: %s", fn))
}
cat("[Setup] All required functions loaded.\n")

# ── 6. Utilities ──────────────────────────────────────────────────────────────

xcols <- function(dat) sort(grep("^x", colnames(dat), value=TRUE), method="radix")
zcols <- function(dat) sort(grep("^z", colnames(dat), value=TRUE), method="radix")

transitive_closure <- function(A) {
  A[is.na(A)] <- 0; diag(A) <- 0
  d <- nrow(A); R <- (A>0)*1; Ap <- A
  for (k in seq_len(d-1L)) { Ap <- ((Ap%*%A)>0)*1; R <- ((R+Ap)>0)*1 }
  diag(R) <- NA; R
}

# defined AFTER source() so it shadows ascend.R's true_ancestral(adj_xx)
true_ancestral <- function(sim_obj) {
  A <- t(sim_obj$adj_xx); A[is.na(A)] <- 0; diag(A) <- 0
  nms <- colnames(A); R <- transitive_closure(A)
  rownames(R) <- colnames(R) <- nms; R
}

cpdag_to_ancestral <- function(cpdag, nms) {
  d <- nrow(cpdag)
  cpdag[is.na(cpdag)] <- 0            # NA-safe: PC's amat can carry NA/conflict marks,
  # which made `if (cpdag[i,j]==1)` throw below.
  A <- matrix(0, d, d, dimnames=list(nms, nms))
  for (i in seq_len(d)) for (j in seq_len(d))
    if (i!=j && cpdag[i,j]==1 && cpdag[j,i]==0) A[i,j] <- 1
  R <- transitive_closure(A); rownames(R) <- colnames(R) <- nms
  R[is.na(R)] <- 0
  for (i in seq_len(d)) for (j in seq_len(d))
    if (i!=j && cpdag[i,j]==1 && cpdag[j,i]==1 && R[i,j]==0 && R[j,i]==0)
      R[i,j] <- NA
  diag(R) <- NA; R
}

make_suffstat <- function(mat) {
  mat <- mat[complete.cases(mat) & apply(mat,1,function(r) all(is.finite(r))),
             , drop=FALSE]
  col_sd <- apply(mat, 2, sd, na.rm=TRUE)
  mat    <- mat[, col_sd > 1e-10, drop=FALSE]
  n <- nrow(mat); p <- ncol(mat)
  if (p < 2L) stop("Fewer than 2 non-constant columns after cleaning")
  C <- cor(mat)
  C[is.nan(C) | is.na(C)] <- 0
  diag(C) <- 1
  ev <- min(eigen(C, only.values=TRUE)$values)
  if (ev <= 1e-10) {
    C  <- C + (abs(ev)+1e-4)*diag(p)
    dv <- sqrt(diag(C)); C <- C/outer(dv,dv)
  }
  list(C=C, n=n, labels=colnames(mat))
}

# ── 7. Evaluation ─────────────────────────────────────────────────────────────

ancestral_metrics <- function(M_est, M_true) {
  nms <- rownames(M_true); M_est <- M_est[nms,nms]; d <- nrow(M_true)
  tp<-0; fp<-0; fn<-0; tn<-0; unres_tp<-0; unres_tn<-0
  for (i in seq_len(d-1L)) for (j in (i+1L):d) {
    tv <- M_true[i,j]; ev <- M_est[i,j]; if (is.na(tv)) next
    tpos <- (tv==1)
    if (is.na(ev)) {
      if (tpos) unres_tp<-unres_tp+1 else unres_tn<-unres_tn+1; next
    }
    epos <- (ev==1||ev==0.5)
    if( tpos && epos)  tp<-tp+1
    if(!tpos && epos)  fp<-fp+1
    if( tpos && !epos) fn<-fn+1
    if(!tpos && !epos) tn<-tn+1
  }
  tot <- tp+fp+fn+tn+unres_tp+unres_tn; res <- tp+fp+fn+tn
  pr <- if ((tp+fp)>0)          tp/(tp+fp)          else NA_real_
  re <- if ((tp+fn+unres_tp)>0) tp/(tp+fn+unres_tp) else NA_real_
  f1 <- if (!is.na(pr)&&!is.na(re)&&(pr+re)>0) 2*pr*re/(pr+re) else NA_real_
  list(precision=pr, recall=re, f1=f1,
       accuracy =if(res>0)(tp+tn)/res else NA_real_,
       coverage =if(tot>0)res/tot     else NA_real_,
       tp=tp, fp=fp, fn=fn, tn=tn,
       unres_tp=unres_tp, unres_tn=unres_tn)
}

# ── 8. Method runners ─────────────────────────────────────────────────────────

# returns list(value=..., elapsed=<sec>, status="ok"|"timeout"|"error")
run_timed <- function(label, expr, d_x, d_z, n_val) {
  cat(sprintf("    [%-8s] d_x=%d d_z=%d n=%d ... ", label, d_x, d_z, n_val))
  flush.console()
  t0 <- proc.time()["elapsed"]
  status <- "ok"
  value <- tryCatch(
    withTimeout(expr, timeout=METHOD_TIMEOUT, onTimeout="error"),
    TimeoutException = function(e) {
      cat(sprintf("TIMEOUT (>%dh)\n", METHOD_TIMEOUT %/% 3600L))
      flush.console(); status <<- "timeout"; NULL
    },
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      flush.console(); status <<- "error"; NULL
    }
  )
  list(value=value, elapsed=as.numeric(proc.time()["elapsed"]-t0), status=status)
}

run_ascend <- function(sim_obj, d_x, d_z, n_val) {
  run_timed("ASCEND", {
    invisible(capture.output(
      M <- ascend(sim_obj, maxiter=12L, alpha=0.05,
                  alpha_mb=0.05, fdr=TRUE, min_votes=1L,
                  prescreen=0.30, verbose=FALSE)
    )); M
  }, d_x, d_z, n_val)
}

run_cbl <- function(sim_obj, d_x, d_z, n_val) {
  run_timed("CBL", {
    M_raw <- cbl_fn(sim_obj, gamma=0.5, maxiter=10L, B=20L)
    if (is.null(M_raw)) { message("CBL returned NULL"); return(NULL) }
    dat    <- as.data.frame(sim_obj$dat)
    x_cols <- xcols(dat)
    if (is.null(rownames(M_raw))) rownames(M_raw) <- x_cols
    if (is.null(colnames(M_raw))) colnames(M_raw) <- x_cols
    M_raw
  }, d_x, d_z, n_val)
}

run_ges <- function(sim_obj, d_x, d_z, n_val) {
  run_timed("GES", {
    dat      <- as.data.frame(sim_obj$dat)
    x_cols   <- xcols(dat); z_cols <- zcols(dat)
    all_cols <- c(z_cols, x_cols)
    mat      <- as.matrix(dat[, all_cols])
    mat      <- mat[complete.cases(mat) & apply(mat,1,function(r) all(is.finite(r))),]
    if (nrow(mat) < ncol(mat)+3L) stop("Insufficient rows for GES")
    score <- new("GaussL0penObsScore", mat)
    fit   <- ges(score, labels=all_cols, iterate=TRUE, verbose=FALSE)
    ess   <- as(fit$essgraph, "matrix")
    rownames(ess) <- colnames(ess) <- all_cols; diag(ess) <- 0
    cpdag_to_ancestral(ess[x_cols, x_cols, drop=FALSE], x_cols)
  }, d_x, d_z, n_val)
}

run_lingam <- function(sim_obj, d_x, d_z, n_val) {
  run_timed("LiNGAM", {
    dat      <- as.data.frame(sim_obj$dat)
    x_cols   <- xcols(dat); z_cols <- zcols(dat)
    all_cols <- c(z_cols, x_cols); total <- length(all_cols)
    fit      <- pcalg::lingam(as.matrix(dat[, all_cols]))
    B        <- if (!is.null(fit$Bpruned)) fit$Bpruned else fit$B
    A        <- matrix(0, total, total, dimnames=list(all_cols,all_cols))
    for (i in seq_len(total)) for (j in seq_len(total))
      if (i!=j && abs(B[i,j])>1e-9) A[j,i] <- 1
    A_xx <- A[x_cols, x_cols, drop=FALSE]
    R    <- transitive_closure(A_xx)
    rownames(R) <- colnames(R) <- x_cols
    R[is.na(R)] <- 0; diag(R) <- NA; R
  }, d_x, d_z, n_val)
}

run_pc <- function(sim_obj, d_x, d_z, n_val, alpha_pc=0.05) {
  run_timed("PC", {
    dat      <- as.data.frame(sim_obj$dat)
    x_cols   <- xcols(dat); z_cols <- zcols(dat)
    all_cols <- c(z_cols, x_cols); total <- length(all_cols)
    ss       <- make_suffstat(as.matrix(dat[, all_cols]))
    kept_labels <- if (!is.null(ss$labels)) ss$labels else all_cols
    if (ss$n < length(kept_labels)+3L) stop("Insufficient observations for PC")
    x_kept <- intersect(x_cols, kept_labels)
    if (length(x_kept) < 2L) stop("Fewer than 2 X variables after filtering")
    fit <- pc(suffStat=ss, indepTest=gaussCItest, labels=kept_labels,
              alpha=alpha_pc, verbose=FALSE, maj.rule=TRUE, solve.confl=TRUE)
    amat <- as(fit, "amat")
    rownames(amat) <- colnames(amat) <- kept_labels
    cpdag_to_ancestral(amat[x_kept, x_kept, drop=FALSE], x_kept)
  }, d_x, d_z, n_val)
}

# ── 9. Record result ──────────────────────────────────────────────────────────

make_row <- function(method, M, M_true, x_cols, n_val, rep_id,
                     d_x, d_z, r2, sp, elapsed_sec=NA_real_, status="ok") {
  n_over_dz <- n_val / d_z
  make_stub <- function(st) data.frame(
    job=COMBO, task_id=TASK_ID, method=method, n=n_val, rep=rep_id,
    d_x=d_x, d_z=d_z, r2=r2, sp=sp, n_over_dz=n_over_dz,
    precision=NA_real_, recall=NA_real_, f1=NA_real_,
    accuracy=NA_real_, coverage=NA_real_,
    tp=NA_integer_, fp=NA_integer_, fn=NA_integer_, tn=NA_integer_,
    unres_tp=NA_integer_, unres_tn=NA_integer_,
    elapsed_sec=elapsed_sec, status=st, stringsAsFactors=FALSE
  )
  if (is.null(M)) return(make_stub(status))             # status carries timeout/error/skip
  M <- tryCatch(M[x_cols, x_cols], error=function(e) NULL)
  if (is.null(M)) return(make_stub("failed_subset"))
  mt  <- ancestral_metrics(M, M_true)
  pr  <- if (is.na(mt$precision)) 0 else mt$precision
  re  <- if (is.na(mt$recall))    0 else mt$recall
  f1  <- if (is.na(mt$f1))        0 else mt$f1
  cov <- if (is.na(mt$coverage))  0 else mt$coverage
  cat(sprintf("done (%.1fs)  pr=%.3f  re=%.3f  f1=%.3f  cov=%.3f\n",
              elapsed_sec, pr, re, f1, cov))
  flush.console()
  data.frame(
    job=COMBO, task_id=TASK_ID, method=method, n=n_val, rep=rep_id,
    d_x=d_x, d_z=d_z, r2=r2, sp=sp, n_over_dz=n_over_dz,
    precision=mt$precision, recall=mt$recall, f1=mt$f1,
    accuracy=mt$accuracy, coverage=mt$coverage,
    tp=mt$tp, fp=mt$fp, fn=mt$fn, tn=mt$tn,
    unres_tp=mt$unres_tp, unres_tn=mt$unres_tn,
    elapsed_sec=elapsed_sec, status="ok", stringsAsFactors=FALSE
  )
}

# ── 10. Single replicate ──────────────────────────────────────────────────────
# `dead` (job-scoped env) = methods that timed out earlier in THIS task.
# Marker files extend the skip to larger-n tasks of the same combo (best effort).

run_one_rep <- function(sim_obj, n_val, rep_id, dead) {
  dat    <- as.data.frame(sim_obj$dat)
  x_cols <- xcols(dat)
  M_true <- true_ancestral(sim_obj)[x_cols, x_cols]
  dx     <- params$d_x; dz <- params$d_z
  r2     <- params$r2;  sp <- params$sp
  
  # CBL reads sim_obj$params$lin_pr; the rewritten sim_dat no longer stores it.
  if (is.null(sim_obj$params$lin_pr)) sim_obj$params$lin_pr <- 1
  
  marker_path <- function(label) file.path(OUT_DIR, sprintf("timeout_%s.rds", label))
  is_dead <- function(label) {
    if (label %in% dead$methods) return(TRUE)
    mp <- marker_path(label)
    if (file.exists(mp)) {
      n0 <- tryCatch(readRDS(mp), error=function(e) NA_real_)
      if (!is.na(n0) && n_val >= n0) return(TRUE)
    }
    FALSE
  }
  mark_dead <- function(label) {
    dead$methods <- union(dead$methods, label)
    mp   <- marker_path(label)
    prev <- if (file.exists(mp)) tryCatch(readRDS(mp), error=function(e) Inf) else Inf
    saveRDS(min(prev, n_val), mp)
  }
  
  run_or_skip <- function(label, runner) {
    if (is_dead(label)) {
      cat(sprintf("    [%-8s] d_x=%d d_z=%d n=%d ... SKIPPED (timed out at <= this n)\n",
                  label, dx, dz, n_val))
      return(make_row(label, NULL, M_true, x_cols, n_val, rep_id,
                      dx, dz, r2, sp, elapsed_sec=0, status="skipped_timeout"))
    }
    tr <- runner(sim_obj, dx, dz, n_val)
    if (identical(tr$status, "timeout")) mark_dead(label)
    make_row(label, tr$value, M_true, x_cols, n_val, rep_id,
             dx, dz, r2, sp, elapsed_sec=tr$elapsed, status=tr$status)
  }
  
  rows <- list(
    run_or_skip("ASCEND", run_ascend),
    run_or_skip("CBL",    run_cbl),
    run_or_skip("GES",    run_ges),
    run_or_skip("LiNGAM", run_lingam),
    run_or_skip("PC",     run_pc)
  )
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)==0L) return(NULL)
  do.call(rbind, rows)
}

# ── 11. Main loop (single n; loop over replicates) ────────────────────────────

cat(sprintf("[Task %d] Output dir : %s\n", TASK_ID, OUT_DIR))
cat(sprintf("[Task %d] n          : %d  (n/d_z = %.1f)\n", TASK_ID, N_VAL, N_VAL/params$d_z))
cat(sprintf("[Task %d] n_rep      : %d\n", TASK_ID, N_REP))
cat(sprintf("[Task %d] timeout    : %ds per method\n", TASK_ID, METHOD_TIMEOUT))

dead <- new.env(parent=emptyenv()); dead$methods <- character(0)

for (rep in seq_len(N_REP)) {
  out_file <- file.path(OUT_DIR, sprintf("n%d_rep%02d.rds", N_VAL, rep))
  
  if (file.exists(out_file)) {
    cat(sprintf("  rep %02d — already saved, skipping\n", rep))
    next
  }
  
  cat(sprintf("\n  rep %02d/%02d\n", rep, N_REP))
  t0 <- proc.time()["elapsed"]
  
  sim <- tryCatch(
    sim_dat(n=N_VAL, d_z=params$d_z, d_x=params$d_x,
            r2=params$r2, lin_pr=1,  sp=params$sp,
            p_cross=P_CROSS, x_effect=X_EFFECT,
            seed=COMBO*1000000L + N_VAL*1000L + rep),
    error=function(e) { message("sim_dat: ", e$message); NULL }
  )
  if (is.null(sim)) { cat("    sim_dat failed — skipping rep\n"); next }
  
  result <- tryCatch(
    run_one_rep(sim, N_VAL, rep, dead),
    error=function(e) { message("run_one_rep: ", e$message); NULL }
  )
  
  if (!is.null(result) && nrow(result) > 0L) {
    saveRDS(result, out_file)
    cat(sprintf("  rep %02d done (%.1fs total) — %d/%d methods ok\n",
                rep, proc.time()["elapsed"]-t0,
                sum(result$status=="ok", na.rm=TRUE), nrow(result)))
  } else {
    cat(sprintf("  rep %02d — no results (%.1fs)\n",
                rep, proc.time()["elapsed"]-t0))
  }
  
  rm(sim, result); gc(verbose=FALSE)
}

saveRDS(list(combo=COMBO, n=N_VAL, params=params, completed=Sys.time()),
        file.path(OUT_DIR, sprintf("done_n%d.rds", N_VAL)))
cat(sprintf("\n[Task %d] Complete (combo=%d, n=%d).\n", TASK_ID, COMBO, N_VAL))