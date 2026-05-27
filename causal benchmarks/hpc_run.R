#!/usr/bin/env Rscript
# =============================================================================
# hpc_run.R — ASCEND benchmark, KCL CREATE HPC
#
# Usage (SLURM array, job index 1-81):
#   Rscript hpc_run.R <job_index>
#
# Parameter grid (81 jobs):
#   r2      : 1/2, 2/3, 3/4
#   d_x     : 20, 40, 80
#   d_z     : d_x, 2*d_x, 3*d_x
#   sp      : 1/2, 2/3, 3/4
#
# Per job: n in 2^(9:17), N_REP=20 replicates each
# Methods: ASCEND, CBL, GES, LiNGAM, PC  (Framework 1: ancestral recovery)
#
# Each method has a hard per-replicate timeout (METHOD_TIMEOUT seconds).
# If a method times out or errors it is recorded as NA and skipped.
# Results are saved immediately after each replicate — safe to restart.
#
# Output:
#   results/job_<idx>/params.rds        — parameters for this job
#   results/job_<idx>/n<n>_rep<r>.rds   — data.frame, one row per method
#   results/job_<idx>/done.rds          — written on clean completion
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
N_VEC          <- 2L^(9:17)     # 512 … 131072
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
GRID <- make_grid()
stopifnot(nrow(GRID) == 81L)

# ── 3. Job index from command line ────────────────────────────────────────────

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1L) stop("Usage: Rscript hpc_run.R <job_index>  [1-81]")
JOB_IDX <- as.integer(args[1L])
if (is.na(JOB_IDX) || JOB_IDX < 1L || JOB_IDX > 81L)
  stop(sprintf("job_index must be 1-81, got: %s", args[1L]))

params <- as.list(GRID[JOB_IDX, ])
cat(sprintf("\n[Job %d] r2=%.4f  d_x=%d  d_z=%d  sp=%.4f\n",
            JOB_IDX, params$r2, params$d_x, params$d_z, params$sp))

# ── 4. Output directory ───────────────────────────────────────────────────────

OUT_DIR <- file.path("results", sprintf("job_%03d", JOB_IDX))
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
saveRDS(params, file.path(OUT_DIR, "params.rds"))

# ── 5. Source algorithm files ─────────────────────────────────────────────────

source("ascend.R")     # sim_dat, ascend_fn, get_true_ancestral_matrix,
# evaluate_ancestral, normalize_to_sorted_ancestral
source("cbl_fixed.R")  # cbl_fn (shah_ss.R embedded)

for (fn in c("sim_dat", "ascend_fn", "cbl_fn", "minD")) {
  if (!exists(fn)) stop(sprintf("Required function not found: %s", fn))
}
cat("[Setup] All required functions loaded.\n")

# ── 6. Utilities ──────────────────────────────────────────────────────────────

# Sort column names numerically: x1,x2,...,x10 not x1,x10,x11
xcols <- function(dat) sort(grep("^x", colnames(dat), value=TRUE), method="radix")
zcols <- function(dat) sort(grep("^z", colnames(dat), value=TRUE), method="radix")

transitive_closure <- function(A) {
  A[is.na(A)] <- 0; diag(A) <- 0
  d <- nrow(A); R <- (A>0)*1; Ap <- A
  for (k in seq_len(d-1L)) { Ap <- ((Ap%*%A)>0)*1; R <- ((R+Ap)>0)*1 }
  diag(R) <- NA; R
}

true_ancestral <- function(sim_obj) {
  A <- t(sim_obj$adj_xx); A[is.na(A)] <- 0; diag(A) <- 0
  nms <- colnames(A); R <- transitive_closure(A)
  rownames(R) <- colnames(R) <- nms; R
}

# CPDAG -> ancestral: directed edges -> transitive closure, undirected -> NA
cpdag_to_ancestral <- function(cpdag, nms) {
  d <- nrow(cpdag)
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
  # Remove rows with any NA or non-finite value
  mat <- mat[complete.cases(mat) & apply(mat,1,function(r) all(is.finite(r))),
             , drop=FALSE]
  # Remove columns with zero or near-zero variance (causes NA in cor())
  col_sd <- apply(mat, 2, sd, na.rm=TRUE)
  mat    <- mat[, col_sd > 1e-10, drop=FALSE]
  n <- nrow(mat); p <- ncol(mat)
  if (p < 2L) stop("Fewer than 2 non-constant columns after cleaning")
  C <- cor(mat)
  # Replace any remaining NaN/NA (from perfectly collinear columns) with 0
  C[is.nan(C) | is.na(C)] <- 0
  diag(C) <- 1
  # Regularise if not positive definite
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

# Wrapper: run with hard timeout.
# Prints header, returns list(value=..., elapsed=...) or NULL on failure.
run_timed <- function(label, expr, d_x, d_z, n_val) {
  cat(sprintf("    [%-8s] d_x=%d d_z=%d n=%d ... ", label, d_x, d_z, n_val))
  flush.console()
  t0 <- proc.time()["elapsed"]
  value <- tryCatch(
    withTimeout(expr, timeout=METHOD_TIMEOUT, onTimeout="error"),
    TimeoutException = function(e) {
      cat(sprintf("TIMEOUT (>%dh)\n", METHOD_TIMEOUT %/% 3600L))
      flush.console(); NULL
    },
    error = function(e) {
      cat(sprintf("ERROR: %s\n", conditionMessage(e)))
      flush.console(); NULL
    }
  )
  list(value=value, elapsed=proc.time()["elapsed"]-t0)
}

# ASCEND
run_ascend <- function(sim_obj, d_x, d_z, n_val) {
  run_timed("ASCEND", {
    invisible(capture.output(
      M <- ascend_fn(sim_obj, maxiter=12L, alpha=0.05,
                     alpha_mb_start=0.15, alpha_mb_floor=0.05,
                     alpha_decay=0.80, fdr_correction=TRUE, min_votes=1L)
    )); M
  }, d_x, d_z, n_val)
}

# CBL
run_cbl <- function(sim_obj, d_x, d_z, n_val) {
  run_timed("CBL", {
    M_raw <- cbl_fn(sim_obj, gamma=0.5, maxiter=10L, B=20L)
    if (is.null(M_raw)) { message("CBL returned NULL"); return(NULL) }
    dat    <- as.data.frame(sim_obj$dat)
    x_cols <- xcols(dat)
    if (is.null(rownames(M_raw))) rownames(M_raw) <- x_cols
    if (is.null(colnames(M_raw))) colnames(M_raw) <- x_cols
    # CBL uses same convention as ASCEND: M[i,j]=1 means i is ancestor of j
    # No transposition needed — return as-is
    M_raw
  }, d_x, d_z, n_val)
}

# GES — no fixedGaps (avoids NA covariance errors); extracts X block
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

# LiNGAM
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

# PC — no fixedGaps (symmetric fixedGaps causes pcalg errors in some configs).
# Two-tier structure is enforced by the data: Z does not respond to X.
run_pc <- function(sim_obj, d_x, d_z, n_val, alpha_pc=0.05) {
  run_timed("PC", {
    dat      <- as.data.frame(sim_obj$dat)
    x_cols   <- xcols(dat); z_cols <- zcols(dat)
    all_cols <- c(z_cols, x_cols); total <- length(all_cols)
    ss       <- make_suffstat(as.matrix(dat[, all_cols]))
    # Use the labels that survived variance filtering in make_suffstat
    kept_labels <- if (!is.null(ss$labels)) ss$labels else all_cols
    if (ss$n < length(kept_labels)+3L) stop("Insufficient observations for PC")
    # Only proceed if enough X columns survived filtering
    x_kept <- intersect(x_cols, kept_labels)
    if (length(x_kept) < 2L) stop("Fewer than 2 X variables after filtering")
    fit <- pc(suffStat=ss, indepTest=gaussCItest, labels=kept_labels,
              alpha=alpha_pc, verbose=FALSE, maj.rule=TRUE, solve.confl=TRUE)
    amat <- as(fit, "amat")
    rownames(amat) <- colnames(amat) <- kept_labels
    # Extract X block using only columns that survived
    cpdag_to_ancestral(amat[x_kept, x_kept, drop=FALSE], x_kept)
  }, d_x, d_z, n_val)
}

# ── 9. Record result ──────────────────────────────────────────────────────────

make_row <- function(method, M, M_true, x_cols, n_val, rep_id, job_idx,
                     d_x, d_z, r2, sp, elapsed=NA) {
  make_failed <- function(status) data.frame(
    job=job_idx, method=method, n=n_val, rep=rep_id,
    d_x=d_x, d_z=d_z, r2=r2, sp=sp,
    precision=NA_real_, recall=NA_real_, f1=NA_real_,
    accuracy=NA_real_, coverage=NA_real_,
    tp=NA_integer_, fp=NA_integer_, fn=NA_integer_, tn=NA_integer_,
    unres_tp=NA_integer_, unres_tn=NA_integer_,
    status=status, stringsAsFactors=FALSE
  )
  if (is.null(M)) return(make_failed("failed"))
  M <- tryCatch(M[x_cols, x_cols], error=function(e) NULL)
  if (is.null(M)) return(make_failed("failed_subset"))
  mt  <- ancestral_metrics(M, M_true)
  pr  <- if (is.na(mt$precision)) 0 else mt$precision
  re  <- if (is.na(mt$recall))    0 else mt$recall
  f1  <- if (is.na(mt$f1))        0 else mt$f1
  cov <- if (is.na(mt$coverage))  0 else mt$coverage
  dt_str <- if (!is.na(elapsed)) sprintf("done (%.1fs)", elapsed) else "done"
  cat(sprintf("%s  pr=%.3f  re=%.3f  f1=%.3f  cov=%.3f\n",
              dt_str, pr, re, f1, cov))
  flush.console()
  data.frame(
    job=job_idx, method=method, n=n_val, rep=rep_id,
    d_x=d_x, d_z=d_z, r2=r2, sp=sp,
    precision=mt$precision, recall=mt$recall, f1=mt$f1,
    accuracy=mt$accuracy, coverage=mt$coverage,
    tp=mt$tp, fp=mt$fp, fn=mt$fn, tn=mt$tn,
    unres_tp=mt$unres_tp, unres_tn=mt$unres_tn,
    status="ok", stringsAsFactors=FALSE
  )
}

# ── 10. Single replicate ──────────────────────────────────────────────────────

run_one_rep <- function(sim_obj, n_val, rep_id) {
  dat    <- as.data.frame(sim_obj$dat)
  x_cols <- xcols(dat)
  M_true <- true_ancestral(sim_obj)[x_cols, x_cols]
  dx     <- params$d_x; dz <- params$d_z
  r2     <- params$r2;  sp <- params$sp
  
  mk <- function(label, timed_result) {
    M       <- timed_result$value
    elapsed <- timed_result$elapsed
    make_row(label, M, M_true, x_cols, n_val, rep_id, JOB_IDX,
             dx, dz, r2, sp, elapsed)
  }
  
  rows <- list(
    mk("ASCEND",  run_ascend( sim_obj, dx, dz, n_val)),
    mk("CBL",     run_cbl(    sim_obj, dx, dz, n_val)),
    mk("GES",     run_ges(    sim_obj, dx, dz, n_val)),
    mk("LiNGAM",  run_lingam( sim_obj, dx, dz, n_val)),
    mk("PC",      run_pc(     sim_obj, dx, dz, n_val))
  )
  rows <- Filter(Negate(is.null), rows)
  if (length(rows)==0L) return(NULL)
  do.call(rbind, rows)
}

# ── 11. Main loop ─────────────────────────────────────────────────────────────

cat(sprintf("[Job %d] Output dir : %s\n", JOB_IDX, OUT_DIR))
cat(sprintf("[Job %d] n_vec      : %s\n", JOB_IDX, paste(N_VEC, collapse=",")))
cat(sprintf("[Job %d] n_rep      : %d\n", JOB_IDX, N_REP))
cat(sprintf("[Job %d] timeout    : %ds per method\n", JOB_IDX, METHOD_TIMEOUT))

for (n_val in N_VEC) {
  cat(sprintf("\n[Job %d] === n = %d  (n/d_z = %.1f) ===\n",
              JOB_IDX, n_val, n_val/params$d_z))
  
  for (rep in seq_len(N_REP)) {
    out_file <- file.path(OUT_DIR, sprintf("n%d_rep%02d.rds", n_val, rep))
    
    # Safe restart: skip completed replicates
    if (file.exists(out_file)) {
      cat(sprintf("  rep %02d — already saved, skipping\n", rep))
      next
    }
    
    cat(sprintf("\n  rep %02d/%02d\n", rep, N_REP))
    t0 <- proc.time()["elapsed"]
    
    sim <- tryCatch(
      sim_dat(n=n_val, d_z=params$d_z, d_x=params$d_x,
              r2=params$r2, lin_pr=1,  sp=params$sp,
              p_cross=P_CROSS, x_effect=X_EFFECT,
              seed=JOB_IDX*1000000L + n_val*1000L + rep),
      error=function(e) { message("sim_dat: ", e$message); NULL }
    )
    if (is.null(sim)) { cat("    sim_dat failed — skipping rep\n"); next }
    
    result <- tryCatch(
      run_one_rep(sim, n_val, rep),
      error=function(e) { message("run_one_rep: ", e$message); NULL }
    )
    
    # Save immediately — partial results (some methods failed) are still saved
    if (!is.null(result) && nrow(result) > 0L) {
      saveRDS(result, out_file)
      cat(sprintf("  rep %02d done (%.1fs total) — %d/%d methods saved\n",
                  rep, proc.time()["elapsed"]-t0,
                  sum(result$status=="ok", na.rm=TRUE), nrow(result)))
    } else {
      cat(sprintf("  rep %02d — no results (%.1fs)\n",
                  rep, proc.time()["elapsed"]-t0))
    }
    
    rm(sim, result); gc(verbose=FALSE)
  }
}

saveRDS(list(job=JOB_IDX, params=params, completed=Sys.time()),
        file.path(OUT_DIR, "done.rds"))
cat(sprintf("\n[Job %d] Complete.\n", JOB_IDX))