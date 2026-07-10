## ======================================================================
## Synthetic benchmark for ASCEND against GENIE3, ARACNe, and WGCNA.
##
## Requires sim_dat() and ascend() to be defined in the global
## environment; source ascend.R before sourcing this file.
##
## Outputs:
##   benchmark_v2_main_raw.csv   per-replicate metrics
##   benchmark_v2_wilcoxon.csv   paired Wilcoxon tests, BH corrected
##
## Adjust N_REP to control replicate count.
##
## ----------------------------------------------------------------------
## Changes vs the previous version, to match the updated ascend.R:
##   (1) wrap_aracne() indexed sim_obj$dat with data.table syntax
##       (`[, xlabs, with = FALSE]`). The updated sim_dat() returns $dat
##       as a plain data.frame, so that line errored with
##       "unused argument (with = FALSE)" on every replicate. Now uses
##       base-R data.frame indexing (`[, xlabs, drop = FALSE]`), wrapped
##       in as.data.frame() so it is correct for both a data.frame and a
##       data.table.
##   (2) The parallel worker export list named the OLD ascend internals
##       (topo_sort_adj, ci_test_pval, ...) which no longer exist and were
##       silently dropped. The updated ascend() depends on topo_order,
##       is_dag, try_edge, ci_pval and iamb; these are now exported so the
##       %dopar% workers can actually find them.
##   (3) wrap_ascend() now calls ascend(..., verbose = FALSE); the updated
##       ascend() reports progress via message(), which capture.output()
##       does not intercept, so this keeps the benchmark log clean.
## ======================================================================

N_REP <- 50L

suppressPackageStartupMessages({
  library(data.table); library(igraph)
  library(matrixStats); library(dplyr)
  library(minet);       library(GENIE3)
  library(WGCNA);       library(bnlearn)
  library(parallel);    library(doParallel); library(foreach)
  library(PRROC)
})

N_CORES <- max(1L, min(detectCores() - 1L, 20L))

if (!exists("ascend") || !exists("sim_dat")) {
  stop("ascend or sim_dat not found. Source ascend.R first.")
}

## --- Ground truth ---------------------------------------------------------

transitive_closure <- function(A) {
  A <- as.matrix(A); A[is.na(A)] <- 0; diag(A) <- 0
  R <- (A > 0) * 1L; Ap <- R; n <- nrow(A)
  for (k in seq_len(n - 1L)) {
    Ap <- ((Ap %*% R) > 0) * 1L
    R  <- ((R  + Ap) > 0) * 1L
  }
  diag(R) <- 0L; R
}

# Returns directed adjacency, ancestral closure, undirected skeleton, and
# edge count, all with consistent x-labels.
prepare_gt <- function(sim_obj) {
  adj   <- sim_obj$adj_xx; adj[is.na(adj)] <- 0
  p     <- nrow(adj); xlabs <- paste0("x", seq_len(p))
  A_dir <- t(adj); diag(A_dir) <- 0
  rownames(A_dir) <- colnames(A_dir) <- xlabs
  A_anc  <- transitive_closure(A_dir)
  rownames(A_anc) <- colnames(A_anc) <- xlabs
  A_skel <- ((A_anc + t(A_anc)) > 0L) * 1L; diag(A_skel) <- 0L
  rownames(A_skel) <- colnames(A_skel) <- xlabs
  K <- sum(A_skel[upper.tri(A_skel)])
  list(A_dir = A_dir, A_anc = A_anc, A_skel = A_skel,
       K = K, p = p, xlabs = xlabs)
}

## --- Evaluation metrics ---------------------------------------------------

compute_auroc <- function(score_sym, A_skel) {
  S  <- pmax(as.matrix(score_sym), t(as.matrix(score_sym)))
  ut <- upper.tri(A_skel)
  sc <- as.numeric(S[ut]); lb <- as.integer(A_skel[ut])
  np <- sum(lb == 1L); nn <- sum(lb == 0L)
  if (np == 0L || nn == 0L) return(NA_real_)
  r   <- rank(sc, ties.method = "average")
  auc <- (sum(r[lb == 1L]) - np * (np + 1L) / 2.0) / (np * nn)
  max(auc, 1.0 - auc)
}

compute_aupr <- function(score_sym, A_skel) {
  S  <- pmax(as.matrix(score_sym), t(as.matrix(score_sym)))
  ut <- upper.tri(A_skel)
  sc <- as.numeric(S[ut]); lb <- as.integer(A_skel[ut])
  if (sum(lb == 1L) == 0L || sum(lb == 0L) == 0L) return(NA_real_)
  pr <- PRROC::pr.curve(scores.class0 = sc[lb == 1L],
                        scores.class1 = sc[lb == 0L],
                        curve = FALSE)
  pr$auc.integral
}

aupr_baseline <- function(A_skel) {
  ut <- upper.tri(A_skel)
  sum(A_skel[ut] == 1L) / sum(ut)
}

prf_from_binary <- function(A_pred, A_skel) {
  A_pred <- (as.matrix(A_pred) > 0) * 1L
  ut <- upper.tri(A_skel)
  tp <- sum(A_skel[ut] == 1L & A_pred[ut] == 1L)
  fp <- sum(A_skel[ut] == 0L & A_pred[ut] == 1L)
  fn <- sum(A_skel[ut] == 1L & A_pred[ut] == 0L)
  pr <- if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_
  re <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  f1 <- if (!is.na(pr) && !is.na(re) && (pr + re) > 0)
    2.0 * pr * re / (pr + re) else NA_real_
  list(precision = pr, recall = re, f1 = f1, tp = tp, fp = fp, fn = fn)
}

# Threshold a continuous weight matrix to the top K_match symmetric edges.
topK_skeleton <- function(W, K_match) {
  W  <- as.matrix(W); W[is.na(W)] <- 0; diag(W) <- 0
  Ws <- pmax(W, t(W))
  d  <- nrow(Ws); A <- matrix(0L, d, d, dimnames = dimnames(Ws))
  if (K_match <= 0L) return(A)
  ut_idx <- which(upper.tri(Ws), arr.ind = TRUE)
  if (nrow(ut_idx) == 0L) return(A)
  vals <- Ws[ut_idx]
  ord  <- order(-vals, seq_len(length(vals)))
  keep <- ut_idx[ord[seq_len(min(K_match, nrow(ut_idx)))], , drop = FALSE]
  for (r in seq_len(nrow(keep))) {
    i <- keep[r, 1]; j <- keep[r, 2]
    A[i, j] <- 1L; A[j, i] <- 1L
  }
  A
}

# ASCEND outputs values in {0, 0.5, 1, NA}. We convert these to a ranked
# score (NA -> 1, confirmed/oriented -> 2, independent -> 0) so they can
# be evaluated on the same AUPR/AUROC scale as the continuous methods.
ascend_score_matrix <- function(M, xlabs) {
  p <- length(xlabs)
  if (!is.null(rownames(M)) && all(xlabs %in% rownames(M)))
    M <- M[xlabs, xlabs, drop = FALSE]
  else dimnames(M) <- list(xlabs, xlabs)
  S <- matrix(0L, p, p, dimnames = list(xlabs, xlabs))
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      v <- M[i, j]
      score <- if (is.na(v)) 1L
      else if (v == 1 || v == 0.5) 2L
      else 0L
      S[i, j] <- score; S[j, i] <- score
    }
  }
  S
}

ascend_binary_skeleton <- function(M, xlabs) {
  p <- length(xlabs)
  if (!is.null(rownames(M)) && all(xlabs %in% rownames(M)))
    M <- M[xlabs, xlabs, drop = FALSE]
  else dimnames(M) <- list(xlabs, xlabs)
  A <- matrix(0L, p, p, dimnames = list(xlabs, xlabs))
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    v <- M[i, j]
    if (!is.na(v) && (v == 1 || v == 0.5)) { A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

ascend_coverage <- function(M, xlabs) {
  p <- length(xlabs)
  if (!is.null(rownames(M)) && all(xlabs %in% rownames(M)))
    M <- M[xlabs, xlabs, drop = FALSE]
  else dimnames(M) <- list(xlabs, xlabs)
  n_tot <- p * (p - 1L) / 2L; n_res <- 0L
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p)
    if (!is.na(M[i, j])) n_res <- n_res + 1L
  n_res / n_tot
}

# Of the claimed ancestral edges, what fraction agree with the true
# ancestral closure (in direction i -> j)?
ascend_direction_acc <- function(M, xlabs, A_anc) {
  p <- length(xlabs)
  if (!is.null(rownames(M)) && all(xlabs %in% rownames(M)))
    M <- M[xlabs, xlabs, drop = FALSE]
  else dimnames(M) <- list(xlabs, xlabs)
  nc <- 0L; nt <- 0L
  for (i in seq_len(p - 1L)) for (j in (i + 1L):p) {
    v <- M[i, j]
    if (is.na(v) || !(v == 1 || v == 0.5)) next
    nt <- nt + 1L
    if (A_anc[i, j] == 1L) nc <- nc + 1L
  }
  if (nt == 0L) return(NA_real_)
  nc / nt
}

## --- Method wrappers ------------------------------------------------------

wrap_ascend <- function(sim_obj) {
  tryCatch({
    invisible(capture.output(
      M <- suppressMessages(ascend(sim_obj, verbose = FALSE)),
      type = "output"))
    M
  }, error = function(e) {
    p <- sim_obj$params$d_x; labs <- paste0("x", seq_len(p))
    matrix(NA_real_, p, p, dimnames = list(labs, labs))
  })
}

wrap_genie3 <- function(sim_obj, nTrees = 500L) {
  all_labs <- colnames(sim_obj$dat)
  xlabs    <- paste0("x", seq_len(sim_obj$params$d_x))
  Xg       <- t(as.matrix(sim_obj$dat)); rownames(Xg) <- all_labs
  tryCatch({
    W <- as.matrix(GENIE3::GENIE3(Xg, nTrees = nTrees, nCores = 1L))
    W <- W[xlabs, xlabs]; W <- pmax(W, t(W)); diag(W) <- 0; W
  }, error = function(e) {
    p <- length(xlabs); matrix(0, p, p, dimnames = list(xlabs, xlabs))
  })
}

# ARACNe via minet, with fallbacks for degenerate MI estimates: try
# equal-frequency binning, then equal-width binning, and finally fall back
# to absolute Spearman correlation if both fail.
wrap_aracne <- function(sim_obj) {
  all_labs <- colnames(sim_obj$dat)
  xlabs    <- paste0("x", seq_len(sim_obj$params$d_x)); p <- length(xlabs)
  Xg_all   <- t(as.matrix(sim_obj$dat)); rownames(Xg_all) <- all_labs
  # sim_dat() returns $dat as a data.frame; as.data.frame() keeps this
  # correct whether $dat is a data.frame or a data.table.
  X_only   <- as.matrix(as.data.frame(sim_obj$dat)[, xlabs, drop = FALSE])
  n_samp   <- nrow(sim_obj$dat)
  n_bins   <- max(3L, min(10L, round(sqrt(n_samp) / 2)))
  
  extract_xx <- function(W_full) {
    ok <- intersect(xlabs, rownames(W_full))
    W  <- matrix(0, p, p, dimnames = list(xlabs, xlabs))
    if (length(ok)) W[ok, ok] <- W_full[ok, ok]
    W[is.na(W)] <- 0; diag(W) <- 0; pmax(W, t(W))
  }
  is_degenerate <- function(W) {
    if (is.null(W)) return(TRUE)
    W <- as.matrix(W); ut <- W[upper.tri(W)]
    if (max(abs(ut), na.rm = TRUE) < 1e-10) return(TRUE)
    length(unique(round(ut, 8))) < 3L
  }
  
  W_full <- tryCatch({
    mim <- minet::build.mim(Xg_all, estimator = "mi.mm",
                            disc = "equalfreq", nbins = n_bins)
    W <- as.matrix(minet::aracne(mim, eps = 0))
    if (is_degenerate(W)) as.matrix(mim) else W
  }, error = function(e) NULL)
  
  if (is.null(W_full) || is_degenerate(W_full)) {
    W_full <- tryCatch({
      mim <- minet::build.mim(Xg_all, estimator = "mi.mm",
                              disc = "equalwidth", nbins = n_bins)
      W <- as.matrix(minet::aracne(mim, eps = 0))
      if (is_degenerate(W)) as.matrix(mim) else W
    }, error = function(e) NULL)
  }
  
  if (is.null(W_full) || is_degenerate(W_full)) {
    W_xx <- abs(cor(X_only, method = "spearman",
                    use = "pairwise.complete.obs"))
    W_xx[is.na(W_xx)] <- 0; diag(W_xx) <- 0
    rownames(W_xx) <- colnames(W_xx) <- xlabs
    return(W_xx)
  }
  extract_xx(W_full)
}

# WGCNA with a data-driven choice of soft-threshold power (R^2 >= 0.80
# scale-free fit), falling back to power 6 if the criterion is not met.
wrap_wgcna <- function(sim_obj) {
  all_labs <- colnames(sim_obj$dat)
  xlabs    <- paste0("x", seq_len(sim_obj$params$d_x)); p <- length(xlabs)
  X_all    <- as.matrix(sim_obj$dat)
  tryCatch({
    powers <- c(1:10, seq(12, 20, 2))
    capture.output(
      sft <- WGCNA::pickSoftThreshold(
        X_all, powerVector = powers, networkType = "signed hybrid",
        RsquaredCut = 0.80, verbose = 0
      ), type = "output"
    )
    beta <- sft$powerEstimate
    if (is.na(beta) || length(beta) == 0L)
      beta <- powers[which.max(sft$fitIndices[, "SFT.R.sq"])]
    if (is.na(beta) || !is.finite(beta)) beta <- 6L
    adj_w <- suppressMessages(suppressWarnings(
      WGCNA::adjacency(X_all, power = beta, type = "signed hybrid",
                       corFnc = "cor",
                       corOptions = list(use = "pairwise.complete.obs"))
    ))
    TOM_full <- suppressMessages(suppressWarnings(
      WGCNA::TOMsimilarity(adj_w, TOMType = "signed", verbose = 0)
    ))
    rownames(TOM_full) <- colnames(TOM_full) <- all_labs
    W <- TOM_full[xlabs, xlabs, drop = FALSE]
    diag(W) <- 0; W
  }, error = function(e) {
    matrix(0, p, p, dimnames = list(xlabs, xlabs))
  })
}

## --- Single replicate -----------------------------------------------------

# Runs all four methods on a single simulated dataset and returns a long
# data.table of metrics. F1, precision and recall are computed at matched
# edge count: competitor weights are thresholded to retain the same number
# of edges that ASCEND claims, so methods are compared at the same operating
# point.
run_one_rep <- function(sim_obj, rep_id, n_val, sp_val, r2_val) {
  gt    <- prepare_gt(sim_obj)
  xlabs <- gt$xlabs; A_skel <- gt$A_skel; A_anc <- gt$A_anc
  
  t0 <- Sys.time(); M_asc <- wrap_ascend(sim_obj)
  t_asc <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time(); W_gen <- wrap_genie3(sim_obj)
  t_gen <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time(); W_ara <- wrap_aracne(sim_obj)
  t_ara <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time(); W_wgc <- wrap_wgcna(sim_obj)
  t_wgc <- as.numeric(Sys.time() - t0, units = "secs")
  
  sc_asc <- ascend_score_matrix(M_asc, xlabs)
  sc_gen <- { W <- pmax(W_gen, t(W_gen)); diag(W) <- 0; W }
  sc_ara <- { W <- pmax(W_ara, t(W_ara)); diag(W) <- 0; W }
  sc_wgc <- { W <- pmax(W_wgc, t(W_wgc)); diag(W) <- 0; W }
  
  au_asc <- compute_auroc(sc_asc, A_skel); ap_asc <- compute_aupr(sc_asc, A_skel)
  au_gen <- compute_auroc(sc_gen, A_skel); ap_gen <- compute_aupr(sc_gen, A_skel)
  au_ara <- compute_auroc(sc_ara, A_skel); ap_ara <- compute_aupr(sc_ara, A_skel)
  au_wgc <- compute_auroc(sc_wgc, A_skel); ap_wgc <- compute_aupr(sc_wgc, A_skel)
  ap_base <- aupr_baseline(A_skel)
  
  bin_asc <- ascend_binary_skeleton(M_asc, xlabs)
  K_match <- sum(bin_asc[upper.tri(bin_asc)])
  
  bin_gen <- topK_skeleton(sc_gen, K_match)
  bin_ara <- topK_skeleton(sc_ara, K_match)
  bin_wgc <- topK_skeleton(sc_wgc, K_match)
  
  m_asc <- prf_from_binary(bin_asc, A_skel)
  m_gen <- prf_from_binary(bin_gen, A_skel)
  m_ara <- prf_from_binary(bin_ara, A_skel)
  m_wgc <- prf_from_binary(bin_wgc, A_skel)
  
  cov_asc <- ascend_coverage(M_asc, xlabs)
  dir_asc <- ascend_direction_acc(M_asc, xlabs, A_anc)
  
  v <- function(x) if (is.na(x)) "  NA" else sprintf("%.2f", x)
  cat(sprintf("rep %2d | n=%d sp=%.1f r2=%.1f K_true=%d K_match=%d\n",
              rep_id, n_val, sp_val, r2_val, gt$K, K_match))
  cat(sprintf("  ASCEND  AUPR=%s AUROC=%s F1=%s P=%s R=%s dir=%s t=%.1fs\n",
              v(ap_asc), v(au_asc), v(m_asc$f1), v(m_asc$precision),
              v(m_asc$recall), v(dir_asc), t_asc))
  cat(sprintf("  GENIE3  AUPR=%s AUROC=%s F1=%s P=%s R=%s         t=%.1fs\n",
              v(ap_gen), v(au_gen), v(m_gen$f1), v(m_gen$precision),
              v(m_gen$recall), t_gen))
  cat(sprintf("  ARACNE  AUPR=%s AUROC=%s F1=%s P=%s R=%s         t=%.1fs\n",
              v(ap_ara), v(au_ara), v(m_ara$f1), v(m_ara$precision),
              v(m_ara$recall), t_ara))
  cat(sprintf("  WGCNA   AUPR=%s AUROC=%s F1=%s P=%s R=%s         t=%.1fs\n\n",
              v(ap_wgc), v(au_wgc), v(m_wgc$f1), v(m_wgc$precision),
              v(m_wgc$recall), t_wgc))
  
  data.table(
    method        = c("ASCEND", "GENIE3", "ARACNE", "WGCNA"),
    n             = n_val,
    sp            = sp_val,
    r2            = r2_val,
    rep           = rep_id,
    K_true        = gt$K,
    K_match       = K_match,
    aupr_base     = ap_base,
    coverage      = c(cov_asc, 1.0, 1.0, 1.0),
    aupr          = c(ap_asc, ap_gen, ap_ara, ap_wgc),
    auroc         = c(au_asc, au_gen, au_ara, au_wgc),
    f1            = c(m_asc$f1, m_gen$f1, m_ara$f1, m_wgc$f1),
    precision     = c(m_asc$precision, m_gen$precision, m_ara$precision, m_wgc$precision),
    recall        = c(m_asc$recall, m_gen$recall, m_ara$recall, m_wgc$recall),
    direction_acc = c(dir_asc, NA_real_, NA_real_, NA_real_),
    runtime_s     = c(t_asc, t_gen, t_ara, t_wgc)
  )
}

## --- Runner ---------------------------------------------------------------

# Iterates the parameter grid; for each cell runs n_rep replicates in
# parallel (or sequentially if parallel=FALSE) and collects results.
# Failed replicates are logged but do not halt the run.
run_full <- function(conditions, n_rep, d_z = 20L, d_x = 15L,
                     p_cross = 0.20, x_effect = 0.8, seed_base = 42L,
                     label = "run", parallel = TRUE) {
  if (parallel) {
    registerDoParallel(N_CORES)
    on.exit(stopImplicitCluster())
  }
  
  # Functions the workers need in scope. These are the CURRENT ascend.R
  # internals (topo_order/is_dag/try_edge/ci_pval/iamb) plus the benchmark
  # helpers; only those that actually exist are exported.
  worker_fns <- c(
    "run_one_rep", "prepare_gt", "transitive_closure",
    "compute_auroc", "compute_aupr", "aupr_baseline",
    "prf_from_binary", "topK_skeleton",
    "ascend_score_matrix", "ascend_binary_skeleton",
    "ascend_coverage", "ascend_direction_acc",
    "wrap_ascend", "wrap_genie3", "wrap_aracne", "wrap_wgcna",
    "ascend", "sim_dat",
    "topo_order", "is_dag", "try_edge", "ci_pval", "iamb",
    "true_ancestral"
  )
  worker_fns <- worker_fns[sapply(worker_fns, exists)]
  
  all_res <- vector("list", nrow(conditions))
  for (ci in seq_len(nrow(conditions))) {
    n_val  <- conditions$n[ci]
    sp_val <- conditions$sp[ci]
    r2_val <- conditions$r2[ci]
    cat(sprintf("\n[%s] n=%d sp=%.1f r2=%.1f -- %d reps on %d cores\n",
                label, n_val, sp_val, r2_val, n_rep,
                if (parallel) N_CORES else 1L))
    
    if (parallel) {
      rep_list <- foreach(
        r              = seq_len(n_rep),
        .packages      = c("data.table", "minet", "GENIE3", "WGCNA",
                           "bnlearn", "igraph", "matrixStats", "PRROC", "dplyr"),
        .export        = worker_fns,
        .errorhandling = "pass"
      ) %dopar% {
        set.seed(seed_base + n_val * 1000L + round(sp_val * 100L) * 10L +
                   round(r2_val * 100L) * 100L + r)
        tryCatch({
          sim <- sim_dat(n = n_val, d_z = d_z, d_x = d_x, r2 = r2_val,
                         lin_pr = 1, sp = sp_val,
                         p_cross = p_cross, x_effect = x_effect)
          run_one_rep(sim, r, n_val, sp_val, r2_val)
        }, error = function(e) {
          list(error = conditionMessage(e), rep = r)
        })
      }
    } else {
      rep_list <- vector("list", n_rep)
      for (r in seq_len(n_rep)) {
        set.seed(seed_base + n_val * 1000L + round(sp_val * 100L) * 10L +
                   round(r2_val * 100L) * 100L + r)
        rep_list[[r]] <- tryCatch({
          sim <- sim_dat(n = n_val, d_z = d_z, d_x = d_x, r2 = r2_val,
                         lin_pr = 1, sp = sp_val,
                         p_cross = p_cross, x_effect = x_effect)
          run_one_rep(sim, r, n_val, sp_val, r2_val)
        }, error = function(e) {
          cat(sprintf("  rep %d failed: %s\n", r, conditionMessage(e)))
          list(error = conditionMessage(e), rep = r)
        })
      }
    }
    
    good_reps <- list(); bad_msgs <- character(0)
    for (rr in seq_along(rep_list)) {
      x <- rep_list[[rr]]
      if (inherits(x, "data.table") && nrow(x) > 0L) {
        good_reps[[length(good_reps) + 1L]] <- x
      } else if (inherits(x, "error") || inherits(x, "simpleError")) {
        bad_msgs <- c(bad_msgs, conditionMessage(x))
      } else if (is.list(x) && !is.null(x$error)) {
        bad_msgs <- c(bad_msgs, as.character(x$error))
      } else {
        bad_msgs <- c(bad_msgs, sprintf("unexpected class: %s",
                                        paste(class(x), collapse = ",")))
      }
    }
    n_ok <- length(good_reps); n_err <- length(bad_msgs)
    cat(sprintf("  condition done: %d/%d reps succeeded, %d errored\n",
                n_ok, n_rep, n_err))
    if (n_err > 0L) {
      for (m in unique(bad_msgs)[1:min(3L, length(unique(bad_msgs)))])
        cat(sprintf("    error: %s\n", m))
    }
    if (n_ok > 0L) {
      all_res[[ci]] <- rbindlist(good_reps, use.names = TRUE, fill = TRUE)
    }
  }
  rbindlist(Filter(Negate(is.null), all_res), use.names = TRUE, fill = TRUE)
}

## --- Significance testing -------------------------------------------------

# Paired one-sided Wilcoxon signed-rank tests of ASCEND vs each competitor,
# with Benjamini-Hochberg correction across the family.
wilcoxon_tests <- function(dt, metric_name) {
  conds <- unique(dt[, .(n, sp, r2)])
  results <- list()
  for (ci in seq_len(nrow(conds))) {
    sub <- dt[n == conds$n[ci] & sp == conds$sp[ci] & r2 == conds$r2[ci]]
    for (other in c("GENIE3", "ARACNE", "WGCNA")) {
      asc_o <- sub[method == "ASCEND"][order(rep), get(metric_name)]
      oth   <- sub[method == other][order(rep), get(metric_name)]
      n_paired <- sum(!is.na(asc_o) & !is.na(oth))
      if (n_paired < 5L) next
      tst <- tryCatch(
        wilcox.test(asc_o, oth, paired = TRUE, alternative = "greater",
                    exact = FALSE),
        error = function(e) NULL
      )
      if (is.null(tst)) next
      results[[length(results) + 1L]] <- data.table(
        n = conds$n[ci], sp = conds$sp[ci], r2 = conds$r2[ci],
        metric = metric_name, comparison = paste0("ASCEND>", other),
        n_paired = n_paired,
        median_diff = median(asc_o - oth, na.rm = TRUE),
        W = as.numeric(tst$statistic),
        p_raw = tst$p.value
      )
    }
  }
  res <- rbindlist(results)
  if (nrow(res) == 0L) return(res)
  res[, p_bh := p.adjust(p_raw, method = "BH")]
  res[, signif := ifelse(p_bh < 0.001, "***",
                         ifelse(p_bh < 0.01,  "**",
                                ifelse(p_bh < 0.05,  "*", "ns")))]
  res
}

## --- Summary tables -------------------------------------------------------

METHOD_ORDER <- c("ASCEND", "GENIE3", "ARACNE", "WGCNA")

make_summary <- function(dt, n_sel, sp_sel, r2_sel = 0.5) {
  dt %>%
    filter(n == n_sel, sp == sp_sel, r2 == r2_sel) %>%
    group_by(method) %>%
    summarise(
      nrep    = n(),
      AUPR_m  = mean(aupr, na.rm = TRUE),
      AUPR_se = sd(aupr,  na.rm = TRUE) / sqrt(n()),
      AU_m    = mean(auroc, na.rm = TRUE),
      AU_se   = sd(auroc,  na.rm = TRUE) / sqrt(n()),
      F1_m    = mean(f1, na.rm = TRUE),
      F1_se   = sd(f1,  na.rm = TRUE) / sqrt(n()),
      PR_m    = mean(precision, na.rm = TRUE),
      RE_m    = mean(recall,    na.rm = TRUE),
      DIR_m   = mean(direction_acc, na.rm = TRUE),
      DIR_se  = sd(direction_acc,   na.rm = TRUE) / sqrt(n()),
      RT_m    = mean(runtime_s, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(match(method, METHOD_ORDER))
}

print_table <- function(summ, n_sel, sp_sel, r2_sel = 0.5, base_aupr = NA) {
  cat(sprintf("\nn=%d  sp=%.1f  r2=%.1f  (%d reps)\n",
              n_sel, sp_sel, r2_sel,
              if (nrow(summ) > 0) summ$nrep[1] else 0L))
  if (nrow(summ) == 0L) {
    cat("  no results for this cell\n"); return(invisible(NULL))
  }
  if (!is.na(base_aupr))
    cat(sprintf("  random-baseline AUPR = %.3f (edge prevalence)\n", base_aupr))
  cat(sprintf("  %-8s %18s %14s %18s %10s %10s\n",
              "Method", "AUPR primary", "AUROC", "F1@matched-K",
              "Precision", "Recall"))
  for (i in seq_len(nrow(summ))) {
    r <- summ[i, ]
    cat(sprintf("  %-8s  %6.3f +/- %.3f    %5.3f +/- %.3f    %6.3f +/- %.3f    %5.3f      %5.3f\n",
                r$method, r$AUPR_m, r$AUPR_se,
                r$AU_m, r$AU_se, r$F1_m, r$F1_se, r$PR_m, r$RE_m))
  }
  asc <- summ %>% filter(method == "ASCEND")
  if (nrow(asc) > 0 && !is.na(asc$DIR_m))
    cat(sprintf("  ASCEND direction accuracy: %.3f +/- %.3f (chance = 0.5)\n",
                asc$DIR_m, asc$DIR_se))
  cat(sprintf("  Runtimes (s): %s\n",
              paste(sprintf("%s=%.1f", summ$method, summ$RT_m), collapse = ", ")))
  cat("\n")
}

print_all_cells <- function(dt) {
  if (nrow(dt) == 0L) {
    cat("no data to print\n"); return(invisible(NULL))
  }
  conds <- unique(dt[, .(n, sp, r2)])
  setorder(conds, sp, n, r2)
  for (i in seq_len(nrow(conds))) {
    base_aupr <- mean(dt[n == conds$n[i] & sp == conds$sp[i] &
                           r2 == conds$r2[i] & method == "ASCEND", aupr_base],
                      na.rm = TRUE)
    s <- make_summary(dt, conds$n[i], conds$sp[i], conds$r2[i])
    print_table(s, conds$n[i], conds$sp[i], conds$r2[i], base_aupr = base_aupr)
  }
}

## --- Run ------------------------------------------------------------------

cat(sprintf("\nReps per condition: %d   Cores: %d\n\n", N_REP, N_CORES))

# Sequential dry run to confirm wrappers, oracle and metrics all work
# before launching the parallel main grid.
cat("Step 1: sanity check (sequential, 1 cell, 2 reps)\n")
sanity_conds <- data.table(n = 500L, sp = 0.9, r2 = 0.5)
sanity_res <- run_full(
  conditions = sanity_conds,
  n_rep      = 2L,
  d_z = 20L, d_x = 15L,
  p_cross = 0.20, x_effect = 0.8,
  seed_base = 42L,
  label     = "sanity",
  parallel  = FALSE
)
if (nrow(sanity_res) == 0L) {
  stop("Sanity check returned no rows; halting before main run.")
}

# Main parameter grid. Edit as needed.
cat("\nStep 2: main grid\n")
main_conds <- CJ(n  = c(1000L, 2000L),
                 sp = c(0.5, 0.7, 0.9),
                 r2 = c(0.5, 0.7))
main_res <- run_full(
  conditions = main_conds,
  n_rep      = N_REP,
  d_z = 20L, d_x = 15L,
  p_cross = 0.20, x_effect = 0.8,
  seed_base = 42L,
  label     = "main",
  parallel  = TRUE
)

fwrite(main_res, "benchmark_v2_main_raw.csv")
print_all_cells(main_res)

wilcox_all <- rbindlist(lapply(
  c("aupr", "auroc", "f1"),
  function(mn) wilcoxon_tests(main_res, mn)
), use.names = TRUE, fill = TRUE)
fwrite(wilcox_all, "benchmark_v2_wilcoxon.csv")
cat("\nWrote benchmark_v2_main_raw.csv and benchmark_v2_wilcoxon.csv\n")