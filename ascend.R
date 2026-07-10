# ======================================================================
# ASCEND : Ancestral Scalable Causal discovEry via iNherited Descent
# ----------------------------------------------------------------------
# Constraint-based ancestral discovery for two-tier systems, where a set
# of background variables Z is known to causally precede a set of
# foreground variables X. ASCEND recovers the ancestral order among the
# foreground variables by testing conditional independencies while
# conditioning only on small, dynamically maintained sets of nearest
# ancestors.
#
# Output: an ancestrality matrix M with M[a, b] in {1, 0.5, 0, NA}
#   1   a is an ancestor of b            (a < b)
#   0.5 a is not a descendant of b       (a <= b)
#   0   neither is an ancestor of the other (a ~ b)
#   NA  undetermined
#
# Dependencies: base R + stats only.
# ======================================================================


# ----------------------------------------------------------------------
# 1. Two-tier linear-Gaussian simulator
#    Returns $dat (z1..z_dz, x1..x_dx) and $adj_xx with adj_xx[child, parent] = 1.
# ----------------------------------------------------------------------

sim_dat <- function(n, d_z, d_x,
                    r2       = 0.5,   # signal-to-noise (variance explained by parents)
                    lin_pr   = 1,     # fraction of linear parents (1 = fully linear)
                    sp       = 0.9,   # sparsity (larger => sparser graph)
                    p_cross  = 0.05,  # P(Z -> X edge)
                    x_effect = 0.8,   # X -> X effect size
                    seed     = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # optionally transform a fraction of parent columns to inject nonlinearity
  prep <- function(mat, pr) {
    mat <- as.matrix(mat); p <- ncol(mat)
    if (p == 0 || pr >= 1) return(mat)
    n_nl <- round((1 - pr) * p); if (n_nl <= 0) return(mat)
    cols  <- sample.int(p, n_nl)
    types <- sample(c("sq", "sqrt", "softplus", "relu"), n_nl, replace = TRUE)
    for (ii in seq_along(cols)) {
      v <- mat[, cols[ii]]
      mat[, cols[ii]] <- switch(types[ii],
                                sq = v^2, sqrt = sqrt(abs(v)), softplus = log1p(exp(v)), relu = pmax(v, 0))
    }
    mat
  }
  
  noise_sd <- function(signal, r2) {
    v <- var(signal)
    if (is.na(v) || v == 0) return(1)
    if (r2 <= 0) return(sqrt(v))
    sqrt(max(0, (v - r2 * v) / r2))
  }
  
  rand_dag <- function(k) {                       # random DAG (upper-triangular, permuted)
    p_edge <- min(max(1 - sp, 0.01), 0.4)
    A <- matrix(0, k, k)
    if (k > 1) for (i in 1:(k - 1)) for (j in (i + 1):k)
      if (runif(1) < p_edge) A[i, j] <- 1
    perm <- sample.int(k); A[perm, perm]
  }
  
  Azz <- rand_dag(d_z); Axx <- rand_dag(d_x)
  Azx <- matrix(rbinom(d_z * d_x, 1, p_cross), d_z, d_x)
  colSds <- function(M) apply(M, 2, sd)
  
  # sample background Z in topological order
  Z <- matrix(0, n, d_z)
  for (idx in topo_order(Azz)) {
    par <- which(Azz[, idx] == 1)
    Z[, idx] <- if (!length(par)) rnorm(n)
    else as.numeric(prep(Z[, par, drop = FALSE], lin_pr) %*% rnorm(length(par))) + rnorm(n)
  }
  
  # sample foreground X in topological order; record the true X -> X edges
  X <- matrix(0, n, d_x)
  adj_xx <- matrix(0, d_x, d_x); diag(adj_xx) <- NA_real_
  ord_x <- topo_order(Axx)
  for (s in seq_along(ord_x)) {
    j  <- ord_x[s]
    zp <- which(Azx[, j] == 1)
    sig_z <- if (length(zp))
      as.numeric(prep(Z[, zp, drop = FALSE], lin_pr) %*% sample(c(1, -1), length(zp), TRUE)) else rep(0, n)
    xp <- intersect(which(Axx[, j] == 1), if (s > 1) ord_x[1:(s - 1)] else integer(0))
    sig_x <- rep(0, n)
    if (length(xp)) {
      Pax <- prep(X[, xp, drop = FALSE], lin_pr)
      adj_xx[j, xp] <- 1
      csd <- colSds(Pax); csd[csd == 0] <- 1
      sig_x <- as.numeric(Pax %*% ((x_effect / csd) * sample(c(1, -1), length(xp), TRUE)))
    }
    sig <- sig_z + sig_x
    if (var(sig) < 1e-6) sig <- rnorm(n)
    X[, j] <- sig + rnorm(n, sd = noise_sd(sig, r2))
  }
  
  dat <- as.data.frame(cbind(Z, X))
  colnames(dat) <- c(paste0("z", seq_len(d_z)), paste0("x", seq_len(d_x)))
  rownames(adj_xx) <- colnames(adj_xx) <- paste0("x", seq_len(d_x))
  list(dat = dat, adj_xx = adj_xx,
       params = list(n = n, d_z = d_z, d_x = d_x, r2 = r2, sp = sp))
}


# ----------------------------------------------------------------------
# 2. Graph utilities (adjacency convention: A[i, j] = 1 means i -> j)
# ----------------------------------------------------------------------

# topological order (sources first); NULL if the graph has a cycle
topo_order <- function(A) {
  A <- (A != 0) * 1; A[is.na(A)] <- 0
  k <- nrow(A); if (k == 0) return(integer(0))
  indeg <- colSums(A); ord <- integer(0)
  zeros <- which(indeg == 0)
  while (length(zeros)) {
    v <- zeros[1]; zeros <- zeros[-1]; ord <- c(ord, v)
    nb <- which(A[v, ] != 0)
    for (u in nb) { indeg[u] <- indeg[u] - 1; if (indeg[u] == 0) zeros <- c(zeros, u) }
    A[v, nb] <- 0
  }
  if (length(ord) != k) NULL else ord
}

is_dag_custom <- function(A) !is.null(topo_order(A))

# write M[i,j]=v_ij, M[j,i]=v_ji only if the precedence graph stays acyclic
try_edge <- function(M, i, j, v_ij, v_ji) {
  T <- M; T[i, j] <- v_ij; T[j, i] <- v_ji
  b <- ((T == 1 | T == 0.5) * 1); b[is.na(b)] <- 0
  if (is_dag_custom(b)) list(M = T, ok = TRUE) else list(M = M, ok = FALSE)
}


# ----------------------------------------------------------------------
# 3. Conditional independence test
#    Fisher's z on the partial correlation rho(x, y | S), read from the
#    precomputed correlation matrix C. Returns a p-value (small => dependent).
# ----------------------------------------------------------------------

ci_pval <- function(x, y, S, C, n) {
  S  <- S[S != x & S != y]
  df <- n - length(S) - 3
  if (df < 1) return(NA_real_)
  m <- C[c(x, y, S), c(x, y, S), drop = FALSE]
  P <- tryCatch(solve(m),
                error = function(e) tryCatch(solve(m + diag(1e-8, nrow(m))),
                                             error = function(e2) NULL))
  if (is.null(P)) return(NA_real_)
  den <- sqrt(P[1, 1] * P[2, 2])
  if (!is.finite(den) || den <= 0) return(NA_real_)
  rho <- max(min(-P[1, 2] / den, 1 - 1e-10), -1 + 1e-10)
  2 * pnorm(-abs(sqrt(df) * atanh(rho)))
}


# ----------------------------------------------------------------------
# 4. Markov-blanket learner (IAMB) restricted to a non-descendant pool.
#    Returns the nearest ancestors pa(target; pool).
# ----------------------------------------------------------------------

iamb <- function(target, pool, C, n, alpha, max_cond) {
  pool <- setdiff(pool, target)
  mb <- character(0)
  repeat {                                        # grow: add the strongest associate
    rest <- setdiff(pool, mb)
    if (!length(rest) || length(mb) >= max_cond) break
    pv <- vapply(rest, function(w) ci_pval(target, w, mb, C, n), numeric(1))
    pv[is.na(pv)] <- 1
    k <- which.min(pv)
    if (pv[k] <= alpha) mb <- c(mb, rest[k]) else break
  }
  changed <- TRUE                                 # shrink: drop now-redundant members
  while (changed && length(mb)) {
    changed <- FALSE
    for (w in mb) {
      p <- ci_pval(target, w, setdiff(mb, w), C, n)
      if (!is.na(p) && p > alpha) { mb <- setdiff(mb, w); changed <- TRUE; break }
    }
  }
  mb
}


# ----------------------------------------------------------------------
# 5. ASCEND
#
# Arguments
#   sim_obj    output of sim_dat() (or any list with $dat holding z*/x* columns)
#   maxiter    maximum outer iterations
#   alpha      CI threshold for the orientation rules R1/R2/R3
#   alpha_mb   CI threshold for Markov-blanket discovery
#   fdr        Benjamini-Hochberg correction on the pairwise R3 tests per sweep
#   min_votes  witnesses that must agree before committing an orientation
#              (1 = first valid witness; raise to >=2 to trade recall for precision)
#   prescreen  keep Z marginally associated to X (BH p < prescreen) as the
#              IAMB candidate pool; bounds cost when d_z is large
# ----------------------------------------------------------------------

ascend <- function(sim_obj,
                   maxiter   = 10,
                   alpha     = 0.05,
                   alpha_mb  = 0.05,
                   fdr       = TRUE,
                   min_votes = 1,
                   prescreen = 0.30,
                   verbose   = TRUE) {
  
  dat   <- as.data.frame(sim_obj$dat)
  xlabs <- grep("^x", colnames(dat), value = TRUE)
  zlabs <- grep("^z", colnames(dat), value = TRUE)
  d_x   <- length(xlabs)
  n     <- nrow(dat)
  is_z  <- setNames(colnames(dat) %in% zlabs, colnames(dat))
  max_cond <- max(1, floor((n - 4) / 3))
  
  # one correlation matrix powers every CI test
  Xall <- as.matrix(dat); Xall[!is.finite(Xall)] <- NA
  C <- suppressWarnings(cor(Xall, use = "pairwise.complete.obs"))
  C[is.na(C)] <- 0; diag(C) <- 1
  
  # per-X background pool: Z marginally associated to X
  prescreen_z <- function(xi) {
    pv <- vapply(zlabs, function(zc) ci_pval(xi, zc, character(0), C, n), numeric(1))
    pv[is.na(pv)] <- 1
    keep <- which(p.adjust(pv, "BH") < prescreen)
    if (length(keep) < min(3, length(zlabs))) keep <- order(pv)[seq_len(min(3, length(pv)))]
    zlabs[keep]
  }
  z_pool <- setNames(lapply(xlabs, prescreen_z), xlabs)
  
  nearest_anc <- function(xi, pool) iamb(xi, intersect(pool, colnames(dat)), C, n, alpha_mb, max_cond)
  
  pa <- setNames(lapply(xlabs, function(xi) nearest_anc(xi, z_pool[[xi]])), xlabs)  # t = 1: over Z
  
  M  <- matrix(NA_real_, d_x, d_x, dimnames = list(xlabs, xlabs))
  
  # Guarded conditioning set: union of nearest ancestors, keeping a foreground
  # W only if it is a known non-descendant of BOTH endpoints. Background Z always
  # qualifies. This excludes mediators of either endpoint, which is what keeps R3 sound.
  build_S <- function(i, j) {
    cand <- setdiff(union(pa[[xlabs[i]]], pa[[xlabs[j]]]), xlabs[c(i, j)])
    keep <- vapply(cand, function(w) {
      if (is_z[w]) return(TRUE)
      wi <- match(w, xlabs)
      (!is.na(M[wi, i]) && M[wi, i] %in% c(0.5, 1)) &&
        (!is.na(M[wi, j]) && M[wi, j] %in% c(0.5, 1))
    }, logical(1))
    cand[keep]
  }
  
  converged <- FALSE; iter <- 0
  while (!converged && iter < maxiter) {
    iter <- iter + 1; converged <- TRUE
    if (verbose) message(sprintf("  [ascend] iter %d  |  %d pairs unresolved",
                                 iter, sum(is.na(M[lower.tri(M)]))))
    
    # --- R3 over all unresolved pairs, with one BH correction per sweep ---
    info <- list()
    for (i in 2:d_x) for (j in 1:(i - 1)) {
      if (!is.na(M[i, j])) next
      S  <- build_S(i, j)
      pv <- ci_pval(xlabs[i], xlabs[j], S, C, n)
      if (is.na(pv)) { converged <- FALSE; next }
      info[[length(info) + 1]] <- list(i = i, j = j, S = S, pv = pv)
    }
    if (length(info)) {
      adj <- if (fdr) p.adjust(vapply(info, function(e) e$pv, numeric(1)), "BH")
      else vapply(info, function(e) e$pv, numeric(1))
      
      for (idx in seq_along(info)) {
        e <- info[[idx]]; i <- e$i; j <- e$j; S <- e$S
        xi <- xlabs[i]; xj <- xlabs[j]
        
        if (adj[idx] > alpha) { M[i, j] <- 0; M[j, i] <- 0; converged <- FALSE; next }  # R3: ~
        if (!length(S)) next
        
        # --- R1 / R2 orientation, tallying votes across all witnesses W ---
        votes <- c(r1 = 0L, r1r = 0L, r2 = 0L, r2r = 0L)
        evid  <- c(r1 = 0,  r1r = 0,  r2 = 0,  r2r = 0)
        for (W in S) {
          Sw  <- setdiff(S, W)
          pj0 <- ci_pval(W, xj, Sw, C, n);  pjI <- ci_pval(W, xj, c(Sw, xi), C, n)
          pi0 <- ci_pval(W, xi, Sw, C, n);  piJ <- ci_pval(W, xi, c(Sw, xj), C, n)
          if (!is.na(pj0) && !is.na(pjI) && pj0 <= alpha && pjI > alpha) {        # R1  => Xi < Xj
            votes["r1"]  <- votes["r1"]  + 1L; evid["r1"]  <- evid["r1"]  - log(pj0 + 1e-300) }
          if (!is.na(pi0) && !is.na(piJ) && pi0 <= alpha && piJ > alpha) {        # R1r => Xj < Xi
            votes["r1r"] <- votes["r1r"] + 1L; evid["r1r"] <- evid["r1r"] - log(pi0 + 1e-300) }
          if (!is.na(pi0) && !is.na(piJ) && pi0 > alpha && piJ <= alpha) {        # R2  => Xi <= Xj
            votes["r2"]  <- votes["r2"]  + 1L; evid["r2"]  <- evid["r2"]  - log(piJ + 1e-300) }
          if (!is.na(pj0) && !is.na(pjI) && pj0 > alpha && pjI <= alpha) {        # R2r => Xj <= Xi
            votes["r2r"] <- votes["r2r"] + 1L; evid["r2r"] <- evid["r2r"] - log(pjI + 1e-300) }
        }
        
        # opposite directions: keep the one with stronger aggregate evidence
        if (votes["r1"]  > 0 && votes["r1r"] > 0) {
          if (evid["r1"]  > evid["r1r"]) votes["r1r"] <- 0L
          else if (evid["r1r"] > evid["r1"])  votes["r1"]  <- 0L
          else { votes["r1"] <- 0L; votes["r1r"] <- 0L } }
        if (votes["r2"]  > 0 && votes["r2r"] > 0) {
          if (evid["r2"]  > evid["r2r"]) votes["r2r"] <- 0L
          else if (evid["r2r"] > evid["r2"])  votes["r2"]  <- 0L
          else { votes["r2"] <- 0L; votes["r2r"] <- 0L } }
        if (votes["r1"]  > 0) votes["r2"]  <- 0L          # strict supersedes weak
        if (votes["r1r"] > 0) votes["r2r"] <- 0L
        
        if (max(votes) < min_votes) next
        res <- NULL
        if (votes["r1"]  >= min_votes && votes["r1"] >= votes["r1r"]) res <- try_edge(M, i, j, 1,   0)
        if ((is.null(res) || !res$ok) && votes["r1r"] >= min_votes)   res <- try_edge(M, i, j, 0,   1)
        if ((is.null(res) || !res$ok) && votes["r2"]  >= min_votes &&
            votes["r2"] >= votes["r2r"])                              res <- try_edge(M, i, j, 0.5, 0)
        if ((is.null(res) || !res$ok) && votes["r2r"] >= min_votes)   res <- try_edge(M, i, j, 0,   0.5)
        if (!is.null(res) && res$ok) { M <- res$M; converged <- FALSE }
      }
    }
    
    # --- transitive closure: a < k < b  =>  a < b ---
    repeat {
      done <- TRUE
      for (k in 1:d_x) {
        ancs <- which(M[, k] == 1); descs <- which(M[k, ] == 1)
        for (a in ancs) for (b in descs) if (a != b && (is.na(M[a, b]) || M[a, b] != 1)) {
          r <- try_edge(M, a, b, 1, 0)
          if (r$ok) { M <- r$M; done <- FALSE; converged <- FALSE }
        }
      }
      if (done) break
    }
    
    # --- symmetry closure: a <= b AND b <= a  =>  a ~ b ---
    for (i in 1:d_x) for (j in 1:d_x) if (i != j &&
                                          !is.na(M[i, j]) && !is.na(M[j, i]) && M[i, j] == 0.5 && M[j, i] == 0.5) {
      M[i, j] <- 0; M[j, i] <- 0; converged <- FALSE
    }
    
    # --- refresh nearest ancestors over the updated non-descendant sets ---
    for (i in 1:d_x) {
      xi <- xlabs[i]; old <- pa[[xi]]
      pool <- setdiff(unique(c(z_pool[[xi]], xlabs[which(M[, i] %in% c(0.5, 1))])), xi)
      pa[[xi]] <- nearest_anc(xi, pool)
      if (!setequal(pa[[xi]], old)) converged <- FALSE
    }
  }
  
  # final R3 sweep on still-unresolved pairs (independent => ~; else leave NA)
  for (i in 2:d_x) for (j in 1:(i - 1)) if (is.na(M[i, j])) {
    pv <- ci_pval(xlabs[i], xlabs[j], build_S(i, j), C, n)
    if (!is.na(pv) && pv > alpha) { M[i, j] <- 0; M[j, i] <- 0 }
  }
  for (i in 1:d_x) for (j in 1:d_x) if (i != j &&
                                        !is.na(M[i, j]) && !is.na(M[j, i]) && M[i, j] == 0.5 && M[j, i] == 0.5) {
    M[i, j] <- 0; M[j, i] <- 0
  }
  
  # reorder to a topological order so strict edges sit above the diagonal
  # (cosmetic: no entry is zeroed)
  sup <- matrix(0, d_x, d_x)
  for (a in 1:d_x) for (b in 1:d_x)
    if (a != b && !is.na(M[a, b]) && M[a, b] %in% c(0.5, 1)) sup[a, b] <- 1
  ord <- topo_order(sup)
  if (!is.null(ord)) M <- M[ord, ord, drop = FALSE]
  diag(M) <- NA
  if (verbose) message(sprintf("  [ascend] converged after %d iteration(s)", iter))
  M
}


# ----------------------------------------------------------------------
# 6. Ground truth and evaluation
# ----------------------------------------------------------------------

# true ancestral matrix: reach[a, b] = 1 iff a is an ancestor of b
# (ancestors ordered first, same convention as ascend()'s output)
true_ancestral <- function(adj_xx) {
  d <- nrow(adj_xx); A <- t(adj_xx); diag(A) <- 0; A[is.na(A)] <- 0  # A[parent, child] = 1
  reach <- (A > 0) * 1; P <- A
  for (k in seq_len(d - 1)) { P <- ((P %*% A) > 0) * 1; reach <- ((reach + P) > 0) * 1 }
  dimnames(reach) <- dimnames(adj_xx); diag(reach) <- NA
  ord <- topo_order(reach); if (!is.null(ord)) reach <- reach[ord, ord]
  reach
}

# Evaluate the directed ancestral relation over ALL ordered pairs (a, b):
# each cell is the claim "a precedes b", so direction is scored, and the
# metric does not depend on either matrix being in any particular order.
evaluate <- function(est, truth, verbose = TRUE) {
  nms <- rownames(truth); est <- est[nms, nms, drop = FALSE]
  d <- nrow(truth)
  
  tp <- fp <- fn <- tn <- unresolved <- 0
  for (a in 1:d) for (b in 1:d) {
    if (a == b) next
    tpos <- !is.na(truth[a, b]) && truth[a, b] == 1
    e <- est[a, b]
    if (is.na(e)) { unresolved <- unresolved + 1; next }
    ppos <- (e == 1 || e == 0.5)
    if      ( tpos &&  ppos) tp <- tp + 1
    else if (!tpos &&  ppos) fp <- fp + 1
    else if ( tpos && !ppos) fn <- fn + 1
    else                     tn <- tn + 1
  }
  
  # direction accuracy and coverage over unordered pairs
  dir_ok <- dir_tot <- res_pairs <- tot_pairs <- 0
  for (a in 1:(d - 1)) for (b in (a + 1):d) {
    tot_pairs <- tot_pairs + 1
    if (!is.na(est[a, b]) || !is.na(est[b, a])) res_pairs <- res_pairs + 1
    t_ab <- !is.na(truth[a, b]) && truth[a, b] == 1
    t_ba <- !is.na(truth[b, a]) && truth[b, a] == 1
    if (!(t_ab || t_ba)) next
    c_ab <- !is.na(est[a, b]) && est[a, b] %in% c(0.5, 1)
    c_ba <- !is.na(est[b, a]) && est[b, a] %in% c(0.5, 1)
    if (!(c_ab || c_ba)) next
    dir_tot <- dir_tot + 1
    if ((t_ab && c_ab && !c_ba) || (t_ba && c_ba && !c_ab)) dir_ok <- dir_ok + 1
  }
  
  pr <- if (tp + fp) tp / (tp + fp) else NA
  re <- if (tp + fn) tp / (tp + fn) else NA
  f1 <- if (!is.na(pr) && !is.na(re) && pr + re > 0) 2 * pr * re / (pr + re) else NA
  out <- list(tp = tp, fp = fp, fn = fn, tn = tn, unresolved = unresolved,
              precision = pr, recall = re, f1 = f1,
              dir_acc = if (dir_tot) dir_ok / dir_tot else NA,
              coverage = res_pairs / tot_pairs)
  
  if (verbose) {
    line <- strrep("-", 46)
    cat(line, "\n")
    cat("  ASCEND evaluation (directed ancestral relation)\n")
    cat(line, "\n")
    cat(sprintf("  TP %-4d  FP %-4d  FN %-4d  TN %-4d  NA %-4d\n",
                tp, fp, fn, tn, unresolved))
    cat(sprintf("  precision %5.3f   recall %5.3f   F1 %5.3f\n", pr, re, f1))
    cat(sprintf("  direction accuracy %5.3f   coverage %5.3f\n", out$dir_acc, out$coverage))
    cat(line, "\n")
  }
  invisible(out)
}


# ----------------------------------------------------------------------
# 7. Pretty-print an ancestrality matrix with numeric values
#     1  a is an ancestor of b
#     0.5  a is a non-descendant of b
#     0  unrelated
#     NA  unresolved
# ----------------------------------------------------------------------

print_ancestral <- function(M, title = NULL) {
  if (!is.null(title)) cat(title, "\n")
  print(M)
  cat("  legend:  1 = ancestor   0.5 = non-descendant   0 = unrelated   NA = unresolved\n")
}


# ----------------------------------------------------------------------
# 8. Example
# ----------------------------------------------------------------------

if (sys.nframe() == 0) {
  set.seed(123)
  sim <- sim_dat(n = 1000, d_z = 20, d_x = 8,
                 r2 = 0.7, sp = 0.7, p_cross = 0.10, x_effect = 1.5, seed = 123)
  
  t0  <- Sys.time()
  est <- ascend(sim, alpha = 0.05, alpha_mb = 0.05, fdr = TRUE, min_votes = 1)
  rt  <- as.numeric(Sys.time() - t0, units = "secs")
  
  truth <- true_ancestral(sim$adj_xx)
  cat(sprintf("\nruntime: %.2f s\n\n", rt))
  print_ancestral(truth, "True ancestral matrix")
  cat("\n")
  print_ancestral(est, "ASCEND estimate")
  cat("\n")
  evaluate(est, truth)
}