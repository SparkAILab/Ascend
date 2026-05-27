## ASCEND: Ancestral Scalable Causal discovEry via iNherited Descent.
#######################################################################

## Constraint-based ancestral causal discovery for two-tier systems
## with a known causal ordering: a background tier Z precedes a
## foreground tier X.
##
## Input:
##   sim_obj : list with element $dat (a data.table whose columns
##             follow the convention z1, ..., z_dz, x1, ..., x_dx).
##
## Output:
##   d_x x d_x matrix M over the foreground variables with entries in
##   {0, 0.5, 1, NA}:
##     M[i, j] = 1    Xi is a strict ancestor of Xj
##     M[i, j] = 0.5  Xi is a non-descendant of Xj (weak precedence)
##     M[i, j] = 0    Xi and Xj are independent given S
##     M[i, j] = NA   pair could not be resolved given the available
##                    conditioning sets and data
##
## Orientation rules applied to a pair (Xi, Xj) with conditioning
## set S and helper W in S, Sw = S \ {W}:
##
##   R3 (independence):
##     Xi _|_ Xj | S  =>  Xi ~ Xj
##
##   R1 (deactivation):
##     W _|_/ Xj | Sw  and  W _|_ Xj | Sw u {Xi}  =>  Xi -< Xj
##     Conditioning on Xi screens W from Xj, so Xi mediates the path
##     from W to Xj.
##
##   R2 (activation):
##     W _|_ Xi | Sw  and  W _|_/ Xi | Sw u {Xj}  =>  Xi <= Xj
##     Conditioning on Xj opens a path between W and Xi, so Xj is a
##     collider or descendant of a collider on that path.
##
## Markov-blanket alpha schedule:
##   The MB oracle's significance level alpha_mb is annealed from
##   alpha_mb_start (liberal) to alpha_mb_floor (conservative) by a
##   geometric factor alpha_decay per iteration. Liberal alpha early
##   minimises omissions from the conditioning set; commission
##   errors are tolerated because a slightly oversized valid
##   conditioning set does not invalidate the rules under
##   faithfulness.
##
## Vote aggregation:
##   The oracle pseudocode commits orientation on the first valid W.
##   In finite samples a single W may fire spuriously. We collect
##   votes from every W in S, take the plurality winner subject to
##   min_votes, and break direction conflicts (R1 vs R1_rev, R2 vs
##   R2_rev) using the total -log(p) evidence across the significant
##   tests. Under perfect tests this recovers the oracle behaviour.

suppressPackageStartupMessages({
  library(bnlearn)
  library(data.table)
  library(matrixStats)
  library(igraph)
})

## --- Graph utilities ------------------------------------------------------

clean_matrix_adj <- function(mat) {
  m <- as.matrix(mat); m[is.na(m)] <- 0; m[m != 0] <- 1; m
}

topo_sort_adj <- function(adj) {
  adj <- clean_matrix_adj(adj); n <- nrow(adj); if (n == 0) return(integer(0))
  A <- (adj != 0) * 1; indeg <- colSums(A); indeg[is.na(indeg)] <- 0
  ord <- integer(0); zeros <- sort(which(indeg == 0))
  while (length(zeros) > 0) {
    v <- zeros[1]; zeros <- zeros[-1]; ord <- c(ord, v)
    nbrs <- which(A[v, ] != 0)
    for (u in nbrs) {
      indeg[u] <- indeg[u] - 1
      if (indeg[u] == 0) zeros <- sort(unique(c(zeros, u)))
    }
    A[v, nbrs] <- 0
  }
  if (length(ord) != n) return(NULL); ord
}

is_dag_adj <- function(adj) !is.null(topo_sort_adj(adj))

sort_ancestral_matrix <- function(m) {
  if (is.null(m) || nrow(m) == 0) return(m)
  sup <- m; sup[is.na(sup)] <- 0
  sup[!(sup %in% c(0.5, 1))] <- 0; sup[sup != 0] <- 1
  ord <- topo_sort_adj(sup)
  if (is.null(ord)) { warning("Not a DAG"); return(m) }
  s <- m[ord, ord, drop = FALSE]; diag(s) <- NA
  s[lower.tri(s)] <- 0; s
}

build_constraint_adj <- function(m) {
  n <- nrow(m); A <- matrix(0, n, n)
  rownames(A) <- rownames(m); colnames(A) <- colnames(m)
  for (j in seq_len(n)) for (i in seq_len(n)) {
    if (i == j || is.na(m[j, i])) next
    if (m[j, i] == 1 || m[j, i] == 0.5) A[i, j] <- 1
  }
  A
}

permute_to_lower <- function(m) {
  ord <- topo_sort_adj(build_constraint_adj(m))
  if (is.null(ord)) stop("Ancestry matrix contains cycles")
  m[as.numeric(ord), as.numeric(ord)]
}

# Attempt to update entry (i, j) of the ancestral matrix. The edge is
# accepted only if the resulting graph is acyclic.
try_edge_update <- function(adj, i, j, val_ij, val_ji) {
  t <- adj; t[i, j] <- val_ij; t[j, i] <- val_ji
  b <- ((t == 1 | t == 0.5) * 1); b[is.na(b)] <- 0
  if (is_dag_adj(b)) return(list(adj = t, success = TRUE))
  list(adj = adj, success = FALSE)
}

## --- Conditional independence test ---------------------------------------
##
## Partial F-test for H0: target _|_ query | cond_set, via a
## likelihood-ratio comparison between nested linear models. Small
## p-value rejects independence.

ci_test_pval <- function(target, query, cond_set, data, min_n = 10) {
  vars <- unique(c(target, query, cond_set))
  miss <- setdiff(vars, colnames(data))
  if (length(miss) > 0) return(NA)
  if (nrow(data) < min_n) return(NA)
  
  sub <- data[, vars, drop = FALSE]
  good <- sapply(sub, function(col) {
    v <- na.omit(col); length(unique(v)) > 1 && sd(v) > 0
  })
  sub <- sub[, as.logical(good), drop = FALSE]
  if (!(target %in% colnames(sub)) || !(query %in% colnames(sub))) return(NA)
  
  cs <- intersect(cond_set, colnames(sub))
  
  rhs0 <- if (length(cs) > 0) paste(cs, collapse = "+") else "1"
  rhs1 <- if (length(cs) > 0) paste(c(cs, query), collapse = "+") else query
  
  f0 <- tryCatch(lm(as.formula(paste(target, "~", rhs0)), data = sub),
                 error = function(e) NULL)
  f1 <- tryCatch(lm(as.formula(paste(target, "~", rhs1)), data = sub),
                 error = function(e) NULL)
  if (is.null(f0) || is.null(f1)) return(NA)
  
  res <- tryCatch(anova(f0, f1), error = function(e) NULL)
  if (is.null(res)) return(NA)
  res$`Pr(>F)`[2]
}

## --- Main routine --------------------------------------------------------

ascend_fn <- function(sim_obj,
                      maxiter        = 10,
                      alpha          = 0.05,
                      alpha_mb_start = 0.20,
                      alpha_mb_floor = 0.05,
                      alpha_decay    = 0.70,
                      fdr_correction = TRUE,
                      min_votes      = 1) {
  
  dat   <- sim_obj$dat
  d_x   <- length(grep("^x", colnames(dat)))
  xlabs <- paste0('x', seq_len(d_x))
  
  # Median/MAD scaling, robust to outliers; falls back to mean/sd if the
  # MAD is zero, and adds tiny jitter if both are degenerate.
  scale_robust <- function(x) {
    if (all(is.na(x))) return(x)
    m <- mad(x, na.rm = TRUE)
    if (is.na(m) || m == 0) {
      s <- sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) {
        x <- x + rnorm(length(x), 0, 1e-6); s <- sd(x, na.rm = TRUE)
      }
      return((x - mean(x, na.rm = TRUE)) / s)
    }
    (x - median(x, na.rm = TRUE)) / m
  }
  
  dataAll <- as.data.frame(dat)
  for (col in colnames(dataAll)) dataAll[[col]] <- scale_robust(dataAll[[col]])
  dataAll[is.infinite(as.matrix(dataAll))] <- 0
  dataAll[is.na(dataAll)] <- 0
  
  z_cols <- grep("^z", colnames(dataAll), value = TRUE)
  
  # Marginal Z screening. The two-tier assumption guarantees any subset
  # of Z is a valid set of non-descendants, so we retain Z with any
  # marginal signal (BH-corrected threshold of 0.30) plus a minimum
  # count so that the initial conditioning set is never empty.
  min_z_keep <- max(3, ceiling(length(z_cols) * 0.10))
  
  prescreen_z <- function(node_i, liberal_thresh = 0.30) {
    if (sd(dataAll[[node_i]], na.rm = TRUE) == 0) return(character(0))
    pv <- sapply(z_cols, function(zc) {
      ct <- tryCatch(cor.test(dataAll[[node_i]], dataAll[[zc]]),
                     error = function(e) NULL)
      if (is.null(ct)) 1 else ct$p.value
    })
    adj_pv <- p.adjust(pv, method = "BH")
    keep   <- which(adj_pv < liberal_thresh)
    if (length(keep) < min_z_keep)
      keep <- order(pv)[1:min(min_z_keep, length(pv))]
    z_cols[keep]
  }
  
  # IAMB Markov-blanket learning at significance alpha_mb, with guards
  # against degenerate columns.
  learn_mb_safe <- function(node_i, candidates, alpha_mb) {
    if (length(candidates) == 0 || sd(dataAll[[node_i]], na.rm = TRUE) == 0)
      return(character(0))
    sub <- dataAll[, c(candidates, node_i), drop = FALSE]
    good <- sapply(sub, function(col) {
      v <- na.omit(col); length(unique(v)) > 1 && sd(v) > 0
    })
    sub <- sub[, as.logical(good), drop = FALSE]
    if (ncol(sub) < 2 || !(node_i %in% colnames(sub))) return(character(0))
    tryCatch(
      learn.mb(sub, node = node_i, method = "iamb", test = "zf",
               alpha = alpha_mb),
      error = function(e) character(0)
    )
  }
  
  ## Initial MB at iteration 0. T_X^(0) = Z; alpha_mb is at its most
  ## liberal value.
  alpha_mb_current <- alpha_mb_start
  
  mb_list <- vector("list", d_x); names(mb_list) <- xlabs
  for (i in seq_len(d_x)) {
    node_i  <- xlabs[i]
    z_cands <- prescreen_z(node_i)
    mb_list[[node_i]] <- learn_mb_safe(node_i, z_cands, alpha_mb_current)
  }
  
  M <- matrix(NA, d_x, d_x); diag(M) <- NA
  rownames(M) <- colnames(M) <- xlabs
  
  converged <- FALSE; iter <- 0
  
  while (!converged && iter <= maxiter) {
    iter <- iter + 1; converged <- TRUE
    
    alpha_mb_current <- max(alpha_mb_floor,
                            alpha_mb_start * alpha_decay ^ (iter - 1))
    
    cat(sprintf("Iteration %d | alpha_mb=%.3f | NA pairs remaining: %d\n",
                iter, alpha_mb_current,
                sum(is.na(M[upper.tri(M)]))))
    
    ## --- R3: pairwise CI tests on unresolved pairs ---------------------
    ## S_ij = mb[xi] u mb[xj], excluding xi and xj themselves.
    pval_info <- list(); pidx <- 1
    
    for (i in 2:d_x) {
      for (j in 1:(i - 1)) {
        if (!is.na(M[i, j])) next
        
        xi <- xlabs[i]; xj <- xlabs[j]
        S <- setdiff(union(mb_list[[xi]], mb_list[[xj]]), c(xi, xj))
        S <- intersect(S, colnames(dataAll))
        
        pv <- ci_test_pval(xj, xi, S, dataAll)
        if (!is.na(pv)) {
          pval_info[[pidx]] <- list(i = i, j = j, xi = xi, xj = xj,
                                    pval = pv, S = S)
          pidx <- pidx + 1
        } else {
          converged <- FALSE
        }
      }
    }
    
    if (length(pval_info) > 0) {
      pvals     <- sapply(pval_info, function(x) x$pval)
      adj_pvals <- if (fdr_correction) p.adjust(pvals, method = "BH") else pvals
      
      for (idx in seq_along(pval_info)) {
        info <- pval_info[[idx]]
        i <- info$i; j <- info$j; xi <- info$xi; xj <- info$xj; S <- info$S
        adj_pval <- adj_pvals[idx]
        
        # R3: independence -> Xi ~ Xj.
        if (adj_pval > alpha) {
          M[i, j] <- 0; M[j, i] <- 0; converged <- FALSE; next
        }
        
        # Dependent: collect R1 and R2 votes across all W in S.
        if (length(S) == 0) next
        
        votes    <- c(r1 = 0L,  r1_rev = 0L,  r2 = 0L,  r2_rev = 0L)
        evidence <- c(r1 = 0.0, r1_rev = 0.0, r2 = 0.0, r2_rev = 0.0)
        
        for (W in S) {
          Sw <- setdiff(S, W)
          p_wj_no_i <- ci_test_pval(W, xj, Sw,         dataAll)
          p_wj_wi_i <- ci_test_pval(W, xj, c(Sw, xi),  dataAll)
          p_wi_no_j <- ci_test_pval(W, xi, Sw,         dataAll)
          p_wi_wi_j <- ci_test_pval(W, xi, c(Sw, xj),  dataAll)
          
          r1f  <- !is.na(p_wj_no_i) && !is.na(p_wj_wi_i) &&
            p_wj_no_i <= alpha && p_wj_wi_i > alpha
          r1rf <- !is.na(p_wi_no_j) && !is.na(p_wi_wi_j) &&
            p_wi_no_j <= alpha && p_wi_wi_j > alpha
          r2f  <- !is.na(p_wi_no_j) && !is.na(p_wi_wi_j) &&
            p_wi_no_j > alpha  && p_wi_wi_j <= alpha
          r2rf <- !is.na(p_wj_no_i) && !is.na(p_wj_wi_i) &&
            p_wj_no_i > alpha  && p_wj_wi_i <= alpha
          
          votes["r1"]     <- votes["r1"]     + r1f
          votes["r1_rev"] <- votes["r1_rev"] + r1rf
          votes["r2"]     <- votes["r2"]     + r2f
          votes["r2_rev"] <- votes["r2_rev"] + r2rf
          
          if (r1f  && !is.na(p_wj_no_i)) evidence["r1"]     <- evidence["r1"]     - log(p_wj_no_i + 1e-300)
          if (r1rf && !is.na(p_wi_no_j)) evidence["r1_rev"] <- evidence["r1_rev"] - log(p_wi_no_j + 1e-300)
          if (r2f  && !is.na(p_wi_wi_j)) evidence["r2"]     <- evidence["r2"]     - log(p_wi_wi_j + 1e-300)
          if (r2rf && !is.na(p_wj_wi_i)) evidence["r2_rev"] <- evidence["r2_rev"] - log(p_wj_wi_i + 1e-300)
        }
        
        # Resolve direction conflicts by total -log(p) evidence; exact
        # ties leave the pair unresolved.
        if (votes["r1"] > 0 && votes["r1_rev"] > 0) {
          if      (evidence["r1"]     > evidence["r1_rev"]) votes["r1_rev"] <- 0L
          else if (evidence["r1_rev"] > evidence["r1"])     votes["r1"]     <- 0L
          else { votes["r1"] <- 0L; votes["r1_rev"] <- 0L }
        }
        if (votes["r2"] > 0 && votes["r2_rev"] > 0) {
          if      (evidence["r2"]     > evidence["r2_rev"]) votes["r2_rev"] <- 0L
          else if (evidence["r2_rev"] > evidence["r2"])     votes["r2"]     <- 0L
          else { votes["r2"] <- 0L; votes["r2_rev"] <- 0L }
        }
        # Strict orientation supersedes weak in the same direction.
        if (votes["r1"]     > 0) votes["r2"]     <- 0L
        if (votes["r1_rev"] > 0) votes["r2_rev"] <- 0L
        
        best <- max(votes)
        if (best < min_votes) next
        
        # Try strict first, then weak; accept the first update that
        # leaves M acyclic.
        result <- NULL
        if (votes["r1"] >= min_votes && votes["r1"] >= votes["r1_rev"])
          result <- try_edge_update(M, i, j, 1, 0)
        if ((is.null(result) || !result$success) &&
            votes["r1_rev"] >= min_votes)
          result <- try_edge_update(M, i, j, 0, 1)
        if ((is.null(result) || !result$success) &&
            votes["r2"] >= min_votes && votes["r2"] >= votes["r2_rev"])
          result <- try_edge_update(M, i, j, 0.5, 0)
        if ((is.null(result) || !result$success) &&
            votes["r2_rev"] >= min_votes)
          result <- try_edge_update(M, i, j, 0, 0.5)
        
        if (!is.null(result) && result$success) {
          M <- result$adj; converged <- FALSE
        }
      }
    }
    
    ## --- Transitive closure --------------------------------------------
    ## Xi -< Xk and Xk -< Xj imply Xi -< Xj.
    done <- FALSE
    while (!done) {
      done <- TRUE
      for (k in seq_len(d_x)) {
        ancs  <- which(M[, k] == 1)
        descs <- which(M[k, ] == 1)
        if (!length(ancs) || !length(descs)) next
        for (anc in ancs) for (desc in descs) {
          if (anc == desc) next
          if (is.na(M[anc, desc]) || M[anc, desc] != 1) {
            res <- try_edge_update(M, anc, desc, 1, 0)
            if (res$success) { M <- res$adj; done <- FALSE; converged <- FALSE }
          }
        }
      }
    }
    
    ## --- Symmetry closure ---------------------------------------------
    ## Xi <= Xj together with Xj <= Xi contradicts acyclicity, so the
    ## pair is reset to independent.
    for (i in seq_len(d_x)) for (j in seq_len(d_x)) {
      if (i == j) next
      if (!is.na(M[i, j]) && !is.na(M[j, i]) &&
          M[i, j] == 0.5 && M[j, i] == 0.5) {
        M[i, j] <- 0; M[j, i] <- 0; converged <- FALSE
      }
    }
    
    bin <- clean_matrix_adj(M); bin[is.na(bin)] <- 0
    if (!is_dag_adj(bin)) stop(sprintf("Cycle detected at iteration %d", iter))
    
    ## --- Update non-descendant sets and recompute MBs ------------------
    ## T_Xi shrinks each iteration toward the true parent set: the Z
    ## variables currently in Xi's MB plus any X confirmed as
    ## ancestors of Xi. Fallback to the marginal Z prescreen if no Z
    ## parents have been found yet.
    for (i in seq_len(d_x)) {
      node_i  <- xlabs[i]; old_mb <- mb_list[[node_i]]
      
      z_parents_Xi <- intersect(mb_list[[node_i]], z_cols)
      if (length(z_parents_Xi) == 0) z_parents_Xi <- prescreen_z(node_i)
      
      confirmed_x_anc <- xlabs[which(M[, i] == 1 | M[, i] == 0.5)]
      T_Xi <- setdiff(unique(c(z_parents_Xi, confirmed_x_anc)), node_i)
      T_Xi <- intersect(T_Xi, colnames(dataAll))
      
      new_mb <- learn_mb_safe(node_i, T_Xi, alpha_mb_current)
      
      if (!setequal(new_mb, old_mb)) converged <- FALSE
      mb_list[[node_i]] <- new_mb
    }
    
  }
  
  ## --- Post-loop R3 pass --------------------------------------------
  ## After the iterative loop some pairs may remain NA. R3 is applied
  ## once more using the final MBs: a pair is set to 0 only if the CI
  ## test now finds independence. Dependent but unoriented pairs are
  ## left as NA: the algorithm does not force an orientation that the
  ## rules did not establish.
  for (i in 2:d_x) {
    for (j in 1:(i - 1)) {
      if (!is.na(M[i, j])) next
      xi <- xlabs[i]; xj <- xlabs[j]
      S  <- setdiff(union(mb_list[[xi]], mb_list[[xj]]), c(xi, xj))
      S  <- intersect(S, colnames(dataAll))
      pv <- ci_test_pval(xj, xi, S, dataAll)
      if (!is.na(pv) && pv > alpha) {
        M[i, j] <- 0; M[j, i] <- 0
      }
    }
  }
  
  for (i in seq_len(d_x)) for (j in seq_len(d_x)) {
    if (i == j) next
    if (!is.na(M[i, j]) && !is.na(M[j, i]) &&
        M[i, j] == 0.5 && M[j, i] == 0.5) {
      M[i, j] <- 0; M[j, i] <- 0
    }
  }
  
  # Topologically sort the final matrix for presentation.
  if (!all(is.na(M))) {
    tryCatch({
      M <- permute_to_lower(M)
      M <- sort_ancestral_matrix(M)
    }, error = function(e) warning("Final normalisation failed: ", e$message))
  }
  
  M
}

## --- Example run ---------------------------------------------
##
## # source("simulation.R")
## # sim_obj <- sim_dat(n = 500, d_z = 15, d_x = 6, r2 = 0.5,
## #                    sp = 0.3, p_cross = 0.15, x_effect = 0.9, seed = 123)
## # M <- ascend_fn(sim_obj)