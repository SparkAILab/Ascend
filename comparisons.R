# ==========================================================
# ASCEND Framework 
# ==========================================================

# Load required packages 
library(bnlearn)
library(dplyr)
library(qvalue)
library(foreach)
library(data.table)
library(pcalg)
library(RBGL)
library(matrixStats)
library(glmnet)
library(lightgbm)
library(tidyverse)
library(doMC)
registerDoMC(16)
set.seed(123, kind = "L'Ecuyer-CMRG")

# ======================================================================
# Helper functions replicating igraph behavior 
# ======================================================================

# Clean matrix: treat NA as 0, and any nonzero value (like 0.5 or 2) as 1
clean_matrix_adj <- function(mat) {
  mat2 <- as.matrix(mat)
  mat2[is.na(mat2)] <- 0
  mat2[mat2 != 0] <- 1
  return(mat2)
}

# Recursive DFS reachability, replicates igraph::subcomponent(mode = "out")
reachable_nodes_adj <- function(adj, start) {
  adj <- clean_matrix_adj(adj)
  n <- nrow(adj)
  visited <- rep(FALSE, n)
  dfs <- function(v) {
    visited[v] <<- TRUE
    nbrs <- which(adj[v, ] != 0)
    for (u in nbrs) {
      if (!visited[u]) dfs(u)
    }
  }
  dfs(start)
  which(visited)
}

# Kahn topological sort — deterministic, smallest index first (matches igraph)
topo_sort_adj <- function(adj) {
  adj <- clean_matrix_adj(adj)
  n <- nrow(adj)
  if (n == 0) return(integer(0))
  A <- (adj != 0) * 1
  indeg <- colSums(A)
  indeg[is.na(indeg)] <- 0
  
  order <- integer(0)
  zeros <- sort(which(indeg == 0))  # deterministic
  
  while (length(zeros) > 0) {
    v <- zeros[1]
    zeros <- zeros[-1]
    order <- c(order, v)
    nbrs <- which(A[v, ] != 0)
    if (length(nbrs) > 0) {
      for (u in nbrs) {
        indeg[u] <- indeg[u] - 1
        if (indeg[u] == 0) zeros <- sort(unique(c(zeros, u)))
      }
    }
    A[v, nbrs] <- 0
  }
  
  if (length(order) != n) {
    return(NULL)  # cycle detected
  }
  return(order)
}

# DAG check: true if topo-sort succeeds
is_dag_adj <- function(adj) {
  !is.null(topo_sort_adj(adj))
}
# ======================================================================

# find_simple_cycles_adj: find simple cycles up to max_len (returns list of node index vectors)
find_simple_cycles_adj <- function(adj, max_len = 4) {
  adj <- as.matrix(adj)
  n <- nrow(adj)
  cycles <- list()
  for (start in seq_len(n)) {
    dfs <- function(path) {
      last <- tail(path, 1)
      if (length(path) > 1 && last == start) {
        cycles <<- append(cycles, list(path))
        return()
      }
      if (length(path) >= max_len) return()
      nbrs <- which(adj[last, ] != 0)
      for (u in nbrs) {
        if (u == start) {
          cycles <<- append(cycles, list(c(path, u)))
        } else if (!(u %in% path)) {
          dfs(c(path, u))
        }
      }
    }
    dfs(c(start))
  }
  return(cycles)
}

# Build constraint adjacency matrix from ancestry-like matrix 'm'
# In your earlier logic: if m[j,i] == 1 or 0.5 then add edge i -> j
build_constraint_adj <- function(m) {
  n <- nrow(m)
  A <- matrix(0, n, n)
  rownames(A) <- rownames(m)
  colnames(A) <- colnames(m)
  for (j in seq_len(n)) {
    for (i in seq_len(n)) {
      if (i == j) next
      if (is.na(m[j, i])) next
      if (m[j, i] == 1 || m[j, i] == 0.5) {
        A[i, j] <- 1
      }
    }
  }
  return(A)
}

# Permute matrix to lower triangular form using constraints (throws if cycle)
permute_to_lower <- function(m) {
  get_valid_permutation <- function(m) {
    A <- build_constraint_adj(m)
    ord <- topo_sort_adj(A)
    if (is.null(ord)) stop("Ancestry matrix contains cycles")
    return(as.numeric(ord))
  }
  perm <- get_valid_permutation(m)
  return(m[perm, perm])
}

# try_edge_update: test DAG property after proposed change
try_edge_update <- function(adj, i, j, val_ij, val_ji, iter) {
  test_adj <- adj
  test_adj[i, j] <- val_ij
  test_adj[j, i] <- val_ji
  bin <- (test_adj == 1 | test_adj == 0.5) * 1
  if (is_dag_adj(bin)) {
    return(list(adj = test_adj, success = TRUE))
  } else {
    return(list(adj = adj, success = FALSE))
  }
}

# thin wrapper to match earlier function name
find_simple_cycles <- function(graph_adj, max_len = 4) {
  find_simple_cycles_adj(graph_adj, max_len = max_len)
}

# ------------------------------------------------------------------
# get_mb helper: compute MB(Xi, Ai) safely (only over Ai ∪ {xi})
# ------------------------------------------------------------------
get_mb <- function(xi, Ai, dataAll, alpha = 0.05, method = "gs", test = "mi-g") {
  # xi: single variable name (character)
  # Ai: character vector of candidate variables (may include z* and some x*)
  if (!is.character(xi) || length(xi) != 1) stop("xi must be a single variable name")
  if (!(xi %in% colnames(dataAll))) stop(paste("xi not found in data:", xi))
  
  # ensure Ai is subset of available columns and remove xi if present
  Ai <- unique(intersect(Ai, colnames(dataAll)))
  Ai <- setdiff(Ai, xi)
  
  # select columns to pass to learn.mb
  cols <- intersect(c(Ai, xi), colnames(dataAll))
  if (length(cols) == 0) return(character(0))
  sub_df <- dataAll[, cols, drop = FALSE]
  
  # Drop constant columns
  keep <- sapply(sub_df, function(col) length(unique(col)) > 1)
  sub_df <- sub_df[, keep, drop = FALSE]
  
  # If after dropping constants we only have xi, return empty MB
  if (ncol(sub_df) <= 1) return(character(0))
  
  mb_res <- try(learn.mb(sub_df, node = xi, method = method, test = test, alpha = alpha), silent = TRUE)
  if (inherits(mb_res, "try-error")) {
    warning("learn.mb failed for ", xi, " with Ai = {", paste(Ai, collapse = ","), "}; returning character(0)")
    return(character(0))
  }
  mb <- as.character(mb_res)
  mb <- intersect(mb, colnames(sub_df))
  return(mb)
}
# ------------------------------------------------------------------

# normalize_to_sorted_ancestral: replaces igraph-based normalization
normalize_to_sorted_ancestral <- function(adj_matrix, is_ancestral = NULL) {
  clean_matrix <- function(mat) {
    mat_clean <- mat
    mat_clean[is.na(mat_clean)] <- 0
    return(mat_clean)
  }
  
  is_ancestral_matrix <- function(mat) {
    mat_clean <- clean_matrix(mat)
    d <- nrow(mat_clean)
    for (i in 1:d) {
      true_descendants <- setdiff(reachable_nodes_adj(mat_clean, i), i)
      matrix_descendants <- which(mat_clean[i, ] == 1)
      if (!setequal(true_descendants, matrix_descendants)) {
        return(FALSE)
      }
    }
    return(TRUE)
  }
  
  dag_to_ancestral <- function(mat) {
    mat_clean <- clean_matrix(mat)
    d <- nrow(mat_clean)
    ancestral_mat <- matrix(0, nrow = d, ncol = d)
    if (!is.null(rownames(mat))) rownames(ancestral_mat) <- rownames(mat)
    if (!is.null(colnames(mat))) colnames(ancestral_mat) <- colnames(mat)
    for (i in 1:d) {
      reachable <- reachable_nodes_adj(mat_clean, i)
      reachable <- setdiff(reachable, i)
      if (length(reachable) > 0) ancestral_mat[i, reachable] <- 1
    }
    diag(ancestral_mat) <- NA
    return(ancestral_mat)
  }
  
  topological_sort_matrix <- function(mat) {
    mat_clean <- clean_matrix(mat)
    if (!is_dag_adj(mat_clean)) {
      warning("Matrix is not a DAG, cannot perform topological sort. Returning original order.")
      return(list(order = 1:nrow(mat), sorted_mat = mat))
    }
    topo_order <- topo_sort_adj(mat_clean)
    sorted_mat <- mat[topo_order, topo_order, drop = FALSE]
    return(list(order = topo_order, sorted_mat = sorted_mat))
  }
  
  adj_clean <- clean_matrix(adj_matrix)
  if (is.null(is_ancestral)) {
    is_ancestral_input <- is_ancestral_matrix(adj_clean)
  } else {
    is_ancestral_input <- is_ancestral
  }
  if (!is_ancestral_input) {
    ancestral_mat <- dag_to_ancestral(adj_clean)
  } else {
    ancestral_mat <- adj_clean
  }
  sort_result <- topological_sort_matrix(ancestral_mat)
  diag(sort_result$sorted_mat) <- NA
  return(sort_result$sorted_mat)
}

# ==========================================================
# 3. CI & edge tests
# ==========================================================

test_ci <- function(xi, xj, df, alpha = 0.05) {
  vars <- intersect(c(xi, xj), names(df))
  if (length(vars) < 2) return(FALSE)
  if (ncol(df) <= 2) return(FALSE)
  cond <- setdiff(names(df), c(xi, xj))
  if (length(cond) == 0) {
    fit_null <- lm(as.formula(paste(xj, "~ 1")), data = df)
    fit_alt <- lm(as.formula(paste(xj, "~", xi)), data = df)
  } else {
    cond_str <- paste(cond, collapse = " + ")
    fit_null <- lm(as.formula(paste(xj, "~", cond_str)), data = df)
    fit_alt  <- lm(as.formula(paste(xj, "~", paste(c(cond, xi), collapse = " + "))), data = df)
  }
  test_res <- anova(fit_null, fit_alt)
  pval <- test_res$`Pr(>F)`[2]
  return(pval > alpha)
}

test_r1r2 <- function(W, xj, S_w, xi, data, alpha = 0.05) {
  vars <- unique(c(xj, W, S_w, xi))
  miss <- setdiff(vars, colnames(data))
  if (length(miss)) {
    return(list(beta_0 = NA, beta_1 = NA, pval_0 = NA, pval_1 = NA,
                p_value_diff = NA, r1 = FALSE, r2 = FALSE))
  }
  form0 <- as.formula(paste(xj, "~", paste(c(W, S_w), collapse = "+")))
  form1 <- as.formula(paste(xj, "~", paste(c(W, S_w, xi), collapse = "+")))
  f0 <- lm(form0, data)
  f1 <- lm(form1, data)
  coef0 <- summary(f0)$coefficients
  coef1 <- summary(f1)$coefficients
  if (!(W %in% rownames(coef0)) || !(W %in% rownames(coef1))) {
    return(list(beta_0 = NA, beta_1 = NA, pval_0 = NA, pval_1 = NA,
                p_value_diff = NA, r1 = FALSE, r2 = FALSE))
  }
  beta_0 <- coef0[W, 1]; se_0 <- coef0[W, 2]; pval_0 <- coef0[W, 4]
  beta_1 <- coef1[W, 1]; se_1 <- coef1[W, 2]; pval_1 <- coef1[W, 4]
  se_diff <- sqrt(se_0^2 + se_1^2)
  nu_combined <- se_diff^4 / (se_0^4/f0$df.residual + se_1^4/f1$df.residual)
  t_stat <- (beta_0 - beta_1) / se_diff
  p_value_diff <- 2 * pt(-abs(t_stat), df = nu_combined)
  r1 <- (pval_0 <= alpha && pval_1 > alpha)
  r2 <- (pval_0 > alpha && pval_1 <= alpha)
  return(list(beta_0 = beta_0, beta_1 = beta_1, pval_0 = pval_0,
              pval_1 = pval_1, p_value_diff = p_value_diff, r1 = r1, r2 = r2))
}

# ==========================================================
# 4. Main ASCEND function 
# ==========================================================

ascend_fn <- function(sim_obj, maxiter = 9, alpha = 0.05, verbose = TRUE) {
  dat <- sim_obj$dat
  z <- as.matrix(select(dat, starts_with('z')))
  x <- as.matrix(select(dat, starts_with('x')))
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  dataAll <- as.data.frame(scale(dat))
  z_cols <- grep("^z", colnames(dataAll), value = TRUE)
  x_cols <- grep("^x", colnames(dataAll), value = TRUE)
  
  # Initialize A_i = Z for each X_i, and learn Anc^*(X_i, A_i) = MB(X_i; A_i)
  A_list <- vector("list", length = d_x); names(A_list) <- xlabs
  mb_list <- vector("list", length = d_x); names(mb_list) <- xlabs
  for (i in seq_len(d_x)) {
    xi <- xlabs[i]
    A_list[[xi]] <- z_cols
    mb_list[[xi]] <- get_mb(xi, A_list[[xi]], dataAll, alpha = alpha)
  }
  
  if (verbose) {
    cat("\n==============================\n")
    cat("Iteration 0: Markov Blankets and A_i\n")
    cat("==============================\n")
    for (v in xlabs) {
      cat("A_", v, " = {", paste(A_list[[v]], collapse = ", "), "}\n", sep = "")
      cat("MB(", v, ") = { ", paste(mb_list[[v]], collapse = ", "), " }\n", sep = "")
    }
    cat("==============================\n\n")
  }
  
  adj_new <- matrix(NA, d_x, d_x); diag(adj_new) <- NA
  converged <- FALSE
  iter <- 0
  while (!converged && iter <= maxiter) {
    iter <- iter + 1
    converged <- TRUE
    
    # Pairwise tests
    for (i in 2:d_x) {
      for (j in 1:(i - 1)) {
        if (is.na(adj_new[i, j])) {
          xi <- xlabs[i]; xj <- xlabs[j]
          
          # S is union of current Anc^*(xi,Ai) and Anc^*(xj,Aj)
          Si <- setdiff(mb_list[[xi]], c(xi, xj))
          Sj <- setdiff(mb_list[[xj]], c(xi, xj))
          S <- setdiff(unique(c(Si, Sj)), c(xi, xj))
          
          S_data_cols <- intersect(c(xi, xj, S), colnames(dataAll))
          S_data <- dataAll[, S_data_cols, drop = FALSE]
          keep <- sapply(S_data, function(col) length(unique(col)) > 1)
          S_data <- S_data[, keep, drop = FALSE]
          
          if (test_ci(xi, xj, S_data)) {
            adj_new[i, j] <- 0; adj_new[j, i] <- 0
            converged <- FALSE
          } else {
            # try to orient edges using r1/r2 tests
            for (W in S) {
              Sw <- setdiff(S, W)
              act <- test_r1r2(W, xj, Sw, xi, dataAll, alpha)
              deact <- test_r1r2(W, xi, Sw, xj, dataAll, alpha)
              
              result <- NULL
              if (!is.null(deact) && is.list(deact) && deact$r1) {
                result <- try_edge_update(adj_new, i, j, 1, 0, iter)
              } else if (!is.null(act) && is.list(act) && act$r1) {
                result <- try_edge_update(adj_new, i, j, 0, 1, iter)
              } else if (!is.null(deact) && is.list(deact) && deact$r2) {
                result <- try_edge_update(adj_new, i, j, 0, 0.5, iter)
              } else if (!is.null(act) && is.list(act) && act$r2) {
                result <- try_edge_update(adj_new, i, j, 0.5, 0, iter)
              }
              
              if (!is.null(result) && result$success) {
                adj_new <- result$adj
                converged <- FALSE
                break
              }
            } # end W loop
          } # end else not CI
        } # end if is.na
      } # end j
    } # end i
    
    # Check DAG property of the current binary adjacency
    bin_current <- (clean_matrix_adj(adj_new) == 1 | clean_matrix_adj(adj_new) == 0.5) * 1
    if (!is_dag_adj(bin_current)) stop("Cycle detected")
    
    # Compute ancestral closure (propagate 1s)
    closure <- FALSE
    while (!closure) {
      closure <- TRUE
      for (k in seq_len(ncol(adj_new))) {
        parents_k <- which(adj_new[, k] == 1)
        if (length(parents_k) > 0) {
          submat <- adj_new[, parents_k, drop = FALSE]
          if (any(submat == 1, na.rm = TRUE)) {
            e_idx <- unique(which(submat == 1, arr.ind = TRUE)[, 1])
            for (row in e_idx) {
              if (is.na(adj_new[row, k]) || adj_new[row, k] != 1) {
                adj_new[row, k] <- 1
                closure <- FALSE
              }
            }
          }
        }
      }
    } # end closure
    
    # Recompute A_i and relearn Anc^*(X_i, A_i) for each node
    for (i in seq_len(d_x)) {
      xi <- xlabs[i]
      old_mb <- mb_list[[xi]]
      preds_idx <- which((adj_new[, i] == 1) | (adj_new[, i] == 0.5))
      preds <- xlabs[preds_idx]
      # A_i <- Z ∪ Anc^*(X_i, A_i) ∪ preds
      A_i_new <- unique(c(z_cols, old_mb, preds))
      A_list[[xi]] <- A_i_new
      # Relearn MB restricted to A_i_new ∪ {xi}
      new_mb <- get_mb(xi, A_list[[xi]], dataAll, alpha = alpha)
      if (!setequal(new_mb, old_mb)) converged <- FALSE
      mb_list[[xi]] <- new_mb
    }
    
    if (verbose) {
      cat("\n==============================\n")
      cat("Iteration", iter, ": A_i and MBs\n")
      cat("==============================\n")
      for (v in xlabs) {
        cat("A_", v, " = {", paste(A_list[[v]], collapse = ", "), "}\n", sep = "")
        cat("MB(", v, ") = { ", paste(mb_list[[v]], collapse = ", "), " }\n", sep = "")
      }
      cat("==============================\n\n")
    }
    
  } # end while !converged
  
  # Final permutation and normalization
  adj_new <- permute_to_lower(adj_new)
  adj_new <- normalize_to_sorted_ancestral(adj_new)
  return(adj_new)
}

# ==========================================================
# 5. Example run 
# ==========================================================

# Simulate one dataset 
sim_obj <- sim_dat(
  n = 1000,        # sample size
  d_z = 10,       # number of z
  d_x = 4,         # number of observed variables
  rho = 0.8,
  r2 = 0.8,        # variance explained by instruments
  lin_pr = 1,      # linearity
  sp = 0.5,       # sparsity
  method = 'er',   # Erdős–Rényi structure
  pref = 1
)

# Run ascend algorithm
amat_ascend <- ascend_fn(sim_obj)

# Ground truth adjacency matrix
amat_true <- normalize_to_sorted_ancestral(sim_obj$adj_mat )  

# Print results
cat("Ground truth adjacency matrix (amat_true):\n")
print(amat_true)
cat("\nEstimated adjacency matrix from ASCEND (amat_ascend):\n")
print(amat_ascend)
