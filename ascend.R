# ======================================================================
# Load required packages
# ======================================================================
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
# Helper functions 
# ======================================================================

clean_matrix_adj <- function(mat) {
  mat2 <- as.matrix(mat)
  mat2[is.na(mat2)] <- 0
  mat2[mat2 != 0] <- 1
  return(mat2)
}

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

topo_sort_adj <- function(adj) {
  adj <- clean_matrix_adj(adj)
  n <- nrow(adj)
  if (n == 0) return(integer(0))
  A <- (adj != 0) * 1
  indeg <- colSums(A)
  indeg[is.na(indeg)] <- 0
  
  order <- integer(0)
  zeros <- sort(which(indeg == 0))
  
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
  
  if (length(order) != n) return(NULL)
  return(order)
}

is_dag_adj <- function(adj) {
  !is.null(topo_sort_adj(adj))
}

sort_ancestral_matrix <- function(adj_matrix) {
  if (is.null(adj_matrix) || nrow(adj_matrix) == 0) return(adj_matrix)
  
  support <- adj_matrix
  support[is.na(support)] <- 0
  support[!(support %in% c(0.5, 1))] <- 0
  support[support != 0] <- 1
  
  ord <- topo_sort_adj(support)
  if (is.null(ord)) {
    warning("Graph is not a DAG; returning original matrix")
    return(adj_matrix)
  }
  
  sorted <- adj_matrix[ord, ord, drop = FALSE]
  diag(sorted) <- NA
  sorted[lower.tri(sorted)] <- 0
  return(sorted)
}

build_constraint_adj <- function(m) {
  n <- nrow(m)
  A <- matrix(0, n, n)
  rownames(A) <- rownames(m)
  colnames(A) <- colnames(m)
  for (j in seq_len(n)) {
    for (i in seq_len(n)) {
      if (i == j) next
      if (is.na(m[j, i])) next
      if (m[j, i] == 1 || m[j, i] == 0.5) A[i, j] <- 1
    }
  }
  return(A)
}

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

test_r1r2 <- function(W, xj, S_w, xi, data, alpha = 0.05) {
  vars <- unique(c(xj, W, S_w, xi))
  miss <- setdiff(vars, colnames(data))
  if (length(miss)) return(list(beta_0 = NA, beta_1 = NA, pval_0 = NA, pval_1 = NA,
                                p_value_diff = NA, r1 = FALSE, r2 = FALSE))
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

# ======================================================================
# ASCEND function 
# ======================================================================

ascend_fn <- function(sim_obj, maxiter = 9, alpha = 0.05, fdr_correction = TRUE) {
  dat <- sim_obj$dat
  z <- as.matrix(dplyr::select(dat, starts_with('z')))
  x <- as.matrix(dplyr::select(dat, starts_with('x')))
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  
  robust_scale_safe <- function(x) {
    if (all(is.na(x))) return(x)
    med <- median(x, na.rm = TRUE)
    mad_val <- mad(x, na.rm = TRUE)
    if (mad_val == 0 || is.na(mad_val)) {
      sd_val <- sd(x, na.rm = TRUE)
      if (sd_val == 0 || is.na(sd_val)) {
        x <- x + rnorm(length(x), 0, 1e-6)
        sd_val <- sd(x, na.rm = TRUE)
      }
      return((x - mean(x, na.rm = TRUE)) / sd_val)
    }
    return((x - med) / mad_val)
  }
  
  dataAll <- as.data.frame(dat)
  for (col in colnames(dataAll)) dataAll[[col]] <- robust_scale_safe(dataAll[[col]])
  dataAll[is.infinite(as.matrix(dataAll))] <- 0
  dataAll[is.na(dataAll)] <- 0
  
  z_cols <- grep("^z", colnames(dataAll), value = TRUE)
  x_cols <- grep("^x", colnames(dataAll), value = TRUE)
  mb_list <- vector("list", length = d_x); names(mb_list) <- xlabs
  
  set.seed(42)
  for (i in seq_len(d_x)) {
    node_i <- xlabs[i]
    node_sd <- sd(dataAll[[node_i]], na.rm = TRUE)
    if (is.na(node_sd) || node_sd == 0) {
      mb_list[[node_i]] <- character(0)
      next
    }
    sub_df <- dataAll[, c(z_cols, node_i), drop = FALSE]
    col_vars <- sapply(sub_df, function(col) sd(col, na.rm = TRUE))
    if (any(col_vars == 0)) sub_df <- sub_df[, col_vars > 0, drop = FALSE]
    if (ncol(sub_df) < 2) {
      mb_list[[node_i]] <- character(0)
      next
    }
    mb_list[[node_i]] <- tryCatch({
      learn.mb(sub_df, node = node_i, method = "iamb", test = "zf", alpha = alpha)
    }, error = function(e) character(0))
  }
  
  adj_new <- matrix(NA, d_x, d_x); diag(adj_new) <- NA
  rownames(adj_new) <- colnames(adj_new) <- xlabs
  
  converged <- FALSE; iter <- 0
  while (!converged && iter <= maxiter) {
    iter <- iter + 1
    converged <- TRUE
    pval_info <- list(); pval_index <- 1
    
    for (i in 2:d_x) {
      for (j in 1:(i-1)) {
        if (is.na(adj_new[i, j])) {
          xi <- xlabs[i]; xj <- xlabs[j]
          mb_i <- mb_list[[xi]]; mb_j <- mb_list[[xj]]
          S <- setdiff(union(mb_i, mb_j), c(xi, xj))
          S_data <- dataAll[, intersect(c(xi, xj, S), colnames(dataAll)), drop = FALSE]
          
          keep <- sapply(S_data, function(col) {
            if(all(is.na(col))) return(FALSE)
            vals <- unique(na.omit(col))
            length(vals) > 1
          })
          keep <- as.logical(keep)
          if (length(keep) == ncol(S_data)) S_data <- S_data[, keep, drop = FALSE]
          
          # Sync S with remaining columns
          S <- intersect(S, colnames(S_data))
          
          if (!(xi %in% colnames(S_data)) || !(xj %in% colnames(S_data))) { converged <- FALSE; next }
          if (nrow(S_data) < 10) { converged <- FALSE; next }
          
          # ---------- SAFE FORMULA CONSTRUCTION ----------
          if (length(S) == 0) {
            fit_null <- lm(as.formula(paste(xj, "~ 1")), data = S_data)
            fit_alt  <- lm(as.formula(paste(xj, "~", xi)), data = S_data)
          } else {
            rhs0 <- paste(S, collapse = " + ")
            rhs1 <- paste(c(S, xi), collapse = " + ")
            fit_null <- lm(as.formula(paste(xj, "~", rhs0)), data = S_data)
            fit_alt  <- lm(as.formula(paste(xj, "~", rhs1)), data = S_data)
          }
          
          test_res <- tryCatch(anova(fit_null, fit_alt), error = function(e) NULL)
          pval <- if(!is.null(test_res)) test_res$`Pr(>F)`[2] else NA
          if (!is.na(pval)) {
            pval_info[[pval_index]] <- list(i=i,j=j,xi=xi,xj=xj,pval=pval,S=S,S_data=S_data)
            pval_index <- pval_index + 1
          }
        }
      }
    }
    
    if (length(pval_info) > 0) {
      pvals <- sapply(pval_info, function(x) x$pval)
      adj_pvals <- if(fdr_correction) p.adjust(pvals, method = "BH") else pvals
      
      for (idx in seq_along(pval_info)) {
        info <- pval_info[[idx]]; i <- info$i; j <- info$j
        xi <- info$xi; xj <- info$xj; S <- info$S; S_data <- info$S_data
        adj_pval <- adj_pvals[idx]
        
        if (adj_pval > alpha) {
          adj_new[i,j] <- 0; adj_new[j,i] <- 0
          converged <- FALSE
        } else {
          for (W in S) {
            Sw <- setdiff(S, W)
            act <- test_r1r2(W,xj,Sw,xi,dataAll,alpha)
            deact <- test_r1r2(W,xi,Sw,xj,dataAll,alpha)
            result <- NULL
            if (!is.null(deact) && is.list(deact) && deact$r1) result <- try_edge_update(adj_new,i,j,1,0,iter)
            else if (!is.null(act) && is.list(act) && act$r1) result <- try_edge_update(adj_new,i,j,0,1,iter)
            else if (!is.null(deact) && is.list(deact) && deact$r2) result <- try_edge_update(adj_new,i,j,0,0.5,iter)
            else if (!is.null(act) && is.list(act) && act$r2) result <- try_edge_update(adj_new,i,j,0.5,0,iter)
            if (!is.null(result) && result$success) { adj_new <- result$adj; converged <- FALSE; break }
          }
        }
      }
    }
    
    bin_current <- (clean_matrix_adj(adj_new) == 1 | clean_matrix_adj(adj_new) == 0.5) * 1
    if (!is_dag_adj(bin_current)) stop("Cycle detected")
    # closure (keep upper-triangular)
    closure <- FALSE
    while(!closure) {
      closure <- TRUE
      for(i in seq_len(ncol(adj_new))) {
        m <- which(adj_new[,i]==1)
        if(length(m)>0) {
          submat <- adj_new[, m, drop=FALSE]
          if(any(submat==1, na.rm=TRUE)) {
            e_idx <- unique(which(submat==1, arr.ind=TRUE)[,1])
            for(row in e_idx) {
              if(is.na(adj_new[row,i]) || adj_new[row,i]!=1) { adj_new[row,i]<-1; closure<-FALSE }
            }
          }
        }
      }
    }
    
    for(i in seq_len(d_x)) {
      node_i <- xlabs[i]; old_mb <- mb_list[[node_i]]
      non_desc <- which((adj_new[,i]==1) | (adj_new[,i]==0.5))
      A_i <- intersect(union(old_mb, xlabs[non_desc]), colnames(dataAll))
      sub_df2 <- dataAll[, intersect(c(A_i, node_i), colnames(dataAll)), drop=FALSE]
      keep2 <- sapply(sub_df2,function(col){ vals<-unique(na.omit(col)); length(vals)>1 })
      keep2 <- as.logical(keep2)
      if(length(keep2)==ncol(sub_df2)) sub_df2 <- sub_df2[,keep2,drop=FALSE]
      if(ncol(sub_df2)>=2 && node_i %in% colnames(sub_df2)) {
        new_mb <- tryCatch({ learn.mb(sub_df2, node=node_i, method="iamb", test="zf", alpha=alpha) }, error=function(e) old_mb)
      } else { new_mb <- old_mb }
      if(!setequal(new_mb, old_mb)) converged<-FALSE
      mb_list[[node_i]] <- new_mb
    }
  }
  
  if(!all(is.na(adj_new))) {
    tryCatch({
      adj_new <- permute_to_lower(adj_new)
      adj_new <- sort_ancestral_matrix(adj_new)
    }, error=function(e){ warning("Final normalization failed:", e$message) })
  }
  return(adj_new)
}

# ======================================================================
# Example run
# ======================================================================

sim_obj <- sim_dat(
  n = 200,
  d_z = 30,
  d_x = 7,
  rho = 0.5,
  r2 = 0.7,
  lin_pr = 1,
  sp = 0.5,
  method = 'er',
  pref = NA
)

amat_ascend <- ascend_fn(sim_obj)
amat_true <- sim_obj$adj_mat

cat("Ground truth adjacency matrix (amat_true):\n")
print(amat_true)
cat("\nEstimated adjacency matrix from ASCEND (amat_ascend):\n")
print(amat_ascend)
