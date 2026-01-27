library(data.table)
library(matrixStats)

# -----------------------------------------------------------------
# Utility: adjacency and DAG helpers
# -----------------------------------------------------------------
clean_matrix_adj <- function(mat) {
  mat2 <- as.matrix(mat)
  mat2[is.na(mat2)] <- 0
  mat2[mat2 != 0] <- 1
  return(mat2)
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
    v <- zeros[1]; zeros <- zeros[-1]
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

sort_ancestral_matrix <- function(adj_matrix) {
  if (is.null(adj_matrix) || nrow(adj_matrix) == 0) return(adj_matrix)
  support <- adj_matrix
  support[is.na(support)] <- 0
  support[!(support %in% c(0.5, 1))] <- 0
  support[support != 0] <- 1
  
  ord <- topo_sort_adj(support)
  if (is.null(ord)) return(adj_matrix)
  
  sorted <- adj_matrix[ord, ord, drop = FALSE]
  diag(sorted) <- NA
  sorted[lower.tri(sorted)] <- 0
  return(sorted)
}

# -----------------------------------------------------------------
# Simulator
# -----------------------------------------------------------------
sim_dat <- function(n, d_z, d_x,
                    rho = 0.5, r2 = 0.5, lin_pr = 1,
                    sp = 0.9, method = "er", pref = NA,
                    p_cross = 0.05, p_z = NULL, p_x = NULL,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  prep <- function(dat, pr) {
    dat <- as.matrix(dat)
    p <- ncol(dat)
    if (p == 0 || pr >= 1) return(dat)
    n_nl <- round((1 - pr) * p)
    if (n_nl <= 0) return(dat)
    cols <- sample.int(p, size = n_nl)
    nl_types <- sample(c('sq', 'sqrt', 'sftpls', 'relu'), size = n_nl, replace = TRUE)
    out <- dat
    for (ii in seq_along(cols)) {
      idx <- cols[ii]; ttype <- nl_types[ii]; v <- dat[, idx]
      out[, idx] <- switch(ttype,
                           sq = v^2,
                           sqrt = sqrt(abs(v)),
                           sftpls = log1p(exp(v)),
                           relu = pmax(v, 0))
    }
    return(out)
  }
  
  sim_noise <- function(signal, r2) {
    signal <- as.numeric(signal)
    var_mu <- var(signal, na.rm = TRUE)
    if (is.na(var_mu) || var_mu == 0) return(rnorm(n, sd = 0))
    var_noise <- if (r2 > 0) (var_mu - r2 * var_mu) / r2 else var_mu
    noise <- rnorm(n, sd = sqrt(pmax(0, var_noise)))
    return(noise)
  }
  
  make_random_dag_upper <- function(k, p_edge = NULL) {
    if (is.null(p_edge)) p_edge <- min(max(1 - sp, 0.01), 0.4)
    A <- matrix(0L, k, k)
    if (k > 1) {
      for (i in 1:(k-1)) for (j in (i+1):k) if (runif(1) < p_edge) A[i,j] <- 1L
    }
    perm <- sample.int(k)
    A_perm <- A[perm, perm, drop = FALSE]
    rownames(A_perm) <- colnames(A_perm) <- paste0("n", seq_len(k))
    list(A = A_perm, perm = perm)
  }
  
  # --- Build DAGs ---
  Azz <- if (d_z > 0) {
    z_dag <- make_random_dag_upper(d_z, p_edge = p_z)
    rownames(z_dag$A) <- colnames(z_dag$A) <- paste0("z", seq_len(d_z))
    z_dag$A
  } else matrix(0L, 0, 0)
  
  x_dag <- make_random_dag_upper(d_x, p_edge = p_x)
  Axx <- x_dag$A
  rownames(Axx) <- colnames(Axx) <- paste0("x", seq_len(d_x))
  
  Azx <- matrix(0L, d_z, d_x, dimnames = list(paste0("z", seq_len(d_z)), paste0("x", seq_len(d_x))))
  if (d_z > 0 && d_x > 0) {
    for (i in seq_len(d_z)) for (j in seq_len(d_x)) if (runif(1) < p_cross) Azx[i,j] <- 1L
  }
  
  full_nodes <- c(if(d_z>0) paste0("z", seq_len(d_z)) else character(0), paste0("x", seq_len(d_x)))
  A_full <- matrix(0L, length(full_nodes), length(full_nodes),
                   dimnames = list(full_nodes, full_nodes))
  if(d_z>0) A_full[1:d_z,1:d_z] <- Azz
  if(d_z>0 && d_x>0) A_full[1:d_z,(d_z+1):(d_z+d_x)] <- Azx
  if(d_x>0) A_full[(d_z+1):(d_z+d_x),(d_z+1):(d_z+d_x)] <- Axx
  
  # --- Simulate Z ---
  Zmat <- if (d_z>0) {
    mat <- matrix(0, n, d_z, dimnames = list(NULL, paste0("z", seq_len(d_z))))
    topoZ <- as.integer(igraph::topo_sort(igraph::graph_from_adjacency_matrix(Azz, mode="directed"), mode="out"))
    for(idx in topoZ) {
      parents_idx <- which(Azz[,idx]==1L)
      signal <- if(length(parents_idx)>0) as.numeric(prep(mat[,parents_idx, drop=FALSE], lin_pr) %*% rnorm(length(parents_idx))) else rep(0,n)
      mat[,idx] <- signal + rnorm(n)
    }
    mat
  } else matrix(nrow=n, ncol=0)
  
  # --- Simulate X ---
  Xmat <- matrix(0, n, d_x, dimnames = list(NULL, paste0("x", seq_len(d_x))))
  adj_xx_out <- matrix(0L, d_x, d_x, dimnames=list(paste0("x", seq_len(d_x)), paste0("x", seq_len(d_x))))
  if(d_x>0) {
    topoX <- as.integer(igraph::topo_sort(igraph::graph_from_adjacency_matrix(Axx, mode="directed"), mode="out"))
    for(j_pos in seq_along(topoX)) {
      j <- topoX[j_pos]
      col_full <- d_z + j
      z_par_idx <- if(d_z>0) which(A_full[1:d_z, col_full]==1L) else integer(0)
      signal_z <- if(length(z_par_idx)>0) as.numeric(prep(Zmat[, z_par_idx, drop=FALSE], lin_pr) %*% sample(c(1,-1), length(z_par_idx), replace=TRUE)) else 0
      x_par_idx <- which(A_full[(d_z+1):(d_z+d_x), col_full]==1L)
      already_sim <- topoX[1:(j_pos-1)]
      if(length(already_sim)==0) x_par_idx <- integer(0) else x_par_idx <- intersect(x_par_idx, already_sim)
      signal_x <- if(length(x_par_idx)>0) {
        PaX <- prep(Xmat[, x_par_idx, drop=FALSE], lin_pr)
        adj_xx_out[j, x_par_idx] <- 1L
        beta_x <- rep(1, ncol(PaX))
        as.numeric(PaX %*% beta_x)
      } else 0
      Xmat[,j] <- signal_z + signal_x + sim_noise(signal_z+signal_x, r2)
    }
  }
  
  dat_dt <- data.table()
  if(d_z>0) dat_dt <- cbind(dat_dt, as.data.table(Zmat))
  if(d_x>0) dat_dt <- cbind(dat_dt, as.data.table(Xmat))
  
  # Normalize adjacency for X->X
  adj_xx_out_norm <- sort_ancestral_matrix(adj_xx_out)
  
  params <- list(n=n, d_z=d_z, d_x=d_x, rho=rho, r2=r2, lin_pr=lin_pr, sp=sp,
                 method=method, pref=pref, p_cross=p_cross, p_z=p_z, p_x=p_x)
  
  return(list(dat=dat_dt, adj_full=A_full, adj_mat=adj_xx_out_norm, params=params))
}

# -------------------------------
# Example usage
# -------------------------------
sim <- sim_dat(n=2000, d_z=10, d_x=4, r2=0.5, lin_pr=0.7, sp=0.5, p_cross=0.08, seed=123)

# Basic checks
A <- sim$adj_full; d_z <- sim$params$d_z; d_x <- sim$params$d_x
cat("Dimensions full adjacency:", dim(A), "\n")
cat("Sum of X->Z edges:", sum(A[(d_z+1):(d_z+d_x), 1:d_z]), "\n")
cat("Z->X edges count:", sum(A[1:d_z, (d_z+1):(d_z+d_x)]), "\n")
cat("X->X direct adj:\n"); print(sim$adj_mat)

# Plot
library(igraph)
g <- graph_from_adjacency_matrix(A, mode = "directed")
V(g)$color <- ifelse(grepl("^z", V(g)$name), "skyblue", "salmon")
plot(g, layout = layout_as_tree, vertex.size = 20, edge.arrow.size = 0.4)

