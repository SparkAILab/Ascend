library(data.table)
library(matrixStats)
library(igraph)

sim_dat <- function(n,
                                 d_z,
                                 d_x,
                                 rho = 0.5,
                                 r2 = 0.5,
                                 lin_pr = 1,
                                 sp = 0.9,
                                 method = "er",
                                 pref = NA,
                                 p_cross = 0.05,
                                 p_z = NULL,
                                 p_x = NULL,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # -----------------------------------
  # Helpers (kept close to original)
  # -----------------------------------
  prep <- function(dat, pr) {
    if (is.null(dim(dat))) dat <- matrix(dat, ncol = 1) else dat <- as.matrix(dat)
    p <- ncol(dat)
    if (p == 0 || pr >= 1) return(dat)
    n_nl <- round((1 - pr) * p)
    if (n_nl <= 0) return(dat)
    cols <- sample.int(p, size = n_nl)
    nl_types <- sample(c('sq', 'sqrt', 'sftpls', 'relu'), size = n_nl, replace = TRUE)
    out <- dat
    for (ii in seq_along(cols)) {
      idx <- cols[ii]; ttype <- nl_types[ii]; v <- dat[, idx]
      if (ttype == 'sq') out[, idx] <- v^2
      else if (ttype == 'sqrt') out[, idx] <- sqrt(abs(v))
      else if (ttype == 'sftpls') out[, idx] <- log1p(exp(v))
      else if (ttype == 'relu') out[, idx] <- pmax(v, 0)
    }
    return(out)
  }
  sim_noise <- function(signal, r2) {
    signal <- as.numeric(signal)
    var_mu <- stats::var(signal, na.rm = TRUE)
    if (is.na(var_mu) || var_mu == 0) return(rnorm(n, sd = 0))
    if (r2 <= 0) return(rnorm(n, sd = sqrt(var_mu)))
    var_noise <- (var_mu - r2 * var_mu) / r2
    var_noise <- pmax(0, var_noise)
    noise <- rnorm(n, sd = sqrt(var_noise))
    return(noise)
  }
  # -----------------------------------
  # Build internal random DAGs (upper triangular construction ensures acyclicity)
  # -----------------------------------
  make_random_dag_upper <- function(k, p_edge = NULL) {
    if (is.null(p_edge)) {
      # default density derived from sp (higher sp -> sparser). convert sp to p_edge roughly
      p_edge <- 1 - sp
      p_edge <- min(max(p_edge, 0.01), 0.4)
    }
    A <- matrix(0, nrow = k, ncol = k)
    # allow edges only i -> j for i < j (upper triangular) to ensure acyclic
    for (i in 1:(k - 1)) {
      for (j in (i + 1):k) {
        if (runif(1) < p_edge) A[i, j] <- 1
      }
    }
    # optionally permute node order to randomize topology while preserving acyclicity:
    perm <- sample.int(k)
    A_perm <- A[perm, perm]
    # reorder back to a canonical labeling (so A_perm is a DAG but nodes are permuted)
    return(list(A = A_perm, perm = perm))
  }
  # -----------------------------------
  # 1) Create Z-DAG and X-DAG (internal)
  # -----------------------------------
  z_dag <- make_random_dag_upper(d_z, p_edge = p_z)
  Azz <- z_dag$A
  # label rows/cols z1..z_d
  rownames(Azz) <- colnames(Azz) <- paste0("z", seq_len(d_z))
  
  x_dag <- make_random_dag_upper(d_x, p_edge = p_x)
  Axx <- x_dag$A
  rownames(Axx) <- colnames(Axx) <- paste0("x", seq_len(d_x))
  
  # -----------------------------------
  # 2) Cross-layer Z -> X edges sampled independently with probability p_cross
  #    (no X->Z edges by construction)
  # -----------------------------------
  Azx <- matrix(0, nrow = d_z, ncol = d_x, dimnames = list(paste0("z", 1:d_z), paste0("x", 1:d_x)))
  for (i in seq_len(d_z)) {
    for (j in seq_len(d_x)) {
      if (runif(1) < p_cross) Azx[i, j] <- 1
    }
  }
  
  # -----------------------------------
  # 3) Full adjacency (Z then X)
  # -----------------------------------
  full_nodes <- c(rownames(Azz), colnames(Axx))
  A_full <- matrix(0, nrow = d_z + d_x, ncol = d_z + d_x, dimnames = list(full_nodes, full_nodes))
  A_full[1:d_z, 1:d_z] <- Azz
  A_full[1:d_z, (d_z + 1):(d_z + d_x)] <- Azx
  A_full[(d_z + 1):(d_z + d_x), (d_z + 1):(d_z + d_x)] <- Axx
  # bottom-left block remains zero -> guarantees no X->Z
  
  # -----------------------------------
  # 4) Simulate Z recursively according to Z-DAG (topological order)
  # -----------------------------------
  # find topological order for the Z-DAG (Azz is permuted; use igraph to topo sort)
  gZ <- igraph::graph_from_adjacency_matrix(Azz, mode = "directed")
  topoZ <- igraph::topo_sort(gZ, mode = "out")
  topoZ <- as.integer(topoZ) # indices in 1:d_z referencing rows/cols of Azz
  # create Z matrix
  Zmat <- matrix(0, nrow = n, ncol = d_z, dimnames = list(NULL, paste0("z", seq_len(d_z))))
  for (idx in topoZ) {
    parents_idx <- which(Azz[, idx] == 1)
    if (length(parents_idx) == 0) {
      # exogenous
      Zmat[, idx] <- rnorm(n)
    } else {
      Pa <- prep(Zmat[, parents_idx, drop = FALSE], lin_pr)
      # sample Rademacher or normalized weights
      if (ncol(Pa) > 0) {
        betas <- rnorm(ncol(Pa), sd = 1)
        signal <- as.numeric(Pa %*% betas)
      } else signal <- rep(0, n)
      Zmat[, idx] <- signal + rnorm(n, sd = 1)
    }
  }
  # apply a mild correlation scaling if rho>0 to preserve earlier behaviour (optional)
  if (!is.null(rho) && rho > 0) {
    # small mixing with Toeplitz structure to preserve some marginal correlation if desired
    # but keep causal structure primary; this step is optional and commented out
    # Sigma <- toeplitz(rho^(0:(d_z - 1)))
    # Zmat <- scale(Zmat %*% chol(Sigma))
    NULL
  }
  
  # -----------------------------------
  # 5) Simulate X recursively using A_full parents (Z->X and X->X)
  # -----------------------------------
  Xmat <- matrix(0, nrow = n, ncol = d_x, dimnames = list(NULL, paste0("x", seq_len(d_x))))
  adj_xx_out <- matrix(0, nrow = d_x, ncol = d_x, dimnames = list(paste0("x", seq_len(d_x)), paste0("x", seq_len(d_x))))
  diag(adj_xx_out) <- NA_real_
  
  # find topo order within X-subgraph
  gX <- igraph::graph_from_adjacency_matrix(Axx, mode = "directed")
  topoX <- as.integer(igraph::topo_sort(gX, mode = "out"))
  # topoX gives indices 1..d_x in an order where parents come before children
  # iterate through topoX order and simulate each x
  for (j_idx in seq_along(topoX)) {
    j <- topoX[j_idx] # j is index among 1..d_x
    col_index_full <- d_z + j
    xname <- paste0("x", j)
    # Z parents
    z_par_idx <- which(A_full[1:d_z, col_index_full] == 1)
    if (length(z_par_idx) > 0) {
      PaZ <- prep(Zmat[, z_par_idx, drop = FALSE], lin_pr)
      if (ncol(PaZ) > 0) {
        beta_z <- sample(c(1, -1), size = ncol(PaZ), replace = TRUE)
        signal_z <- as.numeric(PaZ %*% beta_z)
      } else signal_z <- rep(0, n)
    } else {
      signal_z <- rep(0, n)
    }
    # X parents (only those among X that are already simulated according to topoX order and that have edge -> j)
    x_par_idx_full <- which(A_full[(d_z + 1):(d_z + d_x), col_index_full] == 1)
    # convert to 1..d_x index
    x_par_idx <- x_par_idx_full
    # keep only those parents that have been simulated already (i.e., appear earlier in topoX)
    if (length(x_par_idx) > 0) {
      already_simulated <- topoX[1:(j_idx - 1)]
      x_par_idx <- intersect(x_par_idx, already_simulated)
    }
    signal_x <- rep(0, n)
    if (length(x_par_idx) > 0) {
      PaX <- prep(Xmat[, x_par_idx, drop = FALSE], lin_pr)
      if (ncol(PaX) > 0) {
        adj_xx_out[j, x_par_idx] <- 1
        causal_wt <- 1 / max(1, sqrt(length(z_par_idx) + length(x_par_idx)))
        sigma_xij <- sqrt(pmax(0, causal_wt * stats::var(signal_z, na.rm = TRUE)))
        csd <- matrixStats::colSds(PaX); csd[csd == 0] <- 1
        beta_x <- rep(sigma_xij, length(csd)) / csd
        signal_x <- as.numeric(PaX %*% beta_x)
      }
    }
    signal_total <- signal_z + signal_x
    Xmat[, j] <- signal_total + sim_noise(signal_total, r2)
  }
  
  # -----------------------------------
  # 6) return tidy list
  # -----------------------------------
  dat_dt <- data.table(Zmat, Xmat)
  colnames(dat_dt) <- c(paste0("z", 1:d_z), paste0("x", 1:d_x))
  params <- list(n = n, d_z = d_z, d_x = d_x, rho = rho, r2 = r2,
                 lin_pr = lin_pr, sp = sp, method = method, pref = pref,
                 p_cross = p_cross, p_z = p_z, p_x = p_x)
  return(list(dat = dat_dt, adj_full = A_full, adj_xx = adj_xx_out, params = params))
}



# Example usage
sim <- sim_dat(n = 2000, d_z = 10, d_x = 4,
                            rho = 0.5, r2 = 0.5, lin_pr = 1,
                            sp = 0.5, method = "er", pref = NA,
                            p_cross = 0.08, seed = 123)

# Basic checks
A <- sim$adj_full
d_z <- sim$params$d_z; d_x <- sim$params$d_x

cat("Dimensions full adjacency:", dim(A), "\n")
cat("Sum of X->Z edges (should be 0):", sum(A[(d_z + 1):(d_z + d_x), 1:d_z]), "\n")
cat("Z->X edges count:", sum(A[1:d_z, (d_z + 1):(d_z + d_x)]), "\n")
cat("X->X direct adj (adj_xx):\n"); print(sim$adj_xx)

# Plot
library(igraph)
g <- graph_from_adjacency_matrix(A, mode = "directed")
V(g)$color <- ifelse(grepl("^z", V(g)$name), "skyblue", "salmon")
plot(g, layout = layout_as_tree, vertex.size = 20, edge.arrow.size = 0.4)

