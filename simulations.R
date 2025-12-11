# Minimal safe sim_dat with defensive checks
library(data.table)
library(matrixStats)
library(igraph)

# -----------------------------------------------------------------
# Utility: transitive closure and normalization to "ancestral" form
# -----------------------------------------------------------------
transitive_closure <- function(A) {
  A <- (A != 0) * 1L
  n <- nrow(A)
  R <- A
  for (k in seq_len(n)) {
    for (i in seq_len(n)) {
      if (R[i, k] == 1L) R[i, ] <- (R[i, ] | R[k, ]) * 1L
    }
  }
  diag(R) <- 0L
  return(R)
}



# -----------------------------------------------------------------
# Simulator
# -----------------------------------------------------------------
sim_dat <- function(n, d_z,
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
  # Helpers
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
    # produce noise so that Var(signal) : Var(noise) set by r2
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
  # Make random upper-triangular DAG and permute nodes
  # -----------------------------------
  make_random_dag_upper <- function(k, p_edge = NULL) {
    if (is.null(p_edge)) {
      p_edge <- 1 - sp
      p_edge <- min(max(p_edge, 0.01), 0.4)
    }
    A <- matrix(0L, nrow = k, ncol = k)
    if (k > 1) {
      for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
          if (runif(1) < p_edge) A[i, j] <- 1L
        }
      }
    }
    perm <- sample.int(k)
    A_perm <- A[perm, perm, drop = FALSE]
    rownames(A_perm) <- colnames(A_perm) <- paste0("n", seq_len(k))
    return(list(A = A_perm, perm = perm))
  }
  # -----------------------------------
  # 1) build Z and X subgraphs
  # -----------------------------------
  if (d_z < 0 || d_x < 0) stop("d_z and d_x must be nonnegative integers")
  if (d_z == 0) {
    Azz <- matrix(0L, 0, 0)
    rownames(Azz) <- colnames(Azz) <- character(0)
  } else {
    z_dag <- make_random_dag_upper(d_z, p_edge = p_z)
    Azz <- z_dag$A
    rownames(Azz) <- colnames(Azz) <- paste0("z", seq_len(d_z))
  }
  x_dag <- make_random_dag_upper(d_x, p_edge = p_x)
  Axx <- x_dag$A
  rownames(Axx) <- colnames(Axx) <- paste0("x", seq_len(d_x))
  # -----------------------------------
  # 2) cross-layer Z->X edges
  # -----------------------------------
  Azx <- matrix(0L, nrow = d_z, ncol = d_x,
                dimnames = list(paste0("z", seq_len(d_z)), paste0("x", seq_len(d_x))))
  if (d_z > 0 && d_x > 0) {
    for (i in seq_len(d_z)) for (j in seq_len(d_x)) if (runif(1) < p_cross) Azx[i, j] <- 1L
  }
  # -----------------------------------
  # 3) full adjacency (Z then X)
  # -----------------------------------
  full_nodes <- c(if (d_z > 0) paste0("z", seq_len(d_z)) else character(0),
                  paste0("x", seq_len(d_x)))
  A_full <- matrix(0L, nrow = length(full_nodes), ncol = length(full_nodes),
                   dimnames = list(full_nodes, full_nodes))
  if (d_z > 0) A_full[1:d_z, 1:d_z] <- Azz
  if (d_z > 0 && d_x > 0) A_full[1:d_z, (d_z + 1):(d_z + d_x)] <- Azx
  if (d_x > 0) A_full[(d_z + 1):(d_z + d_x), (d_z + 1):(d_z + d_x)] <- Axx
  # -----------------------------------
  # 4) simulate Z according to Azz
  # -----------------------------------
  if (d_z > 0) {
    gZ <- igraph::graph_from_adjacency_matrix(Azz, mode = "directed")
    topoZ <- as.integer(igraph::topo_sort(gZ, mode = "out"))
    Zmat <- matrix(0, nrow = n, ncol = d_z, dimnames = list(NULL, paste0("z", seq_len(d_z))))
    for (idx in topoZ) {
      parents_idx <- which(Azz[, idx] == 1L)
      if (length(parents_idx) == 0) {
        Zmat[, idx] <- rnorm(n)
      } else {
        Pa <- prep(Zmat[, parents_idx, drop = FALSE], lin_pr)
        if (ncol(Pa) > 0) {
          betas <- rnorm(ncol(Pa), sd = 1)
          signal <- as.numeric(Pa %*% betas)
        } else signal <- rep(0, n)
        Zmat[, idx] <- signal + rnorm(n, sd = 1)
      }
    }
  } else {
    Zmat <- matrix(nrow = n, ncol = 0)
  }
  # -----------------------------------
  # 5) simulate X according to Z->X and X->X (respect topo in X-subgraph)
  # -----------------------------------
  Xmat <- matrix(0, nrow = n, ncol = d_x, dimnames = list(NULL, paste0("x", seq_len(d_x))))
  adj_xx_out <- matrix(0L, nrow = d_x, ncol = d_x, dimnames = list(paste0("x", seq_len(d_x)), paste0("x", seq_len(d_x))))
  if (d_x > 0) {
    gX <- igraph::graph_from_adjacency_matrix(Axx, mode = "directed")
    topoX <- as.integer(igraph::topo_sort(gX, mode = "out"))
    # topoX is ordering of 1:d_x where parents appear before children
    for (j_pos in seq_along(topoX)) {
      j <- topoX[j_pos]            # column index in 1..d_x
      col_full <- d_z + j         # column index in A_full
      # Z parents for this X_j (indices among 1..d_z)
      z_par_idx <- integer(0)
      if (d_z > 0) z_par_idx <- which(A_full[1:d_z, col_full] == 1L)
      signal_z <- rep(0, n)
      if (length(z_par_idx) > 0) {
        PaZ <- prep(Zmat[, z_par_idx, drop = FALSE], lin_pr)
        if (ncol(PaZ) > 0) {
          beta_z <- sample(c(1, -1), size = ncol(PaZ), replace = TRUE)
          signal_z <- as.numeric(PaZ %*% beta_z)
        }
      }
      # X parents: indices within X-submatrix (1..d_x)
      x_par_idx <- which(A_full[(d_z + 1):(d_z + d_x), col_full] == 1L)
      # keep only those X-parents that were already simulated (appear earlier in topoX)
      if (length(x_par_idx) > 0) {
        already_sim <- topoX[1:(j_pos - 1)]
        if (length(already_sim) == 0) x_par_idx <- integer(0) else x_par_idx <- intersect(x_par_idx, already_sim)
      }
      signal_x <- rep(0, n)
      if (length(x_par_idx) > 0) {
        PaX <- prep(Xmat[, x_par_idx, drop = FALSE], lin_pr)
        if (ncol(PaX) > 0) {
          adj_xx_out[j, x_par_idx] <- 1L
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
  }
  # -----------------------------------
  # 6) tidy up and return
  # -----------------------------------
  dat_dt <- data.table()
  if (d_z > 0) dat_dt <- cbind(dat_dt, as.data.table(Zmat))
  if (d_x > 0) dat_dt <- cbind(dat_dt, as.data.table(Xmat))
  # guarantee numeric columns and no NA
  for (nm in colnames(dat_dt)) dat_dt[[nm]] <- as.numeric(dat_dt[[nm]])
  # normalize adjacency for X->X to ancestral form (NA diag)
  adj_xx_out_norm <- normalize_to_sorted_ancestral(adj_xx_out)
  params <- list(n = n, d_z = d_z, d_x = d_x, rho = rho, r2 = r2,
                 lin_pr = lin_pr, sp = sp, method = method, pref = pref,
                 p_cross = p_cross, p_z = p_z, p_x = p_x)
  return(list(dat = dat_dt, adj_full = A_full, adj_mat = adj_xx_out_norm, params = params))
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
cat("X->X direct adj (adj_xx):\n"); print(sim$adj_mat)

# Plot
library(igraph)
g <- graph_from_adjacency_matrix(A, mode = "directed")
V(g)$color <- ifelse(grepl("^z", V(g)$name), "skyblue", "salmon")
plot(g, layout = layout_as_tree, vertex.size = 20, edge.arrow.size = 0.4)

