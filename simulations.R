## Two-tier linear-Gaussian simulation for benchmarking ancestral
## causal-discovery methods.
##
## sim_dat() generates synthetic data from a directed acyclic graph
## with two strata of variables:
##
##   Z : background (e.g. genotype). Edges may exist Z -> Z and Z -> X.
##   X : foreground (e.g. transcriptome). Edges may exist X -> X.
##
## No X -> Z edges are generated; this enforces the two-tier
## assumption that all Z causally precede all X.
##
## Each variable is drawn from a linear structural equation with
## optional element-wise nonlinear transforms (controlled by lin_pr),
## additive Gaussian noise scaled to achieve a target r^2 on each
## node, and signed random coefficients on parents.
##
## Returns a list with:
##   $dat       data.table of generated samples, columns z1..z_dz,
##              x1..x_dx
##   $adj_full  full (d_z + d_x) x (d_z + d_x) direct-edge adjacency
##   $adj_xx    d_x x d_x direct X -> X adjacency, where
##              adj_xx[child, parent] = 1
##   $params    list of parameters used

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(matrixStats)
})

sim_dat <- function(n,
                    d_z,
                    d_x,
                    rho      = 0.5,
                    r2       = 0.5,
                    lin_pr   = 1,
                    sp       = 0.9,
                    method   = "er",
                    pref     = NA,
                    p_cross  = 0.05,
                    p_z      = NULL,
                    p_x      = NULL,
                    x_effect = 0.8,
                    seed     = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Apply element-wise nonlinear transforms to a random subset of columns,
  # controlled by lin_pr (1 = fully linear, 0 = fully nonlinear).
  prep <- function(dat, pr) {
    if (is.null(dim(dat))) dat <- matrix(dat, ncol = 1) else dat <- as.matrix(dat)
    p <- ncol(dat); if (p == 0 || pr >= 1) return(dat)
    n_nl <- round((1 - pr) * p); if (n_nl <= 0) return(dat)
    cols <- sample.int(p, size = n_nl)
    nl_types <- sample(c('sq', 'sqrt', 'sftpls', 'relu'), n_nl, replace = TRUE)
    out <- dat
    for (ii in seq_along(cols)) {
      idx <- cols[ii]; v <- dat[, idx]
      if      (nl_types[ii] == 'sq')     out[, idx] <- v^2
      else if (nl_types[ii] == 'sqrt')   out[, idx] <- sqrt(abs(v))
      else if (nl_types[ii] == 'sftpls') out[, idx] <- log1p(exp(v))
      else if (nl_types[ii] == 'relu')   out[, idx] <- pmax(v, 0)
    }
    out
  }
  
  # Generate additive noise calibrated to a target signal-to-total
  # variance ratio of r2.
  sim_noise <- function(signal, r2) {
    signal <- as.numeric(signal); var_mu <- var(signal, na.rm = TRUE)
    if (is.na(var_mu) || var_mu == 0) return(rnorm(n, sd = 1))
    if (r2 <= 0) return(rnorm(n, sd = sqrt(var_mu)))
    rnorm(n, sd = sqrt(pmax(0, (var_mu - r2 * var_mu) / r2)))
  }
  
  # Random upper-triangular DAG on k nodes with edge probability p_edge,
  # then permute to a random node order. sp controls density (higher sp
  # means sparser).
  make_dag_upper <- function(k, p_edge = NULL) {
    if (is.null(p_edge)) p_edge <- min(max(1 - sp, 0.01), 0.4)
    A <- matrix(0, k, k)
    for (i in 1:(k - 1)) for (j in (i + 1):k) if (runif(1) < p_edge) A[i, j] <- 1
    perm <- sample.int(k); list(A = A[perm, perm], perm = perm)
  }
  
  # Build the three blocks of the full adjacency matrix.
  z_dag <- make_dag_upper(d_z, p_z); Azz <- z_dag$A
  rownames(Azz) <- colnames(Azz) <- paste0("z", seq_len(d_z))
  x_dag <- make_dag_upper(d_x, p_x); Axx <- x_dag$A
  rownames(Axx) <- colnames(Axx) <- paste0("x", seq_len(d_x))
  
  Azx <- matrix(0, d_z, d_x,
                dimnames = list(paste0("z", 1:d_z), paste0("x", 1:d_x)))
  for (i in seq_len(d_z)) for (j in seq_len(d_x))
    if (runif(1) < p_cross) Azx[i, j] <- 1
  
  fn <- c(rownames(Azz), colnames(Axx))
  A_full <- matrix(0, d_z + d_x, d_z + d_x, dimnames = list(fn, fn))
  A_full[1:d_z, 1:d_z] <- Azz
  A_full[1:d_z, (d_z + 1):(d_z + d_x)] <- Azx
  A_full[(d_z + 1):(d_z + d_x), (d_z + 1):(d_z + d_x)] <- Axx
  
  # Sample Z variables in topological order.
  gZ <- igraph::graph_from_adjacency_matrix(Azz, mode = "directed")
  topoZ <- as.integer(igraph::topo_sort(gZ, mode = "out"))
  Zmat <- matrix(0, n, d_z,
                 dimnames = list(NULL, paste0("z", seq_len(d_z))))
  for (idx in topoZ) {
    par <- which(Azz[, idx] == 1)
    if (!length(par)) {
      Zmat[, idx] <- rnorm(n)
    } else {
      Pa <- prep(Zmat[, par, drop = FALSE], lin_pr)
      Zmat[, idx] <- as.numeric(Pa %*% rnorm(ncol(Pa), sd = 1)) +
        rnorm(n, sd = 1)
    }
  }
  
  # Sample X variables in topological order, mixing Z-parent and
  # X-parent contributions with a target signal-to-noise via sim_noise.
  Xmat <- matrix(0, n, d_x,
                 dimnames = list(NULL, paste0("x", seq_len(d_x))))
  adj_xx_out <- matrix(0, d_x, d_x,
                       dimnames = list(paste0("x", seq_len(d_x)),
                                       paste0("x", seq_len(d_x))))
  diag(adj_xx_out) <- NA_real_
  
  gX <- igraph::graph_from_adjacency_matrix(Axx, mode = "directed")
  topoX <- as.integer(igraph::topo_sort(gX, mode = "out"))
  
  for (j_idx in seq_along(topoX)) {
    j <- topoX[j_idx]; cif <- d_z + j
    
    zp <- which(A_full[1:d_z, cif] == 1)
    sig_z <- rep(0, n)
    if (length(zp)) {
      PaZ <- prep(Zmat[, zp, drop = FALSE], lin_pr)
      sig_z <- as.numeric(PaZ %*% sample(c(1, -1), ncol(PaZ), replace = TRUE))
    }
    
    xp <- intersect(which(A_full[(d_z + 1):(d_z + d_x), cif] == 1),
                    if (j_idx > 1) topoX[1:(j_idx - 1)] else integer(0))
    sig_x <- rep(0, n)
    if (length(xp)) {
      PaX <- prep(Xmat[, xp, drop = FALSE], lin_pr)
      adj_xx_out[j, xp] <- 1
      csd <- matrixStats::colSds(PaX); csd[csd == 0] <- 1
      beta_x <- (x_effect / csd) * sample(c(1, -1), length(csd), replace = TRUE)
      sig_x <- as.numeric(PaX %*% beta_x)
    }
    
    sig_total <- sig_z + sig_x
    if (var(sig_total) < 1e-6) sig_total <- rnorm(n)
    Xmat[, j] <- sig_total + sim_noise(sig_total, r2)
  }
  
  dat_dt <- data.table(Zmat, Xmat)
  colnames(dat_dt) <- c(paste0("z", 1:d_z), paste0("x", 1:d_x))
  list(dat = dat_dt, adj_full = A_full, adj_xx = adj_xx_out,
       params = list(n = n, d_z = d_z, d_x = d_x, r2 = r2, lin_pr = lin_pr,
                     sp = sp, p_cross = p_cross, x_effect = x_effect))
}