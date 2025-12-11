# -------------------------
# Minimal, consistent pipeline (IMMEDIATE-EDGE evaluation)
# -------------------------
library(data.table)
library(matrixStats)
library(ggplot2)
library(foreach)
library(doMC)
library(pcalg)
library(igraph)

registerDoMC(cores = max(1, parallel::detectCores() - 1))


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
  adj_xx_out <- normalize_to_sorted_ancestral(adj_xx_out)
  params <- list(n = n, d_z = d_z, d_x = d_x, rho = rho, r2 = r2,
                 lin_pr = lin_pr, sp = sp, method = method, pref = pref,
                 p_cross = p_cross, p_z = p_z, p_x = p_x)
  return(list(dat = dat_dt, adj_full = A_full, adj_mat = adj_xx_out, params = params))
}




# -------------------------
# 2) Small utilities for immediate-adjacency convention
# -------------------------
clean_adj_immediate <- function(A, d = NULL) {
  if (is.null(A)) {
    if (is.null(d)) stop("need dimension")
    M <- matrix(0L, d, d); return(M)
  }
  A <- as.matrix(A)
  A[is.na(A)] <- 0
  mode(A) <- "integer"
  A <- (A != 0L) * 1L
  if (!is.null(d)) {
    if (nrow(A) != d || ncol(A) != d) {
      B <- matrix(0L, d, d)
      rn <- rownames(A); cn <- colnames(A)
      if (!is.null(rn) && !is.null(cn)) {
        rIx <- pmin(length(rn), d); cIx <- pmin(length(cn), d)
        # best-effort: try to place by names if they look like x1..xk
        # else just return zero-padded
      }
      A <- B
    }
  }
  diag(A) <- 0L
  A
}

# -------------------------
# 3) Method wrappers returning immediate X×X adjacency (0/1)
# -------------------------
# GES: extract CPDAG block for Xs; keep symmetric entries for undirected edges (both directions 1)
ges_fn_immediate <- function(sim_obj) {
  dat <- sim_obj$dat
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total <- d_z + d_x
  mat <- as.matrix(dat)
  colnames(mat) <- colnames(dat)
  
  # fixedGaps: forbid outgoing from Z and Z->Z
  fixedGaps <- matrix(FALSE, total, total)
  if (d_z > 0) {
    fixedGaps[1:d_z, ] <- TRUE
    fixedGaps[1:d_z, 1:d_z] <- TRUE
  }
  maxDeg <- max(1, round(sim_obj$params$sp * total))
  score <- new("GaussL0penObsScore", mat)
  fit <- tryCatch(ges(score, labels = colnames(mat), maxDegree = maxDeg,
                      fixedGaps = fixedGaps, iterate = TRUE, verbose = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(clean_adj_immediate(NULL, d_x))
  
  ess <- as(fit$essgraph, "matrix")
  diag(ess) <- 0
  # extract X block
  if (d_z > 0) {
    x_idx <- (d_z + 1):ncol(ess)
    A_x <- ess[x_idx, x_idx, drop = FALSE]
  } else A_x <- ess
  # ensure 0/1 ints
  A_x <- clean_adj_immediate(A_x, d_x)
  A_x
}

# LiNGAM (pcalg::lingam) -> B matrix threshold -> immediate
lingam_fn_immediate <- function(sim_obj, thr = 1e-6) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  mat <- as.matrix(sim_obj$dat)
  fit <- tryCatch(pcalg::lingam(mat), error = function(e) NULL)
  if (is.null(fit)) return(clean_adj_immediate(NULL, d_x))
  B <- fit$Bpruned
  if (is.null(B)) B <- fit$B
  if (is.null(B)) return(clean_adj_immediate(NULL, d_x))
  if (nrow(B) == (d_z + d_x)) {
    Bx <- B[(d_z+1):(d_z+d_x), (d_z+1):(d_z+d_x), drop = FALSE]
  } else Bx <- B
  A <- (abs(Bx) > thr) * 1L
  clean_adj_immediate(A, d_x)
}

# CBL and ASCEND: user-defined functions; wrap defensively
wrapper_cbl_immediate <- function(sim_obj) {
  A <- tryCatch(cbl_fn(sim_obj, B = 20, maxiter = 10, gamma = 0.5),
                error = function(e) NULL)
  #clean_adj_immediate(A, sim_obj$params$d_x)
}
wrapper_ascend_immediate <- function(sim_obj) {
  A <- tryCatch(ascend_fn(sim_obj), error = function(e) NULL)
  clean_adj_immediate(A, sim_obj$params$d_x)
}

# -------------------------
# 4) Immediate-edge evaluation metrics
# -------------------------
compare_immediate <- function(Gtrue, Gest) {
  # Both are immediate adjacency matrices (0/1)
  Gtrue <- clean_adj_immediate(Gtrue)
  Gest  <- clean_adj_immediate(Gest, nrow(Gtrue))
  
  # skeletons: undirected presence via OR
  S_true <- ((Gtrue + t(Gtrue)) > 0) * 1L
  S_est  <- ((Gest  + t(Gest))  > 0) * 1L
  
  # upper-triangle counts for skeleton (undirected comparisons)
  ut <- upper.tri(S_true)
  tp_skel <- sum(S_true[ut] & S_est[ut])
  fp_skel <- sum((!S_true[ut]) & S_est[ut])
  fn_skel <- sum(S_true[ut] & (!S_est[ut]))
  tn_skel <- sum((!S_true[ut]) & (!S_est[ut]))
  
  precision_skel <- if ((tp_skel + fp_skel) == 0) NA_real_ else tp_skel / (tp_skel + fp_skel)
  recall_skel    <- if ((tp_skel + fn_skel) == 0) NA_real_ else tp_skel / (tp_skel + fn_skel)
  f1_skel        <- if (is.na(precision_skel) || is.na(recall_skel) || (precision_skel+recall_skel)==0) NA_real_ else 2 * precision_skel * recall_skel / (precision_skel + recall_skel)
  
  # orientation accuracy: only where skeleton is correct (undirected or directed)
  idx_pairs <- which(S_true == 1 & S_est == 1, arr.ind = TRUE)
  # consider unordered pairs only once
  # keep only upper-tri to avoid double counting
  idx_pairs <- idx_pairs[idx_pairs[,1] < idx_pairs[,2], , drop = FALSE]
  orient_total <- nrow(idx_pairs)
  orient_correct <- 0
  for (k in seq_len(orient_total)) {
    i <- idx_pairs[k,1]; j <- idx_pairs[k,2]
    # true orientation possibilities
    true_ij <- Gtrue[i,j] == 1 && Gtrue[j,i] == 0
    true_ji <- Gtrue[j,i] == 1 && Gtrue[i,j] == 0
    # est orientation
    est_ij <- Gest[i,j] == 1 && Gest[j,i] == 0
    est_ji <- Gest[j,i] == 1 && Gest[i,j] == 0
    if ( (true_ij && est_ij) || (true_ji && est_ji) ) orient_correct <- orient_correct + 1
    # note: if either true or est are undirected (both directions 1) then this counts as not oriented
  }
  orient_accuracy_cond_skel <- if (orient_total == 0) NA_real_ else orient_correct / orient_total
  
  # SHD (immediate): # undirected edge disagreements + orientation mismatches on shared skeleton
  shd_skel <- fp_skel + fn_skel
  orient_mismatch <- orient_total - orient_correct
  shd <- shd_skel + orient_mismatch
  
  list(
    tp_skel = tp_skel, fp_skel = fp_skel, fn_skel = fn_skel, tn_skel = tn_skel,
    precision_skel = precision_skel, recall_skel = recall_skel, f1_skel = f1_skel,
    orient_total = orient_total, orient_correct = orient_correct,
    orient_accuracy_cond_skel = orient_accuracy_cond_skel,
    shd = shd,
    true_directed_count = sum(Gtrue),
    ordered_TP = sum((Gest == 1) & (Gtrue == 1)),
    ordered_FP = sum((Gest == 1) & (Gtrue == 0)),
    ordered_FN = sum((Gest == 0) & (Gtrue == 1))
  )
}

# -------------------------
# 5) Experiment runner (immediate-edge evaluation)
# -------------------------
run_sample_size_experiment_immediate <- function(
    n_vec = c(100,200,500,1000),
    methods = c("ges","lingam","cbl","ASCEND"),
    n_rep = 10,
    d_z = 50, d_x = 6,
    rho = 0.5, r2 = 0.1, lin_pr = 1, sp = 0.8,
    seed_base = 1, cores = max(1, parallel::detectCores()-1),
    verbose = TRUE
) {
  registerDoMC(cores)
  results <- list()
  
  for (ni in seq_along(n_vec)) {
    n_val <- n_vec[ni]
    if (verbose) message("[*] sample size = ", n_val)
    
    runs <- foreach(rep_i = 1:n_rep, .combine = rbind, .packages = c("pcalg","igraph","matrixStats")) %dopar% {
      set.seed(seed_base + 10000*ni + rep_i)
      
      sim <- sim_dat(n = n_val, d_z = d_z, d_x = d_x, rho = rho, r2 = r2, lin_pr = lin_pr, sp = sp)
      
      # ground truth immediate adjacency (X->X)
      Gtrue <- clean_adj_immediate(sim$adj_mat, d = d_x)
      
      # compute estimators (immediate)
      ests <- list(
        ges    = ges_fn_immediate(sim),
        lingam = lingam_fn_immediate(sim),
        cbl    = wrapper_cbl_immediate(sim),
        ASCEND = wrapper_ascend_immediate(sim)
      )
      
      out_rows <- list()
      for (m in methods) {
        Gest <- ests[[m]]
        sc <- compare_immediate(Gtrue, Gest)
        out_rows[[m]] <- data.frame(
          method = m, n = n_val, rep = rep_i, seed = seed_base + 10000*ni + rep_i,
          tp_skel = sc$tp_skel, fp_skel = sc$fp_skel, fn_skel = sc$fn_skel, tn_skel = sc$tn_skel,
          precision_skel = sc$precision_skel, recall_skel = sc$recall_skel, f1_skel = sc$f1_skel,
          orient_total = sc$orient_total, orient_correct = sc$orient_correct,
          orient_accuracy_cond_skel = sc$orient_accuracy_cond_skel,
          true_directed_count = sc$true_directed_count,
          ordered_TP = sc$ordered_TP, ordered_FP = sc$ordered_FP, ordered_FN = sc$ordered_FN,
          precision_dir = if ((sc$ordered_TP + sc$ordered_FP)==0) NA_real_ else sc$ordered_TP/(sc$ordered_TP+sc$ordered_FP),
          recall_dir = if (sc$true_directed_count==0) NA_real_ else sc$ordered_TP / sc$true_directed_count,
          f1_dir = if ((sc$ordered_TP + sc$ordered_FP + sc$true_directed_count)==0) NA_real_ else 2*sc$ordered_TP/(sc$ordered_TP + sc$ordered_FP + sc$true_directed_count),
          shd = sc$shd, true_edges = sc$true_directed_count,
          stringsAsFactors = FALSE
        )
      }
      do.call(rbind, out_rows)
    } # foreach rep
    
    results[[as.character(n_val)]] <- runs
  }
  
  raw_dt <- as.data.table(do.call(rbind, results))
  agg_dt <- raw_dt[, .(
    n_rep = .N,
    mean_f1_skel = mean(f1_skel, na.rm = TRUE),
    sd_f1_skel = sd(f1_skel, na.rm = TRUE),
    mean_orient_acc = mean(orient_accuracy_cond_skel, na.rm = TRUE),
    sd_orient_acc = sd(orient_accuracy_cond_skel, na.rm = TRUE),
    mean_shd = mean(shd, na.rm = TRUE),
    sd_shd = sd(shd, na.rm = TRUE),
    mean_true_edges = mean(true_edges, na.rm = TRUE)
  ), by = .(method, n)]
  
  list(raw = raw_dt, agg = agg_dt)
}

# -------------------------
# 6) Plotting (same custom theme you like)
# -------------------------
seaborn_theme <- function(base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel/Plot Background: A very light grey, like the Seaborn default
      panel.background = element_rect(fill = "#EEEEEE", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      
      # Grid Lines: White/light grey and solid, the most distinctive feature
      panel.grid.major = element_line(colour = "white", linewidth = 0.5),
      panel.grid.minor = element_line(colour = "white", linewidth = 0.25),
      
      # Axis Lines: Keep axis lines visible
      #axis.line = element_line(colour = "black"),
      axis.line = element_blank(),
      
      # Panel Border: No explicit border is usually drawn in Seaborn
      panel.border = element_blank(),
      
      # Additional theme elements from your original base_theme for consistency
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.8, "lines"),
      
      # Reduced margins for less empty space
      plot.margin = margin(0.2, 0.5, 0.2, 0.5, "cm")
    )
}

plot_sample_size_results <- function(agg_dt) {
  method_colors <- c("ASCEND" = "#4DAF4A", "ges" = "#984EA3",
                     "lingam" = "#377EB8", "cbl" = "#E41A1C")
  agg_dt$method <- factor(agg_dt$method, levels = c("ASCEND", "ges", "lingam", "cbl"))
  custom_theme <- seaborn_theme()
  
  p_f1 <- ggplot(agg_dt, aes(x = n, y = mean_f1_skel, color = method, group = method)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
    geom_ribbon(aes(ymin = pmax(0, mean_f1_skel - sd_f1_skel/sqrt(n_rep)),
                    ymax = mean_f1_skel + sd_f1_skel/sqrt(n_rep), fill = method),
                alpha = 0.2, color = NA) +
    scale_x_log10(breaks = unique(agg_dt$n)) +
    scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
    labs(x = "Sample Size", y = "Skeleton F1 (immediate)") +
    custom_theme
  
  p_orient <- ggplot(agg_dt, aes(x = n, y = mean_orient_acc, color = method, group = method)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
    geom_ribbon(aes(ymin = pmax(0, mean_orient_acc - sd_orient_acc/sqrt(n_rep)),
                    ymax = mean_orient_acc + sd_orient_acc/sqrt(n_rep), fill = method),
                alpha = 0.2, color = NA) +
    scale_x_log10(breaks = unique(agg_dt$n)) +
    scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
    labs(x = "Sample Size", y = "Orientation accuracy (conditional on skeleton)") +
    custom_theme
  
  p_shd <- ggplot(agg_dt, aes(x = n, y = mean_shd, color = method, group = method)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.5) +
    geom_ribbon(aes(ymin = pmax(0, mean_shd - sd_shd/sqrt(n_rep)),
                    ymax = mean_shd + sd_shd/sqrt(n_rep), fill = method),
                alpha = 0.2, color = NA) +
    scale_x_log10(breaks = unique(agg_dt$n)) +
    scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
    labs(x = "Sample Size", y = "SHD (immediate)") +
    custom_theme
  
  list(f1 = p_f1, orient = p_orient, shd = p_shd)
}

# -------------------------
# 7) Example run (smaller N_rep for quick test)
# -------------------------
# Adjust n_vec and n_rep for a full experiment; this is an example
# Make sure you have cbl_fn() and ascend_fn() in your environment.
n_vec <- c(500, 1000, 2000)
res_n <- run_sample_size_experiment_immediate(
  n_vec = n_vec,
  methods = c("ges", "lingam", "cbl", "ASCEND"),
  n_rep = 8,
  d_z = 60, 
  d_x = 10,
  rho = 0.5, 
  r2 = 0.5, 
  lin_pr = 1, 
  sp = 0.5,
  seed_base = 123, cores = max(1, parallel::detectCores()-1),
  verbose = TRUE
)

print(res_n$agg)
plots <- plot_sample_size_results(res_n$agg)
print(plots$f1)
print(plots$orient)
print(plots$shd)
