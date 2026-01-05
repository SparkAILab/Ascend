# ============================================================================
# ASCEND: COMPREHENSIVE COMPARISON WITH ESTABLISHED NETWORK INFERENCE METHODS
# ============================================================================

# Load required packages
library(ppcor)      # Partial correlation (Kim 2015)
library(pcalg)      # PC algorithm (Spirtes et al. 2000)
library(glasso)     # Graphical Lasso (Friedman et al. 2008)
library(GENIE3)     # Tree-based GRN inference (Huynh-Thu et al. 2010)
library(dplyr)      # Data manipulation
library(ggplot2)    # Visualization
library(foreach)    # Parallel processing

# Set up parallel processing (optional)
registerDoMC(cores = max(1, parallel::detectCores() - 1))

# ============================================================================
# 1. DATA PREPROCESSING
# ============================================================================

robust_scale <- function(x) {
  if (all(is.na(x))) return(x)
  
  # Winsorize extreme outliers (5% tails)
  q05 <- quantile(x, 0.05, na.rm = TRUE)
  q95 <- quantile(x, 0.95, na.rm = TRUE)
  x[x < q05] <- q05
  x[x > q95] <- q95
  
  # Standardize using median and MAD (more robust than mean/SD)
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  
  if (mad_val == 0 || is.na(mad_val)) {
    # Fallback to standard scaling if MAD fails
    x <- scale(x)
    return(x[, 1])
  }
  
  return((x - med) / mad_val)
}

preprocess_data <- function(sim_obj, method_type = "standard") {
  # Different preprocessing for different method families
  dat <- as.data.frame(sim_obj$dat)
  
  if (method_type == "genie3") {
    # GENIE3 prefers data in original scale or gently normalized
    for (col in colnames(dat)) {
      dat[[col]] <- scale(dat[[col]])  # Simple z-scoring
    }
  } else {
    # Robust scaling for other methods
    for (col in colnames(dat)) {
      dat[[col]] <- robust_scale(dat[[col]])
    }
  }
  
  # Handle any NA/Inf values
  dat[is.na(dat) | is.infinite(as.matrix(dat))] <- 0
  
  # Ensure no zero-variance columns
  col_vars <- apply(dat, 2, var, na.rm = TRUE)
  zero_var <- which(col_vars == 0)
  if (length(zero_var) > 0) {
    n <- nrow(dat)
    for (col in zero_var) {
      dat[, col] <- dat[, col] + rnorm(n, 0, 1e-8)
    }
  }
  
  sim_obj$dat <- dat
  return(sim_obj)
}

# ============================================================================
# 2. EVALUATION METRICS (SKELETON COMPARISON)
# ============================================================================

extract_x_submatrix <- function(adj_matrix, x_cols) {
  # Extract the X->X submatrix from a full adjacency matrix
  if (is.null(adj_matrix) || nrow(adj_matrix) == 0) {
    return(NULL)
  }
  
  # Get indices of X variables
  x_indices <- which(colnames(adj_matrix) %in% x_cols)
  
  if (length(x_indices) == 0) {
    return(matrix(0, nrow = length(x_cols), ncol = length(x_cols)))
  }
  
  # Extract submatrix
  x_submatrix <- adj_matrix[x_indices, x_indices, drop = FALSE]
  
  # Ensure correct dimensions
  if (nrow(x_submatrix) != length(x_cols)) {
    x_submatrix <- matrix(0, nrow = length(x_cols), ncol = length(x_cols))
  }
  
  return(x_submatrix)
}

compare_skeletons <- function(Gtrue, Gest) {
  # Ensure both are matrices
  Gtrue <- as.matrix(Gtrue)
  Gest <- as.matrix(Gest)
  
  if (nrow(Gtrue) != nrow(Gest) || ncol(Gtrue) != ncol(Gest)) {
    stop(paste("Dimension mismatch: Gtrue =", dim(Gtrue), "Gest =", dim(Gest)))
  }
  
  # Convert to undirected skeletons
  S_true <- ((Gtrue + t(Gtrue)) > 0) * 1L
  S_est <- ((Gest + t(Gest)) > 0) * 1L
  
  # Compare only upper triangle (undirected edges)
  ut <- upper.tri(S_true)
  tp <- sum(S_true[ut] & S_est[ut])
  fp <- sum((!S_true[ut]) & S_est[ut])
  fn <- sum(S_true[ut] & (!S_est[ut]))
  
  # Calculate metrics (handle edge cases)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  
  f1 <- if (is.na(precision) || is.na(recall) || 
            (precision + recall) == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }
  
  return(list(
    precision = precision,
    recall = recall,
    f1 = f1,
    tp = tp,
    fp = fp,
    fn = fn,
    shd = fp + fn  # Structural Hamming Distance
  ))
}

# ============================================================================
# 3. METHOD IMPLEMENTATIONS (UPDATED TO USE Z+X)
# ============================================================================

# 3.1 Correlation (baseline) - Pearson correlation thresholding
run_correlation <- function(sim_obj, top_pct = 0.3) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Use all data (Z+X)
  full_data <- as.matrix(sim_obj$dat)
  
  # Compute correlation matrix
  cor_mat <- abs(cor(full_data))
  
  # Initialize full adjacency matrix
  adj_full <- matrix(0, total_vars, total_vars)
  rownames(adj_full) <- colnames(adj_full) <- colnames(full_data)
  
  # For evaluation, we only care about X->X edges
  # But we still need to compute correlations on full data
  x_indices <- (d_z + 1):total_vars
  x_cor_mat <- cor_mat[x_indices, x_indices, drop = FALSE]
  
  # Threshold only the X submatrix
  score_vec <- x_cor_mat[upper.tri(x_cor_mat)]
  n_to_keep <- round(length(score_vec) * top_pct)
  thresh <- sort(score_vec, decreasing = TRUE)[n_to_keep]
  
  # Set X->X edges
  x_adj <- (x_cor_mat >= thresh) * 1
  diag(x_adj) <- 0
  
  # Insert into full adjacency
  adj_full[x_indices, x_indices] <- x_adj
  
  return(adj_full)
}

# 3.2 ARACNE - Updated to use Z+X
run_aracne <- function(sim_obj, top_pct = 0.3) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Use all data (Z+X)
  full_data <- as.matrix(sim_obj$dat)
  full_data[is.na(full_data)] <- 0
  full_data <- full_data + matrix(rnorm(prod(dim(full_data)), 0, 1e-10), 
                                  nrow = nrow(full_data))
  
  tryCatch({
    # Compute mutual information on full data
    mim <- minet::build.mim(full_data, estimator = "spearman")
    
    # Apply ARACNE to full matrix
    net_full <- minet::aracne(mim, eps = 0)
    
    # Extract X submatrix
    x_indices <- (d_z + 1):total_vars
    x_net <- net_full[x_indices, x_indices, drop = FALSE]
    
    # Threshold X submatrix
    score_vec <- x_net[upper.tri(x_net)]
    survivors <- score_vec[score_vec > 0]
    
    if (length(survivors) < 2) {
      x_adj <- matrix(0, d_x, d_x)
    } else {
      n_to_keep <- max(1, round(length(score_vec) * top_pct))
      thresh <- quantile(score_vec, 1 - top_pct, na.rm = TRUE)
      x_adj <- (x_net >= thresh & x_net > 0) * 1
    }
    
    diag(x_adj) <- 0
    
    # Create full adjacency
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(full_data)
    adj_full[x_indices, x_indices] <- x_adj
    
    return(adj_full)
  }, error = function(e) {
    message("ARACNE error: ", e$message)
    # Return empty full adjacency
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(sim_obj$dat)
    return(adj_full)
  })
}

# 3.3 CLR - Updated to use Z+X
run_clr <- function(sim_obj, top_pct = 0.3) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Use all data (Z+X)
  full_data <- as.matrix(sim_obj$dat)
  full_data[is.na(full_data)] <- 0
  full_data <- full_data + matrix(rnorm(prod(dim(full_data)), 0, 1e-10), 
                                  nrow = nrow(full_data))
  
  tryCatch({
    # Compute mutual information on full data
    mim <- minet::build.mim(full_data, estimator = "pearson")
    net_full <- minet::clr(mim)
    
    # Extract X submatrix
    x_indices <- (d_z + 1):total_vars
    x_net <- net_full[x_indices, x_indices, drop = FALSE]
    
    # Threshold X submatrix
    score_vec <- x_net[upper.tri(x_net)]
    n_to_keep <- max(1, round(length(score_vec) * top_pct))
    thresh <- quantile(score_vec, 1 - top_pct, na.rm = TRUE)
    
    x_adj <- (x_net >= thresh) * 1
    diag(x_adj) <- 0
    
    # Create full adjacency
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(full_data)
    adj_full[x_indices, x_indices] <- x_adj
    
    return(adj_full)
  }, error = function(e) {
    message("CLR error: ", e$message)
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(sim_obj$dat)
    return(adj_full)
  })
}

# 3.4 PPCOR - Partial correlation with FDR correction (uses Z+X)
run_ppcor <- function(sim_obj, alpha = 0.05) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Use all data (Z+X)
  full_data <- as.matrix(sim_obj$dat)
  
  if (ncol(full_data) < 3) {
    return(run_correlation(sim_obj))
  }
  
  tryCatch({
    # Compute partial correlations on full data
    pcor_result <- ppcor::pcor(full_data, method = "pearson")
    
    # Extract p-values
    pvals <- pcor_result$p.value
    diag(pvals) <- 1  # No self-edges
    
    # Apply Benjamini-Hochberg FDR correction
    ut <- upper.tri(pvals)
    pvals_vec <- pvals[ut]
    pvals_adj <- p.adjust(pvals_vec, method = "BH")
    
    # Create full adjacency matrix
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(full_data)
    adj_full[ut] <- (pvals_adj < alpha) * 1
    adj_full <- adj_full + t(adj_full)  # Symmetrize
    
    # For evaluation, we only care about X->X edges
    # But we return the full matrix for consistency
    
    return(adj_full)
  }, error = function(e) {
    return(run_correlation(sim_obj))
  })
}

# 3.5 GLASSO - Sparse inverse covariance estimation (uses Z+X)
run_glasso <- function(sim_obj, rho = 0.1) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Use all data (Z+X)
  full_data <- as.matrix(sim_obj$dat)
  
  tryCatch({
    # Standardize data
    full_scaled <- scale(full_data)
    
    # Compute sample covariance
    S <- cov(full_scaled)
    
    # Graphical Lasso with given regularization
    result <- glasso::glasso(S, rho = rho)
    
    # Non-zero entries in precision matrix = edges
    adj_full <- (result$wi != 0) * 1
    diag(adj_full) <- 0
    adj_full <- pmax(adj_full, t(adj_full))  # Symmetrize
    
    rownames(adj_full) <- colnames(adj_full) <- colnames(full_data)
    
    return(adj_full)
  }, error = function(e) {
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(sim_obj$dat)
    return(adj_full)
  })
}

# 3.6 PC Algorithm - Constraint-based causal discovery (uses Z+X)
run_pc <- function(sim_obj, alpha = 0.05) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Use all data (Z+X)
  full_data <- as.matrix(sim_obj$dat)
  
  if (ncol(full_data) < 3) {
    return(run_correlation(sim_obj))
  }
  
  tryCatch({
    # Prepare sufficient statistics for PC
    full_scaled <- scale(full_data)
    suffStat <- list(C = cor(full_scaled), n = nrow(full_scaled))
    
    # Run PC algorithm on all variables
    pc_fit <- pcalg::pc(
      suffStat = suffStat,
      indepTest = pcalg::gaussCItest,
      p = ncol(full_scaled),
      alpha = alpha,
      skel.method = "stable",
      u2pd = "retry"
    )
    
    # Extract skeleton (undirected edges)
    pc_amat <- as(pc_fit, "amat")
    adj_full <- matrix(as.numeric(pc_amat), nrow = nrow(pc_amat))
    rownames(adj_full) <- colnames(adj_full) <- colnames(full_data)
    
    # Convert to undirected skeleton
    adj_full <- (adj_full | t(adj_full)) * 1
    
    return(adj_full)
  }, error = function(e) {
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(sim_obj$dat)
    return(adj_full)
  })
}

# 3.7 GENIE3 - Tree-based ensemble (uses Z+X)
run_genie3 <- function(sim_obj, n_trees = 500, top_pct = 0.3) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Use all data (Z+X)
  full_data <- as.matrix(sim_obj$dat)
  
  tryCatch({
    full_transposed <- t(full_data)
    var_names <- colnames(full_data)
    rownames(full_transposed) <- var_names
    
    # Run GENIE3 on all variables
    weight_matrix <- GENIE3::GENIE3(
      exprMatrix = full_transposed,
      treeMethod = "RF",
      nTrees = n_trees,
      nCores = 1,
      verbose = FALSE
    )
    
    # Convert to symmetric scores
    sym_weights <- (weight_matrix + t(weight_matrix)) / 2
    
    # Extract X submatrix
    x_indices <- (d_z + 1):total_vars
    x_weights <- sym_weights[x_indices, x_indices, drop = FALSE]
    
    # Threshold X submatrix
    thresh <- quantile(x_weights[upper.tri(x_weights)], 1 - top_pct, na.rm = TRUE)
    x_adj <- (x_weights >= thresh) * 1
    diag(x_adj) <- 0
    
    # Create full adjacency
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- var_names
    adj_full[x_indices, x_indices] <- x_adj
    
    return(adj_full)
  }, error = function(e) {
    message("GENIE3 error: ", e$message)
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- colnames(sim_obj$dat)
    return(adj_full)
  })
}

# 3.8 ASCEND - Your method (already uses Z+X by design)
run_ascend <- function(sim_obj, alpha = 0.05) {
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  total_vars <- d_z + d_x
  
  # Check if ascend_fn exists
  if (!exists("ascend_fn")) {
    stop("ASCEND implementation (ascend_fn) not found. Please load your ASCEND code.")
  }
  
  tryCatch({
    # Run ASCEND
    adj_directed <- ascend_fn(sim_obj, alpha = alpha)
    
    # ASCEND returns X-to-X adjacency directly
    # But we need to embed it in full adjacency for consistency
    
    # Create full adjacency matrix
    adj_full <- matrix(0, total_vars, total_vars)
    all_vars <- c(paste0("z", 1:d_z), paste0("x", 1:d_x))
    rownames(adj_full) <- colnames(adj_full) <- all_vars
    
    # Insert ASCEND's X adjacency
    x_indices <- (d_z + 1):total_vars
    if (nrow(adj_directed) == d_x) {
      # Convert to undirected skeleton for comparison
      adj_skeleton <- (adj_directed | t(adj_directed)) * 1
      adj_full[x_indices, x_indices] <- adj_skeleton
    }
    
    return(adj_full)
  }, error = function(e) {
    message("ASCEND error: ", e$message)
    adj_full <- matrix(0, total_vars, total_vars)
    rownames(adj_full) <- colnames(adj_full) <- c(paste0("z", 1:d_z), paste0("x", 1:d_x))
    return(adj_full)
  })
}

# ============================================================================
# 4. UPDATED COMPARISON FRAMEWORK
# ============================================================================

compare_all_methods <- function(sim_obj) {
  # Define methods to compare
  methods <- c("ASCEND", "PC", "ARACNE", "CLR", "GLASSO", "PPCOR", "GENIE3", "Correlation")
  
  # Get ground truth (X-to-X subgraph)
  Gtrue <- sim_obj$adj_mat
  d_z <- sim_obj$params$d_z
  d_x <- sim_obj$params$d_x
  
  # Get X variable names for extraction
  x_cols <- paste0("x", 1:d_x)
  
  results <- list()
  
  cat("\n" , strrep("=", 70))
  cat("\nMETHOD COMPARISON: ", d_z, " Z + ", d_x, " X variables\n", sep = "")
  cat(strrep("=", 70), "\n")
  
  # Preprocess data appropriately for each method
  sim_standard <- preprocess_data(sim_obj, "standard")
  sim_genie3 <- preprocess_data(sim_obj, "genie3")
  
  for (method in methods) {
    cat("\nRunning ", method, "... ", sep = "")
    
    # Select appropriate preprocessed data
    if (method == "GENIE3") {
      sim_used <- sim_genie3
    } else {
      sim_used <- sim_standard
    }
    
    # Run method (now on Z+X data)
    start_time <- Sys.time()
    
    adj_full <- switch(method,
                       "ASCEND" = run_ascend(sim_used),
                       "PC" = run_pc(sim_used),
                       "ARACNE" = run_aracne(sim_used),
                       "CLR" = run_clr(sim_used),
                       "GLASSO" = run_glasso(sim_used),
                       "PPCOR" = run_ppcor(sim_used),
                       "GENIE3" = run_genie3(sim_used),
                       "Correlation" = run_correlation(sim_used)
    )
    
    runtime <- difftime(Sys.time(), start_time, units = "secs")
    
    # Extract X-to-X submatrix for evaluation
    adj_x_submatrix <- extract_x_submatrix(adj_full, x_cols)
    
    # Ensure dimensions match
    if (!is.null(adj_x_submatrix) && 
        nrow(adj_x_submatrix) == d_x && 
        ncol(adj_x_submatrix) == d_x) {
      
      # Evaluate only on X-to-X subgraph
      metrics <- compare_skeletons(Gtrue, adj_x_submatrix)
      
      results[[method]] <- list(
        full_adjacency = adj_full,
        x_adjacency = adj_x_submatrix,
        metrics = metrics,
        runtime = as.numeric(runtime),
        edges_found = sum(adj_x_submatrix[upper.tri(adj_x_submatrix)])
      )
      
      cat(sprintf("F1=%.3f (%.1fs)", metrics$f1, runtime))
    } else {
      results[[method]] <- list(
        full_adjacency = adj_full,
        x_adjacency = NULL,
        metrics = list(f1 = NA, precision = NA, recall = NA),
        runtime = as.numeric(runtime),
        edges_found = 0
      )
      cat("Failed")
    }
  }
  
  cat("\n", strrep("=", 70), "\n")
  
  return(list(
    ground_truth = Gtrue,
    method_results = results,
    simulation_params = sim_obj$params
  ))
}

# ============================================================================
# 5. UPDATED BATCH EXPERIMENT
# ============================================================================

run_comprehensive_experiment <- function(n_reps = 20,
                                         n = 1000, d_z = 100, d_x = 25,
                                         sparsity = 0.3, signal = 0.3) {
  
  methods <- c("ASCEND", "PC", "ARACNE", "CLR", "GLASSO", "PPCOR", "GENIE3", "Correlation")
  
  cat("\n" , strrep("=", 70))
  cat("\nCOMPREHENSIVE EVALUATION: ", n_reps, " REPLICATIONS\n", sep = "")
  cat("Parameters: n=", n, ", d_z=", d_z, ", d_x=", d_x, 
      ", sparsity=", sparsity, ", signal=", signal, "\n", sep = "")
  cat(strrep("=", 70), "\n")
  
  # Initialize results storage
  all_results <- data.frame()
  
  for (rep in 1:n_reps) {
    cat("\nReplication ", rep, "/", n_reps, ": ", sep = "")
    
    # Generate simulation
    sim_obj <- sim_dat(
      n = n, d_z = d_z, d_x = d_x,
      sp = sparsity, r2 = signal,
      lin_pr = 1.0,
      seed = 12345 + rep * 100
    )
    
    true_edges <- sum(sim_obj$adj_mat[upper.tri(sim_obj$adj_mat)])
    cat(true_edges, " true X-to-X edges\n", sep = "")
    
    # Compare methods
    comparison <- compare_all_methods(sim_obj)
    
    # Collect results
    for (method in methods) {
      if (!is.null(comparison$method_results[[method]]$metrics)) {
        m <- comparison$method_results[[method]]$metrics
        all_results <- rbind(all_results, data.frame(
          Rep = rep,
          Method = method,
          F1 = m$f1,
          Precision = m$precision,
          Recall = m$recall,
          Runtime = comparison$method_results[[method]]$runtime,
          Edges_Found = comparison$method_results[[method]]$edges_found,
          True_Edges = true_edges,
          Method_Class = ifelse(method == "ASCEND", "Tiered Causal",
                                ifelse(method %in% c("PC"), "Standard Causal",
                                       ifelse(method %in% c("GLASSO", "PPCOR"), 
                                              "Graphical Model", "Association"))),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # ==========================================================================
  # 6. STATISTICAL ANALYSIS AND REPORTING
  # ==========================================================================
  
  if (nrow(all_results) > 0) {
    # Summary statistics
    summary_stats <- all_results %>%
      group_by(Method, Method_Class) %>%
      summarise(
        N = sum(!is.na(F1)),
        F1_mean = mean(F1, na.rm = TRUE),
        F1_sd = sd(F1, na.rm = TRUE),
        F1_se = F1_sd / sqrt(N),
        Precision_mean = mean(Precision, na.rm = TRUE),
        Recall_mean = mean(Recall, na.rm = TRUE),
        Runtime_mean = mean(Runtime, na.rm = TRUE),
        Runtime_sd = sd(Runtime, na.rm = TRUE),
        Success_Rate = round(N / n_reps * 100, 1),
        .groups = 'drop'
      ) %>%
      arrange(desc(F1_mean))
    
    # Print comprehensive summary
    cat("\n" , strrep("=", 70))
    cat("\nEXPERIMENT SUMMARY\n")
    cat(strrep("=", 70), "\n\n")
    
    cat(sprintf("%-12s %-15s %8s %8s %8s %8s %8s %6s\n", 
                "Method", "Class", "F1", "±SE", "Prec", "Rec", "Time(s)", "Success%"))
    cat(strrep("-", 80), "\n")
    
    for (i in 1:nrow(summary_stats)) {
      cat(sprintf("%-12s %-15s %8.3f %8.3f %8.3f %8.3f %8.1f %8.0f%%\n",
                  summary_stats$Method[i],
                  summary_stats$Method_Class[i],
                  summary_stats$F1_mean[i],
                  summary_stats$F1_se[i],
                  summary_stats$Precision_mean[i],
                  summary_stats$Recall_mean[i],
                  summary_stats$Runtime_mean[i],
                  summary_stats$Success_Rate[i]))
    }
    
    # Statistical significance testing
    cat("\n" , strrep("=", 70))
    cat("\nSTATISTICAL SIGNIFICANCE TESTS (Paired t-tests)\n")
    cat(strrep("=", 70), "\n")
    
    # Compare ASCEND against each baseline
    ascend_f1 <- all_results$F1[all_results$Method == "ASCEND"]
    
    for (baseline in setdiff(methods, "ASCEND")) {
      baseline_f1 <- all_results$F1[all_results$Method == baseline]
      
      # Paired t-test (same replications)
      if (length(ascend_f1) == length(baseline_f1) && 
          length(ascend_f1) > 1) {
        test <- t.test(ascend_f1, baseline_f1, paired = TRUE)
        improvement <- (mean(ascend_f1, na.rm = TRUE) - 
                          mean(baseline_f1, na.rm = TRUE)) / 
          mean(baseline_f1, na.rm = TRUE) * 100
        
        stars <- ifelse(test$p.value < 0.001, "***",
                        ifelse(test$p.value < 0.01, "**",
                               ifelse(test$p.value < 0.05, "*", "NS")))
        
        cat(sprintf("ASCEND vs %-10s: ΔF1 = %+6.3f (%+6.1f%%) p = %.4f %s\n",
                    baseline,
                    mean(ascend_f1, na.rm = TRUE) - mean(baseline_f1, na.rm = TRUE),
                    improvement,
                    test$p.value,
                    stars))
      }
    }
    
    # ======================================================================
    # 7. VISUALIZATION
    # ======================================================================
    
    if (require(ggplot2)) {
      # F1 Score comparison
      p1 <- ggplot(all_results, aes(x = Method, y = F1, fill = Method_Class)) +
        geom_boxplot(alpha = 0.8, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        labs(title = "X-to-X Skeleton Recovery Performance (F1 Score)",
             subtitle = paste(n_reps, "replications | Using Z+X data"),
             x = "Method", y = "F1 Score",
             fill = "Method Class") +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        scale_fill_brewer(palette = "Set2")
      
      # Precision-Recall tradeoff
      p2 <- ggplot(summary_stats, 
                   aes(x = Recall_mean, y = Precision_mean, 
                       color = Method, shape = Method_Class)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = Precision_mean - F1_sd/2, 
                          ymax = Precision_mean + F1_sd/2), 
                      width = 0.01) +
        geom_errorbarh(aes(xmin = Recall_mean - F1_sd/2, 
                           xmax = Recall_mean + F1_sd/2), 
                       height = 0.01) +
        labs(title = "Precision-Recall Tradeoff (X-to-X edges)",
             x = "Recall", y = "Precision",
             color = "Method", shape = "Method Class") +
        theme_bw(base_size = 11) +
        theme(legend.position = "right") +
        scale_color_brewer(palette = "Set1")
      
      # Runtime comparison
      p3 <- ggplot(all_results, aes(x = Method, y = Runtime, fill = Method_Class)) +
        geom_boxplot(alpha = 0.8) +
        scale_y_log10() +
        labs(title = "Computational Runtime (log scale)",
             subtitle = "All methods use Z+X data",
             x = "Method", y = "Runtime (seconds, log10)",
             fill = "Method Class") +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        scale_fill_brewer(palette = "Set3")
      
      # Display plots
      print(p1)
      print(p2)
      print(p3)
    }
  }
  
  return(list(
    raw_results = all_results,
    summary = summary_stats,
    n_replications = n_reps,
    parameters = list(n = n, d_z = d_z, d_x = d_x, sparsity = sparsity, signal = signal)
  ))
}

# ============================================================================
# 8. MAIN EXECUTION
# ============================================================================

# Run the comprehensive evaluation
cat("\nStarting comprehensive evaluation of ASCEND with Z+X data...\n")

results <- run_comprehensive_experiment(
  n_reps = 5,      # Start with 5, increase to 20-30 for paper
  n = 1000,
  d_z = 50,        
  d_x = 15,       
  sparsity = 0.3,
  signal = 0.3
)

# Final summary
if (!is.null(results$summary)) {
  cat("\n" , strrep("=", 70))
  cat("\nFINAL RANKING BY F1 SCORE (X-to-X skeleton recovery)\n")
  cat(strrep("=", 70), "\n")
  
  for (i in 1:nrow(results$summary)) {
    cat(sprintf("%d. %-12s: F1 = %.3f (±%.3f), Class: %s\n", 
                i, 
                results$summary$Method[i],
                results$summary$F1_mean[i],
                results$summary$F1_se[i],
                results$summary$Method_Class[i]))
  }
}

# Save results for reproducibility
saveRDS(results, file = "ascend_comparison_results_zx.rds")
cat("\nResults saved to 'ascend_comparison_results_zx.rds'\n")