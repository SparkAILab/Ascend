## ASCEND-PC vs plain PC on synthetic two-tier data.
##
## Sweep over n_samples at fixed d_x = 20, d_z = 40. Three methods
## compared across 10 replicates per sample size:
##
##   ASCEND-PC      ASCEND uses Z+X, PC uses X with ASCEND's constraints
##   Plain PC (X)   PC sees only X (standard practice)
##   Plain PC (Z+X) PC sees Z and X as flat node set; evaluated on
##                  X-X induced subgraph
##
## Ground truth: sim_obj$adj_xx (true X->X direct edges, plus the
## transitive closure for the ancestral evaluation).
##
## Outputs:
##   sim_sweep_results.rds       per-rep metric records (full detail)
##   sim_sweep_results.csv       flat table for plotting
##   sim_sweep_checkpoint.rds    auto-saved every rep for crash recovery
##
## Inputs sourced (in working directory):
##   ascend.R       canonical ASCEND
##   simulation.R   sim_dat()

# source("ascend.R")
# source("simulation.R")

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(pcalg)
  library(PRROC)
  library(dplyr)
})

## --- Configuration -------------------------------------------------------

D_X <- 20
D_Z <- 40
N_LEVELS <- c(256, 1024, 4096, 16384)
N_REPS   <- 10

# Simulation parameters held constant. These are sim_dat's own
# defaults and match the regime used by the existing GRN benchmark
# where ASCEND outperforms baselines. Diagnostic confirmed:
# - Direct X-X edge density ~5-8%
# - Ancestral prevalence ~5-9%
# - Ancestor/non-ancestor correlation ratio ~1.3-1.8 (sufficient signal)
SIM_PARAMS <- list(
  r2       = 0.5,
  lin_pr   = 1,
  sp       = 0.9,         # ER density factor (higher = sparser)
  p_cross  = 0.05,        # P(Z -> X) edge
  x_effect = 0.8
)

ALPHA_PC <- 0.05

CHECKPOINT_PATH <- "sim_sweep_checkpoint.rds"
RESULTS_RDS     <- "sim_sweep_results.rds"
RESULTS_CSV     <- "sim_sweep_results.csv"

## --- Utility: get_true_ancestral_matrix (transitive closure) ------------

get_true_ancestral_matrix <- function(adj_xx) {
  d <- nrow(adj_xx); nms <- rownames(adj_xx)
  A <- t(adj_xx); diag(A) <- 0; A[is.na(A)] <- 0
  reach <- (A > 0) * 1; Ap <- A
  for (k in seq_len(d - 1)) {
    Ap <- ((Ap %*% A) > 0) * 1
    reach <- ((reach + Ap) > 0) * 1
  }
  rownames(reach) <- colnames(reach) <- nms
  diag(reach) <- NA
  reach
}

## --- Methods (each returns d_x x d_x matrix on labels x1..xN) -----------

run_ascend_pc <- function(sim_obj, alpha = ALPHA_PC) {
  amat_ascend <- ascend_fn(sim_obj)
  
  d_x <- length(grep("^x", colnames(sim_obj$dat)))
  fg_names <- paste0("x", seq_len(d_x))
  amat_fg <- amat_ascend[fg_names, fg_names]
  n_fg <- length(fg_names)
  
  # Skeleton constraints.
  fixed_gaps <- matrix(FALSE, n_fg, n_fg,
                       dimnames = list(fg_names, fg_names))
  n_gaps <- 0
  for (i in 1:(n_fg - 1)) for (j in (i + 1):n_fg) {
    if (!is.na(amat_fg[i, j]) && !is.na(amat_fg[j, i]) &&
        amat_fg[i, j] == 0 && amat_fg[j, i] == 0) {
      fixed_gaps[i, j] <- fixed_gaps[j, i] <- TRUE
      n_gaps <- n_gaps + 1
    }
  }
  
  # Build correlation matrix on X only.
  X <- as.data.frame(sim_obj$dat[, fg_names, with = FALSE])
  X <- X[complete.cases(X), ]
  cor_mat <- cor(X, method = "spearman")
  min_eig <- min(eigen(cor_mat, only.values = TRUE)$values)
  if (min_eig <= 1e-10) {
    lambda  <- abs(min_eig) + 1e-4
    cor_mat <- cor_mat + lambda * diag(n_fg)
    d_sc    <- sqrt(diag(cor_mat))
    cor_mat <- cor_mat / outer(d_sc, d_sc)
  }
  suffStat <- list(C = cor_mat, n = nrow(X))
  
  pc_fit <- pc(suffStat = suffStat, indepTest = gaussCItest,
               labels = fg_names, alpha = alpha,
               fixedGaps = fixed_gaps, verbose = FALSE,
               maj.rule = TRUE, solve.confl = TRUE)
  amat_pc <- as(pc_fit, "amat")
  rownames(amat_pc) <- colnames(amat_pc) <- fg_names
  
  # Apply ASCEND's directional constraints post hoc.
  amat_constrained <- amat_pc
  for (i in 1:n_fg) for (j in 1:n_fg) {
    if (i == j) next
    if (amat_constrained[i, j] == 1) {
      if (!is.na(amat_fg[j, i]) && amat_fg[j, i] == 1 &&
          !is.na(amat_fg[i, j]) && amat_fg[i, j] == 0) {
        if (amat_constrained[j, i] == 1) {
          amat_constrained[i, j] <- 0
        } else {
          amat_constrained[i, j] <- 0
          amat_constrained[j, i] <- 0
        }
      }
    }
  }
  for (i in 1:(n_fg - 1)) for (j in (i + 1):n_fg) {
    if (amat_constrained[i, j] == 1 && amat_constrained[j, i] == 1) {
      if (!is.na(amat_fg[i, j]) && amat_fg[i, j] == 1 &&
          !is.na(amat_fg[j, i]) && amat_fg[j, i] == 0) {
        amat_constrained[i, j] <- 1; amat_constrained[j, i] <- 0
      } else if (!is.na(amat_fg[j, i]) && amat_fg[j, i] == 1 &&
                 !is.na(amat_fg[i, j]) && amat_fg[i, j] == 0) {
        amat_constrained[j, i] <- 1; amat_constrained[i, j] <- 0
      }
    }
  }
  
  list(amat = amat_constrained,
       amat_fg_ascend = amat_fg,
       n_gaps = n_gaps,
       n_pairs_total = n_fg * (n_fg - 1) / 2,
       prune_pct = 100 * n_gaps / (n_fg * (n_fg - 1) / 2))
}

run_plain_pc_x <- function(sim_obj, alpha = ALPHA_PC) {
  d_x <- length(grep("^x", colnames(sim_obj$dat)))
  fg_names <- paste0("x", seq_len(d_x))
  
  X <- as.data.frame(sim_obj$dat[, fg_names, with = FALSE])
  X <- X[complete.cases(X), ]
  n_fg <- ncol(X)
  cor_mat <- cor(X, method = "spearman")
  min_eig <- min(eigen(cor_mat, only.values = TRUE)$values)
  if (min_eig <= 1e-10) {
    lambda  <- abs(min_eig) + 1e-4
    cor_mat <- cor_mat + lambda * diag(n_fg)
    d_sc    <- sqrt(diag(cor_mat))
    cor_mat <- cor_mat / outer(d_sc, d_sc)
  }
  suffStat <- list(C = cor_mat, n = nrow(X))
  
  pc_fit <- pc(suffStat = suffStat, indepTest = gaussCItest,
               labels = fg_names, alpha = alpha,
               verbose = FALSE, maj.rule = TRUE, solve.confl = TRUE)
  amat_p <- as(pc_fit, "amat")
  rownames(amat_p) <- colnames(amat_p) <- fg_names
  amat_p
}

run_plain_pc_zx <- function(sim_obj, alpha = ALPHA_PC) {
  d_x <- length(grep("^x", colnames(sim_obj$dat)))
  fg_names <- paste0("x", seq_len(d_x))
  
  all_data <- as.data.frame(sim_obj$dat)
  all_data <- all_data[complete.cases(all_data), ]
  n_all <- ncol(all_data)
  cor_full <- cor(all_data, method = "spearman")
  min_eig <- min(eigen(cor_full, only.values = TRUE)$values)
  if (min_eig <= 1e-10) {
    lambda  <- abs(min_eig) + 1e-4
    cor_full <- cor_full + lambda * diag(n_all)
    d_sc    <- sqrt(diag(cor_full))
    cor_full <- cor_full / outer(d_sc, d_sc)
  }
  suffStat <- list(C = cor_full, n = nrow(all_data))
  labels   <- colnames(all_data)
  
  pc_fit <- pc(suffStat = suffStat, indepTest = gaussCItest,
               labels = labels, alpha = alpha,
               verbose = FALSE, maj.rule = TRUE, solve.confl = TRUE)
  amat_full <- as(pc_fit, "amat")
  rownames(amat_full) <- colnames(amat_full) <- labels
  
  # Induced subgraph on X only.
  amat_full[fg_names, fg_names]
}

## --- Metrics -------------------------------------------------------------

compute_metrics <- function(amat_est, amat_true) {
  # amat_true: d_x x d_x ancestral matrix (transitive closure of direct
  # edges), upper-triangular by construction. amat_est: ASCEND-style
  # adjacency from a PC variant.
  d <- nrow(amat_true)
  fg_names <- rownames(amat_true)
  amat_est <- amat_est[fg_names, fg_names]
  
  pairs_records <- list()
  for (i in 1:(d - 1)) for (j in (i + 1):d) {
    t_ij <- amat_true[i, j]  # 1 if i is ancestor of j, else 0
    t_ji <- amat_true[j, i]  # always 0 by upper-triangular structure
    if (is.na(t_ij)) next
    a_ij <- amat_est[i, j]
    a_ji <- amat_est[j, i]
    pred_dir_ij <- as.integer(!is.na(a_ij) && !is.na(a_ji) &&
                                a_ij == 1 && a_ji == 0)
    pred_dir_ji <- as.integer(!is.na(a_ij) && !is.na(a_ji) &&
                                a_ji == 1 && a_ij == 0)
    pred_und    <- as.integer(!is.na(a_ij) && !is.na(a_ji) &&
                                a_ij == 1 && a_ji == 1)
    pred_edge   <- pred_dir_ij | pred_dir_ji | pred_und
    edge_score  <- pred_dir_ij * 2L + pred_und * 1L
    
    pairs_records[[length(pairs_records) + 1]] <- data.frame(
      i = i, j = j,
      true_label = as.integer(t_ij == 1 || t_ji == 1),
      true_dir_ij = as.integer(t_ij == 1),
      true_dir_ji = as.integer(t_ji == 1),
      pred_dir_ij = pred_dir_ij,
      pred_dir_ji = pred_dir_ji,
      pred_und    = pred_und,
      pred_edge   = pred_edge,
      edge_score  = edge_score
    )
  }
  pairs <- bind_rows(pairs_records)
  
  n_true <- sum(pairs$true_label)
  prev   <- if (nrow(pairs) > 0) n_true / nrow(pairs) else 0
  
  if (n_true == 0 || nrow(pairs) == 0) {
    return(list(n_true = n_true, prev = prev,
                auprc_ratio = NA, epr = NA, dir_acc = NA,
                f1 = NA, precision = NA, recall = NA,
                fisher_or = NA, fisher_p = NA,
                TP = 0, FP = 0, FN = 0, TN = 0,
                n_predicted = 0, n_directed = 0))
  }
  
  pr <- pr.curve(scores.class0 = pairs$edge_score[pairs$true_label == 1],
                 scores.class1 = pairs$edge_score[pairs$true_label == 0],
                 curve = FALSE)
  auprc_ratio <- pr$auc.integral / prev
  
  k <- n_true
  top_k <- order(pairs$edge_score, decreasing = TRUE)[1:k]
  epr <- (sum(pairs$true_label[top_k]) / k) / prev
  
  # Direction accuracy among true positives.
  tp_df <- pairs[pairs$pred_edge == 1 & pairs$true_label == 1, , drop = FALSE]
  if (nrow(tp_df) > 0) {
    correct <- (tp_df$pred_dir_ij == 1 & tp_df$true_dir_ij == 1) |
      (tp_df$pred_dir_ji == 1 & tp_df$true_dir_ji == 1)
    dir_acc <- mean(correct)
  } else {
    dir_acc <- NA
  }
  
  # Binary classification on edge presence.
  TP <- sum(pairs$pred_edge == 1 & pairs$true_label == 1)
  FP <- sum(pairs$pred_edge == 1 & pairs$true_label == 0)
  FN <- sum(pairs$pred_edge == 0 & pairs$true_label == 1)
  TN <- sum(pairs$pred_edge == 0 & pairs$true_label == 0)
  precision <- if (TP + FP > 0) TP / (TP + FP) else NA
  recall    <- if (TP + FN > 0) TP / (TP + FN) else NA
  f1 <- if (!is.na(precision) && !is.na(recall) && precision + recall > 0)
    2 * precision * recall / (precision + recall) else NA
  
  fr <- tryCatch(
    fisher.test(matrix(c(TP, FP, FN, TN), nrow = 2),
                alternative = "greater"),
    error = function(e) list(estimate = NA, p.value = NA)
  )
  
  n_predicted <- sum(pairs$pred_edge)
  n_directed  <- sum(pairs$pred_dir_ij | pairs$pred_dir_ji)
  
  list(n_true = n_true, prev = prev,
       auprc_ratio = auprc_ratio, epr = epr, dir_acc = dir_acc,
       f1 = f1, precision = precision, recall = recall,
       fisher_or = unname(fr$estimate), fisher_p = fr$p.value,
       TP = TP, FP = FP, FN = FN, TN = TN,
       n_predicted = n_predicted, n_directed = n_directed)
}

## --- Per-rep runner ------------------------------------------------------

run_one_rep <- function(n_samp, rep_id, seed) {
  cat(sprintf("  [n=%d rep=%d] sim_dat\n", n_samp, rep_id))
  set.seed(seed)
  sim_obj <- sim_dat(n        = n_samp,
                     d_z      = D_Z,
                     d_x      = D_X,
                     r2       = SIM_PARAMS$r2,
                     lin_pr   = SIM_PARAMS$lin_pr,
                     sp       = SIM_PARAMS$sp,
                     p_cross  = SIM_PARAMS$p_cross,
                     x_effect = SIM_PARAMS$x_effect,
                     seed     = seed)
  amat_true <- get_true_ancestral_matrix(sim_obj$adj_xx)
  
  # ASCEND-PC.
  cat(sprintf("  [n=%d rep=%d] ASCEND-PC\n", n_samp, rep_id))
  t_a <- system.time({
    out_a <- tryCatch(run_ascend_pc(sim_obj),
                      error = function(e) {
                        message("ASCEND-PC error: ", conditionMessage(e))
                        NULL
                      })
  })
  
  # Plain PC X-only.
  cat(sprintf("  [n=%d rep=%d] Plain PC (X)\n", n_samp, rep_id))
  t_px <- system.time({
    amat_px <- tryCatch(run_plain_pc_x(sim_obj),
                        error = function(e) {
                          message("Plain PC (X) error: ", conditionMessage(e))
                          NULL
                        })
  })
  
  # Plain PC Z+X.
  cat(sprintf("  [n=%d rep=%d] Plain PC (Z+X)\n", n_samp, rep_id))
  t_pzx <- system.time({
    amat_pzx <- tryCatch(run_plain_pc_zx(sim_obj),
                         error = function(e) {
                           message("Plain PC (Z+X) error: ", conditionMessage(e))
                           NULL
                         })
  })
  
  # Metrics for each method that ran.
  m_a <- if (!is.null(out_a))   compute_metrics(out_a$amat, amat_true) else NULL
  m_px <- if (!is.null(amat_px))  compute_metrics(amat_px,   amat_true) else NULL
  m_pzx <- if (!is.null(amat_pzx)) compute_metrics(amat_pzx, amat_true) else NULL
  
  make_row <- function(method, m, t_sec, extra = list()) {
    if (is.null(m)) return(data.frame(
      n = n_samp, rep = rep_id, method = method,
      auprc_ratio = NA, epr = NA, dir_acc = NA, f1 = NA,
      precision = NA, recall = NA,
      fisher_or = NA, fisher_p = NA,
      TP = NA, FP = NA, FN = NA, TN = NA,
      n_predicted = NA, n_directed = NA,
      n_true_pairs = NA, prevalence = NA,
      wall_clock_s = t_sec,
      n_pruned = NA, n_pairs_total = NA, prune_pct = NA,
      stringsAsFactors = FALSE
    ))
    data.frame(
      n = n_samp, rep = rep_id, method = method,
      auprc_ratio = m$auprc_ratio, epr = m$epr, dir_acc = m$dir_acc,
      f1 = m$f1, precision = m$precision, recall = m$recall,
      fisher_or = m$fisher_or, fisher_p = m$fisher_p,
      TP = m$TP, FP = m$FP, FN = m$FN, TN = m$TN,
      n_predicted = m$n_predicted, n_directed = m$n_directed,
      n_true_pairs = m$n_true, prevalence = m$prev,
      wall_clock_s = t_sec,
      n_pruned     = if (!is.null(extra$n_pruned))      extra$n_pruned      else NA,
      n_pairs_total= if (!is.null(extra$n_pairs_total)) extra$n_pairs_total else NA,
      prune_pct    = if (!is.null(extra$prune_pct))     extra$prune_pct     else NA,
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(
    make_row("ASCEND-PC", m_a, as.numeric(t_a["elapsed"]),
             extra = list(
               n_pruned      = if (!is.null(out_a)) out_a$n_gaps        else NA,
               n_pairs_total = if (!is.null(out_a)) out_a$n_pairs_total else NA,
               prune_pct     = if (!is.null(out_a)) out_a$prune_pct     else NA)),
    make_row("Plain PC (X)",   m_px,  as.numeric(t_px["elapsed"])),
    make_row("Plain PC (Z+X)", m_pzx, as.numeric(t_pzx["elapsed"]))
  )
}

## --- Main loop with checkpointing ---------------------------------------

# Resume from checkpoint if present.
if (file.exists(CHECKPOINT_PATH)) {
  cp <- readRDS(CHECKPOINT_PATH)
  results <- cp$results
  done_keys <- cp$done_keys
  cat(sprintf("Resuming from checkpoint: %d reps already completed\n",
              length(done_keys)))
} else {
  results <- list()
  done_keys <- character(0)
}

# Master seed for reproducibility; per-rep seeds derived from it.
MASTER_SEED <- 20240101

total_cells <- length(N_LEVELS) * N_REPS
cell_idx <- 0

for (n_samp in N_LEVELS) {
  for (rep_id in seq_len(N_REPS)) {
    cell_idx <- cell_idx + 1
    key <- sprintf("n=%d_rep=%d", n_samp, rep_id)
    if (key %in% done_keys) {
      cat(sprintf("[%d/%d] %s already done, skipping\n",
                  cell_idx, total_cells, key))
      next
    }
    
    seed <- MASTER_SEED + 1000 * which(N_LEVELS == n_samp) + rep_id
    cat(sprintf("\n[%d/%d] %s (seed=%d)\n",
                cell_idx, total_cells, key, seed))
    
    rep_df <- tryCatch(
      run_one_rep(n_samp, rep_id, seed),
      error = function(e) {
        message("Rep error: ", conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(rep_df)) {
      results[[length(results) + 1]] <- rep_df
      done_keys <- c(done_keys, key)
      # Save checkpoint after each rep.
      saveRDS(list(results = results, done_keys = done_keys),
              CHECKPOINT_PATH)
      cat(sprintf("  checkpoint saved (%d reps complete)\n",
                  length(done_keys)))
    } else {
      cat(sprintf("  rep failed, not saved\n"))
    }
  }
}

## --- Finalise outputs ----------------------------------------------------

all_results <- bind_rows(results)
saveRDS(all_results, RESULTS_RDS)
write.csv(all_results, RESULTS_CSV, row.names = FALSE)

cat(sprintf("\nWrote %s (%d rows)\n", RESULTS_CSV, nrow(all_results)))

# Quick summary.
summ <- all_results %>%
  group_by(n, method) %>%
  summarise(
    auprc_ratio_mean = mean(auprc_ratio, na.rm = TRUE),
    auprc_ratio_sd   = sd(auprc_ratio,   na.rm = TRUE),
    epr_mean         = mean(epr, na.rm = TRUE),
    dir_acc_mean     = mean(dir_acc, na.rm = TRUE),
    f1_mean          = mean(f1, na.rm = TRUE),
    wall_mean_s      = mean(wall_clock_s, na.rm = TRUE),
    n_pruned_mean    = mean(n_pruned, na.rm = TRUE),
    n_pairs_total    = mean(n_pairs_total, na.rm = TRUE),
    prune_mean       = mean(prune_pct, na.rm = TRUE),
    n_reps_ok        = sum(!is.na(auprc_ratio)),
    .groups = "drop"
  )
cat("\n=== Summary ===\n")
print(as.data.frame(summ))

# Headline pruning statement for ASCEND-PC.
asc <- summ %>% filter(method == "ASCEND-PC")
if (nrow(asc) > 0) {
  cat("\n=== ASCEND pruning summary ===\n")
  for (k in seq_len(nrow(asc))) {
    r <- asc[k, ]
    cat(sprintf("n=%d:  ASCEND pruned %d of %d pairs (%.1f%%) before PC\n",
                r$n,
                as.integer(round(r$n_pruned_mean)),
                as.integer(round(r$n_pairs_total)),
                r$prune_mean))
  }
}

cat("\nDone.\n")