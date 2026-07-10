### Simulations for subgraph discovery with many ancestors ###
# https://github.com/dswatson/cbl
# Load libraries and Shah's code, register cores
source('shah_ss.R')
library(data.table)
library(pcalg)
library(RBGL)
library(matrixStats)
library(glmnet)
library(lightgbm)
library(tidyverse)
library(doMC)
registerDoMC(16)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

################################################################################

### SIMULATION ###

#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param d_x Dimensionality of X.
#' @param rho Auto-correlation of the Toeplitz matrix for Z.
#' @param r2 Proportion of variance explained for all foreground variables. 
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' @param method Method used for generating the graph structure. Options are
#'   \code{"er"} for Erdós-Rényi and \code{"barabasi"} for Barabási-Albert.
#' @param sp Average sparsity of the graph. Note that this must be high for 
#'   \code{method = "barabasi"} or else you'll run into errors.
#' @param pref Strength of preferential attachment if \code{method = "barabasi"}.
#' 



################################################################################

### CONFOUNDER BLANKET LEARNER ###


#' @param x Design matrix.
#' @param y Outcome vector.
#' @param f Regression method, either \code{"lasso"} or \code{"gbm"}.
#' @param prms List of parameters to use when \code{f = "gbm"}.
#' 

# Fit regressions, return bit vector for feature selection.
l0 <- function(x, y, f, prms) {
  n <- nrow(x)
  trn <- sample(n, round(0.8 * n))
  tst <- seq_len(n)[-trn]
  if (f == 'lasso') {
    fit <- glmnet(x[trn, ], y[trn], intercept = FALSE)
    y_hat <- predict(fit, newx = x[tst, ], s = fit$lambda)
    eps <- y_hat - y[tst]
    mse <- colMeans(eps^2)
    betas <- coef(fit, s = fit$lambda)[-1, which.min(mse)]
    out <- ifelse(betas == 0, 0, 1)
  } else if (f == 'gbm') {
    d_trn <- lgb.Dataset(x[trn, ], label = y[trn])
    d_tst <- lgb.Dataset.create.valid(d_trn, x[tst, ], label = y[tst])
    fit <- lgb.train(params = prms, data = d_trn, valids = list(tst = d_tst), 
                     nrounds = 3500, early_stopping_rounds = 10, verbose = 0)
    vimp <- lgb.importance(fit)
    out <- as.numeric(colnames(x) %in% vimp$Feature)
  } 
  return(out)
}


#' @param df Table of (de)activation rates.
#' @param B Number of complementary pairs to draw for stability selection.

# Compute consistency lower bound
epsilon_fn <- function(df, B) {
  # Nullify 
  dji <- drji <- aji <- arji <- dij <- drij <- aij <- arij <- tau <- tt <-
    int_err <- ext_err <- NULL
  # Loop through thresholds in search of inconsistencies
  err_check <- function(tau) {
    # Inferences at this value of tau
    df[, dji := ifelse(drji >= tau, 1, 0)]
    df[, aji := ifelse(arji >= tau, 1, 0)]
    df[, dij := ifelse(drij >= tau, 1, 0)]
    df[, aij := ifelse(arij >= tau, 1, 0)]
    # Internal consistency (for a single Z)
    df[, int_err := ifelse(dji + aji + dij + aij > 1, 1, 0)]
    int_err <- ifelse(sum(df$int_err) > 0, 1, 0)
    # External consistency (across multiple Z's)
    if (df[, sum(dji) > 0] & df[, sum(dij + aij) > 0]) {
      ext_err <- 1
    } else if (df[, sum(dij) > 0] & df[, sum(dji + aji) > 0]) {
      ext_err <- 1
    } else {
      ext_err <- 0
    }
    # Export
    out <- data.table('tau' = tau, 'int_err' = int_err, 'ext_err' = ext_err)
    return(out)
  }
  err_df <- foreach(tt = seq_len(2 * B) / (2 * B), .combine = rbind) %do% 
    err_check(tt)
  # Compute minimal thresholds, export
  epsilon <- err_df[int_err == 0 & ext_err == 0, min(tau)]
  return(epsilon)
}


#' @param df Table of (de)activation rates.
#' @param eps Consistency lower bound, as computed by \code{epsilon_fn}.
#' @param order Causal order of interest, either \code{"ij"} or \code{"ji"}.
#' @param rule Inference rule, either \code{"R1"} or \code{"R2"}.
#' @param B Number of complementary pairs to draw for stability selection.

# Infer causal direction using stability selection
ss_fn <- function(df, eps, order, rule, B) {
  # Find the right rate
  if (order == 'ji' & rule == 'R1') {
    r <- df[, drji]
  } else if (order == 'ji' & rule == 'R2') {
    r <- df[, arji]
  } else if (order == 'ij' & rule == 'R1') {
    r <- df[, drij]
  } else if (order == 'ij' & rule == 'R2') {
    r <- df[, arij]
  } 
  # Stability selection parameters
  theta <- mean(r)
  ub <- minD(theta, B) * sum(r <= theta)
  tau <- seq_len(2 * B) / (2 * B)
  # Do any features exceed the upper bound?
  dat <- data.frame(tau, err_bound = ub) %>%
    filter(tau >= eps) %>%
    rowwise() %>%
    mutate(detected = sum(r >= tau)) %>% 
    ungroup(.) %>%
    mutate(surplus = ifelse(detected > err_bound, 1, 0))
  # Export
  out <- data.table(
    'order' = order, 'rule' = rule, 
    'decision' = ifelse(sum(dat$surplus) > 0, 1, 0)
  )
  return(out)
}

#' @param b Subsample index.
#' @param i First foreground variable index.
#' @param j Second foreground variable index.
#' @param x Matrix of foreground variables.
#' @param z_t Intersection of iteration-t known non-descendants for foreground
#'   variables \code{i} and \code{j}.
#' @param s Regression method. 
#' @param params Optional list of parameters to use when \code{s = "boost"}.

# Complementary pairs subsampling loop
sub_loop <- function(b, i, j, x, z_t, s, params) {
  # Prelimz
  n <- nrow(x) 
  d_zt <- ncol(z_t)
  # Take complementary subsets
  a_set <- sample(n, round(0.5 * n))
  b_set <- seq_len(n)[-a_set]
  # Fit reduced models
  s0 <- sapply(c(i, j), function(k) {
    c(l0(z_t[a_set, ], x[a_set, k], s, params), 
      l0(z_t[b_set, ], x[b_set, k], s, params))
  })
  # Fit expanded models
  s1 <- sapply(c(i, j), function(k) {
    not_k <- setdiff(c(i, j), k)
    c(l0(cbind(z_t[a_set, ], x[a_set, not_k]), x[a_set, k], s, params),
      l0(cbind(z_t[b_set, ], x[b_set, not_k]), x[b_set, k], s, params))
  })
  # Record disconnections and (de)activations
  dis_a <- any(s1[d_zt + 1, ] == 0)
  dis_b <- any(s1[2 * (d_zt + 1), ] == 0)
  dis <- rep(c(dis_a, dis_b), each = d_zt)
  d_ji <- s0[, 1] == 1 & s1[seq_len(d_zt), 1] == 0
  a_ji <- s0[, 2] == 0 & s1[seq_len(d_zt), 2] == 1
  d_ij <- s0[, 2] == 1 & s1[seq_len(d_zt), 2] == 0
  a_ij <- s0[, 1] == 0 & s1[seq_len(d_zt), 1] == 1
  # Export
  out <- data.table(b = rep(c(2 * b - 1, 2 * b), each = d_zt), i, j,
                    z = rep(colnames(z_t), times = 2),
                    dis, d_ji, a_ji, d_ij, a_ij)
  return(out)
}


#' @param sim_obj Simulation object as computed by \code{sim_dat}.
#' @param maxiter Maximum number of iterations to loop through if convergence
#'   is elusive.
#' @param gamma Omission threshold.
#' @param B Number of complementary pairs to draw for stability selection.

# Subdag discovery via confounder blanket regression
cbl_fn <- function(sim_obj, gamma = 0.5, maxiter = 100, B = 50) {
  ### PRELIMINARIES ###
  # Get data, hyperparameters, train/test split
  dat <- sim_obj$dat
  n <- nrow(dat)
  z <- as.matrix(select(dat, starts_with('z')))
  d_z <- ncol(z)
  x <- as.matrix(select(dat, starts_with('x')))
  d_x <- ncol(x)
  xlabs <- paste0('x', seq_len(d_x))
  linear <- sim_obj$params$lin_pr == 1
  if (linear) {
    f <- 'lasso'
    prms <- NULL
  } else {
    f <- 'gbm'
    prms <- list(
      objective = 'regression', max_depth = 1, 
      bagging.fraction = 0.5, feature_fraction = 0.8, 
      num_threads = 1, force_col_wise = TRUE
    )
  }
  # Initialize
  adj_list <- list(
    matrix(NA_real_, nrow = d_x, ncol = d_x, 
           dimnames = list(colnames(x), colnames(x)))
  )
  adj0 <- adj1 <- adj_list[[1]]
  converged <- FALSE
  iter <- 0
  ### BIG LOOP ###
  while(converged == FALSE & iter <= maxiter) {
    # Pairwise test loop
    for (i in 2:d_x) {
      for (j in 1:(i - 1)) {
        # Only continue if relationship is unknown
        if (is.na(adj1[i, j]) & is.na(adj1[j, i])) { 
          preceq_i <- which(adj0[i, ] > 0)
          preceq_j <- which(adj0[j, ] > 0)
          a0 <- intersect(preceq_i, preceq_j) 
          preceq_i <- which(adj1[i, ] > 0)
          preceq_j <- which(adj1[j, ] > 0)
          a1 <- intersect(preceq_i, preceq_j) 
          # Only continue if the set of non-descendants has increased since last 
          # iteration (i.e., have we learned anything new?)
          if (iter == 0 | length(a1) > length(a0)) {
            z_t <- cbind(z, x[, a1])
            df <- foreach(bb = seq_len(B), .combine = rbind) %do%
              sub_loop(bb, i, j, x, z_t, f, prms)
            # Compute rates
            df[, disr := sum(dis) / .N]
            if (df$disr[1] > gamma) { 
              adj1[i, j] <- adj1[j, i] <- 0
            } else {
              df[, drji := sum(d_ji) / .N, by = z]
              df[, arji := sum(a_ji) / .N, by = z]
              df[, drij := sum(d_ij) / .N, by = z]
              df[, arij := sum(a_ij) / .N, by = z]
              df <- unique(df[, .(i, j, z, disr, drji, arji, drij, arij)])
              # Consistent lower bound
              eps <- epsilon_fn(df, B)
              # Stable upper bound
              out <- foreach(oo = c('ji', 'ij'), .combine = rbind) %:%
                foreach(rr = c('R1', 'R2'), .combine = rbind) %do%
                ss_fn(df, eps, oo, rr, B)
              # Update adjacency matrix
              if (out[rule == 'R1' & order == 'ji', decision == 1]) {
                adj1[i, j] <- 1
                adj1[j, i] <- 0
              } else if (out[rule == 'R1' & order == 'ij', decision == 1]) {
                adj1[j, i] <- 1
                adj1[i, j] <- 0
              } else if (out[rule == 'R2', sum(decision) == 2]) {
                adj1[i, j] <- adj1[j, i] <- 0
              } else if (out[rule == 'R2' & order == 'ji', decision == 1]) {
                adj1[i, j] <- 0.5 
              } else if (out[rule == 'R2' & order == 'ij', decision == 1]) {
                adj1[j, i] <- 0.5
              } 
            }
          }
        } 
      }
    }
    # Check for transitivity - FIXED VERSION
    # Check for transitivity - SIMPLE FIX
    closure <- FALSE
    while(closure == FALSE) {
      closure <- TRUE
      for (i in seq_len(d_x)) {
        m <- which(adj1[, i] == 1)  # direct parents of i
        if (length(m) > 0) {
          # For each parent, get its parents
          for (parent in m) {
            e <- which(adj1[, parent] == 1)  # grandparents through this parent
            if (length(e) > 0) {
              for (grandparent in e) {
                if (is.na(adj1[grandparent, i]) || adj1[grandparent, i] != 1) {
                  adj1[grandparent, i] <- 1
                  closure <- FALSE
                }
              }
            }
          }
        }
      }
    }
    # Iterate, check for convergence
    iter <- iter + 1
    if (identical(adj0, adj1)) {
      converged <- TRUE
    } else {
      adj_list <- append(adj_list, list(adj1))
      adj0 <- adj_list[[iter]]
      adj1 <- adj_list[[iter + 1]]
      adj_list[[iter]] <- 0
    }
  }
  # Export final adjacency matrix
  #adj1 <- permute_to_lower(adj1)
  
  return(adj1)
}

################################################################################

### EXAMPLE ###
#
# Guarded the same way as ascend.R's example block: sys.nframe() == 0 only
# holds when this file is run as the top-level script (e.g. `Rscript cbl.R`),
# not when it's source()'d from another script such as the benchmark. This
# also means cbl_fn() is always exercised against whatever sim_dat() is
# in scope at call time -- when sourced from the benchmark, that's the same
# sim_dat() (and hence the same simulated dataset per cell) that ASCEND runs
# against, since sim_dat() is defined once in ascend.R.
#
# NOTE: this demo call uses sim_dat()'s CURRENT signature (n, d_z, d_x, r2,
# lin_pr, sp, p_cross, x_effect, seed) as defined in ascend.R. The previous
# version of this block called an older simulator with rho/method/pref
# arguments that no longer exist.

if (sys.nframe() == 0) {
  sim_obj <- sim_dat(
    n        = 1000,  # sample size
    d_z      = 20,    # number of z
    d_x      = 6,     # number of observed variables
    r2       = 0.4,   # variance explained by instruments
    lin_pr   = 1,      # linearity
    sp       = 0.5,   # sparsity
    p_cross  = 0.05,  # P(Z -> X edge)
    x_effect = 0.8,   # X -> X effect size
    seed     = 123
  )
  
  # Run cbl algorithm
  amat_cbl <- cbl_fn(sim_obj)
  
  print("amat_cbl")
  print(amat_cbl)
}