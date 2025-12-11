###############################################################################
# FULL EXPERIMENT PIPELINE FOR:
#   ASCEND, GENIE3, TIGRESS, ARACNE, CLR
#   → Immediate-edge DAG evaluation
#   → Repeated runs across sample sizes
#   → F1 (skeleton), orientation accuracy, SHD
###############################################################################

library(data.table)
library(igraph)
library(ggplot2)
library(matrixStats)

try(library(dplyr), silent = TRUE)

#### ================================================================
#### Utilities
#### ================================================================
clean_adj <- function(A) {
  A <- as.matrix(A)
  A[is.na(A)] <- 0
  (A != 0) * 1L
}

transitive_reduction <- function(M) {
  M <- clean_adj(M)
  g <- graph_from_adjacency_matrix(M, mode="directed")
  red <- igraph::transitive.reduction(g)
  A <- as.matrix(as_adj(red))
  (A != 0) * 1L
}

topk_from_weighted <- function(W, K, xlabs) {
  W <- as.matrix(W); diag(W) <- 0
  idx <- which(W != 0, arr.ind = TRUE)
  
  if (length(idx) == 0)
    return(matrix(0, length(xlabs), length(xlabs), dimnames=list(xlabs,xlabs)))
  
  vals <- W[idx]
  ord <- order(-vals)
  keep <- head(ord, K)
  
  A <- matrix(0, length(xlabs), length(xlabs), dimnames=list(xlabs,xlabs))
  for (k in keep) {
    i <- idx[k,1]; j <- idx[k,2]
    A[i,j] <- 1
  }
  A
}

#### ================================================================
#### Evaluation — immediate edges
#### ================================================================
evaluate_immediate <- function(A_true, A_est) {
  A_true <- clean_adj(A_true)
  A_est  <- clean_adj(A_est)
  
  S_true <- (A_true + t(A_true) > 0) * 1L
  S_est  <- (A_est  + t(A_est)  > 0) * 1L
  
  tp <- sum(S_true==1 & S_est==1)
  fp <- sum(S_true==0 & S_est==1)
  fn <- sum(S_true==1 & S_est==0)
  
  precision <- if ((tp+fp)==0) NA else tp/(tp+fp)
  recall    <- if ((tp+fn)==0) NA else tp/(tp+fn)
  f1        <- if (is.na(precision) || is.na(recall) || precision+recall==0)
    NA else 2*precision*recall/(precision+recall)
  
  idx <- which(S_true==1 & S_est==1, arr.ind=TRUE)
  good <- 0
  for (k in seq_len(nrow(idx))) {
    i <- idx[k,1]; j <- idx[k,2]
    if (A_true[i,j] == A_est[i,j] && A_true[j,i]==A_est[j,i])
      good <- good + 1
  }
  orient_acc <- if (nrow(idx)==0) NA else good/nrow(idx)
  
  list(
    f1_skel = f1,
    orient_accuracy = orient_acc,
    shd = sum(A_true != A_est)
  )
}



safe_run <- function(expr, p) {
  tryCatch(
    expr,
    error = function(e) {
      message("[METHOD ERROR] ", e$message)
      return(matrix(0, nrow=p, ncol=p))
    }
  )
}


#### ================================================================
#### ASCEND wrapper  (YOU MUST PROVIDE ascend_fn)
#### ================================================================
wrapper_ascend <- function(sim) {
  p <- sim$params$d_x
  xlabs <- paste0("x", 1:p)
  X <- as.matrix(sim$dat[, ..xlabs])
  
  safe_run({
    out <- ascend::ascend(X)
    
    if (!is.list(out) || is.null(out$dag))
      stop("ASCEND returned invalid structure")
    
    A <- out$dag
    A[is.na(A)] <- 0
    A
  }, p)
}



#### ================================================================
#### GENIE3
#### ================================================================
wrapper_genie3 <- function(sim) {
  if (!requireNamespace("GENIE3", quietly=TRUE)) stop("GENIE3 not installed")
  
  p <- sim$params$d_x
  xlabs <- paste0("x", 1:p)
  X <- as.matrix(sim$dat[, ..xlabs])
  
  safe_run({
    w <- GENIE3::GENIE3(X)
    A <- GENIE3::getLinkMatrix(w)
    A
  }, p)
}



#### ================================================================
#### TIGRESS
#### ================================================================
wrapper_tigress <- function(sim) {
  if (!requireNamespace("tigress", quietly=TRUE)) stop("tigress not installed")
  
  p <- sim$params$d_x
  xlabs <- paste0("x", 1:p)
  X <- as.matrix(sim$dat[, ..xlabs])
  Xt <- t(X)
  
  safe_run({
    f <- tigress::tigress
    argnames <- names(formals(f))
    
    if (all(c("nt","R") %in% argnames)) {
      res <- f(Xt, nt=100, R=100)
    } else {
      message("[TIGRESS] old API, using default parameters")
      res <- f(Xt)
    }
    
    A <- res$score
    rownames(A) <- colnames(A) <- xlabs
    A
  }, p)
}




#### ================================================================
#### ARACNE + CLR
#### ================================================================
wrapper_aracne <- function(sim) {
  if (!requireNamespace("minet", quietly=TRUE)) stop("minet not installed")
  
  p <- sim$params$d_x
  xlabs <- paste0("x", 1:p)
  X <- as.matrix(sim$dat[, ..xlabs])
  
  safe_run({
    A <- minet::minet(X, method="aracne")
    A
  }, p)
}



wrapper_clr <- function(sim) {
  if (!requireNamespace("minet", quietly=TRUE)) stop("minet not installed")
  
  p <- sim$params$d_x
  xlabs <- paste0("x", 1:p)
  X <- as.matrix(sim$dat[, ..xlabs])
  
  safe_run({
    A <- minet::minet(X, method="clr")
    A
  }, p)
}



#### ================================================================
#### Run a single method-set on a single dataset
#### ================================================================
run_all_methods <- function(sim) {
  p <- sim$params$d_x
  Gtrue <- clean_adj(sim$adj_mat)
  
  methods <- list(
    ASCEND = wrapper_ascend(sim),
    GENIE3 = wrapper_genie3(sim),
    TIGRESS = wrapper_tigress(sim),
    ARACNE = wrapper_aracne(sim),
    CLR = wrapper_clr(sim)
  )
  
  out <- list()
  for (m in names(methods)) {
    out[[m]] <- evaluate_immediate_edges(Gtrue, methods[[m]])
  }
  
  out
}


#### ================================================================
#### Experiment: repeated runs at multiple sample sizes
#### ================================================================
run_sample_size_experiment <- function(
    n_vec,
    n_rep,
    d_z=10, d_x=5,
    rho=0.5, r2=0.5,
    lin_pr=1, sp=0.5,
    method="er", pref=NA, p_cross=0.05,
    seed_base=100
) {
  
  results <- list()
  
  for (n in n_vec) {
    message("Sample size: ", n)
    
    mets <- list()
    
    for (r in 1:n_rep) {
      set.seed(seed_base + r + n)
      sim <- sim_dat(n, d_z, d_x, rho, r2, lin_pr, sp, method, pref, p_cross)
      
      out <- run_all_methods(sim)
      mets[[r]] <- out
    }
    
    results[[as.character(n)]] <- mets
  }
  res <- results
  return(res)
}


#### ================================================================
#### Plotting
#### ================================================================
plot_results <- function(res) {
  
  p1 <- ggplot(res, aes(x=n, y=f1_skel, color=method)) +
    stat_summary(fun=mean, geom="line") +
    stat_summary(fun=mean, geom="point") +
    theme_bw() + ggtitle("Skeleton F1 vs sample size")
  
  p2 <- ggplot(res, aes(x=n, y=orient_accuracy, color=method)) +
    stat_summary(fun=mean, geom="line") +
    stat_summary(fun=mean, geom="point") +
    theme_bw() + ggtitle("Orientation accuracy vs sample size")
  
  p3 <- ggplot(res, aes(x=n, y=shd, color=method)) +
    stat_summary(fun=mean, geom="line") +
    stat_summary(fun=mean, geom="point") +
    theme_bw() + ggtitle("SHD vs sample size")
  
  list(F1=p1, ORIENT=p2, SHD=p3)
}

###############################################################################
# END PIPELINE
###############################################################################


n_vec <- c(500, 1000, 2000)
res <- run_sample_size_experiment(
  n_vec = n_vec,
  n_rep = 10,
  d_z = 20, d_x = 5,
  rho = 0.5, r2 = 0.5,
  lin_pr = 1, sp = 0.5,
  method = "er",
  pref = NA,
  p_cross = 0.08,
  seed_base = 123
)

plots <- plot_results(res$raw)
plots$F1
plots$ORIENT
plots$SHD


