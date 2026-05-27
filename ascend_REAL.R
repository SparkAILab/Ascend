# ======================================================================
# ASCEND - Complete Self-Contained Implementation
# Ancestral Scalable Causal discovEry via iNherited Descent
#
# Provides: sim_dat(), ascend_fn(), get_true_ancestral_matrix(),
#           evaluate_ancestral()
#
# No top-level side effects. Safe to source() in any environment.
# ======================================================================


# ======================================================================
# 1. Libraries
# ======================================================================
suppressPackageStartupMessages({
  library(bnlearn)
  library(dplyr)
  library(data.table)
  library(matrixStats)
  library(igraph)
  library(doMC)
})


# ======================================================================
# 2. Simulation function
#
# Returns:
#   $dat      - data.table, columns z1..z_dz, x1..x_dx
#   $adj_full - full (dz+dx)x(dz+dx) direct-edge adjacency
#   $adj_xx   - dx x dx direct X->X edges, adj_xx[child,parent]=1
#   $params   - parameters used
# ======================================================================

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
  
  prep <- function(dat, pr) {
    if (is.null(dim(dat))) dat <- matrix(dat, ncol=1) else dat <- as.matrix(dat)
    p <- ncol(dat); if (p==0 || pr>=1) return(dat)
    n_nl <- round((1-pr)*p); if (n_nl<=0) return(dat)
    cols <- sample.int(p, size=n_nl)
    nl_types <- sample(c('sq','sqrt','sftpls','relu'), n_nl, replace=TRUE)
    out <- dat
    for (ii in seq_along(cols)) {
      idx <- cols[ii]; v <- dat[,idx]
      if      (nl_types[ii]=='sq')     out[,idx] <- v^2
      else if (nl_types[ii]=='sqrt')   out[,idx] <- sqrt(abs(v))
      else if (nl_types[ii]=='sftpls') out[,idx] <- log1p(exp(v))
      else if (nl_types[ii]=='relu')   out[,idx] <- pmax(v, 0)
    }
    out
  }
  
  sim_noise <- function(signal, r2) {
    signal <- as.numeric(signal); var_mu <- var(signal, na.rm=TRUE)
    if (is.na(var_mu) || var_mu==0) return(rnorm(n, sd=1))
    if (r2<=0) return(rnorm(n, sd=sqrt(var_mu)))
    rnorm(n, sd=sqrt(pmax(0, (var_mu - r2*var_mu) / r2)))
  }
  
  make_dag_upper <- function(k, p_edge=NULL) {
    if (is.null(p_edge)) p_edge <- min(max(1-sp, 0.01), 0.4)
    A <- matrix(0, k, k)
    for (i in 1:(k-1)) for (j in (i+1):k) if (runif(1) < p_edge) A[i,j] <- 1
    perm <- sample.int(k); list(A=A[perm,perm], perm=perm)
  }
  
  z_dag <- make_dag_upper(d_z, p_z); Azz <- z_dag$A
  rownames(Azz) <- colnames(Azz) <- paste0("z", seq_len(d_z))
  x_dag <- make_dag_upper(d_x, p_x); Axx <- x_dag$A
  rownames(Axx) <- colnames(Axx) <- paste0("x", seq_len(d_x))
  
  Azx <- matrix(0, d_z, d_x,
                dimnames=list(paste0("z",1:d_z), paste0("x",1:d_x)))
  for (i in seq_len(d_z)) for (j in seq_len(d_x))
    if (runif(1) < p_cross) Azx[i,j] <- 1
  
  fn <- c(rownames(Azz), colnames(Axx))
  A_full <- matrix(0, d_z+d_x, d_z+d_x, dimnames=list(fn,fn))
  A_full[1:d_z, 1:d_z] <- Azz
  A_full[1:d_z, (d_z+1):(d_z+d_x)] <- Azx
  A_full[(d_z+1):(d_z+d_x), (d_z+1):(d_z+d_x)] <- Axx
  
  gZ    <- igraph::graph_from_adjacency_matrix(Azz, mode="directed")
  topoZ <- as.integer(igraph::topo_sort(gZ, mode="out"))
  Zmat  <- matrix(0, n, d_z, dimnames=list(NULL, paste0("z",seq_len(d_z))))
  for (idx in topoZ) {
    par <- which(Azz[,idx]==1)
    if (!length(par)) {
      Zmat[,idx] <- rnorm(n)
    } else {
      Pa <- prep(Zmat[,par,drop=FALSE], lin_pr)
      Zmat[,idx] <- as.numeric(Pa %*% rnorm(ncol(Pa), sd=1)) + rnorm(n, sd=1)
    }
  }
  
  Xmat       <- matrix(0, n, d_x, dimnames=list(NULL, paste0("x",seq_len(d_x))))
  adj_xx_out <- matrix(0, d_x, d_x,
                       dimnames=list(paste0("x",seq_len(d_x)),
                                     paste0("x",seq_len(d_x))))
  diag(adj_xx_out) <- NA_real_
  
  gX    <- igraph::graph_from_adjacency_matrix(Axx, mode="directed")
  topoX <- as.integer(igraph::topo_sort(gX, mode="out"))
  
  for (j_idx in seq_along(topoX)) {
    j   <- topoX[j_idx]; cif <- d_z + j
    zp  <- which(A_full[1:d_z, cif]==1)
    sig_z <- rep(0, n)
    if (length(zp)) {
      PaZ   <- prep(Zmat[,zp,drop=FALSE], lin_pr)
      sig_z <- as.numeric(PaZ %*% sample(c(1,-1), ncol(PaZ), replace=TRUE))
    }
    xp <- intersect(which(A_full[(d_z+1):(d_z+d_x), cif]==1),
                    if (j_idx>1) topoX[1:(j_idx-1)] else integer(0))
    sig_x <- rep(0, n)
    if (length(xp)) {
      PaX          <- prep(Xmat[,xp,drop=FALSE], lin_pr)
      adj_xx_out[j,xp] <- 1
      csd          <- matrixStats::colSds(PaX); csd[csd==0] <- 1
      beta_x       <- (x_effect/csd) * sample(c(1,-1), length(csd), replace=TRUE)
      sig_x        <- as.numeric(PaX %*% beta_x)
    }
    sig_total <- sig_z + sig_x
    if (var(sig_total) < 1e-6) sig_total <- rnorm(n)
    Xmat[,j] <- sig_total + sim_noise(sig_total, r2)
  }
  
  dat_dt <- data.table(Zmat, Xmat)
  colnames(dat_dt) <- c(paste0("z",1:d_z), paste0("x",1:d_x))
  list(dat     = dat_dt,
       adj_full = A_full,
       adj_xx   = adj_xx_out,
       params   = list(n=n, d_z=d_z, d_x=d_x, r2=r2, lin_pr=lin_pr,
                       sp=sp, p_cross=p_cross, x_effect=x_effect))
}


# ======================================================================
# 3. Graph utilities
# ======================================================================

clean_matrix_adj <- function(mat) {
  m <- as.matrix(mat); m[is.na(m)] <- 0; m[m!=0] <- 1; m
}

topo_sort_adj <- function(adj) {
  adj <- clean_matrix_adj(adj); n <- nrow(adj)
  if (n==0) return(integer(0))
  A <- (adj!=0)*1; indeg <- colSums(A); indeg[is.na(indeg)] <- 0
  ord <- integer(0); zeros <- sort(which(indeg==0))
  while (length(zeros)>0) {
    v <- zeros[1]; zeros <- zeros[-1]; ord <- c(ord, v)
    nbrs <- which(A[v,]!=0)
    for (u in nbrs) {
      indeg[u] <- indeg[u]-1
      if (indeg[u]==0) zeros <- sort(unique(c(zeros,u)))
    }
    A[v, nbrs] <- 0
  }
  if (length(ord)!=n) return(NULL); ord
}

is_dag_adj <- function(adj) !is.null(topo_sort_adj(adj))

sort_ancestral_matrix <- function(m) {
  if (is.null(m) || nrow(m)==0) return(m)
  sup <- m; sup[is.na(sup)] <- 0
  sup[!(sup %in% c(0.5,1))] <- 0; sup[sup!=0] <- 1
  ord <- topo_sort_adj(sup)
  if (is.null(ord)) { warning("Not a DAG"); return(m) }
  s <- m[ord, ord, drop=FALSE]; diag(s) <- NA; s[lower.tri(s)] <- 0; s
}

build_constraint_adj <- function(m) {
  n <- nrow(m); A <- matrix(0, n, n)
  rownames(A) <- rownames(m); colnames(A) <- colnames(m)
  for (j in seq_len(n)) for (i in seq_len(n)) {
    if (i==j || is.na(m[j,i])) next
    if (m[j,i]==1 || m[j,i]==0.5) A[i,j] <- 1
  }; A
}

permute_to_lower <- function(m) {
  ord <- topo_sort_adj(build_constraint_adj(m))
  if (is.null(ord)) stop("Ancestry matrix contains cycles")
  m[as.numeric(ord), as.numeric(ord)]
}

try_edge_update <- function(adj, i, j, val_ij, val_ji) {
  t <- adj; t[i,j] <- val_ij; t[j,i] <- val_ji
  b <- ((t==1 | t==0.5)*1); b[is.na(b)] <- 0
  if (is_dag_adj(b)) return(list(adj=t, success=TRUE))
  list(adj=adj, success=FALSE)
}

# Stub required by cbl_fn if cbl_fixed.R is also sourced.
normalize_to_sorted_ancestral <- function(m) m


# ======================================================================
# 4. CI test helper
#
# F-test for H0: target _|_ query | conditioning_set
# Returns p-value. Small p => reject independence => dependent.
# ======================================================================

ci_test_pval <- function(target, query, cond_set, data, min_n=10) {
  vars <- unique(c(target, query, cond_set))
  miss <- setdiff(vars, colnames(data))
  if (length(miss)>0) return(NA)
  if (nrow(data)<min_n) return(NA)
  
  sub  <- data[, vars, drop=FALSE]
  good <- sapply(sub, function(col) {
    v <- na.omit(col); length(unique(v))>1 && sd(v)>0
  })
  sub <- sub[, as.logical(good), drop=FALSE]
  if (!(target %in% colnames(sub)) || !(query %in% colnames(sub))) return(NA)
  
  cs   <- intersect(cond_set, colnames(sub))
  rhs0 <- if (length(cs)>0) paste(cs, collapse="+") else "1"
  rhs1 <- if (length(cs)>0) paste(c(cs, query), collapse="+") else query
  
  f0   <- tryCatch(lm(as.formula(paste(target,"~",rhs0)), data=sub), error=function(e) NULL)
  f1   <- tryCatch(lm(as.formula(paste(target,"~",rhs1)), data=sub), error=function(e) NULL)
  if (is.null(f0) || is.null(f1)) return(NA)
  
  res  <- tryCatch(anova(f0, f1), error=function(e) NULL)
  if (is.null(res)) return(NA)
  res$`Pr(>F)`[2]
}


# ======================================================================
# 5. ASCEND main function
# ======================================================================

ascend_fn <- function(sim_obj,
                      maxiter        = 10,
                      alpha          = 0.05,
                      alpha_mb_start = 0.20,
                      alpha_mb_floor = 0.05,
                      alpha_decay    = 0.70,
                      fdr_correction = TRUE,
                      min_votes      = 1) {
  
  dat   <- sim_obj$dat
  d_x   <- length(grep("^x", colnames(dat)))
  xlabs <- paste0('x', seq_len(d_x))
  
  scale_robust <- function(x) {
    if (all(is.na(x))) return(x)
    m <- mad(x, na.rm=TRUE)
    if (is.na(m) || m==0) {
      s <- sd(x, na.rm=TRUE)
      if (is.na(s) || s==0) { x <- x+rnorm(length(x),0,1e-6); s <- sd(x,na.rm=TRUE) }
      return((x - mean(x,na.rm=TRUE)) / s)
    }
    (x - median(x,na.rm=TRUE)) / m
  }
  
  dataAll <- as.data.frame(dat)
  for (col in colnames(dataAll)) dataAll[[col]] <- scale_robust(dataAll[[col]])
  dataAll[is.infinite(as.matrix(dataAll))] <- 0
  dataAll[is.na(dataAll)] <- 0
  
  z_cols <- grep("^z", colnames(dataAll), value=TRUE)
  
  min_z_keep <- max(3, ceiling(length(z_cols) * 0.10))
  
  prescreen_z <- function(node_i, liberal_thresh=0.30) {
    if (sd(dataAll[[node_i]], na.rm=TRUE)==0) return(character(0))
    pv <- sapply(z_cols, function(zc) {
      ct <- tryCatch(cor.test(dataAll[[node_i]], dataAll[[zc]]), error=function(e) NULL)
      if (is.null(ct)) 1 else ct$p.value
    })
    adj_pv <- p.adjust(pv, method="BH")
    keep   <- which(adj_pv < liberal_thresh)
    if (length(keep) < min_z_keep)
      keep <- order(pv)[1:min(min_z_keep, length(pv))]
    z_cols[keep]
  }
  
  learn_mb_safe <- function(node_i, candidates, alpha_mb) {
    if (length(candidates)==0 || sd(dataAll[[node_i]],na.rm=TRUE)==0)
      return(character(0))
    sub  <- dataAll[, c(candidates, node_i), drop=FALSE]
    good <- sapply(sub, function(col) {
      v <- na.omit(col); length(unique(v))>1 && sd(v)>0
    })
    sub <- sub[, as.logical(good), drop=FALSE]
    if (ncol(sub)<2 || !(node_i %in% colnames(sub))) return(character(0))
    tryCatch(
      learn.mb(sub, node=node_i, method="iamb", test="zf", alpha=alpha_mb),
      error=function(e) character(0)
    )
  }
  
  alpha_mb_current <- alpha_mb_start
  mb_list <- vector("list", d_x); names(mb_list) <- xlabs
  
  for (i in seq_len(d_x)) {
    node_i <- xlabs[i]
    z_cands <- prescreen_z(node_i)
    mb_list[[node_i]] <- learn_mb_safe(node_i, z_cands, alpha_mb_current)
  }
  
  M <- matrix(NA, d_x, d_x); diag(M) <- NA
  rownames(M) <- colnames(M) <- xlabs
  
  converged <- FALSE; iter <- 0
  
  while (!converged && iter <= maxiter) {
    iter <- iter+1; converged <- TRUE
    alpha_mb_current <- max(alpha_mb_floor,
                            alpha_mb_start * alpha_decay^(iter-1))
    cat(sprintf("Iteration %d | alpha_mb=%.3f | NA pairs remaining: %d\n",
                iter, alpha_mb_current, sum(is.na(M[upper.tri(M)]))))
    
    pval_info <- list(); pidx <- 1
    for (i in 2:d_x) {
      for (j in 1:(i-1)) {
        if (!is.na(M[i,j])) next
        xi <- xlabs[i]; xj <- xlabs[j]
        S  <- setdiff(union(mb_list[[xi]], mb_list[[xj]]), c(xi,xj))
        S  <- intersect(S, colnames(dataAll))
        pv <- ci_test_pval(xj, xi, S, dataAll)
        if (!is.na(pv)) {
          pval_info[[pidx]] <- list(i=i,j=j,xi=xi,xj=xj,pval=pv,S=S)
          pidx <- pidx+1
        } else {
          converged <- FALSE
        }
      }
    }
    
    if (length(pval_info)>0) {
      pvals     <- sapply(pval_info, function(x) x$pval)
      adj_pvals <- if (fdr_correction) p.adjust(pvals, method="BH") else pvals
      
      for (idx in seq_along(pval_info)) {
        info     <- pval_info[[idx]]
        i<-info$i; j<-info$j; xi<-info$xi; xj<-info$xj; S<-info$S
        adj_pval <- adj_pvals[idx]
        
        if (adj_pval > alpha) {
          M[i,j] <- 0; M[j,i] <- 0; converged <- FALSE; next
        }
        if (length(S)==0) next
        
        votes    <- c(r1=0L, r1_rev=0L, r2=0L, r2_rev=0L)
        evidence <- c(r1=0.0, r1_rev=0.0, r2=0.0, r2_rev=0.0)
        
        for (W in S) {
          Sw <- setdiff(S, W)
          p_wj_no_i <- ci_test_pval(W, xj, Sw,        dataAll)
          p_wj_wi_i <- ci_test_pval(W, xj, c(Sw,xi),  dataAll)
          p_wi_no_j <- ci_test_pval(W, xi, Sw,        dataAll)
          p_wi_wi_j <- ci_test_pval(W, xi, c(Sw,xj),  dataAll)
          
          r1f  <- !is.na(p_wj_no_i)&&!is.na(p_wj_wi_i) && p_wj_no_i<=alpha && p_wj_wi_i>alpha
          r1rf <- !is.na(p_wi_no_j)&&!is.na(p_wi_wi_j) && p_wi_no_j<=alpha && p_wi_wi_j>alpha
          r2f  <- !is.na(p_wi_no_j)&&!is.na(p_wi_wi_j) && p_wi_no_j>alpha  && p_wi_wi_j<=alpha
          r2rf <- !is.na(p_wj_no_i)&&!is.na(p_wj_wi_i) && p_wj_no_i>alpha  && p_wj_wi_i<=alpha
          
          votes["r1"]     <- votes["r1"]     + r1f
          votes["r1_rev"] <- votes["r1_rev"] + r1rf
          votes["r2"]     <- votes["r2"]     + r2f
          votes["r2_rev"] <- votes["r2_rev"] + r2rf
          
          if (r1f  && !is.na(p_wj_no_i)) evidence["r1"]     <- evidence["r1"]     - log(p_wj_no_i  + 1e-300)
          if (r1rf && !is.na(p_wi_no_j)) evidence["r1_rev"] <- evidence["r1_rev"] - log(p_wi_no_j  + 1e-300)
          if (r2f  && !is.na(p_wi_wi_j)) evidence["r2"]     <- evidence["r2"]     - log(p_wi_wi_j  + 1e-300)
          if (r2rf && !is.na(p_wj_wi_i)) evidence["r2_rev"] <- evidence["r2_rev"] - log(p_wj_wi_i  + 1e-300)
        }
        
        if (votes["r1"]>0 && votes["r1_rev"]>0) {
          if      (evidence["r1"] > evidence["r1_rev"]) votes["r1_rev"] <- 0L
          else if (evidence["r1_rev"] > evidence["r1"]) votes["r1"]     <- 0L
          else { votes["r1"] <- 0L; votes["r1_rev"] <- 0L }
        }
        if (votes["r2"]>0 && votes["r2_rev"]>0) {
          if      (evidence["r2"] > evidence["r2_rev"]) votes["r2_rev"] <- 0L
          else if (evidence["r2_rev"] > evidence["r2"]) votes["r2"]     <- 0L
          else { votes["r2"] <- 0L; votes["r2_rev"] <- 0L }
        }
        if (votes["r1"]>0)     votes["r2"]     <- 0L
        if (votes["r1_rev"]>0) votes["r2_rev"] <- 0L
        
        best <- max(votes)
        if (best < min_votes) next
        
        result <- NULL
        if (votes["r1"]>=min_votes && votes["r1"]>=votes["r1_rev"])
          result <- try_edge_update(M,i,j,1,0)
        if ((is.null(result)||!result$success) && votes["r1_rev"]>=min_votes)
          result <- try_edge_update(M,i,j,0,1)
        if ((is.null(result)||!result$success) && votes["r2"]>=min_votes && votes["r2"]>=votes["r2_rev"])
          result <- try_edge_update(M,i,j,0.5,0)
        if ((is.null(result)||!result$success) && votes["r2_rev"]>=min_votes)
          result <- try_edge_update(M,i,j,0,0.5)
        
        if (!is.null(result) && result$success) {
          M <- result$adj; converged <- FALSE
        }
      }
    }
    
    # Transitive closure
    done <- FALSE
    while (!done) {
      done <- TRUE
      for (k in seq_len(d_x)) {
        ancs  <- which(M[,k]==1)
        descs <- which(M[k,]==1)
        if (!length(ancs) || !length(descs)) next
        for (anc in ancs) for (desc in descs) {
          if (anc==desc) next
          if (is.na(M[anc,desc]) || M[anc,desc]!=1) {
            res <- try_edge_update(M, anc, desc, 1, 0)
            if (res$success) { M <- res$adj; done <- FALSE; converged <- FALSE }
          }
        }
      }
    }
    
    # Symmetry closure: Xi <= Xj AND Xj <= Xi => Xi ~ Xj
    for (i in seq_len(d_x)) for (j in seq_len(d_x)) {
      if (i==j) next
      if (!is.na(M[i,j]) && !is.na(M[j,i]) && M[i,j]==0.5 && M[j,i]==0.5) {
        M[i,j] <- 0; M[j,i] <- 0; converged <- FALSE
      }
    }
    
    bin <- clean_matrix_adj(M); bin[is.na(bin)] <- 0
    if (!is_dag_adj(bin)) stop(sprintf("Cycle detected at iteration %d", iter))
    
    # Update non-descendant sets and recompute MB
    for (i in seq_len(d_x)) {
      node_i       <- xlabs[i]
      z_parents_Xi <- intersect(mb_list[[node_i]], z_cols)
      if (length(z_parents_Xi)==0) z_parents_Xi <- prescreen_z(node_i)
      confirmed_x_anc <- xlabs[which(M[,i]==1 | M[,i]==0.5)]
      T_Xi <- setdiff(unique(c(z_parents_Xi, confirmed_x_anc)), node_i)
      T_Xi <- intersect(T_Xi, colnames(dataAll))
      new_mb <- learn_mb_safe(node_i, T_Xi, alpha_mb_current)
      if (!setequal(new_mb, mb_list[[node_i]])) converged <- FALSE
      mb_list[[node_i]] <- new_mb
    }
  }
  
  # Final R3 pass on remaining NA pairs
  for (i in 2:d_x) {
    for (j in 1:(i-1)) {
      if (!is.na(M[i,j])) next
      xi <- xlabs[i]; xj <- xlabs[j]
      S  <- setdiff(union(mb_list[[xi]], mb_list[[xj]]), c(xi,xj))
      S  <- intersect(S, colnames(dataAll))
      pv <- ci_test_pval(xj, xi, S, dataAll)
      if (!is.na(pv) && pv > alpha) { M[i,j] <- 0; M[j,i] <- 0 }
    }
  }
  
  # Final symmetry closure
  for (i in seq_len(d_x)) for (j in seq_len(d_x)) {
    if (i==j) next
    if (!is.na(M[i,j]) && !is.na(M[j,i]) && M[i,j]==0.5 && M[j,i]==0.5) {
      M[i,j] <- 0; M[j,i] <- 0
    }
  }
  
  if (!all(is.na(M))) {
    tryCatch({
      M <- permute_to_lower(M)
      M <- sort_ancestral_matrix(M)
    }, error=function(e) warning("Final normalisation failed: ", e$message))
  }
  
  M
}


# ======================================================================
# 6. Evaluation utilities
# ======================================================================

get_true_ancestral_matrix <- function(adj_xx) {
  d <- nrow(adj_xx); nms <- rownames(adj_xx)
  A <- t(adj_xx); diag(A) <- 0; A[is.na(A)] <- 0
  reach <- (A>0)*1; Ap <- A
  for (k in seq_len(d-1)) { Ap <- ((Ap%*%A)>0)*1; reach <- ((reach+Ap)>0)*1 }
  rownames(reach) <- colnames(reach) <- nms; diag(reach) <- NA
  ord <- topo_sort_adj(reach); if (!is.null(ord)) reach <- reach[ord,ord]
  reach
}

evaluate_ancestral <- function(estimated, truth, verbose=TRUE) {
  nms       <- rownames(truth)
  estimated <- estimated[nms, nms]
  d         <- nrow(truth)
  tp<-0; fp<-0; fn<-0; tn<-0; unresolved_tp<-0; unresolved_tn<-0
  
  for (i in seq_len(d-1)) for (j in (i+1):d) {
    t_ij <- truth[i,j]; e_ij <- estimated[i,j]
    if (is.na(t_ij)) next
    true_pos <- (t_ij==1)
    pred_pos <- !is.na(e_ij) && (e_ij==1 || e_ij==0.5)
    resolved <- !is.na(e_ij)
    if (!resolved) {
      if (true_pos) unresolved_tp <- unresolved_tp+1
      else          unresolved_tn <- unresolved_tn+1
      next
    }
    if ( true_pos &&  pred_pos) tp <- tp+1
    if (!true_pos &&  pred_pos) fp <- fp+1
    if ( true_pos && !pred_pos) fn <- fn+1
    if (!true_pos && !pred_pos) tn <- tn+1
  }
  
  total_pairs    <- tp+fp+fn+tn+unresolved_tp+unresolved_tn
  resolved_pairs <- tp+fp+fn+tn
  coverage       <- if (total_pairs>0)    resolved_pairs/total_pairs else NA
  acc            <- if (resolved_pairs>0) (tp+tn)/resolved_pairs     else NA
  pr             <- if ((tp+fp)>0)        tp/(tp+fp)                 else NA
  re_res         <- if ((tp+fn)>0)        tp/(tp+fn)                 else NA
  all_true_pos   <- tp+fn+unresolved_tp
  re_all         <- if (all_true_pos>0)   tp/all_true_pos            else NA
  f1             <- if (!is.na(pr)&&!is.na(re_res)&&(pr+re_res)>0)
    2*pr*re_res/(pr+re_res) else NA
  
  if (verbose) {
    cat("=== ASCEND Evaluation ===\n")
    cat(sprintf("  Resolved pairs : %d / %d  (Coverage=%.3f)\n",
                resolved_pairs, total_pairs, coverage))
    cat(sprintf("  Unresolved (NA): %d  [%d true ancestors, %d true non-ancestors]\n",
                unresolved_tp+unresolved_tn, unresolved_tp, unresolved_tn))
    cat(sprintf("  TP=%-4d FP=%-4d FN=%-4d TN=%-4d\n", tp, fp, fn, tn))
    cat(sprintf("  Precision         = %.3f\n", pr))
    cat(sprintf("  Recall (resolved) = %.3f\n", re_res))
    cat(sprintf("  Recall (overall)  = %.3f\n", re_all))
    cat(sprintf("  F1 (resolved)     = %.3f\n", f1))
    cat(sprintf("  Accuracy          = %.3f\n", acc))
    cat("=========================\n")
  }
  invisible(list(tp=tp, fp=fp, fn=fn, tn=tn,
                 unresolved_tp=unresolved_tp, unresolved_tn=unresolved_tn,
                 coverage=coverage, precision=pr,
                 recall_resolved=re_res, recall_overall=re_all,
                 f1=f1, accuracy=acc))
}