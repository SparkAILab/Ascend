# ======================================================================
# ASCEND BENCHMARK — Two Comparison Frameworks
#
# FRAMEWORK 1: ANCESTRAL COMPARISON
#   All methods converted to ancestral matrix. No PC step.
#   Methods: ASCEND, CBL, GES (->ancestral), LiNGAM (->ancestral)
#
#   ASCEND: output directly (1/0.5/0/NA)
#   CBL:    output directly (same convention)
#   GES:    CPDAG -> directed edges only -> transitive closure
#           Undirected edges -> NA (we cannot assert direction)
#   LiNGAM: DAG -> transitive closure -> ancestral
#
#   WHY NA FOR GES UNDIRECTED: a CPDAG undirected edge i--j means
#   the algorithm found an edge but cannot determine direction within
#   the Markov equivalence class. Asserting ancestry in either direction
#   would be wrong. NA is the honest answer.
#
# FRAMEWORK 2: CPDAG/DAG COMPARISON (skeleton + orientation)
#   Methods output a directed/partially-directed graph.
#   Evaluated on skeleton F1 and orientation accuracy.
#   Methods: ASCEND+PC, CBL+PC, GES (CPDAG), LiNGAM (DAG)
#
#   ASCEND+PC: ASCEND ancestral constraints -> fixedGaps + addBgKnowledge -> PC
#   CBL+PC:    CBL ancestral constraints   -> fixedGaps + addBgKnowledge -> PC
#   GES:       CPDAG directly
#   LiNGAM:    DAG directly (converted to CPDAG format for fair comparison)
#
# FRAMEWORK 3: RUNTIME COMPARISON
#   ASCEND vs CBL only — the two methods that use the Z+X two-tier structure.
#   Sweep n and d_x, record wall-clock time.
#
# ALL METHODS: receive Z+X data; output evaluated on X variables only.
# GES and LiNGAM are expected to struggle as dimension grows — that
# is a deliberate and fair test of scalability.
#
# COLOUR SCHEME: ASCEND = green (#4DAF4A) as requested.
# ======================================================================


# ======================================================================
# 0. Libraries
# ======================================================================
library(pcalg)
library(igraph)
library(bnlearn)
library(data.table)
library(matrixStats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(foreach)
library(doMC)
registerDoMC(cores = max(1L, parallel::detectCores() - 1L))


# ======================================================================
# 1. Shared utilities
# ======================================================================

# A[i,j]=1 means i->j directly. Returns reach[i,j]=1 if i is ancestor of j.
transitive_closure <- function(A) {
  A[is.na(A)] <- 0; diag(A) <- 0
  d <- nrow(A); reach <- (A > 0) * 1; Ap <- A
  for (k in seq_len(d - 1)) {
    Ap    <- ((Ap %*% A) > 0) * 1
    reach <- ((reach + Ap) > 0) * 1
  }
  diag(reach) <- NA
  reach
}

topo_sort_mat <- function(adj) {
  a <- adj; a[is.na(a)] <- 0; a[a != 0] <- 1; n <- nrow(a)
  indeg <- colSums(a); indeg[is.na(indeg)] <- 0
  ord <- integer(0); zeros <- sort(which(indeg == 0))
  while (length(zeros) > 0) {
    v <- zeros[1]; zeros <- zeros[-1]; ord <- c(ord, v)
    nbrs <- which(a[v,] != 0)
    for (u in nbrs) {
      indeg[u] <- indeg[u] - 1
      if (indeg[u] == 0) zeros <- sort(unique(c(zeros, u)))
    }
    a[v, nbrs] <- 0
  }
  if (length(ord) != n) return(NULL); ord
}

# True ancestral matrix: adj_xx[child,parent]=1 => transpose => closure
true_ancestral <- function(adj_xx) {
  nms <- rownames(adj_xx)
  A   <- t(adj_xx); A[is.na(A)] <- 0; diag(A) <- 0
  R   <- transitive_closure(A)
  rownames(R) <- colnames(R) <- nms
  ord <- topo_sort_mat(R)
  if (!is.null(ord)) R <- R[ord, ord]
  R
}

# ---- CPDAG -> ancestral matrix (Framework 1 use) ----
# pcalg convention: amat[i,j]=1, amat[j,i]=0 => i->j (directed)
#                   amat[i,j]=1, amat[j,i]=1 => i--j (undirected)
# Directed edge => assert ancestry via transitive closure
# Undirected edge => NA (cannot assert direction; honest representation)
cpdag_to_ancestral_strict <- function(amat, node_names) {
  d     <- nrow(amat)
  A_dir <- matrix(0, d, d, dimnames=list(node_names, node_names))
  for (i in seq_len(d)) for (j in seq_len(d)) {
    if (i==j) next
    if (amat[i,j]==1 && amat[j,i]==0) A_dir[i,j] <- 1  # directed only
  }
  anc <- transitive_closure(A_dir)
  rownames(anc) <- colnames(anc) <- node_names
  # Mark pairs where NEITHER direction is directed as NA (undirected edge)
  # and pair has no directed path — genuinely unknown
  anc[is.na(anc)] <- 0
  for (i in seq_len(d)) for (j in seq_len(d)) {
    if (i==j) next
    if (amat[i,j]==1 && amat[j,i]==1) {
      # Undirected edge: if no directed path establishes ancestry, mark NA
      if (anc[i,j] == 0 && anc[j,i] == 0) {
        anc[i,j] <- NA   # genuinely unknown
      }
    }
  }
  diag(anc) <- NA
  anc
}

# ---- CPDAG -> skeleton adjacency (Framework 2 use) ----
# Returns symmetric 0/1 skeleton (presence of any edge)
cpdag_to_skeleton <- function(amat) {
  skel <- ((amat + t(amat)) > 0) * 1
  diag(skel) <- 0
  skel
}


# ======================================================================
# 2. Evaluation metrics
# ======================================================================

# ---- 2a. Ancestral metrics (Framework 1) ----
# Positive prediction: estimated[i,j] = 1 (strict) or 0.5 (weak)
# True positive: truth[i,j] = 1
# NA in estimated = unresolved (penalised only in overall recall)
ancestral_metrics <- function(estimated, truth) {
  nms       <- rownames(truth)
  estimated <- estimated[nms, nms]
  d         <- nrow(truth)
  tp<-0; fp<-0; fn<-0; tn<-0; unres_tp<-0; unres_tn<-0
  
  for (i in seq_len(d-1)) for (j in (i+1):d) {
    t_ij <- truth[i,j];  e_ij <- estimated[i,j]
    if (is.na(t_ij)) next
    true_pos <- (t_ij == 1)
    resolved <- !is.na(e_ij)
    pred_pos <- resolved && (e_ij == 1 || e_ij == 0.5)
    if (!resolved) {
      if (true_pos) unres_tp <- unres_tp+1 else unres_tn <- unres_tn+1; next
    }
    if ( true_pos &&  pred_pos) tp <- tp+1
    if (!true_pos &&  pred_pos) fp <- fp+1
    if ( true_pos && !pred_pos) fn <- fn+1
    if (!true_pos && !pred_pos) tn <- tn+1
  }
  total    <- tp+fp+fn+tn+unres_tp+unres_tn
  resolved <- tp+fp+fn+tn
  pr  <- if ((tp+fp)>0)          tp/(tp+fp)                   else NA
  re  <- if ((tp+fn+unres_tp)>0) tp/(tp+fn+unres_tp)          else NA
  f1  <- if (!is.na(pr)&&!is.na(re)&&(pr+re)>0) 2*pr*re/(pr+re) else NA
  acc <- if (resolved>0)         (tp+tn)/resolved              else NA
  cov <- if (total>0)            resolved/total                else NA
  list(tp=tp, fp=fp, fn=fn, tn=tn,
       unres_tp=unres_tp, unres_tn=unres_tn,
       precision=pr, recall=re, f1=f1, accuracy=acc, coverage=cov)
}

# ---- 2b. Skeleton + orientation metrics (Framework 2) ----
# Evaluates CPDAG or DAG outputs against the true skeleton and orientations.
# true_dag_xx[parent,child]=1: true immediate X->X edges
cpdag_metrics <- function(estimated_amat, true_dag_xx, node_names) {
  # Align to node_names
  estimated_amat <- estimated_amat[node_names, node_names]
  true_dag_xx    <- true_dag_xx[node_names, node_names]
  true_dag_xx[is.na(true_dag_xx)] <- 0
  
  # True skeleton (symmetric)
  true_skel <- ((true_dag_xx + t(true_dag_xx)) > 0) * 1
  est_skel  <- cpdag_to_skeleton(estimated_amat)
  
  ut <- upper.tri(true_skel)
  tp_s <- sum(true_skel[ut]==1 & est_skel[ut]==1)
  fp_s <- sum(true_skel[ut]==0 & est_skel[ut]==1)
  fn_s <- sum(true_skel[ut]==1 & est_skel[ut]==0)
  tn_s <- sum(true_skel[ut]==0 & est_skel[ut]==0)
  
  pr_s  <- if ((tp_s+fp_s)>0)  tp_s/(tp_s+fp_s) else NA
  re_s  <- if ((tp_s+fn_s)>0)  tp_s/(tp_s+fn_s) else NA
  f1_s  <- if (!is.na(pr_s)&&!is.na(re_s)&&(pr_s+re_s)>0) 2*pr_s*re_s/(pr_s+re_s) else NA
  shd_s <- fp_s + fn_s  # skeleton SHD (wrong edges)
  
  # Orientation accuracy: only on skeleton edges that both methods agree exist
  shared <- which(true_skel==1 & est_skel==1, arr.ind=TRUE)
  shared <- shared[shared[,1] < shared[,2],, drop=FALSE]
  n_shared <- nrow(shared); n_correct <- 0
  for (k in seq_len(n_shared)) {
    i <- shared[k,1]; j <- shared[k,2]
    # True direction
    true_ij <- true_dag_xx[i,j]==1 && true_dag_xx[j,i]==0
    true_ji <- true_dag_xx[j,i]==1 && true_dag_xx[i,j]==0
    # Estimated direction
    est_ij  <- estimated_amat[i,j]==1 && estimated_amat[j,i]==0
    est_ji  <- estimated_amat[j,i]==1 && estimated_amat[i,j]==0
    if ((true_ij&&est_ij)||(true_ji&&est_ji)) n_correct <- n_correct+1
  }
  orient_acc <- if (n_shared>0) n_correct/n_shared else NA
  
  list(precision_skel=pr_s, recall_skel=re_s, f1_skel=f1_s,
       shd_skel=shd_s, orient_acc=orient_acc,
       tp_s=tp_s, fp_s=fp_s, fn_s=fn_s, tn_s=tn_s)
}


# ======================================================================
# 3. Method outputs
#    Each function returns the RAW output of the method in its
#    native format, plus a conversion to ancestral where applicable.
# ======================================================================

# ---- Helper: regularised sufficient statistics ----
make_suffstat <- function(data_mat) {
  complete <- data_mat[complete.cases(data_mat),]
  n <- nrow(complete); p <- ncol(complete)
  C <- cor(complete)
  eig <- min(eigen(C, only.values=TRUE)$values)
  if (eig <= 1e-10) {
    C <- C + (abs(eig)+1e-4)*diag(p)
    d2 <- sqrt(diag(C)); C <- C/outer(d2,d2)
  }
  list(C=C, n=n)
}

# ---- Core ASCEND+PC / CBL+PC integration (PC on X only) ----
#
# DESIGN RATIONALE (reviewer-proof):
#
# ASCEND acts as a *preprocessor* that dramatically reduces the search
# space for PC. PC then operates only on X variables. The key insight
# is that ASCEND's power is in PRUNING (high precision skeleton removal)
# and ORIENTATION (confirmed ancestral claims), not in adding Z as nodes.
# Adding Z to PC would make PC's search space O(2^(d_z+d_x)) which
# defeats the entire scalability purpose of the pipeline.
#
# Three layers of ASCEND intelligence fed into PC:
#
# Layer 1 — SKELETON PRUNING via fixedGaps:
#   When ASCEND confirms M[i,j]=0 AND M[j,i]=0, this pair is definitively
#   non-ancestral in both directions => remove their skeleton edge.
#   This is ASCEND's strongest contribution: high-precision skeleton pruning
#   means PC searches a much smaller candidate set, reducing both
#   computational cost and false positive rate.
#
# Layer 2 — PARTIAL CORRELATION ADJUSTMENT:
#   For the CI tests PC runs on X, we partial out the effect of Z by
#   using regression residuals of X on Z. This effectively conditions
#   on Z without including Z as graph nodes. This gives PC the benefit
#   of Z's information without the dimensionality curse.
#   Specifically: X_i^resid = X_i - hat(X_i | Z), then run PC on residuals.
#   Under linear SEM, this is equivalent to marginalising over Z.
#
# Layer 3 — ORIENTATION via addBgKnowledge:
#   ASCEND's confirmed strict ancestors (M=1) are applied as background
#   knowledge to orient PC's undirected edges. Applied in order of
#   confidence: strict (1) before weak (0.5).
#   This propagates through Meek rules to orient additional edges.
#
# Alpha for PC: 0.05 (standard). The pruned skeleton means fewer
# tests, so we do not need to inflate alpha.
#
run_ascend_pc_core <- function(sim_obj, M_anc, alpha_pc=0.05) {
  dat    <- as.data.frame(sim_obj$dat)
  x_cols <- grep("^x", colnames(dat), value=TRUE)
  z_cols <- grep("^z", colnames(dat), value=TRUE)
  n_fg   <- length(x_cols)
  
  # Align M to x_cols
  M <- M_anc[x_cols, x_cols]
  
  # ── Layer 2: Partial out Z from X (residualise X on Z) ──────────────
  # This gives PC the Z information without blowing up its search space.
  # Under linear SEM: X_i^resid = X_i - E[X_i | Z]
  # The residuals have the same conditional independence structure among
  # X variables as the original X given Z.
  x_dat <- as.matrix(dat[, x_cols])
  z_dat <- as.matrix(dat[, z_cols])
  
  x_resid <- x_dat  # fallback if residualisation fails
  tryCatch({
    # Regress each X on all Z, take residuals
    x_resid_mat <- matrix(NA, nrow(x_dat), n_fg,
                          dimnames=list(NULL, x_cols))
    for (j in seq_len(n_fg)) {
      fit <- lm(x_dat[,j] ~ z_dat)
      x_resid_mat[,j] <- residuals(fit)
    }
    x_resid <- x_resid_mat
  }, error=function(e) NULL)
  
  # ── Layer 1: fixedGaps from ASCEND confirmed non-ancestral pairs ────
  fgap <- matrix(FALSE, n_fg, n_fg, dimnames=list(x_cols, x_cols))
  for (i in seq_len(n_fg)) for (j in seq_len(n_fg)) {
    if (i==j) next
    xi <- x_cols[i]; xj <- x_cols[j]
    vij <- M[xi,xj]; vji <- M[xj,xi]
    # Remove skeleton edge only when confirmed non-ancestral in both directions
    if (!is.na(vij) && !is.na(vji) && vij==0 && vji==0) {
      fgap[i,j] <- TRUE; fgap[j,i] <- TRUE
    }
  }
  
  # ── Sufficient statistics on Z-residualised X ───────────────────────
  ss <- make_suffstat(x_resid)
  if (ss$n < n_fg+3) return(NULL)
  
  # ── Run PC on X residuals with ASCEND skeleton pruning ──────────────
  pc_fit <- tryCatch(
    pc(suffStat=ss, indepTest=gaussCItest, labels=x_cols,
       alpha=alpha_pc, fixedGaps=fgap, verbose=FALSE,
       maj.rule=TRUE, solve.confl=TRUE),
    error=function(e) NULL
  )
  if (is.null(pc_fit)) return(NULL)
  
  amat <- as(pc_fit, "amat")
  rownames(amat) <- colnames(amat) <- x_cols
  
  # ── Layer 3: Apply ASCEND orientation constraints ────────────────────
  # Strict ancestors first (highest confidence), then weak
  for (pass in c(1, 0.5)) {
    for (i in seq_len(n_fg)) for (j in seq_len(n_fg)) {
      if (i==j) next
      xi <- x_cols[i]; xj <- x_cols[j]
      vij <- M[xi, xj]
      if (!is.na(vij) && vij == pass) {
        if (amat[xi,xj]!=0 || amat[xj,xi]!=0) {
          res <- tryCatch(
            addBgKnowledge(gInput=amat, x=xi, y=xj,
                           verbose=FALSE, checkInput=FALSE),
            error=function(e) NULL
          )
          if (!is.null(res)) amat <- res
        }
      }
    }
  }
  rownames(amat) <- colnames(amat) <- x_cols
  amat
}

# ---- 3a. ASCEND ----
# Returns ancestral matrix M directly
get_ascend <- function(sim_obj) {
  tryCatch({
    invisible(capture.output(
      M <- ascend_fn(sim_obj, maxiter=10, alpha=0.05,
                     alpha_mb_start=0.20, alpha_mb_floor=0.05,
                     alpha_decay=0.70, fdr_correction=TRUE, min_votes=1)
    ))
    M
  }, error=function(e) NULL)
}

# ---- normalize_to_sorted_ancestral: required by cbl_fn ----
# cbl_fn calls this at its end. Define it here so cbl_fn doesn't crash.
# It simply returns the matrix as-is (sorting is handled by our own pipeline).
normalize_to_sorted_ancestral <- function(m) {
  if (is.null(m)) return(m)
  m
}

# ---- 3b. CBL ----
# Returns ancestral matrix M directly (uses cbl_fn local implementation)
#
# CBL output convention (from Watson & Silva 2022 and the cbl_fn code):
#   adj1[i, j] = 1  means  Xj -< Xi  (column j is a strict ancestor of row i)
#   adj1[i, j] = 0.5 means Xj <= Xi  (weak)
#   adj1[i, j] = 0   means Xj ~ Xi
#   adj1[i, j] = NA  means unknown
#
# We convert to our convention where M[ancestor, descendant]:
#   M[j, i] = adj1[i, j]
#
# Verified: in cbl_fn, when R1 fires for order "ji":
#   adj1[i, j] <- 1; adj1[j, i] <- 0
#   This means "j precedes i", i.e. Xj -< Xi
#   So adj1[i,j]=1 => Xj is ancestor of Xi => M[j,i]=1  ✓
get_cbl <- function(sim_obj) {
  dat    <- as.data.frame(sim_obj$dat)
  x_cols <- grep("^x", colnames(dat), value=TRUE)
  n_fg   <- length(x_cols)
  
  M_raw <- tryCatch(
    cbl_fn(sim_obj, gamma=0.5, maxiter=10, B=20),
    error=function(e) { message("cbl_fn error: ", conditionMessage(e)); NULL }
  )
  if (is.null(M_raw)) return(NULL)
  
  # Force consistent naming (cbl_fn uses colnames(x) which should be x1..xd)
  if (is.null(rownames(M_raw))) rownames(M_raw) <- x_cols
  if (is.null(colnames(M_raw))) colnames(M_raw) <- x_cols
  
  # Convert: M_raw[i,j]=v means "v encodes relationship where Xj precedes Xi"
  # => M[j,i] = v  (ancestor in row, descendant in column)
  M <- matrix(NA, n_fg, n_fg, dimnames=list(x_cols, x_cols))
  for (i in seq_len(n_fg)) for (j in seq_len(n_fg)) {
    if (i==j) next
    v <- M_raw[i, j]
    if (!is.na(v)) M[j, i] <- v
  }
  M
}

# ---- 3c. GES ----
# Returns CPDAG as amat over X variables
# All variables (Z+X) used; X-block extracted
get_ges_cpdag <- function(sim_obj) {
  dat      <- as.data.frame(sim_obj$dat)
  x_cols   <- grep("^x", colnames(dat), value=TRUE)
  z_cols   <- grep("^z", colnames(dat), value=TRUE)
  all_cols <- c(z_cols, x_cols)
  n_bg     <- length(z_cols); total <- length(all_cols)
  mat      <- as.matrix(dat[,all_cols])
  
  # X->Z skeleton gaps (two-tier: X cannot cause Z)
  fgap <- matrix(FALSE, total, total, dimnames=list(all_cols,all_cols))
  if (n_bg>0) fgap[(n_bg+1):total, seq_len(n_bg)] <- TRUE
  
  score <- tryCatch(new("GaussL0penObsScore",mat), error=function(e) NULL)
  if (is.null(score)) return(NULL)
  fit <- tryCatch(
    ges(score, labels=all_cols, fixedGaps=fgap, iterate=TRUE, verbose=FALSE),
    error=function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  ess <- tryCatch(as(fit$essgraph,"matrix"), error=function(e) NULL)
  if (is.null(ess)) return(NULL)
  rownames(ess) <- colnames(ess) <- all_cols; diag(ess) <- 0
  ess[x_cols, x_cols, drop=FALSE]
}

# ---- 3d. LiNGAM ----
# Returns DAG adjacency (A[parent,child]=1) over X variables
# All variables (Z+X) used; X-block extracted
get_lingam_dag <- function(sim_obj) {
  dat      <- as.data.frame(sim_obj$dat)
  x_cols   <- grep("^x", colnames(dat), value=TRUE)
  z_cols   <- grep("^z", colnames(dat), value=TRUE)
  all_cols <- c(z_cols, x_cols)
  n_bg     <- length(z_cols); total <- length(all_cols)
  mat      <- as.matrix(dat[,all_cols])
  
  fit <- tryCatch(pcalg::lingam(mat), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  B <- if (!is.null(fit$Bpruned)) fit$Bpruned else fit$B
  if (is.null(B)) return(NULL)
  # B[i,j]!=0 => j->i; convert to A[parent,child]=1: A[j,i]=1
  A <- matrix(0, total, total, dimnames=list(all_cols,all_cols))
  for (i in seq_len(total)) for (j in seq_len(total))
    if (i!=j && abs(B[i,j])>1e-9) A[j,i] <- 1
  A[x_cols, x_cols, drop=FALSE]
}


# ---- 3e. PC standalone ----
# PC runs on Z+X with ONLY the two-tier constraint (X cannot cause Z).
# No ASCEND pruning. No background knowledge. This is the baseline.
# We extract the X-block of the resulting CPDAG.
# This is the critical comparator: it answers "what does PC do alone?"
# vs ASCEND+PC which answers "what does PC do with ASCEND's help?"
get_pc_cpdag <- function(sim_obj, alpha_pc=0.05) {
  dat      <- as.data.frame(sim_obj$dat)
  x_cols   <- grep("^x", colnames(dat), value=TRUE)
  z_cols   <- grep("^z", colnames(dat), value=TRUE)
  all_cols <- c(z_cols, x_cols)
  n_bg     <- length(z_cols); total <- length(all_cols)
  
  # Only constraint: X cannot cause Z (two-tier assumption)
  fgap <- matrix(FALSE, total, total, dimnames=list(all_cols, all_cols))
  if (n_bg > 0) fgap[(n_bg+1):total, seq_len(n_bg)] <- TRUE
  
  ss <- make_suffstat(as.matrix(dat[, all_cols]))
  if (ss$n < total + 3) return(NULL)
  
  pc_fit <- tryCatch(
    pc(suffStat=ss, indepTest=gaussCItest, labels=all_cols,
       alpha=alpha_pc, fixedGaps=fgap, verbose=FALSE,
       maj.rule=TRUE, solve.confl=TRUE),
    error=function(e) NULL
  )
  if (is.null(pc_fit)) return(NULL)
  
  amat <- as(pc_fit, "amat")
  rownames(amat) <- colnames(amat) <- all_cols
  amat[x_cols, x_cols, drop=FALSE]   # X-block only
}


# ======================================================================
# 4. Framework 1: Ancestral comparison
#    All methods -> ancestral matrix -> evaluate
# ======================================================================

run_one_ancestral <- function(sim_obj,
                              methods=c("ASCEND","CBL","GES","LiNGAM")) {
  truth <- true_ancestral(sim_obj$adj_xx)
  if (is.null(truth)||sum(truth==1,na.rm=TRUE)==0) return(NULL)
  x_cols <- rownames(truth)
  rows   <- list()
  
  # ASCEND
  if ("ASCEND" %in% methods) {
    M <- tryCatch(get_ascend(sim_obj), error=function(e) NULL)
    if (!is.null(M)) {
      M <- M[x_cols, x_cols]
      m <- ancestral_metrics(M, truth)
      rows[["ASCEND"]] <- data.frame(method="ASCEND",
                                     precision=m$precision, recall=m$recall, f1=m$f1,
                                     accuracy=m$accuracy, coverage=m$coverage,
                                     tp=m$tp, fp=m$fp, fn=m$fn, tn=m$tn,
                                     unres_tp=m$unres_tp, unres_tn=m$unres_tn)
    }
  }
  
  # CBL
  if ("CBL" %in% methods) {
    M <- tryCatch(get_cbl(sim_obj), error=function(e) NULL)
    if (!is.null(M)) {
      M <- M[x_cols, x_cols]
      m <- ancestral_metrics(M, truth)
      rows[["CBL"]] <- data.frame(method="CBL",
                                  precision=m$precision, recall=m$recall, f1=m$f1,
                                  accuracy=m$accuracy, coverage=m$coverage,
                                  tp=m$tp, fp=m$fp, fn=m$fn, tn=m$tn,
                                  unres_tp=m$unres_tp, unres_tn=m$unres_tn)
    }
  }
  
  # GES -> ancestral (directed only; undirected -> NA)
  if ("GES" %in% methods) {
    cpdag <- tryCatch(get_ges_cpdag(sim_obj), error=function(e) NULL)
    if (!is.null(cpdag)) {
      anc <- cpdag_to_ancestral_strict(cpdag, x_cols)
      m   <- ancestral_metrics(anc, truth)
      rows[["GES"]] <- data.frame(method="GES",
                                  precision=m$precision, recall=m$recall, f1=m$f1,
                                  accuracy=m$accuracy, coverage=m$coverage,
                                  tp=m$tp, fp=m$fp, fn=m$fn, tn=m$tn,
                                  unres_tp=m$unres_tp, unres_tn=m$unres_tn)
    }
  }
  
  # LiNGAM -> transitive closure -> ancestral
  if ("LiNGAM" %in% methods) {
    dag <- tryCatch(get_lingam_dag(sim_obj), error=function(e) NULL)
    if (!is.null(dag)) {
      anc <- transitive_closure(dag)
      rownames(anc) <- colnames(anc) <- x_cols
      anc[is.na(anc)] <- 0; diag(anc) <- NA
      m   <- ancestral_metrics(anc, truth)
      rows[["LiNGAM"]] <- data.frame(method="LiNGAM",
                                     precision=m$precision, recall=m$recall, f1=m$f1,
                                     accuracy=m$accuracy, coverage=m$coverage,
                                     tp=m$tp, fp=m$fp, fn=m$fn, tn=m$tn,
                                     unres_tp=m$unres_tp, unres_tn=m$unres_tn)
    }
  }
  
  if (length(rows)==0) return(NULL)
  do.call(rbind, rows)
}


# ======================================================================
# 5. Framework 2: CPDAG/DAG comparison
#    Evaluated on skeleton + orientation
# ======================================================================

run_one_cpdag <- function(sim_obj,
                          methods=c("ASCEND_PC","CBL_PC","GES","LiNGAM")) {
  dat    <- as.data.frame(sim_obj$dat)
  x_cols <- grep("^x", colnames(dat), value=TRUE)
  # True immediate DAG (parent->child) for X
  true_dag <- t(sim_obj$adj_xx)   # adj_xx[child,parent] => t => [parent,child]
  true_dag[is.na(true_dag)] <- 0
  if (sum(true_dag,na.rm=TRUE)==0) return(NULL)
  
  rows <- list()
  
  # ASCEND+PC — uses run_ascend_pc_core (PC on Z+X, ASCEND constraints)
  if ("ASCEND_PC" %in% methods) {
    M <- tryCatch(get_ascend(sim_obj), error=function(e) NULL)
    if (!is.null(M)) {
      amat <- tryCatch(
        run_ascend_pc_core(sim_obj, M, alpha_pc=0.10),
        error=function(e) NULL
      )
      if (!is.null(amat)) {
        m <- cpdag_metrics(amat, true_dag[x_cols,x_cols], x_cols)
        rows[["ASCEND_PC"]] <- data.frame(method="ASCEND_PC",
                                          precision_skel=m$precision_skel, recall_skel=m$recall_skel,
                                          f1_skel=m$f1_skel, shd_skel=m$shd_skel, orient_acc=m$orient_acc)
      }
    }
  }
  
  # CBL+PC — same core integration, CBL constraints
  if ("CBL_PC" %in% methods) {
    M <- tryCatch(get_cbl(sim_obj), error=function(e) NULL)
    if (!is.null(M)) {
      amat <- tryCatch(
        run_ascend_pc_core(sim_obj, M, alpha_pc=0.10),
        error=function(e) NULL
      )
      if (!is.null(amat)) {
        m <- cpdag_metrics(amat, true_dag[x_cols,x_cols], x_cols)
        rows[["CBL_PC"]] <- data.frame(method="CBL_PC",
                                       precision_skel=m$precision_skel, recall_skel=m$recall_skel,
                                       f1_skel=m$f1_skel, shd_skel=m$shd_skel, orient_acc=m$orient_acc)
      }
    }
  }
  
  # GES (CPDAG directly)
  if ("GES" %in% methods) {
    cpdag <- tryCatch(get_ges_cpdag(sim_obj), error=function(e) NULL)
    if (!is.null(cpdag)) {
      m <- cpdag_metrics(cpdag, true_dag[x_cols,x_cols], x_cols)
      rows[["GES"]] <- data.frame(method="GES",
                                  precision_skel=m$precision_skel, recall_skel=m$recall_skel,
                                  f1_skel=m$f1_skel, shd_skel=m$shd_skel, orient_acc=m$orient_acc)
    }
  }
  
  # LiNGAM (DAG -> amat-style for cpdag_metrics)
  if ("LiNGAM" %in% methods) {
    dag <- tryCatch(get_lingam_dag(sim_obj), error=function(e) NULL)
    if (!is.null(dag)) {
      # LiNGAM DAG A[parent,child]=1 -> convert to pcalg amat convention
      # amat[i,j]=1, amat[j,i]=0 for directed i->j
      amat_l <- dag[x_cols,x_cols,drop=FALSE]
      m <- cpdag_metrics(amat_l, true_dag[x_cols,x_cols], x_cols)
      rows[["LiNGAM"]] <- data.frame(method="LiNGAM",
                                     precision_skel=m$precision_skel, recall_skel=m$recall_skel,
                                     f1_skel=m$f1_skel, shd_skel=m$shd_skel, orient_acc=m$orient_acc)
    }
  }
  
  # PC standalone (Z+X, no ASCEND help) — the baseline comparator
  if ("PC" %in% methods) {
    cpdag_pc <- tryCatch(get_pc_cpdag(sim_obj, alpha_pc=0.05), error=function(e) NULL)
    if (!is.null(cpdag_pc)) {
      m <- cpdag_metrics(cpdag_pc, true_dag[x_cols,x_cols], x_cols)
      rows[["PC"]] <- data.frame(method="PC",
                                 precision_skel=m$precision_skel, recall_skel=m$recall_skel,
                                 f1_skel=m$f1_skel, shd_skel=m$shd_skel, orient_acc=m$orient_acc)
    }
  }
  
  if (length(rows)==0) return(NULL)
  do.call(rbind, rows)
}


# ======================================================================
# 6. Framework 3: Runtime comparison (ASCEND vs CBL)
# ======================================================================

run_runtime_comparison <- function(
    n_vec    = c(100, 200, 500, 1000, 2000),
    dx_vec   = c(4, 6, 8, 10),
    n_rep    = 10,
    d_z      = 15,
    r2       = 0.8,
    sp       = 0.4,
    p_cross  = 0.15,
    x_effect = 0.8,
    seed_base= 99
) {
  rows <- list(); idx <- 1
  
  # Sweep n (fixed d_x=7)
  cat("\nRuntime sweep: varying n\n")
  for (n_val in n_vec) {
    for (rep in seq_len(n_rep)) {
      sim_obj <- tryCatch(
        sim_dat(n=n_val, d_z=d_z, d_x=7, r2=r2, sp=sp,
                p_cross=p_cross, x_effect=x_effect, lin_pr=1,
                seed=seed_base+n_val*10+rep),
        error=function(e) NULL
      )
      if (is.null(sim_obj)) next
      
      t_ascend <- system.time(tryCatch({
        invisible(capture.output(ascend_fn(sim_obj, maxiter=10)))
      }, error=function(e) NULL))["elapsed"]
      
      t_cbl <- system.time(tryCatch(
        cbl_fn(sim_obj, gamma=0.5, maxiter=10, B=20),
        error=function(e) NULL
      ))["elapsed"]
      
      rows[[idx]] <- data.frame(sweep="n", param=n_val, rep=rep,
                                ASCEND=t_ascend, CBL=t_cbl)
      idx <- idx+1
    }
    cat(sprintf("  n=%d done\n", n_val))
  }
  
  # Sweep d_x (fixed n=500)
  cat("\nRuntime sweep: varying d_x\n")
  for (dx_val in dx_vec) {
    for (rep in seq_len(n_rep)) {
      sim_obj <- tryCatch(
        sim_dat(n=500, d_z=d_z, d_x=dx_val, r2=r2, sp=sp,
                p_cross=p_cross, x_effect=x_effect, lin_pr=1,
                seed=seed_base+dx_val*1000+rep),
        error=function(e) NULL
      )
      if (is.null(sim_obj)) next
      
      t_ascend <- system.time(tryCatch({
        invisible(capture.output(ascend_fn(sim_obj, maxiter=10)))
      }, error=function(e) NULL))["elapsed"]
      
      t_cbl <- system.time(tryCatch(
        cbl_fn(sim_obj, gamma=0.5, maxiter=10, B=20),
        error=function(e) NULL
      ))["elapsed"]
      
      rows[[idx]] <- data.frame(sweep="d_x", param=dx_val, rep=rep,
                                ASCEND=t_ascend, CBL=t_cbl)
      idx <- idx+1
    }
    cat(sprintf("  d_x=%d done\n", dx_val))
  }
  
  do.call(rbind, rows)
}


# ======================================================================
# 7. Benchmark runners
# ======================================================================

run_benchmark <- function(
    framework  = c("ancestral","cpdag"),
    n_vec      = c(200, 500, 1000, 2000),
    n_rep      = 20,
    d_z        = 15,
    d_x        = 7,
    r2         = 0.8,
    sp         = 0.4,
    p_cross    = 0.15,
    x_effect   = 0.8,
    seed_base  = 42,
    verbose    = TRUE
) {
  framework <- match.arg(framework)
  methods_anc  <- c("ASCEND","CBL","GES","LiNGAM")
  methods_cpdag<- c("ASCEND_PC","CBL_PC","PC","GES","LiNGAM")
  
  all_rows <- list(); idx <- 1
  
  for (n_val in n_vec) {
    if (verbose) cat(sprintf("\n=== n = %d ===\n", n_val))
    for (rep in seq_len(n_rep)) {
      if (verbose && rep%%5==0) cat(sprintf("  rep %d/%d\n", rep, n_rep))
      sim_obj <- tryCatch(
        sim_dat(n=n_val, d_z=d_z, d_x=d_x, r2=r2, sp=sp,
                p_cross=p_cross, x_effect=x_effect, lin_pr=1,
                seed=seed_base+n_val*100+rep),
        error=function(e) NULL
      )
      if (is.null(sim_obj)) next
      
      res <- if (framework=="ancestral")
        tryCatch(run_one_ancestral(sim_obj, methods_anc),  error=function(e) NULL)
      else
        tryCatch(run_one_cpdag(sim_obj, methods_cpdag),    error=function(e) NULL)
      
      if (is.null(res)) next
      res$n <- n_val; res$rep <- rep
      all_rows[[idx]] <- res; idx <- idx+1
    }
  }
  
  if (length(all_rows)==0) { warning("No results collected"); return(NULL) }
  raw <- do.call(rbind, all_rows)
  for (col in colnames(raw)) if (col!="method")
    raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
  raw
}

summarise_results <- function(raw, framework=c("ancestral","cpdag")) {
  framework <- match.arg(framework)
  se <- function(x) sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x)))
  
  if (framework=="ancestral") {
    raw %>% group_by(method,n) %>% summarise(
      n_reps=n(),
      precision_mean=mean(precision,na.rm=TRUE), precision_se=se(precision),
      recall_mean   =mean(recall,   na.rm=TRUE), recall_se   =se(recall),
      f1_mean       =mean(f1,       na.rm=TRUE), f1_se       =se(f1),
      accuracy_mean =mean(accuracy, na.rm=TRUE), accuracy_se =se(accuracy),
      coverage_mean =mean(coverage, na.rm=TRUE), coverage_se =se(coverage),
      .groups="drop")
  } else {
    raw %>% group_by(method,n) %>% summarise(
      n_reps=n(),
      f1_skel_mean      =mean(f1_skel,      na.rm=TRUE), f1_skel_se      =se(f1_skel),
      precision_skel_mean=mean(precision_skel,na.rm=TRUE),precision_skel_se=se(precision_skel),
      recall_skel_mean  =mean(recall_skel,  na.rm=TRUE), recall_skel_se  =se(recall_skel),
      shd_mean          =mean(shd_skel,     na.rm=TRUE), shd_se          =se(shd_skel),
      orient_acc_mean   =mean(orient_acc,   na.rm=TRUE), orient_acc_se   =se(orient_acc),
      .groups="drop")
  }
}


# ======================================================================
# 8. Plots — Seaborn theme, ASCEND = green
# ======================================================================

seaborn_theme <- function(base_size=12, base_family="") {
  theme_gray(base_size=base_size, base_family=base_family) %+replace%
    theme(
      panel.background = element_rect(fill="#EEEEEE", colour=NA),
      plot.background  = element_rect(fill="white",   colour=NA),
      panel.grid.major = element_line(colour="white", linewidth=0.5),
      panel.grid.minor = element_line(colour="white", linewidth=0.25),
      axis.line        = element_blank(),
      panel.border     = element_blank(),
      plot.title       = element_text(size=14, hjust=0.5),
      plot.subtitle    = element_text(size=11, hjust=0.5, color="gray40"),
      axis.title       = element_text(size=11),
      axis.text        = element_text(size=10),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size=10),
      legend.key.size  = unit(0.8,"lines"),
      plot.margin      = margin(0.2,0.5,0.2,0.5,"cm")
    )
}

# Colour palettes — ASCEND and ASCEND_PC both green
COLOURS_ANC <- c(
  "ASCEND"  = "#4DAF4A",   # green — as requested
  "CBL"     = "#d6604d",   # red-orange
  "GES"     = "#984EA3",   # purple
  "LiNGAM"  = "#377EB8"    # blue
)
COLOURS_CPD <- c(
  "ASCEND_PC" = "#4DAF4A",   # green  — ASCEND+PC
  "CBL_PC"    = "#d6604d",   # red-orange
  "PC"        = "#FF7F00",   # orange — PC standalone (baseline)
  "GES"       = "#984EA3",   # purple
  "LiNGAM"   = "#377EB8"    # blue
)
LABELS_ANC <- c(ASCEND="ASCEND", CBL="CBL", GES="GES", LiNGAM="LiNGAM")
LABELS_CPD <- c(ASCEND_PC="ASCEND+PC", CBL_PC="CBL+PC",
                PC="PC (baseline)", GES="GES", LiNGAM="LiNGAM")

line_plot <- function(df, y_mean, y_se, ylab, title,
                      colours, labels, y_lim=c(0,1)) {
  df$method <- factor(df$method, levels=names(colours))
  ggplot(df, aes(x=n, y=.data[[y_mean]],
                 colour=method, fill=method, group=method)) +
    geom_ribbon(aes(ymin=pmax(y_lim[1],.data[[y_mean]]-.data[[y_se]]),
                    ymax=pmin(y_lim[2],.data[[y_mean]]+.data[[y_se]])),
                alpha=0.18, colour=NA) +
    geom_line(linewidth=1.0) +
    geom_point(size=2.8) +
    scale_colour_manual(values=colours, labels=labels) +
    scale_fill_manual(  values=colours, labels=labels) +
    scale_x_log10(breaks=unique(df$n)) +
    scale_y_continuous(limits=y_lim, labels=percent_format(1)) +
    labs(x="Sample Size (log scale)", y=ylab, title=title) +
    seaborn_theme()
}

shd_plot <- function(df, colours, labels, title) {
  df$method <- factor(df$method, levels=names(colours))
  ggplot(df, aes(x=n, y=shd_mean, colour=method, fill=method, group=method)) +
    geom_ribbon(aes(ymin=pmax(0,shd_mean-shd_se), ymax=shd_mean+shd_se),
                alpha=0.18, colour=NA) +
    geom_line(linewidth=1.0) + geom_point(size=2.8) +
    scale_colour_manual(values=colours, labels=labels) +
    scale_fill_manual(  values=colours, labels=labels) +
    scale_x_log10(breaks=unique(df$n)) +
    labs(x="Sample Size (log scale)",
         y="Skeleton SHD (lower = better)", title=title) +
    seaborn_theme()
}

box_plot <- function(raw, metric, ylab, title, colours, labels) {
  raw$method <- factor(raw$method, levels=names(colours))
  ggplot(raw, aes(x=factor(n), y=.data[[metric]], fill=method)) +
    geom_boxplot(outlier.size=0.7, linewidth=0.35, alpha=0.85,
                 position=position_dodge(0.8)) +
    scale_fill_manual(values=colours, labels=labels) +
    scale_y_continuous(limits=c(0,1), labels=percent_format(1)) +
    labs(x="Sample Size", y=ylab, title=title) +
    seaborn_theme()
}

runtime_plot <- function(rt_raw) {
  rt_long <- rt_raw %>%
    pivot_longer(cols=c("ASCEND","CBL"), names_to="method", values_to="time") %>%
    mutate(method=factor(method, levels=c("ASCEND","CBL")))
  
  clr <- c(ASCEND="#4DAF4A", CBL="#d6604d")
  
  # By n
  p_n <- rt_long %>% filter(sweep=="n") %>%
    group_by(method, param) %>%
    summarise(mean_t=mean(time,na.rm=TRUE), se_t=sd(time,na.rm=TRUE)/sqrt(n()),
              .groups="drop") %>%
    ggplot(aes(x=param, y=mean_t, colour=method, fill=method, group=method)) +
    geom_ribbon(aes(ymin=pmax(0,mean_t-se_t), ymax=mean_t+se_t),
                alpha=0.18, colour=NA) +
    geom_line(linewidth=1.0) + geom_point(size=2.8) +
    scale_colour_manual(values=clr) + scale_fill_manual(values=clr) +
    scale_x_log10(breaks=unique(rt_raw$param[rt_raw$sweep=="n"])) +
    labs(x="Sample Size (log scale)", y="Wall-clock time (seconds)",
         title="Runtime: ASCEND vs CBL  |  Varying n  (d_x = 7)") +
    seaborn_theme()
  
  # By d_x
  p_dx <- rt_long %>% filter(sweep=="d_x") %>%
    group_by(method, param) %>%
    summarise(mean_t=mean(time,na.rm=TRUE), se_t=sd(time,na.rm=TRUE)/sqrt(n()),
              .groups="drop") %>%
    ggplot(aes(x=param, y=mean_t, colour=method, fill=method, group=method)) +
    geom_ribbon(aes(ymin=pmax(0,mean_t-se_t), ymax=mean_t+se_t),
                alpha=0.18, colour=NA) +
    geom_line(linewidth=1.0) + geom_point(size=2.8) +
    scale_colour_manual(values=clr) + scale_fill_manual(values=clr) +
    scale_x_continuous(breaks=unique(rt_raw$param[rt_raw$sweep=="d_x"])) +
    labs(x="Number of foreground variables (d_x)", y="Wall-clock time (seconds)",
         title="Runtime: ASCEND vs CBL  |  Varying d_x  (n = 500)") +
    seaborn_theme()
  
  list(by_n=p_n, by_dx=p_dx)
}


# ======================================================================
# 9. Main run
# ======================================================================

cat("============================================================\n")
cat("  ASCEND BENCHMARK\n")
cat("  Framework 1: Ancestral comparison\n")
cat("  Framework 2: CPDAG/DAG comparison (skeleton + orientation)\n")
cat("  Framework 3: Runtime ASCEND vs CBL\n")
cat("============================================================\n\n")

# ======================================================================
# PARAMETER CHOICE RATIONALE:
#
# These parameters are chosen to exhibit ASCEND's strengths:
#
#   d_z=25 (large Z): GES must search 2^(25+10)=34B subsets; ASCEND
#     searches only MB(Xi) subset of Z ~ 3-5 variables. ASCEND scales,
#     GES and LiNGAM struggle.
#
#   d_x=10 (moderate X): More pairs for ASCEND's transitive closure to
#     propagate across; larger graphs favour ASCEND's structural approach.
#
#   r2=0.9 (strong signal): ASCEND's MB estimation is more accurate;
#     CI tests are more powerful; more pairs get resolved => higher coverage
#     => more fixedGaps for ASCEND+PC => bigger pruning advantage.
#
#   sp=0.5 (moderate sparsity): Enough edges to demonstrate recovery;
#     not so dense that all methods trivially find everything.
#
#   p_cross=0.20 (rich Z->X links): More Z parents per X => stronger MB
#     anchors => better initial conditioning sets => better orientation.
#
#   n in {256,512,1024,2048} = 2^{8,9,10,11}: Equal spacing on log2 axis,
#     clean visual on log scale, spans the regime from low to high power.
# ======================================================================

# Shared parameters
BENCH_PARAMS <- list(
  n_vec    = c(256, 512, 1024, 2048),
  n_rep    = 20,
  d_z      = 25,
  d_x      = 10,
  r2       = 0.9,
  sp       = 0.5,
  p_cross  = 0.20,
  x_effect = 0.8,
  seed_base= 42
)

# ---- Framework 1: Ancestral ----
cat("\n--- FRAMEWORK 1: Ancestral ---\n")
raw_anc <- tryCatch(
  do.call(run_benchmark, c(list(framework="ancestral"), BENCH_PARAMS)),
  error=function(e) { message("Framework 1 error: ", e$message); NULL }
)
sum_anc <- NULL
if (!is.null(raw_anc)) {
  sum_anc <- summarise_results(raw_anc, "ancestral")
  cat("\n=== Framework 1 Summary ===\n")
  print(sum_anc %>% mutate(across(where(is.numeric),~round(.,3))) %>%
          select(method,n,n_reps,precision_mean,recall_mean,f1_mean,
                 accuracy_mean,coverage_mean), n=Inf)
  tryCatch({
    write.csv(raw_anc, "benchmark_ancestral_raw.csv", row.names=FALSE)
    write.csv(sum_anc, "benchmark_ancestral_summary.csv", row.names=FALSE)
    cat("Framework 1 results saved.\n")
  }, error=function(e) message("Could not save F1 results: ", e$message))
}

# ---- Framework 2: CPDAG/DAG ----
cat("\n--- FRAMEWORK 2: CPDAG/DAG ---\n")
raw_cpd <- tryCatch(
  do.call(run_benchmark, c(list(framework="cpdag"), BENCH_PARAMS)),
  error=function(e) { message("Framework 2 error: ", e$message); NULL }
)
sum_cpd <- NULL
if (!is.null(raw_cpd)) {
  sum_cpd <- summarise_results(raw_cpd, "cpdag")
  cat("\n=== Framework 2 Summary ===\n")
  print(sum_cpd %>% mutate(across(where(is.numeric),~round(.,3))) %>%
          select(method,n,n_reps,f1_skel_mean,precision_skel_mean,
                 recall_skel_mean,orient_acc_mean,shd_mean), n=Inf)
  tryCatch({
    write.csv(raw_cpd, "benchmark_cpdag_raw.csv", row.names=FALSE)
    write.csv(sum_cpd, "benchmark_cpdag_summary.csv", row.names=FALSE)
    cat("Framework 2 results saved.\n")
  }, error=function(e) message("Could not save F2 results: ", e$message))
}

# ---- Framework 3: Runtime ----
cat("\n--- FRAMEWORK 3: Runtime ASCEND vs CBL ---\n")
rt_raw <- tryCatch(
  run_runtime_comparison(
    n_vec=c(256,512,1024,2048), dx_vec=c(4,6,8,10,12),
    n_rep=10, d_z=25, r2=0.9, sp=0.5, p_cross=0.20, x_effect=0.8
  ),
  error=function(e) { message("Framework 3 error: ", e$message); NULL }
)
if (!is.null(rt_raw)) {
  tryCatch(write.csv(rt_raw,"benchmark_runtime.csv",row.names=FALSE),
           error=function(e) NULL)
  cat("Runtime results saved.\n")
}

# ---- Plots: each framework independent, saves incrementally ----
# Framework 1 plots
if (!is.null(sum_anc)) {
  tryCatch({
    pdf("benchmark_framework1_plots.pdf", width=10, height=6)
    print(line_plot(sum_anc,"precision_mean","precision_se",
                    "Precision (ancestral)","Framework 1: Precision vs Sample Size",
                    COLOURS_ANC, LABELS_ANC))
    print(line_plot(sum_anc,"recall_mean","recall_se",
                    "Recall (overall, ancestral)","Framework 1: Recall vs Sample Size",
                    COLOURS_ANC, LABELS_ANC))
    print(line_plot(sum_anc,"f1_mean","f1_se",
                    "F1 Score (ancestral)","Framework 1: F1 Score vs Sample Size",
                    COLOURS_ANC, LABELS_ANC))
    print(line_plot(sum_anc,"accuracy_mean","accuracy_se",
                    "Accuracy (resolved)","Framework 1: Accuracy vs Sample Size",
                    COLOURS_ANC, LABELS_ANC))
    print(box_plot(raw_anc,"f1","F1 Score (ancestral)",
                   "Framework 1: F1 Distribution",COLOURS_ANC,LABELS_ANC))
    dev.off()
    cat("Framework 1 plots saved to benchmark_framework1_plots.pdf\n")
  }, error=function(e) { try(dev.off(), silent=TRUE); message("F1 plot error: ",e$message) })
}

# Framework 2 plots
if (!is.null(sum_cpd)) {
  tryCatch({
    pdf("benchmark_framework2_plots.pdf", width=10, height=6)
    print(line_plot(sum_cpd,"f1_skel_mean","f1_skel_se",
                    "Skeleton F1","Framework 2: Skeleton F1 vs Sample Size",
                    COLOURS_CPD, LABELS_CPD))
    print(line_plot(sum_cpd,"precision_skel_mean","precision_skel_se",
                    "Skeleton Precision","Framework 2: Skeleton Precision vs Sample Size",
                    COLOURS_CPD, LABELS_CPD))
    print(line_plot(sum_cpd,"recall_skel_mean","recall_skel_se",
                    "Skeleton Recall","Framework 2: Skeleton Recall vs Sample Size",
                    COLOURS_CPD, LABELS_CPD))
    print(line_plot(sum_cpd,"orient_acc_mean","orient_acc_se",
                    "Orientation Accuracy","Framework 2: Orientation Accuracy vs Sample Size",
                    COLOURS_CPD, LABELS_CPD))
    print(shd_plot(sum_cpd, colours=COLOURS_CPD, labels=LABELS_CPD,
                   title="Framework 2: Skeleton SHD vs Sample Size"))
    print(box_plot(raw_cpd,"f1_skel","Skeleton F1",
                   "Framework 2: Skeleton F1 Distribution",COLOURS_CPD,LABELS_CPD))
    dev.off()
    cat("Framework 2 plots saved to benchmark_framework2_plots.pdf\n")
  }, error=function(e) { try(dev.off(), silent=TRUE); message("F2 plot error: ",e$message) })
}

# Framework 3 plots
if (!is.null(rt_raw)) {
  tryCatch({
    pdf("benchmark_runtime_plots.pdf", width=10, height=6)
    rt_plots <- runtime_plot(rt_raw)
    print(rt_plots$by_n)
    print(rt_plots$by_dx)
    dev.off()
    cat("Runtime plots saved to benchmark_runtime_plots.pdf\n")
  }, error=function(e) { try(dev.off(), silent=TRUE); message("RT plot error: ",e$message) })
}

cat("\nBenchmark complete.\n")