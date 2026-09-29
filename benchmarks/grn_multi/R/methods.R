## ======================================================================
## benchmarks/grn_multi/R/methods.R
## One function per R method. Each takes a dataset (read_dataset()) and the
## method's settings, and returns a list with some of
##   score  X x X confidence matrix [regulator, target]
##   adj    X x X native graph (a -> b is adj[a, b] = 1; undirected = both)
##   anc    X x X ancestral matrix (eval_metrics.R convention)
##   stats  named list of extra numbers to keep (tests, sizes, ...)
## Every method sees the same input: all background AND foreground
## variables. Only the X x X block is kept for scoring.
##
## Reviewer comments: 1 (baselines), 10 (tier-aware PC and score-based search)
## Packages are loaded inside each function, so a method whose package is
## missing fails on its own and the others still run.
## ======================================================================

## --- ASCEND -----------------------------------------------------------------
m_ascend <- function(D, s) {
  source(file.path(ROOT, "R", "ascend.R"), local = TRUE)
  sim <- list(dat = D$dat[, c(D$zlabs, D$xlabs)], params = list(d_x = length(D$xlabs)))
  M <- suppressMessages(ascend(sim, maxiter = s$maxiter, alpha = s$alpha,
                               alpha_mb = s$alpha_mb, fdr = TRUE,
                               min_votes = s$min_votes, prescreen = s$prescreen,
                               verbose = FALSE))
  xl <- D$xlabs
  anc <- unclass(M)[xl, xl]; attributes(anc) <- list(dim = dim(anc), dimnames = list(xl, xl))
  # Ranking score for AUPR/AUROC: the p-value of the last pairwise test
  # (small = strong dependence), stored as -log10 p; NA if never tested.
  pp <- attr(M, "pair_p")
  score <- if (!is.null(pp)) -log10(pmax(pp[xl, xl], 1e-300)) else NULL
  # Native direct-edge graph: the minimal edges implied by the ancestor claims
  adj <- transitive_reduction((!is.na(anc) & anc == 1) * 1)
  st <- attr(M, "stats")
  list(anc = anc, adj = adj, score = score,
       stats = if (is.null(st)) list() else lapply(unlist(st), identity))
}

## --- tiered PC (tpc: PC-stable with tiered background knowledge) -------------
m_pc_tiered <- function(D, s) {
  suppressPackageStartupMessages(library(tpc))
  all <- c(D$zlabs, D$xlabs)
  X <- as.matrix(D$dat[, all])
  fit <- tpc::tpc(suffStat = list(C = cor(X), n = nrow(X)), indepTest = pcalg::gaussCItest,
                  alpha = s$alpha, labels = all,
                  tiers = c(rep(1, length(D$zlabs)), rep(2, length(D$xlabs))))
  # graph slot: [a, b] = 1 means a -> b (both = undirected). NOT as(fit, "amat"),
  # which pcalg stores transposed ([b, a] = 1 for a -> b).
  G <- as(fit@graph, "matrix"); dimnames(G) <- list(all, all)
  list(adj = G[D$xlabs, D$xlabs], anc = graph_to_ancestral(G, D$xlabs))
}

## --- tiered hill climbing (score-based search, BIC, tier blacklist) ---------
m_hc_tiered <- function(D, s) {
  suppressPackageStartupMessages(library(bnlearn))
  all <- c(D$zlabs, D$xlabs)
  df  <- as.data.frame(D$dat[, all])
  bl  <- bnlearn::tiers2blacklist(list(D$zlabs, D$xlabs))    # forbids every X -> Z edge
  fit <- bnlearn::hc(df, score = s$score, blacklist = bl)
  G <- bnlearn::amat(fit)[all, all]                            # [a, b] = 1: a -> b
  list(adj = G[D$xlabs, D$xlabs], anc = graph_to_ancestral(G, D$xlabs))
}

## --- GES (pcalg, no tier constraint) ------------------------------------------
m_ges <- function(D, s) {
  suppressPackageStartupMessages(library(pcalg))
  all <- c(D$zlabs, D$xlabs)
  score <- new("GaussL0penObsScore", as.matrix(D$dat[, all]))
  fit <- pcalg::ges(score, labels = all, verbose = FALSE)
  G <- as(fit$essgraph, "matrix") * 1; dimnames(G) <- list(all, all)   # [a, b]: a -> b
  list(adj = G[D$xlabs, D$xlabs], anc = graph_to_ancestral(G, D$xlabs))
}

## --- CBL (Watson & Silva 2022) -----------------------------------------------
m_cbl <- function(D, s) {
  suppressPackageStartupMessages({ library(data.table); library(dplyr); library(foreach) })
  e <- new.env()
  assign("ROOT", ROOT, envir = e)
  sys.source(file.path(ROOT, "R", "cbl.R"), envir = e)
  # lin_pr = 1: CBL uses its lasso learner, the same model class as ASCEND's test
  sim <- list(dat = D$dat[, c(D$zlabs, D$xlabs)], params = list(lin_pr = 1))
  m <- e$cbl_fn(sim, gamma = s$gamma, maxiter = s$maxiter, B = s$B)
  # cbl_fn returns m[descendant, ancestor]; convert to [ancestor, descendant]
  anc <- e$cbl_as_ancestral(m)[D$xlabs, D$xlabs]
  diag(anc) <- NA
  list(anc = anc, adj = transitive_reduction((!is.na(anc) & anc == 1) * 1))
}

## --- GENIE3 --------------------------------------------------------------------
m_genie3 <- function(D, s) {
  suppressPackageStartupMessages(library(GENIE3))
  all <- c(D$zlabs, D$xlabs)
  E <- t(as.matrix(D$dat[, all])); rownames(E) <- all
  W <- GENIE3::GENIE3(E, targets = D$xlabs, nTrees = s$nTrees, K = s$K, nCores = 1)
  list(score = W[D$xlabs, D$xlabs])                         # [regulator, target]
}

## --- minet: CLR, MRNET, ARACNe -------------------------------------------------
minet_mim <- function(D, s) {
  suppressPackageStartupMessages(library(minet))
  all <- c(D$zlabs, D$xlabs)
  minet::build.mim(as.data.frame(D$dat[, all]), estimator = s$estimator)
}
m_clr    <- function(D, s) list(score = minet::clr(minet_mim(D, s))[D$xlabs, D$xlabs])
m_mrnet  <- function(D, s) list(score = minet::mrnet(minet_mim(D, s))[D$xlabs, D$xlabs])
m_aracne <- function(D, s) list(score = minet::aracne(minet_mim(D, s), eps = 0)[D$xlabs, D$xlabs])

## --- WGCNA (topological overlap of the soft-thresholded adjacency) ------------
m_wgcna <- function(D, s) {
  suppressPackageStartupMessages(library(WGCNA))
  all <- c(D$zlabs, D$xlabs)
  X <- as.matrix(D$dat[, all])
  powers <- c(1:10, seq(12, 20, 2))
  invisible(capture.output(sft <- WGCNA::pickSoftThreshold(
    X, powerVector = powers, networkType = s$type, RsquaredCut = s$RsquaredCut, verbose = 0)))
  beta <- sft$powerEstimate
  if (length(beta) == 0 || is.na(beta)) beta <- powers[which.max(sft$fitIndices[, "SFT.R.sq"])]
  if (length(beta) == 0 || !is.finite(beta)) beta <- 6
  A <- WGCNA::adjacency(X, power = beta, type = s$type)
  TOM <- WGCNA::TOMsimilarity(A, TOMType = if (s$type == "unsigned") "unsigned" else "signed", verbose = 0)
  dimnames(TOM) <- list(all, all)
  list(score = TOM[D$xlabs, D$xlabs], stats = list(power = beta))
}

## --- PPCOR (full-order partial correlation, as in BEELINE) -----------------------
m_ppcor <- function(D, s) {
  all <- c(D$zlabs, D$xlabs)
  X <- as.matrix(D$dat[, all])
  P <- tryCatch(ppcor::pcor(X)$estimate, error = function(e) NULL)
  shrunk <- is.null(P) || any(!is.finite(P))
  if (shrunk) {           # singular covariance: shrinkage estimate instead
    Om <- solve(cov(X) + diag(1e-6 * mean(diag(cov(X))), ncol(X)))
    P <- -cov2cor(Om)
  }
  dimnames(P) <- list(all, all)
  list(score = abs(P[D$xlabs, D$xlabs]), stats = list(shrunk = shrunk))
}

## --- absolute Pearson correlation --------------------------------------------------
m_pearson <- function(D, s) {
  C <- abs(cor(as.matrix(D$dat[, D$xlabs])))
  list(score = C)
}

R_METHODS <- list(ascend = m_ascend, pc_tiered = m_pc_tiered, hc_tiered = m_hc_tiered,
                  ges = m_ges, cbl = m_cbl, genie3 = m_genie3, clr = m_clr,
                  mrnet = m_mrnet, aracne = m_aracne, wgcna = m_wgcna,
                  ppcor = m_ppcor, pearson = m_pearson)
