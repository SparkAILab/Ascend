## ======================================================================
## benchmarks/grn_multi/R/scoring.R
## Every metric of the benchmark, for every method class, against both
## ground truths. Sourced by 03_score.R; tested by tests/test_grn_multi.R.
##
## Reviewer comments: 2 (AUPR/AUROC everywhere; oriented precision/recall,
##   ancestor-descendant accuracy, SHD), 3 (BEELINE's AUPR ratio and early
##   precision ratio), 11 (orientation), minor 9 (coverage)
##
## Two ground truths (both X x X, [a, b] = 1 means a -> b / a ancestor of b):
##   direct : the true regulatory edges (what GRN methods estimate)
##   anc    : its transitive closure (what ASCEND estimates)
## Causal methods are scored on their ancestral matrix against "anc" and on
## their native graph against "direct" (ASCEND's native graph is the
## transitive reduction of its ancestor claims).
##
## Metric families (one row per metric in the long output)
##   Ranking, all methods (unordered pairs; score = max of both directions)
##     aupr        area under the precision-recall curve (Davis-Goadrich
##                 interpolation, correct with tied scores)
##     auroc       area under the ROC curve (ties count one half)
##     aupr_ratio  aupr / edge density (BEELINE; 1 = random)
##     epr         early precision ratio: precision in the top-k pairs,
##                 k = number of true pairs, / edge density (BEELINE)
##     prec_Ktrue, f1_Ktrue           top-k with k = number of true pairs
##     prec_Kasc, rec_Kasc, f1_Kasc   top-k with k = number of pairs ASCEND
##                 calls related on the same dataset (matched budget)
##     Ties at the cut-off are split in expectation (random tie-breaking).
##   Direction, directed and causal methods (ordered pairs)
##     dir_aupr, dir_auroc     score[a, b] ranks the ordered pair a -> b
##     orient_acc_rank         among true a -> b pairs: P(score[a,b] > score[b,a])
##   Native graph, causal methods (eval_metrics.R::eval_ancestral, strict:
##     only 1 counts as a claim; len_* also count 0.5 "non-descendant")
##     skel_precision/recall/f1, dir_precision/recall/f1, orient_acc,
##     orient_undet, ad_acc, shd, coverage, n_claims, n_reversed
##   Directed GRN methods at the matched budget
##     *_Kasc versions of the native-graph metrics, after keeping the top
##     K_asc pairs and orienting each towards the larger score
## ======================================================================

## --- ranking helpers -------------------------------------------------------
auroc_ties <- function(sc, lb) {
  np <- sum(lb == 1); nn <- sum(lb == 0)
  if (np == 0 || nn == 0) return(NA_real_)
  r <- rank(sc, ties.method = "average")
  (sum(r[lb == 1]) - np * (np + 1) / 2) / (np * nn)
}

aupr_dg <- function(sc, lb) {
  if (sum(lb == 1) == 0 || sum(lb == 0) == 0) return(NA_real_)
  if (length(unique(sc)) == 1) return(mean(lb))          # one tie group = random
  PRROC::pr.curve(scores.class0 = sc[lb == 1], scores.class1 = sc[lb == 0],
                  curve = FALSE, dg.compute = TRUE)$auc.davis.goadrich
}

# Expected TP among the top k of `sc` under random tie-breaking
expected_tp_topk <- function(sc, lb, k) {
  if (k <= 0) return(0)
  k <- min(k, length(sc))
  thr <- sort(sc, decreasing = TRUE)[k]
  above <- sc > thr; tie <- sc == thr
  sum(lb[above]) + (k - sum(above)) * mean(lb[tie])
}

prf_topk <- function(sc, lb, k) {
  if (k <= 0) return(c(prec = NA, rec = NA, f1 = NA))
  tp <- expected_tp_topk(sc, lb, k); k <- min(k, length(sc))
  pr <- tp / k; re <- tp / max(1, sum(lb))
  c(prec = pr, rec = re, f1 = if (pr + re > 0) 2 * pr * re / (pr + re) else 0)
}

upper_pairs <- function(M) {
  ut <- which(upper.tri(M), arr.ind = TRUE); ut
}

# Deterministic top-k of a symmetric score (ties broken by pair index), each
# kept pair oriented towards the larger directed score (ties: undirected)
topk_oriented <- function(S, k) {
  d <- nrow(S); Ssym <- pmax(S, t(S)); ut <- upper_pairs(S)
  v <- Ssym[ut]; ord <- order(-v, seq_along(v))
  keep <- ut[ord[seq_len(min(k, length(v)))], , drop = FALSE]
  E <- matrix(0, d, d, dimnames = dimnames(S))
  for (r in seq_len(nrow(keep))) {
    a <- keep[r, 1]; b <- keep[r, 2]
    if (S[a, b] > S[b, a]) E[a, b] <- 1
    else if (S[b, a] > S[a, b]) E[b, a] <- 1
    else E[a, b] <- E[b, a] <- 1
  }
  E
}

## --- the metric set for one (method output, truth) -------------------------
# out: list(class, score, anc, adj) as read from disk
# Tr : truth matrix for this target; K_asc: ASCEND's related-pair count (or NA)
# target: "direct" or "anc"
score_one <- function(out, Tr, target, K_asc) {
  d <- nrow(Tr); xl <- rownames(Tr)
  Tr[is.na(Tr)] <- 0; diag(Tr) <- 0
  Tskel <- ((Tr + t(Tr)) > 0) * 1
  ut <- upper_pairs(Tr); lb <- Tskel[ut]
  K_true <- sum(lb); density <- K_true / length(lb)
  res <- list(K_true = K_true, K_asc = K_asc, density = density)

  # the graph a causal method is judged on for this target
  est <- NULL
  if (out$class == "causal") {
    est <- if (target == "anc") out$anc else out$adj
    est <- est[xl, xl]
  }

  # ranking score over unordered pairs
  if (out$class == "causal") {
    rel <- (!is.na(est) & est == 1) | (!is.na(t(est)) & t(est) == 1)
    na2 <- is.na(est) & is.na(t(est))
    lvl <- ifelse(rel, 2, ifelse(na2, 1, 0))
    if (!is.null(out$score)) {                    # ASCEND: break ties by -log10 p
      s <- out$score[xl, xl]; s[!is.finite(s)] <- 0; s <- pmax(s, t(s))
      lvl <- lvl + 0.9 * s / (max(s) + 1)
    }
    Ssym <- lvl
  } else {
    S <- abs(out$score[xl, xl]); S[!is.finite(S)] <- 0; diag(S) <- 0
    Ssym <- pmax(S, t(S))
  }
  sc <- Ssym[ut]
  res$aupr  <- aupr_dg(sc, lb)
  res$auroc <- auroc_ties(sc, lb)
  res$aupr_ratio <- res$aupr / density
  pk <- prf_topk(sc, lb, K_true)
  res$prec_Ktrue <- pk[["prec"]]; res$f1_Ktrue <- pk[["f1"]]
  res$epr <- pk[["prec"]] / density
  if (is.finite(K_asc)) {
    pa <- prf_topk(sc, lb, K_asc)
    res$prec_Kasc <- pa[["prec"]]; res$rec_Kasc <- pa[["rec"]]; res$f1_Kasc <- pa[["f1"]]
  }

  # direction
  off <- which(row(Tr) != col(Tr))
  if (out$class == "causal") {
    Sdir <- (!is.na(est) & est == 1) * 2 + (is.na(est)) * 1     # claim > unknown > no
  } else if (out$class == "directed") {
    Sdir <- abs(out$score[xl, xl]); Sdir[!is.finite(Sdir)] <- 0
  } else Sdir <- NULL
  if (!is.null(Sdir)) {
    res$dir_aupr  <- aupr_dg(Sdir[off], Tr[off])
    res$dir_auroc <- auroc_ties(Sdir[off], Tr[off])
    tp <- which(Tr == 1, arr.ind = TRUE)
    if (nrow(tp)) {
      f <- Sdir[tp]; b <- Sdir[tp[, 2:1, drop = FALSE]]
      res$orient_acc_rank <- mean((f > b) + 0.5 * (f == b))
    }
  }

  # native graph (causal) or matched-budget graph (directed)
  add_eval <- function(E, pos, suffix = "", prefix = "") {
    m <- eval_ancestral(E, Tr, pos_values = pos)
    keep <- c("skel_precision", "skel_recall", "skel_f1", "dir_precision", "dir_recall",
              "dir_f1", "orient_acc", "orient_undet", "ad_acc", "shd", "coverage",
              "n_reversed")
    v <- m[keep]; v$n_claims <- m$dir_claims
    names(v) <- paste0(prefix, names(v), suffix)
    v
  }
  if (out$class == "causal") {
    res <- c(res, add_eval(est, 1))
    if (any(est == 0.5, na.rm = TRUE)) res <- c(res, add_eval(est, c(1, 0.5), prefix = "len_"))
  } else if (out$class == "directed" && is.finite(K_asc) && K_asc > 0) {
    E <- topk_oriented(abs(out$score[xl, xl]), K_asc)
    v <- add_eval(E, 1, suffix = "_Kasc")
    res <- c(res, v[c("dir_precision_Kasc", "dir_recall_Kasc", "dir_f1_Kasc",
                      "orient_acc_Kasc", "ad_acc_Kasc", "shd_Kasc")])
  }
  res
}

# Number of pairs ASCEND calls related (either direction), per target
ascend_K <- function(out, target, xl) {
  est <- if (target == "anc") out$anc else out$adj
  est <- est[xl, xl]
  rel <- (!is.na(est) & est == 1) | (!is.na(t(est)) & t(est) == 1)
  sum(rel[upper.tri(rel)])
}
