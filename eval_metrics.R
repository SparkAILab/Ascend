# ======================================================================
# Shared evaluation for ancestral / directed-graph output
# ----------------------------------------------------------------------
# One set of definitions used by every benchmark script, so precision,
# recall, coverage and orientation are computed identically everywhere.
#
# Inputs
#   est    d x d matrix over the foreground, est[a, b] in {1, 0.5, 0, NA}
#            1   a is an ancestor of b            (a < b)
#            0.5 a is not a descendant of b       (a <= b)
#            0   no ancestral claim in this direction
#            NA  undetermined
#          Methods that only output ancestors (GES/PC/LiNGAM closures) use
#          1 / 0 / NA. An undirected CPDAG edge may be coded as 1 in both
#          directions or NA in both directions.
#   truth  d x d 0/1 matrix, truth[a, b] = 1 iff a is an ancestor of b.
#          Both matrices must carry the same row/column names; est is
#          reordered to truth's names, so neither needs any particular order.
#
# Pair states (unordered pair {a, b})
#   truth: "a<b", "b<a" or "none"   (none = neither is an ancestor)
#   est  : "a<b", "b<a", "none", "undirected" (related, direction not
#          committed: positive claims in both directions) or "NA"
#          (unresolved: both directions NA)
#
# Definitions (all reported by eval_ancestral())
#   coverage        fraction of unordered pairs whose est state is not "NA"
#   undetermined    1 - coverage
#   Skeleton level (is the pair ancestrally related, either direction?)
#     skel_precision/recall/f1   NA pairs count as "not related" for recall
#   Directed level (ordered pairs; a claim a->b is est[a, b] in pos_values)
#     dir_precision   TP / (all directed claims)
#     dir_recall      TP / (all true ordered ancestral pairs); unresolved
#                     true pairs count as misses
#     dir_f1
#     n_reversed      claims a->b where truth is b->a
#   Orientation (among truly related pairs the method also calls related)
#     orient_acc      correct direction / pairs with a committed direction
#     orient_undet    fraction of those pairs left "undirected"
#   3-class ancestor-descendant accuracy over resolved pairs
#     ad_acc          est state == truth state, over pairs with est state
#                     in {a<b, b<a, none}
#   Structural Hamming distance on the ancestral graph
#     shd             one unit per pair whose est state != truth state;
#                     reversal = 1, undirected on a true edge = 1,
#                     NA treated as "none"
#     shd_resolved    same, over resolved pairs only
#
# pos_values: which est values count as a positive directed claim. The
#   paper's convention counts both 1 and 0.5; use pos_values = 1 for the
#   strict "ancestor" claims only.
# ======================================================================

pair_states <- function(est, truth, pos_values = c(1, 0.5)) {
  nms <- rownames(truth)
  stopifnot(!is.null(nms), all(nms %in% rownames(est)), all(nms %in% colnames(est)))
  est <- est[nms, nms, drop = FALSE]
  d   <- length(nms)
  ut  <- which(upper.tri(matrix(0, d, d)), arr.ind = TRUE)
  a <- ut[, 1]; b <- ut[, 2]

  t_ab <- truth[cbind(a, b)]; t_ba <- truth[cbind(b, a)]
  t_ab[is.na(t_ab)] <- 0; t_ba[is.na(t_ba)] <- 0
  t_state <- ifelse(t_ab == 1, "a<b", ifelse(t_ba == 1, "b<a", "none"))

  e_ab <- est[cbind(a, b)]; e_ba <- est[cbind(b, a)]
  p_ab <- !is.na(e_ab) & e_ab %in% pos_values
  p_ba <- !is.na(e_ba) & e_ba %in% pos_values
  na_pair <- is.na(e_ab) & is.na(e_ba)
  e_state <- ifelse(na_pair, "NA",
                    ifelse(p_ab & p_ba, "undirected",
                           ifelse(p_ab, "a<b",
                                  ifelse(p_ba, "b<a", "none"))))

  data.frame(a = nms[a], b = nms[b], truth = t_state, est = e_state,
             stringsAsFactors = FALSE)
}

eval_ancestral <- function(est, truth, pos_values = c(1, 0.5)) {
  ps <- pair_states(est, truth, pos_values)
  safe_div <- function(x, y) if (y > 0) x / y else NA_real_
  f1_of    <- function(p, r) if (!is.na(p) && !is.na(r) && p + r > 0) 2 * p * r / (p + r) else NA_real_

  n_pairs  <- nrow(ps)
  resolved <- ps$est != "NA"
  t_rel    <- ps$truth != "none"
  e_rel    <- ps$est %in% c("a<b", "b<a", "undirected")

  # skeleton (relatedness) level
  s_tp <- sum(t_rel & e_rel); s_fp <- sum(!t_rel & e_rel); s_fn <- sum(t_rel & !e_rel)
  skel_precision <- safe_div(s_tp, s_tp + s_fp)
  skel_recall    <- safe_div(s_tp, s_tp + s_fn)

  # directed level over ordered pairs; an undirected pair is a claim in both directions
  claim_ab <- ps$est %in% c("a<b", "undirected")
  claim_ba <- ps$est %in% c("b<a", "undirected")
  d_tp <- sum(claim_ab & ps$truth == "a<b") + sum(claim_ba & ps$truth == "b<a")
  d_claims <- sum(claim_ab) + sum(claim_ba)
  d_true   <- sum(t_rel)
  n_reversed <- sum(ps$est == "a<b" & ps$truth == "b<a") + sum(ps$est == "b<a" & ps$truth == "a<b")
  dir_precision <- safe_div(d_tp, d_claims)
  dir_recall    <- safe_div(d_tp, d_true)

  # orientation among truly related pairs called related
  both_rel   <- t_rel & e_rel
  committed  <- both_rel & ps$est != "undirected"
  n_correct  <- sum(committed & ps$est == ps$truth)
  orient_acc   <- safe_div(n_correct, sum(committed))
  orient_undet <- safe_div(sum(both_rel & ps$est == "undirected"), sum(both_rel))

  # 3-class accuracy over resolved, direction-committed pairs
  cls   <- ps$est %in% c("a<b", "b<a", "none")
  ad_acc <- safe_div(sum(cls & ps$est == ps$truth), sum(cls))

  # SHD on the ancestral graph (NA treated as "none")
  e_for_shd <- ifelse(ps$est == "NA", "none", ps$est)
  shd          <- sum(e_for_shd != ps$truth)
  shd_resolved <- sum(resolved & ps$est != ps$truth)

  list(
    n_pairs = n_pairs, n_true_rel = d_true,
    coverage = mean(resolved), undetermined = 1 - mean(resolved),
    skel_tp = s_tp, skel_fp = s_fp, skel_fn = s_fn,
    skel_precision = skel_precision, skel_recall = skel_recall,
    skel_f1 = f1_of(skel_precision, skel_recall),
    dir_tp = d_tp, dir_claims = d_claims, n_reversed = n_reversed,
    dir_precision = dir_precision, dir_recall = dir_recall,
    dir_f1 = f1_of(dir_precision, dir_recall),
    orient_acc = orient_acc, orient_undet = orient_undet,
    orient_n = sum(committed),
    ad_acc = ad_acc,
    shd = shd, shd_resolved = shd_resolved,
    n_undirected = sum(ps$est == "undirected"),
    n_unresolved_true = sum(!resolved & t_rel)
  )
}

# Flat one-row data.frame, convenient for rbind-ing benchmark rows
eval_ancestral_row <- function(est, truth, pos_values = c(1, 0.5), prefix = "") {
  m <- eval_ancestral(est, truth, pos_values)
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  if (nzchar(prefix)) names(df) <- paste0(prefix, names(df))
  df
}

# Skeleton-only scoring for methods that output an undirected weight
# matrix (GENIE3, ARACNe, WGCNA, CLR, ...): thresholded binary skeleton ->
# the orientation metrics are NA by construction.
eval_undirected <- function(A_pred, truth) {
  A_pred <- (as.matrix(A_pred) > 0) * 1
  A_pred <- pmax(A_pred, t(A_pred))
  m <- eval_ancestral(A_pred, truth, pos_values = 1)
  m$dir_precision <- m$dir_recall <- m$dir_f1 <- NA_real_
  m$orient_acc <- m$ad_acc <- NA_real_; m$n_reversed <- NA_integer_
  m
}

# Ground truth from a simulator adjacency (adj_xx[child, parent] = 1)
truth_from_adj <- function(adj_xx) {
  A <- t(adj_xx); A[is.na(A)] <- 0; diag(A) <- 0
  d <- nrow(A); R <- (A > 0) * 1; P <- R
  for (k in seq_len(max(d - 1, 0))) { P <- ((P %*% A) > 0) * 1; R <- ((R + P) > 0) * 1 }
  dimnames(R) <- dimnames(A); diag(R) <- 0
  R
}
