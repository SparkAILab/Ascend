## ======================================================================
## Reviewer comment 11: what drives ASCEND's orientation accuracy?
## ----------------------------------------------------------------------
## For every foreground pair of every simulated run, record the truth, the
## estimate, the step that committed it (R1 / R2 / R3 / transitivity /
## symmetry / final R3, from attr(M, "rule")), the number of agreeing
## witnesses, and structural features of the true graph. Then tabulate:
##   1. contribution and accuracy of each rule
##   2. how often orientation is left undetermined, and why
##   3. orientation accuracy by in-degree of the descendant, by confounding
##      (a shared direct parent; a foreground common ancestor reaching both
##      without going through either), by the number of parallel paths, by
##      adjacency (direct edge vs ancestral only), and by dependence
##      strength (|partial correlation| given the true background parents,
##      a proxy for near-unfaithfulness)
##
## Run from the repo root (base R + R.utils):
##   Rscript analysis/orientation_analysis.R           # 20 seeds x 6 cells
##   SEEDS=5 Rscript analysis/orientation_analysis.R   # quick look
## Writes analysis/out/orientation_pairs.csv (one row per pair per run) and
##        analysis/out/orientation_*.csv summary tables.
## ======================================================================

suppressPackageStartupMessages(library(R.utils))
source("ascend.R"); source("eval_metrics.R"); source("sim_tiers.R")

SEEDS <- 100L + seq_len(as.integer(Sys.getenv("SEEDS", "20")))
GRID  <- expand.grid(n = c(1024, 4096), sp = c(0.5, 0.7, 0.9))
DZ <- 50; DX <- 15; R2 <- 0.5

closure <- function(A) {                          # A[parent, child] -> reach[a, b]
  R <- (A > 0) * 1; P <- R
  for (k in seq_len(nrow(A) - 1)) { P <- ((P %*% A) > 0) * 1; if (!any(P > 0)) break; R <- ((R + P) > 0) * 1 }
  R
}
n_paths <- function(A) {                          # number of directed paths a -> b
  N <- A; P <- A
  for (k in seq_len(nrow(A) - 1)) { P <- P %*% A; if (!any(P > 0)) break; N <- N + P }
  N
}

pair_features <- function(sim) {
  A  <- sim$A_full; V <- rownames(A)
  xs <- unname(sim$labels[grep("^x", names(sim$labels))])      # true names of observed X
  xl <- names(sim$labels)[grep("^x", names(sim$labels))]
  R  <- closure(A); NP <- n_paths(A)
  indeg_x <- colSums(A[grepl("^X", V), , drop = FALSE])[xs]
  indeg_z <- colSums(A[grepl("^Z", V), , drop = FALSE])[xs]
  # R_minus[[v]]: reachability after deleting node v
  R_minus <- lapply(setNames(xs, xs), function(v) { B <- A; B[v, ] <- 0; B[, v] <- 0; closure(B) })
  dat <- as.matrix(sim$dat)
  pc_given_z <- function(a, b) {                  # |partial cor| of a, b given their Z parents
    zp <- intersect(names(which(A[, xs[a]] == 1 | A[, xs[b]] == 1)), grep("^Z", V, value = TRUE))
    zc <- names(sim$labels)[match(zp, sim$labels)]; zc <- zc[!is.na(zc)]
    ya <- dat[, xl[a]]; yb <- dat[, xl[b]]
    if (length(zc)) { Zm <- dat[, zc, drop = FALSE]; ya <- resid(lm(ya ~ Zm)); yb <- resid(lm(yb ~ Zm)) }
    abs(cor(ya, yb))
  }
  d <- length(xs); out <- list()
  for (a in 1:(d - 1)) for (b in (a + 1):d) {
    va <- xs[a]; vb <- xs[b]
    anc_ab <- R[va, vb] == 1; anc_ba <- R[vb, va] == 1
    # confounding, two strengths: (i) the pair shares a direct parent;
    # (ii) some foreground c reaches a while avoiding b and b while avoiding a
    # (background common ancestors are near-universal in this simulator, so
    # they are not informative)
    shared_parent <- any(A[, va] == 1 & A[, vb] == 1)
    conf_x <- any(R_minus[[vb]][, va] == 1 & R_minus[[va]][, vb] == 1 &
                    !(V %in% c(va, vb)) & grepl("^X", V))
    child <- if (anc_ab) b else if (anc_ba) a else NA
    out[[length(out) + 1]] <- data.frame(
      a = xl[a], b = xl[b],
      adjacent = A[va, vb] == 1 || A[vb, va] == 1,
      n_paths = if (anc_ab) NP[va, vb] else if (anc_ba) NP[vb, va] else 0,
      shared_parent = shared_parent, confounded_x = conf_x,
      child_indeg_x = if (is.na(child)) NA else indeg_x[[child]],
      child_indeg_z = if (is.na(child)) NA else indeg_z[[child]],
      pcor = pc_given_z(a, b))
  }
  do.call(rbind, out)
}

rows <- list()
for (g in seq_len(nrow(GRID))) for (seed in SEEDS) {
  sim <- sim_tiers(n = GRID$n[g], d_z = DZ, d_x = DX, r2 = R2, sp = GRID$sp[g], seed = seed)
  M <- tryCatch(withTimeout(ascend(sim, verbose = FALSE), timeout = 600, onTimeout = "error"),
                error = function(e) NULL)
  if (is.null(M)) next
  ps <- pair_states(M, sim$truth)
  rule <- attr(M, "rule"); votes <- attr(M, "rule_votes"); sweep <- attr(M, "rule_sweep")
  ps$rule  <- rule[cbind(ps$a, ps$b)]
  ps$votes <- votes[cbind(ps$a, ps$b)]
  ps$sweep <- sweep[cbind(ps$a, ps$b)]
  ps <- merge(ps, pair_features(sim), by = c("a", "b"))
  rows[[length(rows) + 1]] <- cbind(n = GRID$n[g], sp = GRID$sp[g], seed = seed, ps)
  cat(sprintf("n=%d sp=%.1f seed=%d done\n", GRID$n[g], GRID$sp[g], seed))
}
P <- do.call(rbind, rows)
P$rule[is.na(P$rule)] <- "unresolved"
P$related   <- P$truth != "none"
P$oriented  <- P$est %in% c("a<b", "b<a")
P$correct   <- P$est == P$truth
dir.create("analysis/out", recursive = TRUE, showWarnings = FALSE)
write.csv(P, "analysis/out/orientation_pairs.csv", row.names = FALSE)

rate <- function(x) if (length(x)) mean(x) else NA_real_

## 1. what each rule commits, and how often it is right
T1 <- do.call(rbind, lapply(split(P, P$rule), function(d) data.frame(
  rule = d$rule[1], pairs = nrow(d), share = nrow(d) / nrow(P),
  correct_state = rate(d$correct),
  oriented = sum(d$oriented),
  orient_correct_when_related = rate(d$correct[d$oriented & d$related]),
  median_votes = median(d$votes, na.rm = TRUE))))
write.csv(T1, "analysis/out/orientation_by_rule.csv", row.names = FALSE)

## 2. truly related pairs: oriented correctly / reversed / called unrelated / undetermined
rel <- P[P$related, ]
rel$outcome <- ifelse(rel$est == rel$truth, "correct",
               ifelse(rel$est %in% c("a<b", "b<a"), "reversed",
               ifelse(rel$est == "none", "missed (called unrelated)",
               ifelse(rel$est == "undirected", "undetermined (undirected)", "undetermined (NA)"))))
T2 <- as.data.frame(prop.table(table(sp = rel$sp, n = rel$n, outcome = rel$outcome), c(1, 2)))
write.csv(T2, "analysis/out/orientation_outcomes.csv", row.names = FALSE)

## 3. orientation accuracy among related & oriented pairs, by structural feature
ro <- rel[rel$oriented, ]
by_feat <- function(f, lab) {
  sp <- split(ro$correct, f)
  data.frame(feature = lab, level = names(sp), pairs = lengths(sp),
             orient_acc = vapply(sp, rate, numeric(1)))
}
T3 <- rbind(
  by_feat(cut(ro$child_indeg_x, c(-Inf, 1, 2, 3, Inf), c("1", "2", "3", "4+")), "X-parents of descendant"),
  by_feat(cut(ro$child_indeg_z, c(-Inf, 0, 2, 5, Inf), c("0", "1-2", "3-5", "6+")), "Z-parents of descendant"),
  by_feat(ifelse(ro$shared_parent, "yes", "no"), "shares a direct parent"),
  by_feat(ifelse(ro$confounded_x, "yes", "no"), "foreground common ancestor"),
  by_feat(cut(ro$n_paths, c(0, 1, 2, Inf), c("1", "2", "3+")), "directed paths"),
  by_feat(ifelse(ro$adjacent, "direct edge", "ancestral only"), "adjacency"),
  by_feat(cut(ro$pcor, c(0, 0.05, 0.1, 0.2, 1), c("<0.05", "0.05-0.1", "0.1-0.2", ">0.2"), include.lowest = TRUE),
          "|partial cor| given Z parents"),
  by_feat(ifelse(is.na(ro$rule), "?", ro$rule), "committing rule"))
write.csv(T3, "analysis/out/orientation_by_feature.csv", row.names = FALSE)

options(width = 200)
cat("\n== 1. Rules ==\n");                 print(T1, digits = 3, row.names = FALSE)
cat("\n== 2. Outcomes for truly related pairs (all cells pooled) ==\n")
print(round(prop.table(table(rel$outcome)), 3))
cat("\n== 3. Orientation accuracy by feature ==\n"); print(T3, digits = 3, row.names = FALSE)
