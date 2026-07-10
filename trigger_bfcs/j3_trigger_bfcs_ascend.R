# ==========================================================================
# FAIR COMPARISON (3-way): ASCEND vs TRIGGER vs BFCS on curated YEASTRACT truth
# ==========================================================================
# This is a corrected redo of the ASCEND-vs-TRIGGER harness. Four issues in the
# previous version each happened to penalise ASCEND specifically; all four are
# fixed here, and BFCS (Bucur et al., PGM 2018) is added as a genuine third
# method computed from the SAME yeast data (no fabricated numbers).
#
# WHAT CHANGED, AND WHY EACH ONE MATTERED
# ---------------------------------------------------------------------------
# FIX 1  Edge definition.  Old code: ascend_binary <- (M==0.5 | M==1).
#        M==0.5 means "not a descendant" -- a weak, high-frequency, largely
#        uninformative claim, NOT a directed regulatory edge. Counting it as a
#        predicted edge floods ASCEND's prediction set with near-meaningless
#        pairs that almost never match curated truth, which tanks ASCEND's
#        precision. Now: a predicted ASCEND edge is M==1 ONLY (ancestor).
#
# FIX 2  Depth-symmetric ancestral truth.  ASCEND's M==1 is transitively closed
#        with NO depth cap (a 5-hop chain is written as 1). The old "ancestral"
#        panel capped truth at <=2 hops, so a correct long-range ASCEND ancestor
#        was scored as a false positive against a truth that could not credit it.
#        TRIGGER/BFCS never hit this (each is a single pairwise claim). Fix: the
#        ancestral panel's truth is now the FULL transitive closure of curated
#        edges, restricted to paths *within g* (the induced subgraph). That is
#        exactly the set of ancestral relations ASCEND could possibly recover
#        given only genes in g, so both estimands are depth-consistent, and
#        confining paths to g keeps the closure from exploding to "everything".
#
# FIX 3  Matched-K is computed and spent on the SAME region that is scored.
#        Scoring only ever looks at rows of defined_regulators (~20-40 genes),
#        but old K summed ascend_binary over the whole 150x150 matrix, and
#        competitors picked their top-K from anywhere. Now K counts only
#        regulator-sourced ASCEND edges, and every competitor's top-K is drawn
#        from the same regulator-sourced candidate pool. K is also panel-specific
#        (each estimand matches ASCEND's edge count for that estimand).
#
# FIX 4  Ground-truth density is reported.  Building g from curated TFs + their
#        curated targets makes the induced graph denser than a typical GRN, and
#        ASCEND's advantage is documented to widen as the graph gets sparser
#        (paper Fig. 2 / 7f). This is not a bug, but the reader must be able to
#        judge it, so we print the density of the induced truth alongside results.
#
# BFCS   Implemented from Bucur et al. (2018), eq. (6), DMAG + background-
#        knowledge prior (Table 1), scoring the causal chain L -> Ti -> Tj as the
#        posterior of model M6 (Tj _||_ L | Ti). Instruments are the strongest-
#        linked markers to Ti (cis-linkage in spirit; the causality-equivalence
#        theorem only needs L linked to Ti). Validated on synthetic chains /
#        colliders / confounders before use. This is a real run on the same data.
# ==========================================================================



#Attached Ground Truth data was downloaded from https://yeastract-plus.org/yeastract/scerevisiae/formregmatrix.php

suppressMessages({
  library(trigger)
  library(org.Sc.sgd.db)
  library(AnnotationDbi)
  library(ggplot2)
})

# ======================= EDIT THESE PATHS ================================
BFCS_DIR     <- "~/Phd/gen2-r/BFCS"   # holds the .RData with the TRIGGER matrix
ASCEND_SRC   <- "~/Phd/gen2-r/ascend.R"
GROUND_TRUTH <- "~/Phd/gen2-r/June29_last/Ascend/trigger_bfcs/yeasttract999.tsv"
OUT_DIR      <- "~/Phd/gen2-r/June29_last/Ascend/trigger_bfcs/fair_out_3way"
# knobs (not tuned to any method's output)
N_TF_SEED      <- 40    # seed TFs to attempt
MAX_GENES      <- 150   # cap on gene-set size
MIN_TARGETS    <- 3     # a TF must have >= this many measurable targets to seed
MIN_REGULATORS <- 5     # abort scoring if fewer defined regulators than this
N_INSTRUMENTS  <- 5     # BFCS: markers per source gene to maximise the chain over
ASCEND_DIRECT_MODE <- "reduction"  # "reduction" (fair) or "closure" (strict)
#   ASCEND returns a transitively-closed ANCESTRAL matrix, not direct edges.
#   The ANCESTRAL panel is the like-for-like estimand for ASCEND (and, since
#   TRIGGER/BFCS also estimate causal reachability rather than direct binding,
#   for all three methods) -- report it as primary.
#   On the DIRECT (exact-edge) panel, "reduction" scores ASCEND's transitive
#   reduction, i.e. its minimal direct-cause edges, which is the like-for-like
#   direct comparison. "closure" scores raw M==1 and is deliberately strict:
#   it penalises every indirect-but-correct ancestor as a false positive.
# ==========================================================================

dir.create(path.expand(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
source(path.expand(ASCEND_SRC))
set.seed(1)

## =========================================================================
## BFCS: Bayes Factors of Covariance Structures  (Bucur et al., PGM 2018)
## p(M6 | D) = p(X1 -> X2 -> X3 | D), the causal chain through the mediator X2.
## X1 = instrument L, X2 = Ti (mediator/regulator), X3 = Tj (target).
## Correlations: r12 = cor(L,Ti), r13 = cor(L,Tj), r23 = cor(Ti,Tj).
## Vectorised over the target index j (r13, r23 are vectors; r12 is scalar).
## =========================================================================
# DMAG-with-BK model counts from Table 1 (M0..M10); M2 has prior 0 under BK.
BFCS_PRIOR <- c(3, 2, 0, 2, 1, 1, 1, 3, 1, 1, 1)

bfcs_chain_prob_vec <- function(r12, r13, r23, n, nu = 4, prior = BFCS_PRIOR) {
  eps <- 1e-12
  a12 <- max(1 - r12^2, eps)
  a13 <- pmax(1 - r13^2, eps)
  a23 <- pmax(1 - r23^2, eps)
  D   <- pmax(1 - r12^2 - r13^2 - r23^2 + 2 * r12 * r13 * r23, eps)
  e1  <- (n + nu) / 2; e2 <- (n + nu - 1) / 2
  lf  <- log((n + nu - 2) / (nu - 2))
  lg  <- lgamma((n + nu) / 2) + lgamma((nu - 1) / 2) -
    lgamma((n + nu - 1) / 2) - lgamma(nu / 2)
  lD  <- log(D); l12 <- log(a12); l13 <- log(a13); l23 <- log(a23)
  logB <- cbind(
    M0  = rep(0, length(r13)),                  # full (reference)
    M1  = lf - lg + e2 * l12,                    # X1 _||_ X2         (marginal)
    M2  = lf - lg + e2 * l23,                    # X2 _||_ X3
    M3  = lf - lg + e2 * l13,                    # X3 _||_ X1
    M4  = lg + e1 * (lD - l13 - l23),            # X1 _||_ X2 | X3
    M5  = lg + e1 * (lD - l12 - l13),            # X2 _||_ X3 | X1
    M6  = lg + e1 * (lD - l12 - l23),            # X3 _||_ X1 | X2  <- causal chain
    M7  = lf + e1 * (lD - l23),                  # X1 _||_ (X2,X3)
    M8  = lf + e1 * (lD - l13),                  # X2 _||_ (X3,X1)
    M9  = lf + e1 * (lD - l12),                  # X3 _||_ (X1,X2)
    M10 = lf + lg + e1 * lD                      # empty
  )
  lp <- log(prior / sum(prior))                  # prior 0 -> -Inf -> weight 0
  W  <- sweep(logB, 2, lp, `+`)
  W  <- W - apply(W, 1, function(v) max(v[is.finite(v)]))
  P  <- exp(W); P[!is.finite(P)] <- 0
  P  <- P / rowSums(P)
  P[, "M6"]
}

# Full BFCS matrix over g: B[i,j] = P(gene i -> gene j).
build_bfcs_matrix <- function(z, x, g, n, n_instruments = N_INSTRUMENTS, nu = 4) {
  xg  <- x[, g, drop = FALSE]
  CzX <- suppressWarnings(cor(z, xg)); CzX[is.na(CzX)] <- 0   # d_z x |g|
  Cxx <- suppressWarnings(cor(xg));    Cxx[is.na(Cxx)] <- 0   # |g| x |g|
  m   <- length(g)
  B   <- matrix(0, m, m, dimnames = list(g, g))
  for (i in seq_len(m)) {
    inst <- order(CzX[, i]^2, decreasing = TRUE)[seq_len(min(n_instruments, nrow(CzX)))]
    best <- numeric(m)
    for (k in inst) {
      r12 <- CzX[k, i]
      if (r12^2 < 1e-8) next                     # instrument not linked to Ti -> no chain
      pv  <- bfcs_chain_prob_vec(r12, CzX[k, ], Cxx[i, ], n, nu)
      best <- pmax(best, pv, na.rm = TRUE)
    }
    best[i] <- 0
    B[i, ]  <- best
  }
  B
}

## graph helpers ------------------------------------------------------------
close_reach_mat <- function(A) {                 # boolean transitive closure (NA-safe)
  R <- (A > 0); R[is.na(R)] <- FALSE; R <- R * 1
  repeat { R2 <- ((R + R %*% R) > 0) * 1; if (all(R2 == R)) break; R <- R2 }
  diag(R) <- 0; R > 0
}
transitive_reduction <- function(Rlogical) {     # minimal edges whose closure = R
  R <- (Rlogical > 0) * 1; diag(R) <- 0
  red <- (R > 0) & !((R %*% R) > 0)
  diag(red) <- FALSE; red
}
# FIX 2: ancestral truth = full closure of curated edges, confined to paths in g.
reachable_closure_within <- function(gt_regulons, nodes) {
  m <- length(nodes)
  A <- matrix(0, m, m, dimnames = list(nodes, nodes))
  for (r in nodes) { tt <- intersect(gt_regulons[[r]], nodes); if (length(tt)) A[r, tt] <- 1 }
  R <- close_reach_mat(A)
  setNames(lapply(nodes, function(r) nodes[R[r, ]]), nodes)
}

## -------------------------------------------------------------------------
## STAGE 0a - load + normalise ground truth (YEASTRACT ";"-sep, no header)
## -------------------------------------------------------------------------
gt_raw <- read.delim(path.expand(GROUND_TRUTH), sep = ";", header = FALSE,
                     stringsAsFactors = FALSE, check.names = FALSE,
                     quote = "", strip.white = TRUE)[, 1:2]
if (length(unique(gt_raw[[1]])) <= length(unique(gt_raw[[2]]))) {
  reg_v <- gt_raw[[1]]; tgt_v <- gt_raw[[2]]
} else { reg_v <- gt_raw[[2]]; tgt_v <- gt_raw[[1]] }

build_name2orf <- function() {
  valid_orf <- keys(org.Sc.sgd.db, keytype = "ORF"); dict <- character(0)
  add_col <- function(col) {
    df <- tryCatch(suppressMessages(AnnotationDbi::select(
      org.Sc.sgd.db, keys = valid_orf, keytype = "ORF", columns = col)),
      error = function(e) NULL)
    if (is.null(df) || !col %in% colnames(df)) return(invisible())
    nm <- toupper(trimws(df[[col]])); orf <- df[["ORF"]]
    ok <- !is.na(nm) & nm != "" & !is.na(orf)
    d <- orf[ok]; names(d) <- nm[ok]; d <- d[!duplicated(names(d))]
    dict <<- c(dict, d[setdiff(names(d), names(dict))])
  }
  for (col in c("GENENAME", "COMMON", "ALIAS")) add_col(col)
  self <- valid_orf; names(self) <- toupper(valid_orf)
  c(dict, self[setdiff(names(self), names(dict))])
}
name2orf <- build_name2orf()
to_orf <- function(v) unname(name2orf[toupper(trimws(v))])

gt <- data.frame(regulator = to_orf(reg_v), target = to_orf(tgt_v),
                 stringsAsFactors = FALSE)
gt <- unique(gt[!is.na(gt$regulator) & !is.na(gt$target) & gt$regulator != gt$target, ])
gt_regulons <- split(gt$target, gt$regulator)
curated_tfs <- names(gt_regulons)
cat(sprintf("Curated: %d direct pairs, %d regulators\n", nrow(gt), length(curated_tfs)))

## -------------------------------------------------------------------------
## STAGE 0b - data + TRIGGER matrix; define the SHARED method gene-space
## -------------------------------------------------------------------------
data(yeast)
z <- t(yeast$marker - 1)
x <- t(yeast$exp)
n_samples <- nrow(x)
load(file.path(path.expand(BFCS_DIR), "inst/extdata/yeast_regulatory_probabilities.RData"))
stopifnot(exists("yeast_trigger_w50k"))
method_space <- intersect(colnames(x), rownames(yeast_trigger_w50k))
cat(sprintf("Samples n = %d | shared gene-space (x AND trigger matrix): %d genes\n",
            n_samples, length(method_space)))
cat(sprintf("Curated TFs inside shared space: %d\n",
            length(intersect(curated_tfs, method_space))))

## -------------------------------------------------------------------------
## STAGE 0b' - DIAGNOSTIC: why do so few curated TFs land in the method space?
## If regulators-in-space is tiny but targets-in-space is large, the regulator
## column is a different identifier type (or the orientation is flipped). If the
## RAW tokens match method_space better than the mapped ORFs, then method_space
## is NOT ORF-based and to_orf() is mapping you off it -- match on raw instead.
## -------------------------------------------------------------------------
## -------------------------------------------------------------------------
## STAGE 0b - CORRECTED: data + TRIGGER matrix; map space to systematic ORFs
## -------------------------------------------------------------------------
data(yeast)
z <- t(yeast$marker - 1)
x <- t(yeast$exp)
n_samples <- nrow(x)

load(file.path(path.expand(BFCS_DIR), "inst/extdata/yeast_regulatory_probabilities.RData"))
stopifnot(exists("yeast_trigger_w50k"))

# 1. Identify raw shared space names
raw_method_space <- intersect(colnames(x), rownames(yeast_trigger_w50k))

# 2. Map raw names to ORFs using your existing name2orf dictionary
mapped_orfs <- name2orf[toupper(trimws(raw_method_space))]
valid_idx <- !is.na(mapped_orfs)

# 3. Rename columns of x and rows/cols of trigger matrix to clean ORFs
# (This filters down to elements that successfully map to a yeast ORF)
method_space <- unname(mapped_orfs[valid_idx])
raw_names_mapped <- raw_method_space[valid_idx]

# Subset and update names to ORFs
x <- x[, raw_names_mapped, drop = FALSE]
colnames(x) <- name2orf[toupper(trimws(colnames(x)))]

yeast_trigger_w50k <- yeast_trigger_w50k[raw_names_mapped, raw_names_mapped]
rownames(yeast_trigger_w50k) <- colnames(yeast_trigger_w50k) <- colnames(x)

# Re-verify method space is now purely ORF-based
method_space <- colnames(x)

cat(sprintf("Samples n = %d | shared gene-space successfully mapped to ORFs: %d genes\n",
            n_samples, length(method_space)))
cat(sprintf("Curated TFs inside shared space: %d\n",
            length(intersect(curated_tfs, method_space))))

## -------------------------------------------------------------------------
## STAGE 0c - DIAGNOSTIC: how many TFs can this dataset actually support?
## -------------------------------------------------------------------------
tf_usable <- intersect(curated_tfs, method_space)
tf_ntar   <- sapply(tf_usable, function(tf) length(intersect(gt_regulons[[tf]], method_space)))
tf_usable <- tf_usable[tf_ntar >= MIN_TARGETS]
cat(sprintf("Usable TFs (>=%d measurable targets): %d\n", MIN_TARGETS, length(tf_usable)))
if (length(tf_usable) < MIN_REGULATORS) {
  stop(sprintf(paste0("Only %d usable TFs in the shared gene-space (need >= %d).\n",
                      "This dataset/matrix cannot support a fair TF-level comparison.\n",
                      "-> report the synthetic-data comparison instead, or widen method_space."),
               length(tf_usable), MIN_REGULATORS))
}

## -------------------------------------------------------------------------
## STAGE 0d - build TF-centred gene set within the shared space
## -------------------------------------------------------------------------
tf_usable <- tf_usable[order(-sapply(tf_usable,
                                     function(tf) length(intersect(gt_regulons[[tf]], method_space))))]
seed_tfs     <- head(tf_usable, N_TF_SEED)
seed_targets <- unique(unlist(lapply(seed_tfs,
                                     function(tf) intersect(gt_regulons[[tf]], method_space))))
g <- intersect(unique(c(seed_tfs, seed_targets)), method_space)
if (length(g) > MAX_GENES) {
  tgt_only <- setdiff(g, seed_tfs)
  deg <- sapply(tgt_only, function(gene)
    sum(vapply(seed_tfs, function(tf) gene %in% gt_regulons[[tf]], logical(1))))
  keep <- names(sort(deg, decreasing = TRUE))[seq_len(MAX_GENES - length(seed_tfs))]
  g <- unique(c(seed_tfs, keep))
}
gt_in_g <- gt[gt$regulator %in% g & gt$target %in% g, ]
defined_regulators <- intersect(unique(gt_in_g$regulator), g)
cat(sprintf("Gene set: %d genes | defined regulators: %d | direct edges in g: %d\n",
            length(g), length(defined_regulators), nrow(gt_in_g)))
if (length(defined_regulators) < MIN_REGULATORS)
  stop(sprintf("Only %d defined regulators landed in g (need >= %d). Raise N_TF_SEED/MAX_GENES.",
               length(defined_regulators), MIN_REGULATORS))

## truth sets (both panels), on the SAME regulator set
direct_truth <- setNames(lapply(defined_regulators,
                                function(r) intersect(gt_regulons[[r]], g)), defined_regulators)
reach_all    <- reachable_closure_within(gt_regulons, g)          # FIX 2: within-g closure
reach_truth  <- reach_all[defined_regulators]

## -------------------------------------------------------------------------
## STAGE 0d' - FIX 4: report induced ground-truth density so the reader can judge
## -------------------------------------------------------------------------
poss_pairs_reg <- length(defined_regulators) * (length(g) - 1)
dens_direct    <- nrow(gt_in_g) / poss_pairs_reg
n_reach_edges  <- sum(sapply(reach_truth, length))
dens_reach     <- n_reach_edges / poss_pairs_reg
cat("\n----- ground-truth density in g (FIX 4) -----\n")
cat(sprintf("  direct  : %d edges | avg %.2f targets/reg | density(reg-sourced) %.3f\n",
            nrow(gt_in_g), mean(sapply(direct_truth, length)), dens_direct))
cat(sprintf("  ancestral(within-g closure): %d edges | avg %.2f desc/reg | density %.3f\n",
            n_reach_edges, mean(sapply(reach_truth, length)), dens_reach))
cat("  (Real GRNs are sparse; ASCEND's advantage over other methods widens as the\n")
cat("   true graph gets sparser. A TF+targets construction is comparatively dense,\n")
cat("   which is ASCEND's least-favourable regime -- interpret accordingly.)\n")

trigger_subset <- yeast_trigger_w50k[g, g]

## -------------------------------------------------------------------------
## STAGE 0e - DIAGNOSTICS: does each competitor even have signal on these TFs?
## If ~0, a later score of 0 is a starved harness, not a fair loss.
## -------------------------------------------------------------------------
trig_reg_rows <- trigger_subset[defined_regulators, , drop = FALSE]
cat("\n----- competitor signal on defined regulators (STAGE 0e) -----\n")
cat(sprintf("  TRIGGER: %d nonzero entries in reg rows (max=%.3f); %d/%d regs have >=1 nonzero\n",
            sum(trig_reg_rows > 0, na.rm = TRUE), max(trigger_subset, na.rm = TRUE),
            sum(rowSums(trig_reg_rows > 0, na.rm = TRUE) > 0), length(defined_regulators)))
# BFCS instrument availability on the regulators
CzX_reg <- suppressWarnings(cor(z, x[, defined_regulators, drop = FALSE]))
best_r2 <- apply(CzX_reg^2, 2, max, na.rm = TRUE)
cat(sprintf("  BFCS   : %d/%d regs have an instrument with cor^2 > 0.05 (max cor^2 = %.3f)\n",
            sum(best_r2 > 0.05), length(defined_regulators), max(best_r2)))

## -------------------------------------------------------------------------
## STAGE 1 - run ASCEND
## -------------------------------------------------------------------------
ad <- as.data.frame(cbind(z, x[, g]))
colnames(ad) <- c(paste0("z", seq_len(ncol(z))), paste0("x", seq_along(g)))
for (cc in colnames(ad)) ad[[cc]] <- as.numeric(scale(ad[[cc]]))
ad[is.na(ad)] <- 0
sim_obj <- list(dat = ad, adj_xx = NULL, params = list(d_z = ncol(z), d_x = length(g)))
ascend_network <- ascend(sim_obj, alpha = 0.05, alpha_mb = 0.05,
                         fdr = TRUE, min_votes = 1, verbose = TRUE)
idx <- as.integer(sub("^x", "", rownames(ascend_network)))
rownames(ascend_network) <- colnames(ascend_network) <- g[idx]
ascend_subset <- ascend_network[g, g]

## -------------------------------------------------------------------------
## STAGE 1b - build BFCS matrix on the same data / gene set
## -------------------------------------------------------------------------
cat("\n[bfcs] scoring triplets ...\n")
bfcs_subset <- build_bfcs_matrix(z, x, g, n_samples)
bfcs_reg_rows <- bfcs_subset[defined_regulators, , drop = FALSE]
cat(sprintf("  BFCS: %d nonzero entries in reg rows (max=%.3f); %d/%d regs have >=1 nonzero\n",
            sum(bfcs_reg_rows > 1e-6), max(bfcs_subset),
            sum(rowSums(bfcs_reg_rows > 1e-6) > 0), length(defined_regulators)))

## -------------------------------------------------------------------------
## STAGE 2 - binarise; matched-K restricted to the SCORED region (FIX 1 + FIX 3)
## -------------------------------------------------------------------------
# FIX 1: an ASCEND edge is M==1 only (ancestor). 0.5 ("non-descendant") is NOT an edge.
ascend_closed <- close_reach_mat(ascend_subset == 1)         # idempotent; guards closure
ascend_bin_anc    <- ascend_closed                            # ancestral panel
ascend_bin_direct <- if (ASCEND_DIRECT_MODE == "reduction")   # direct panel
  transitive_reduction(ascend_closed) else ascend_closed

# FIX 3: K counts ONLY regulator-sourced ASCEND edges (the region that is scored),
# and is computed per panel (each estimand matches ASCEND's own edge count).
K_anc    <- sum(ascend_bin_anc[defined_regulators, , drop = FALSE])
K_direct <- sum(ascend_bin_direct[defined_regulators, , drop = FALSE])
cat(sprintf("\nMatched K (regulator-sourced ASCEND edges): direct=%d  ancestral=%d\n",
            K_direct, K_anc))

# competitors: top-K drawn ONLY from regulator-sourced candidate edges (FIX 3)
top_scored_edges_mat <- function(M, K, allowed_rows) {
  out <- matrix(FALSE, nrow(M), ncol(M), dimnames = dimnames(M))
  if (K <= 0) return(out)
  M2 <- M; diag(M2) <- NA
  rowmask <- rownames(M2) %in% allowed_rows              # length-nrow, recycles over cols
  cand <- which(is.finite(M2) & M2 > 0 & rowmask)
  if (!length(cand)) return(out)
  ord <- cand[order(M2[cand], decreasing = TRUE)][seq_len(min(K, length(cand)))]
  out[ord] <- TRUE; out
}
trigger_bin_direct <- top_scored_edges_mat(trigger_subset, K_direct, defined_regulators)
trigger_bin_anc    <- top_scored_edges_mat(trigger_subset, K_anc,    defined_regulators)
bfcs_bin_direct    <- top_scored_edges_mat(bfcs_subset,    K_direct, defined_regulators)
bfcs_bin_anc       <- top_scored_edges_mat(bfcs_subset,    K_anc,    defined_regulators)
cat(sprintf("TRIGGER claimed %d/%d (direct) and %d/%d (anc) edges at matched K\n",
            sum(trigger_bin_direct), K_direct, sum(trigger_bin_anc), K_anc))
cat(sprintf("BFCS    claimed %d/%d (direct) and %d/%d (anc) edges at matched K\n",
            sum(bfcs_bin_direct), K_direct, sum(bfcs_bin_anc), K_anc))

get_targets <- function(r, m) if (!r %in% rownames(m)) character(0) else
  setdiff(names(which(m[r, ])), r)

## -------------------------------------------------------------------------
## STAGE 3 - edge-level scoring, both panels, three methods
## -------------------------------------------------------------------------
edge_metrics <- function(binary_mat, method_name, truth_sets, panel) {
  tp <- fp <- fn <- 0
  for (r in defined_regulators) {
    true_t <- truth_sets[[r]]; pred_t <- get_targets(r, binary_mat)
    tp <- tp + length(intersect(pred_t, true_t))
    fp <- fp + length(setdiff(pred_t, true_t))
    fn <- fn + length(setdiff(true_t, pred_t))
  }
  precision <- if (tp + fp > 0) tp / (tp + fp) else NA
  recall    <- if (tp + fn > 0) tp / (tp + fn) else NA
  f1 <- if (!is.na(precision) && !is.na(recall) && precision + recall > 0)
    2 * precision * recall / (precision + recall) else NA
  data.frame(method = method_name, panel = panel, precision = precision,
             recall = recall, f1 = f1, tp = tp, fp = fp, fn = fn)
}
edge_results <- rbind(
  edge_metrics(ascend_bin_direct,  "ASCEND",  direct_truth, "direct"),
  edge_metrics(trigger_bin_direct, "TRIGGER", direct_truth, "direct"),
  edge_metrics(bfcs_bin_direct,    "BFCS",    direct_truth, "direct"),
  edge_metrics(ascend_bin_anc,     "ASCEND",  reach_truth,  "ancestral"),
  edge_metrics(trigger_bin_anc,    "TRIGGER", reach_truth,  "ancestral"),
  edge_metrics(bfcs_bin_anc,       "BFCS",    reach_truth,  "ancestral"))
cat("\n===== EDGE-LEVEL (both panels, three methods) =====\n"); print(edge_results, row.names = FALSE)

## -------------------------------------------------------------------------
## STAGE 4 - per-regulator Jaccard + paired Wilcoxon (all method pairs), BH
## -------------------------------------------------------------------------
jaccard <- function(a, b) if (length(a) == 0 && length(b) == 0) NA else
  length(intersect(a, b)) / length(union(a, b))
regulon_tbl <- function(binary_mat, method_name, truth_sets, panel)
  do.call(rbind, lapply(defined_regulators, function(r) {
    pred_t <- get_targets(r, binary_mat)
    data.frame(method = method_name, panel = panel, regulator = r,
               jaccard = jaccard(pred_t, truth_sets[[r]]),
               n_true = length(truth_sets[[r]]), n_pred = length(pred_t))
  }))
regulon_results <- rbind(
  regulon_tbl(ascend_bin_direct,  "ASCEND",  direct_truth, "direct"),
  regulon_tbl(trigger_bin_direct, "TRIGGER", direct_truth, "direct"),
  regulon_tbl(bfcs_bin_direct,    "BFCS",    direct_truth, "direct"),
  regulon_tbl(ascend_bin_anc,     "ASCEND",  reach_truth,  "ancestral"),
  regulon_tbl(trigger_bin_anc,    "TRIGGER", reach_truth,  "ancestral"),
  regulon_tbl(bfcs_bin_anc,       "BFCS",    reach_truth,  "ancestral"))
cat("\n===== MEAN Jaccard =====\n")
print(aggregate(jaccard ~ method + panel, regulon_results,
                function(v) round(mean(v, na.rm = TRUE), 4)))

pairs_to_test <- list(c("ASCEND", "TRIGGER"), c("ASCEND", "BFCS"), c("TRIGGER", "BFCS"))
sig_tbl <- do.call(rbind, lapply(c("direct", "ancestral"), function(pn)
  do.call(rbind, lapply(pairs_to_test, function(pr) {
    a <- regulon_results$jaccard[regulon_results$method == pr[1] & regulon_results$panel == pn]
    b <- regulon_results$jaccard[regulon_results$method == pr[2] & regulon_results$panel == pn]
    ok <- is.finite(a) & is.finite(b)
    wt <- tryCatch(wilcox.test(a[ok], b[ok], paired = TRUE), error = function(e) NULL)
    data.frame(panel = pn, pair = paste(pr, collapse = " vs "), n_pairs = sum(ok),
               median_diff = if (any(ok)) median(a[ok] - b[ok]) else NA,
               p_raw = if (!is.null(wt)) wt$p.value else NA)
  }))))
sig_tbl$q <- p.adjust(sig_tbl$p_raw, method = "BH")
cat("\n===== Paired Wilcoxon (per-regulator Jaccard), BH =====\n"); print(sig_tbl, row.names = FALSE)
if (any(sig_tbl$n_pairs < 8, na.rm = TRUE))
  cat("NOTE: <8 paired regulators -> Wilcoxon is underpowered; treat p-values as descriptive.\n")

## -------------------------------------------------------------------------
## STAGE 5 - pairwise concordance (ancestral panel, at matched K_anc)
## -------------------------------------------------------------------------
edge_set <- function(m) { w <- which(m, arr.ind = TRUE)
if (!nrow(w)) character(0) else paste(rownames(m)[w[, 1]], colnames(m)[w[, 2]], sep = "->") }
es <- list(ASCEND = edge_set(ascend_bin_anc), TRIGGER = edge_set(trigger_bin_anc),
           BFCS = edge_set(bfcs_bin_anc))
cat("\nConcordance (ancestral panel, matched K):\n")
for (pr in pairs_to_test)
  cat(sprintf("  %-16s %d shared edges\n", paste(pr, collapse = "-"),
              length(intersect(es[[pr[1]]], es[[pr[2]]]))))

## -------------------------------------------------------------------------
## FIGURES
## -------------------------------------------------------------------------
cols <- c(ASCEND = "#0072B2", TRIGGER = "#D55E00", BFCS = "#009E73")
edge_results$method <- factor(edge_results$method, levels = names(cols))
regulon_results$method <- factor(regulon_results$method, levels = names(cols))

p_f1 <- ggplot(edge_results, aes(method, f1, fill = method)) +
  geom_col(alpha = .85) +
  geom_text(aes(label = ifelse(is.na(f1), "NA", sprintf("F1=%.2f\nP=%.2f R=%.2f",
                                                        f1, precision, recall))),
            vjust = -0.2, size = 2.7) +
  facet_wrap(~panel) + scale_fill_manual(values = cols) +
  expand_limits(y = max(edge_results$f1, na.rm = TRUE) * 1.25) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Edge-level F1 by estimand", y = "F1", x = NULL)
p_jac <- ggplot(regulon_results[is.finite(regulon_results$jaccard), ],
                aes(method, jaccard, colour = method)) +
  geom_boxplot(outlier.shape = NA, alpha = .3) +
  geom_jitter(width = .15, size = 1.5, alpha = .6) +
  facet_wrap(~panel) + scale_colour_manual(values = cols) +
  theme_classic() + theme(legend.position = "none") +
  labs(title = "Per-regulator Jaccard by estimand", y = "Jaccard", x = NULL)
ggsave(file.path(path.expand(OUT_DIR), "edge_f1_3way.pdf"),        p_f1, width = 7.5, height = 4)
ggsave(file.path(path.expand(OUT_DIR), "regulon_jaccard_3way.pdf"), p_jac, width = 7.5, height = 4)

save(g, K_direct, K_anc, gt, gt_in_g, defined_regulators, direct_truth, reach_truth,
     dens_direct, dens_reach, ascend_subset, trigger_subset, bfcs_subset,
     ascend_bin_direct, ascend_bin_anc, edge_results, regulon_results, sig_tbl,
     file = file.path(path.expand(OUT_DIR), "fair_comparison_3way.RData"))
cat("\n=== DONE. Read STAGE 0c/0d'/0e diagnostics before trusting any panel. ===\n")
cat("=== The ANCESTRAL panel is the like-for-like estimand for ASCEND; the DIRECT\n")
cat("=== panel is inherently unfavourable to any transitively-closed method. ===\n")






library(patchwork)

combined <- p_f1 + p_jac +
  plot_layout(ncol = 2, widths = c(1,1))

ggsave(
  file.path(path.expand(OUT_DIR), "comparison_3way.pdf"),
  combined,
  width = 10,
  height = 4.5
)