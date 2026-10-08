## ======================================================================
## tests/test_grn_multi.R
## Tests for the multi-method benchmark (benchmarks/grn_multi) and for the
## orientation conventions of every comparator.
##
## Reviewer comments: 1, 2, 11 (the numbers behind them)
## What it checks
##   A. scoring.R: AUPR/AUROC with ties, top-k in expectation, orientation
##      of directed scores, the causal-method path through eval_ancestral()
##   B. graph_to_ancestral(): paths through background variables, undirected edges
##   C. orientation conventions on a simulated chain x1 -> x2 -> x3
##      (skipped per method when its package is not installed):
##        pcalg::pc via as(fit, "amat") is transposed -> t() needed
##        pcalg::ges essgraph matrix, tpc graph, bnlearn amat are [from, to]
##        cbl_fn() is [descendant, ancestor] -> cbl_as_ancestral()
## How to run: Rscript tests/test_grn_multi.R
## ======================================================================

ROOT <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  d
})
setwd(ROOT)
source("R/eval_metrics.R")
source("benchmarks/grn_multi/R/common.R")
source("benchmarks/grn_multi/R/scoring.R")
has <- function(p) requireNamespace(p, quietly = TRUE)
ok <- function(cond, msg) { if (!isTRUE(cond)) stop("FAILED: ", msg); cat("ok  ", msg, "\n") }

## --- A. scoring ---------------------------------------------------------------
ok(abs(auroc_ties(c(3, 2, 1, 0), c(1, 1, 0, 0)) - 1) < 1e-12, "AUROC perfect ranking = 1")
ok(abs(auroc_ties(c(1, 1, 1, 1), c(1, 0, 1, 0)) - 0.5) < 1e-12, "AUROC all tied = 0.5")
ok(abs(expected_tp_topk(c(2, 1, 1, 1), c(0, 1, 1, 0), 2) - 2 / 3) < 1e-12,
   "top-k splits ties in expectation")
if (has("PRROC")) {
  ok(abs(aupr_dg(c(1, 1, 1, 1), c(1, 0, 0, 0)) - 0.25) < 1e-12, "AUPR of a constant score = density")
  ok(aupr_dg(c(4, 3, 2, 1), c(1, 1, 0, 0)) > 0.99, "AUPR perfect ranking ~ 1")
}

nm <- paste0("x", 1:4)
Tr <- matrix(0, 4, 4, dimnames = list(nm, nm)); Tr["x1", "x2"] <- Tr["x2", "x3"] <- 1
S  <- matrix(0, 4, 4, dimnames = list(nm, nm))
S["x1", "x2"] <- 0.9; S["x2", "x1"] <- 0.1          # right way round
S["x3", "x2"] <- 0.8; S["x2", "x3"] <- 0.2          # wrong way round
E <- topk_oriented(S, 2)
ok(E["x1", "x2"] == 1 && E["x3", "x2"] == 1 && sum(E) == 2, "top-k keeps 2 pairs, oriented to the larger score")
if (has("PRROC")) {
  r <- score_one(list(class = "directed", score = S), Tr, "direct", K_asc = 2)
  ok(r$orient_acc_rank == 0.5, "orient_acc_rank: one of two true edges ranked the right way")
  ok(r$f1_Kasc == 1 && r$dir_f1_Kasc == 0.5, "matched-K skeleton F1 = 1, directed F1 = 0.5")
  anc <- Tr; diag(anc) <- NA
  rc <- score_one(list(class = "causal", anc = anc, adj = Tr), Tr, "anc", K_asc = 2)
  ok(rc$dir_f1 == 1 && rc$orient_acc == 1 && rc$shd == 0 && rc$aupr > 0.99, "a perfect causal estimate scores perfectly")
  rr <- score_one(list(class = "causal", anc = t(Tr), adj = t(Tr)), Tr, "anc", K_asc = 2)
  ok(rr$dir_f1 == 0 && rr$orient_acc == 0 && rr$skel_f1 == 1, "a fully reversed estimate: skeleton right, direction wrong")
}

## --- B. graph_to_ancestral -------------------------------------------------------
v <- c("z1", "x1", "x2", "x3")
G <- matrix(0, 4, 4, dimnames = list(v, v))
G["x1", "z1"] <- 1; G["z1", "x2"] <- 1                  # x1 -> z1 -> x2 (only possible without tiers)
G["x2", "x3"] <- G["x3", "x2"] <- 1                     # x2 - x3 undirected
A <- graph_to_ancestral(G, c("x1", "x2", "x3"))
ok(A["x1", "x2"] == 1 && A["x2", "x1"] == 0, "a directed path through a background node counts")
ok(A["x2", "x3"] == 1 && A["x3", "x2"] == 1, "an undirected edge is related in both directions")
ok(transitive_reduction(transitive_closure(Tr))["x1", "x3"] == 0, "transitive reduction drops the implied edge")

## --- C. orientation conventions on a chain -----------------------------------------
set.seed(7); n <- 3000
z <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, paste0("z", 1:3)))
x1 <- z[, 1] + rnorm(n); x2 <- 0.8 * x1 + z[, 2] + rnorm(n); x3 <- 0.8 * x2 + z[, 3] + rnorm(n)
X <- cbind(z, x1 = x1, x2 = x2, x3 = x3)
if (has("pcalg")) {
  suppressPackageStartupMessages(library(pcalg))
  fit <- pc(list(C = cor(X), n = n), gaussCItest, labels = colnames(X), alpha = 0.01)
  am <- as(fit, "amat")
  ok(am["x2", "x1"] == 1 && am["x1", "x2"] == 0, "pcalg as(fit, 'amat') is transposed ([to, from])")
  gm <- as(fit@graph, "matrix")
  ok(gm["x1", "x2"] == 1 && gm["x2", "x1"] == 0, "pcalg fit@graph is [from, to]")
  es <- as(ges(new("GaussL0penObsScore", X))$essgraph, "matrix"); dimnames(es) <- list(colnames(X), colnames(X))
  ok(es["x1", "x2"] && !es["x2", "x1"], "pcalg GES essgraph matrix is [from, to]")
}
if (has("tpc")) {
  ft <- tpc::tpc(list(C = cor(X), n = n), pcalg::gaussCItest, alpha = 0.01, labels = colnames(X),
                 tiers = c(1, 1, 1, 2, 2, 2))
  tg <- as(ft@graph, "matrix")
  ok(tg["x1", "x2"] == 1 && tg["x2", "x1"] == 0, "tpc graph is [from, to]")
}
if (has("bnlearn")) {
  fb <- bnlearn::hc(as.data.frame(X), blacklist = bnlearn::tiers2blacklist(list(colnames(z), c("x1", "x2", "x3"))))
  ab <- bnlearn::amat(fb)
  ok(ab["x1", "x2"] == 1 && ab["x2", "x1"] == 0, "bnlearn amat is [from, to]")
}
if (all(sapply(c("data.table", "dplyr", "foreach", "glmnet", "lightgbm", "pcalg", "doMC"), has))) {
  suppressPackageStartupMessages({ library(data.table); library(dplyr); library(foreach) })
  e <- new.env(); assign("ROOT", ROOT, envir = e)
  suppressMessages(sys.source(file.path(ROOT, "R", "ascend.R"), envir = e))
  suppressMessages(sys.source(file.path(ROOT, "R", "cbl.R"), envir = e))
  sim <- e$sim_dat(2000, 20, 6, r2 = .7, sp = .6, p_cross = .2, z_scale = TRUE, seed = 5)
  At <- t(sim$adj_xx); At[is.na(At)] <- 0
  truth <- transitive_closure(At)
  m <- e$cbl_fn(sim, B = 50)
  good <- eval_ancestral(e$cbl_as_ancestral(m), truth, 1)$orient_acc
  bad  <- eval_ancestral(m, truth, 1)$orient_acc
  ok(good > bad, sprintf("CBL output must be transposed (orientation %.2f vs %.2f raw)", good, bad))
}
cat("grn_multi tests OK\n")
