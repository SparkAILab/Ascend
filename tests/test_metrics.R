## ======================================================================
## tests/test_metrics.R
## Unit tests for R/eval_metrics.R and R/stats_utils.R on hand-built graphs
## whose correct scores are known (reviewer comments 2, 9, 11, minor 9).
## Run: Rscript tests/test_metrics.R      (base R only, < 5 s)
## ======================================================================
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run the tests from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
source(file.path(ROOT, "R", "eval_metrics.R"))
n <- c("x1","x2","x3")
T <- matrix(0,3,3,dimnames=list(n,n)); T["x1","x2"] <- 1; T["x2","x3"] <- 1; T["x1","x3"] <- 1  # chain x1->x2->x3
# perfect estimate
m <- eval_ancestral(T, T); stopifnot(m$dir_precision==1, m$dir_recall==1, m$orient_acc==1, m$shd==0, m$coverage==1, m$ad_acc==1)
# reverse x1-x2, NA on x1-x3, x2-x3 correct
E <- matrix(0,3,3,dimnames=list(n,n)); E["x2","x1"] <- 1; E["x2","x3"] <- 1; E["x1","x3"] <- NA; E["x3","x1"] <- NA
m <- eval_ancestral(E, T)
stopifnot(m$dir_tp==1, m$dir_claims==2, m$n_reversed==1, abs(m$dir_recall-1/3)<1e-12,
          m$orient_acc==0.5, abs(m$coverage-2/3)<1e-12, m$shd==2, m$shd_resolved==1)
# undirected (both directions 1) on x1-x2
E2 <- T; E2["x2","x1"] <- 1
m <- eval_ancestral(E2, T); stopifnot(m$orient_undet==1/3, m$n_undirected==1, m$shd==1, m$dir_claims==4, m$dir_tp==3)
# name order independence
p <- c(3,1,2); m2 <- eval_ancestral(E2[p,p], T); stopifnot(identical(m, m2))
# 0.5 counted by default, not under strict
E3 <- T; E3["x1","x2"] <- 0.5
stopifnot(eval_ancestral(E3,T)$dir_tp==3, eval_ancestral(E3,T,pos_values=1)$dir_tp==2)
cat("eval_metrics tests OK\n")

source(file.path(ROOT, "R", "stats_utils.R"))
set.seed(2); x <- rnorm(30, .1); y <- x - 0.05 + rnorm(30, 0, .02)
pc <- paired_compare(x, y); stopifnot(pc$ci_lo < pc$mean_diff, pc$mean_diff < pc$ci_hi, pc$n_pairs==30)
dt <- data.frame(method=rep(c("ASCEND","G","W"),each=20), cell=rep(1:2,30), rep=rep(1:20,3), f1=runif(60), aupr=runif(60))
pt <- paired_table(dt, c("f1","aupr"), "cell"); stopifnot(nrow(pt)==8, all(pt$family_size==4))
print(head(pt[,c("cell","metric","competitor","mean_diff","ci_lo","ci_hi","p_wilcox","q_bh","family")],3))
st <- summary_table(dt, c("f1"), "cell"); stopifnot(nrow(st)==6)
cat("stats_utils tests OK\n")
