## ======================================================================
## tests/test_ascend_regression.R
## 1. The revised R/ascend.R returns exactly the same ancestrality matrix
##    as the submitted ascend.R (tests/fixtures/ascend_submitted.R) on the
##    submitted simulator, so none of the revision instrumentation changes
##    the algorithm.
## 2. The new attributes are present and consistent: work counters
##    ("stats"), pairwise p-values ("pair_p") and rule provenance ("rule",
##    "rule_sweep", "rule_votes": every resolved pair has a rule).
## 3. sim_dat(z_scale = TRUE) keeps every background variable at unit
##    variance; z_scale = FALSE reproduces the submitted data exactly.
## Reviewer comments: 5, 6, 11 (the attributes), and the simulator fix.
## Run: Rscript tests/test_ascend_regression.R     (base R only, ~30 s)
## ======================================================================
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run the tests from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
source(file.path(ROOT, "R", "ascend.R"))
old <- new.env(); sys.source(file.path(ROOT, "tests", "fixtures", "ascend_submitted.R"), old)
new_attrs <- c("stats", "pair_p", "rule", "rule_sweep", "rule_votes")

for (cfg in list(c(1000, 20, 8, .7), c(500, 30, 12, .5), c(2000, 10, 15, .9))) {
  args <- list(n = cfg[1], d_z = cfg[2], d_x = cfg[3], r2 = .6, sp = cfg[4],
               p_cross = .15, x_effect = 1, seed = 42)
  sim_old <- do.call(old$sim_dat, args)
  sim_new <- do.call(sim_dat, args)                       # z_scale = FALSE by default
  stopifnot(identical(sim_old$dat, sim_new$dat), identical(sim_old$adj_xx, sim_new$adj_xx))

  M0 <- old$ascend(sim_old, verbose = FALSE)
  M1 <- ascend(sim_old, verbose = FALSE)
  M1s <- M1; for (a in new_attrs) attr(M1s, a) <- NULL
  stopifnot(identical(M0, M1s))

  st <- attr(M1, "stats")
  stopifnot(st$n_ci_tests == st$n_pair_tests + st$n_witness_tests + st$n_mb_tests,
            st$n_ci_tests > 0, st$ci_time_sec >= 0, st$n_sweeps >= 1)
  P <- attr(M1, "pair_p"); stopifnot(all(dim(P) == dim(M1)), all(P[!is.na(P)] >= 0 & P[!is.na(P)] <= 1))
  R <- attr(M1, "rule"); off <- row(M1) != col(M1)
  stopifnot(all(!is.na(R[off & !is.na(M1)])),              # every resolved pair has a rule
            all(R[!is.na(R)] %in% c("R1", "R2", "R3", "transitivity", "symmetry", "R3_final")))
  cat(sprintf("n=%d d_z=%d d_x=%d: identical to submitted; %d CI tests (%d pair, %d witness, %d blanket)\n",
              cfg[1], cfg[2], cfg[3], st$n_ci_tests, st$n_pair_tests, st$n_witness_tests, st$n_mb_tests))
}

s <- sim_dat(n = 2000, d_z = 60, d_x = 10, r2 = .5, sp = .7, p_cross = .1, x_effect = .9,
             z_scale = TRUE, seed = 1)
zv <- apply(s$dat[, grep("^z", names(s$dat))], 2, var)
stopifnot(all(abs(zv - 1) < 1e-8))
s0 <- sim_dat(n = 2000, d_z = 60, d_x = 10, r2 = .5, sp = .7, p_cross = .1, x_effect = .9,
              z_scale = FALSE, seed = 1)
cat(sprintf("z_scale=TRUE: var(Z) all 1; z_scale=FALSE: max var(Z) = %.3g\n",
            max(apply(s0$dat[, grep("^z", names(s0$dat))], 2, var))))
cat("ascend regression tests OK\n")
