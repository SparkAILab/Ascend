## ======================================================================
## tests/test_sim_tiers.R
## Checks R/sim_tiers.R, the tier-violation simulator (reviewer comment 4):
##   - with no violations, the truth equals the transitive closure of the
##     foreground DAG and background variables have unit variance;
##   - each violation type removes / relabels the expected number of
##     variables and is recorded in $violations;
##   - feedback variables really are caused by a foreground variable;
##   - the truth always refers to the analyst's labels.
## Run: Rscript tests/test_sim_tiers.R      (base R only, < 5 s)
## ======================================================================
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run the tests from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
source(file.path(ROOT, "R", "sim_tiers.R")); source(file.path(ROOT, "R", "eval_metrics.R"))

s <- sim_tiers(n = 3000, d_z = 30, d_x = 12, seed = 7)
stopifnot(all(abs(apply(s$dat[, grep("^z", names(s$dat))], 2, var) - 1) < 1e-8))
stopifnot(identical(unname(s$truth), unname(truth_from_adj(s$adj_xx))))   # no hidden paths without violations
stopifnot(all(lengths(s$violations) == 0))

h <- sim_tiers(n = 500, d_z = 40, d_x = 10, p_hide_z = 0.5, seed = 7)
stopifnot(length(h$violations$hidden_z) == 20, h$params$d_z == 20)

m <- sim_tiers(n = 500, d_z = 20, d_x = 10, p_x_as_z = 0.2, seed = 7)
stopifnot(length(m$violations$x_as_z) == 2, m$params$d_x == 8, m$params$d_z == 22)

zx <- sim_tiers(n = 500, d_z = 20, d_x = 10, p_z_as_x = 0.1, seed = 7)
stopifnot(length(zx$violations$z_as_x) == 2, zx$params$d_x == 12)

f <- sim_tiers(n = 500, d_z = 20, d_x = 10, p_feedback = 0.2, seed = 7)
fb <- f$violations$feedback_z; stopifnot(length(fb) == 4)
for (z in fb) {
  xpar <- rownames(f$A_full)[f$A_full[, z] == 1 & grepl("^X", rownames(f$A_full))]
  stopifnot(length(xpar) >= 1, match(z, f$order) > match(xpar[1], f$order))
}
stopifnot(identical(rownames(f$truth), paste0("x", seq_len(f$params$d_x))))
cat("sim_tiers tests OK\n")
