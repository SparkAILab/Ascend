## ======================================================================
## benchmarks/grn_multi/01_make_datasets.R
## Simulates every two-tier dataset of the benchmark (arms lin, nonlin,
## large in config.R) and writes each one, with its ground truth, to
## WORK/data/<dataset_id>/. The SERGIO arm is written by py/make_sergio.py.
##
## Reviewer comments: 1, 2, 9 (the data behind every comparison), minor 8
## How to run: Rscript benchmarks/grn_multi/01_make_datasets.R
##   Skips datasets that already exist, so it is safe to re-run.
##   ~1 min for the full design; GRNM_TEST=1 writes the small test design.
## Output per dataset (all methods read exactly these files):
##   data.csv.gz       n x (z1..z_dz, x1..x_dx)
##   truth_direct.csv  X x X, [a, b] = 1 iff a -> b in the true graph
##   truth_anc.csv     X x X, [a, b] = 1 iff a is an ancestor of b
##   meta.json         arm, cell settings, replicate, seed, labels
## Full write-up: benchmarks/grn_multi/README.md
## ======================================================================

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  setwd(d)
})
source("benchmarks/grn_multi/config.R")
source(file.path(GRNM_DIR, "R", "common.R"))
source(file.path(ROOT, "R", "ascend.R"))          # sim_dat()
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x) trimws(formatC(x, format = "fg", digits = 3))

n_new <- 0L
arm_names <- names(ARMS)
for (ai in seq_along(arm_names)) {
  arm <- arm_names[ai]; A <- ARMS[[arm]]
  if (A$generator != "sim_dat") next
  for (ci in seq_len(nrow(A$grid))) {
    cell <- as.list(A$grid[ci, , drop = FALSE])
    for (rep in seq_len(A$reps)) {
      id <- sprintf("%s_n%d_sp%s_r2%s_rep%03d", arm, cell$n, fmt(cell$sp), fmt(cell$r2), rep)
      dir <- dataset_dir(id)
      if (file.exists(file.path(dir, "meta.json"))) next
      # one seed per dataset, stable across machines and design edits
      seed <- 1e6 * ai + 1e3 * ci + rep
      p <- c(cell, A$fixed)
      sim <- sim_dat(n = p$n, d_z = p$d_z, d_x = p$d_x, r2 = p$r2, sp = p$sp,
                     lin_pr = p$lin_pr, p_cross = p$p_cross, x_effect = p$x_effect,
                     z_scale = TRUE, seed = seed)
      xl <- paste0("x", seq_len(p$d_x)); zl <- paste0("z", seq_len(p$d_z))
      A_dir <- t(sim$adj_xx); A_dir[is.na(A_dir)] <- 0; diag(A_dir) <- 0
      dimnames(A_dir) <- list(xl, xl)
      A_anc <- transitive_closure(A_dir); dimnames(A_anc) <- list(xl, xl)
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      con <- gzfile(file.path(dir, "data.csv.gz"), "w")
      write.csv(signif(as.matrix(sim$dat[, c(zl, xl)]), 7), con, row.names = FALSE)
      close(con)
      write_matrix(A_dir, file.path(dir, "truth_direct.csv"))
      write_matrix(A_anc, file.path(dir, "truth_anc.csv"))
      meta <- list(id = id, arm = arm, rep = rep, seed = seed,
                   cell = cell, params = p, n = p$n, d_z = p$d_z, d_x = p$d_x,
                   zlabs = zl, xlabs = xl,
                   n_direct = sum(A_dir), n_anc = sum(A_anc),
                   generator = "sim_dat (R/ascend.R), z_scale = TRUE")
      write_json(meta, file.path(dir, "meta.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
      n_new <- n_new + 1L
    }
  }
}
cat(sprintf("%d new datasets written to %s\n", n_new, DATA_DIR))
