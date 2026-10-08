## ======================================================================
## realdata/drosophila/dgrp_controls.R
## One ASCEND run of the DGRP statistical controls. Run with no argument it
## prints the list of runs; with a run number it does that run and saves a
## small summary (edge counts, degrees, edge list, CI-test counts).
##
## Reviewer comments: 7 (negative controls with permuted genotypes and with
##   randomised tier labels, stability under subsampling, selection bias of
##   the top-250 genes, number of tests), minor 7 (filtering counts)
## The runs (seeds are the run's replicate number):
##   observed      1 run    the paper's analysis, unchanged
##   perm_geno   100 runs   genotypes permuted across lines (rows of Z
##                          shuffled together): breaks every SNP-gene link but
##                          keeps linkage among SNPs and correlation among genes
##   swap_tiers   20 runs   250 of the 750 variables drawn at random as the
##                          foreground, the rest as the background
##   subsample   100 runs   80% of lines, drawn without replacement
##   random_genes 20 runs   250 genes drawn at random among all genes passing
##                          the expression filter (instead of the top 250 by
##                          variance)
## How to run:
##   one run:  Rscript realdata/drosophila/dgrp_controls.R 1
##   all runs on CREATE: bash realdata/drosophila/slurm/submit_dgrp.sh
##   then:     Rscript realdata/drosophila/dgrp_controls_summarise.R
##   DGRP_DATA = folder with the two data files (default realdata/drosophila)
##   DGRP_OUT  = output folder (default realdata/drosophila/results)
## Full write-up: docs/REVISION_REPORT.md
## ======================================================================

ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
source(file.path(ROOT, "R", "ascend.R"))
source(file.path(ROOT, "realdata", "drosophila", "dgrp_data.R"))
DATA <- Sys.getenv("DGRP_DATA", file.path(ROOT, "realdata", "drosophila"))
OUT  <- Sys.getenv("DGRP_OUT",  file.path(ROOT, "realdata", "drosophila", "results"))
dir.create(file.path(OUT, "runs"), recursive = TRUE, showWarnings = FALSE)

N_REPS <- c(observed = 1, perm_geno = 100, swap_tiers = 20, subsample = 100, random_genes = 20)
runs <- data.frame(type = rep(names(N_REPS), N_REPS),
                   rep  = unlist(lapply(N_REPS, seq_len)), stringsAsFactors = FALSE)
runs$id <- seq_len(nrow(runs))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat(nrow(runs), "runs:\n"); print(table(factor(runs$type, names(N_REPS))))
  quit(save = "no")
}
r <- runs[as.integer(args[1]), ]
out_file <- file.path(OUT, "runs", sprintf("%03d_%s_%03d.rds", r$id, r$type, r$rep))
if (file.exists(out_file)) { cat("already done:", out_file, "\n"); quit(save = "no") }

inp <- dgrp_input(DATA)
z <- inp$z; x <- inp$x
set.seed(1000 * match(r$type, names(N_REPS)) + r$rep)

if (r$type == "perm_geno") {
  z <- z[sample(nrow(z)), , drop = FALSE]
} else if (r$type == "subsample") {
  keep <- sort(sample(nrow(z), round(0.8 * nrow(z))))
  z <- z[keep, , drop = FALSE]; x <- x[keep, , drop = FALSE]
} else if (r$type == "random_genes") {
  x <- inp$x_pool[, sort(sample(ncol(inp$x_pool), ncol(inp$x))), drop = FALSE]
} else if (r$type == "swap_tiers") {
  all <- cbind(z, x)
  fg <- sort(sample(ncol(all), ncol(x)))
  z <- all[, -fg, drop = FALSE]; x <- all[, fg, drop = FALSE]
}

obj <- dgrp_ascend_obj(z, x)
t0 <- Sys.time()
M <- ascend(obj, alpha = 0.05, alpha_mb = 0.05, fdr = TRUE, min_votes = 1, verbose = FALSE)
secs <- as.numeric(Sys.time() - t0, units = "secs")

res <- c(list(type = r$type, rep = r$rep, n_lines = nrow(z), seconds = secs,
              genes = colnames(x)),
         dgrp_summary(M, colnames(x)))
if (r$type == "observed") { res$M <- M; res$flow <- inp$flow }
tmp <- paste0(out_file, ".tmp"); saveRDS(res, tmp); invisible(file.rename(tmp, out_file))
cat(sprintf("%s %d: %d directed, %d undirected, max out-degree %d, %.0f s\n",
            r$type, r$rep, res$n_directed, res$n_undirected, max(res$out_degree), secs))
