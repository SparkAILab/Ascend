#!/usr/bin/env Rscript
## ======================================================================
## benchmarks/causal_grid/hpc_test.R
## Smoke test for hpc_run.R (task 1, 2 replicates).
##
## Reviewer comments: none (test helper)
## Revision notes: Unchanged.
## How to run: cd benchmarks/causal_grid && Rscript hpc_test.R
## Full write-up: docs/REVISION_REPORT.md
## ======================================================================

# =============================================================================
# hpc_test.R — Smoke test: task 1 (combo 1, smallest n), 2 replicates.
# Verifies all 5 methods run, produce output, time correctly, and save.
# Usage: Rscript hpc_test.R
# =============================================================================
cat("=== SMOKE TEST: verifying hpc_run.R on task 1 (combo 1, n=N_VEC[1]) ===\n\n")

lines <- readLines("hpc_run.R")

# Override: 2 reps, and force task id 1 (decodes to combo 1 at the smallest n).
lines <- sub("^N_REP\\s*<-.*$", "N_REP <- 2L",   lines)
lines <- sub("^args\\s*<-.*$",  "args <- c('1')", lines)

tmp <- tempfile(fileext=".R")
writeLines(lines, tmp)
source(tmp, local=FALSE)
unlink(tmp)

cat("\n=== Smoke test complete: check results/job_001/ for n*_rep*.rds (5 rows each) ===\n")