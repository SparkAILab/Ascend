#!/usr/bin/env Rscript
# =============================================================================
# hpc_test.R — Smoke test: job 1, n=512 only, 2 replicates
# Verifies all 5 methods run, produce output, and save correctly.
# Usage: Rscript hpc_test.R
# =============================================================================
cat("=== SMOKE TEST: verifying hpc_run.R on job 1 ===\n\n")

# Read hpc_run.R and patch constants for a quick test run
lines <- readLines("hpc_run.R")

# Override: 2 reps, single n value, job index 1
lines <- sub("^N_REP\\s*<-.*$",  "N_REP <- 2L",          lines)
lines <- sub("^N_VEC\\s*<-.*$",  "N_VEC <- c(512L)",     lines)
lines <- sub("^args\\s*<-.*$",   "args <- c('1')",        lines)

tmp <- tempfile(fileext=".R")
writeLines(lines, tmp)
source(tmp, local=FALSE)
unlink(tmp)