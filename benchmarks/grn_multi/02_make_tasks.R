## ======================================================================
## benchmarks/grn_multi/02_make_tasks.R
## Turns (datasets x methods) into a task list and packs it into chunks,
## one chunk per SLURM array job, so that no job comes near the 48 h limit.
##
## How the packing works
##   Each task (one method on one dataset) has
##     budget_s  hard wall-clock limit from config.R (timeout -> "timeout")
##     est_s     expected run time: from the test run's timings if present
##               (GRNM_PILOT, default work_test/results/runs.csv), otherwise
##               from the rough defaults below. Always capped at budget_s.
##   Tasks are assigned longest-first to the least-loaded chunk. A chunk runs
##   its tasks on GRNM_CPUS parallel workers; its expected length
##   (sum est_s / workers, plus its longest task) is kept below
##   GRNM_CHUNK_HOURS (default 20 h), far inside the 48 h limit.
##   Only tasks without a finished run are listed, so re-running this after
##   a partial run packs just what is left.
##
## How to run: Rscript benchmarks/grn_multi/02_make_tasks.R
##   env: GRNM_CPUS (8), GRNM_CHUNK_HOURS (20), GRNM_RETRY_ERRORS=1 to
##   delete failed/timed-out runs first so they are attempted again.
## Output: WORK/tasks.csv (task, chunk, dataset, method, lang, budget_s, est_s)
##         WORK/chunks.csv (chunk, n_tasks, expected hours)
## ======================================================================

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  setwd(d)
})
source("benchmarks/grn_multi/config.R")
source(file.path(GRNM_DIR, "R", "common.R"))
suppressPackageStartupMessages(library(data.table))

CPUS  <- as.integer(Sys.getenv("GRNM_CPUS", "8"))
CHUNK_H <- as.numeric(Sys.getenv("GRNM_CHUNK_HOURS", "20"))
RETRY <- Sys.getenv("GRNM_RETRY_ERRORS", "0") == "1"

ids <- sort(list.dirs(DATA_DIR, full.names = FALSE, recursive = FALSE))
ids <- ids[file.exists(file.path(DATA_DIR, ids, "meta.json"))]
if (!length(ids)) stop("No datasets in ", DATA_DIR, ": run 01_make_datasets.R and py/make_sergio.py first")
dsm <- rbindlist(lapply(ids, function(id) {
  m <- fromJSON(file.path(DATA_DIR, id, "meta.json"))
  data.table(dataset = id, arm = m$arm, arm_key = if (!is.null(m$arm_key)) m$arm_key else m$arm,
             n = m$n, d_x = m$d_x, d_z = m$d_z)
}))

tasks <- CJ(dataset = dsm$dataset, method = METHODS$method)
tasks <- merge(tasks, dsm, by = "dataset")
tasks <- merge(tasks, as.data.table(METHODS)[, .(method, lang, class, budget_h, skip_arms)], by = "method")
skip <- mapply(function(sa, a, k) { s <- strsplit(sa, ";")[[1]]; a %in% s || k %in% s },
               tasks$skip_arms, tasks$arm, tasks$arm_key)
tasks <- tasks[!skip]
mult <- ifelse(tasks$arm_key %in% names(BUDGET_MULT), BUDGET_MULT[tasks$arm_key], 1)
tasks[, budget_s := as.integer(pmin(budget_h * mult, MAX_BUDGET_H) * 3600)]

## finished runs are not listed again
meta_file <- file.path(RAW_DIR, tasks$method, paste0(tasks$dataset, ".meta.json"))
if (RETRY) {
  st <- vapply(meta_file, function(f) if (file.exists(f)) fromJSON(f)$status else "none", "")
  bad <- st %in% c("error", "timeout")
  if (any(bad)) { file.remove(meta_file[bad]); cat(sum(bad), "failed runs cleared for retry\n") }
}
done <- file.exists(meta_file)
cat(sprintf("%d tasks in the design, %d already finished, %d to run\n",
            nrow(tasks), sum(done), sum(!done)))
tasks <- tasks[!done]

## expected run times
# rough defaults (seconds) as a function of the number of foreground genes;
# replaced by the test run's measurements whenever those exist
default_est <- function(method, d_x, n) {
  base <- switch(method,
                 ascend = 2, pc_tiered = 5, hc_tiered = 3, ges = 3, cbl = 300,
                 genie3 = 120, grnboost2 = 10, regdiffusion = 30, pidc = 10,
                 clr = 1, mrnet = 1, aracne = 1, wgcna = 10, ppcor = 1, pearson = 1, 60)
  expo <- switch(method, cbl = 2.5, ascend = 2, pc_tiered = 2, ges = 2, hc_tiered = 2,
                 genie3 = 1.2, pidc = 3, 1.5)
  base * (d_x / 15)^expo * (n / 2000)
}
tasks[, est_s := mapply(default_est, method, d_x, n)]
pilot <- Sys.getenv("GRNM_PILOT", file.path(GRNM_DIR, "work_test", "results", "runs.csv"))
if (file.exists(pilot) && !GRNM_TEST) {
  pr <- fread(pilot)[status == "ok"]
  # scale each pilot time to the task's size (d_x, n) and keep the largest
  pr <- pr[, .(pil_s = max(runtime_s), pil_dx = d_x[which.max(runtime_s)], pil_n = n[which.max(runtime_s)]),
           by = .(method, arm_key)]
  tasks <- merge(tasks, pr, by = c("method", "arm_key"), all.x = TRUE)
  tasks[!is.na(pil_s), est_s := 1.5 * pil_s * pmax(1, (d_x / pil_dx))^2 * pmax(1, n / pil_n)]
  # timed-out pilot runs: assume the full budget
  pto <- fread(pilot)[status == "timeout", unique(paste(method, arm_key))]
  tasks[paste(method, arm_key) %in% pto, est_s := budget_s]
  tasks[, c("pil_s", "pil_dx", "pil_n") := NULL]
  cat("Run-time estimates from the test run:", pilot, "\n")
} else cat("No test-run timings found: using rough default run-time estimates\n")
tasks[, est_s := pmin(est_s, budget_s)]

## pack: longest first into the chunk with the smallest expected length
cap_s <- CHUNK_H * 3600
setorder(tasks, -est_s)
load <- numeric(0); longest <- numeric(0); chunk <- integer(nrow(tasks))
for (i in seq_len(nrow(tasks))) {
  e <- tasks$est_s[i]
  len_if <- (load + e) / CPUS + pmax(longest, e)
  j <- if (length(load)) which.min(len_if) else integer(0)
  if (!length(j) || len_if[j] > cap_s) { load <- c(load, 0); longest <- c(longest, 0); j <- length(load) }
  load[j] <- load[j] + e; longest[j] <- max(longest[j], e); chunk[i] <- j
}
tasks[, chunk := chunk]
setorder(tasks, chunk, -est_s)
tasks[, task := seq_len(.N)]
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
fwrite(tasks[, .(task, chunk, dataset, method, lang, budget_s, est_s = round(est_s))],
       file.path(WORK, "tasks.csv"))
ch <- tasks[, .(n_tasks = .N, cpu_hours = round(sum(est_s) / 3600, 1),
                expected_hours = round(sum(est_s) / CPUS / 3600 + max(est_s) / 3600, 1)), by = chunk]
fwrite(ch, file.path(WORK, "chunks.csv"))
cat(sprintf("%d chunks for %d workers each; expected chunk length %.1f-%.1f h; total %.0f CPU hours\n",
            nrow(ch), CPUS, min(ch$expected_hours), max(ch$expected_hours), sum(ch$cpu_hours)))
cat("Submit with: sbatch --array=1-", nrow(ch), " benchmarks/grn_multi/slurm/methods.sbatch\n", sep = "")
