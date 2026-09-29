## ======================================================================
## benchmarks/grn_multi/05_test_report.R
## After the test run: did every method run on every arm, how long did it
## take, and what will the full run cost?
##
## How to run: GRNM_TEST=1 Rscript benchmarks/grn_multi/05_test_report.R
##   (slurm/test_pipeline.sbatch runs it). Writes WORK/results/TEST_REPORT.txt
## ======================================================================

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  setwd(d)
})
source("benchmarks/grn_multi/config.R")
suppressPackageStartupMessages(library(data.table))
runs <- fread(file.path(RES_DIR, "runs.csv"))
test_res <- RES_DIR

# the FULL design, for the projection
local({
  Sys.setenv(GRNM_TEST = "0"); on.exit(Sys.setenv(GRNM_TEST = "1"))
  source("benchmarks/grn_multi/config.R", local = TRUE)
  full <<- rbindlist(lapply(names(ARMS), function(a) {
    A <- ARMS[[a]]
    if (a == "sergio") data.table(arm_key = paste0("sergio_", A$grid$ds), k = A$reps)[, .(n_datasets = sum(k)), by = arm_key]
    else data.table(arm_key = a, n_datasets = nrow(A$grid) * A$reps)
  }))
  full_methods <<- METHODS
})

# runs a method is not meant to do (skip_arms in config.R) are not failures
skips <- setNames(strsplit(as.character(full_methods$skip_arms), ";"), full_methods$method)
runs <- runs[!mapply(function(m, k) k %in% trimws(skips[[m]]) || sub("_.*", "", k) %in% trimws(skips[[m]]),
                     method, arm_key)]

lines <- c(sprintf("TEST REPORT  %s", format(Sys.time(), "%Y-%m-%d %H:%M")), "")
st <- dcast(runs, method ~ arm_key, value.var = "status", fun.aggregate = function(s) paste(sort(unique(s)), collapse = "/"))
lines <- c(lines, "1. Status per method and arm (every cell should read 'ok'):", "",
           capture.output(print(st, row.names = FALSE)), "")
bad <- runs[status != "ok"]
if (nrow(bad)) {
  lines <- c(lines, "   Not ok:", capture.output(print(bad[, .(dataset, method, status, error = substr(error, 1, 90))], row.names = FALSE)),
             "   timeout = needs longer than the test cap (2 h); it gets its full budget in the real run.",
             "   error / crashed = look at WORK/logs/<method>/<dataset>.log before the full run.", "")
} else lines <- c(lines, "   All runs ok.", "")

rt <- runs[status %in% c("ok", "timeout"), .(sec = max(runtime_s)), by = .(method, arm_key)]
rt <- merge(rt, full, by = "arm_key")
rt <- merge(rt, as.data.table(full_methods)[, .(method, skip_arms, budget_h)], by = "method")
rt <- rt[!mapply(function(s, k) k %in% strsplit(s, ";")[[1]] || sub("_.*", "", k) %in% strsplit(s, ";")[[1]],
                 skip_arms, arm_key)]
# a run that hit the test cap may need its whole budget in the real run
tmo <- runs[status == "timeout", .(timed_out = TRUE), by = .(method, arm_key)]
rt <- merge(rt, tmo, by = c("method", "arm_key"), all.x = TRUE)
rt[, full_budget_s := budget_h * 3600 * ifelse(arm_key %in% names(BUDGET_MULT), BUDGET_MULT[arm_key], 1)]
rt[!is.na(timed_out), sec := full_budget_s]
rt[, cpu_h := sec * n_datasets / 3600]
lines <- c(lines, "2. Longest test run (s) per method and arm, and projected CPU hours for the full design:", "",
           capture.output(print(dcast(rt, method ~ arm_key, value.var = "sec"), row.names = FALSE)), "",
           "   (A method that timed out in the test is counted at its full budget.)",
           sprintf("   Projected total: about %.0f CPU hours (upper estimate: the slowest test run times the number of datasets).",
                   sum(rt$cpu_h)),
           sprintf("   With %s workers per job that is about %.0f job-hours.", Sys.getenv("GRNM_CPUS", "8"),
                   sum(rt$cpu_h) / as.numeric(Sys.getenv("GRNM_CPUS", "8"))), "",
           "3. Next: bash benchmarks/grn_multi/slurm/submit_all.sh", "")
writeLines(lines, file.path(test_res, "TEST_REPORT.txt"))
cat(lines, sep = "\n")
