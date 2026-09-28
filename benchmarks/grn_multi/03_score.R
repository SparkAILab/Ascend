## ======================================================================
## benchmarks/grn_multi/03_score.R
## Scores every finished run against both ground truths and writes two
## tables. No method is re-run: this reads WORK/raw only, so it can be
## repeated (e.g. after adding a metric) at no cost.
##
## Reviewer comments: 2, 3, 9, 11, minors 6 and 9
## How to run: Rscript benchmarks/grn_multi/03_score.R
##   (uses SLURM_CPUS_PER_TASK or GRNM_CORES cores; a few minutes)
## Output (WORK/results/):
##   runs.csv             one row per dataset x method: status, run time,
##                        CPU time, error message, method counters
##   metrics_long.csv.gz  one row per dataset x method x target x metric
## Metric definitions: R/scoring.R. Full write-up: README.md
## ======================================================================

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  setwd(d)
})
source("benchmarks/grn_multi/config.R")
source(file.path(GRNM_DIR, "R", "common.R"))
source(file.path(GRNM_DIR, "R", "scoring.R"))
source(file.path(ROOT, "R", "eval_metrics.R"))
suppressPackageStartupMessages({ library(parallel); library(data.table) })
stopifnot(requireNamespace("PRROC", quietly = TRUE))
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("GRNM_CORES", "1")))
ids <- sort(list.dirs(DATA_DIR, full.names = FALSE, recursive = FALSE))
ids <- ids[file.exists(file.path(DATA_DIR, ids, "meta.json"))]
cat(sprintf("Scoring %d datasets x %d methods on %d cores\n", length(ids), nrow(METHODS), cores))

read_run <- function(method, id, class) {
  mp <- raw_path(method, id, "meta.json")
  if (!file.exists(mp)) return(list(meta = list(status = "missing")))
  meta <- fromJSON(mp, simplifyVector = TRUE)
  out <- list(meta = meta, class = class)
  if (identical(meta$status, "ok"))
    for (w in c("score", "anc", "adj")) {
      f <- raw_path(method, id, paste0(w, ".csv.gz"))
      if (file.exists(f)) out[[w]] <- read_matrix(f)
    }
  out
}

score_dataset <- function(id) {
  D <- read_dataset(dataset_dir(id))
  m <- D$meta; xl <- D$xlabs
  cell <- m$cell
  base <- data.table(dataset = id, arm = m$arm,
                     arm_key = if (!is.null(m$arm_key)) m$arm_key else m$arm,
                     n = m$n, sp = if (!is.null(cell$sp)) cell$sp else NA_real_,
                     r2 = if (!is.null(cell$r2)) cell$r2 else NA_real_,
                     ds = if (!is.null(cell$ds)) cell$ds else NA_character_,
                     noise = if (!is.null(cell$noise)) cell$noise else NA_character_,
                     rep = m$rep, d_z = m$d_z, d_x = m$d_x)
  runs <- list(); mets <- list()
  outs <- setNames(lapply(seq_len(nrow(METHODS)), function(i)
    read_run(METHODS$method[i], id, METHODS$class[i])), METHODS$method)
  truths <- list(direct = D$truth_direct[xl, xl], anc = D$truth_anc[xl, xl])
  K_asc <- sapply(names(truths), function(tg)
    if (identical(outs$ascend$meta$status, "ok")) ascend_K(outs$ascend, tg, xl) else NA_real_)
  for (i in seq_len(nrow(METHODS))) {
    me <- METHODS$method[i]; o <- outs[[me]]; mt <- o$meta
    st <- if (!is.null(mt$stats)) unlist(mt$stats) else NULL
    runs[[me]] <- cbind(base, data.table(
      method = me, label = METHODS$label[i], class = METHODS$class[i],
      status = mt$status,
      runtime_s = if (!is.null(mt$runtime_s)) mt$runtime_s else NA_real_,
      cpu_s = if (!is.null(mt$cpu_s)) mt$cpu_s else NA_real_,
      error = if (!is.null(mt$error)) substr(mt$error, 1, 300) else NA_character_,
      stats = if (length(st)) paste(names(st), signif(st, 6), sep = "=", collapse = ";") else NA_character_))
    if (!identical(mt$status, "ok")) next
    for (tg in names(truths)) {
      r <- tryCatch(score_one(o, truths[[tg]], tg, K_asc[[tg]]),
                    error = function(e) { message(id, " ", me, " ", tg, ": ", conditionMessage(e)); NULL })
      if (is.null(r)) next
      r <- unlist(r)
      mets[[paste(me, tg)]] <- data.table(dataset = id, method = me, target = tg,
                                          metric = names(r), value = as.numeric(r))
    }
  }
  list(runs = rbindlist(runs, fill = TRUE), mets = rbindlist(mets))
}

res <- mclapply(ids, function(id) tryCatch(score_dataset(id), error = function(e) {
  message("dataset ", id, " failed: ", conditionMessage(e)); NULL }), mc.cores = cores)
runs <- rbindlist(lapply(res, `[[`, "runs"), fill = TRUE)
mets <- rbindlist(lapply(res, `[[`, "mets"))
fwrite(runs, file.path(RES_DIR, "runs.csv"))
fwrite(mets, file.path(RES_DIR, "metrics_long.csv.gz"))

cat("\nRun status by method:\n")
print(dcast(runs, method ~ status, fun.aggregate = length, value.var = "dataset"))
cat(sprintf("\nWrote %s and metrics_long.csv.gz (%d rows)\n", file.path(RES_DIR, "runs.csv"), nrow(mets)))
