## ======================================================================
## benchmarks/grn_multi/run_task.R
## Runs ONE R method on ONE dataset and writes its raw output to
## WORK/raw/<method>/<dataset_id>.* (see R/common.R). Everything is kept:
## scores, graphs, run time and the method's own counters, so metrics
## can be recomputed later without re-running any method.
##
## Reviewer comments: 1, 2, 6 (run time and test counts per run)
## How to run: normally called by run_chunk.sh under a wall-clock limit;
##   by hand: Rscript benchmarks/grn_multi/run_task.R <dataset_id> <method>
## Full write-up: benchmarks/grn_multi/README.md
## ======================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript run_task.R <dataset_id> <method>")
id <- args[1]; method <- args[2]

# single-threaded everywhere, so run times are comparable across methods
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  setwd(d)
})
source("benchmarks/grn_multi/config.R")
source(file.path(GRNM_DIR, "R", "common.R"))
source(file.path(GRNM_DIR, "R", "methods.R"))
if (!method %in% names(R_METHODS)) stop("unknown R method: ", method)

dir.create(file.path(RAW_DIR, method), recursive = TRUE, showWarnings = FALSE)
D <- read_dataset(dataset_dir(id))
s <- SETTINGS[[if (method %in% c("clr", "mrnet", "aracne")) "minet" else method]]
set.seed(D$meta$seed)

t0  <- proc.time()
out <- tryCatch(R_METHODS[[method]](D, s), error = function(e) e)
tt  <- proc.time() - t0

meta <- list(dataset = id, method = method, lang = "R",
             runtime_s = unname(tt[["elapsed"]]),
             cpu_s = unname(tt[["user.self"]] + tt[["sys.self"]]),
             settings = s, host = Sys.info()[["nodename"]],
             finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             r_version = R.version.string)
if (inherits(out, "error")) {
  meta$status <- "error"; meta$error <- conditionMessage(out)
} else {
  meta$status <- "ok"
  for (w in c("score", "anc", "adj")) if (!is.null(out[[w]])) {
    M <- as.matrix(out[[w]]); storage.mode(M) <- "double"
    stopifnot(identical(rownames(M), D$xlabs), identical(colnames(M), D$xlabs))
    write_matrix(M, raw_path(method, id, paste0(w, ".csv.gz")))
  }
  meta$outputs <- intersect(c("score", "anc", "adj"), names(out))
  if (length(out$stats)) meta$stats <- out$stats
}
# meta.json is written last: its presence marks the run as finished
write_json(meta, raw_path(method, id, "meta.json"), auto_unbox = TRUE, pretty = TRUE,
           digits = NA, null = "null", na = "null")
cat(sprintf("%s %s %s %.1fs%s\n", id, method, meta$status, meta$runtime_s,
            if (!is.null(meta$error)) paste0(" : ", meta$error) else ""))
