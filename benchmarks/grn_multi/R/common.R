## ======================================================================
## benchmarks/grn_multi/R/common.R
## Shared helpers: dataset I/O, graph utilities (closure, reduction,
## CPDAG -> ancestral), and the file layout of one method run.
##
## Conventions used everywhere in this benchmark
##   Every graph or score matrix is indexed [from, to]:
##     adj[a, b] = 1           edge a -> b
##     adj[a, b] = adj[b, a]   = 1  undirected edge a - b
##     anc[a, b] = 1           a is an ancestor of b   (eval_metrics.R)
##     score[a, b]             confidence that a regulates b
##   pcalg's as(fit, "amat") uses the TRANSPOSED convention and CBL's
##   output is [descendant, ancestor]; both are converted on the way in.
## ======================================================================

suppressPackageStartupMessages(library(jsonlite))

## --- dataset I/O ---------------------------------------------------------
dataset_dir <- function(id) file.path(DATA_DIR, id)

read_dataset <- function(dir) {
  meta <- fromJSON(file.path(dir, "meta.json"))
  dat  <- read.csv(file.path(dir, "data.csv.gz"), check.names = FALSE)
  list(meta = meta, dat = dat,
       zlabs = meta$zlabs, xlabs = meta$xlabs,
       truth_direct = read_matrix(file.path(dir, "truth_direct.csv")),
       truth_anc    = read_matrix(file.path(dir, "truth_anc.csv")))
}

write_matrix <- function(M, path) {
  df <- data.frame(row = rownames(M), as.data.frame(M, check.names = FALSE),
                   check.names = FALSE)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "w") else file(path, "w")
  on.exit(close(con))
  write.csv(df, con, row.names = FALSE, na = "NA")
}

read_matrix <- function(path) {
  df <- read.csv(path, check.names = FALSE, na.strings = "NA")
  M  <- as.matrix(df[, -1, drop = FALSE]); rownames(M) <- df[[1]]
  storage.mode(M) <- "double"; M
}

## --- graph utilities -----------------------------------------------------
transitive_closure <- function(A) {
  A <- as.matrix(A); A[is.na(A)] <- 0; diag(A) <- 0
  R <- (A > 0) * 1
  repeat {
    R2 <- ((R + (R %*% R)) > 0) * 1
    if (identical(R2, R)) break
    R <- R2
  }
  diag(R) <- 0; R
}

# Minimal edge set with the same reachability (for acyclic A)
transitive_reduction <- function(A) {
  A <- (as.matrix(A) > 0) * 1; A[is.na(A)] <- 0; diag(A) <- 0
  R <- transitive_closure(A)
  redundant <- ((A %*% R) > 0) * 1        # a -> c -> ... -> b exists
  out <- A * (1 - redundant); dimnames(out) <- dimnames(A); out
}

# Any mixed graph over ALL variables (adj[a,b]=1 a->b, both = undirected)
# -> ancestral matrix over `keep` in the eval_metrics convention.
#   Directed paths are closed over the full graph, so a path through a
#   background variable counts. An undirected edge between two foreground
#   variables that no directed path already orders is "related, direction
#   not committed": 1 in both directions.
graph_to_ancestral <- function(adj, keep) {
  adj <- as.matrix(adj); adj[is.na(adj)] <- 0; diag(adj) <- 0
  dir <- (adj > 0 & t(adj) == 0) * 1
  und <- (adj > 0 & t(adj) > 0)
  R <- transitive_closure(dir)
  dimnames(R) <- dimnames(adj)
  R <- R[keep, keep, drop = FALSE]
  U <- und[keep, keep, drop = FALSE] & R == 0 & t(R) == 0
  R[U] <- 1
  diag(R) <- NA; R
}

## --- one method run on disk ---------------------------------------------
## RAW_DIR/<method>/<dataset_id>.meta.json   status, runtime, settings, stats
## RAW_DIR/<method>/<dataset_id>.score.csv.gz  score matrix over X
## RAW_DIR/<method>/<dataset_id>.anc.csv.gz    ancestral matrix over X (causal)
## RAW_DIR/<method>/<dataset_id>.adj.csv.gz    native graph over X (causal)
raw_path <- function(method, id, what) file.path(RAW_DIR, method, paste0(id, ".", what))

run_done <- function(method, id) file.exists(raw_path(method, id, "meta.json"))
