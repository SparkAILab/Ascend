#!/usr/bin/env Rscript
# =============================================================================
# hpc_merge.R — Merge all job results into a single data.frame
#
# Usage: Rscript hpc_merge.R
#
# Run from the same directory as results/
# Output: ascend_benchmark_merged.rds
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

res_dir <- "results"
if (!dir.exists(res_dir)) stop("No results/ directory found")

job_dirs <- list.dirs(res_dir, recursive=FALSE, full.names=TRUE)
job_dirs <- job_dirs[grepl("job_\\d+$", job_dirs)]
cat(sprintf("Found %d job directories\n", length(job_dirs)))

all_rows <- list()
n_files  <- 0L
n_failed <- 0L

for (jd in job_dirs) {
  # Load params for this job (fallback: extract from row data)
  param_file <- file.path(jd, "params.rds")
  params <- if (file.exists(param_file)) readRDS(param_file) else NULL
  
  rds_files <- list.files(jd, pattern="^n\\d+_rep\\d+\\.rds$", full.names=TRUE)
  
  for (f in rds_files) {
    df <- tryCatch(readRDS(f), error=function(e) {
      message(sprintf("  Cannot read %s: %s", f, e$message)); NULL
    })
    if (is.null(df) || nrow(df)==0) { n_failed <- n_failed+1; next }
    
    # Back-fill params if columns missing (older format compatibility)
    if (!is.null(params)) {
      if (!"d_x" %in% colnames(df)) df$d_x <- params$d_x
      if (!"d_z" %in% colnames(df)) df$d_z <- params$d_z
      if (!"r2"  %in% colnames(df)) df$r2  <- params$r2
      if (!"sp"  %in% colnames(df)) df$sp  <- params$sp
    }
    if (!"status" %in% colnames(df)) df$status <- "ok"
    
    all_rows[[length(all_rows)+1]] <- df
    n_files <- n_files+1
  }
}

cat(sprintf("Loaded %d replicate files (%d unreadable)\n", n_files, n_failed))

if (length(all_rows)==0) stop("No data found")

merged <- bind_rows(all_rows)

# Ensure numeric columns are numeric
num_cols <- c("precision","recall","f1","accuracy","coverage",
              "tp","fp","fn","tn","unres_tp","unres_tn",
              "n","d_x","d_z","r2","sp","rep","job")
for (col in intersect(num_cols, colnames(merged)))
  merged[[col]] <- suppressWarnings(as.numeric(merged[[col]]))

# Report completeness
cat(sprintf("\nTotal rows     : %d\n", nrow(merged)))
cat(sprintf("Unique jobs    : %d\n", length(unique(merged$job))))
cat(sprintf("Methods        : %s\n", paste(sort(unique(merged$method)), collapse=", ")))
cat(sprintf("n values       : %s\n", paste(sort(unique(merged$n)), collapse=", ")))
cat(sprintf("Status ok/fail : %d / %d\n",
            sum(merged$status=="ok", na.rm=TRUE),
            sum(merged$status!="ok", na.rm=TRUE)))

# Quick summary: mean F1 by method and n
cat("\n=== Mean F1 by method × n (status=ok only) ===\n")
sum_df <- merged %>%
  filter(status=="ok") %>%
  group_by(method, n) %>%
  summarise(
    n_reps    = n(),
    f1_mean   = round(mean(f1,        na.rm=TRUE), 3),
    pr_mean   = round(mean(precision, na.rm=TRUE), 3),
    re_mean   = round(mean(recall,    na.rm=TRUE), 3),
    cov_mean  = round(mean(coverage,  na.rm=TRUE), 3),
    .groups="drop"
  ) %>%
  arrange(method, n)
print(sum_df, n=Inf)

out_file <- "ascend_benchmark_merged.rds"
saveRDS(merged, out_file)
cat(sprintf("\nSaved: %s\n", out_file))

# Also save CSV for convenience
write.csv(merged, sub("\\.rds$", ".csv", out_file), row.names=FALSE)
cat(sprintf("Saved: %s\n", sub("\\.rds$", ".csv", out_file)))