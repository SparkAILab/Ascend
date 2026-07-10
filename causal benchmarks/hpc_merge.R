#!/usr/bin/env Rscript
# =============================================================================
# hpc_merge.R — Merge all task results into a single data.frame
#
# Usage: Rscript hpc_merge.R       (run from the directory holding results/)
#
# Output:
#   ascend_benchmark_merged.rds
#   ascend_benchmark_merged.csv
#
# Prints three summaries useful for the paper:
#   1. Mean F1 / precision / recall / coverage by method x n  (status==ok)
#   2. Runtime (median / mean / max seconds) by method x n     (status==ok)
#   3. Completion & failure mode (ok / timeout / error / skipped) by method x n
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

res_dir <- "results"
if (!dir.exists(res_dir)) stop("No results/ directory found")

job_dirs <- list.dirs(res_dir, recursive=FALSE, full.names=TRUE)
job_dirs <- job_dirs[grepl("job_\\d+$", job_dirs)]
cat(sprintf("Found %d job (combo) directories\n", length(job_dirs)))

all_rows <- list(); n_files <- 0L; n_failed <- 0L

for (jd in job_dirs) {
  param_file <- file.path(jd, "params.rds")
  params <- if (file.exists(param_file)) readRDS(param_file) else NULL
  
  rds_files <- list.files(jd, pattern="^n\\d+_rep\\d+\\.rds$", full.names=TRUE)
  
  for (f in rds_files) {
    df <- tryCatch(readRDS(f), error=function(e) {
      message(sprintf("  Cannot read %s: %s", f, e$message)); NULL
    })
    if (is.null(df) || nrow(df)==0) { n_failed <- n_failed+1; next }
    
    if (!is.null(params)) {
      if (!"d_x" %in% colnames(df)) df$d_x <- params$d_x
      if (!"d_z" %in% colnames(df)) df$d_z <- params$d_z
      if (!"r2"  %in% colnames(df)) df$r2  <- params$r2
      if (!"sp"  %in% colnames(df)) df$sp  <- params$sp
    }
    if (!"status"      %in% colnames(df)) df$status      <- "ok"
    if (!"elapsed_sec" %in% colnames(df)) df$elapsed_sec <- NA_real_
    if (!"n_over_dz"   %in% colnames(df)) df$n_over_dz   <- df$n / df$d_z
    
    all_rows[[length(all_rows)+1]] <- df
    n_files <- n_files+1
  }
}

cat(sprintf("Loaded %d replicate files (%d unreadable)\n", n_files, n_failed))
if (length(all_rows)==0) stop("No data found")

merged <- bind_rows(all_rows)

num_cols <- c("precision","recall","f1","accuracy","coverage",
              "tp","fp","fn","tn","unres_tp","unres_tn",
              "elapsed_sec","n_over_dz",
              "n","d_x","d_z","r2","sp","rep","job","task_id")
for (col in intersect(num_cols, colnames(merged)))
  merged[[col]] <- suppressWarnings(as.numeric(merged[[col]]))

# ── Overview ──────────────────────────────────────────────────────────────────
cat(sprintf("\nTotal rows     : %d\n", nrow(merged)))
cat(sprintf("Unique combos  : %d\n", length(unique(merged$job))))
cat(sprintf("Methods        : %s\n", paste(sort(unique(merged$method)), collapse=", ")))
cat(sprintf("n values       : %s\n", paste(sort(unique(merged$n)), collapse=", ")))
status_tab <- table(merged$status)
cat("Status counts  : ",
    paste(sprintf("%s=%d", names(status_tab), as.integer(status_tab)), collapse="  "),
    "\n", sep="")

# ── 1. Accuracy by method x n (status==ok only) ──────────────────────────────
cat("\n=== Mean F1 / precision / recall / coverage by method x n (status==ok) ===\n")
acc_df <- merged %>%
  filter(status=="ok") %>%
  group_by(method, n) %>%
  summarise(n_reps   = n(),
            f1_mean  = round(mean(f1,        na.rm=TRUE), 3),
            pr_mean  = round(mean(precision, na.rm=TRUE), 3),
            re_mean  = round(mean(recall,    na.rm=TRUE), 3),
            cov_mean = round(mean(coverage,  na.rm=TRUE), 3),
            .groups="drop") %>%
  arrange(method, n)
print(acc_df, n=Inf)

# ── 2. Runtime by method x n (status==ok only) ───────────────────────────────
cat("\n=== Runtime seconds by method x n (status==ok) ===\n")
time_df <- merged %>%
  filter(status=="ok", !is.na(elapsed_sec)) %>%
  group_by(method, n) %>%
  summarise(n_reps     = n(),
            t_median_s = round(median(elapsed_sec), 2),
            t_mean_s   = round(mean(elapsed_sec),   2),
            t_max_s    = round(max(elapsed_sec),    2),
            .groups="drop") %>%
  arrange(method, n)
print(time_df, n=Inf)

# ── 3. Completion & failure mode by method x n ───────────────────────────────
# Fraction of attempts that finished within the per-method budget — the core
# scalability statistic (e.g. "PC completes 0% of runs at n>=65536, d_x=80").
cat("\n=== Completion / failure mode by method x n (fraction of attempts) ===\n")
comp_df <- merged %>%
  group_by(method, n) %>%
  summarise(attempts   = n(),
            frac_ok      = round(mean(status=="ok"),              3),
            frac_timeout = round(mean(status=="timeout"),         3),
            frac_error   = round(mean(status=="error"),           3),
            frac_skipped = round(mean(status=="skipped_timeout"), 3),
            .groups="drop") %>%
  arrange(method, n)
print(comp_df, n=Inf)

# ── Save ──────────────────────────────────────────────────────────────────────
out_file <- "ascend_benchmark_merged.rds"
saveRDS(merged, out_file)
cat(sprintf("\nSaved: %s\n", out_file))
write.csv(merged, sub("\\.rds$", ".csv", out_file), row.names=FALSE)
cat(sprintf("Saved: %s\n", sub("\\.rds$", ".csv", out_file)))