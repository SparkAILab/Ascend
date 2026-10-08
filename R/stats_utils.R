## ======================================================================
## R/stats_utils.R
## Paired statistics for benchmark tables.
##
## Reviewer comments: 9 (paired mean differences with CIs, explicit comparison family, SE on precision), minor 5
## Revision notes: paired_table() gives mean paired difference + 95% bootstrap CI, Wilcoxon, BH within a named family; summary_table() gives mean +/- SE.
## How to run: Tested by tests/test_metrics.R.
## Full write-up: docs/REVISION_REPORT.md
## ======================================================================

# ======================================================================
# Paired-comparison statistics for benchmark tables
# ----------------------------------------------------------------------
# For every (cell, competitor, metric) the reported quantities are:
#   mean_diff      mean of per-replicate (ASCEND - competitor)
#   ci_lo, ci_hi   95% CI of the mean paired difference (percentile
#                  bootstrap over replicates; t-interval also returned)
#   median_diff    for reference only
#   wins/ties/losses
#   p_wilcox       paired Wilcoxon signed-rank p-value
#   q_bh           Benjamini-Hochberg q within the stated family
#
# The comparison family is explicit: by default one family per metric,
# containing every (cell x competitor) comparison of that metric. Pass
# `family_by` to change it and the column `family` records what was used.
#
# Dependencies: base R + stats only.
# ======================================================================

mean_se <- function(x) {
  x <- x[is.finite(x)]
  c(mean = if (length(x)) mean(x) else NA_real_,
    se   = if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_,
    n    = length(x))
}

boot_ci_mean <- function(d, B = 5000L, level = 0.95, seed = 1L) {
  d <- d[is.finite(d)]
  if (length(d) < 2) return(c(NA_real_, NA_real_))
  if (all(d == d[1])) return(c(d[1], d[1]))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  bm <- vapply(seq_len(B), function(i) mean(sample(d, replace = TRUE)), numeric(1))
  a  <- (1 - level) / 2
  unname(quantile(bm, c(a, 1 - a)))
}

paired_compare <- function(x, y, alternative = "two.sided", B = 5000L) {
  ok <- is.finite(x) & is.finite(y)
  d  <- x[ok] - y[ok]
  n  <- length(d)
  out <- list(n_pairs = n, mean_diff = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
              t_lo = NA_real_, t_hi = NA_real_, median_diff = NA_real_,
              wins = sum(d > 0), ties = sum(d == 0), losses = sum(d < 0),
              p_wilcox = NA_real_)
  if (n < 2) return(out)
  out$mean_diff   <- mean(d)
  out$median_diff <- median(d)
  ci <- boot_ci_mean(d, B = B); out$ci_lo <- ci[1]; out$ci_hi <- ci[2]
  # t.test refuses differences that are constant up to rounding (sd > 0 but
  # tiny next to the mean, e.g. two methods tying in every replicate)
  tt <- if (sd(d) > 0) tryCatch(t.test(d), error = function(e) NULL) else NULL
  if (!is.null(tt)) { out$t_lo <- tt$conf.int[1]; out$t_hi <- tt$conf.int[2] }
  else out$t_lo <- out$t_hi <- mean(d)
  out$p_wilcox <- if (any(d != 0))
    suppressWarnings(wilcox.test(x[ok], y[ok], paired = TRUE,
                                 alternative = alternative, exact = FALSE)$p.value)
  else 1
  out
}

# Long-format benchmark table -> paired comparisons of `ref` vs every other method.
#   dt        data.frame with columns method, <cell_cols>, <rep_col>, <metrics>
#   metrics   character vector of metric columns
#   family_by columns defining a BH family (default: metric only)
paired_table <- function(dt, metrics, cell_cols, rep_col = "rep", ref = "ASCEND",
                         alternative = "two.sided", family_by = "metric", B = 5000L) {
  dt <- as.data.frame(dt)
  comps <- setdiff(unique(dt$method), ref)
  cells <- unique(dt[, cell_cols, drop = FALSE])
  rows <- list()
  for (ci in seq_len(nrow(cells))) {
    key <- cells[ci, , drop = FALSE]
    sel <- Reduce(`&`, lapply(cell_cols, function(cc) dt[[cc]] == key[[cc]]))
    sub <- dt[sel, , drop = FALSE]
    r   <- sub[sub$method == ref, , drop = FALSE]
    for (cm in comps) {
      c_ <- sub[sub$method == cm, , drop = FALSE]
      m  <- merge(r, c_, by = rep_col, suffixes = c(".ref", ".cmp"))
      for (mt in metrics) {
        pc <- paired_compare(m[[paste0(mt, ".ref")]], m[[paste0(mt, ".cmp")]],
                             alternative = alternative, B = B)
        rows[[length(rows) + 1]] <- cbind(key, data.frame(
          metric = mt, reference = ref, competitor = cm, as.data.frame(pc),
          stringsAsFactors = FALSE))
      }
    }
  }
  out <- do.call(rbind, rows)
  fam <- interaction(out[, family_by, drop = FALSE], drop = TRUE, sep = "|")
  out$family <- paste0(paste(family_by, collapse = "+"), "=", as.character(fam))
  out$family_size <- ave(seq_len(nrow(out)), fam, FUN = length)
  out$q_bh <- ave(out$p_wilcox, fam, FUN = function(p) p.adjust(p, "BH"))
  out
}

# Per-cell mean +/- SE for every metric and method (used for Tables 1-2)
summary_table <- function(dt, metrics, cell_cols) {
  dt <- as.data.frame(dt)
  grp <- unique(dt[, c(cell_cols, "method"), drop = FALSE])
  rows <- lapply(seq_len(nrow(grp)), function(i) {
    key <- grp[i, , drop = FALSE]
    sel <- Reduce(`&`, lapply(names(key), function(cc) dt[[cc]] == key[[cc]]))
    vals <- unlist(lapply(metrics, function(mt) {
      s <- mean_se(dt[sel, mt]); setNames(s, paste0(mt, c("_mean", "_se", "_n")))
    }))
    cbind(key, as.data.frame(as.list(vals)))
  })
  do.call(rbind, rows)
}

fmt_mean_se <- function(m, se, digits = 3) {
  ifelse(is.na(m), "NA",
         ifelse(is.na(se), formatC(m, digits = digits, format = "f"),
                sprintf(paste0("%.", digits, "f ± %.", digits, "f"), m, se)))
}

fmt_diff_ci <- function(md, lo, hi, q, digits = 3) {
  star <- ifelse(is.na(q), "", ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", ""))))
  sprintf(paste0("%+.", digits, "f [%+.", digits, "f, %+.", digits, "f]%s"), md, lo, hi, star)
}
