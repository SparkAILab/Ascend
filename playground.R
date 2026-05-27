## Statistical analysis of sim_sweep_results.csv.
##
## For each n level, runs paired Wilcoxon signed-rank tests of
## ASCEND-PC vs each baseline on three headline metrics:
##   - AUPRC ratio
##   - Direction accuracy
##   - F1 score
##
## Reports:
##   - Per-rep distributions (violin plots, saved as PDF/PNG)
##   - Wilcoxon p-values, paired by rep within each n
##   - Cliff's delta (non-parametric effect size)
##   - Bootstrap 95% CIs on means
##   - Win-rate table (how many reps ASCEND-PC beat each baseline)
##
## Outputs:
##   sweep_stats.csv          per-(n, metric, comparison) test stats
##   sweep_winrates.csv       per-(n, metric) win counts
##   sweep_distributions.{pdf,png}  violin/box figure

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

if (!file.exists("sim_sweep_results.csv"))
  stop("sim_sweep_results.csv not found in working directory")

dat <- read.csv("sim_sweep_results.csv", stringsAsFactors = FALSE)

# Sanity check.
cat(sprintf("Loaded %d rows: %d n-levels x %d methods x %d reps\n",
            nrow(dat),
            length(unique(dat$n)),
            length(unique(dat$method)),
            length(unique(dat$rep))))

dat$method <- factor(dat$method,
                     levels = c("Plain PC (X)", "Plain PC (Z+X)", "ASCEND-PC"))

## --- Cliff's delta -------------------------------------------------------
##
## Non-parametric effect size for two paired samples. Range [-1, 1].
##   |d| < 0.147 : negligible
##   |d| < 0.33  : small
##   |d| < 0.474 : medium
##   |d| >= 0.474: large

cliffs_delta <- function(x, y) {
  # x and y need not be the same length; treat as independent samples.
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) return(NA)
  more <- sum(outer(x, y, ">"))
  less <- sum(outer(x, y, "<"))
  (more - less) / (length(x) * length(y))
}

cliffs_label <- function(d) {
  if (is.na(d)) return("NA")
  a <- abs(d)
  if (a < 0.147)  "negligible"
  else if (a < 0.33)   "small"
  else if (a < 0.474)  "medium"
  else                 "large"
}

## --- Paired comparison runner -------------------------------------------

run_paired_test <- function(df, metric, baseline) {
  # Pair on rep id within a single n level.
  w <- df %>%
    select(rep, method, all_of(metric)) %>%
    pivot_wider(names_from = method, values_from = all_of(metric))
  asc  <- w[["ASCEND-PC"]]
  base <- w[[baseline]]
  
  keep <- !is.na(asc) & !is.na(base)
  asc  <- asc[keep]; base <- base[keep]
  if (length(asc) < 3) return(NULL)
  
  diff <- asc - base
  test <- tryCatch(
    wilcox.test(asc, base, paired = TRUE, alternative = "greater",
                exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(test)) return(NULL)
  
  cd <- cliffs_delta(asc, base)
  
  list(
    n_pairs    = length(asc),
    asc_mean   = mean(asc),
    asc_sd     = sd(asc),
    base_mean  = mean(base),
    base_sd    = sd(base),
    diff_mean  = mean(diff),
    diff_sd    = sd(diff),
    win_count  = sum(diff > 0),
    tie_count  = sum(diff == 0),
    loss_count = sum(diff < 0),
    wilcox_p   = test$p.value,
    wilcox_W   = unname(test$statistic),
    cliffs_d   = cd,
    cliffs_lab = cliffs_label(cd)
  )
}

## --- Run all comparisons -------------------------------------------------

metrics <- c("auprc_ratio", "dir_acc", "f1")
baselines <- c("Plain PC (X)", "Plain PC (Z+X)")
n_levels <- sort(unique(dat$n))

stats_rows <- list()
for (n_val in n_levels) {
  sub <- dat %>% filter(n == n_val)
  for (m in metrics) {
    for (b in baselines) {
      r <- run_paired_test(sub, m, b)
      if (is.null(r)) next
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        n           = n_val,
        metric      = m,
        comparison  = sprintf("ASCEND-PC vs %s", b),
        n_pairs     = r$n_pairs,
        ascend_mean = round(r$asc_mean, 3),
        ascend_sd   = round(r$asc_sd, 3),
        baseline_mean = round(r$base_mean, 3),
        baseline_sd   = round(r$base_sd, 3),
        diff_mean   = round(r$diff_mean, 3),
        diff_sd     = round(r$diff_sd, 3),
        wins        = r$win_count,
        ties        = r$tie_count,
        losses      = r$loss_count,
        wilcox_p    = signif(r$wilcox_p, 3),
        cliffs_d    = round(r$cliffs_d, 3),
        cliffs_lab  = r$cliffs_lab,
        stringsAsFactors = FALSE
      )
    }
  }
}
stats_df <- bind_rows(stats_rows)

# Add a holm correction over the 24 tests (4 n x 3 metrics x 2 baselines).
stats_df$wilcox_p_holm <- signif(p.adjust(stats_df$wilcox_p, method = "holm"), 3)

write.csv(stats_df, "sweep_stats.csv", row.names = FALSE)
cat(sprintf("\nWrote sweep_stats.csv (%d rows)\n", nrow(stats_df)))

# Print a focused view: just direction accuracy, since that's the headline.
cat("\n=== Direction accuracy (headline metric) ===\n")
da <- stats_df %>% filter(metric == "dir_acc") %>%
  select(n, comparison, ascend_mean, baseline_mean, diff_mean,
         wins, ties, losses, wilcox_p, wilcox_p_holm, cliffs_d, cliffs_lab)
print(as.data.frame(da), row.names = FALSE)

cat("\n=== AUPRC ratio ===\n")
ar <- stats_df %>% filter(metric == "auprc_ratio") %>%
  select(n, comparison, ascend_mean, baseline_mean, diff_mean,
         wins, ties, losses, wilcox_p, wilcox_p_holm, cliffs_d, cliffs_lab)
print(as.data.frame(ar), row.names = FALSE)

cat("\n=== F1 ===\n")
f1 <- stats_df %>% filter(metric == "f1") %>%
  select(n, comparison, ascend_mean, baseline_mean, diff_mean,
         wins, ties, losses, wilcox_p, wilcox_p_holm, cliffs_d, cliffs_lab)
print(as.data.frame(f1), row.names = FALSE)

## --- Bootstrap CIs on means ---------------------------------------------

boot_ci <- function(x, B = 1000, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(c(NA, NA))
  set.seed(42)
  boot_means <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  quantile(boot_means, probs = c((1 - conf) / 2, 1 - (1 - conf) / 2))
}

ci_rows <- list()
for (n_val in n_levels) {
  for (meth in levels(dat$method)) {
    for (m in metrics) {
      vals <- dat %>% filter(n == n_val, method == meth) %>% pull(.data[[m]])
      ci <- boot_ci(vals)
      ci_rows[[length(ci_rows) + 1]] <- data.frame(
        n = n_val, method = meth, metric = m,
        mean = round(mean(vals, na.rm = TRUE), 3),
        ci_lo = round(unname(ci[1]), 3),
        ci_hi = round(unname(ci[2]), 3),
        stringsAsFactors = FALSE
      )
    }
  }
}
ci_df <- bind_rows(ci_rows)
write.csv(ci_df, "sweep_ci.csv", row.names = FALSE)
cat("\nWrote sweep_ci.csv\n")

## --- Win-rate table ------------------------------------------------------

winrates <- stats_df %>%
  select(n, metric, comparison, wins, ties, losses) %>%
  arrange(metric, n, comparison)
write.csv(winrates, "sweep_winrates.csv", row.names = FALSE)
cat("Wrote sweep_winrates.csv\n")

## --- Figure: distributions ----------------------------------------------

theme_pub <- function(base = 10) {
  theme_bw(base_size = base) +
    theme(
      plot.title = element_text(face = "bold", size = base + 1,
                                colour = "#1A237E"),
      plot.subtitle = element_text(size = base - 1, colour = "#546E7A"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.35),
      panel.border     = element_rect(colour = "grey70", linewidth = 0.5),
      panel.background = element_rect(fill = "#F8F9FA"),
      strip.background = element_rect(fill = "#E3F2FD", colour = "grey70"),
      strip.text       = element_text(face = "bold", size = base - 1),
      legend.position  = "top",
      legend.title     = element_blank()
    )
}

cols <- c("Plain PC (X)" = "#C62828",
          "Plain PC (Z+X)" = "#EF9A9A",
          "ASCEND-PC" = "#1565C0")

dat_long <- dat %>%
  pivot_longer(c(auprc_ratio, dir_acc, f1),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         auprc_ratio = "AUPRC ratio",
                         dir_acc = "Direction accuracy",
                         f1 = "F1"),
         metric = factor(metric,
                         levels = c("AUPRC ratio", "Direction accuracy", "F1")),
         n_lab = sprintf("n = %d", n),
         n_lab = factor(n_lab, levels = paste0("n = ", n_levels)))

# Per-metric panel with all n on x-axis, method colour.
build_panel <- function(metric_name) {
  d <- dat_long %>% filter(metric == metric_name)
  ggplot(d, aes(x = n_lab, y = value, fill = method, colour = method)) +
    geom_violin(alpha = 0.30, position = position_dodge(width = 0.8),
                width = 0.7, scale = "width", colour = NA) +
    geom_boxplot(alpha = 0.85, position = position_dodge(width = 0.8),
                 width = 0.25, outlier.shape = NA, colour = "black",
                 linewidth = 0.4) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                               dodge.width = 0.8),
               size = 1.2, alpha = 0.75) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    labs(title = metric_name, x = NULL, y = metric_name) +
    theme_pub() +
    theme(axis.text.x = element_text(size = 9))
}

p_auprc <- build_panel("AUPRC ratio") +
  geom_hline(yintercept = 1.0, linetype = "22",
             colour = "grey50", linewidth = 0.5)
p_dir   <- build_panel("Direction accuracy") +
  geom_hline(yintercept = 0.5, linetype = "22",
             colour = "grey50", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 1),
                     labels = percent_format(accuracy = 1))
p_f1    <- build_panel("F1")

fig <- (p_auprc / p_dir / p_f1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Per-replicate distributions across sample sizes",
    subtitle = sprintf(
      "10 replicates per cell; d_x = %d, d_z = %d. Dashed line: random baseline.",
      20, 40),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, colour = "#0D47A1"),
      plot.subtitle = element_text(size = 9, colour = "#546E7A"),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(12, 12, 12, 12),
      legend.position = "top"
    )
  ) &
  theme(legend.position = "top")

ggsave("sweep_distributions.pdf", fig,
       width = 11, height = 12, device = cairo_pdf, bg = "white")
ggsave("sweep_distributions.png", fig,
       width = 11, height = 12, dpi = 300, bg = "white")
cat("Wrote sweep_distributions.{pdf,png}\n")

cat("\nDone.\n")