## ======================================================================
## analysis/make_report_figures.R
## Draws every figure in docs/REVISION_REPORT.md from the CSV outputs of
## the analysis and benchmark scripts. No simulation happens here.
##
## Reviewer comments: see each figure below.
## Inputs: for each file, analysis/out/<file> if it exists (your own run),
##   otherwise docs/pilot_results/<file> (the pilot results in the report).
##   The ASCEND vs CBL figure reads benchmarks/cbl/results/ first.
## Run: Rscript analysis/make_report_figures.R   (ggplot2, data.table; ~10 s)
## Writes docs/figures/fig*.png
## ======================================================================
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
FIG <- file.path(ROOT, "docs", "figures"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
inp <- function(f, extra = NULL) {
  p <- Filter(file.exists, c(extra, file.path(ROOT, "analysis", "out", f), file.path(ROOT, "docs", "pilot_results", f)))
  if (!length(p)) { message("skip: ", f, " not found"); return(NULL) }
  message("read ", p[1]); fread(p[1])
}
save <- function(p, name, w = 8, h = 4.5) { ggsave(file.path(FIG, name), p, width = w, height = h, dpi = 150); message("wrote ", name) }
th <- theme_bw(base_size = 10) + theme(legend.position = "bottom")
se <- function(x) { x <- x[is.finite(x)]; if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_ }

## Fig 1. Scoring audit (section 1; comments 1, 2) --------------------------
a <- inp("scoring_audit_grn.csv")
if (!is.null(a)) {
  l <- melt(a, id.vars = c("sp", "r2", "seed"),
            measure.vars = c("f1_submitted", "f1_fixed", "direction_acc_submitted", "orient_acc_fixed"))
  l[, metric := ifelse(grepl("^f1", variable), "ASCEND F1 (GRN benchmark)", "ASCEND orientation accuracy")]
  l[, scoring := ifelse(grepl("submitted", variable), "submitted scoring", "fixed scoring")]
  s <- l[, .(mean = mean(value, na.rm = TRUE), se = se(value)), by = .(sp, r2, metric, scoring)]
  s[, cell := sprintf("sp=%.1f\nR2=%.1f", sp, r2)]
  save(ggplot(s, aes(cell, mean, fill = scoring)) +
         geom_col(position = position_dodge(0.8), width = 0.75) +
         geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.8), width = 0.2) +
         facet_wrap(~metric) + scale_fill_manual(values = c("#1b9e77", "#bbbbbb")) +
         labs(x = NULL, y = "Mean over replicates (+/- SE)", fill = NULL,
              title = "Same ASCEND output, scored with the submitted and the fixed functions",
              subtitle = "n = 2000, d_z = 20, d_x = 15, submitted simulator") + th,
       "fig1_scoring_audit.png", 9, 4.5)
}
b <- inp("scoring_audit_grid_status.csv")
if (!is.null(b)) {
  b <- melt(b, id.vars = "method", measure.vars = intersect(c("ok", "error", "timeout", "skipped_timeout"), names(b)))
  save(ggplot(b, aes(method, value, fill = variable)) + geom_col(position = "fill") +
         scale_fill_manual(values = c(ok = "#1b9e77", error = "#d95f02", timeout = "#7570b3", skipped_timeout = "#bbbbbb")) +
         labs(x = NULL, y = "Share of runs", fill = NULL,
              title = "Submitted causal grid (Fig. 6): run status by method",
              subtitle = "GES/PC errors: crash on undirected CPDAG edges, now fixed") + th,
       "fig1b_grid_run_status.png", 6, 4)
}

## Fig 2. Simulator diagnostic (section 2) -----------------------------------
d <- inp("simulator_diagnostic.csv")
if (!is.null(d)) {
  d[, max_var_z := as.numeric(max_var_z)]
  s <- d[, .(max_var_z = log10(max(max_var_z)), pcor = mean(median_edge_pcor, na.rm = TRUE),
             tp = mean(dir_tp, na.rm = TRUE)), by = .(d_z, simulator)]
  l <- melt(s, id.vars = c("d_z", "simulator"))
  l[, variable := factor(variable, c("max_var_z", "pcor", "tp"),
                         c("log10 of the largest background variance", "Median |partial cor| of true X->X edges given Z",
                           "ASCEND true positives (directed)"))]
  save(ggplot(l, aes(d_z, value, colour = simulator)) + geom_line() + geom_point() +
         facet_wrap(~variable, scales = "free_y") + scale_x_log10() +
         labs(x = "d_z (background variables)", y = NULL, colour = NULL,
              title = "Submitted simulator: background variance explodes and the foreground signal vanishes",
              subtitle = "n = 2000, d_x = 20; same seeds for both simulators") + th,
       "fig2_simulator.png", 10, 4)
}
g <- inp("submitted_grid_tp_by_dz.csv")
if (!is.null(g)) {
  save(ggplot(g, aes(factor(dz_over_dx), tp, fill = factor(d_x))) + geom_col(position = "dodge") +
         labs(x = "d_z / d_x", y = "Mean ASCEND true positives", fill = "d_x",
              title = "Submitted 81-cell grid: ASCEND finds almost nothing once d_z >= 2 d_x") + th,
       "fig2b_submitted_grid_tp.png", 6, 4)
}

## Fig 3. Cost and accuracy vs d_z, prescreen (sections 4, 7; comments 5, 6) --
p <- inp("prescreen_scaling.csv")
if (!is.null(p)) {
  p[, prescreen := factor(ifelse(prescreen < 1, "prescreen 0.30 (default)", "prescreen off"))]
  s <- p[, .(`Wall time (s)` = mean(time_sec, na.rm = TRUE), `CI tests` = mean(n_ci_tests, na.rm = TRUE),
             `Time per CI test (microseconds)` = 1e6 * mean(time_per_test, na.rm = TRUE),
             `Mean |S| in pairwise tests` = mean(mean_pair_cond, na.rm = TRUE),
             `Directed F1` = mean(dir_f1, na.rm = TRUE), `Orientation accuracy` = mean(orient_acc, na.rm = TRUE)),
         by = .(n, d_z, prescreen)]
  l <- melt(s, id.vars = c("n", "d_z", "prescreen"))
  l[, n := factor(paste0("n = ", n))]
  save(ggplot(l, aes(d_z, value, colour = n, linetype = prescreen)) + geom_line() + geom_point(size = 1) +
         facet_wrap(~variable, scales = "free_y", ncol = 3) + scale_x_log10() +
         labs(x = "d_z", y = NULL, colour = NULL, linetype = NULL,
              title = "ASCEND (guarded conditioning set): cost and accuracy as the background grows",
              subtitle = "d_x = 10, fixed simulator") + th,
       "fig3_scaling_prescreen.png", 10, 6)
}

## Fig 4. ASCEND vs CBL cost per test (section 4; comment 6) -----------------
c4 <- inp("results_ascend_vs_cbl_v2.csv", file.path(ROOT, "benchmarks", "cbl", "results", "results_ascend_vs_cbl_v2.csv"))
if (!is.null(c4)) {
  c4 <- c4[status == "ok"]
  sw <- rbind(c4[n == 1024 & d_x == 5, .(sweep = "d_z (n = 1024, d_x = 5)", x = d_z, method, time_sec, n_ci_tests, time_per_test, dir_f1, coverage)],
              c4[d_z == 10 & d_x == 5, .(sweep = "n (d_x = 5, d_z = 10)", x = n, method, time_sec, n_ci_tests, time_per_test, dir_f1, coverage)],
              c4[n == 1024 & d_z == 10, .(sweep = "d_x (n = 1024, d_z = 10)", x = d_x, method, time_sec, n_ci_tests, time_per_test, dir_f1, coverage)])
  s <- sw[, .(`Wall time (s)` = mean(time_sec), `Tests (CBL: l0 verdicts)` = mean(n_ci_tests),
              `Time per test (ms)` = 1e3 * mean(time_per_test, na.rm = TRUE), Coverage = mean(coverage, na.rm = TRUE)),
          by = .(sweep, x, method)]
  l <- melt(s, id.vars = c("sweep", "x", "method"))
  save(ggplot(l, aes(x, value, colour = method)) + geom_line() + geom_point() +
         facet_grid(variable ~ sweep, scales = "free") + scale_y_log10() +
         scale_colour_manual(values = c(ascend = "#1b9e77", cbl = "#d95f02")) +
         labs(x = NULL, y = NULL, colour = NULL, title = "ASCEND vs CBL: work, cost per test and coverage") + th,
       "fig4_ascend_vs_cbl_cost.png", 9, 8)
}

## Fig 5. Tier sensitivity (section 5; comment 4) -----------------------------
t5 <- inp("tier_sensitivity.csv")
if (!is.null(t5)) {
  t5 <- t5[status == "ok"]
  l <- melt(t5, id.vars = c("type", "level", "seed"), measure.vars = c("dir_precision", "dir_recall", "dir_f1", "orient_acc"))
  s <- l[, .(mean = mean(value, na.rm = TRUE), se = se(value)), by = .(type, level, variable)]
  s[, type := factor(type, c("x_as_z", "feedback", "hide_z", "z_as_x"),
                     c("Foreground labelled as background", "Feedback X -> Z", "Unmeasured background",
                       "Background labelled as foreground"))]
  save(ggplot(s, aes(level, mean, colour = variable)) + geom_line() + geom_point() +
         geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.01) +
         facet_wrap(~type, scales = "free_x") + coord_cartesian(ylim = c(0, 1)) +
         labs(x = "Violation rate", y = "Mean (+/- SE)", colour = NULL,
              title = "ASCEND under violations of the two-tier assumption",
              subtitle = "n = 1024, d_z = 50, d_x = 15, sp = 0.7") + th,
       "fig5_tier_sensitivity.png", 9, 6)
}

## Fig 6. Orientation (section 6; comment 11) ---------------------------------
r6 <- inp("orientation_by_rule.csv")
if (!is.null(r6)) {
  r6 <- r6[rule != "unresolved"]
  l <- melt(r6, id.vars = "rule", measure.vars = c("share", "correct_state", "orient_correct_when_related"))
  l[, variable := factor(variable, c("share", "correct_state", "orient_correct_when_related"),
                         c("Share of all pairs committed", "Committed state correct",
                           "Orientation correct (truly related pairs)"))]
  save(ggplot(l[is.finite(value)], aes(rule, value, fill = rule)) + geom_col(show.legend = FALSE) +
         geom_hline(data = data.frame(variable = factor("Orientation correct (truly related pairs)", levels(l$variable)), y = 0.5),
                    aes(yintercept = y), linetype = 2) +
         facet_wrap(~variable) + labs(x = NULL, y = NULL, title = "Which rule commits each pair, and how often it is right",
                                      subtitle = "Dashed line: a coin flip for orientation") + th,
       "fig6_orientation_rules.png", 9, 3.8)
}
f6 <- inp("orientation_by_feature.csv")
if (!is.null(f6)) {
  f6 <- f6[feature != "committing rule"]
  f6[, level := factor(level, unique(level))]
  save(ggplot(f6, aes(level, orient_acc)) + geom_col(fill = "#7570b3") +
         geom_text(aes(label = pairs), vjust = -0.3, size = 2.5) +
         facet_wrap(~feature, scales = "free_x", ncol = 4) + coord_cartesian(ylim = c(0, 1)) +
         labs(x = NULL, y = "Orientation accuracy", title = "Orientation accuracy by graph feature (numbers: pairs)") + th,
       "fig6b_orientation_features.png", 11, 5)
}
o6 <- inp("orientation_outcomes.csv")
if (!is.null(o6)) {
  save(ggplot(o6, aes(factor(sp), Freq, fill = outcome)) + geom_col() + facet_wrap(~paste0("n = ", n)) +
         labs(x = "sp (sparsity)", y = "Share of truly related pairs", fill = NULL,
              title = "What ASCEND does with truly related pairs") + th + guides(fill = guide_legend(nrow = 2)),
       "fig6c_orientation_outcomes.png", 7, 4.5)
}

## Fig 7. Runtime vs discovered edges (minor 6), submitted grid -------------
cd <- fread(Filter(file.exists, c(file.path(ROOT, "benchmarks", "causal_grid", "ascend_benchmark_v3_merged.csv"),
                                  file.path(ROOT, "benchmarks", "causal_grid", "results_submitted", "ascend_benchmark_merged.csv")))[1])
cd <- cd[status == "ok" & is.finite(elapsed_sec) & elapsed_sec > 0]
cd[, edges := if ("dir_claims" %in% names(cd)) dir_claims else tp + fp]
save(ggplot(cd, aes(pmax(edges, 0.5), elapsed_sec, colour = method)) + geom_point(alpha = 0.15, size = 0.5) +
       geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + scale_x_log10() + scale_y_log10() +
       labs(x = "Discovered (claimed) ancestral edges", y = "Wall-clock time (s)", colour = NULL,
            title = "Runtime against discovered edges (submitted grid until re-run)") + th,
     "fig7_runtime_vs_edges.png", 6, 4.5)
