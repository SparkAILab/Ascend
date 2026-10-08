## ======================================================================
## benchmarks/grn_multi/04_summarise.R
## Turns the scored runs into the tables and figures for the paper and the
## response letter. Reads WORK/results/{runs.csv, metrics_long.csv.gz}.
##
## Reviewer comments: 1 (all baselines), 2 (every metric in every cell),
##   3 (SERGIO arm), 9 (mean +/- SE for every metric; paired mean
##   differences with 95% bootstrap CIs, wins/ties/losses, Wilcoxon, BH
##   within an explicit family), 11 (orientation), minor 6 (run time)
## How to run: Rscript benchmarks/grn_multi/04_summarise.R  (a few minutes)
## Output (WORK/results/):
##   summary.csv        mean, SE, n_reps per arm x cell x method x target x metric
##   paired.csv         ASCEND minus each competitor, per cell; BH family =
##                      arm x target x metric (all cells and competitors)
##   ranks.csv          ASCEND's rank among all methods per cell and metric
##   status.csv         ok / error / timeout counts per arm x method
##   runtime.csv        run-time summary per arm x method
##   TABLES.md          the main tables, ready to read
##   fig_*.pdf / .png   heatmaps (method x cell), orientation, run time
## Two comparison metrics are added for methods of different classes:
##   dirF1  = dir_f1 (causal: own graph) or dir_f1_Kasc (directed GRN:
##            top-K_asc pairs oriented by the larger score)
##   orient = orient_acc or orient_acc_Kasc in the same way
## ======================================================================

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  setwd(d)
})
source("benchmarks/grn_multi/config.R")
source(file.path(ROOT, "R", "stats_utils.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

runs <- fread(file.path(RES_DIR, "runs.csv"))
mets <- fread(file.path(RES_DIR, "metrics_long.csv.gz"))
cellcols <- c("arm", "n", "sp", "r2", "ds", "noise")
info <- unique(runs[, c("dataset", cellcols, "rep", "d_x", "d_z"), with = FALSE])
## SERGIO cells are ds x noise; n varies by a few cells per noisy replicate
## (cells dropped by the noise model) and must not split the cell
info[arm == "sergio", n := NA]
lab  <- setNames(METHODS$label, METHODS$method)
cls  <- setNames(METHODS$class, METHODS$method)

## derived cross-class metrics
w <- dcast(mets, dataset + method + target ~ metric, value.var = "value")
num <- function(x) if (is.null(x)) NA_real_ else x
w[, dirF1  := fifelse(!is.na(num(w$dir_f1)), num(w$dir_f1), num(w$dir_f1_Kasc))]
w[, orient := fifelse(!is.na(num(w$orient_acc)), num(w$orient_acc), num(w$orient_acc_Kasc))]
mets <- melt(w, id.vars = c("dataset", "method", "target"), variable.name = "metric",
             value.name = "value", na.rm = TRUE)
mets <- merge(mets, info, by = "dataset")
mets[, label := lab[method]]
rt <- merge(runs[status == "ok", .(dataset, method, value = runtime_s)], info, by = "dataset")
rt[, `:=`(target = "anc", metric = "runtime_s", label = lab[method])]
rt2 <- copy(rt)[, target := "direct"]
mets <- rbind(mets, rt, rt2, use.names = TRUE)
cellkey <- function(d) d[, cell := fifelse(arm == "sergio", paste(ds, noise),
                                           sprintf("n=%d sp=%s r2=%s", n, sp, r2))]
cellkey(mets)

## 1. mean +/- SE
summ <- mets[, .(mean = mean(value), se = if (.N > 1) sd(value) / sqrt(.N) else NA_real_, n_reps = .N),
             by = c(cellcols, "cell", "method", "label", "target", "metric")]
fwrite(summ, file.path(RES_DIR, "summary.csv"))

## 2. paired differences vs ASCEND
PAIR_METRICS <- c("aupr", "auroc", "aupr_ratio", "epr", "f1_Kasc", "prec_Kasc", "rec_Kasc",
                  "f1_Ktrue", "dir_aupr", "dir_auroc", "orient_acc_rank", "dirF1", "orient", "runtime_s")
pm <- mets[metric %in% PAIR_METRICS]
ref <- pm[method == "ascend", .(dataset, target, metric, ref = value)]
cmp <- merge(pm[method != "ascend"], ref, by = c("dataset", "target", "metric"))
set.seed(1)
paired <- cmp[, {
  pc <- paired_compare(ref, value, B = 2000L)
  as.list(unlist(pc))
}, by = c(cellcols, "cell", "target", "metric", "method", "label")]
paired[, family := paste(arm, target, metric, sep = "|")]
paired[, q_bh := p.adjust(p_wilcox, "BH"), by = family]
paired[, family_size := .N, by = family]
setnames(paired, c("method", "label"), c("competitor", "competitor_label"))
fwrite(paired, file.path(RES_DIR, "paired.csv"))

## 3. ASCEND's rank per cell (1 = best; lower is better only for runtime and shd)
lower_better <- c("runtime_s", "shd", "shd_Kasc")
ranks <- summ[metric %in% c(PAIR_METRICS, "dir_f1", "orient_acc")]
ranks[, rank := frank(if (metric[1] %in% lower_better) mean else -mean, ties.method = "min"),
      by = c(cellcols, "cell", "target", "metric")]
ranks[, n_methods := .N, by = c(cellcols, "cell", "target", "metric")]
fwrite(ranks[method == "ascend", c(cellcols, "cell", "target", "metric", "mean", "rank", "n_methods"), with = FALSE],
       file.path(RES_DIR, "ranks.csv"))

## 4. status and run time
status <- dcast(merge(runs, unique(runs[, .(dataset)]), by = "dataset"),
                arm + method ~ status, fun.aggregate = length, value.var = "dataset")
fwrite(status, file.path(RES_DIR, "status.csv"))
runtime <- runs[status == "ok", .(median_s = median(runtime_s), mean_s = mean(runtime_s),
                                  max_s = max(runtime_s), cpu_hours = sum(cpu_s, na.rm = TRUE) / 3600),
                by = .(arm, method)]
fwrite(runtime, file.path(RES_DIR, "runtime.csv"))

## 5. TABLES.md
f3 <- function(m, s) ifelse(is.na(m), "", ifelse(is.na(s), sprintf("%.3f", m), sprintf("%.3f ± %.3f", m, s)))
md_table <- function(d, rowcol, colcol, valcol) {
  wide <- dcast(d, as.formula(paste(rowcol, "~", colcol)), value.var = valcol)
  hdr <- paste0("| ", paste(names(wide), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(wide)), collapse = "|"), "|")
  body <- apply(wide, 1, function(r) paste0("| ", paste(ifelse(is.na(r), "", r), collapse = " | "), " |"))
  c(hdr, sep, body, "")
}
order_methods <- function(d) d[, label := factor(label, levels = METHODS$label)]
out <- c("# Multi-method benchmark: tables", "",
         sprintf("Generated %s from %d runs (%d ok). Metrics: see R/scoring.R. Values: mean ± SE over replicates.",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(runs), sum(runs$status == "ok")), "")
KEY <- c("aupr", "auroc", "epr", "f1_Kasc", "dirF1", "orient", "runtime_s")
pc <- PRIMARY
for (tg in c("anc", "direct")) {
  d <- summ[arm == pc$arm & n == pc$n & sp == pc$sp & r2 == pc$r2 & target == tg & metric %in% KEY]
  if (!nrow(d)) next
  d[, v := f3(mean, se)]; order_methods(d)
  out <- c(out, sprintf("## Primary cell (%s, n=%d, sp=%s, R²=%s), truth = %s", pc$arm, pc$n, pc$sp, pc$r2,
                        if (tg == "anc") "ancestral relations" else "direct edges"), "",
           md_table(d, "label", "metric", "v"))
}
for (a in names(ARMS)) for (tg in c("anc", "direct")) {
  d <- mets[arm == a & target == tg & metric %in% KEY, .(m = mean(value), s = sd(value) / sqrt(.N)),
            by = .(label, metric)]
  if (!nrow(d)) next
  d[, v := f3(m, s)]; order_methods(d)
  out <- c(out, sprintf("## Arm %s, all cells pooled, truth = %s", a, tg), "", md_table(d, "label", "metric", "v"))
}
for (mt in c("f1_Kasc", "aupr", "dirF1")) {
  d <- paired[arm == "lin" & target == "anc" & metric == mt]
  if (!nrow(d)) next
  d[, v := fmt_diff_ci(mean_diff, ci_lo, ci_hi, q_bh)]
  d[, competitor_label := factor(competitor_label, levels = METHODS$label)]
  out <- c(out, sprintf("## Paired %s, ASCEND minus competitor, arm lin, ancestral truth", mt),
           "Mean difference [95% bootstrap CI]; * q<0.05, ** q<0.01, *** q<0.001 (Wilcoxon, BH within arm × target × metric).", "",
           md_table(d, "cell", "competitor_label", "v"))
}
st <- copy(status); out <- c(out, "## Run status", "", md_table(melt(st, id.vars = c("arm", "method"))[, v := as.character(value)], "method", "arm+variable", "v"))
writeLines(out, file.path(RES_DIR, "TABLES.md"))

## 6. figures
theme_set(theme_bw(base_size = 9))
plot_heat <- function(mt, tg, file) {
  d <- summ[metric == mt & target == tg]
  if (!nrow(d)) return(invisible())
  order_methods(d)
  g <- ggplot(d, aes(cell, label, fill = mean)) + geom_tile() +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 2.2) +
    scale_fill_viridis_c(name = mt) + facet_grid(~ arm, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, title = sprintf("%s (truth: %s)", mt, tg)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  for (ext in c("pdf", "png")) ggsave(file.path(RES_DIR, paste0(file, ".", ext)), g, width = 12, height = 4.5, dpi = 200)
}
for (tg in c("anc", "direct")) {
  plot_heat("f1_Kasc", tg, paste0("fig_f1K_", tg))
  plot_heat("aupr", tg, paste0("fig_aupr_", tg))
  plot_heat("orient", tg, paste0("fig_orient_", tg))
}
d <- mets[metric == "runtime_s" & target == "anc"]; order_methods(d)
g <- ggplot(d, aes(label, value)) + geom_boxplot(outlier.size = 0.4) + scale_y_log10() + coord_flip() +
  facet_grid(~ arm) + labs(x = NULL, y = "run time per dataset (s, log scale)")
for (ext in c("pdf", "png")) ggsave(file.path(RES_DIR, paste0("fig_runtime.", ext)), g, width = 12, height = 4, dpi = 200)
cat("Wrote summary.csv, paired.csv, ranks.csv, status.csv, runtime.csv, TABLES.md and figures to", RES_DIR, "\n")
