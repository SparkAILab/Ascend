## ======================================================================
## benchmarks/grn_multi/06_paper_tables.R
## LaTeX tables and figures for the revised manuscript, drawn from the
## output of 04_summarise.R (summary.csv, paired.csv, status.csv).
##
## Reviewer comments: 1 (baselines; main-text panel vs supplement),
##   2 (AUPR/AUROC and directed metrics for every method and cell),
##   3 (SERGIO arm), 9 (mean +/- SE everywhere; paired mean differences
##   with 95% CIs and BH q-values; stated comparison family), minor 5
##   (each significance statement names its metric)
## How to run: Rscript benchmarks/grn_multi/06_paper_tables.R  (seconds)
##   Reads and writes under GRNM_WORK as set in config.R.
## Output (WORK/results/paper/):
##   tab_primary.tex       Table 1: main panel at the primary cell
##   tab_sweep.tex         Table 2: ASCEND vs best main-panel competitor, 12 cells
##   tab_arms.tex          Table 3: main panel pooled over the other arms
##   tabS_primary_all.tex  all 16 methods at the primary cell (supplement)
##   tabS_arms_all.tex     all 16 methods pooled per arm (supplement)
##   tabS_direct.tex       all 16 methods, direct-edge truth (supplement)
##   tabS_aupr_cells.tex   AUPR and AUROC, every method and lin cell (supplement)
##   tabS_paired.tex       paired AUPR / F1@K / directed F1, every competitor (supplement)
##   tabS_paired_arms.tex  paired AUPR in the non-linear, large and SERGIO arms (supplement)
##   tabS_coverage.tex     coverage / directed precision of the causal methods (supplement)
##   tabS_status.tex       run status per method and arm (supplement)
##   fig_main_heat.pdf/png Figure 2: main panel, 12 lin cells
## The main-text panel is MAIN below. The four methods outside it (tiered
## HC, tiered NOTEARS, PPCOR, Pearson) and GES are reported in full in the
## supplement; q-values are BH-corrected over ALL 15 competitors, so
## reporting a subset never loosens the correction.
## ======================================================================

local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) d <- dirname(d)
  setwd(d)
})
source("benchmarks/grn_multi/config.R")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

OUT <- file.path(RES_DIR, "paper"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
summ   <- fread(file.path(RES_DIR, "summary.csv"))
paired <- fread(file.path(RES_DIR, "paired.csv"))
status <- fread(file.path(RES_DIR, "status.csv"))

## Main-text panel: ASCEND, its closest relative (CBL), PC with the same
## tier knowledge, the three methods of the submitted paper, and the
## methods reviewer 1 named (CLR, MRNET, GRNBoost2, PIDC, deep learning).
MAIN <- c("ascend", "cbl", "pc_tiered", "genie3", "grnboost2", "regdiffusion",
          "pidc", "clr", "mrnet", "aracne", "wgcna")
ALL  <- METHODS$method
lab  <- setNames(METHODS$label, METHODS$method)
cls  <- setNames(METHODS$class, METHODS$method)

lin_cells <- unique(summ[arm == "lin", .(n, sp, r2, cell)])[order(r2, sp, n)]
pc <- PRIMARY

f3  <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
f2  <- function(x) ifelse(is.na(x), "--", sprintf("%.2f", x))
fse <- function(m, s, d = 3) ifelse(is.na(m), "--",
  ifelse(is.na(s), sprintf(paste0("%.", d, "f"), m), sprintf(paste0("%.", d, "f $\\pm$ %.", d, "f"), m, s)))
frt <- function(x) ifelse(is.na(x), "--", ifelse(x < 0.1, "$<$0.1", ifelse(x < 10, sprintf("%.1f", x), sprintf("%.0f", x))))
stars <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, "$^{***}$", ifelse(q < 0.01, "$^{**}$", ifelse(q < 0.05, "$^{*}$", ""))))
fdiff <- function(m, lo, hi, q) sprintf("%+.3f [%+.3f, %+.3f]%s", m, lo, hi, stars(q))
bold_best <- function(txt, val, higher = TRUE) {
  ok <- !is.na(val); if (!any(ok)) return(txt)
  b <- if (higher) max(val[ok]) else min(val[ok])
  ifelse(ok & abs(val - b) < 1e-12, paste0("\\textbf{", txt, "}"), txt)
}
write_tex <- function(lines, f) { writeLines(lines, file.path(OUT, f)); cat("wrote", f, "\n") }
tabular <- function(cols, header, rows) c(sprintf("\\begin{tabular}{@{}%s@{}}", cols), "\\toprule",
  paste(header, collapse = " & "), "\\\\\\midrule", paste0(rows, " \\\\"), "\\bottomrule", "\\end{tabular}")

## ---- one-cell table (mean +/- SE) ------------------------------------------
cell_table <- function(d, methods, tg = "anc") {
  M <- c("aupr", "auroc", "f1_Kasc", "prec_Kasc", "dirF1", "orient", "runtime_s")
  d <- d[target == tg & metric %in% M & method %in% methods]
  w <- function(mt, m) { r <- d[metric == mt & method == m]; if (nrow(r)) r else data.table(mean = NA_real_, se = NA_real_) }
  cols <- lapply(M, function(mt) {
    v <- sapply(methods, function(m) w(mt, m)$mean); s <- sapply(methods, function(m) w(mt, m)$se)
    if (mt == "runtime_s") return(frt(v))
    if (mt %in% c("aupr", "f1_Kasc", "dirF1", "auroc")) bold_best(fse(v, s), v) else if (mt == "prec_Kasc") bold_best(f3(v), v) else f3(v)
  })
  rows <- sapply(seq_along(methods), function(i) paste(c(lab[methods[i]], sapply(cols, `[`, i)), collapse = " & "))
  tabular("lccccccr", c("Method", "AUPR", "AUROC", "$F_1$@$K$", "Prec.@$K$", "Dir.\\ $F_1$", "Orient.\\ acc.", "Time (s)"), rows)
}
prim <- summ[arm == pc$arm & n == pc$n & sp == pc$sp & r2 == pc$r2]
write_tex(cell_table(prim, MAIN), "tab_primary.tex")
write_tex(cell_table(prim, ALL), "tabS_primary_all.tex")
write_tex(cell_table(prim, ALL, "direct"), "tabS_direct.tex")

## ---- Table 2: ASCEND vs the best main-panel competitor, every lin cell ------
sweep_rows <- c()
for (i in seq_len(nrow(lin_cells))) {
  cl <- lin_cells[i]; parts <- c()
  for (mt in c("aupr", "f1_Kasc", "dirF1")) {
    s <- summ[arm == "lin" & cell == cl$cell & target == "anc" & metric == mt & method %in% MAIN]
    a <- s[method == "ascend", mean]
    best <- s[method != "ascend"][which.max(mean)]
    p <- paired[arm == "lin" & cell == cl$cell & target == "anc" & metric == mt & competitor == best$method]
    parts <- c(parts, f3(a), sprintf("%s %s", f3(best$mean), paste0("\\scriptsize(", lab[best$method], ")")),
               fdiff(p$mean_diff, p$ci_lo, p$ci_hi, p$q_bh))
  }
  sweep_rows <- c(sweep_rows, paste(c(cl$r2, sprintf("%s / %d", cl$sp, cl$n), parts), collapse = " & "))
}
hdr <- c("$R^2$", "$sp$ / $n$",
         "\\multicolumn{3}{c}{AUPR}", "\\multicolumn{3}{c}{$F_1$ at matched $K$}", "\\multicolumn{3}{c}{Directed $F_1$}")
sub <- c("", "", rep(c("ASCEND", "Best other", "Difference [95\\% CI]"), 3))
write_tex(c("\\begin{tabular}{@{}cc|ccc|ccc|ccc@{}}", "\\toprule", paste(hdr, collapse = " & "), "\\\\",
            paste(sub, collapse = " & "), "\\\\\\midrule", paste0(sweep_rows, " \\\\"), "\\bottomrule", "\\end{tabular}"),
          "tab_sweep.tex")

## ---- pooled per arm; SERGIO per dataset as BEELINE's AUPR ratio ----------------
## SEM arms: mean over cells of AUPR, F1@K and directed F1. SERGIO: AUPR ratio
## (AUPR / edge density; 1 = random) for each of the four datasets, because
## pooling clean and noisy data hides that methods recover signal only from
## the clean ones.
SERGIO_CELLS <- c("DS1 clean", "DS2 clean", "DS1 noisy", "DS2 noisy")
arm_table <- function(methods, arms = c("nonlin", "large")) {
  M <- c("aupr", "f1_Kasc", "dirF1")
  d <- summ[target == "anc" & metric %in% M & method %in% methods & arm %in% arms,
            .(m = mean(mean)), by = .(arm, method, metric)]
  sg <- summ[target == "anc" & metric == "aupr_ratio" & method %in% methods & arm == "sergio",
             .(m = mean(mean)), by = .(cell, method)]
  cols <- list()
  for (a in arms) for (mt in M) {
    v <- sapply(methods, function(x) { r <- d[arm == a & method == x & metric == mt, m]; if (length(r)) r else NA_real_ })
    cols[[length(cols) + 1]] <- bold_best(f3(v), v)
  }
  for (cc in SERGIO_CELLS) {
    v <- sapply(methods, function(x) { r <- sg[cell == cc & method == x, m]; if (length(r) && is.finite(r)) r else NA_real_ })
    cols[[length(cols) + 1]] <- bold_best(f2(v), v)
  }
  rows <- sapply(seq_along(methods), function(i) paste(c(lab[methods[i]], sapply(cols, `[`, i)), collapse = " & "))
  armlab <- c(lin = "Linear, $d_X=15$", nonlin = "Non-linear, $d_X=15$", large = "Large, $d_X=50$")
  c(sprintf("\\begin{tabular}{@{}l%s|cccc@{}}", paste(rep("|ccc", length(arms)), collapse = "")), "\\toprule",
    paste(c("", sprintf("\\multicolumn{3}{c}{%s}", armlab[arms]), "\\multicolumn{4}{c}{SERGIO, AUPR ratio}"), collapse = " & "), "\\\\",
    paste(c("Method", rep(c("AUPR", "$F_1$@$K$", "Dir.\\ $F_1$"), length(arms)),
            "DS1", "DS2", "DS1$^{n}$", "DS2$^{n}$"), collapse = " & "),
    "\\\\\\midrule", paste0(rows, " \\\\"), "\\bottomrule", "\\end{tabular}")
}
write_tex(arm_table(MAIN), "tab_arms.tex")
write_tex(arm_table(ALL, c("lin", "nonlin", "large")), "tabS_arms_all.tex")

## ---- AUPR and AUROC in every lin cell, every method (comment 2) --------------
rows <- c()
for (m in ALL) {
  for (mt in c("aupr", "auroc")) {
    v <- sapply(lin_cells$cell, function(cc) { r <- summ[arm == "lin" & cell == cc & method == m & target == "anc" & metric == mt, mean]; if (length(r)) r else NA })
    rows <- c(rows, paste(c(if (mt == "aupr") lab[m] else "", toupper(sub("aupr", "AUPR", sub("auroc", "AUROC", mt))), f2(v)), collapse = " & "))
  }
}
ch <- sprintf("%s/%s/%d", lin_cells$r2, lin_cells$sp, lin_cells$n)
write_tex(tabular(paste0("ll", strrep("c", nrow(lin_cells))),
                  c("Method", "", sprintf("\\rotatebox{90}{%s}", ch)), rows), "tabS_aupr_cells.tex")

## ---- paired differences, every competitor, lin arm --------------------------
prow <- c()
for (m in setdiff(ALL, "ascend")) {
  for (mt in c("aupr", "f1_Kasc", "dirF1")) {
    p <- paired[arm == "lin" & target == "anc" & metric == mt & competitor == m]
    if (!nrow(p)) next
    prow <- c(prow, paste(c(if (mt == "aupr") lab[m] else "", c(aupr = "AUPR", f1_Kasc = "$F_1$@$K$", dirF1 = "Dir.\\ $F_1$")[mt],
      sprintf("%+.3f", mean(p$mean_diff)), sprintf("%+.3f to %+.3f", min(p$mean_diff), max(p$mean_diff)),
      sum(p$mean_diff > 0 & p$q_bh < 0.05), sum(p$mean_diff < 0 & p$q_bh < 0.05),
      sprintf("%.0f\\%%", 100 * sum(p$wins) / sum(p$n_pairs))), collapse = " & "))
  }
}
write_tex(tabular("llcccccc", c("Competitor", "Metric", "Mean diff.", "Range over cells", "ASCEND sig.\\ better", "ASCEND sig.\\ worse", "ASCEND wins"), prow),
          "tabS_paired.tex")

## ---- paired AUPR in the other arms (count of cells ASCEND is sig. better/worse) ---
pa <- paired[target == "anc" & metric == "aupr" & arm != "lin"]
pa[, grp := fifelse(arm == "sergio", cell, arm)]
grps <- c("nonlin", "large", SERGIO_CELLS)
arow <- c()
for (m in setdiff(ALL, "ascend")) {
  v <- sapply(grps, function(g) {
    p <- pa[grp == g & competitor == m]
    if (!nrow(p)) return("--")
    sprintf("%+.3f (%d/%d/%d)", mean(p$mean_diff), sum(p$mean_diff > 0 & p$q_bh < 0.05),
            sum(p$q_bh >= 0.05 | is.na(p$q_bh)), sum(p$mean_diff < 0 & p$q_bh < 0.05))
  })
  arow <- c(arow, paste(c(lab[m], v), collapse = " & "))
}
write_tex(tabular(paste0("l", strrep("c", length(grps))),
                  c("Competitor", "Non-linear", "Large", "DS1", "DS2", "DS1 noisy", "DS2 noisy"), arow),
          "tabS_paired_arms.tex")

## ---- coverage of the causal methods, every lin cell (minor 9) -----------------
CAUSAL <- METHODS$method[METHODS$class == "causal"]
crow <- c()
for (i in seq_len(nrow(lin_cells))) {
  cl <- lin_cells[i]
  v <- sapply(CAUSAL, function(m) {
    cv <- summ[arm == "lin" & cell == cl$cell & method == m & target == "anc" & metric == "coverage", mean]
    pr <- summ[arm == "lin" & cell == cl$cell & method == m & target == "anc" & metric == "dir_precision", mean]
    if (!length(cv)) "--" else sprintf("%.2f / %.2f", cv, if (length(pr)) pr else NA)
  })
  crow <- c(crow, paste(c(cl$r2, sprintf("%s / %d", cl$sp, cl$n), v), collapse = " & "))
}
write_tex(tabular(paste0("cc", strrep("c", length(CAUSAL))), c("$R^2$", "$sp$ / $n$", lab[CAUSAL]), crow),
          "tabS_coverage.tex")

## ---- run status ----------------------------------------------------------------
st <- copy(status); for (k in c("ok", "error", "timeout")) if (!k %in% names(st)) st[, (k) := 0L]
st[, txt := ifelse(error + timeout == 0, as.character(ok), sprintf("%d / %d / %d", ok, error, timeout))]
sw <- dcast(st, method ~ arm, value.var = "txt"); sw <- sw[match(ALL, sw$method)]
srow <- apply(sw, 1, function(r) paste(c(lab[r[["method"]]], ifelse(is.na(r[-1]), "skipped", r[-1])), collapse = " & "))
write_tex(tabular(paste0("l", strrep("c", ncol(sw) - 1)), c("Method", names(sw)[-1]), srow), "tabS_status.tex")

## ---- Figure 2: main panel, 12 lin cells -----------------------------------------
d <- summ[arm == "lin" & target == "anc" & metric %in% c("aupr", "f1_Kasc", "dirF1") & method %in% MAIN]
d[, label := factor(lab[method], levels = rev(lab[MAIN]))]
d[, cellf := factor(sprintf("R²=%s\nsp=%s\nn=%d", r2, sp, n),
                    levels = sprintf("R²=%s\nsp=%s\nn=%d", lin_cells$r2, lin_cells$sp, lin_cells$n))]
d[, panel := factor(c(aupr = "(a) AUPR", f1_Kasc = "(b) F1 at matched edge count", dirF1 = "(c) Directed F1")[metric],
                    levels = c("(a) AUPR", "(b) F1 at matched edge count", "(c) Directed F1"))]
d[, best := mean == max(mean), by = .(cell, metric)]
g <- ggplot(d, aes(cellf, label, fill = mean)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", mean), fontface = ifelse(best, "bold", "plain")), size = 2.1) +
  scale_fill_gradient(low = "#f7fbff", high = "#2171b5", limits = c(0, 1), name = NULL) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        legend.position = "right", legend.key.height = unit(1.2, "cm"))
for (ext in c("pdf", "png")) ggsave(file.path(OUT, paste0("fig_main_heat.", ext)), g, width = 7.2, height = 5.4, dpi = 300)
cat("Wrote paper tables and Figure 2 to", OUT, "\n")
