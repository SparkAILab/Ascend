## ======================================================================
## realdata/drosophila/dgrp_controls_summarise.R
## Collects the DGRP control runs and writes every number the paper's DGRP
## paragraph and Supplementary Section S10 need, with two figures.
##
## Reviewer comments: 7 (negative controls, stability, selection bias,
##   number of tests and false discovery proportion, formal enrichment of
##   immune genes among the hubs), minor 3 (gene symbols from FlyBase IDs),
##   minor 7 (filtering counts)
## How to run (after dgrp_controls.R has done its runs):
##   Rscript realdata/drosophila/dgrp_controls_summarise.R
##   DGRP_OUT     = folder of dgrp_controls.R (default realdata/drosophila/results)
##   DGRP_GENESET = optional file of FlyBase IDs (one per line) to test for
##                  enrichment; without it the script uses GO through
##                  org.Dm.eg.db (BiocManager::install("org.Dm.eg.db"))
## Output (in DGRP_OUT): DGRP_CONTROLS.md (read this first), dgrp_flow.csv,
##   dgrp_hubs.csv, dgrp_null.csv, dgrp_stability.csv, dgrp_enrichment.csv,
##   fig_dgrp_null.pdf/png, fig_dgrp_stability.pdf/png
## Full write-up: docs/REVISION_REPORT.md
## ======================================================================

suppressMessages({ library(data.table); library(ggplot2) })
ROOT <- local({                      # repo root = folder holding R/ascend.R
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  d
})
OUT <- Sys.getenv("DGRP_OUT", file.path(ROOT, "realdata", "drosophila", "results"))
TOP_HUBS <- 20   # hubs tested for enrichment
files <- list.files(file.path(OUT, "runs"), pattern = "\\.rds$", full.names = TRUE)
runs  <- lapply(files, readRDS)
type  <- vapply(runs, `[[`, "", "type")
if (!"observed" %in% type) stop("The observed run (run 1) is missing in ", OUT)
obs <- runs[[which(type == "observed")]]
md <- c("# DGRP statistical controls", "",
        sprintf("Runs found: %s.", paste(sprintf("%s %d", names(table(type)), table(type)), collapse = ", ")), "")

## ---- gene symbols (minor 3): org.Dm.eg.db if installed --------------------
have_org <- requireNamespace("org.Dm.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)
sym <- setNames(rep(NA_character_, length(obs$genes)), obs$genes)
if (have_org) {
  m <- suppressMessages(AnnotationDbi::select(org.Dm.eg.db::org.Dm.eg.db, keys = obs$genes,
                                              keytype = "FLYBASE", columns = "SYMBOL"))
  m <- m[!duplicated(m$FLYBASE), ]
  sym[m$FLYBASE] <- m$SYMBOL
}

## ---- filtering flow (minor 7) ---------------------------------------------
flow <- data.table(step = names(obs$flow), count = as.numeric(obs$flow))
fwrite(flow, file.path(OUT, "dgrp_flow.csv"))
md <- c(md, "## Filtering flow (Supplementary Fig. S8)", "",
        "| step | count |", "|---|---|", sprintf("| %s | %s |", flow$step, format(flow$count, big.mark = ",")), "")

## ---- observed network and hubs ---------------------------------------------
hubs <- data.table(generic_id = paste0("x", seq_along(obs$genes)), gene_id = obs$genes,
                   symbol = sym[obs$genes], out_degree = obs$out_degree, in_degree = obs$in_degree)
setorder(hubs, -out_degree)
hubs[, rank := seq_len(.N)]
fwrite(hubs, file.path(OUT, "dgrp_hubs.csv"))
md <- c(md, "## Observed network", "",
        sprintf("- %d directed edges, %d undirected entries, density %.2f%% (paper: 880, 61, 3.02%%)",
                obs$n_directed, obs$n_undirected, 100 * obs$density),
        sprintf("- genes with out-degree > 0: %d of %d; with out-degree >= 5: %d",
                sum(hubs$out_degree > 0), nrow(hubs), sum(hubs$out_degree >= 5)),
        sprintf("- top hubs: %s", paste(sprintf("%s (%s, %d)", hubs$gene_id[1:7],
                                                 ifelse(is.na(hubs$symbol[1:7]), "?", hubs$symbol[1:7]),
                                                 hubs$out_degree[1:7]), collapse = "; ")),
        sprintf("- ASCEND run time %.0f s", obs$seconds), "")

## ---- multiplicity -----------------------------------------------------------
st <- obs$stats
md <- c(md, "## Number of tests", "",
        sprintf("- %s CI tests in total (%s pairwise, %s witness, %s Markov blanket), %d sweeps",
                format(st$n_ci_tests, big.mark = ","), format(st$n_pair_tests, big.mark = ","),
                format(st$n_witness_tests, big.mark = ","), format(st$n_mb_tests, big.mark = ","), st$n_sweeps),
        "- Pairwise dependence calls are BH-controlled per sweep at q = 0.05, so the expected false discovery proportion among them is at most 5%.", "")

## ---- negative controls and random gene sets ---------------------------------
one <- function(r) data.table(type = r$type, rep = r$rep, n_directed = r$n_directed,
                              n_undirected = r$n_undirected, max_out = max(r$out_degree),
                              n_hubs5 = sum(r$out_degree >= 5), seconds = r$seconds)
all <- rbindlist(lapply(runs, one))
fwrite(all, file.path(OUT, "dgrp_null.csv"))
emp_p <- function(null, o) (1 + sum(null >= o)) / (1 + length(null))
null_rows <- list()
for (tp in intersect(c("perm_geno", "swap_tiers", "random_genes"), all$type)) {
  d <- all[type == tp]
  for (stat in c("max_out", "n_directed", "n_hubs5")) {
    o <- all[type == "observed"][[stat]]
    null_rows[[length(null_rows) + 1]] <- data.table(
      control = tp, statistic = stat, observed = o, runs = nrow(d),
      null_mean = mean(d[[stat]]), null_q025 = quantile(d[[stat]], 0.025),
      null_q975 = quantile(d[[stat]], 0.975), null_max = max(d[[stat]]),
      empirical_p = emp_p(d[[stat]], o))
  }
}
nt <- rbindlist(null_rows)
fwrite(nt, file.path(OUT, "dgrp_null_summary.csv"))
lab <- c(perm_geno = "Genotypes permuted across lines", swap_tiers = "Tier labels randomised",
         random_genes = "250 random genes passing the filter")
md <- c(md, "## Negative controls and selection bias", "",
        "Empirical p = (1 + number of control runs at or above the observed value) / (1 + runs).", "",
        "| control | statistic | observed | control mean [2.5%, 97.5%] | control max | runs | empirical p |",
        "|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %g | %.1f [%.1f, %.1f] | %g | %d | %.3f |", lab[nt$control], nt$statistic,
                nt$observed, nt$null_mean, nt$null_q025, nt$null_q975, nt$null_max, nt$runs, nt$empirical_p), "",
        "max_out = largest out-degree; n_hubs5 = genes with out-degree >= 5.", "")

if (nrow(nt)) {
  pd <- all[type %in% names(lab)]
  pd[, control := factor(lab[type], lab)]
  ov <- all[type == "observed"]
  pl <- melt(pd, id.vars = c("control"), measure.vars = c("max_out", "n_directed"))
  ol <- melt(ov, measure.vars = c("max_out", "n_directed"))
  g <- ggplot(pl, aes(value)) + geom_histogram(bins = 20, fill = "grey70") +
    geom_vline(data = ol, aes(xintercept = value), colour = "#1b9e77", linewidth = 0.8) +
    facet_grid(control ~ variable, scales = "free",
               labeller = labeller(variable = c(max_out = "Largest out-degree", n_directed = "Directed edges"))) +
    labs(x = NULL, y = "Control runs", caption = "Green line: observed DGRP network") +
    theme_bw(base_size = 9)
  for (ext in c("pdf", "png")) ggsave(file.path(OUT, paste0("fig_dgrp_null.", ext)), g, width = 6.5, height = 5.5, dpi = 200)
}

## ---- each top hub against its own out-degree with permuted genotypes ----------
## The largest out-degree can be matched by chance hubs elsewhere in the
## network; the specific claim is that THESE genes are hubs, so each top hub
## is also compared with its own out-degree across the permuted-genotype runs.
perm <- runs[type == "perm_geno"]
if (length(perm)) {
  top <- hubs[1:10]
  po <- sapply(perm, function(r) r$out_degree[match(top$gene_id, r$genes)])
  ph <- data.table(top[, .(gene_id, symbol, out_degree)],
                   perm_mean = rowMeans(po), perm_q975 = apply(po, 1, quantile, 0.975),
                   perm_max = apply(po, 1, max),
                   empirical_p = (1 + rowSums(po >= top$out_degree)) / (1 + length(perm)))
  fwrite(ph, file.path(OUT, "dgrp_hub_perm.csv"))
  md <- c(md, sprintf("## Each top hub against its own out-degree in %d permuted-genotype runs", length(perm)), "",
          "| hub | symbol | observed out-degree | permuted mean | permuted 97.5% | permuted max | empirical p |",
          "|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %d | %.1f | %.1f | %d | %.3f |", ph$gene_id, ifelse(is.na(ph$symbol), "", ph$symbol),
                  ph$out_degree, ph$perm_mean, ph$perm_q975, ph$perm_max, ph$empirical_p), "")
}

## ---- stability under subsampling ----------------------------------------------
sub <- runs[type == "subsample"]
if (length(sub)) {
  top <- hubs[1:10]
  rk <- sapply(sub, function(r) {
    o <- r$out_degree; rank(-o, ties.method = "min")[match(top$gene_id, r$genes)]
  })
  stab <- data.table(top[, .(gene_id, symbol, out_degree, rank)],
                     in_top5 = rowMeans(rk <= 5), in_top10 = rowMeans(rk <= 10),
                     in_top20 = rowMeans(rk <= 20), median_rank = apply(rk, 1, median))
  fwrite(stab, file.path(OUT, "dgrp_stability.csv"))
  ef <- table(unlist(lapply(sub, `[[`, "edges")))
  freq <- as.numeric(ef[obs$edges]); freq[is.na(freq)] <- 0; freq <- freq / length(sub)
  md <- c(md, sprintf("## Stability over %d subsamples of 80%% of lines", length(sub)), "",
          "| hub | symbol | observed out-degree (rank) | in top 5 | in top 10 | in top 20 | median rank |",
          "|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %d (%d) | %.0f%% | %.0f%% | %.0f%% | %g |", stab$gene_id,
                  ifelse(is.na(stab$symbol), "", stab$symbol), stab$out_degree, stab$rank,
                  100 * stab$in_top5, 100 * stab$in_top10, 100 * stab$in_top20, stab$median_rank), "",
          sprintf("- Directed edges of the observed network: median selection frequency %.2f; %.0f%% are found in at least half of the subsamples.",
                  median(freq), 100 * mean(freq >= 0.5)), "")
  sd <- melt(stab[, .(hub = factor(ifelse(is.na(symbol), gene_id, symbol), rev(ifelse(is.na(symbol), gene_id, symbol))),
                      in_top5, in_top10, in_top20)], id.vars = "hub")
  g <- ggplot(sd, aes(value, hub, fill = variable)) + geom_col(position = "dodge") +
    scale_fill_manual(values = c(in_top5 = "#1b9e77", in_top10 = "#7570b3", in_top20 = "#bbbbbb"),
                      labels = c("top 5", "top 10", "top 20")) +
    labs(x = "Share of 80% subsamples", y = NULL, fill = "Hub ranked in") + theme_bw(base_size = 9)
  for (ext in c("pdf", "png")) ggsave(file.path(OUT, paste0("fig_dgrp_stability.", ext)), g, width = 5.5, height = 3.8, dpi = 200)
}

## ---- enrichment of immune genes among the top hubs ------------------------------
sets <- list()
gs_file <- Sys.getenv("DGRP_GENESET", "")
if (nzchar(gs_file)) {
  sets[[basename(gs_file)]] <- unique(trimws(readLines(gs_file)))
} else if (have_org) {
  go <- c(`defense response (GO:0006952)` = "GO:0006952",
          `Toll signalling (GO:0008063)` = "GO:0008063",
          `Imd / PGRP signalling (GO:0061057)` = "GO:0061057",
          `antimicrobial humoral response (GO:0019730)` = "GO:0019730")
  for (nm in names(go)) {
    m <- tryCatch(suppressMessages(AnnotationDbi::select(org.Dm.eg.db::org.Dm.eg.db, keys = go[[nm]],
                                                          keytype = "GOALL", columns = "FLYBASE")),
                  error = function(e) NULL)
    if (!is.null(m)) sets[[nm]] <- unique(m$FLYBASE)
  }
  tollimd <- unique(unlist(sets[grep("Toll|Imd", names(sets))]))
  sets[["Toll or Imd (union of the two above)"]] <- tollimd
}
if (length(sets)) {
  universe <- obs$genes; top_ids <- hubs$gene_id[1:TOP_HUBS]
  en <- rbindlist(lapply(names(sets), function(nm) {
    K <- sum(universe %in% sets[[nm]]); k <- sum(top_ids %in% sets[[nm]])
    data.table(gene_set = nm, in_top = k, top = TOP_HUBS, in_universe = K, universe = length(universe),
               expected = TOP_HUBS * K / length(universe),
               p_hypergeometric = phyper(k - 1, K, length(universe) - K, TOP_HUBS, lower.tail = FALSE),
               genes = paste(ifelse(is.na(sym[top_ids[top_ids %in% sets[[nm]]]]), top_ids[top_ids %in% sets[[nm]]],
                                    sym[top_ids[top_ids %in% sets[[nm]]]]), collapse = ", "))
  }))
  fwrite(en, file.path(OUT, "dgrp_enrichment.csv"))
  md <- c(md, sprintf("## Enrichment among the top %d hubs (universe: the %d analysed genes)", TOP_HUBS, length(universe)), "",
          "| gene set | in top hubs | expected | in the 250 genes | hypergeometric p | genes |", "|---|---|---|---|---|---|",
          sprintf("| %s | %d / %d | %.1f | %d | %.2g | %s |", en$gene_set, en$in_top, en$top, en$expected,
                  en$in_universe, en$p_hypergeometric, en$genes), "")
} else {
  md <- c(md, "## Enrichment", "", "Not run: install org.Dm.eg.db or set DGRP_GENESET to a file of FlyBase IDs.", "")
}

writeLines(md, file.path(OUT, "DGRP_CONTROLS.md"))
cat(md, sep = "\n")
