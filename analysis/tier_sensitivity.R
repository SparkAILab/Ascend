## ======================================================================
## Reviewer comment 4: sensitivity to violations of the two-tier assumption
## ----------------------------------------------------------------------
## One violation type at a time, at increasing rates, everything else at
## the default setting. See sim_tiers.R for the definition of each type.
##
## Run from the repo root (base R + R.utils + ggplot2):
##   Rscript analysis/tier_sensitivity.R               # 10 seeds, ~15-30 min
##   SEEDS=3 Rscript analysis/tier_sensitivity.R       # quick look
##   N=2000 DZ=100 DX=20 Rscript analysis/tier_sensitivity.R
## Writes analysis/out/tier_sensitivity.csv (one row per run),
##        analysis/out/tier_sensitivity_summary.csv and .pdf
## ======================================================================

suppressPackageStartupMessages({ library(R.utils); library(ggplot2) })
source("ascend.R"); source("eval_metrics.R"); source("sim_tiers.R"); source("stats_utils.R")

env_num <- function(k, d) as.numeric(Sys.getenv(k, d))
SEEDS   <- 100L + seq_len(env_num("SEEDS", "10"))
N  <- env_num("N", "1024"); DZ <- env_num("DZ", "50"); DX <- env_num("DX", "15")
SP <- env_num("SP", "0.7"); R2 <- env_num("R2", "0.5")
TIMEOUT <- env_num("TIMEOUT", "600")

LEVELS <- list(
  x_as_z   = c(0, 0.05, 0.10, 0.20, 0.30),   # foreground mislabelled as background (violates A1)
  feedback = c(0, 0.05, 0.10, 0.20, 0.30),   # X -> Z reverse causation         (violates A1)
  hide_z   = c(0, 0.25, 0.50, 0.75, 0.90),   # unmeasured background (latent confounding)
  z_as_x   = c(0, 0.05, 0.10, 0.20, 0.30)    # background mislabelled as foreground
)

run_one <- function(type, level, seed) {
  args <- list(n = N, d_z = DZ, d_x = DX, r2 = R2, sp = SP, seed = seed)
  args[[switch(type, x_as_z = "p_x_as_z", feedback = "p_feedback",
               hide_z = "p_hide_z", z_as_x = "p_z_as_x")]] <- level
  sim <- do.call(sim_tiers, args)
  t0 <- Sys.time()
  M <- tryCatch(withTimeout(ascend(sim, verbose = FALSE), timeout = TIMEOUT, onTimeout = "error"),
                error = function(e) NULL)
  el <- as.numeric(Sys.time() - t0, units = "secs")
  base <- data.frame(type = type, level = level, seed = seed,
                     d_z_obs = sim$params$d_z, d_x_obs = sim$params$d_x, time_s = el)
  if (is.null(M)) return(cbind(base, status = "timeout/error"))
  e <- eval_ancestral(M, sim$truth)
  # errors that touch a violating variable (mislabelled Z -> X, feedback path)
  cbind(base, status = "ok",
        as.data.frame(e[c("dir_precision", "dir_recall", "dir_f1", "orient_acc",
                          "orient_undet", "ad_acc", "coverage", "shd", "n_true_rel",
                          "dir_claims", "n_reversed")]))
}

rows <- list()
for (type in names(LEVELS)) for (lv in LEVELS[[type]]) for (s in SEEDS) {
  rows[[length(rows) + 1]] <- run_one(type, lv, s)
  cat(sprintf("%-9s level=%.2f seed=%d done\n", type, lv, s))
}
res <- do.call(dplyr::bind_rows, rows)
dir.create("analysis/out", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "analysis/out/tier_sensitivity.csv", row.names = FALSE)

mets <- c("dir_precision", "dir_recall", "dir_f1", "orient_acc", "coverage", "shd")
ok   <- res[res$status == "ok", ]; ok$method <- "ASCEND"
summ <- summary_table(ok, mets, c("type", "level"))
write.csv(summ, "analysis/out/tier_sensitivity_summary.csv", row.names = FALSE)

long <- do.call(rbind, lapply(setdiff(mets, "shd"), function(m)
  data.frame(type = summ$type, level = summ$level, metric = m,
             mean = summ[[paste0(m, "_mean")]], se = summ[[paste0(m, "_se")]])))
p <- ggplot(long, aes(level, mean, colour = metric)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.01) +
  facet_wrap(~type, scales = "free_x",
             labeller = as_labeller(c(x_as_z = "Foreground labelled as background",
                                      feedback = "Feedback X -> Z",
                                      hide_z = "Unmeasured background",
                                      z_as_x = "Background labelled as foreground"))) +
  labs(x = "Violation rate", y = "Mean (+/- SE)", colour = NULL) + theme_bw()
ggsave("analysis/out/tier_sensitivity.pdf", p, width = 9, height = 6)

options(width = 200)
print(summ[, c("type", "level", paste0(mets, "_mean"))], digits = 3, row.names = FALSE)
