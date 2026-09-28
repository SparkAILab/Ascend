## ======================================================================
## benchmarks/grn_multi/config.R
## The whole design of the multi-method benchmark in one place: which
## datasets are simulated, which methods run on them, and the time budget
## of every run.
##
## Reviewer comments: 1 (baselines), 3 (BEELINE-style SERGIO arm),
##          4/minor 8 (non-linear arm), 2 and 9 (every metric, every cell)
## How to run: sourced by 01_make_datasets.R, 02_make_tasks.R, 03_score.R
##   and 04_summarise.R. Edit here, nowhere else. Set GRNM_TEST=1 for the
##   small test design used by slurm/test_pipeline.sh.
## Full write-up: benchmarks/grn_multi/README.md
## ======================================================================

GRNM_TEST <- Sys.getenv("GRNM_TEST", "0") == "1"

## --------------------------------------------------------------------
## 1. Dataset arms
##    Each arm is a grid of simulation settings x replicates. Every
##    dataset gets its own seed, so data are identical on any machine.
## --------------------------------------------------------------------
ARMS <- list(

  # (a) The paper's GRN benchmark design (Tables 1-3, Fig. 2), unchanged
  #     except for the unit-variance background (z_scale = TRUE).
  lin = list(
    generator = "sim_dat",
    grid  = expand.grid(n = c(1000L, 2000L), sp = c(0.5, 0.7, 0.9), r2 = c(0.5, 0.7),
                        KEEP.OUT.ATTRS = FALSE),
    fixed = list(d_z = 20L, d_x = 15L, p_cross = 0.20, x_effect = 0.8, lin_pr = 1),
    reps  = 50L
  ),

  # (b) Non-linear parents (half of each variable's parents enter through
  #     x^2, sqrt|x|, softplus or ReLU). Tests the linear-Gaussian
  #     assumption (minor 8) and is the setting where tree-based GRN
  #     methods should do best.
  nonlin = list(
    generator = "sim_dat",
    grid  = expand.grid(n = 2000L, sp = c(0.7, 0.9), r2 = 0.7, KEEP.OUT.ATTRS = FALSE),
    fixed = list(d_z = 20L, d_x = 15L, p_cross = 0.20, x_effect = 0.8, lin_pr = 0.5),
    reps  = 50L
  ),

  # (c) Larger networks: 50 genes, 100 background variables.
  large = list(
    generator = "sim_dat",
    grid  = expand.grid(n = c(1000L, 2000L), sp = c(0.9, 0.95), r2 = 0.7,
                        KEEP.OUT.ATTRS = FALSE),
    fixed = list(d_z = 100L, d_x = 50L, p_cross = 0.05, x_effect = 0.8, lin_pr = 1),
    reps  = 20L
  ),

  # (d) BEELINE-style arm: SERGIO's published steady-state datasets
  #     (Dibaeinia & Sinha 2020). These are GeneNetWeaver sub-networks of
  #     E. coli (DS1, 100 genes) and yeast (DS2, 400 genes), 9 cell types x
  #     300 cells, 15 replicates each. The master regulators (genes with no
  #     regulator) are the background tier; every other gene is foreground.
  #     "clean" = SERGIO's de-noised expression, log1p.
  #     "noisy" = SERGIO's 10x-like technical noise (outlier genes, library
  #     size, dropout, UMI counts; parameters from SERGIO's tutorial),
  #     library-size normalised, log1p.
  #     Built by py/make_sergio.py, not by 01_make_datasets.R.
  sergio = list(
    generator = "sergio",
    grid  = expand.grid(ds = c("DS1", "DS2"), noise = c("clean", "noisy"),
                        stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE),
    fixed = list(),
    reps  = 15L
  )
)

## Small design for the test run: one replicate of the smallest and the
## largest setting of every arm, and one SERGIO DS1 replicate of each kind.
if (GRNM_TEST) {
  ARMS$lin$grid    <- ARMS$lin$grid[c(1, nrow(ARMS$lin$grid)), ]
  ARMS$nonlin$grid <- ARMS$nonlin$grid[1, , drop = FALSE]
  ARMS$large$grid  <- ARMS$large$grid[1, , drop = FALSE]
  ARMS$sergio$grid <- ARMS$sergio$grid[ARMS$sergio$grid$ds == "DS1", ]
  for (a in names(ARMS)) ARMS[[a]]$reps <- 1L
}

## Primary cell, used for the one-table summary (paper Table 2)
PRIMARY <- list(arm = "lin", n = 2000L, sp = 0.9, r2 = 0.7)

## --------------------------------------------------------------------
## 2. Methods
##    class:
##      causal     outputs a (partially) directed graph; scored on every
##                 directed metric (ASCEND and the tier-aware comparators)
##      directed   GRN method with a regulator -> target score; scored on
##                 the skeleton AND on direction (which way is stronger)
##      undirected symmetric association score; skeleton metrics only
##    lang: which runner executes it (run_task.R or run_task.py)
##    budget_h: wall-clock limit of ONE run (one dataset); a run that hits
##      it is recorded as "timeout" and never retried automatically
##    skip_arms: arms the method is not run on (only where it cannot finish)
## --------------------------------------------------------------------
METHODS <- read.csv(text = "
method,       label,        lang, class,      budget_h, skip_arms,  reference
ascend,       ASCEND,       R,    causal,     6,        ,           this paper
pc_tiered,    Tiered PC,    R,    causal,     6,        ,           Spirtes et al. 2001; tpc: Witte et al. 2022
hc_tiered,    Tiered HC,    R,    causal,     6,        ,           score-based (BIC) hill climbing with tier blacklist; bnlearn (Scutari 2010)
ges,          GES,          R,    causal,     6,        ,           Chickering 2002; pcalg (no tier constraint)
cbl,          CBL,          R,    causal,     24,       large;sergio_DS2, Watson & Silva 2022
notears_tiered, Tiered NOTEARS, py, causal,   12,       ,           Zheng et al. 2018; tier constraint as bounds
genie3,       GENIE3,       R,    directed,   6,        ,           Huynh-Thu et al. 2010
grnboost2,    GRNBoost2,    py,   directed,   6,        ,           Moerman et al. 2019 (network step of SCENIC)
regdiffusion, RegDiffusion, py,   directed,   6,        ,           Zhu & Slonim 2024 (deep learning; successor of DeepSEM/DAZZLE)
pidc,         PIDC,         py,   undirected, 6,        ,           Chan et al. 2017
clr,          CLR,          R,    undirected, 2,        ,           Faith et al. 2007; minet
mrnet,        MRNET,        R,    undirected, 2,        ,           Meyer et al. 2007; minet
aracne,       ARACNe,       R,    undirected, 2,        ,           Margolin et al. 2006; minet
wgcna,        WGCNA,        R,    undirected, 2,        ,           Langfelder & Horvath 2008
ppcor,        PPCOR,        R,    undirected, 2,        ,           Kim 2015 (partial correlation; in BEELINE)
pearson,      Pearson,      R,    undirected, 1,        ,           absolute Pearson correlation (naive baseline)
", strip.white = TRUE, stringsAsFactors = FALSE)
METHODS$skip_arms[is.na(METHODS$skip_arms)] <- ""

## Optional subset: GRNM_METHODS=ascend,genie3,... restricts every script
.sel <- Sys.getenv("GRNM_METHODS", "")
if (nzchar(.sel)) METHODS <- METHODS[METHODS$method %in% strsplit(.sel, ",")[[1]], ]

## The test run caps every run at 2 h: long runs show up as "timeout" and
## 02_make_tasks.R then plans the full budget for them
if (GRNM_TEST) METHODS$budget_h <- pmin(METHODS$budget_h, 2)

## SERGIO DS2 has 400 genes: give every method longer
BUDGET_MULT <- c(sergio_DS2 = 3)
MAX_BUDGET_H <- 40   # must stay below the 48 h partition limit

## Method settings (fixed for every dataset; recorded in each run's meta)
SETTINGS <- list(
  ascend    = list(alpha = 0.05, alpha_mb = 0.05, maxiter = 10L, prescreen = 0.30, min_votes = 1L),
  pc_tiered = list(alpha = 0.05),
  hc_tiered = list(score = "bic-g"),
  ges       = list(score = "GaussL0penObsScore"),
  cbl       = list(gamma = 0.5, B = 50L, maxiter = 10L),
  genie3    = list(nTrees = 1000L, K = "sqrt"),
  minet     = list(estimator = "spearman"),     # minet default
  wgcna     = list(type = "unsigned", RsquaredCut = 0.80)   # WGCNA default; effects have both signs
)

## --------------------------------------------------------------------
## 3. Where things live
## --------------------------------------------------------------------
GRNM_DIR <- local({
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "ascend.R"))) {
    if (dirname(d) == d) stop("Run this script from inside the Ascend repository")
    d <- dirname(d)
  }
  file.path(d, "benchmarks", "grn_multi")
})
ROOT <- dirname(dirname(GRNM_DIR))
# All generated files go under WORK (default benchmarks/grn_multi/work).
# On the cluster, point it at scratch:  export GRNM_WORK=/scratch/users/$USER/grn_multi
WORK <- Sys.getenv("GRNM_WORK", file.path(GRNM_DIR, if (GRNM_TEST) "work_test" else "work"))
DATA_DIR <- file.path(WORK, "data")
RAW_DIR  <- file.path(WORK, "raw")
RES_DIR  <- file.path(WORK, "results")
