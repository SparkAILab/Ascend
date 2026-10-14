# ASCEND: Ancestral Scalable Causal discovEry via iNherited Descent

This repository holds the code, data and results behind the paper *Causal ASCEND: scalable two-tier causal discovery on high-dimensional multi-omics data* (Asiedu & Watson).

ASCEND is a constraint-based method for ancestral causal discovery in **two-tier systems**: a set of background variables Z (for example genotypes or methylation sites) is known to causally precede a set of foreground variables X (for example gene expression). ASCEND recovers the ancestral (cause-before-effect) relations among the foreground variables. Each conditional-independence test conditions only on a small, iteratively refined set of nearest ancestors instead of the whole background layer, so the number of tests grows polynomially with the number of variables.

Output: an ancestrality matrix `M` over the foreground, with `M[a, b]` equal to
`1` (a is an ancestor of b), `0.5` (a is not a descendant of b), `0` (no ancestral claim) or `NA` (undetermined).

## Using ASCEND: the R package

If you want to *use* ASCEND on your own data, install the R package in **[`bioc/causalASCEND`](bioc/causalASCEND)** (prepared for Bioconductor). It accepts matrices, `SummarizedExperiment` objects or a `MultiAssayExperiment`, and comes with a vignette:

```r
# from a clone of this repository
install.packages("bioc/causalASCEND", repos = NULL, type = "source")
library(causalASCEND)
sim <- simulateTwoTier(n = 2000, nBackground = 20, nForeground = 15,
    r2 = 0.7, sparsity = 0.9, pCross = 0.2, effect = 0.9, seed = 1)
res <- ascend(sim$background, sim$foreground)
ancestralEdges(res)
vignette("causalASCEND")
```

The rest of this repository is the research code used for the paper's simulations, benchmarks and real-data analyses. It works from plain R scripts (`R/ascend.R`) so every number can be reproduced without installing the package.

## Repository layout

```
R/                        the research library; every script sources from here
  ascend.R                  ASCEND and the two-tier simulator sim_dat()
  sim_tiers.R               simulator with controlled violations of the two-tier assumption
  eval_metrics.R            shared scoring: directed P/R/F1, orientation, SHD, coverage
  stats_utils.R             paired differences with CIs, Wilcoxon, BH within a family
  cbl.R, shah_ss.R          CBL comparator (Watson & Silva)
benchmarks/
  grn/                      ASCEND vs GENIE3 / ARACNe / WGCNA (paper Tables 1-3, Fig. 2)
  grn_multi/                ASCEND vs 15 methods: SEM arms + SERGIO arm, SLURM pipeline
  cbl/                      ASCEND vs CBL sweeps over n, d_x, d_z (paper Fig. 5)
  causal_grid/              81-cell grid vs CBL / GES / LiNGAM / PC on SLURM (paper Fig. 6)
  */results_submitted/      archived benchmark results used by the paper's figures
realdata/
  drosophila/               DGRP analysis and its statistical controls (data not in git)
  yeast/                    ASCEND vs TRIGGER vs BFCS on YEASTRACT ground truth
analysis/                 supplementary analyses; each writes CSVs to analysis/out/
  prescreen_scaling.R       cost and accuracy vs d_z, with and without the background prescreen
  tier_sensitivity.R        robustness to violations of the two-tier assumption
  orientation_analysis.R    which rules drive orientation, and when it fails
  benchmark_tables.R        summary and paired tables (with CIs) from the benchmark CSVs
  make_figures.R            draws docs/figures/ from the CSVs
docs/
  figures/                  figures of the supplementary analyses
  results_local/            full-seed results of the three analyses above
tests/                    base-R tests; run these first
bioc/causalASCEND/        the R package for users
```

Every code file starts with a short header saying what it does, its inputs and outputs, and how to run it. Every script finds the repository root by itself, so you can run it from the root or from any sub-folder.

## Quick start

R >= 4.3. From the repository root:

```bash
Rscript tests/test_metrics.R              # scorer and paired statistics on hand-built graphs
Rscript tests/test_sim_tiers.R            # tier-violation simulator
Rscript tests/test_ascend_regression.R    # R/ascend.R matches a frozen reference implementation
Rscript tests/test_benchmark_smoke.R      # benchmark scoring pipelines (stand-ins if packages are missing)
Rscript tests/test_grn_multi.R            # multi-method benchmark metrics and orientation conventions
Rscript R/ascend.R                        # the 8-gene demo
```

The tests need base R, plus `data.table`, `dplyr` and `R.utils` for the smoke tests. The full benchmarks need:

```r
install.packages(c("data.table", "dplyr", "igraph", "matrixStats", "R.utils", "glmnet",
                   "lightgbm", "doParallel", "doMC", "foreach", "PRROC", "bnlearn",
                   "WGCNA", "ggplot2", "pcalg"))
install.packages("BiocManager")
BiocManager::install(c("graph", "RBGL", "GENIE3", "minet"))
```

`benchmarks/grn_multi/install_R_packages.R` and `benchmarks/grn_multi/requirements.txt` install everything the multi-method benchmark needs (R and Python).

## Reproducing the paper

| Result | How to reproduce |
|---|---|
| GRN benchmark (Tables 1-3, Fig. 2) | `Rscript benchmarks/grn/run_grn_benchmark.R`, then `Rscript benchmarks/grn/plot_f1_heatmap.R` (see the script headers) |
| Multi-method benchmark (16 methods, SEM and SERGIO arms) | [`benchmarks/grn_multi/README.md`](benchmarks/grn_multi/README.md): SLURM pipeline, then `04_summarise.R` and `06_paper_tables.R` |
| ASCEND vs CBL (Fig. 5) | `sbatch benchmarks/cbl/run_ascend_vs_cbl.sbatch` or `Rscript benchmarks/cbl/run_ascend_vs_cbl.R`; plot with `benchmarks/cbl/plot_ascend_vs_cbl.py` |
| Causal-discovery grid (Fig. 6) | `cd benchmarks/causal_grid && bash submit_all.sh`, then `Rscript hpc_merge.R` (see `hpc_run.R`) |
| DGRP analysis and controls | [`realdata/drosophila/README.md`](realdata/drosophila/README.md) |
| Yeast comparison (ASCEND vs TRIGGER vs BFCS) | `BFCS_DIR=/path/to/BFCS Rscript realdata/yeast/yeast_trigger_bfcs_ascend.R` |
| Scaling, prescreen, tier sensitivity, orientation | `Rscript analysis/prescreen_scaling.R`, `analysis/tier_sensitivity.R`, `analysis/orientation_analysis.R` |
| Benchmark summary tables | `Rscript analysis/benchmark_tables.R` (writes `analysis/out/`) |
| Figures in `docs/figures/` | `Rscript analysis/make_figures.R` (reads `analysis/out/`, else `docs/results_local/`) |

Simulations use fixed seeds, so they give the same data on any machine. The scripts that read benchmark results prefer a fresh run (`benchmarks/*/results/`) and fall back to the archived results in `benchmarks/*/results_submitted/`.

## Data availability

- **Simulated data** are generated by the scripts (`R/ascend.R::sim_dat()`, `R/sim_tiers.R`, `benchmarks/grn_multi/01_make_datasets.R`).
- **SERGIO** datasets DS1 and DS2 come from https://github.com/PayamDiba/SERGIO (cloned by `benchmarks/grn_multi/slurm/setup.sh`).
- **DGRP** (Drosophila Genetic Reference Panel): expression from GEO series GSE117850 (`GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz`) and genotypes `dgrp2.tgeno.txt` from https://zenodo.org/records/14871341. These files are too large for git; see `realdata/drosophila/README.md`.
- **Yeast**: expression and genotypes from the `trigger` Bioconductor package (`data(yeast)`), the TRIGGER probability matrix from the BFCS code (`BFCS_DIR`), and the YEASTRACT regulator-target ground truth in `realdata/yeast/yeasttract999.tsv`.

## Citation

Asiedu S, Watson D. Causal ASCEND: scalable two-tier causal discovery on high-dimensional multi-omics data.

Archived code and results: [10.5281/zenodo.23377773](https://doi.org/10.5281/zenodo.23377773)

The CBL comparator in `R/cbl.R` follows https://github.com/dswatson/cbl (Watson & Silva).
