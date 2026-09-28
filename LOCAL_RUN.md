# Running the revision code locally

Everything runs from a checkout of the branch `revision/metrics-stats-ablation`. There are four pipelines. Each one has a quick local mode, and the full versions go on the cluster.

## 0. One-time setup

```bash
git fetch origin
git checkout revision/metrics-stats-ablation
```

Install the R packages once (R ≥ 4.3):

```r
install.packages(c("data.table", "dplyr", "tidyverse", "igraph", "matrixStats",
                   "R.utils", "glmnet", "lightgbm", "doParallel", "foreach",
                   "PRROC", "bnlearn", "WGCNA", "doMC"))       # doMC: macOS/Linux only
install.packages("BiocManager")
BiocManager::install(c("graph", "RBGL", "GENIE3", "minet"))
install.packages("pcalg")
```

On a laptop, change `registerDoMC(16)` in `cbl.R` (line 13) to your number of cores.

## 1. Sanity checks (under a minute)

```bash
Rscript tests/test_metrics.R          # prints "eval_metrics tests OK" and "stats_utils tests OK"
Rscript ascend.R                      # the original 8-gene example, unchanged output
```

## 2. Guarded vs full conditioning set (base R only, ~5–10 min)

```bash
Rscript analysis/compare_condsets.R   # SEEDS=10 for more seeds
```

This runs the same algorithm, test and data twice, and only the conditioning set changes. It prints time, CI tests (including how many are Markov-blanket tests), time per test, |S|, and directed P/R/F1 and orientation accuracy for each (n, d_z). Per-run rows go to `analysis/out/compare_condsets.csv`.

## 3. The four pipelines

| Pipeline | Quick local run | Full run | Output |
|---|---|---|---|
| **A. GRN comparison** (Tables 1–3, Fig. 2) | `cd "GRN comparison"`<br>`GRN_QUICK=1 N_REP=3 Rscript -e 'source("../ascend.R"); source("ascend_vs_GRN.R")'`<br>Primary cell only. GENIE3 takes ~2–3 min per rep at n = 2000. | same command without `GRN_QUICK` and `N_REP` (12 cells × 50 reps) | `benchmark_v3_main_raw.csv`, `_wilcoxon.csv`, `_paired.csv`, `_summary.csv` |
| **B. ASCEND vs CBL** (Fig. 5, comment 6 ablation) | `cd ascendCBL`<br>`BENCH_QUICK=1 Rscript ASCEND_vs_CBL.R`<br>6 runs, ~1–2 min | `Rscript ASCEND_vs_CBL.R` (resumable; `BENCH_TIMEOUT` caps each run) | `results_ascend_vs_cbl_v2.csv` |
| **C. Causal-discovery grid** (Fig. 6) | from the **repo root**:<br>`N_REP=2 METHOD_TIMEOUT=300 Rscript "causal benchmarks/hpc_run.R" 1`<br>Task 1 = first combo at n = 512. Tasks 1–81 are all combos at n = 512. | `cd "causal benchmarks" && bash submit_all.sh` (copy `ascend.R`, `cbl.R`, `shah_ss.R`, `eval_metrics.R` next to `hpc_run.R`, as before) | `results/job_XXX/n*_rep*.rds`, then `MERGE_OUT=ascend_benchmark_v3_merged.rds Rscript hpc_merge.R` |
| **D. Tables and figures** | from the repo root:<br>`Rscript analysis/revision_tables.R` | same | `analysis/out/`: GRN summary + paired CIs, Table 2/3 markdown, coverage table, ablation cost plot, runtime vs edges plot |

Pipeline D uses whatever results exist. It prefers the `v3`/`v2` re-run files and falls back to the submitted CSVs, so you can run it right away to see the new tables built from the old data.

## What changed and what stayed the same

**Stays the same**
- The ASCEND algorithm with default arguments. Same output matrix on the same data (checked against the submitted `ascend.R`).
- The CI test (Fisher-z), IAMB, the R1/R2/R3 rules, and closure.
- The parameter grids, seeds, and the SLURM scripts (`submit.sh`, `submit_all.sh`).
- The competitor wrappers (GENIE3, ARACNe, WGCNA, CBL, GES, LiNGAM, PC).

**Changes**
- **Simulator:** background variables are kept at unit variance (`z_scale = TRUE`, your choice). Every re-run gets new data, so all numbers will move. `Z_SCALE=0` reproduces the submitted data.
- **Scoring:** `eval_metrics.R` everywhere. Reversed orientations now count as errors, "direction accuracy" is a real orientation accuracy, and SHD and coverage have one definition. The CSV columns `precision/recall/f1` in pipeline C are now the directed versions.
- **GRN edge counting:** ASCEND claims in either direction are counted. Before, only i < j was counted, which also changes K_match for the competitors. The old functions are kept as `*_legacy`.
- **New outputs:** mean paired differences with CIs; `aupr_cont`/`auroc_cont` (ASCEND ranked by its pairwise test p-value); an `ascend_full` arm and per-test timing in pipeline B; saved estimated matrices and ASCEND work counts in pipeline C.
- **New output file names** (`v3`, `_v2`), so the submitted results are never overwritten.
