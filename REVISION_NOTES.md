# Revision changes (reviewer round 1)

This branch adds the evaluation, statistics and ablation code needed for the major revision. `ascend()` in its default mode returns exactly the same matrices as before. The only change is extra attributes on the result.

## New files

| File | What it does | Reviewer item |
|---|---|---|
| `eval_metrics.R` | One shared scorer, `eval_ancestral(est, truth)`. It returns directed precision, recall and F1 over all ordered pairs (reversed edges count as errors), orientation accuracy, the undetermined-orientation rate, 3-class ancestor–descendant accuracy, SHD on the ancestral graph, and coverage. `eval_undirected()` scores skeleton-only methods. | 2, 11, minor 9 |
| `stats_utils.R` | `paired_table()`: for each comparison, the mean paired difference with a 95% bootstrap CI, a t-interval, wins/ties/losses, a paired Wilcoxon test, and BH correction within an explicitly named family. Also `summary_table()` (mean ± SE for every metric, precision included). | 9, minor 5 |
| `analysis/revision_tables.R` | Builds the GRN summary for every cell (AUROC/AUPR everywhere), a Table 2/3 replacement with paired CIs, the ASCEND vs CBL coverage table, the conditioning-set ablation cost table and plot, and runtime against discovered edge count. | 2, 5, 6, 9, minor 6, 9 |
| `tests/test_metrics.R` | Unit tests for the two modules above. Run with `Rscript tests/test_metrics.R` from the repo root. | |

## Changes to existing code

- **`ascend.R`**
  - New `cond_set = "full"` ablation. It runs the same algorithm and the same Fisher-z test, but conditions on CBL's full valid set: Z plus the known non-descendants of both endpoints. No Markov-blanket step is run. (Reviewer comment 6.)
  - The result now carries `attr(M, "stats")`: CI-test counts split into pairwise R3, R1/R2 witness, and Markov-blanket tests; time spent inside the CI test; mean and max conditioning-set size; and the number of sweeps. (Reviewer comments 5 and 6.)
  - Markov-blanket results are cached: IAMB is skipped for a variable whose non-descendant pool hasn't changed since its last search. The output is identical and it takes about half the tests.
  - The result also carries `attr(M, "pair_p")`, the p-value of each pair's last R3 test. It serves as a continuous score for AUPR/AUROC.
  - `sim_dat()` has a new `z_scale` argument, which keeps each background variable at unit variance. With the default (`FALSE`), Z variances grow roughly geometrically along the background DAG, and at d_z of about 30 or more the X→X signal vanishes. In the submitted 81-cell grid, ASCEND averaged 0.1 true positives at d_z = 3·d_x. Re-runs use `z_scale = TRUE`. `FALSE` reproduces the submitted data exactly.
- **`GRN comparison/ascend_vs_GRN.R`**
  - `ascend_score_matrix()` and `ascend_binary_skeleton()` used to read only `M[i, j]` with i < j, so ASCEND claims oriented j → i were scored as "no edge". This affected K_match, F1, precision and recall. Both now read either direction. The old versions are kept as `*_legacy`.
  - `direction_acc` is now `orient_acc`. The submitted definition only counted upper-triangle claims and counted false positives as wrong directions; it is kept as `direction_acc_legacy`.
  - Output adds `aupr_cont`/`auroc_cont`, the shared metrics, and `benchmark_v3_paired.csv`/`benchmark_v3_summary.csv`.
- **`causal benchmarks/hpc_run.R`**
  - Scored with `eval_ancestral()`. The old scorer only read the upper triangle in label order, so reversed orientations were never counted as errors.
  - Undirected CPDAG edges are now coded as related-but-unoriented instead of NA.
  - Estimated matrices are saved in each replicate `.rds` (`SAVE_MATRICES=1`), so results can be re-scored without re-running.
  - ASCEND's work statistics are recorded.
  - Uses `Z_SCALE=1` by default.
- **`ascendCBL/ASCEND_vs_CBL.R`**
  - Adds a third arm, `ascend_full`.
  - Records time per test for all methods, including CBL's `l0` fits, conditioning-set sizes, the true maximum in-degree, and the shared metrics.
  - Writes to `results_ascend_vs_cbl_v2.csv`.
  - Environment variables `BENCH_METHODS`, `BENCH_OUT` and `Z_SCALE`.

## Re-running

```
# GRN benchmark (needs GENIE3, minet, WGCNA, PRROC)
cd "GRN comparison" && Rscript -e 'source("../ascend.R"); source("ascend_vs_GRN.R")'
# ASCEND vs ASCEND-full vs CBL
cd ascendCBL && Rscript ASCEND_vs_CBL.R
# 81-cell causal-discovery grid (SLURM, unchanged submission scripts)
cd "causal benchmarks" && bash submit_all.sh
# tables and figures
Rscript analysis/revision_tables.R
```
