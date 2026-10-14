# Local results

Local runs of the analyses below, each on one CPU core. `analysis/make_figures.R` reads these files when `analysis/out/` does not hold your own run.

| File | Script | Design |
|---|---|---|
| `tier_sensitivity*.csv` | `SEEDS=30 Rscript analysis/tier_sensitivity.R` | 4 violation types x 5 levels x 30 seeds |
| `orientation_*.csv` | `SEEDS=20 Rscript analysis/orientation_analysis.R` | 6 cells x 20 seeds (12,600 pairs); the per-pair file is not kept here (regenerate it) |
| `results_ascend_vs_cbl_v2.csv` | `BENCH_SEEDS=2 Rscript benchmarks/cbl/run_ascend_vs_cbl.R` (subset of cells) | ASCEND vs CBL, 9 cells (n, d_x and d_z sweeps) x 2 seeds; the full sweep is `benchmarks/cbl/run_ascend_vs_cbl.sbatch` |
| `prescreen_scaling.csv` | `SEEDS=20 Rscript analysis/prescreen_scaling.R` | n in {300, 1024} x d_z in {25, 50, 100, 200} x 20 seeds, prescreen on/off |

`Rscript analysis/make_figures.R` redraws `docs/figures/` from `analysis/out/` (falling back to this folder).
