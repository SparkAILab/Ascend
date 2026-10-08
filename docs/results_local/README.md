# Full-seed local runs (8 October 2026)

These replace the pilot values in `docs/pilot_results/` for three analyses, all run on one CPU core each:

| File | Script | Design |
|---|---|---|
| `tier_sensitivity*.csv` | `SEEDS=30 Rscript analysis/tier_sensitivity.R` | 4 violation types x 5 levels x 30 seeds |
| `orientation_*.csv` | `SEEDS=20 Rscript analysis/orientation_analysis.R` | 6 cells x 20 seeds (12,600 pairs); the per-pair file is not kept here (regenerate it) |
| `prescreen_scaling.csv` | `SEEDS=20 Rscript analysis/prescreen_scaling.R` | n in {300, 1024} x d_z in {25, 50, 100, 200} x 20 seeds, prescreen on/off |

`Rscript analysis/make_report_figures.R` redraws docs/figures from analysis/out/.
