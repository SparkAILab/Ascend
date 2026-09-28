# ASCEND: Ancestral Scalable Causal discovEry via iNherited Descent

This repository holds the code for *Causal ASCEND* (Asiedu & Watson). The paper covers constraint-based ancestral discovery for two-tier systems, where background variables Z (for example genotypes) are known to precede foreground variables X (for example expression).

Branch `revision/metrics-stats-ablation` holds the code for the major revision. **[docs/REVISION_REPORT.md](docs/REVISION_REPORT.md)** covers:

- what was changed and why
- what each analysis found
- how to reproduce every number and figure

## Layout

```
R/                      the library: every script sources from here
  ascend.R                ASCEND + two-tier simulator sim_dat()
  sim_tiers.R             simulator with tier violations (comment 4)
  eval_metrics.R          shared scoring: directed P/R/F1, orientation, SHD, coverage
  stats_utils.R           paired differences with CIs, Wilcoxon, BH families
  cbl.R, shah_ss.R        CBL comparator (Watson & Silva)
benchmarks/
  grn/                    ASCEND vs GENIE3 / ARACNe / WGCNA   (Tables 1-3, Fig. 2)
  cbl/                    ASCEND vs CBL sweeps                (Fig. 5)
  causal_grid/            81-cell grid vs CBL / GES / LiNGAM / PC on SLURM (Fig. 6)
  */results_submitted/    the results behind the submitted manuscript (never overwritten)
realdata/
  drosophila/             DGRP analysis (data files are not in git)
  yeast/                  ASCEND vs TRIGGER vs BFCS on YEASTRACT
analysis/               revision analyses; each writes CSVs to analysis/out/
  scoring_audit.R         what the submitted scoring bugs changed
  simulator_diagnostic.R  why the simulator was fixed
  prescreen_scaling.R     cost and accuracy vs d_z, prescreen setting
  tier_sensitivity.R      comment 4
  orientation_analysis.R  comment 11
  revision_tables.R       tables with CIs for the benchmarks
  make_report_figures.R   draws docs/figures/ from the CSVs
tests/                  run these first (base R)
docs/
  REVISION_REPORT.md      the findings report
  figures/                every figure in the report
  pilot_results/          the pilot CSVs the report quotes
```

Every code file starts with a short brief covering what it does, which reviewer comments it addresses, and how to run it. Every script finds the repository root by itself, so you can run it from the root or from any sub-folder.

## Quick start

```bash
git checkout revision/metrics-stats-ablation
Rscript tests/test_metrics.R              # scorer and statistics
Rscript tests/test_sim_tiers.R            # tier-violation simulator
Rscript tests/test_ascend_regression.R    # revised ascend() == submitted ascend()
Rscript tests/test_benchmark_smoke.R      # benchmark scoring pipelines (stubs if packages missing)
Rscript R/ascend.R                        # the 8-gene demo
```

The report's section 9 lists the R packages needed for the full benchmarks.
