# Multi-method benchmark (reviewer comments 1, 2, 3, 9, 11)

This folder holds the benchmark that replaces the submitted GRN comparison, which covered only ASCEND, GENIE3, ARACNe and WGCNA (`benchmarks/grn/`, kept for reference). It compares ASCEND with 15 methods:

- 4 tier-aware or causal comparators
- CBL
- 10 GRN inference methods

The comparison runs on three simulation arms and on SERGIO's published datasets, with every metric computed in every cell. The work is split into small SLURM jobs that resume where they stopped.

**Contents:** [Methods](#methods) · [Direction-aware comparison](#how-direction-is-compared-across-method-types) · [Datasets](#datasets) · [Metrics](#metrics) · [Running on CREATE](#running-it-on-kings-create) · [Outputs](#outputs) · [Files](#files) · [Choosing the method set](#choosing-the-method-set)

## Methods

| Method | Class | Why it is here | Implementation |
|---|---|---|---|
| **ASCEND** | causal | this paper | `R/ascend.R` |
| Tiered PC | causal | constraint-based search with the same tier knowledge (comment 10) | `tpc::tpc` (Witte et al. 2022), Fisher-z, α = 0.05 |
| Tiered HC | causal | score-based (BIC) search with the tier constraint as a blacklist (comment 10) | `bnlearn::hc`, `tiers2blacklist` |
| Tiered NOTEARS | causal | modern continuous-optimisation DAG learner; tiers as bounds | reference code of Zheng et al. 2018, on standardised data |
| GES | causal | classic score-based search, no tier constraint (as in the causal grid) | `pcalg::ges` |
| CBL | causal | the method ASCEND localises | `R/cbl.R`, B = 50, γ = 0.5, lasso |
| GENIE3 | directed GRN | tree ensemble; DREAM4/5 winner | Bioconductor `GENIE3`, 1000 trees |
| GRNBoost2 | directed GRN | the network step of SCENIC (comment 1) | `arboreto`, GRNBoost2 settings |
| RegDiffusion | directed GRN | modern deep learning (comment 1); successor of DeepSEM and DAZZLE by the same group | `regdiffusion` (Zhu & Slonim 2024), CPU |
| PIDC | undirected | information-theoretic; a top performer in BEELINE (comment 1) | port of NetworkInference.jl, Julia not needed |
| CLR | undirected | requested (comment 1) | `minet` |
| MRNET | undirected | requested (comment 1) | `minet` |
| ARACNe | undirected | in the submitted paper | `minet`, ε = 0 |
| WGCNA | undirected | in the submitted paper | `WGCNA`, unsigned TOM |
| PPCOR | undirected | partial correlation, in BEELINE | `ppcor` |
| Pearson | undirected | naive baseline | base R |

Every method receives the same input: all background variables plus all foreground variables. Every run is single-threaded, so run times are comparable.

### Methods left out, and why

These reasons can go in the response letter.

- **SCENIC's motif-pruning step (cisTarget).** It needs motif databases, which simulated genes do not have. GRNBoost2, SCENIC's network-inference step, is included.
- **SINCERITIES, SCODE, LEAP, GRISLI, SCRIBE.** They need pseudotime or time-course data, which neither these simulations nor the DGRP data have.
- **DeepSEM.** It has no maintained package. RegDiffusion comes from the same group and outperformed DeepSEM and DAZZLE on BEELINE, so it stands in as the deep-learning representative.
- **LiNGAM.** It assumes non-Gaussian noise. It is in the causal-discovery grid (`benchmarks/causal_grid`).

## How direction is compared across method types

This is the reviewer's point that ARACNe and WGCNA "do not output oriented edges". Every method is scored on the **skeleton**, i.e. on which pairs are related. Only methods that make directional claims are scored on **direction**, and the scoring depends on the kind of claim.

1. **Causal methods: ASCEND, the three tiered methods, GES and CBL.**
   - They output a directed or partially directed graph.
   - They are scored on every directed metric: directed precision, recall and F1, orientation accuracy, undetermined rate, ancestor–descendant accuracy, SHD and coverage.
   - Tiered PC, Tiered HC and Tiered NOTEARS get exactly ASCEND's background knowledge, so they are the like-for-like comparison.
2. **Directed GRN methods: GENIE3, GRNBoost2, RegDiffusion.**
   - They give a regulator→target score, so they do make a directional claim: they rank a→b against b→a.
   - They are scored in two ways, without adding any threshold of our own:
     - **`orient_acc_rank`**: for every true edge a→b, how often score(a→b) is greater than score(b→a). 0.5 is chance.
     - **`dir_aupr`, `dir_auroc`**: ranking over ordered pairs.
   - They are also scored at ASCEND's edge budget. The top-K pairs are kept, each is oriented towards its larger score, and the result is scored like a causal graph (`dir_f1_Kasc`, `orient_acc_Kasc`).
3. **Undirected methods: PIDC, CLR, MRNET, ARACNe, WGCNA, PPCOR, Pearson.**
   - These are scored on the skeleton only, and their direction metrics are empty.
   - In text, their orientation can be described as chance (0.5).

`04_summarise.R` merges the native and matched-budget versions into two comparison columns, `dirF1` and `orient`, so all direction-capable methods sit in one table.

## Datasets

All datasets are defined in `config.R`.

| Arm | What | Cells × replicates | Genes (d_x) | Background (d_z) | n |
|---|---|---|---|---|---|
| `lin` | the paper's GRN design, linear SEM | 12 × 50 | 15 | 20 | 1000, 2000 |
| `nonlin` | half of all parent effects non-linear (minor 8) | 2 × 50 | 15 | 20 | 2000 |
| `large` | larger networks | 4 × 20 | 50 | 100 | 1000, 2000 |
| `sergio` | BEELINE-style: SERGIO DS1 (100 genes) and DS2 (400 genes), clean and with 10x-like technical noise (comment 3) | 4 × 15 | 93 / 380 | 7 / 20 | 2700 |

All SEM arms use the unit-variance background (`z_scale = TRUE`). Every dataset has a fixed seed, so it is identical on any machine.

**How the SERGIO arm answers comment 3.**
- BEELINE's generators simulate expression only, so no dataset has a genotype layer. SERGIO's datasets are GeneNetWeaver networks (E. coli and yeast sub-networks) simulated with a stochastic model of transcription.
- Master regulators are regulated by nothing, so they satisfy the two-tier order exactly. They are therefore the background tier, and every other gene is foreground.
- These data are non-linear and non-Gaussian, so this arm is the hardest test of ASCEND's linear-Gaussian CI test.
- We follow BEELINE's metrics (AUPR, AUPR ratio, early precision ratio) and state the deviation from its generators with this reason.

## Metrics

Metrics are defined in `R/scoring.R`. Each is computed against two ground truths:

- **`direct`**: the true regulatory edges. This is what GRN methods estimate.
- **`anc`**: the ancestral relations, i.e. the transitive closure of the edges. This is what ASCEND estimates.

Causal methods are scored on their ancestral matrix against `anc`, and on their native graph against `direct`. ASCEND's native graph is the transitive reduction of its ancestor claims.

| Family | Metrics |
|---|---|
| ranking (all methods) | `aupr`, `auroc`, `aupr_ratio`, `epr` (BEELINE), `prec/f1_Ktrue` (top-k, k = number of true pairs), `prec/rec/f1_Kasc` (k = ASCEND's number of claims) |
| direction (causal and directed) | `dir_aupr`, `dir_auroc`, `orient_acc_rank` |
| native graph (causal) | `skel_*`, `dir_precision/recall/f1`, `orient_acc`, `orient_undet`, `ad_acc`, `shd`, `coverage`, `n_claims`; for ASCEND and CBL also `len_*`, which counts ⪯ (0.5) claims as well |
| matched budget (directed) | `dir_f1_Kasc`, `orient_acc_Kasc`, `ad_acc_Kasc`, `shd_Kasc` |
| cost | `runtime_s`, `cpu_s`; ASCEND also records its CI-test counts |

Tied scores are handled exactly:
- AUPR uses the Davis–Goadrich interpolation.
- AUROC counts ties as one half.
- Top-k splits ties at the cut-off in expectation.

For comment 9, `04_summarise.R` reports:
- mean ± SE of every metric;
- the paired difference between ASCEND and every competitor in every cell, with a 95% bootstrap CI, wins/ties/losses and a two-sided Wilcoxon test;
- BH correction within the family *arm × truth × metric*.

## Running it on King's CREATE

Run everything from the repository root.

**1. One-time setup.** First edit the top block of `slurm/env.sh` with your module names (`module avail r`, `module avail python`) and scratch path. Then:
```bash
bash benchmarks/grn_multi/slurm/setup.sh
```
This installs the Python venv (CPU torch) and the R packages, clones SERGIO's datasets (1.7 GB) and runs the tests.

**2. Test run.** It uses 7 datasets (one or two from every arm, including SERGIO DS2) and all 16 methods, with each run capped at 2 h. It takes one job of up to 12 h.
```bash
sbatch benchmarks/grn_multi/slurm/test_pipeline.sbatch
```
When it finishes, read `$GRNM_WORK/../grn_multi_test/results/TEST_REPORT.txt`. It shows:
- **status**: every method should read `ok` on every arm. `timeout` only means the method needs longer than 2 h there. Anything else (`error`, `crashed`) should be fixed before the full run; see `logs/<method>/<dataset>.log`.
- **timings and a projection of the total CPU hours**: the measured times are saved as `benchmarks/grn_multi/pilot_runs.csv` and used to pack the full run into jobs.

The test run's `results/TABLES.md` already shows every table on a handful of datasets.

**3. Full run.**
```bash
bash benchmarks/grn_multi/slurm/submit_all.sh
```
This submits three things in order:
1. **Prepare job:** simulates every dataset (about 900) and packs the tasks into chunks, each expected to take under 20 h on 8 workers.
2. **Array job:** one job per chunk, each with a 48 h limit.
3. **Scoring job:** starts after all array jobs finish.

Each single run also has its own time budget (`budget_h` in `config.R`; ×3 on SERGIO DS2). A run that exceeds it is recorded as `timeout` and does not hold up the rest. Tiered PC did not finish SERGIO DS1 (100 genes) in 6 h in the test, so it is not run on DS2 (400 genes); CBL is not run on `large` or DS2 for the same reason. These are reported as not run, which itself documents the run-time comparison.

**4. If a job is killed or times out.** Run `submit_all.sh` again. Finished runs are never repeated, and only unfinished tasks are re-packed and submitted. To retry the failed runs too, run:
```bash
GRNM_RETRY_ERRORS=1 bash benchmarks/grn_multi/slurm/submit_all.sh
```

**5. Re-scoring.** `sbatch benchmarks/grn_multi/slurm/score.sbatch` works at any time, including mid-run. It never re-runs a method, because each run's raw scores and graphs are kept under `$GRNM_WORK/raw`. That means you can add a metric later without running anything again.

**Running a subset.** For example:
```bash
GRNM_METHODS=ascend,genie3,grnboost2 bash benchmarks/grn_multi/slurm/submit_all.sh
```
To change the design, edit `config.R` (arms, replicates, budgets). Then run `submit_all.sh` again and only the new tasks run.

## Outputs

All outputs are in `$GRNM_WORK/results/`.

| File | Use |
|---|---|
| `TABLES.md` | primary-cell tables (both truths), pooled table per arm, paired F1/AUPR/dirF1 tables per cell, run status |
| `summary.csv` | mean, SE, n_reps for every arm × cell × method × truth × metric (source for any table) |
| `paired.csv` | ASCEND minus each competitor: mean difference, 95% CI, wins/ties/losses, Wilcoxon p, BH q, family |
| `ranks.csv` | ASCEND's rank among all methods, per cell and metric |
| `status.csv`, `runtime.csv` | failures and run times per arm × method (minor 6) |
| `fig_f1K_*`, `fig_aupr_*`, `fig_orient_*`, `fig_runtime` | heatmaps (method × cell, all arms) and run-time boxplots |
| `runs.csv`, `metrics_long.csv.gz` | everything, one row per run / per metric |

## Files

| File | What it does |
|---|---|
| `config.R` | the whole design: arms, methods, budgets, settings |
| `01_make_datasets.R` | simulates the SEM arms |
| `py/make_sergio.py` | builds the SERGIO arm |
| `02_make_tasks.R` | builds the task list and packs it into chunks |
| `run_chunk.sh` | runs one chunk with N workers, each run under `timeout` |
| `run_task.R`, `run_task.py` | run one method on one dataset |
| `R/methods.R`, `py/grn_methods.py` | the method wrappers |
| `R/scoring.R`, `03_score.R` | the metrics, and scoring every run |
| `04_summarise.R` | tables, paired statistics, figures |
| `05_test_report.R` | the test-run report |
| `slurm/` | `env.sh` (edit once), `setup.sh`, `test_pipeline.sbatch`, `submit_all.sh` and the job scripts it uses |

Tests are in `tests/test_grn_multi.R`. They cover the metrics, the graph conversions, and the orientation convention of every comparator.

## Choosing the method set

All 16 methods are run so that nothing needs re-running. However, a reviewer will ask why a method is absent, so the safest course is to fix the set reported in the paper **before** looking at the results, e.g. "all methods in the table above". Any method dropped should then be dropped for a reason stated in the text (cost, assumptions, data type), not because of how it performed. All methods can still be reported in the supplement.
