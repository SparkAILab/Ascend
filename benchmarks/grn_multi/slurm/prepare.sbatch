#!/usr/bin/env bash
## ======================================================================
## benchmarks/grn_multi/slurm/prepare.sbatch
## Full design: simulate every dataset, build the SERGIO arm, and pack the
## task list into chunks (uses the test run's timings in pilot_runs.csv).
## How to run: called by submit_all.sh; by hand: sbatch --wait this file
## ======================================================================
#SBATCH --job-name=grnm_prep
#SBATCH --partition=cpu
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=grnm_prep_%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source benchmarks/grn_multi/slurm/env.sh
Rscript benchmarks/grn_multi/01_make_datasets.R
"$GRNM_PYTHON" benchmarks/grn_multi/py/make_sergio.py --sergio "$GRNM_SERGIO" --ds DS1,DS2 --noise clean,noisy --reps 15
GRNM_PILOT=benchmarks/grn_multi/pilot_runs.csv Rscript benchmarks/grn_multi/02_make_tasks.R
