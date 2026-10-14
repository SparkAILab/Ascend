#!/usr/bin/env bash
## ======================================================================
## benchmarks/grn_multi/slurm/score.sbatch
## Scores every finished run and writes the summary tables and figures.
## Runs after all method jobs (submit_all.sh sets the dependency); can be
## re-run at any time, e.g. while jobs are still going.
## ======================================================================
#SBATCH --job-name=grnm_score
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=grnm_score_%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source benchmarks/grn_multi/slurm/env.sh
GRNM_CORES=${SLURM_CPUS_PER_TASK:-8} Rscript benchmarks/grn_multi/03_score.R
Rscript benchmarks/grn_multi/04_summarise.R
