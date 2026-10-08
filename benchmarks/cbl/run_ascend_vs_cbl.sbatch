#!/usr/bin/env bash
## ======================================================================
## benchmarks/cbl/run_ascend_vs_cbl.sbatch
## The ASCEND vs CBL sweeps on CREATE, one configuration per array task
## (12 configurations: 5 values of n, 4 of d_x, 5 of d_z, sharing one point),
## 10 seeds each, both methods on the same node so times are comparable.
##
## Reviewer comments: 6 (CI tests vs cost per test), 5 (empirical scaling),
##   minor 9 (coverage)
## How to run (from the repository root):
##   sbatch benchmarks/cbl/run_ascend_vs_cbl.sbatch
## Results: benchmarks/cbl/results/shard_*.csv (one per configuration);
##   send Claude the benchmarks/cbl/results folder.
## ======================================================================
#SBATCH --job-name=ascend_cbl
#SBATCH --partition=cpu
#SBATCH --array=1-12
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=ascend_cbl_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source benchmarks/grn_multi/slurm/env.sh
export BENCH_SEEDS="${BENCH_SEEDS:-10}"
export BENCH_SHARD="${SLURM_ARRAY_TASK_ID}/12"
export BENCH_OUT="benchmarks/cbl/results/shard_${SLURM_ARRAY_TASK_ID}.csv"
Rscript benchmarks/cbl/run_ascend_vs_cbl.R
