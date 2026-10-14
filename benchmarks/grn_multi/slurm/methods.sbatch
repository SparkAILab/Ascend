#!/usr/bin/env bash
## ======================================================================
## benchmarks/grn_multi/slurm/methods.sbatch
## One array job = one chunk of tasks (WORK/tasks.csv), run on
## --cpus-per-task parallel workers. Resumable: finished tasks are skipped.
## How to run: called by submit_all.sh
##   (sbatch --array=1-<n_chunks> --cpus-per-task=$GRNM_CPUS this file)
## ======================================================================
#SBATCH --job-name=grnm
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs_grnm/grnm_%A_%a.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source benchmarks/grn_multi/slurm/env.sh
bash benchmarks/grn_multi/run_chunk.sh "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK"
