#!/bin/bash
# =============================================================================
# submit.sh — one SLURM array task = one (combo, n) pair of hpc_run.R
#
# Quick uniform submission (same time/mem for every task):
#     mkdir -p logs results && sbatch submit.sh
#
# Preferred: let submit_all.sh override --array/--time/--mem per n-band:
#     bash submit_all.sh
#
# Command-line flags passed to sbatch OVERRIDE the #SBATCH lines below, which is
# exactly how submit_all.sh sizes each band.
# =============================================================================
#SBATCH --job-name=ascend
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-729           # default grid = 81 combos x 9 n-values;
                                # `Rscript hpc_run.R count` prints the true total,
                                # and submit_all.sh overrides this per band.
#SBATCH --cpus-per-task=1        # methods run serially; 1 core packs best
#SBATCH --mem=8G                 # generic default; large-n bands need more
#SBATCH --time=12:00:00          # generic default; large-n bands need more
# #SBATCH --partition=<your_partition>    # uncomment + set if required
# #SBATCH --account=<your_account>        # uncomment + set if required

set -uo pipefail                  # NB: no -e (module functions can return nonzero)

module load r/4.3.1               # <-- match `module avail r` on the cluster

mkdir -p logs results

echo "Host    : $(hostname)"
echo "Task    : ${SLURM_ARRAY_TASK_ID:-NA}  (job ${SLURM_ARRAY_JOB_ID:-NA})"
echo "Start   : $(date)"
echo "----------------------------------------------------------------------"

Rscript hpc_run.R "${SLURM_ARRAY_TASK_ID}"
rc=$?

echo "----------------------------------------------------------------------"
echo "End     : $(date)   (exit ${rc})"
exit ${rc}