#!/bin/bash
#SBATCH --job-name=ascend_bench
#SBATCH --array=1-81
#SBATCH --output=logs/job_%a_%j.out
#SBATCH --error=logs/job_%a_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

# =============================================================================
# SLURM array job for ASCEND benchmark on KCL CREATE HPC
#
# Submit with:   sbatch submit.sh
# Monitor with:  squeue -u $USER
# Check logs:    tail -f logs/job_1_*.out
#
# Notes:
#   - 81 array jobs, one per parameter combination
#   - 48h wall time (generous; most jobs will finish in 4-12h)
#   - 32GB RAM: sufficient for d_x=80, d_z=240, n=131072
#     (data matrix ~335MB, methods add overhead; 32GB is safe)
#   - 1 CPU per job (all methods are single-threaded here;
#     GENIE3 nCores=1 is set in hpc_run.R)
#   - Jobs are restartable: completed replicates are skipped on rerun
#
# If some jobs run out of time, resubmit the same script — completed
# replicates will be skipped automatically.
# =============================================================================

# Fail immediately on any error in this shell script
set -euo pipefail

# Create log directory if it doesn't exist
mkdir -p logs results

# Load R module — adjust version to what is available on CREATE
# Check available versions with: module avail R
module load r/4.5.1-gcc-13.2.0-withx-rmath-standalone-python-3.11.6

# Print job info for debugging
echo "=============================="
echo "Job array index : ${SLURM_ARRAY_TASK_ID}"
echo "Job ID          : ${SLURM_JOB_ID}"
echo "Node            : $(hostname)"
echo "Start time      : $(date)"
echo "Working dir     : $(pwd)"
echo "=============================="

# Run the R script with this array index
Rscript hpc_run.R ${SLURM_ARRAY_TASK_ID}

echo "=============================="
echo "End time: $(date)"
echo "=============================="