#!/usr/bin/env bash
## ======================================================================
## benchmarks/grn_multi/slurm/submit_all.sh
## The full run in one command (from the repository root, after the test run):
##   bash benchmarks/grn_multi/slurm/submit_all.sh
## 1. prepare.sbatch (waits for it): datasets + task chunks
## 2. methods.sbatch as an array, one job per chunk
## 3. score.sbatch after all of them (afterany: also when some fail)
## Running it again later submits only what has not finished (e.g. after
## a job was killed); GRNM_RETRY_ERRORS=1 also retries errors and timeouts.
## ======================================================================
set -euo pipefail
source benchmarks/grn_multi/slurm/env.sh
mkdir -p logs_grnm
sbatch --wait --partition="$GRNM_PARTITION" benchmarks/grn_multi/slurm/prepare.sbatch
N=$(( $(wc -l < "$GRNM_WORK/chunks.csv") - 1 ))
if (( N < 1 )); then echo "Nothing left to run."; exit 0; fi
JID=$(sbatch --parsable --partition="$GRNM_PARTITION" --array=1-"$N" \
      --cpus-per-task="$GRNM_CPUS" benchmarks/grn_multi/slurm/methods.sbatch)
SID=$(sbatch --parsable --partition="$GRNM_PARTITION" --dependency=afterany:"$JID" \
      benchmarks/grn_multi/slurm/score.sbatch)
echo "methods: array job $JID ($N chunks); scoring: job $SID after it."
echo "Progress:  squeue -u $USER ;  ls $GRNM_WORK/raw/*/ | wc -l"
