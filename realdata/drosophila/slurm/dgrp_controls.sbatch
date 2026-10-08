#!/usr/bin/env bash
## ======================================================================
## realdata/drosophila/slurm/dgrp_controls.sbatch
## One SLURM array task = one or more DGRP control runs of dgrp_controls.R
## (array index i does runs FIRST+i*STEP .. FIRST+i*STEP+STEP-1).
##
## Reviewer comments: 7 (DGRP controls)
## How to run: called by submit_dgrp.sh; nothing to run by hand.
## ======================================================================
#SBATCH --job-name=dgrp_ctrl
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=dgrp_ctrl_%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source benchmarks/grn_multi/slurm/env.sh
STEP="${DGRP_STEP:-1}"
first=$(( ${DGRP_FIRST:-1} + ${SLURM_ARRAY_TASK_ID:-0} * STEP ))
for id in $(seq "$first" $(( first + STEP - 1 ))); do
  [ "$id" -le "${DGRP_NRUNS:-241}" ] || break
  Rscript realdata/drosophila/dgrp_controls.R "$id"
done
