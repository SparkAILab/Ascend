#!/usr/bin/env bash
## ======================================================================
## realdata/drosophila/slurm/submit_dgrp.sh
## Submits the whole DGRP control analysis on CREATE: run 1 (the observed
## network, which also caches the filtered data), then the 240 control runs
## as an array once run 1 has finished, then the summary.
##
## Reviewer comments: 7 (DGRP controls), minor 7 (filtering counts)
## How to run (from the repository root, on a login node):
##   export DGRP_DATA=/scratch/users/$USER/dgrp   # folder with the two data files
##   bash realdata/drosophila/slurm/submit_dgrp.sh
## Results: $DGRP_OUT (default /scratch/users/$USER/dgrp/results); read
##   DGRP_CONTROLS.md there, then send Claude the results folder.
## ======================================================================
set -euo pipefail
source benchmarks/grn_multi/slurm/env.sh
export DGRP_DATA="${DGRP_DATA:-/scratch/users/$USER/dgrp}"
export DGRP_OUT="${DGRP_OUT:-$DGRP_DATA/results}"
for f in dgrp2.tgeno.txt GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz; do
  [ -f "$DGRP_DATA/$f" ] || { echo "Missing $DGRP_DATA/$f (see realdata/drosophila/dgrp_data.R for where to download it)" >&2; exit 1; }
done
mkdir -p "$DGRP_OUT"
NRUNS=$(Rscript realdata/drosophila/dgrp_controls.R | head -1 | awk '{print $1}')
STEP=4                                    # control runs per array task
LAST=$(( (NRUNS - 2) / STEP ))            # array over runs 2..NRUNS
P="--partition=$GRNM_PARTITION"
# run 1 reads the full genotype file once and caches the filtered data
obs=$(sbatch --parsable $P --export=ALL --mem=64G --time=24:00:00 --job-name=dgrp_obs --output=dgrp_obs_%j.out \
        --wrap "source benchmarks/grn_multi/slurm/env.sh && Rscript realdata/drosophila/dgrp_controls.R 1")
arr=$(sbatch --parsable $P --export=ALL,DGRP_FIRST=2,DGRP_STEP=$STEP,DGRP_NRUNS=$NRUNS \
        --dependency=afterok:$obs --array=0-$LAST realdata/drosophila/slurm/dgrp_controls.sbatch)
sumj=$(sbatch --parsable $P --export=ALL --dependency=afterany:$arr --time=2:00:00 --mem=8G \
        --job-name=dgrp_sum --output=dgrp_sum_%j.out \
        --wrap "source benchmarks/grn_multi/slurm/env.sh && Rscript realdata/drosophila/dgrp_controls_summarise.R")
echo "Submitted: observed run $obs, controls array $arr (0-$LAST), summary $sumj"
echo "Watch with: squeue -u $USER ; results in $DGRP_OUT"
