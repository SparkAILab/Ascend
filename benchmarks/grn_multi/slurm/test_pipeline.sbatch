#!/usr/bin/env bash
## ======================================================================
## benchmarks/grn_multi/slurm/test_pipeline.sbatch
## The test run: every method on one or two datasets of EVERY arm (7
## datasets), each run capped at 2 h, then scoring and the summary. It
## checks that everything works on the cluster before the full run, and
## its timings are used to pack the full run into jobs.
## How to run (from the repository root):
##   sbatch benchmarks/grn_multi/slurm/test_pipeline.sbatch
## Then read  $GRNM_WORK/../grn_multi_test/results/TEST_REPORT.txt
## (every method should be "ok" on every arm; "timeout" means it needs
## longer than 2 h there and gets its full budget in the real run).
## ======================================================================
#SBATCH --job-name=grnm_test
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --output=grnm_test_%j.out
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source benchmarks/grn_multi/slurm/env.sh
export GRNM_TEST=1
export GRNM_WORK="$(dirname "$GRNM_WORK")/grn_multi_test"
W=${SLURM_CPUS_PER_TASK:-8}

Rscript benchmarks/grn_multi/01_make_datasets.R
"$GRNM_PYTHON" benchmarks/grn_multi/py/make_sergio.py --sergio "$GRNM_SERGIO" --ds DS1,DS2 --reps 1
GRNM_CPUS=$W Rscript benchmarks/grn_multi/02_make_tasks.R
for c in $(tail -n +2 "$GRNM_WORK/chunks.csv" | cut -d, -f1); do
  bash benchmarks/grn_multi/run_chunk.sh "$c" "$W"
done
GRNM_CORES=$W Rscript benchmarks/grn_multi/03_score.R
Rscript benchmarks/grn_multi/04_summarise.R
Rscript benchmarks/grn_multi/05_test_report.R
# keep the timings next to the code so 02_make_tasks.R finds them for the full run
cp "$GRNM_WORK/results/runs.csv" benchmarks/grn_multi/pilot_runs.csv
