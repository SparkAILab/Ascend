#!/bin/bash
# =============================================================================
# submit_all.sh — submit the array in n-bands, each with its own --time/--mem,
#                 so cheap small-n tasks aren't billed the wall-time and memory
#                 the n=131072 tasks need.
#
# The grid size is read from hpc_run.R, so this keeps working if you edit
# make_grid() or N_VEC. Layout (must match hpc_run.R): the NCOMBO combos at
# N_VEC[k] occupy tasks ((k-1)*NCOMBO + 1) .. (k*NCOMBO).
#
# Usage:
#     bash submit_all.sh
# =============================================================================
set -uo pipefail

mkdir -p logs results

# Ask hpc_run.R for the current grid dimensions: "NCOMBO NVALS NTASKS".
read -r NCOMBO NVALS NTASKS < <(Rscript hpc_run.R count)
if [[ -z "${NCOMBO:-}" ]]; then
  echo "Could not read grid size from 'Rscript hpc_run.R count'. Aborting." >&2
  exit 1
fi
echo "Grid: NCOMBO=${NCOMBO}  n-values=${NVALS}  total tasks=${NTASKS}"

# Per-band resources, indexed by n (band 1 = smallest n). Verify against your
# partition limits. If NVALS exceeds this table, the last entry is reused.
BAND_TIME=(04:00:00 04:00:00 04:00:00 06:00:00 08:00:00 12:00:00 18:00:00 24:00:00 36:00:00)
BAND_MEM=(4G       4G       4G       6G       8G       12G      16G      24G      32G)
# Optional concurrency throttle per band (max simultaneous tasks); "" = unlimited.
BAND_THROTTLE=("" "" "" "" "" "" "" "" "")

last=$(( ${#BAND_TIME[@]} - 1 ))
for (( k=1; k<=NVALS; k++ )); do
  idx=$(( k-1 )); [[ $idx -gt $last ]] && idx=$last
  lo=$(( (k-1)*NCOMBO + 1 ))
  hi=$(( k*NCOMBO ))
  arr="${lo}-${hi}"
  thr="${BAND_THROTTLE[$idx]}"; [[ -n "$thr" ]] && arr="${arr}%${thr}"
  echo "band ${k}/${NVALS}  array=${arr}  time=${BAND_TIME[$idx]}  mem=${BAND_MEM[$idx]}"
  sbatch --array="${arr}" --time="${BAND_TIME[$idx]}" --mem="${BAND_MEM[$idx]}" submit.sh
done

echo "Submitted. Watch with:  squeue -u \$USER"