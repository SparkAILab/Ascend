#!/usr/bin/env bash
## ======================================================================
## benchmarks/grn_multi/run_chunk.sh
## Runs every task of one chunk (see 02_make_tasks.R) on N parallel
## workers. Each task runs under `timeout` with its own budget; a task that
## hits it gets a meta.json with status "timeout". Tasks that already have
## a meta.json are skipped, so a chunk that was killed can simply be
## submitted again and only the unfinished tasks run.
##
## How to run:  bash benchmarks/grn_multi/run_chunk.sh <chunk> [workers]
##   (normally from slurm/methods.sbatch with the array index as <chunk>)
## Logs: WORK/logs/<method>/<dataset>.log
## ======================================================================
set -uo pipefail

CHUNK=${1:?usage: run_chunk.sh <chunk> [workers]}
WORKERS=${2:-${SLURM_CPUS_PER_TASK:-1}}
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${GRNM_TEST:-0}" == "1" ]]; then DEF_WORK="$HERE/work_test"; else DEF_WORK="$HERE/work"; fi
export GRNM_WORK="${GRNM_WORK:-$DEF_WORK}"
export HERE
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYTHON="${GRNM_PYTHON:-python3}"; export PYTHON

TASKS="$GRNM_WORK/tasks.csv"
[[ -f "$TASKS" ]] || { echo "no $TASKS: run 02_make_tasks.R first" >&2; exit 1; }

run_one() {
  # one CSV line: task,chunk,dataset,method,lang,budget_s,est_s
  IFS=, read -r task chunk ds method lang budget est <<< "$1"
  local meta="$GRNM_WORK/raw/$method/$ds.meta.json"
  [[ -f "$meta" ]] && return 0
  mkdir -p "$GRNM_WORK/raw/$method" "$GRNM_WORK/logs/$method"
  local log="$GRNM_WORK/logs/$method/$ds.log"
  local t0=$SECONDS rc
  if [[ "$lang" == "R" ]]; then
    timeout --kill-after=120 "$budget" Rscript "$HERE/run_task.R" "$ds" "$method" > "$log" 2>&1; rc=$?
  else
    timeout --kill-after=120 "$budget" "$PYTHON" "$HERE/run_task.py" "$ds" "$method" > "$log" 2>&1; rc=$?
  fi
  if [[ ! -f "$meta" ]]; then
    # no meta.json: the method was stopped (timeout) or crashed hard (e.g. out of memory)
    local status="crashed"
    [[ $rc -eq 124 || $rc -eq 137 ]] && status="timeout"
    printf '{"dataset":"%s","method":"%s","lang":"%s","status":"%s","runtime_s":%d,"exit_code":%d,"budget_s":%d,"host":"%s"}\n' \
      "$ds" "$method" "$lang" "$status" $((SECONDS - t0)) "$rc" "$budget" "$(hostname)" > "$meta"
  fi
  echo "[$(date +%H:%M:%S)] task $task $method $ds -> $(grep -o '"status": *"[a-z]*"' "$meta" | head -1) ($((SECONDS - t0)) s)"
}
export -f run_one

echo "chunk $CHUNK on $WORKERS workers, WORK=$GRNM_WORK, started $(date)"
awk -F, -v c="$CHUNK" 'NR > 1 && $2 == c' "$TASKS" | \
  xargs -d '\n' -P "$WORKERS" -I{} bash -c 'run_one "$@"' _ {}
echo "chunk $CHUNK finished $(date)"
