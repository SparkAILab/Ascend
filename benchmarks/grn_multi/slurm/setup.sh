#!/usr/bin/env bash
## ======================================================================
## benchmarks/grn_multi/slurm/setup.sh
## One-time setup: Python virtual env, R packages, SERGIO datasets.
## How to run (login node, from the repository root, after editing env.sh):
##   bash benchmarks/grn_multi/slurm/setup.sh
## If the login node forbids long installs, run it inside an interactive
## job instead:  srun -p cpu --cpus-per-task=4 --mem=8G --time=2:00:00 --pty bash
## ======================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cd "$GRNM_REPO"

echo "== Python env in $GRNM_VENV"
[[ -x "$GRNM_PYTHON" ]] || python3 -m venv "$GRNM_VENV"
"$GRNM_PYTHON" -m pip install --upgrade pip
# CPU-only torch keeps the install small; RegDiffusion runs on CPU here
"$GRNM_PYTHON" -m pip install torch --index-url https://download.pytorch.org/whl/cpu
"$GRNM_PYTHON" -m pip install -r benchmarks/grn_multi/requirements.txt
"$GRNM_PYTHON" -c "import arboreto, regdiffusion, astropy, sklearn; print('python ok')"

echo "== R packages in ${R_LIBS_USER%%:*}"
# cmake builds the libuv bundled with the fs package (no system libuv on CREATE)
module load cmake 2>/dev/null || true
Rscript benchmarks/grn_multi/install_R_packages.R || R_FAILED=1

echo "== SERGIO datasets in $GRNM_SERGIO (about 1.7 GB)"
[[ -d "$GRNM_SERGIO/data_sets" ]] || git clone --depth 1 https://github.com/PayamDiba/SERGIO.git "$GRNM_SERGIO"

if [[ "${R_FAILED:-0}" == 1 ]]; then
  echo "SETUP INCOMPLETE: some R packages did not install (see above). Fix those, then run setup.sh again." >&2
  exit 1
fi
echo "== tests"
Rscript tests/test_metrics.R > /dev/null && Rscript tests/test_grn_multi.R
echo "Setup complete. Next: sbatch benchmarks/grn_multi/slurm/test_pipeline.sbatch"
