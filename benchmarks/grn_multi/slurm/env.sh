## ======================================================================
## benchmarks/grn_multi/slurm/env.sh
## Environment for every job of the benchmark. EDIT THE FIRST BLOCK ONCE for
## your account on King's CREATE (or any SLURM cluster); every other script
## sources this file.
## ======================================================================

# --- edit these ---------------------------------------------------------
# Modules: see `module avail r` and `module avail python` on CREATE and
# put the exact names here (R >= 4.3, Python >= 3.9).
module load r 2>/dev/null || true
module load python 2>/dev/null || true
export GRNM_WORK="${GRNM_WORK:-/scratch/users/$USER/grn_multi}"   # all generated files
export GRNM_VENV="${GRNM_VENV:-$HOME/venvs/grn_multi}"             # Python virtual env
export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/grn_multi}"             # R package library
# R's default personal library too, so packages installed by hand with
# install.packages() are found (set the version to your R: R --version)
export R_LIBS_USER="$R_LIBS_USER:$HOME/R/x86_64-pc-linux-gnu-library/4.5"
export GRNM_PARTITION="${GRNM_PARTITION:-cpu}"                     # CREATE partition (48 h limit)
export GRNM_CPUS="${GRNM_CPUS:-8}"                                 # workers per array job
# ------------------------------------------------------------------------

export GRNM_SERGIO="${GRNM_SERGIO:-$GRNM_WORK/SERGIO}"
export GRNM_PYTHON="$GRNM_VENV/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p "${R_LIBS_USER%%:*}" "$GRNM_WORK"
# repository root (this file is benchmarks/grn_multi/slurm/env.sh)
GRNM_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
export GRNM_REPO
