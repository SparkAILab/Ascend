"""
benchmarks/grn_multi/run_task.py
Runs ONE Python method (GRNBoost2, RegDiffusion, PIDC, NOTEARS) on ONE dataset and
writes its raw output to WORK/raw/<method>/<dataset_id>.* in the same format as
run_task.R: <id>.score.csv.gz (foreground x foreground, [regulator, target]),
for NOTEARS also <id>.adj.csv.gz and <id>.anc.csv.gz, and <id>.meta.json
(written last; its presence marks the run as finished).

Reviewer comments: 1, 2, 6
How to run: normally called by run_chunk.sh under a wall-clock limit; by hand:
    python benchmarks/grn_multi/run_task.py <dataset_id> <method>
WORK defaults to benchmarks/grn_multi/work (work_test if GRNM_TEST=1); set GRNM_WORK
to override, exactly as for the R scripts.
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "1")        # single-threaded, like the R methods
os.environ.setdefault("TQDM_DISABLE", "1")   # no progress bars in the logs

import gzip
import json
import socket
import sys
import time
import traceback

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "py"))
from grn_methods import PY_METHODS  # noqa: E402

TEST = os.environ.get("GRNM_TEST", "0") == "1"
WORK = os.environ.get("GRNM_WORK", os.path.join(HERE, "work_test" if TEST else "work"))


def graph_to_ancestral(adj, keep):
    """Same as R/common.R::graph_to_ancestral: directed paths closed over all
    variables, an undirected edge not ordered by a path = 1 both ways, NA diagonal."""
    A = (np.asarray(adj) > 0)
    D = (A & ~A.T).astype(np.int64)
    U = A & A.T
    R = D.copy()
    while True:
        R2 = ((R + R @ R) > 0).astype(np.int64)
        if np.array_equal(R2, R):
            break
        R = R2
    R = R[np.ix_(keep, keep)].astype(float)
    Uk = U[np.ix_(keep, keep)] & (R == 0) & (R.T == 0)
    R[Uk] = 1.0
    np.fill_diagonal(R, np.nan)
    return R


def main(ds_id, method):
    if method not in PY_METHODS:
        sys.exit(f"unknown Python method: {method}")
    ddir = os.path.join(WORK, "data", ds_id)
    with open(os.path.join(ddir, "meta.json")) as f:
        meta_ds = json.load(f)
    zl, xl = list(meta_ds["zlabs"]), list(meta_ds["xlabs"])
    df = pd.read_csv(os.path.join(ddir, "data.csv.gz"))
    names = zl + xl
    data = df[names].to_numpy(dtype=np.float64)
    seed = int(meta_ds["seed"]) % (2 ** 31)
    np.random.seed(seed)

    out_dir = os.path.join(WORK, "raw", method)
    os.makedirs(out_dir, exist_ok=True)
    meta = {"dataset": ds_id, "method": method, "lang": "py",
            "host": socket.gethostname(), "python": sys.version.split()[0]}

    t0, c0 = time.perf_counter(), time.process_time()
    try:
        if method == "grnboost2":
            S = PY_METHODS[method](data, names, seed, targets=xl)
        elif method == "notears_tiered":
            S = PY_METHODS[method](data, names, seed, n_background=len(zl))
        else:
            S = PY_METHODS[method](data, names, seed)
        status, err = "ok", None
    except Exception as e:  # recorded, never re-raised: the chunk moves on
        status, err, S = "error", f"{type(e).__name__}: {e}", None
        traceback.print_exc()
    meta["runtime_s"] = time.perf_counter() - t0
    meta["cpu_s"] = time.process_time() - c0
    meta["status"] = status
    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    if err:
        meta["error"] = err
    else:
        ix = [names.index(g) for g in xl]
        outs = S if isinstance(S, dict) else {"score": S}
        if "adj" in outs:                     # causal output: add the ancestral matrix
            outs["anc"] = graph_to_ancestral(outs["adj"], ix)
        for what, M in outs.items():
            Mx = M if what == "anc" else M[np.ix_(ix, ix)]
            df_out = pd.DataFrame(Mx, index=xl, columns=xl)
            df_out.index.name = "row"
            with gzip.open(os.path.join(out_dir, f"{ds_id}.{what}.csv.gz"), "wt") as f:
                df_out.to_csv(f, float_format="%.8g", na_rep="NA")
        meta["outputs"] = sorted(outs)
    with open(os.path.join(out_dir, f"{ds_id}.meta.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(f"{ds_id} {method} {status} {meta['runtime_s']:.1f}s" + (f" : {err}" if err else ""))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: python run_task.py <dataset_id> <method>")
    main(sys.argv[1], sys.argv[2])
