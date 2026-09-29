"""
benchmarks/grn_multi/py/make_sergio.py
Builds the BEELINE-style arm from SERGIO's published steady-state datasets
(Dibaeinia & Sinha 2020, Cell Systems; https://github.com/PayamDiba/SERGIO).

Reviewer comment: 3 (BEELINE generators and benchmark networks)

Why SERGIO and how the two tiers are defined
    BEELINE's generators simulate expression only, so no dataset has a genotype
    layer. SERGIO's steady-state datasets DS1 (100 genes, E. coli sub-network) and
    DS2 (400 genes, yeast sub-network) are GeneNetWeaver networks simulated with a
    chemical Langevin equation, 9 cell types x 300 cells, 15 replicates each. In
    SERGIO the master regulators (genes nothing regulates) drive the whole network
    and nothing feeds back into them, so they satisfy the two-tier assumption
    exactly: they are the background tier Z, and every other gene is foreground X.
    The data are non-linear (Hill kinetics) and non-Gaussian, unlike the SEM arms.

Versions written for each replicate
    clean : SERGIO's de-noised expression, log1p
    noisy : SERGIO's technical-noise model with the parameters of SERGIO's own
            tutorial (outlier genes p = 0.01, lognormal(0.8, 1); library size
            lognormal(4.6, 0.4); dropout logistic shape 6.5 at the 82nd percentile;
            Poisson UMI counts), then library-size normalisation and log1p.
    Genes that end up constant are dropped (recorded in meta.json).

Usage
    git clone https://github.com/PayamDiba/SERGIO.git   (once)
    python benchmarks/grn_multi/py/make_sergio.py --sergio /path/to/SERGIO \\
        [--ds DS1,DS2] [--noise clean,noisy] [--reps 15]
Writes WORK/data/sergio_<DS>_<noise>_rep<NNN>/ in the same format as
01_make_datasets.R. Skips datasets that already exist.
"""
import argparse
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TEST = os.environ.get("GRNM_TEST", "0") == "1"
WORK = os.environ.get("GRNM_WORK", os.path.join(os.path.dirname(HERE), "work_test" if TEST else "work"))

DS_DIRS = {"DS1": ("De-noised_100G_9T_300cPerT_4_DS1", 4, 100),
           "DS2": ("De-noised_400G_9T_300cPerT_5_DS2", 5, 400),
           "DS3": ("De-noised_1200G_9T_300cPerT_6_DS3", 6, 1200)}
N_BINS, N_SC = 9, 300


def read_grn(ds_path, cid, G):
    """Edges from SERGIO's Interaction file: target, #regs, regs..., K..., coop..."""
    A = np.zeros((G, G), dtype=int)                      # A[reg, target] = 1
    with open(os.path.join(ds_path, f"Interaction_cID_{cid}.txt")) as f:
        for line in f:
            v = [float(x) for x in line.strip().split(",") if x != ""]
            if not v:
                continue
            tgt, k = int(v[0]), int(v[1])
            for r in v[2:2 + k]:
                if int(r) != tgt:
                    A[int(r), tgt] = 1
    mrs = []
    with open(os.path.join(ds_path, f"Regs_cID_{cid}.txt")) as f:
        for line in f:
            if line.strip():
                mrs.append(int(float(line.split(",")[0])))
    return A, sorted(mrs)


def closure(A):
    R = (A > 0).astype(np.int64)
    while True:
        R2 = ((R + R @ R) > 0).astype(np.int64)
        if np.array_equal(R2, R):
            break
        R = R2
    np.fill_diagonal(R, 0)
    return R


def add_noise(expr, G, seed, sergio_dir):
    """SERGIO's technical-noise pipeline (tutorial parameters) on a genes x cells matrix."""
    sys.path.insert(0, sergio_dir)
    from SERGIO.sergio import sergio
    np.random.seed(seed)
    sim = sergio(number_genes=G, number_bins=N_BINS, number_sc=N_SC, noise_params=1,
                 decays=0.8, sampling_state=15, noise_type="dpd")
    e3 = np.array(np.split(expr, N_BINS, axis=1))       # bins x genes x cells
    e_o = sim.outlier_effect(e3, outlier_prob=0.01, mean=0.8, scale=1)
    _, e_ol = sim.lib_size_effect(e_o, mean=4.6, scale=0.4)
    ind = sim.dropout_indicator(e_ol, shape=6.5, percentile=82)
    counts = sim.convert_to_UMIcounts(np.multiply(ind, e_ol))
    counts = np.concatenate(counts, axis=1).astype(np.float64)
    lib = counts.sum(axis=0)
    keep = lib > 0
    counts = counts[:, keep]
    norm = counts / lib[keep] * np.median(lib[keep])
    return np.log1p(norm), int((~keep).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sergio", required=True, help="path to a clone of PayamDiba/SERGIO")
    ap.add_argument("--ds", default="DS1,DS2")
    ap.add_argument("--noise", default="clean,noisy")
    ap.add_argument("--reps", type=int, default=15)
    a = ap.parse_args()
    os.makedirs(os.path.join(WORK, "data"), exist_ok=True)
    n_new = 0
    for di, ds in enumerate(a.ds.split(",")):
        sub, cid, G = DS_DIRS[ds]
        ds_path = os.path.join(a.sergio, "data_sets", sub)
        A, mrs = read_grn(ds_path, cid, G)
        for ni, noise in enumerate(a.noise.split(",")):
            for rep in range(1, a.reps + 1):
                ds_id = f"sergio_{ds}_{noise}_rep{rep:03d}"
                out = os.path.join(WORK, "data", ds_id)
                if os.path.exists(os.path.join(out, "meta.json")):
                    continue
                seed = 9_000_000 + 100_000 * di + 10_000 * ni + rep
                raw = pd.read_csv(os.path.join(ds_path, f"simulated_noNoise_{rep - 1}.csv"),
                                  index_col=0).to_numpy(dtype=np.float64)   # genes x cells
                assert raw.shape == (G, N_BINS * N_SC), raw.shape
                dropped_cells = 0
                if noise == "clean":
                    E = np.log1p(raw)
                else:
                    E, dropped_cells = add_noise(raw, G, seed, a.sergio)
                sd = E.std(axis=1)
                genes = [g for g in range(G) if sd[g] > 1e-8]
                dropped = [g for g in range(G) if sd[g] <= 1e-8]
                z_genes = [g for g in mrs if g in genes]
                x_genes = [g for g in genes if g not in mrs]
                zl = [f"z{i + 1}" for i in range(len(z_genes))]
                xl = [f"x{i + 1}" for i in range(len(x_genes))]
                R = closure(A)                         # ancestry over the FULL network
                ix = np.array(x_genes)
                T_dir = A[np.ix_(ix, ix)]
                T_anc = R[np.ix_(ix, ix)]
                os.makedirs(out, exist_ok=True)
                data = pd.DataFrame(E[z_genes + x_genes, :].T, columns=zl + xl)
                with gzip.open(os.path.join(out, "data.csv.gz"), "wt") as f:
                    data.to_csv(f, index=False, float_format="%.7g")
                for name, M in (("truth_direct.csv", T_dir), ("truth_anc.csv", T_anc)):
                    df = pd.DataFrame(M, index=xl, columns=xl)
                    df.index.name = "row"
                    df.to_csv(os.path.join(out, name))
                meta = {"id": ds_id, "arm": "sergio", "arm_key": f"sergio_{ds}", "rep": rep,
                        "seed": seed, "cell": {"ds": ds, "noise": noise},
                        "n": int(E.shape[1]), "d_z": len(zl), "d_x": len(xl),
                        "zlabs": zl, "xlabs": xl,
                        "sergio_gene_id": {**dict(zip(zl, z_genes)), **dict(zip(xl, x_genes))},
                        "dropped_constant_genes": dropped, "dropped_empty_cells": dropped_cells,
                        "n_direct": int(T_dir.sum()), "n_anc": int(T_anc.sum()),
                        "generator": f"SERGIO {sub}/simulated_noNoise_{rep - 1}.csv"}
                with open(os.path.join(out, "meta.json"), "w") as f:
                    json.dump(meta, f, indent=1)
                n_new += 1
                print(f"{ds_id}: n={meta['n']} d_z={meta['d_z']} d_x={meta['d_x']} "
                      f"direct={meta['n_direct']} ancestral={meta['n_anc']}")
    print(f"{n_new} new SERGIO datasets in {os.path.join(WORK, 'data')}")


if __name__ == "__main__":
    main()
