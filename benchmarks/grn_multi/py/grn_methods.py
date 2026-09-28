"""
benchmarks/grn_multi/py/grn_methods.py
The Python methods of the benchmark: GRNBoost2, RegDiffusion, PIDC, and
NOTEARS with the tier constraint.

Reviewer comments: 1 (GRNBoost2/SCENIC, PIDC, a modern deep-learning method,
a direction-aware comparator), 10 (tier knowledge in a continuous-optimisation
DAG learner)

Every function takes
    data   : n x p float array, columns = all background AND foreground variables
    names  : the p column names
    seed   : integer seed
and returns a p x p score matrix S with S[a, b] = confidence that a regulates b
(rows = regulators). The runner keeps the foreground block.

Implementation notes
- GRNBoost2 calls arboreto's own per-target routine (arboreto.core.infer_partial_network)
  with arboreto's GRNBoost2 settings, one target at a time and without dask, so it
  runs single-threaded and deterministically. The result is identical to
  arboreto.algo.grnboost2 up to the random seed.
- RegDiffusion is the published package (pip install regdiffusion), on CPU with
  its default hyper-parameters.
- NOTEARS (Zheng et al. 2018) is the authors' reference linear implementation
  (github.com/xunzheng/notears, notears/linear.py, Apache-2.0), with the tier
  constraint added as L-BFGS-B bounds that fix every foreground -> background
  weight at 0. The data are standardised first, because unstandardised
  simulated data let NOTEARS read the causal order off the variances
  (Reisach et al. 2021, "Beware of the simulated DAG"); ASCEND's tests are
  scale-free, so this puts both on the same footing. Default settings:
  lambda1 = 0.1, least-squares loss, weight threshold 0.3.
  Returns a dict with the thresholded graph and |W| before thresholding as
  a ranking score.
- PIDC is a faithful port of NetworkInference.jl (Chan et al. 2017; the reference
  implementation used by BEELINE) with its defaults: Bayesian-blocks
  discretisation, maximum-likelihood probabilities, base-2 logarithms, PUC scores
  and the gamma-CDF context step. Julia is not needed.
"""

import numpy as np


# --------------------------------------------------------------------------
# GRNBoost2
# --------------------------------------------------------------------------
def grnboost2(data, names, seed, targets=None):
    from arboreto import core
    names = list(names)
    targets = names if targets is None else list(targets)
    tf_matrix = np.asarray(data, dtype=np.float64)
    S = np.zeros((len(names), len(names)))
    idx = {g: i for i, g in enumerate(names)}
    for t in targets:
        j = idx[t]
        links = core.infer_partial_network(
            regressor_type="GBM", regressor_kwargs=core.SGBM_KWARGS,
            tf_matrix=tf_matrix, tf_matrix_gene_names=names,
            target_gene_name=t, target_gene_expression=tf_matrix[:, j],
            include_meta=False,
            early_stop_window_length=core.EARLY_STOP_WINDOW_LENGTH, seed=seed)
        if len(links) == 0:
            # arboreto's retry() swallows errors and returns an empty frame
            raise RuntimeError(f"GRNBoost2 regression failed for target {t}")
        for tf, imp in zip(links["TF"], links["importance"]):
            S[idx[tf], j] = imp
    return S


# --------------------------------------------------------------------------
# RegDiffusion (deep learning; Zhu & Slonim 2024)
# --------------------------------------------------------------------------
def regdiffusion(data, names, seed):
    import torch
    import regdiffusion as rd
    torch.set_num_threads(1)
    X = np.asarray(data, dtype=np.float32)
    trainer = rd.RegDiffusionTrainer(X, device="cpu", seed=int(seed))
    trainer.train()
    adj = np.asarray(trainer.get_adj(), dtype=np.float64)      # [regulator, target]
    return np.abs(adj)


# --------------------------------------------------------------------------
# PIDC (Chan et al. 2017), port of NetworkInference.jl
# --------------------------------------------------------------------------
def _discretize(x):
    """Bayesian-blocks bins (Scargle 2013, p0 = 0.05), as Discretizers.jl;
    uniform-width 10 bins if that fails, as NetworkInference.jl does."""
    from astropy.stats import bayesian_blocks
    x = np.asarray(x, dtype=np.float64)
    try:
        edges = bayesian_blocks(x, fitness="events", p0=0.05)
        if len(edges) < 2:
            raise ValueError
    except Exception:
        edges = np.linspace(x.min(), x.max(), 11)
    k = len(edges) - 1
    b = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, k - 1)
    return b, k


def _xlogy_ratio(p, q):
    """p * log2(p / q) with 0 for non-finite terms (remove_non_finite)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        v = p * np.log2(p / q)
    v[~np.isfinite(v)] = 0.0
    return v


def pidc(data, names, seed=None):
    X = np.asarray(data, dtype=np.float64)
    n, p = X.shape
    bins, nb = zip(*(_discretize(X[:, i]) for i in range(p)))
    nb = np.array(nb)
    B = int(nb.max())
    marg = np.zeros((p, B))
    for i in range(p):
        marg[i, :nb[i]] = np.bincount(bins[i], minlength=nb[i]) / n

    # pairwise MI and specific information SI[src, tgt, s] = I_spec(tgt = s; src)
    MI = np.zeros((p, p))
    SI = np.zeros((p, p, B))
    for i in range(p):
        for j in range(i + 1, p):
            joint = np.zeros((nb[i], nb[j]))
            np.add.at(joint, (bins[i], bins[j]), 1.0)
            joint /= n
            p1 = joint.sum(axis=1, keepdims=True)      # node i
            p2 = joint.sum(axis=0, keepdims=True)      # node j
            mi = _xlogy_ratio(joint, p1 * p2).sum()
            MI[i, j] = MI[j, i] = mi
            with np.errstate(divide="ignore", invalid="ignore"):
                t1 = (joint / p2) * np.log2(joint / (p1 * p2))
                t2 = (joint / p1) * np.log2(joint / (p1 * p2))
            t1[~np.isfinite(t1)] = 0.0
            t2[~np.isfinite(t2)] = 0.0
            SI[i, j, :nb[j]] = t1.sum(axis=0)          # i about the states of j
            SI[j, i, :nb[i]] = t2.sum(axis=1)          # j about the states of i

    # PUC: for target k and sources i, j:
    #   redundancy R = sum_s p_k(s) min(SI[i,k,s], SI[j,k,s])
    #   u(i -> k | j) = (MI[i,k] - R) / MI[i,k], set to 0 if negative or not finite
    # puc[a, b] sums u over every third variable, with b as target and with a as target.
    U = np.zeros((p, p))
    for k in range(p):
        S = SI[:, k, :]                                 # (p, B)
        R = np.zeros((p, p))
        for s in range(nb[k]):
            if marg[k, s] > 0:
                R += marg[k, s] * np.minimum.outer(S[:, s], S[:, s])
        mik = MI[:, k][:, None]                         # MI[i, k] for row i
        with np.errstate(divide="ignore", invalid="ignore"):
            u = (mik - R) / mik
        u[~np.isfinite(u) | (u < 0)] = 0.0
        u[k, :] = 0.0
        u[:, k] = 0.0
        np.fill_diagonal(u, 0.0)
        U[:, k] = u.sum(axis=1)
    puc = U + U.T

    # context step: gamma CDF of each gene's scores, CLR z-score if the fit fails
    from scipy import stats
    W = np.zeros((p, p))
    fits = []
    for i in range(p):
        sc = np.delete(puc[:, i], i)
        try:
            if np.any(sc <= 0) or np.var(sc) == 0:
                raise ValueError
            a, loc, scale = stats.gamma.fit(sc, floc=0)
            if not (np.isfinite(a) and np.isfinite(scale)):
                raise ValueError
            fits.append((a, scale))
        except Exception:
            fits.append(None)
    for i in range(p):
        for j in range(i + 1, p):
            sc_ij = puc[i, j]
            if fits[i] is not None and fits[j] is not None:
                w = stats.gamma.cdf(sc_ij, fits[i][0], scale=fits[i][1]) + \
                    stats.gamma.cdf(sc_ij, fits[j][0], scale=fits[j][1])
            else:
                si = np.delete(puc[:, i], i)
                sj = np.delete(puc[:, j], j)
                di, dj = sc_ij - si.mean(), sc_ij - sj.mean()
                vi, vj = si.var(ddof=1), sj.var(ddof=1)
                w = np.sqrt((0 if vi == 0 or di < 0 else di ** 2 / vi) +
                            (0 if vj == 0 or dj < 0 else dj ** 2 / vj))
            W[i, j] = W[j, i] = w
    return W


# --------------------------------------------------------------------------
# NOTEARS, linear, with tier bounds (reference implementation + bounds)
# --------------------------------------------------------------------------
def notears_linear(X, lambda1=0.1, max_iter=100, h_tol=1e-8, rho_max=1e16,
                   w_threshold=0.3, forbid=None):
    """Reference NOTEARS (least-squares loss). forbid: d x d bool, True = edge
    i -> j fixed at 0. Returns (W thresholded, W raw); W[i, j] != 0 is i -> j."""
    import scipy.linalg as slin
    import scipy.optimize as sopt
    n, d = X.shape

    def _adj(w):
        return (w[:d * d] - w[d * d:]).reshape([d, d])

    def _loss(W):
        R = X - X @ W
        return 0.5 / n * (R ** 2).sum(), -1.0 / n * X.T @ R

    def _h(W):
        E = slin.expm(W * W)
        return np.trace(E) - d, E.T * W * 2

    def _func(w):
        W = _adj(w)
        loss, G_loss = _loss(W)
        h, G_h = _h(W)
        obj = loss + 0.5 * rho * h * h + alpha * h + lambda1 * w.sum()
        G_smooth = G_loss + (rho * h + alpha) * G_h
        return obj, np.concatenate((G_smooth + lambda1, -G_smooth + lambda1), axis=None)

    if forbid is None:
        forbid = np.zeros((d, d), dtype=bool)
    fb = forbid | np.eye(d, dtype=bool)
    bnds = [(0, 0) if fb[i, j] else (0, None) for _ in range(2) for i in range(d) for j in range(d)]
    X = X - X.mean(axis=0, keepdims=True)
    w_est, rho, alpha, h = np.zeros(2 * d * d), 1.0, 0.0, np.inf
    for _ in range(max_iter):
        w_new, h_new = None, None
        while rho < rho_max:
            sol = sopt.minimize(_func, w_est, method="L-BFGS-B", jac=True, bounds=bnds)
            w_new = sol.x
            h_new, _ = _h(_adj(w_new))
            if h_new > 0.25 * h:
                rho *= 10
            else:
                break
        w_est, h = w_new, h_new
        alpha += rho * h
        if h <= h_tol or rho >= rho_max:
            break
    W_raw = _adj(w_est)
    W = W_raw.copy()
    W[np.abs(W) < w_threshold] = 0
    return W, W_raw


def notears_tiered(data, names, seed, n_background):
    """Background = the first n_background columns. Returns dict(adj, score) over all
    variables; the runner derives the ancestral matrix."""
    X = np.asarray(data, dtype=np.float64)
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    d = X.shape[1]
    forbid = np.zeros((d, d), dtype=bool)
    forbid[n_background:, :n_background] = True      # no foreground -> background edge
    W, W_raw = notears_linear(X, forbid=forbid)
    return {"adj": (W != 0).astype(float), "score": np.abs(W_raw)}


PY_METHODS = {"grnboost2": grnboost2, "regdiffusion": regdiffusion, "pidc": pidc,
              "notears_tiered": notears_tiered}
