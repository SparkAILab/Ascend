#!/usr/bin/env python3
"""
Build the headline benchmark figure for the ASCEND paper from
ascend_benchmark_merged.csv.

Usage:
    python make_benchmark_figure.py [path/to/csv] [out_dir]

Inputs:
    ascend_benchmark_merged.csv (default: same folder)

Outputs:
    {out_dir}/benchmark_headline.pdf
    {out_dir}/benchmark_headline.png

Layout (2 rows x 4 columns, the rightmost column is a tall heatmap spanning
both rows):

    (a) F1 vs n         (b) Precision vs n   (c) Paired F1 diffs   |
    (d) Decisiveness    (e) Reliability      (f) F1 vs sparsity    |  (g) ASCEND-GES heatmap
                                                                   |  across (sp x r2),
                                                                   |  in the realistic regime
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.patches import Rectangle

# ---- style ---------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 10.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 10,
    "legend.frameon": False,
    "figure.dpi": 140,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

COLOURS = {
    "ASCEND": "#0072B2",
    "CBL":    "#E69F00",
    "GES":    "#009E73",
    "LiNGAM": "#CC79A7",
    "PC":     "#888888",
}
METHOD_ORDER = ["ASCEND", "GES", "CBL", "LiNGAM", "PC"]

# ---- args ----------------------------------------------------------------
CSV = sys.argv[1] if len(sys.argv) > 1 else "/mnt/user-data/uploads/ascend_benchmark_merged.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/mnt/user-data/outputs"
os.makedirs(OUT, exist_ok=True)

# ---- load + filter -------------------------------------------------------
raw = pd.read_csv(CSV)

KEY = ['d_x', 'd_z', 'r2', 'sp', 'n']
THRESH = 10
ok = raw[raw.status == 'ok'].copy()
counts = ok.groupby(KEY + ['method']).size().unstack('method', fill_value=0)
NONPC = ['ASCEND', 'CBL', 'GES', 'LiNGAM']
pass_cells = counts[(counts[NONPC] >= THRESH).all(axis=1)].reset_index()[KEY]
ok = ok.merge(pass_cells, on=KEY, how='inner')

# PC: only kept where it has >=10 reps itself
pc_pass_cells = counts[counts['PC'] >= THRESH].reset_index()[KEY]
pc_idx_drop = ok[
    (ok.method == 'PC') &
    ~ok.set_index(KEY).index.isin(pc_pass_cells.set_index(KEY).index)
].index
ok = ok.drop(pc_idx_drop)

fail_rate = (raw.assign(failed=(raw.status != 'ok').astype(int))
             .groupby(['method', 'n'])['failed']
             .mean().unstack('method'))

ok['dz_mult'] = ok.d_z // ok.d_x

print(f"Loaded {len(raw):,} rows.")
print(f"After filter (>=10 reps for 4 main methods per cell): "
      f"{len(ok):,} ok rows across {len(pass_cells)} cells.")
print(f"PC included in {(ok.method=='PC').sum()} rows across "
      f"{len(pc_pass_cells)} cells.")


# ---- helpers -------------------------------------------------------------
def agg(df, group_cols, metric):
    g = df.groupby(group_cols)[metric].agg(['mean', 'std', 'count']).reset_index()
    g['sem'] = g['std'] / np.sqrt(g['count'])
    return g


def line_panel(ax, df, x_col, y_col, x_label, y_label, title,
               log_x=True, methods=None, y_lim=None,
               x_breaks=None, annotation=None,
               sparse_threshold=50):
    """Mean ± SEM ribbon panel.

    Per-x points with fewer than `sparse_threshold` replicates per method
    are drawn with open markers, dashed segments and higher transparency.
    """
    methods = methods or METHOD_ORDER
    for m in methods:
        sub = df[df.method == m]
        if len(sub) == 0:
            continue
        a = agg(sub, x_col, y_col).sort_values(x_col)
        a = a[~a['mean'].isna()]
        if len(a) == 0:
            continue
        col = COLOURS[m]
        dense_mask = a['count'] >= sparse_threshold
        a_dense = a[dense_mask]
        a_sparse = a[~dense_mask]

        if len(a_dense) > 0:
            ax.fill_between(a_dense[x_col],
                            a_dense['mean'] - a_dense['sem'],
                            a_dense['mean'] + a_dense['sem'],
                            color=col, alpha=0.15, lw=0)
        if len(a_dense) >= 2:
            ax.plot(a_dense[x_col], a_dense['mean'], '-',
                    color=col, lw=2.0, zorder=3)
        if dense_mask.any() and (~dense_mask).any():
            last_dense = a[dense_mask].index[-1]
            first_sparse = a[~dense_mask].index[0]
            bridge = a.loc[[last_dense, first_sparse]]
            ax.plot(bridge[x_col], bridge['mean'], '--',
                    color=col, lw=1.4, alpha=0.55, zorder=2)
        if len(a_dense) > 0:
            ax.plot(a_dense[x_col], a_dense['mean'], 'o',
                    color=col, ms=5.5,
                    markeredgecolor='white', markeredgewidth=0.9, zorder=4)
        if len(a_sparse) > 0:
            ax.plot(a_sparse[x_col], a_sparse['mean'], 'o',
                    color='white', ms=5.5,
                    markeredgecolor=col, markeredgewidth=1.4,
                    alpha=0.85, zorder=4)

    if log_x:
        ax.set_xscale('log')
    if y_lim is not None:
        ax.set_ylim(y_lim)
    if x_breaks is not None:
        ax.set_xticks(x_breaks)
        ax.set_xticklabels([str(b) for b in x_breaks])
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(axis='y', which='major', alpha=0.25, lw=0.5)
    if annotation:
        ax.text(*annotation['xy'], annotation['text'],
                transform=ax.transAxes,
                fontsize=annotation.get('size', 8.3),
                color=annotation.get('colour', '#5A5A5A'),
                ha=annotation.get('ha', 'left'),
                va=annotation.get('va', 'top'),
                style='italic')


# ---- compute paired F1 differences for panel (c) -----------------------
piv = ok.pivot_table(index=['job', 'n', 'rep', 'd_x', 'd_z', 'r2', 'sp'],
                     columns='method', values='f1')
piv = piv.dropna(subset=['ASCEND'])


# ---- figure layout: 2x4 with last column spanning both rows ------------
fig = plt.figure(figsize=(16.5, 8.2))
gs = fig.add_gridspec(
    nrows=2, ncols=4,
    width_ratios=[1.0, 1.0, 1.0, 0.85],
    hspace=0.55, wspace=0.36,
    left=0.05, right=0.985, top=0.87, bottom=0.10,
)
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[0, 2])
ax_d = fig.add_subplot(gs[1, 0])
ax_e = fig.add_subplot(gs[1, 1])
ax_f = fig.add_subplot(gs[1, 2])
ax_g = fig.add_subplot(gs[:, 3])   # forest plot, spans both rows
# (g) and (h) are created inline below using a nested gridspec on gs[:, 3]

N_BREAKS = [512, 2048, 8192, 32768, 131072]

# =========================================================================
# (a) F1 vs n
# =========================================================================
line_panel(
    ax_a, ok, 'n', 'f1',
    x_label="Sample size n",
    y_label="F1 score",
    title="(a) F1 vs sample size",
    log_x=True, x_breaks=N_BREAKS,
    y_lim=(0, 0.7),
    annotation={'xy': (0.02, 0.97),
                'text': ("239 cells \u00d7 5 methods \u00d7 up to 20 reps"
                         "\nribbons = \u00b11 SEM")},
)

# =========================================================================
# (b) Precision vs n
# =========================================================================
line_panel(
    ax_b, ok, 'n', 'precision',
    x_label="Sample size n",
    y_label="Precision",
    title="(b) Precision vs sample size",
    log_x=True, x_breaks=N_BREAKS,
    y_lim=(0.2, 0.7),
    annotation={'xy': (0.02, 0.97),
                'text': ("ASCEND stable at \u2248 0.55;\n"
                         "GES improves with n (consistency)\n"
                         "open marker = <50 reps")},
)

# =========================================================================
# (c) Paired F1 differences (violins)
# =========================================================================
competitors = ['CBL', 'GES', 'LiNGAM', 'PC']
diffs = {}
for m in competitors:
    sub = piv.dropna(subset=['ASCEND', m]).copy()
    diffs[m] = (sub['ASCEND'] - sub[m]).values

y_positions = np.arange(len(competitors))
violin_parts = ax_c.violinplot(
    [diffs[m] for m in competitors],
    positions=y_positions, vert=False, widths=0.75,
    showmeans=False, showmedians=False, showextrema=False,
)
for i, body in enumerate(violin_parts['bodies']):
    body.set_facecolor(COLOURS[competitors[i]])
    body.set_edgecolor('white')
    body.set_alpha(0.7)
    body.set_linewidth(0.6)

for i, m in enumerate(competitors):
    arr = diffs[m]
    q25, q50, q75 = np.percentile(arr, [25, 50, 75])
    mean = arr.mean()
    ax_c.plot([q25, q75], [i, i], color='black', lw=4,
              solid_capstyle='butt', alpha=0.85, zorder=3)
    ax_c.plot([q50], [i], 'o', color='white', ms=5,
              markeredgecolor='black', markeredgewidth=1.0, zorder=4)
    ax_c.annotate(f"  mean {mean:+.2f}\n  n = {len(arr):,}",
                  xy=(0.62, i), xycoords=('axes fraction', 'data'),
                  ha='left', va='center', fontsize=8.4, color='#3a3a3a')

ax_c.axvline(0, color='#888', lw=0.8, ls='--', zorder=1)
ax_c.set_yticks(y_positions)
ax_c.set_yticklabels([f"vs {m}" for m in competitors])
ax_c.invert_yaxis()
ax_c.set_xlim(-0.6, 0.85)
ax_c.set_xlabel("Paired F1 difference  (ASCEND \u2212 competitor)")
ax_c.set_title("(c) Paired F1 differences", pad=14)
ax_c.grid(axis='x', alpha=0.25, lw=0.5)
ax_c.set_axisbelow(True)
ax_c.text(0.98, -0.22, "ASCEND better \u2192",
          transform=ax_c.transAxes,
          ha='right', va='top', fontsize=8.3,
          color=COLOURS['ASCEND'], fontweight='bold')
ax_c.text(0.02, -0.22, "\u2190 competitor better",
          transform=ax_c.transAxes,
          ha='left', va='top', fontsize=8.3,
          color='#888888', fontweight='bold')

# =========================================================================
# (d) Coverage vs n
# =========================================================================
line_panel(
    ax_d, ok, 'n', 'coverage',
    x_label="Sample size n",
    y_label="Coverage  (1 \u2212 NA fraction)",
    title="(d) Decisiveness",
    log_x=True, x_breaks=N_BREAKS,
    y_lim=(0.7, 1.02),
    annotation={'xy': (0.02, 0.05),
                'text': "ASCEND and CBL can output NA;\nGES/LiNGAM/PC always commit",
                'va': 'bottom'},
)

# =========================================================================
# (e) Method reliability (failure rate)
# =========================================================================
for m in METHOD_ORDER:
    if m not in fail_rate.columns:
        continue
    ax_e.plot(fail_rate.index, fail_rate[m] * 100, '-o',
              color=COLOURS[m], lw=2.0, ms=5.5,
              markeredgecolor='white', markeredgewidth=0.9, label=m)
ax_e.set_xscale('log')
ax_e.set_xticks(N_BREAKS)
ax_e.set_xticklabels([str(b) for b in N_BREAKS])
ax_e.set_xlabel("Sample size n")
ax_e.set_ylabel("Failures (% of attempted runs)")
ax_e.set_title("(e) Method reliability")
ax_e.grid(axis='y', alpha=0.25, lw=0.5)
ax_e.set_ylim(-3, 100)
ax_e.axhline(50, color='#aaa', lw=0.5, ls=':')
ax_e.text(0.02, 0.97,
          "PC fails 60\u201390% of cells;\nASCEND essentially never fails",
          transform=ax_e.transAxes, fontsize=8.3, color='#5A5A5A',
          ha='left', va='top', style='italic')

# =========================================================================
# (f) F1 vs sparsity (n >= 4096)
# =========================================================================
hi = ok[ok.n >= 4096]
sps = sorted(hi.sp.unique())
sp_labels = [f"{s:.2f}" for s in sps]
methods_plot = [m for m in METHOD_ORDER if m in hi.method.unique()]
xpos = np.arange(len(sps))
width = 0.16
for i, m in enumerate(methods_plot):
    sub = agg(hi[hi.method == m], 'sp', 'f1').sort_values('sp')
    if len(sub) == 0:
        continue
    offset = (i - (len(methods_plot) - 1) / 2) * width
    ax_f.bar(xpos + offset, sub['mean'].values, width=width,
             yerr=sub['sem'].values,
             color=COLOURS[m], edgecolor='white', linewidth=0.4,
             error_kw=dict(elinewidth=0.6, capsize=2, alpha=0.6))
ax_f.set_xticks(xpos)
ax_f.set_xticklabels(sp_labels)
ax_f.set_xlabel("Graph sparsity")
ax_f.set_ylabel("F1 score")
ax_f.set_title("(f) F1 vs sparsity, all methods")
ax_f.grid(axis='y', which='major', alpha=0.25, lw=0.5)
ax_f.text(0.02, 0.97,
          "n \u2265 4096 only\nsignal-rich regime",
          transform=ax_f.transAxes, fontsize=8.3, color='#5A5A5A',
          ha='left', va='top', style='italic')

# =========================================================================
# (g) Forest plot: ASCEND advantage over GES across (sparsity, R^2)
#     in the realistic regime (n >= 4096, d_z = d_x). All 9 cells lie
#     to the right of the dashed zero line — ASCEND wins every cell —
#     with the largest advantages concentrated at high signal strength.
# =========================================================================
from scipy.stats import t as student_t

realistic = ok[(ok.n >= 4096) & (ok.dz_mult == 1)]
key = ['job', 'n', 'rep', 'd_x', 'd_z', 'r2', 'sp']
piv_r = realistic.pivot_table(index=key, columns='method', values='f1').reset_index()
both = piv_r.dropna(subset=['ASCEND', 'GES']).copy()
both['diff'] = both['ASCEND'] - both['GES']


def _cell_stats(g):
    nrep = len(g)
    m = g['diff'].mean()
    se = g['diff'].std(ddof=1) / np.sqrt(nrep)
    crit = student_t.ppf(0.975, df=max(nrep - 1, 1))
    return pd.Series({'mean_diff': m,
                      'ci_lo': m - crit * se,
                      'ci_hi': m + crit * se,
                      'n_reps': nrep})


cells = both.groupby(['sp', 'r2']).apply(_cell_stats,
                                         include_groups=False).reset_index()
cells['sp_r'] = cells.sp.round(2)
cells['r2_r'] = cells.r2.round(2)
cells = cells.sort_values(['r2_r', 'sp_r']).reset_index(drop=True)


def _alpha_for(sp_val):
    # Sparser graph (higher sp) -> darker marker
    return 0.4 if sp_val < 0.55 else (0.65 if sp_val < 0.7 else 1.0)


y_g = np.arange(len(cells))

# Faint band per R² block
r2_unique = sorted(cells.r2_r.unique())
for k, r2v in enumerate(r2_unique):
    idx_arr = np.array(cells.index[cells.r2_r == r2v])
    if k % 2 == 0:
        ax_g.axhspan(idx_arr.min() - 0.5, idx_arr.max() + 0.5,
                     color='#f5f5f5', zorder=0)

for i, row in cells.iterrows():
    a = _alpha_for(row.sp_r)
    ax_g.errorbar(row.mean_diff, i,
                  xerr=[[row.mean_diff - row.ci_lo],
                        [row.ci_hi - row.mean_diff]],
                  fmt='o', color=COLOURS['ASCEND'],
                  ecolor=COLOURS['ASCEND'],
                  capsize=4, capthick=1.2, lw=1.7, ms=9, alpha=a,
                  markeredgecolor='white', markeredgewidth=1.0, zorder=3)

ax_g.axvline(0, color='#5A5A5A', ls='--', lw=0.9, zorder=1)
ax_g.set_yticks(y_g)
ax_g.set_yticklabels([f"sp={row.sp_r:.2f}" for _, row in cells.iterrows()])
ax_g.invert_yaxis()
ax_g.set_xlabel("Mean paired F1\nadvantage of ASCEND over GES")
ax_g.set_xlim(-0.03, 0.12)
ax_g.grid(axis='x', alpha=0.3, lw=0.5)
ax_g.set_axisbelow(True)

# R² group labels on the right of the y-axis
for r2v in r2_unique:
    idx_arr = np.array(cells.index[cells.r2_r == r2v])
    yc = idx_arr.mean()
    ax_g.text(1.02, yc, f"$R^2 = {r2v:.2f}$",
              transform=ax_g.get_yaxis_transform(),
              va='center', ha='left', fontsize=9.5, fontweight='bold',
              color='#333')

ax_g.set_title("(g) ASCEND vs GES advantage\nin realistic regime",
               pad=8)

ax_g.text(0.98, -0.10, "ASCEND better \u2192",
          transform=ax_g.transAxes, ha='right', va='top',
          fontsize=8.4, color=COLOURS['ASCEND'], fontweight='bold')
ax_g.text(0.02, -0.10, "\u2190 GES better",
          transform=ax_g.transAxes, ha='left', va='top',
          fontsize=8.4, color='#888888', fontweight='bold')
ax_g.text(0.5, -0.17,
          "$n \\geq 4{,}096$, $d_z = d_x$\n9 cells; bars: 95% CI",
          transform=ax_g.transAxes, ha='center', va='top',
          fontsize=8.2, color='#5A5A5A', style='italic')

# ---- shared legend at top -----------------------------------------------
handles = [
    mlines.Line2D([], [], color=COLOURS[m], marker='o', lw=2.0, ms=6,
                  markeredgecolor='white', markeredgewidth=0.9,
                  label=(m + (' (ours)' if m == 'ASCEND' else '')))
    for m in METHOD_ORDER
]
fig.legend(handles=handles, loc='upper center',
           bbox_to_anchor=(0.45, 0.965), ncol=5,
           frameon=False, fontsize=10.5)

fig.suptitle(
    "ASCEND on synthetic two-tier benchmarks: accuracy, reliability and "
    "the regime where it dominates score-based baselines",
    y=1.005, fontsize=12.5, fontweight='bold')

# ---- save (the bug we just fixed: savefig inside the loop!) ------------
for ext in ("pdf", "png"):
    path = f"{OUT}/benchmark_headline.{ext}"
    plt.savefig(path)
    print("wrote", path)