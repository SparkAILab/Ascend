#!/usr/bin/env python3
"""
Build the headline ASCEND vs CBL figure from results CSV.

Usage:
    python make_figure.py results.csv outdir/

Output:
    outdir/ascend_vs_cbl.pdf
    outdir/ascend_vs_cbl.png
"""
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
from matplotlib.patches import Patch

# ---- style ----------------------------------------------------------------
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

C_ASCEND = "#0072B2"   # blue
C_CBL    = "#E69F00"   # orange
C_TIMEOUT= "#C44E52"   # red
GREY     = "#5A5A5A"

CSV_PATH = sys.argv[1] if len(sys.argv) > 1 else "/mnt/user-data/uploads/results_12.csv"
OUT_DIR  = sys.argv[2] if len(sys.argv) > 2 else "/mnt/user-data/outputs"

d = pd.read_csv(CSV_PATH)
d["timed_out"] = (d.status == "timeout")

# For metric panels, only ok rows. For runtime/speedup, include timeouts
# (they ARE part of CBL's time cost in any honest comparison).
ok = d[d.status == "ok"].copy()

def agg(df, group_cols, metric):
    """mean and SEM over seeds."""
    g = df.groupby(group_cols)[metric].agg(["mean", "std", "count"]).reset_index()
    g["sem"] = g["std"] / np.sqrt(g["count"])
    return g


def sweep_panel(ax, df, axis, metric, log_y=False, ribbon=True,
                show_timeouts=False, timeout_y_for_mean=False):
    """Plot mean ± SEM ribbon for each method along `axis`."""
    for method, colour in [("cbl", C_CBL), ("ascend", C_ASCEND)]:
        s = df[df.method == method]
        if len(s) == 0:
            continue
        a = agg(s, axis, metric).dropna(subset=["mean"]).sort_values(axis)
        if len(a) == 0:
            continue
        if ribbon and "sem" in a.columns:
            lo = a["mean"] - a["sem"]
            hi = a["mean"] + a["sem"]
            if log_y:
                lo = np.maximum(lo, 1e-9)
            ax.fill_between(a[axis], lo, hi, color=colour, alpha=0.18, lw=0)
        ax.plot(a[axis], a["mean"], "-", color=colour, lw=2.2, zorder=3)
        ax.plot(a[axis], a["mean"], "o", color=colour, ms=6.0,
                markeredgecolor="white", markeredgewidth=1.0, zorder=4)
    # Mark timeouts (CBL only in practice) with a red ✕ at the timeout time
    if show_timeouts:
        t = df[(df.timed_out) & (df.method == "cbl")]
        for xv, sub in t.groupby(axis):
            y = sub.time_sec.mean()
            ax.scatter([xv], [y], marker="x", color=C_TIMEOUT, s=85, lw=2.4,
                       zorder=6)
    if log_y:
        ax.set_yscale("log")


# Slices through default point
slice_dx = d[(d.n == 1024) & (d.d_z == 10)]
slice_dz = d[(d.n == 1024) & (d.d_x == 5)]


# ---- figure --------------------------------------------------------------
fig = plt.figure(figsize=(12.5, 7.4))
gs = fig.add_gridspec(2, 3, hspace=0.50, wspace=0.32,
                      left=0.06, right=0.985, top=0.88, bottom=0.08)
axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]

# =========================== ROW 1 : SPEED ================================
# (a) Runtime vs d_x
ax = axes[0][0]
sweep_panel(ax, slice_dx, "d_x", "time_sec", log_y=True, show_timeouts=True)
ax.set_xlabel(r"$d_x$  (foreground variables)")
ax.set_ylabel("Wall-clock time (s)")
ax.set_title(r"(a) Runtime vs $d_x$")
ax.set_xticks([5, 10, 15, 20])
ax.grid(axis="y", which="major", alpha=0.25, lw=0.5)
ax.text(0.02, 0.97, r"n=1024,  $d_z$=10",
        transform=ax.transAxes, va="top", ha="left",
        fontsize=8.8, color=GREY)

# (b) Runtime vs d_z — the most dramatic scaling story
ax = axes[0][1]
sweep_panel(ax, slice_dz, "d_z", "time_sec", log_y=True, show_timeouts=True)
ax.set_xlabel(r"$d_z$  (background variables)")
ax.set_ylabel("Wall-clock time (s)")
ax.set_title(r"(b) Runtime vs $d_z$")
ax.set_xticks([10, 20, 30, 40, 50])
ax.grid(axis="y", which="major", alpha=0.25, lw=0.5)
ax.text(0.02, 0.97, r"n=1024,  $d_x$=5",
        transform=ax.transAxes, va="top", ha="left",
        fontsize=8.8, color=GREY)
# Scaling exponent annotations, positioned out of the way.
ax.text(0.97, 0.62, r"CBL: $t \propto d_z^{+0.81}$",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9.5, color=C_CBL, fontweight="bold")
ax.text(0.97, 0.20, r"ASCEND: $t \propto d_z^{-0.45}$",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9.5, color=C_ASCEND, fontweight="bold")

# (c) Speedup factor bar chart
ax = axes[0][2]
def speedup_along(df, axis):
    # Treat timeouts as their observed (censored) time — this is a lower
    # bound on the true CBL cost, so the resulting speedup is conservative.
    g = df.groupby([axis, "method"])["time_sec"].mean().unstack("method")
    g = g.reset_index().sort_values(axis)
    g["speedup"] = g["cbl"] / g["ascend"]
    return g

sp_dx = speedup_along(slice_dx, "d_x")
sp_dz = speedup_along(slice_dz, "d_z")

# Track which cells contained CBL timeouts so we can hatch those bars
to_dx = set(d[(d.timed_out) & (d.n == 1024) & (d.d_z == 10)].d_x.unique())
to_dz = set(d[(d.timed_out) & (d.n == 1024) & (d.d_x == 5)].d_z.unique())

bars, labels, colours, hatches = [], [], [], []
for _, r in sp_dx.iterrows():
    bars.append(r["speedup"])
    labels.append(rf"$d_x$={int(r['d_x'])}")
    colours.append("#7FB3D5")
    hatches.append("///" if int(r["d_x"]) in to_dx else "")
for _, r in sp_dz.iterrows():
    bars.append(r["speedup"])
    labels.append(rf"$d_z$={int(r['d_z'])}")
    colours.append("#1f6f9c")
    hatches.append("///" if int(r["d_z"]) in to_dz else "")

x = np.arange(len(bars))
for xi, h, col, hat in zip(x, bars, colours, hatches):
    ax.bar(xi, h, color=col, edgecolor="white" if not hat else C_TIMEOUT,
           linewidth=0.8 if not hat else 1.2, hatch=hat)
for xi, h in zip(x, bars):
    ax.text(xi, h * 1.07, f"{h:,.0f}×".replace(",", ","),
            ha="center", va="bottom", fontsize=9.0)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9.0)
ax.set_yscale("log")
ax.set_ylabel("Speedup  (CBL time / ASCEND time)")
ax.set_title("(c) ASCEND speedup factor")
ax.axhline(1, color=GREY, lw=0.7, linestyle="--")
ax.set_ylim(0.7, max(bars) * 4)
ax.grid(axis="y", which="major", alpha=0.25, lw=0.5)
ax.legend(handles=[
    Patch(facecolor="#7FB3D5", label=r"varying $d_x$"),
    Patch(facecolor="#1f6f9c", label=r"varying $d_z$"),
    Patch(facecolor="white", edgecolor=C_TIMEOUT, hatch="///",
          label="≥1 CBL timeout (lower bound)"),
], loc="upper left", fontsize=8.6)

# =========================== ROW 2 : QUALITY ==============================
# (d) Computational work — CI tests performed
ax = axes[1][0]
sweep_panel(ax, slice_dz, "d_z", "n_ci_tests", log_y=True)
ax.set_xlabel(r"$d_z$")
ax.set_ylabel("Conditional independence\ntests performed")
ax.set_title("(d) Computational work")
ax.set_xticks([10, 20, 30, 40, 50])
ax.grid(axis="y", which="major", alpha=0.25, lw=0.5)
ax.text(0.02, 0.97, r"n=1024,  $d_x$=5",
        transform=ax.transAxes, va="top", ha="left",
        fontsize=8.8, color=GREY)

# (e) Decisiveness — coverage (1 - NA fraction)
ax = axes[1][1]
sweep_panel(ax, slice_dz, "d_z", "coverage")
ax.set_xlabel(r"$d_z$")
ax.set_ylabel("Coverage  (1 − NA fraction)")
ax.set_title("(e) Decisiveness")
ax.set_xticks([10, 20, 30, 40, 50])
ax.set_ylim(0.45, 1.03)
ax.grid(axis="y", which="major", alpha=0.25, lw=0.5)
ax.text(0.02, 0.05, r"n=1024,  $d_x$=5",
        transform=ax.transAxes, va="bottom", ha="left",
        fontsize=8.8, color=GREY)

# (f) Accuracy: F1 vs d_x
ax = axes[1][2]
sweep_panel(ax, slice_dx, "d_x", "f1")
ax.set_xlabel(r"$d_x$")
ax.set_ylabel("F1 score")
ax.set_title(r"(f) Accuracy vs $d_x$")
ax.set_xticks([5, 10, 15, 20])
ax.set_ylim(-0.02, 1.03)
ax.grid(axis="y", which="major", alpha=0.25, lw=0.5)
ax.text(0.02, 0.05, r"n=1024,  $d_z$=10",
        transform=ax.transAxes, va="bottom", ha="left",
        fontsize=8.8, color=GREY)

# ---- shared legend (top) -------------------------------------------------
handles = [
    mlines.Line2D([], [], color=C_ASCEND, marker="o", lw=2.2, ms=7,
                  markeredgecolor="white", markeredgewidth=1.0,
                  label="ASCEND (ours)"),
    mlines.Line2D([], [], color=C_CBL, marker="o", lw=2.2, ms=7,
                  markeredgecolor="white", markeredgewidth=1.0,
                  label="CBL (baseline)"),
    mlines.Line2D([], [], color=C_TIMEOUT, marker="x", lw=0, ms=9,
                  markeredgewidth=2.4, label="CBL exceeded 1 h timeout"),
]
fig.legend(handles=handles, loc="upper center",
           bbox_to_anchor=(0.5, 0.965), ncol=3, frameon=False,
           fontsize=10.5)

fig.suptitle("ASCEND vs CBL  —  scaling, computational cost, and discovery quality "
             "(5 seeds, ribbons = ±1 SEM)",
             y=1.00, fontsize=12.5, fontweight="bold")

for ext in ("pdf", "png"):
    out = f"{OUT_DIR}/ascend_vs_cbl.{ext}"
    plt.savefig(out)
    print("wrote", out)