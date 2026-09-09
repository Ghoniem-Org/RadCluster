#!/usr/bin/env python3
"""
make_size_effect_figures.py — dose dependence of each observable across the
Table 4 system-size ladder, with the EUROFER97 measurements overlaid.

One figure per observable, curves ordered D, C4, C3, C2, C1, C0 (increasing
number of equations after the exact arm).  d_111 is excluded by request.
C4/P=1 is excluded: it is a closure-ORDER column, not a system-size rung, and at
+345% it would set the scale of every panel.

DOSE COVERAGE.  These runs stop at 2 dpa; the EUROFER97 neutron data sit at
14.6-32 dpa (loops) and 2.7-32 dpa (cavities).  The overlay therefore shows
almost no overlap, and that is the honest picture rather than a defect of the
figure: Table 4 is a VERIFICATION at low dose, and nothing in it may be read as
agreement with experiment.  The band is drawn where the data actually are so a
reader cannot mistake the unconstrained region for a validated one.

Output: docs/Formulation/figs/size_effect_<observable>.pdf
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgba

MOD = Path(__file__).resolve().parent.parent
REPO = MOD.parent
sys.path.insert(0, str(REPO))
DB = REPO / "docs" / "Database"
OUT = REPO / "docs" / "Formulation" / "figs"

LABEL_PT, TICK_PT, LEGEND_PT = 22, 20, 15
plt.rcParams.update({
    "axes.labelsize": LABEL_PT, "axes.titlesize": LABEL_PT,
    "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
    "legend.fontsize": LEGEND_PT, "lines.linewidth": 2.6,
})

# The ladder, in the order requested: exact first, then coarsest -> finest.
LADDER = [
    ("D  (exact, 12006 eq)", "20260907_093333_T4_D_4k_*",   "#111111", "-"),
    ("C4  1/800   (84 eq)",  "*_T4_C4_nc_*",                "#d62728", "-"),
    ("C3  1/160  (108 eq)",  "*_T4_C3_nc_*",                "#ff7f0e", "-"),
    ("C2  1/40   (244 eq)",  "*_T4_C2_nc_*",                "#2ca02c", "-"),
    ("C1  1/10   (830 eq)",  "*_T4_C1_nc_*",                "#1f77b4", "-"),
    ("C0  1/2   (3214 eq)",  "*_T4_C0_nc_*",                "#9467bd", "-"),
]

# observable key -> (axis label, log?, experiment source, experiment column)
OBSERVABLES = [
    ("N_loops_111", r"$\frac{1}{2}\langle111\rangle$ density (m$^{-3}$)", True,  "loop111", "Density [m^-3]"),
    ("N_loops_100", r"$\langle100\rangle$ density (m$^{-3}$)",            True,  "loop100", "Density [m^-3]"),
    ("d_100_nm",    r"$\langle100\rangle$ mean diameter (nm)",            False, "loop100", "Diameter [nm]"),
    ("N_voids",     r"Cavity density (m$^{-3}$)",                         True,  "cavity",  "Density [m^-3]"),
    ("d_cavity_nm", r"Cavity mean diameter (nm)",                         False, "cavity",  "Diameter [nm]"),
]
EXP_COLOR = {"loop111": "#1f77b4", "loop100": "#7d3c98", "cavity": "#e8663c"}
EXP_MARK = {"loop111": "o", "loop100": "s", "cavity": "^"}


def euro(df):
    m = (df["Material"].astype(str).str.contains("urofer97", na=False)
         & df["Type of Irradiation"].astype(str).str.strip().str.lower().eq("neutron")
         & df["Irradiation Temp [C]"].between(300, 350))
    return df[m].copy()


def experiment():
    L = euro(pd.read_excel(DB / "InterstitialLoop.xlsx"))
    V = euro(pd.read_excel(DB / "Void.xlsx"))
    lt = L["Loop Type"].astype(str)
    return {"loop111": L[lt.str.contains("111", na=False)],
            "loop100": L[lt.str.contains("100", na=False)],
            "cavity":  V}


def with_diameters(d):
    """Add d_100_nm / d_cavity_nm as TIME SERIES.

    They are not in `results`: post_process stores mean_n_100 and mean_n_v, and
    the diameters are formed only inside the at_dose ladder.  Plotting the raw
    key therefore silently produced panels with the experimental band and no
    model curves at all -- caught by looking at the figure, not by the run
    succeeding.  Same expressions the ladder uses, so figure and table agree.
    """
    r = d["results"]
    idata = d["input_data"]
    Om = float(r.get("Omega", idata.derived["Omega"]))
    b111 = float(r.get("b_111", idata.derived["b_111"]))
    b100 = float(r.get("b_100", idata.energetics.get("b_100", 0.2867) * 1e-9))
    out = dict(r)
    n100 = np.maximum(np.asarray(r.get("mean_n_100", []), float), 0.0)
    nv = np.maximum(np.asarray(r.get("mean_n_v", []), float), 0.0)
    n111 = np.maximum(np.asarray(r.get("mean_n_111", []), float), 0.0)
    out["d_100_nm"] = 2 * np.sqrt(n100 * Om / (np.pi * b100)) * 1e9
    out["d_111_nm"] = 2 * np.sqrt(n111 * Om / (np.pi * b111)) * 1e9
    out["d_cavity_nm"] = 2 * (3 * nv * Om / (4 * np.pi)) ** (1.0 / 3.0) * 1e9
    return out


def trajectories():
    from RadCluster_2_1.py_utils import visualization as viz
    import run_ensemble as re_mod            # noqa: F401  (import parity)
    out = []
    for label, pat, colour, ls in LADDER:
        f = sorted(glob.glob(str(MOD / "output" / pat / "plots" / "plot_data.pkl")))
        if not f:
            print(f"  MISSING trajectory for {label} ({pat})")
            continue
        d = viz.load_plot_data(f[-1])
        r = with_diameters(d)
        out.append((label, colour, ls, np.asarray(r["dose"], float), r))
    return out


def one_figure(key, ylabel, logy, expkey, expcol, traj, EXP):
    fig, ax = plt.subplots(figsize=(11, 7.5))
    for label, colour, ls, dose, r in traj:
        y = np.asarray(r.get(key, []), float)
        if y.size != dose.size:
            continue
        m = dose > 0
        lw = 3.4 if label.startswith("D ") else 2.4
        z = 5 if label.startswith("D ") else 3
        ax.plot(dose[m], y[m], ls, color=colour, lw=lw, label=label, zorder=z)

    df = EXP[expkey].dropna(subset=["Dose [dpa]", expcol])
    if len(df):
        c = EXP_COLOR[expkey]
        x0, x1 = df["Dose [dpa]"].min(), df["Dose [dpa]"].max()
        y0, y1 = df[expcol].min(), df[expcol].max()
        if x1 > x0 and y1 > y0:
            ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                                   facecolor=to_rgba(c, 0.12),
                                   edgecolor=to_rgba(c, 0.85), lw=2.0,
                                   ls=(0, (6, 4)), zorder=1,
                                   label="EUROFER97 measured range"))
        ax.scatter(df["Dose [dpa]"], df[expcol], s=150, marker=EXP_MARK[expkey],
                   facecolors=c, edgecolors="white", linewidths=1.6, zorder=6,
                   label="EUROFER97 neutron, 300-350 $^\\circ$C")

    ax.set_xscale("log")
    # The solution is featureless below ~1e-6 dpa (clusters still nucleating out
    # of the C_floor initial condition) and plotting 15 decades squeezed every
    # curve and the experimental band into the right-hand third.  Cut the left
    # edge; the right edge keeps the measured band in view.
    ax.set_xlim(1e-6, 60)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("Dose (dpa)")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", frameon=True, framealpha=0.92)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / f"size_effect_{key}.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main():
    EXP = experiment()
    traj = trajectories()
    print(f"{len(traj)} trajectories loaded")
    for label, _, _, dose, _ in traj:
        print(f"  {label:26s} {len(dose):3d} points, to {dose.max():.3g} dpa")
    for key, ylabel, logy, ek, ec, in [(o[0], o[1], o[2], o[3], o[4]) for o in OBSERVABLES]:
        p = one_figure(key, ylabel, logy, ek, ec, traj, EXP)
        print(f"  wrote {p.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
