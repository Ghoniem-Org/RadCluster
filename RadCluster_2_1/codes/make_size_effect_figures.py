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
#
# The counts are the TOTAL integrated system -- main block plus the appended
# <100> loop block, which cpp_bridge splits off after the solve so the solver
# only ever reports the main one.  Labelling D with the total (12006) while the
# C rungs carried their main block (830, ...) made the same run read as 8006 in
# Table 4 and 12006 in the legend.  One quantity everywhere; Table 4 now prints
# both columns.
LADDER = [
    ("D  (exact, 12006 eq) - to 2 dpa", "20260907_093333_T4_D_4k_*",   "#111111", "-"),
    # C2-C4 are the 40 dpa reruns: at 2 dpa no model curve reached the dose
    # where the EUROFER97 data sit (14.6-32 dpa for loops), so the overlay had
    # no shared x-range at all.  Same grids as their 2 dpa twins.
    ("C4  1/800   (129 eq)", "*_T4_C4_40_*",                "#d62728", "-"),
    ("C3  1/160   (161 eq)", "*_T4_C3_40_*",                "#ff7f0e", "-"),
    ("C2  1/40    (364 eq)", "*_T4_C2_40_*",                "#2ca02c", "-"),
    ("C1  1/10   (1242 eq) - to 2 dpa", "*_T4_C1_nc_*",                "#1f77b4", "-"),
    ("C0  1/2    (4826 eq)", "*_T4_C0_nc_*",                "#9467bd", "-"),
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

# Lower bound of the ordinate, per observable.  On the log panels the solution
# rises out of the C_floor initial condition through ten empty decades that
# carry no information and compress every curve into the top of the frame;
# these floors put the axis where the populations actually are.
YMIN = {
    "N_loops_111": 1e19,
    "N_loops_100": 1e12,
    "N_voids":     1e19,
}

# One EUROFER97 <100> measurement reports a mean loop diameter of 48 nm at
# 17.4 dpa, 350 C, against 2.8-8.0 nm for the other nine rows of the same
# population.  Carried into the diameter panel it sets the ordinate on its own
# and compresses the model curves and the remaining data into the bottom tenth
# of the frame, so it is excluded from the diameter comparison.  Its density is
# unremarkable (3.7e21 m^-3, interior to the other rows) and is retained, so
# the exclusion is by quantity rather than a discarded measurement.
EXP_EXCLUDE = {
    "d_100_nm": lambda df: df["Diameter [nm]"] < 40.0,
}

# --- provenance ---------------------------------------------------------
# The "Paper" column of both workbooks is merged-cell style: the citation is
# written once at the head of a block and left blank on the rows beneath it, so
# it has to be forward-filled before a row can be attributed.
#
# Forward-filling the loop sheet puts one NEUTRON row (BOR-60, 15 dpa, <100>,
# 6.2 nm) under Boulanger & Serruys, which is an ION-irradiation study and
# therefore cannot be its source.  That row pairs with the 1/2<111> row directly
# below it -- same 6.2 nm, same facility, same dose -- which carries the Chauhan
# label, so the block boundary is off by one and the row belongs to Chauhan.
# Corrected explicitly rather than silently: an ion paper cited as the source of
# a neutron measurement would be a real error in the report.
PAPER_TO_KEY = {
    "Dethloff et al. - 2016":   "dethloff2016",
    "Dethloff et al. - 2018":   "dethloff2018",
    "Klimenkov et al. - 2011":  "klimenkov2011",
    "Klimenkov et al. - 2020":  "klimenkov2020",
    "Weiß et al. - 2012":       "weiss2012",
    "Chauhan et al. - 2021":    "chauhan2021",
    "Coppola et al. - 2019":    "coppola2019",
    "Boulanger and Serruys - 2009": "chauhan2021",   # see note above
}
# Order in which keys are listed in a citation, so captions read consistently.
KEY_ORDER = ["dethloff2016", "dethloff2018", "klimenkov2011", "klimenkov2020",
             "weiss2012", "chauhan2021", "coppola2019"]
SOURCES_TEX = OUT / "sources.tex"

# Plotted dose range.  The solution is featureless below 1e-6 dpa (clusters are
# still nucleating out of the C_floor initial condition); the upper edge keeps
# the measured band in view.
XLIM = (1e-6, 60.0)

# Per-observable overrides.  The <100> diameter spends its first three decades
# in the nucleation transient, which on a linear ordinate costs most of the
# frame and says nothing about the closure; the panel starts where the
# populations are established instead.
XLIM_OVERRIDE = {
    "d_100_nm": (1e-3, 60.0),
}


def euro(df):
    m = (df["Material"].astype(str).str.contains("urofer97", na=False)
         & df["Type of Irradiation"].astype(str).str.strip().str.lower().eq("neutron")
         & df["Irradiation Temp [C]"].between(300, 350))
    return df[m].copy()


def _read(name):
    df = pd.read_excel(DB / name)
    df["Paper"] = df["Paper"].ffill()
    return euro(df)


def experiment():
    L = _read("InterstitialLoop.xlsx")
    V = _read("Void.xlsx")
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


def _linear_extent(key, traj, df, expcol):
    """Range of the model curves and measurements *inside the plotted window*.

    Restricted to the panel's own x-limits: the transient below the left edge
    carries diameters well above anything on screen, and including it set the
    ordinate from data the reader never sees.
    """
    lo_x, hi_x = XLIM_OVERRIDE.get(key, XLIM)
    vals = []
    for _, _, _, dose, r in traj:
        y = np.asarray(r.get(key, []), float)
        if y.size == dose.size:
            y = y[(dose >= lo_x) & (dose <= hi_x) & np.isfinite(y)]
            if y.size:
                vals.append((y.min(), y.max()))
    if len(df):
        vals.append((float(df[expcol].min()), float(df[expcol].max())))
    if not vals:
        return None, None
    return min(v[0] for v in vals), max(v[1] for v in vals)


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
    keep = EXP_EXCLUDE.get(key)
    if keep is not None and len(df):
        df = df[keep(df)]
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
    ax.set_xlim(*XLIM_OVERRIDE.get(key, XLIM))
    if logy:
        ax.set_yscale("log")
        if key in YMIN:
            ax.set_ylim(bottom=YMIN[key])
    else:
        # Linear panels: frame the model curves and the retained measurements
        # together, with a 5 percent margin, rather than letting matplotlib
        # pad around whatever the widest single series happens to be.
        lo, hi = _linear_extent(key, traj, df, expcol)
        if lo is not None:
            pad = 0.05 * (hi - lo) if hi > lo else 0.05 * max(hi, 1.0)
            ax.set_ylim(max(0.0, lo - pad), hi + pad)
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


def _keys(df):
    """Citation keys for the rows actually plotted, in KEY_ORDER."""
    seen = set()
    for paper in df["Paper"].dropna().astype(str):
        head = " - ".join(paper.split(" - ")[:2])
        key = PAPER_TO_KEY.get(head)
        if key is None:
            print(f"  UNMAPPED source: {head}")
        else:
            seen.add(key)
    return [k for k in KEY_ORDER if k in seen]


# LaTeX control sequences may contain letters only, so the digits in the
# observable keys are spelled out.  Left as digits, \newcommand{\srcNLoops111}
# defines \srcNLoops, drops "111" into the preamble, and then typesets the
# \cite -- which is how this first went wrong.
_DIGIT = {"0": "Zero", "1": "One", "2": "Two", "3": "Three", "4": "Four",
          "5": "Five", "6": "Six", "7": "Seven", "8": "Eight", "9": "Nine"}


def _macro_name(key):
    camel = "".join(w.capitalize() for w in key.split("_"))
    return "".join(_DIGIT.get(c, c) for c in camel)


def write_sources(panel_keys):
    """Emit \\srcNloops111 etc. so the captions cite what was actually drawn.

    Hand-copied citations drift the moment a row is added to the database; this
    keeps the caption and the plotted points derived from one source of truth.
    """
    lines = ["% Generated by codes/make_size_effect_figures.py -- do not edit.",
             "% One macro per panel, listing the sources of the points drawn on it."]
    for key, keys in panel_keys.items():
        macro = "src" + _macro_name(key)
        cites = ",".join(keys) if keys else ""
        lines.append(f"\\newcommand{{\\{macro}}}{{\\cite{{{cites}}}}}")
    SOURCES_TEX.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return SOURCES_TEX


def main():
    EXP = experiment()
    traj = trajectories()
    print(f"{len(traj)} trajectories loaded")
    for label, _, _, dose, _ in traj:
        print(f"  {label:26s} {len(dose):3d} points, to {dose.max():.3g} dpa")
    panel_keys = {}
    for key, ylabel, logy, ek, ec, in [(o[0], o[1], o[2], o[3], o[4]) for o in OBSERVABLES]:
        p = one_figure(key, ylabel, logy, ek, ec, traj, EXP)
        df = EXP[ek].dropna(subset=["Dose [dpa]", ec])
        keep = EXP_EXCLUDE.get(key)
        if keep is not None and len(df):
            df = df[keep(df)]
        panel_keys[key] = _keys(df)
        print(f"  wrote {p.relative_to(REPO)}  sources: {', '.join(panel_keys[key])}")
    sp = write_sources(panel_keys)
    print(f"  wrote {sp.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
