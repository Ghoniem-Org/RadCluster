#!/usr/bin/env python3
"""make_M_figures.py — campaign M: binning convergence, drawn.

Two figures:

  M_convergence.pdf   the <100> and cavity observables against N_eq, one panel
                      each (the two <111> panels are omitted), with
                      the exact arm as a horizontal line.  This is the
                      convergence statement: a rung is converged when its
                      panel has gone flat AND flat AT the reference line.
                      Flat alone is self-consistency, which is why the
                      reference is drawn even when it is the only thing on the
                      panel that is not a bin-moment solve.

  M_2d_sweep.pdf      deviation against i_discrete, one curve per I_bin.  Kept
                      as a diagnostic; it is no longer carried in the report,
                      where the M2D tables make the same point.

Reads the campaign board, so it draws whatever has actually completed and says
so, rather than failing when a rung is missing.
"""
from __future__ import annotations

import glob
import json
import sys
import types
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

HERE = Path(__file__).resolve().parent
DT = HERE.parent
MOD = DT.parent
REPO = MOD.parent
for p in (str(REPO), str(DT), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import runs as manifest_mod                      # noqa: E402

SYNC = HERE / ".sync" / "claims"
OUT = REPO / "docs" / "Formulation" / "figs"
SCORE_DOSE = "20"

LABEL_PT, TICK_PT, LEGEND_PT = 22, 20, 15
plt.rcParams.update({
    "axes.labelsize": LABEL_PT, "axes.titlesize": LABEL_PT,
    "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT,
    "legend.fontsize": LEGEND_PT, "lines.linewidth": 3,
})

# (key, axis name, unit) -- the unit is kept separate from the name so that a
# common power of ten can be factored out of the tick labels and carried in the
# axis label instead.  Densities vary by a few parts in a thousand across these
# rungs, and a log axis over that range labels its MINOR ticks, which produced
# ordinates reading "3.74 x 10^21", "3.72 x 10^21" and so on down the panel.
OBS = [("N_loops_111", r"$\frac{1}{2}\langle111\rangle$ density", "m$^{-3}$"),
       ("d_111_nm",    r"$\frac{1}{2}\langle111\rangle$ diameter", "nm"),
       ("N_loops_100", r"$\langle100\rangle$ density",             "m$^{-3}$"),
       ("d_100_nm",    r"$\langle100\rangle$ diameter",            "nm"),
       ("N_voids",     r"Cavity density",                           "m$^{-3}$"),
       ("d_cavity_nm", r"Cavity diameter",                          "nm")]

# Panels of the convergence figure.  The two <111> panels are omitted: the
# mean-size behaviour they showed is set out quantitatively in the text, and
# the density panel adds nothing the <100> and cavity panels do not.
CONV_KEYS = ("N_loops_100", "d_100_nm", "N_voids", "d_cavity_nm")


def _factor(values):
    """Common power of ten for a set of values, and its label fragment.

    Returns (divisor, prefix) so that the panel plots y/divisor and the axis
    label reads e.g. "density ($10^{21}$ m$^{-3}$)".  Values within two decades
    of unity are left alone.
    """
    v = np.abs(np.asarray(list(values), float))
    v = v[np.isfinite(v) & (v > 0)]
    if v.size == 0:
        return 1.0, ""
    e = int(np.floor(np.log10(np.median(v))))
    if -2 <= e <= 2:
        return 1.0, ""
    return 10.0 ** e, rf"$10^{{{e}}}$ "


def board() -> dict:
    out = {}
    for f in sorted(SYNC.glob("*.json")):
        try:
            j = json.loads(f.read_text(encoding="utf-8"))
            out[j.get("run_id", f.stem)] = j
        except Exception:
            pass
    return out


def scored(rec):
    if not rec or rec.get("status") != "done":
        return None
    v = (rec.get("scored") or {}).get(SCORE_DOSE)
    return None if (not v or v.get("missing") or v.get("off_grid")) else v


def ref_record(B):
    """The exact arm, whichever M1 rung is serving as it."""
    rid = getattr(manifest_mod, "M_REF_RUN", None)
    if rid and scored(B.get(rid)):
        return rid, scored(B[rid])
    return None, None


def convergence(B):
    rows = []
    # The REFERENCED ladder (M2R, I = V = 1000).  M2 at 10000 has no exact arm
    # -- the discrete solve there diverges -- so plotting it here would put a
    # convergence figure in front of a reader with nothing to converge TO.
    for nb in manifest_mod.M2R_BINS:
        rec = B.get(f"M2R_B{nb}")
        v = scored(rec)
        if v:
            rows.append((rec.get("N_eq", np.nan), nb, v))
    if not rows:
        print("  no completed M2 rungs; nothing to draw")
        return None
    rid, ref = ref_record(B)

    panels = [o for o in OBS if o[0] in CONV_KEYS]
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    x = [r[0] for r in rows]
    for ax, (key, name, unit) in zip(axes.ravel(), panels):
        y = [r[2][key] for r in rows]
        series = list(y) + ([ref[key]] if ref is not None else [])
        div, pref = _factor(series)
        ax.plot(x, [v / div for v in y], "o-", color="steelblue", ms=10,
                label="bin-moment rungs")
        if ref is not None:
            ax.axhline(ref[key] / div, color="tomato", ls="--", lw=3,
                       label=f"exact ({rid})")
        else:
            # Say it on the figure, not only in the caption: a panel with no
            # reference line shows self-consistency, not accuracy, and a reader
            # scanning figures will not have the caption in view.
            ax.text(0.5, 0.06, "no exact arm — self-convergence only",
                    transform=ax.transAxes, ha="center", color="tomato",
                    fontsize=LEGEND_PT)
        # Both axes linear.  N_eq spans 120 to 250 and the observables a few
        # parts in a thousand; a logarithmic axis over less than one decade
        # labels its minor ticks, and those labels ran into one another.
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5, integer=True))
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
        ax.set_xlabel(r"$N_{\rm eq}$")
        ax.set_ylabel(f"{name} ({pref}{unit})")
        ax.grid(True, which="major", alpha=0.25)
    axes.ravel()[0].legend(loc="best", frameon=True, framealpha=0.92)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "M_convergence.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def two_d(B):
    """The 2-D sweep: deviation against i_discrete, one curve per I_bin.

    This is the figure the M2R ladder could not draw.  M2R swept I_bin at fixed
    i_discrete and found d_111 flat at -5 to -7%; the convergence axis for that
    observable is i_discrete, so the ladder was moving along a contour of the
    error surface rather than down it.  Plotted this way the curves collapse:
    the spread BETWEEN I_bin values at one i_discrete is small, while the fall
    ALONG i_discrete is the whole effect.
    """
    rid, ref = ref_record(B)
    if ref is None:
        print("  no exact arm; skipping the 2-D figure")
        return None
    fig, axes = plt.subplots(3, 2, figsize=(15, 15))
    cmap = plt.get_cmap("viridis")
    drawn = False
    for ax, (key, name, _unit) in zip(axes.ravel(), OBS):
        for c, nb in enumerate(manifest_mod.M2D_BINS):
            xs, ys = [], []
            for i_d in manifest_mod.M2D_IDISC:
                v = scored(B.get(f"M2D_I{i_d}_B{nb}"))
                if v is None or not ref[key]:
                    continue
                xs.append(i_d)
                ys.append((v[key] - ref[key]) / ref[key] * 100.0)
            if xs:
                drawn = True
                ax.plot(xs, ys, "o-", ms=9,
                        color=cmap(c / max(1, len(manifest_mod.M2D_BINS) - 1)),
                        label=f"$I_{{\\rm bin}}$ = {nb}")
        ax.axhline(0.0, color="tomato", ls="--", lw=3)
        ax.set_xscale("log")
        ax.set_xlabel(r"$i_{\rm discrete}$")
        ax.set_ylabel(f"{name} dev (%)")
        ax.grid(True, which="both", alpha=0.25)
    if not drawn:
        plt.close(fig)
        print("  no completed 2-D cells yet")
        return None
    axes.ravel()[0].legend(loc="best", frameon=True, framealpha=0.92)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / "M_2d_sweep.pdf"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main():
    B = board()
    for f in (convergence, two_d):
        p = f(B)
        if p:
            print(f"  wrote {p.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
