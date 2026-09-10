#!/usr/bin/env python3
"""make_M_figures.py — campaign M: binning convergence, drawn.

Two figures:

  M_convergence.pdf   the six observables against N_eq, one panel each, with
                      the exact arm as a horizontal line.  This is the
                      convergence statement: a rung is converged when its
                      panel has gone flat AND flat AT the reference line.
                      Flat alone is self-consistency, which is why the
                      reference is drawn even when it is the only thing on the
                      panel that is not a bin-moment solve.

  M_dose.pdf          the same observables against dose, all rungs overlaid.
                      A rung can agree at the scoring dose and disagree along
                      the way; the table cannot show that and this can.

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

OBS = [("N_loops_111", r"$\frac{1}{2}\langle111\rangle$ density (m$^{-3}$)", True),
       ("d_111_nm",    r"$\frac{1}{2}\langle111\rangle$ diameter (nm)",      False),
       ("N_loops_100", r"$\langle100\rangle$ density (m$^{-3}$)",            True),
       ("d_100_nm",    r"$\langle100\rangle$ diameter (nm)",                 False),
       ("N_voids",     r"Cavity density (m$^{-3}$)",                         True),
       ("d_cavity_nm", r"Cavity diameter (nm)",                              False)]


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

    fig, axes = plt.subplots(3, 2, figsize=(15, 15))
    for ax, (key, ylabel, logy) in zip(axes.ravel(), OBS):
        x = [r[0] for r in rows]
        y = [r[2][key] for r in rows]
        ax.plot(x, y, "o-", color="steelblue", ms=10, label="bin-moment rungs")
        if ref is not None:
            ax.axhline(ref[key], color="tomato", ls="--", lw=3,
                       label=f"exact ({rid})")
        else:
            # Say it on the figure, not only in the caption: a panel with no
            # reference line shows self-consistency, not accuracy, and a reader
            # scanning figures will not have the caption in view.
            ax.text(0.5, 0.06, "no exact arm — self-convergence only",
                    transform=ax.transAxes, ha="center", color="tomato",
                    fontsize=LEGEND_PT)
        ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel(r"$N_{\rm eq}$")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.25)
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
    for ax, (key, ylabel, _logy) in zip(axes.ravel(), OBS):
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
        ax.set_ylabel(ylabel.replace("(m$^{-3}$)", "dev (%)")
                            .replace("(nm)", "dev (%)"))
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
