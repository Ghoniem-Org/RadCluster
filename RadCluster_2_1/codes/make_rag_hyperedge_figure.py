#!/usr/bin/env python
"""make_rag_hyperedge_figure.py — the binary half of the EUROFER-97 RAG.

Companion to ``make_rag_figure.py``, which draws the unary skeleton and
omits COALESCENCE and ANNIHILATION because they are directed *hyperedges*
(two tails, one head) and O(n_max**2) in count -- 360 arcs at n_max=8,
6564 at n_max=30.  A node-link drawing of those is a solid block.

A hyperedge family is a *reaction table*, so it is drawn as one: the
(n, m) plane, with the cell carrying the head the pair maps to.  Every
cell IS one hyperedge, so nothing is sampled or elided.

The pairs are reconstructed from the materialized graph by grouping arcs
on ``hyper_id`` -- each id carries exactly the two tail->head arcs of one
hyperedge -- so this figure, like its companion, is generated from the
declaration rather than from a transcription of the arithmetic.

ENCODING.  Coalescence carries a magnitude (the product size n+m), so it
gets the sequential blue ramp.  Annihilation carries a polarity -- which
sub-space survives the encounter -- so it gets the diverging blue<->red
pair with a neutral gray midpoint at n == m, where the two annihilate
exactly and the content leaves through SINK.  Blue means vacancy in both
panels and in the companion figure.

Usage::

    cd RadCluster_2_1
    python codes/make_rag_hyperedge_figure.py --n-max 8 \
        --out ../docs/Formulation/rag/eurofer_rag_hyperedges.png
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from py_utils.input_data import InputData                      # noqa: E402
from py_utils.reaction_rates import ReactionRates              # noqa: E402
from py_utils.materials import build_eurofer_rag               # noqa: E402
from py_utils.core.to_networkx import to_networkx              # noqa: E402

# ── house style (py_utils/visualization.py) ──────────────────────────────
_LABEL_FONTSIZE = 22
_TICK_FONTSIZE = 20
_PLOT_FONTSIZE = 16
_LW = 3.0
_LW_THIN = 2.4
_DPI = 150

plt.rcParams.update({
    "axes.labelsize": _LABEL_FONTSIZE,
    "axes.titlesize": _LABEL_FONTSIZE,
    "xtick.labelsize": _TICK_FONTSIZE,
    "ytick.labelsize": _TICK_FONTSIZE,
    "legend.fontsize": _PLOT_FONTSIZE,
    "legend.title_fontsize": _PLOT_FONTSIZE,
    "figure.dpi": _DPI,
})

_INK = "#1a1a1a"
_GREY_TEXT = "#6a6a6a"
_OUTSIDE = "#e8e7e3"     # pair declared but head falls outside the domain

# Documented blue ramp, 100 -> 700 (sequential: magnitude).
_SEQ = LinearSegmentedColormap.from_list("radcluster_seq", [
    "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b",
])
# Documented diverging pair: blue <-> red, neutral gray midpoint.
_DIV = LinearSegmentedColormap.from_list("radcluster_div", [
    "#8f2020", "#d03b3b", "#e06a6a", "#eda3a3", "#f0efec",
    "#9ec5f4", "#3987e5", "#184f95", "#0d366b",
])

_POP_LABEL = {
    "bulk-111": r"$\frac{1}{2}\langle 111\rangle$",
    "bulk-100": r"$\langle 100\rangle$",
    "bulk": "vacancy",
}


def _key(node: str):
    if node in ("SOURCE", "SINK"):
        return node, None
    _, pop, n = node.split(":")
    return pop, int(n)


def collect_hyperedges(G):
    """Group materialized arcs by ``hyper_id`` -> (tails, head, label).

    ``to_networkx`` emits each binary reaction as its two tail->head arcs
    sharing one ``hyper_id``, which is exactly enough to recover the
    (n, m) pair and its product without re-deriving the arithmetic here.
    """
    groups = defaultdict(list)
    for u, v, d in G.edges(data=True):
        if d.get("hyperedge") and d.get("hyper_id") is not None:
            groups[d["hyper_id"]].append((u, v, d))
    out = []
    for hid, arcs in groups.items():
        head = arcs[0][1]
        tails = [a[0] for a in arcs]
        out.append((tails, head, arcs[0][2]["edge_class"], arcs[0][2]["label"]))
    return out


def _cell_text_colour(rgba) -> str:
    r, g, b = rgba[:3]
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "white" if lum < 0.55 else _INK


def _draw_matrix(ax, grid, mask_outside, n_max, cmap, norm, fmt):
    for i in range(n_max):
        for j in range(n_max):
            n, m = j + 1, i + 1
            if mask_outside[i, j]:
                ax.add_patch(Rectangle((n - 0.5, m - 0.5), 1, 1,
                                       facecolor=_OUTSIDE, edgecolor="white",
                                       lw=1.5, hatch="///", zorder=1))
                continue
            if np.isnan(grid[i, j]):
                ax.add_patch(Rectangle((n - 0.5, m - 0.5), 1, 1,
                                       facecolor="white", edgecolor="white",
                                       lw=1.5, zorder=1))
                continue
            rgba = cmap(norm(grid[i, j]))
            ax.add_patch(Rectangle((n - 0.5, m - 0.5), 1, 1, facecolor=rgba,
                                   edgecolor="white", lw=1.5, zorder=1))
            ax.text(n, m, fmt(grid[i, j]), ha="center", va="center",
                    fontsize=_PLOT_FONTSIZE, color=_cell_text_colour(rgba),
                    zorder=2)
    ax.set_xlim(0.5, n_max + 0.5)
    ax.set_ylim(0.5, n_max + 0.5)
    ax.set_xticks(range(1, n_max + 1))
    ax.set_yticks(range(1, n_max + 1))
    ax.set_aspect("equal")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#b4b4b4")


def _schematic(ax):
    """The two-tails/one-head convention, drawn once."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect("equal")
    ax.axis("off")

    def reaction(y, t1, t2, head, colour, caption):
        for (cx, cy), lab in (((1.5, y + 1.1), t1), ((1.5, y - 1.1), t2)):
            ax.add_patch(Circle((cx, cy), 0.72, facecolor="white",
                                edgecolor="#3f3f3f", lw=_LW_THIN, zorder=3))
            ax.text(cx, cy, lab, ha="center", va="center", fontsize=_PLOT_FONTSIZE,
                    color=_INK, zorder=4)
        ax.add_patch(Rectangle((4.55, y - 0.45), 0.9, 0.9, facecolor=colour,
                               edgecolor="none", zorder=3))
        for cy in (y + 1.1, y - 1.1):
            ax.add_patch(FancyArrowPatch((2.25, cy), (4.5, y),
                                         arrowstyle="-|>", mutation_scale=18,
                                         lw=_LW_THIN, color=colour, zorder=2))
        ax.add_patch(FancyArrowPatch((5.5, y), (7.55, y), arrowstyle="-|>",
                                     mutation_scale=18, lw=_LW, color=colour,
                                     zorder=2))
        ax.add_patch(Circle((8.3, y), 0.72, facecolor="white",
                            edgecolor="#3f3f3f", lw=_LW_THIN, zorder=3))
        ax.text(8.3, y, head, ha="center", va="center", fontsize=_PLOT_FONTSIZE,
                color=_INK, zorder=4)
        ax.text(5.0, y - 1.85, caption, ha="center", va="top",
                fontsize=_PLOT_FONTSIZE, color=_GREY_TEXT)

    reaction(7.0, "$n$", "$m$", "$n{+}m$", "#256abf", "coalescence")
    reaction(2.4, "$n$", "$m$", "$n{-}m$", "#d03b3b", "annihilation")


def build_figure(n_max: int, out: Path) -> None:
    inp = InputData(I=200, V=200, physics_option="full_CD_fission")
    rag, _ = build_eurofer_rag(inp, ReactionRates(inp))
    G = to_networkx(rag, n_max=n_max, include_binary=True)
    hyper = collect_hyperedges(G)

    coal = [h for h in hyper if h[2] == "coalescence"]
    annih = [h for h in hyper if h[2] == "annihilation"]

    # ── coalescence: head size n+m, NaN where the pair is not declared ───
    c_grid = np.full((n_max, n_max), np.nan)
    for tails, head, _, _ in coal:
        sizes = sorted(_key(t)[1] for t in tails)
        n, m = sizes[0], sizes[-1]
        hp, hn = _key(head)
        if n <= n_max and m <= n_max:
            c_grid[m - 1, n - 1] = hn
            c_grid[n - 1, m - 1] = hn
    # a pair whose product leaves the domain is skipped by the exporter
    c_outside = np.zeros((n_max, n_max), bool)
    for i in range(n_max):
        for j in range(n_max):
            if np.isnan(c_grid[i, j]) and (i + 1) + (j + 1) > n_max:
                c_outside[i, j] = True

    # ── annihilation: signed survivor, + vacancy / - SIA / 0 mutual ─────
    a_grid = np.full((n_max, n_max), np.nan)
    for tails, head, _, _ in annih:
        vac = [t for t in tails if _key(t)[0] == "bulk"]
        sia = [t for t in tails if _key(t)[0] in ("bulk-111", "bulk-100")]
        if not vac or not sia:
            continue
        nv, ns = _key(vac[0])[1], _key(sia[0])[1]
        if nv <= n_max and ns <= n_max:
            a_grid[ns - 1, nv - 1] = nv - ns
    a_outside = np.zeros((n_max, n_max), bool)

    fig = plt.figure(figsize=(26.0, 10.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.0, 1.0], wspace=0.28)
    ax0, ax1, ax2 = (fig.add_subplot(gs[0, k]) for k in range(3))

    _schematic(ax0)

    lo, hi = np.nanmin(c_grid), np.nanmax(c_grid)
    n_seq = Normalize(vmin=lo, vmax=hi)
    _draw_matrix(ax1, c_grid, c_outside, n_max, _SEQ, n_seq,
                 lambda v: f"{int(v)}")
    ax1.set_xlabel("$n$")
    ax1.set_ylabel("$m$")
    fig.colorbar(plt.cm.ScalarMappable(norm=n_seq, cmap=_SEQ), ax=ax1,
                 fraction=0.046, pad=0.03,
                 label="product size $n{+}m$")

    amax = np.nanmax(np.abs(a_grid))
    n_div = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)
    _draw_matrix(ax2, a_grid, a_outside, n_max, _DIV, n_div,
                 lambda v: f"{int(abs(v))}" if v else "0")
    ax2.set_xlabel("vacancy cluster size $n$")
    ax2.set_ylabel("SIA cluster size $m$")
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=n_div, cmap=_DIV), ax=ax2,
                      fraction=0.046, pad=0.03,
                      label="surviving cluster")
    cb.set_ticks([-amax, 0, amax])
    cb.set_ticklabels(["SIA", "mutual\n(SINK)", "vacancy"])

    fams_c = sorted({h[3] for h in coal})
    fams_a = sorted({h[3] for h in annih})
    _wrap = lambda names: "\n".join(textwrap.wrap(", ".join(names), 46))
    ax1.annotate("union over " + str(len(fams_c))
                 + " families, all with head $n{+}m$:\n"
                 + _wrap(fams_c)
                 + "\n\nhatched: pair declared, product outside the domain",
                 xy=(0.5, -0.15), xycoords="axes fraction", ha="center",
                 va="top", fontsize=_PLOT_FONTSIZE, color=_GREY_TEXT,
                 linespacing=1.5)
    ax2.annotate("union over " + str(len(fams_a)) + " families:\n"
                 + _wrap(fams_a),
                 xy=(0.5, -0.15), xycoords="axes fraction", ha="center",
                 va="top", fontsize=_PLOT_FONTSIZE, color=_GREY_TEXT,
                 linespacing=1.5)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=_DPI, bbox_inches="tight", facecolor="white")
    print(f"coalescence hyperedges {len(coal)}   "
          f"annihilation hyperedges {len(annih)}   total arcs "
          f"{G.number_of_edges()}")
    print(f"figure -> {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-max", type=int, default=8)
    ap.add_argument("--out", type=Path,
                    default=Path("../docs/Formulation/rag/eurofer_rag_hyperedges.png"))
    a = ap.parse_args()
    build_figure(a.n_max, a.out)
