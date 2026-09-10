#!/usr/bin/env python
"""make_rag_figure.py — node-link drawing of the EUROFER-97 RAG.

Renders the reaction admissibility graph declared in
``py_utils/materials/eurofer97`` as a three-lane ladder: one lane per
population, cluster size along x, one arc per declared edge family.

The drawing is generated from the MATERIALIZED graph
(``py_utils.core.to_networkx``), not from a hand-transcribed edge list, so
it cannot drift from the declaration.  If a family is added to the host
declaration it appears here on the next run.

Binary classes (COALESCENCE, ANNIHILATION) are omitted: they are directed
hyperedges and O(n_max**2) in count -- 6564 of the 7078 arcs at n_max=30 --
and drawing them turns the figure into a solid block.  Their count is
reported on the figure instead.

ENCODING.  Colour carries the EDGE CLASS here, not the population, because
the edge classes are what must be told apart: several share a vertex pair
and are drawn as nested arcs millimetres from each other.  Population
identity is carried instead by a tinted lane band plus a direct lane label
-- a large, low-saturation mark that cannot be confused with a thin
saturated arc, and which leaves the whole categorical budget for the arcs.
The lane tints are the suite population colours from
py_utils/visualization.py, so the association with the physics figures
survives.

The four arc hues are reference categorical slots 1/2/3/7, validated under
the ALL-PAIRS gate (worst CVD dE 9.2, worst normal-vision dE 16.3).
All-pairs is the right gate here because nested arcs sit side by side
rather than in a legend-ordered stack.  Line style is a redundant second
channel, so the figure survives greyscale printing.

Style otherwise follows py_utils/visualization.py: titles suppressed (the
document captions the figure) and type set once via rcParams at 22/20/16.

Usage::

    cd RadCluster_2_1
    python codes/make_rag_figure.py --n-max 12 \
        --out ../docs/Formulation/rag/eurofer_rag_figure.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from py_utils.input_data import InputData                      # noqa: E402
from py_utils.reaction_rates import ReactionRates              # noqa: E402
from py_utils.materials import build_eurofer_rag               # noqa: E402
from py_utils.core.to_networkx import to_networkx              # noqa: E402

# ── house style (py_utils/visualization.py) ──────────────────────────────
_LABEL_FONTSIZE = 22
_TICK_FONTSIZE = 20
_PLOT_FONTSIZE = 20
_NODE_FONTSIZE = 26   # the size label inside a vertex
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

_INK = "#1a1a1a"        # text tokens, and the cross-lane transfer arrows
_NODE_RING = "#3f3f3f"
_SENTINEL = "#b4b4b4"   # sentinel node outline
_FAN_SIA = "#1a1a1a"    # source/sink fans touching an SIA population
_FAN_VAC = "#1a4f8a"    # ... and a vacancy population
_GREY_TEXT = "#6a6a6a"

# Suite population colours, used here only as a 12%-alpha lane tint.
_POP = {
    "bulk-111": ("steelblue",  r"$\frac{1}{2}\langle 111\rangle$ SIA loops"),
    "bulk-100": ("darkviolet", r"$\langle 100\rangle$ SIA loops"),
    "bulk":     ("tomato",     "vacancy clusters / cavities"),
}

# Validated all-pairs categorical set for the four classes that co-locate
# as nested arcs.  inter_population and solute_trapping are NOT given
# categorical hues: they are geometrically unmistakable (a vertical arrow
# between lanes; a closed loop under a vertex), so spending scarce
# validated slots on them would degrade the four that need separating.
_EDGE = {
    "growth":           ("#2a78d6", "-",  _LW),
    "shrinkage":        ("#eb6834", "--", _LW_THIN),
    "dissociation":     ("#1baf7a", ":",  _LW_THIN),
    "recombination":    ("#4a3aa7", "-.", _LW_THIN),
    "inter_population": (_INK,      "-",  _LW),
    "solute_trapping":  ("#4a4a4a", "-", _LW_THIN),
    "source":           (_FAN_SIA, "-",            1.4),
    "sink":             (_FAN_SIA, (0, (7, 3)), 1.4),
}

# Opened up on both axes: the previous 1.0/1.0 spacing left a 0.4-unit gap
# between circles, far too short for a dash pattern to read.
_X_SCALE = 3.00
_LANE_Y = {"bulk-111": 6.2, "bulk-100": 3.1, "bulk": 0.0}
_R = 0.44


def _node_key(node: str):
    """'sia:bulk-111:7' -> ('bulk-111', 7); sentinels -> (None, None)."""
    if node in ("SOURCE", "SINK"):
        return None, None
    _, pop, n = node.split(":")
    return pop, int(n)


# ── worked examples of the two omitted binary classes ────────────────────
# The lane ladder cannot show COALESCENCE or ANNIHILATION in place (they are
# O(n_max^2) hyperedges), but one instance of each, drawn full size above the
# ladder, states the two-tails/one-head convention that the count in the
# header would otherwise leave abstract.
_EX_R = 0.64
_POP_SHORT = {
    "bulk-111": r"$\frac{1}{2}\langle 111\rangle$",
    "bulk-100": r"$\langle 100\rangle$",
    "bulk": "vacancy",
}
# Compact coloured marks, matching the companion hyperedge figure's
# schematic.  Only the reaction VERTEX is coloured, never the arrows: a
# coloured arrow up here would be read as one of the lane arc classes.
_HYPER = {"coalescence": "#256abf", "annihilation": "#d03b3b"}


def find_hyperedge(rag, edge_class, tails, head):
    """Assert that (tails -> head) is really declared, and name its family.

    The examples are hand-picked sizes, so unlike everything else in this
    figure they could drift from the declaration.  This looks them up in a
    graph materialized wide enough to contain them and returns the family
    label, so a stale example fails loudly instead of drawing a reaction
    the host does not admit.
    """
    from collections import defaultdict
    want_t, want_h = sorted(tails), tuple(head)
    span = max(max(n for _, n in tails), head[1]) + 1
    G = to_networkx(rag, n_max=span, include_binary=True)
    groups = defaultdict(list)
    for u, v, d in G.edges(data=True):
        if d.get("hyper_id") is not None:
            groups[d["hyper_id"]].append((u, v, d))
    for arcs in groups.values():
        d = arcs[0][2]
        if d["edge_class"] != edge_class:
            continue
        if sorted(_node_key(a[0]) for a in arcs) == want_t \
                and _node_key(arcs[0][1]) == want_h:
            return d["label"]
    raise LookupError(
        f"{edge_class} {tails} -> {head} is not declared; the example is "
        f"stale or the host declaration changed")


def _draw_example(ax, x0, y, tails, head, edge_class, label, n_max,
                  note=None):
    """One hyperedge: two tails -> reaction vertex -> one head."""
    colour = _HYPER[edge_class]
    tail_xy = [(x0, y + 1.25), (x0, y - 1.25)]
    sq, hx = x0 + 3.5, x0 + 7.0

    def vertex(cx, cy, pop, n, side):
        ax.add_patch(Circle((cx, cy), _EX_R, facecolor=_POP[pop][0],
                            alpha=0.30, edgecolor="none", zorder=5))
        ax.add_patch(Circle((cx, cy), _EX_R, facecolor="none",
                            edgecolor=_NODE_RING, lw=_LW_THIN, zorder=6))
        ax.text(cx, cy, str(n), ha="center", va="center", color=_INK,
                fontweight="bold", fontsize=_NODE_FONTSIZE, zorder=7)
        dx = -(_EX_R + 0.3) if side == "left" else (_EX_R + 0.3)
        ax.text(cx + dx, cy, _POP_SHORT[pop], ha=side and
                ("right" if side == "left" else "left"), va="center",
                color=_GREY_TEXT, fontsize=_PLOT_FONTSIZE, zorder=7)

    for (tx, ty), (pop, n) in zip(tail_xy, tails):
        vertex(tx, ty, pop, n, "left")
        ax.add_patch(FancyArrowPatch((tx, ty), (sq, y), arrowstyle="-|>",
                                     mutation_scale=17, lw=_LW_THIN,
                                     color=colour, shrinkA=38, shrinkB=14,
                                     zorder=4))
    ax.add_patch(Rectangle((sq - 0.42, y - 0.42), 0.84, 0.84,
                           facecolor=colour, edgecolor="none", zorder=6))
    ax.add_patch(FancyArrowPatch((sq, y), (hx, y), arrowstyle="-|>",
                                 mutation_scale=20, lw=_LW, color=colour,
                                 shrinkA=14, shrinkB=38, zorder=4))
    vertex(hx, y, head[0], head[1], "right")

    ax.text(x0 + 3.5, y + 2.45, f"{edge_class}  ·  {label}", ha="center",
            va="bottom", color=_INK, fontsize=_PLOT_FONTSIZE, zorder=8,
            bbox=dict(facecolor="white", edgecolor="none", pad=2.0))
    if note:
        ax.text(x0 + 3.5, y - 2.45, note, ha="center", va="top",
                color=_GREY_TEXT, fontsize=_PLOT_FONTSIZE, zorder=8,
                bbox=dict(facecolor="white", edgecolor="none", pad=2.0))

    # Leaders tie every example vertex that IS a ladder vertex back to it --
    # the head as much as the tails, and on WHICHEVER lane it lives.  The
    # annihilation example spans polarities, so restricting leaders to the
    # <111> lane (as this did originally) left its vacancy tail floating:
    # the one vertex whose lane the reader most needs pointed out, because
    # it is the reason the reaction cannot be drawn inside a single lane.
    # Its leader crosses the other lanes, which is honest -- the reaction
    # crosses them too.
    #
    # The only vertex that still gets none is a head beyond the cut-off (the
    # coalescence head n+m = 11): there is no rung to point at, and that
    # absence is the truncation the caption describes.
    leader_pts = list(zip(tail_xy, tails)) + [((hx, y), tuple(head))]
    for (cx, cy), (pop, n) in leader_pts:
        if pop not in _LANE_Y or n > n_max:
            continue
        # Strong black dash: this is a callout tying the example to the
        # ladder, and it has to survive crossing the source/sink fans.  The
        # dash period is deliberately longer than the sink arcs' (7,3) so
        # the two do not read as the same class.
        ax.plot([cx, _X_SCALE * n], [cy - _EX_R, _LANE_Y[pop] + _R],
                ls=(0, (9, 4)), lw=2.2, color=_INK, alpha=0.95, zorder=2)


def build_figure(n_max: int, out: Path, examples: bool = True) -> None:
    inp = InputData(I=200, V=200, physics_option="full_CD_fission")
    rag, _ = build_eurofer_rag(inp, ReactionRates(inp))

    G = to_networkx(rag, n_max=n_max, include_binary=False)
    Gfull = to_networkx(rag, n_max=n_max, include_binary=True)
    n_binary = Gfull.number_of_edges() - G.number_of_edges()

    pos = {}
    for node in G.nodes:
        pop, n = _node_key(node)
        if pop is not None:
            pos[node] = (_X_SCALE * n, _LANE_Y[pop])
    x_lo, x_hi = 0.0, _X_SCALE * (n_max + 1)
    pos["SOURCE"] = (x_lo, _LANE_Y["bulk-111"] + 2.7)
    pos["SINK"] = (x_hi, _LANE_Y["bulk"] - 2.7)

    _w = 1.95 * n_max + 8.0
    fig, ax = plt.subplots(figsize=(_w, _w / 1.42 if examples else 12.0))

    # ── lane bands: population identity, large and low-saturation ────────
    for pop, y in _LANE_Y.items():
        ax.add_patch(Rectangle((x_lo - 0.6, y - 0.95),
                               x_hi - x_lo + 1.2, 1.9,
                               facecolor=_POP[pop][0], alpha=0.12,
                               edgecolor="none", zorder=0))

    # ── arcs ─────────────────────────────────────────────────────────────
    # Parallel families between the same vertex pair fan out by increasing
    # |rad| -- the graph is a MultiDiGraph precisely because several
    # classes share endpoints (cavity growth and trap mutation both
    # advance the vacancy lane).
    seen: dict[tuple, int] = {}
    for u, v, d in G.edges(data=True):
        ec = d["edge_class"]
        if u not in pos or v not in pos:
            continue
        upop, un = _node_key(u)
        vpop, vn = _node_key(v)
        colour, ls, lw = _EDGE.get(ec, (_INK, "-", _LW_THIN))
        alpha = 0.95
        zorder = 3
        if ec in ("source", "sink"):
            # Fan colour encodes the POLARITY of the population it touches:
            # the cascade feeds, and the sinks drain, the two sub-spaces
            # independently, and that split is the thing worth seeing.
            # source: head is the population vertex; sink: tail is.
            fan_pop = vpop if ec == "source" else upop
            colour = _FAN_VAC if fan_pop == "bulk" else _FAN_SIA
            alpha, zorder = 0.55, 1

        idx = seen.get((u, v), 0)
        seen[(u, v)] = idx + 1

        if u == v:
            # Solute trapping is a self-loop because the composition axis
            # is NOT expanded here (expanding it would multiply the node
            # count by |C|).  Drawn TANGENT to the vertex and above the
            # arc layer, on a white halo: the loss arcs bow up to 1.3
            # units below the lane, so an untangent loop drawn underneath
            # them either disappears or reads as a second vertex.
            x, y = pos[u]
            cx, cy, r = x, y - _R - 0.27, 0.27
            ax.add_patch(Circle((cx, cy), r, fill=False, ec="white",
                                lw=lw + 3.0, zorder=6))
            ax.add_patch(Circle((cx, cy), r, fill=False, ec=colour,
                                lw=lw, ls=ls, zorder=7))
            ax.plot([cx - r], [cy], marker=(3, 0, 90), markersize=8,
                    color=colour, zorder=7)
            continue

        if un is not None and vn is not None and upop == vpop:
            # growth bows above the lane, loss below; nested by family
            base = 0.24 if vn > un else -0.24
            rad = base * (1.0 + 0.62 * idx)
        else:
            rad = 0.10 * (1 + idx)

        ax.add_patch(FancyArrowPatch(
            pos[u], pos[v], connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>", mutation_scale=19,
            shrinkA=24, shrinkB=24,
            lw=lw, ls=ls, color=colour, alpha=alpha, zorder=zorder))

    # ── vertices ─────────────────────────────────────────────────────────
    # Neutral on purpose: the lane band already says which population this
    # is, and a saturated node would compete with the arc palette.
    for node, (x, y) in pos.items():
        pop, n = _node_key(node)
        if pop is None:
            continue                        # sentinels drawn as label boxes
        ax.add_patch(Circle((x, y), _R, facecolor="white",
                            edgecolor=_NODE_RING, lw=_LW_THIN, zorder=5))
        ax.text(x, y, str(n), ha="center", va="center", color=_INK,
                fontweight="bold", fontsize=_NODE_FONTSIZE, zorder=6)

    for pop, y in _LANE_Y.items():
        ax.text(x_lo - 0.9, y, _POP[pop][1], ha="right", va="center",
                color=_INK, fontsize=_LABEL_FONTSIZE)

    _box = dict(boxstyle="round,pad=0.45", fc="white", ec=_SENTINEL,
                lw=_LW_THIN)
    ax.text(*pos["SOURCE"], "cascade\nsource", ha="center", va="center",
            color=_GREY_TEXT, fontsize=_PLOT_FONTSIZE, bbox=_box, zorder=6)
    ax.text(*pos["SINK"], "sinks\n(dislocations, GBs, precipitates)",
            ha="center", va="center", color=_GREY_TEXT,
            fontsize=_PLOT_FONTSIZE, bbox=_box, zorder=6)

    # ── legend ───────────────────────────────────────────────────────────
    handles = [
        Line2D([], [], color=c, ls=ls, lw=max(lw, 2.6),
               label=name.replace("_", " "))
        for name, (c, ls, lw) in _EDGE.items()
        if name not in ("source", "sink")
    ]
    _dash = (0, (7, 3))
    handles += [
        Line2D([], [], color=_FAN_SIA, ls="-", lw=2.0, label="source → SIA"),
        Line2D([], [], color=_FAN_VAC, ls="-", lw=2.0, label="source → vacancy"),
        Line2D([], [], color=_FAN_SIA, ls=_dash, lw=2.0, label="sink ← SIA"),
        Line2D([], [], color=_FAN_VAC, ls=_dash, lw=2.0, label="sink ← vacancy"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.008), frameon=False,
               title="edge class", ncol=4)

    if examples:
        ey = _LANE_Y["bulk-111"] + 5.6
        lab_a = find_hyperedge(rag, "annihilation",
                               [("bulk", 2), ("bulk-111", 4)], ("bulk-111", 2))
        _draw_example(ax, x_lo + 2.0, ey,
                      [("bulk-111", 4), ("bulk", 2)], ("bulk-111", 2),
                      "annihilation", lab_a, n_max,
                      "the larger loop survives, shrunk by\nthe vacancy "
                      "content it absorbed")
        lab_c = find_hyperedge(rag, "coalescence",
                               [("bulk-111", 3), ("bulk-111", 8)],
                               ("bulk-111", 11))
        _draw_example(ax, x_lo + 15.6, ey,
                      [("bulk-111", 3), ("bulk-111", 8)], ("bulk-111", 11),
                      "coalescence", lab_c, n_max,
                      f"head $n{{+}}m=11$ lies beyond the cut-off\n"
                      f"$n_{{max}}={n_max}$: not an arc in the ladder")
    ax.set_xlabel("cluster size $n$")
    ax.set_xlim(x_lo - 1.4, x_hi + 1.4)
    top = _LANE_Y["bulk-111"] + (9.6 if examples else 4.3)
    ax.set_ylim(_LANE_Y["bulk"] - 4.4, top)
    ax.set_xticks([_X_SCALE * n for n in range(1, n_max + 1)])
    ax.set_xticklabels([str(n) for n in range(1, n_max + 1)])
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(_SENTINEL)
    ax.set_aspect("equal")

    # No title: suppressed suite-wide; the document captions the figure.
    fig.tight_layout(rect=(0.0, 0.15 if examples else 0.23, 1.0, 1.0))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=_DPI, bbox_inches="tight", facecolor="white")
    print(f"nodes {G.number_of_nodes()}  unary arcs {G.number_of_edges()}  "
          f"binary omitted {n_binary}")
    print(f"figure -> {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-max", type=int, default=8)
    ap.add_argument("--no-examples", action="store_true",
                    help="omit the two worked hyperedge examples")
    ap.add_argument("--out", type=Path,
                    default=Path("../docs/Formulation/rag/eurofer_rag_figure.png"))
    a = ap.parse_args()
    build_figure(a.n_max, a.out, examples=not a.no_examples)
