#!/usr/bin/env python3
r"""NetworkX-generated TikZ RAG for the REGIONAL scale (Step 2).

Reads the REAL model structure from the genealogy package
(regional.rstate_index region declarations, populations.child_cluster)
and emits an \input-able TikZ fragment in the style of rag_twosex.tex.

Representative subset: n=2 clusters (C_0..C_2) x 2 sexes x 4 Census regions
(24 state vertices; the full model has 120 states at n=14).

Vertices: 24 state vertices (d, sex, region), the immigration source,
death/emigration sinks, and per-region mating-pair reaction diamonds
M_qr(r) (9 per region, 36 total).

Edges (all real terms of genealogy.regional.rhs_regional):
- mate_in (brown dashed, catalytic): (q,pat,r)->M_qr(r), (r,mat,r)->M_qr(r);
  WITHIN-REGION mating (region assortativity = 1, assumed); parents not
  consumed.
- child (shallowest=orange / deepest=teal / both=dark): M_qr(r) -> (c,pat,r),
  (c,mat,r) with c = populations.child_cluster(q,r,rule) -- actual code.
- imm (dark blue, thick): SRC -> (0,sex,r), labeled with the regional
  allocation sigma*I(t)*lam_r(t) / (1-sigma)*I(t)*lam_r(t).
- death (gray) / emig (gray dotted): (d,sex,r) -> sinks.
- transfer (purple, curved): inter-regional migration (d,sex,r)->(d,sex,r'),
  r != r'. The full model has 6 parallel edges per ordered pair (one per
  (d,sex)); ONE representative curved arrow per ordered pair is drawn,
  labeled m_{r->r'}, with the legend stating the 6-fold multiplicity.

Layout: states pinned on a region-band x cluster grid; mating diamonds via
networkx.spring_layout (fixed seed) clamped to their region band.
"""

import os

os.environ.setdefault("GENEALOGY_MAX_DEPTH", "14")

import numpy as np
import networkx as nx

import genealogy.populations as P
from genealogy import regional as REG

from make_rag_figures import (
    TIKZ_HEAD, EDGE_STYLE, _emit_edges, _legend, _nid, SEED,
)

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "drive_doc", "figs")

SEXES = ("pat", "mat")
RULES = ("shallowest_inheritance", "deepest_inheritance")

REGION_COLORS = {  # (name, RGB)
    "NE": ("Northeast", (31, 119, 180)),
    "MW": ("Midwest", (255, 127, 14)),
    "S": ("South", (44, 160, 44)),
    "W": ("West", (214, 39, 40)),
}
REGION_SHORT = {"NE": "NE", "MW": "MW", "S": "S", "W": "W"}


def _set_depth(n):
    P.MAX_DEPTH = n
    P.N_CLUSTERS = n + 1


def build_regional(n):
    """Regional RAG at depth n: (d,sex,r) states + source/sinks +
    within-region mating diamonds + immigration allocation + death/emig +
    inter-regional transfer edges (all 6*(d,sex) per ordered pair)."""
    _set_depth(n)
    G = nx.DiGraph()
    for r in REG.REGIONS:
        for d in range(n + 1):
            for sex in SEXES:
                G.add_node(("S", d, sex, r), kind="state", d=d, sex=sex, r=r)
    G.add_node("SRC", kind="source")
    G.add_node("DEATH", kind="sink")
    G.add_node("EMIG", kind="sink")

    # within-region mating reaction vertices + child edges (actual model code)
    for r in REG.REGIONS:
        for q in range(n + 1):
            for rr_ in range(n + 1):
                m = ("M", q, rr_, r)
                G.add_node(m, kind="rxn", q=q, r=rr_, reg=r)
                G.add_edge(("S", q, "pat", r), m, kind="mate_in")
                G.add_edge(("S", rr_, "mat", r), m, kind="mate_in")
                for rule in RULES:
                    c = P.child_cluster(q, rr_, rule)
                    G.add_edge(m, ("S", c, "pat", r), kind="child", rule=rule)
                    G.add_edge(m, ("S", c, "mat", r), kind="child", rule=rule)

    # immigration into C_0 of every region and sex (regional allocation)
    for r in REG.REGIONS:
        for sex in SEXES:
            G.add_edge("SRC", ("S", 0, sex, r), kind="imm", reg=r, sex=sex)

    # sinks
    for r in REG.REGIONS:
        for d in range(n + 1):
            for sex in SEXES:
                G.add_edge(("S", d, sex, r), "DEATH", kind="death")
                G.add_edge(("S", d, sex, r), "EMIG", kind="emig")

    # inter-regional migration transfers: (d,sex,r1)->(d,sex,r2), all of them
    for r1 in REG.REGIONS:
        for r2 in REG.REGIONS:
            if r1 == r2:
                continue
            for d in range(n + 1):
                for sex in SEXES:
                    G.add_edge(("S", d, sex, r1), ("S", d, sex, r2),
                               kind="transfer", r1=r1, r2=r2)
    return G


def layout_regional(G, n, width_cm=16.0):
    """Pin states on a 2x2 region grid (cluster along x inside each cell);
    relax mating diamonds inside their cell."""
    dx = 1.55
    cell_w, cell_h, lane = 6.2, 4.6, 0.78
    init, fixed = {}, []
    for ri, r in enumerate(REG.REGIONS):
        gx, gy = ri % 2, ri // 2
        x0c = gx * cell_w
        y0c = (1 - gy) * cell_h
        for d in range(n + 1):
            for sex, y in (("pat", y0c + lane), ("mat", y0c - lane)):
                k = ("S", d, sex, r)
                init[k] = (x0c + d * dx, y)
                fixed.append(k)
    xmax = cell_w + n * dx
    init["SRC"] = (-1.9, cell_h / 2); fixed.append("SRC")
    init["DEATH"] = (xmax + 1.6, cell_h / 2 + 0.9); fixed.append("DEATH")
    init["EMIG"] = (xmax + 1.6, cell_h / 2 - 0.9); fixed.append("EMIG")
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") == "rxn":
            ri = REG.REGION_INDEX[dat["reg"]]
            gx, gy = ri % 2, ri // 2
            init[nd] = (gx * cell_w + n * dx / 2,
                        (1 - gy) * cell_h)
    clamp = {}
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") == "rxn":
            ri = REG.REGION_INDEX[dat["reg"]]
            gx, gy = ri % 2, ri // 2
            x0c = gx * cell_w
            y0c = (1 - gy) * cell_h
            clamp[nd] = (x0c - 0.5, x0c + n * dx + 0.5, y0c - 1.0, y0c + 1.0)
    # spring layout with fixed nodes
    pos = nx.spring_layout(G, pos=init, fixed=fixed, seed=SEED,
                           iterations=150)
    for k, box in clamp.items():
        if k in pos:
            x, y = pos[k]
            (x0, x1, y0, y1) = box
            pos[k] = (min(max(x, x0), x1), min(max(y, y0), y1))
    xs = np.array([p[0] for p in pos.values()])
    ys = np.array([p[1] for p in pos.values()])
    x0, x1 = xs.min(), xs.max()
    y0, y1 = ys.min(), ys.max()
    sx = width_cm / max(x1 - x0, 1e-9)
    out = {k: ((x - x0) * sx, (y - y0) * sx) for k, (x, y) in pos.items()}
    return out, (x1 - x0) * sx, (y1 - y0) * sx


def _cell_bbox(pos, r, n):
    """Bounding box of region r's state vertices (pre-label, in cm)."""
    xs = [pos[("S", d, sex, r)][0] for d in range(n + 1) for sex in SEXES]
    ys = [pos[("S", d, sex, r)][1] for d in range(n + 1) for sex in SEXES]
    return (min(xs), max(xs), min(ys), max(ys))


def _emit_transfer_representatives(pos, n):
    """One curved representative arrow per ordered region pair, labeled
    m_{r1->r2}, routed AROUND the outside of the 2x2 cell grid so no arrow
    crosses a region cell. 8 perimeter arrows + 4 diagonal corner arcs."""
    lines = ["  % inter-regional migration transfers (representative: each",
             "  % drawn arrow stands for 6 parallel (d,sex) transfer edges)"]
    bb = {r: _cell_bbox(pos, r, n) for r in REG.REGIONS}
    # NE top-left, MW top-right, S bottom-left, W bottom-right
    ne, mw, s, w = (bb[r] for r in ("NE", "MW", "S", "W"))

    def arrow(a, b, opts, lab):
        (x1, y1), (x2, y2) = a, b
        lines.append(
            f"  \\draw[->, ragTrans, semithick, {opts}] "
            f"({x1:.3f},{y1:.3f}) to ({x2:.3f},{y2:.3f})"
            f" node[midway, font=\\tiny, fill=white, inner sep=0.5pt] "
            f"{{$m_{{{lab}}}$}};")

    # top row: NE <-> MW above the cells
    yT = max(ne[3], mw[3]) + 0.85
    arrow(((ne[0]+ne[1])/2, yT), ((mw[0]+mw[1])/2, yT),
          "bend left=14", "NE\\to MW")
    arrow(((mw[0]+mw[1])/2, yT+0.42), ((ne[0]+ne[1])/2, yT+0.42),
          "bend right=14", "MW\\to NE")
    # bottom row: S <-> W below the cells
    yB = min(s[2], w[2]) - 0.85
    arrow(((s[0]+s[1])/2, yB), ((w[0]+w[1])/2, yB),
          "bend right=14", "S\\to W")
    arrow(((w[0]+w[1])/2, yB-0.42), ((s[0]+s[1])/2, yB-0.42),
          "bend left=14", "W\\to S")
    # left column: NE <-> S left of the cells
    xL = min(ne[0], s[0]) - 0.85
    arrow((xL, (ne[2]+ne[3])/2), (xL, (s[2]+s[3])/2),
          "bend right=14", "NE\\to S")
    arrow((xL-0.42, (s[2]+s[3])/2), (xL-0.42, (ne[2]+ne[3])/2),
          "bend left=14", "S\\to NE")
    # right column: MW <-> W right of the cells
    xR = max(mw[1], w[1]) + 0.85
    arrow((xR, (mw[2]+mw[3])/2), (xR, (w[2]+w[3])/2),
          "bend left=14", "MW\\to W")
    arrow((xR+0.42, (w[2]+w[3])/2), (xR+0.42, (mw[2]+mw[3])/2),
          "bend right=14", "W\\to MW")
    # diagonals, each around its own outer corner
    arrow((ne[1]+0.35, ne[3]+0.35), (w[1]+0.35, w[2]-0.35),
          "out=25, in=-25, looseness=1.5", "NE\\to W")
    arrow((w[0]-0.35, w[2]-0.35), (ne[0]-0.35, ne[3]+0.35),
          "out=205, in=155, looseness=1.5", "W\\to NE")
    arrow((mw[0]-0.35, mw[3]+0.35), (s[0]-0.35, s[2]-0.35),
          "out=155, in=205, looseness=1.5", "MW\\to S")
    arrow((s[1]+0.35, s[2]-0.35), (mw[1]+0.35, mw[3]+0.35),
          "out=-25, in=25, looseness=1.5", "S\\to MW")
    lowest = min(yB - 0.42, s[2] - 0.35, w[2] - 0.35)
    return lines, lowest


def _emit_imm_stubs(G, pos):
    """Immigration: one short labeled stub per (region, sex) into C_0.
    The model has 8 allocation edges from the single source; the legend
    records the source and the multiplicity."""
    lines = ["  % immigration stubs (regional allocation, labeled)"]
    for r in REG.REGIONS:
        for sex in SEXES:
            x, y = pos[("S", 0, sex, r)]
            lines.append(
                f"  \\draw[->, ragImm, thin] ({x-0.62:.3f},{y:.3f}) -- "
                f"({x-0.06:.3f},{y:.3f});")
            frac = "{\\sigma}" if sex == "pat" else "{(1-\\sigma)}"
            dy = 0.30 if sex == "pat" else -0.30
            anch = "south" if sex == "pat" else "north"
            lines.append(
                f"  \\node[font=\\tiny, text=ragImm, anchor={anch}] "
                f"at ({x-0.34:.3f},{y+dy:.3f}) {{${frac}I\\lambda_{{{r}}}$}};")
    return lines


def _emit_src_key(y):
    return [
        "  % immigration source (all 8 regional allocation edges",
        "  % originate here; one stub per (region, sex) drawn at C_0)",
        f"  \\node[ragSrc] (srcKey) at (1.0,{y:.3f}) {{$\\varnothing$}};",
        f"  \\node[font=\\scriptsize, anchor=west] at (1.7,{y:.3f}) "
        f"{{immigration source $\\varnothing$: "
        f"$\\sigma I\\lambda_r$ (pat), $(1-\\sigma)I\\lambda_r$ (mat) "
        f"into each region's $C_0$ ($\\times 8$ edges)}};",
    ]


def _emit_sink_stubs(G, pos):
    """Death/emigration: one short stub per state (24 parallel edges each
    to the shared sink in the model). Full edge list stays in the graph;
    the legend records the multiplicity."""
    lines = ["  % death / emigration stubs (one per state)"]
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") != "state":
            continue
        x, y = pos[nd]
        lines.append(
            f"  \\draw[->, ragSink!70, thin] ({x:.3f},{y:.3f}) -- "
            f"({x-0.28:.3f},{y-0.34:.3f});")
        lines.append(
            f"  \\draw[->, ragSink!70, thin, dotted] ({x:.3f},{y:.3f}) -- "
            f"({x+0.28:.3f},{y-0.34:.3f});")
    return lines


def _emit_sink_key(y):
    """Shared sink vertices placed by the legend, with one representative
    edge each."""
    return [
        "  % shared sink vertices (all 24 death / 24 emigration edges",
        "  % of the model terminate here; one representative drawn)",
        f"  \\node[ragState, rectangle, fill=ragR0, fill opacity=0.3, "
        f"minimum width=0.5cm, minimum height=0.36cm] (sinkKeyS) "
        f"at (1.0,{y:.3f}) {{$C_d$}};",
        f"  \\node[ragSinkSt] (sinkKeyD) at (2.6,{y+0.42:.3f}) {{$\\mu$}};",
        f"  \\node[ragSinkSt] (sinkKeyE) at (2.6,{y-0.42:.3f}) "
        f"{{$\\varepsilon$}};",
        f"  \\draw[->, ragSink!70, thin] (sinkKeyS) -- (sinkKeyD);",
        f"  \\draw[->, ragSink!70, thin, dotted] (sinkKeyS) -- (sinkKeyE);",
        f"  \\node[font=\\scriptsize, anchor=west] at (3.1,{y:.3f}) "
        f"{{death $\\mu$ / emigration $\\varepsilon$ to shared sinks "
        f"($\\times24$ parallel edges each)}};",
    ]


def _emit_nodes_regional(G, pos, node_cm):
    lines = []
    for nd, dat in G.nodes(data=True):
        x, y = pos[nd]
        kind = dat.get("kind")
        nid = _nid(nd)
        if kind == "state":
            d, sex, r = dat["d"], dat["sex"], dat["r"]
            shape = "rectangle" if sex == "pat" else "ellipse"
            gi = REG.REGION_INDEX[r]
            fill = f"ragR{gi}"
            lines.append(
                f"  \\node[ragState, {shape}, fill={fill}, fill opacity=0.30, "
                f"text opacity=1, minimum width={node_cm:.2f}cm, "
                f"minimum height={node_cm*0.72:.2f}cm] "
                f"({nid}) at ({x:.3f},{y:.3f}) {{$C_{{{d}}}$}};")
        elif kind == "rxn":
            lines.append(f"  \\node[ragRxn] ({nid}) at ({x:.3f},{y:.3f}) {{}};")
        elif kind == "source":
            continue  # drawn beside the legend by _emit_src_key
        elif kind == "sink":
            continue  # drawn beside the legend by _emit_sink_key
    return lines


def _region_backgrounds(G, pos, n, pad=0.42):
    lines = ["  % region band backgrounds"]
    for ri, r in enumerate(REG.REGIONS):
        xs = [pos[("S", d, sex, r)][0]
              for d in range(n + 1) for sex in SEXES]
        ys = [pos[("S", d, sex, r)][1]
              for d in range(n + 1) for sex in SEXES]
        x0, x1 = min(xs) - pad, max(xs) + pad
        y0, y1 = min(ys) - pad, max(ys) + pad
        name, _ = REGION_COLORS[r]
        lines.append(
            f"  \\draw[fill=ragR{ri}, fill opacity=0.08, draw=ragR{ri}, "
            f"draw opacity=0.35, rounded corners=4pt] "
            f"({x0:.3f},{y0:.3f}) rectangle ({x1:.3f},{y1:.3f});")
        lines.append(
            f"  \\node[font=\\tiny\\bfseries, text=ragR{ri}!70!black, anchor=west] "
            f"at ({x0:.3f},{y1+0.30:.3f}) {{{name}}};")
    return lines


def write_regional(path):
    n = 2
    G = build_regional(n)
    pos, w, h = layout_regional(G, n)
    node_cm = 0.55

    groupcolors = "".join(
        f"\\definecolor{{ragR{ri}}}{{RGB}}{{{rgb[0]},{rgb[1]},{rgb[2]}}}\n"
        for ri, r in enumerate(REG.REGIONS)
        for rgb in [REGION_COLORS[r][1]])
    groupcolors += ("\\definecolor{ragTrans}{RGB}{120,40,160}\n")
    extra = (
        "% Regional RAG: 24 state vertices (d,sex,r), d=0..2; full model has\n"
        "% 120 states (n=14). Mating is WITHIN-REGION (region assortativity\n"
        "% lambda_region = 1, assumed). Inter-regional migration: 12 ordered\n"
        f"% region pairs x 6 (d,sex) = 72 transfer edges in the model; one\n"
        "% representative curved arrow per ordered pair is drawn, labeled\n"
        "% m_{r1->r2}.\n")
    head = TIKZ_HEAD.format(
        title=f"Regional RAG, n={n} (4 Census regions; representative subset)",
        seed=SEED, extra_comment=extra, groupcolors=groupcolors)

    body = []
    body.extend(_region_backgrounds(G, pos, n))
    body.extend(_emit_nodes_regional(G, pos, node_cm))
    # edges: draw everything except transfers (representatives instead)
    Gdraw = G.copy()
    for u, v, dat in list(Gdraw.edges(data=True)):
        if dat.get("kind") in ("transfer", "imm", "death", "emig"):
            Gdraw.remove_edge(u, v)
    body.extend(_emit_edges(Gdraw, imm_style="imm_thin"))
    body.extend(_emit_imm_stubs(G, pos))
    body.extend(_emit_sink_stubs(G, pos))
    tlines, tlowest = _emit_transfer_representatives(pos, n)
    body.extend(tlines)
    ymin = min(min(p[1] for p in pos.values()), tlowest)

    # source/sink labels

    legend = _legend([
        ("->, ragImm, thin",
         "immigration into $C_0$: $\\sigma I\\lambda_r$ (pat), "
         "$(1-\\sigma)I\\lambda_r$ (mat)"),
        ("->, ragTrans, semithick",
         "inter-regional migration transfer $m_{r_1\\to r_2}$ "
         "(one drawn per ordered pair = 6 parallel $(d,\\mathrm{sex})$ edges)"),
        ("->, ragMate!75, thin, dashed",
         "mating: WITHIN-REGION, parents catalytic (not consumed)"),
        ("->, ragShallowest, semithick", "child edge, shallowest inheritance"),
        ("->, ragDeepest, semithick", "child edge, deepest inheritance"),
        ("->, ragBoth, semithick", "child edge, both rules agree"),
    ], 0.0, ymin - 1.0)
    body.extend(legend)
    body.extend(_emit_src_key(ymin - 3.8))
    body.extend(_emit_sink_key(ymin - 4.9))
    body.append(
        f"  \\node[font=\\scriptsize, anchor=north west] at (0.0, {ymin - 6.1:.3f}) "
        f"{{\\tikz\\node[ragState, rectangle, fill=ragR0, fill opacity=0.3, "
        f"minimum width=0.5cm, minimum height=0.34cm] {{$C_d$}}; paternal "
        f"$\\;\\;$"
        f"\\tikz\\node[ragState, ellipse, fill=ragR0, fill opacity=0.3, "
        f"minimum width=0.5cm, minimum height=0.34cm] {{$C_d$}}; maternal "
        f"$\\;\\;$"
        f"\\tikz\\node[ragRxn] {{}}; mating-pair vertex $M_{{qr}}(r)$}};")

    tex = head + "\n".join(body) + "\n\\end{tikzpicture}\n"
    with open(path, "w") as f:
        f.write(tex)

    # edge census for the header/log
    from collections import Counter
    kinds = Counter()
    for u, v, dat in G.edges(data=True):
        k = dat.get("kind")
        kinds[f"child[{dat['rule']}]" if k == "child" else k] += 1
    nn = Counter(dat.get("kind") for _, dat in G.nodes(data=True))
    print(f"regional n={n}: {G.number_of_nodes()} vertices "
          f"({dict(nn)}), {G.number_of_edges()} edges ({dict(kinds)})")
    return G


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    p = os.path.join(FIGDIR, "rag_regional.tex")
    write_regional(p)
    print("wrote", p)


if __name__ == "__main__":
    main()
