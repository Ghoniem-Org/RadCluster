#!/usr/bin/env python3
r"""Build NetworkX-generated TikZ figures of the genealogy-dynamics RAG.

Reads the REAL model structure from the genealogy package
(populations.child_cluster, religion group declarations) and emits
\input-able TikZ fragments with \node/\draw commands (no tikzplotlib).

Figures:
  (a) drive_doc/figs/rag_twosex.tex    -- two-sex RAG at n=4
  (b) drive_doc/figs/rag_religious.tex  -- religious RAG at n=2

Mating is drawn as a bipartite reaction graph: each parental pair (q,r)
gets a small diamond "reaction vertex" M_qr; the parents feed it with
dashed catalytic edges (not consumed) and it emits child edges into the
child cluster states under the shallowest and deepest inheritance rules
(shallowest = orange, deepest = teal; one dark edge when both rules agree).

Layout: networkx.spring_layout with fixed seed; state vertices are pinned
in depth lanes (grouped by depth, and by religious group for (b)) while
the reaction vertices relax around them. Colors defined inside each file.
"""

import os

os.environ.setdefault("GENEALOGY_MAX_DEPTH", "12")

import numpy as np
import networkx as nx

import genealogy.populations as P
from genealogy import religion as R

SEED = 20260926
FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "drive_doc", "figs")

SEXES = ("pat", "mat")
RULES = ("shallowest_inheritance", "deepest_inheritance")

# ---------------------------------------------------------------- palette
GROUP_COLORS = {  # (name, RGB)
    "pro": ("Protestant", (31, 119, 180)),
    "cat": ("Catholic", (255, 127, 14)),
    "mor": ("Mormon (LDS)", (44, 160, 44)),
    "jew": ("Jewish", (214, 39, 40)),
    "mus": ("Muslim", (148, 103, 189)),
    "dhr": ("Hindu/Buddhist", (140, 86, 75)),
    "oth": ("Orthodox/other", (227, 119, 194)),
    "una": ("Unaffiliated", (127, 127, 127)),
}


def _set_depth(n):
    """Point the model code at pedigree depth n (clusters C_0..C_n)."""
    P.MAX_DEPTH = n
    P.N_CLUSTERS = n + 1


# ---------------------------------------------------------------- builders
def build_twosex(n):
    """Two-sex RAG: state vertices (d,sex), source, two sinks, mating
    reaction vertices. Returns (DiGraph, edge_kind_counts_info)."""
    _set_depth(n)
    G = nx.DiGraph()
    for d in range(n + 1):
        for sex in SEXES:
            G.add_node(("S", d, sex), kind="state", d=d, sex=sex)
    G.add_node("SRC", kind="source")
    G.add_node("DEATH", kind="sink")
    G.add_node("EMIG", kind="sink")

    for q in range(n + 1):
        for r in range(n + 1):
            m = ("M", q, r)
            G.add_node(m, kind="rxn", q=q, r=r)
            G.add_edge(("S", q, "pat"), m, kind="mate_in")
            G.add_edge(("S", r, "mat"), m, kind="mate_in")
            for rule in RULES:
                c = P.child_cluster(q, r, rule)  # actual model code
                G.add_edge(m, ("S", c, "pat"), kind="child", rule=rule)
                G.add_edge(m, ("S", c, "mat"), kind="child", rule=rule)

    for sex in SEXES:
        G.add_edge("SRC", ("S", 0, sex), kind="imm")
    for d in range(n + 1):
        for sex in SEXES:
            G.add_edge(("S", d, sex), "DEATH", kind="death")
            G.add_edge(("S", d, sex), "EMIG", kind="emig")
    return G


def build_religious(n):
    """Religious RAG at depth n: (d,sex,g) states + source/sinks +
    within-group mating reaction vertices + disaffiliation transfers +
    one faint representative cross-group mating edge per ordered pair."""
    _set_depth(n)
    G = nx.DiGraph()
    for g in R.GROUPS:
        for d in range(n + 1):
            for sex in SEXES:
                G.add_node(("S", d, sex, g), kind="state", d=d, sex=sex, g=g)
    G.add_node("SRC", kind="source")
    G.add_node("DEATH", kind="sink")
    G.add_node("EMIG", kind="sink")

    # within-group mating (child takes the mother's group: here g2 == g1)
    for g in R.GROUPS:
        for q in range(n + 1):
            for r in range(n + 1):
                m = ("M", q, r, g)
                G.add_node(m, kind="rxn", q=q, r=r, g=g)
                G.add_edge(("S", q, "pat", g), m, kind="mate_in")
                G.add_edge(("S", r, "mat", g), m, kind="mate_in")
                for rule in RULES:
                    c = P.child_cluster(q, r, rule)  # actual model code
                    G.add_edge(m, ("S", c, "pat", g), kind="child", rule=rule)
                    G.add_edge(m, ("S", c, "mat", g), kind="child", rule=rule)

    # disaffiliation transfers g -> una (same depth and sex)
    for g in R.DISAFFILIATING_GROUPS:
        for d in range(n + 1):
            for sex in SEXES:
                G.add_edge(("S", d, sex, g), ("S", d, sex, "una"),
                           kind="disaff")

    # immigration into C_0 of every group and sex
    for g in R.GROUPS:
        for sex in SEXES:
            G.add_edge("SRC", ("S", 0, sex, g), kind="imm")

    # sinks
    for g in R.GROUPS:
        for d in range(n + 1):
            for sex in SEXES:
                G.add_edge(("S", d, sex, g), "DEATH", kind="death")
                G.add_edge(("S", d, sex, g), "EMIG", kind="emig")

    # cross-group mating: ONE faint representative edge per ordered pair
    # (g1 -> g2), g1 != g2 — the full set is G*(G-1)*(n+1)^2 combos.
    rng = np.random.default_rng(SEED)
    for g1 in R.GROUPS:
        for g2 in R.GROUPS:
            if g1 == g2:
                continue
            q = int(rng.integers(0, n + 1))
            r = int(rng.integers(0, n + 1))
            G.add_edge(("S", q, "pat", g1), ("S", r, "mat", g2),
                       kind="xmate", q=q, r=r)
    return G


# ---------------------------------------------------------------- layout
def _layout(G, init_pos, fixed, width_cm, clamp=None):
    pos = nx.spring_layout(G, pos=init_pos, fixed=fixed, seed=SEED,
                           iterations=300)
    if clamp is not None:
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
    out = {}
    for k, (x, y) in pos.items():
        out[k] = ((x - x0) * sx, (y - y0) * sx)
    return out, (x1 - x0) * sx, (y1 - y0) * sx


def layout_twosex(G, n, width_cm=15.0):
    dx, yp, ym = 2.6, 2.0, -2.0
    init, fixed = {}, []
    for d in range(n + 1):
        for sex, y in (("pat", yp), ("mat", ym)):
            k = ("S", d, sex)
            init[k] = (d * dx, y)
            fixed.append(k)
    init["SRC"] = (-2.4, 0.0); fixed.append("SRC")
    init["DEATH"] = (n * dx + 2.0, yp); fixed.append("DEATH")
    init["EMIG"] = (n * dx + 2.0, ym); fixed.append("EMIG")
    rng = np.random.default_rng(SEED)
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") == "rxn":
            q, r = dat["q"], dat["r"]
            init[nd] = ((q + r) / 2 * dx + rng.normal(0, 0.25),
                        rng.normal(0, 0.35))
    # keep the mating cloud between the lanes: clamp reaction vertices to
    # the state-node bounding box (plus padding)
    sx0 = -2.4
    sx1 = n * dx + 2.0
    clamp = {nd: (sx0, sx1, ym - 1.2, yp + 1.2)
             for nd, dat in G.nodes(data=True) if dat.get("kind") == "rxn"}
    return _layout(G, init, fixed, width_cm, clamp=clamp)


def layout_religious(G, n, width_cm=16.0):
    gw, dx, yp, ym = 2.0, 0.55, 1.15, -1.15
    init, fixed = {}, []
    for gi, g in enumerate(R.GROUPS):
        for d in range(n + 1):
            for sex, y in (("pat", yp), ("mat", ym)):
                k = ("S", d, sex, g)
                init[k] = (gi * gw + d * dx, y)
                fixed.append(k)
    xmax = (len(R.GROUPS) - 1) * gw + n * dx
    init["SRC"] = (-1.6, 0.0); fixed.append("SRC")
    init["DEATH"] = (xmax + 1.2, yp); fixed.append("DEATH")
    init["EMIG"] = (xmax + 1.2, ym); fixed.append("EMIG")
    rng = np.random.default_rng(SEED)
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") == "rxn":
            gi = R.GROUP_INDEX[dat["g"]]
            init[nd] = (gi * gw + n * dx / 2 + rng.normal(0, 0.3),
                        rng.normal(0, 0.3))
    # keep each group's mating cloud inside its cluster: clamp reaction
    # vertices to their own group's state-node bounding box (plus padding)
    clamp = {}
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") == "rxn":
            gi = R.GROUP_INDEX[dat["g"]]
            clamp[nd] = (gi * gw - 0.55, gi * gw + n * dx + 0.55,
                         ym - 1.0, yp + 1.0)
    return _layout(G, init, fixed, width_cm, clamp=clamp)


# ---------------------------------------------------------------- tikz
def _nid(k):
    if isinstance(k, str):
        return k
    return "-".join(str(x) for x in k)


TIKZ_HEAD = """% {title}
% Generated by make_rag_figures.py from the genealogy package (NetworkX, seed {seed}).
% Vertices: real model state variables + immigration source + death/emigration sinks
% + mating-pair reaction vertices (diamonds). Child cluster from
% genealogy.populations.child_cluster (shallowest / deepest inheritance).
% Requires: \\usepackage{{tikz}} \\usetikzlibrary{{arrows.meta,shapes}}
{extra_comment}\\definecolor{{ragPatFill}}{{RGB}}{{219,232,251}}
\\definecolor{{ragMatFill}}{{RGB}}{{179,205,245}}
\\definecolor{{ragImm}}{{RGB}}{{0,51,153}}
\\definecolor{{ragSink}}{{RGB}}{{110,110,110}}
\\definecolor{{ragMate}}{{RGB}}{{150,90,30}}
\\definecolor{{ragShallowest}}{{RGB}}{{230,110,20}}
\\definecolor{{ragDeepest}}{{RGB}}{{0,150,130}}
\\definecolor{{ragBoth}}{{RGB}}{{60,60,60}}
\\definecolor{{ragDis}}{{RGB}}{{130,60,160}}
\\definecolor{{ragX}}{{RGB}}{{170,170,170}}
{groupcolors}\\begin{{tikzpicture}}[>=Stealth,
  ragState/.style={{draw=ragSink, thick, font=\\tiny\\bfseries, inner sep=1pt}},
  ragRxn/.style={{diamond, draw=ragMate, fill=ragMate!25, inner sep=0pt, minimum size=3.2pt}},
  ragSrc/.style={{circle, draw=ragImm, fill=ragImm!12, thick, font=\\small, minimum size=0.55cm}},
  ragSinkSt/.style={{circle, draw=ragSink, fill=ragSink!12, font=\\scriptsize, minimum size=0.5cm}},
]
"""

EDGE_STYLE = {
    "imm": "->, ragImm, very thick",
    "imm_thin": "->, ragImm, thin",
    "death": "->, ragSink!70, thin",
    "emig": "->, ragSink!70, thin, dotted",
    "mate_in": "->, ragMate!75, thin, dashed",
    "child_shallowest": "->, ragShallowest, semithick",
    "child_deepest": "->, ragDeepest, semithick",
    "child_both": "->, ragBoth, semithick",
    "disaff": "->, ragDis, semithick, dotted",
    "xmate": "->, ragX!60, ultra thin, dashed",
}


def _emit_edges(G, imm_style="imm"):
    """Merge shallowest/deepest child edges to the same target into one."""
    lines = []
    seen_both = set()
    # first pass: find (u,v) pairs carrying both rules
    child_rules = {}
    for u, v, dat in G.edges(data=True):
        if dat.get("kind") == "child":
            child_rules.setdefault((u, v), set()).add(dat["rule"])
    for u, v, dat in G.edges(data=True):
        kind = dat.get("kind")
        if kind == "child":
            rules = child_rules[(u, v)]
            if len(rules) == 2:
                if (u, v) in seen_both:
                    continue
                seen_both.add((u, v))
                style = EDGE_STYLE["child_both"]
                bend = ""
            else:
                rule = dat["rule"]
                style = EDGE_STYLE["child_shallowest" if rule == "shallowest_inheritance"
                                   else "child_deepest"]
                bend = ", bend left=14" if rule == "shallowest_inheritance" else ", bend right=14"
            lines.append(f"  \\draw[{style}{bend}] ({_nid(u)}) -- ({_nid(v)});")
        elif kind == "imm":
            lines.append(f"  \\draw[{EDGE_STYLE[imm_style]}] ({_nid(u)}) -- ({_nid(v)});")
        elif kind in EDGE_STYLE:
            lines.append(f"  \\draw[{EDGE_STYLE[kind]}] ({_nid(u)}) -- ({_nid(v)});")
    return lines


def _emit_nodes_twosex(G, pos, node_cm):
    lines = []
    for nd, dat in G.nodes(data=True):
        x, y = pos[nd]
        kind = dat.get("kind")
        nid = _nid(nd)
        if kind == "state":
            d, sex = dat["d"], dat["sex"]
            shape = "rectangle" if sex == "pat" else "ellipse"
            fill = "ragPatFill" if sex == "pat" else "ragMatFill"
            lines.append(
                f"  \\node[ragState, {shape}, fill={fill}, "
                f"minimum width={node_cm:.2f}cm, minimum height={node_cm*0.72:.2f}cm] "
                f"({nid}) at ({x:.3f},{y:.3f}) {{$C_{{{d}}}$}};")
        elif kind == "rxn":
            lines.append(f"  \\node[ragRxn] ({nid}) at ({x:.3f},{y:.3f}) {{}};")
        elif kind == "source":
            lines.append(f"  \\node[ragSrc] ({nid}) at ({x:.3f},{y:.3f}) "
                         f"{{$\\varnothing$}};")
        elif kind == "sink":
            lab = "$\\mu$" if nd == "DEATH" else "$\\varepsilon$"
            lines.append(f"  \\node[ragSinkSt] ({nid}) at ({x:.3f},{y:.3f}) "
                         f"{{{lab}}};")
    return lines


def _emit_nodes_religious(G, pos, node_cm):
    lines = []
    for nd, dat in G.nodes(data=True):
        x, y = pos[nd]
        kind = dat.get("kind")
        nid = _nid(nd)
        if kind == "state":
            d, sex, g = dat["d"], dat["sex"], dat["g"]
            shape = "rectangle" if sex == "pat" else "ellipse"
            gi = R.GROUP_INDEX[g]
            lines.append(
                f"  \\node[ragState, {shape}, fill=ragG{gi}, fill opacity=0.28, "
                f"text opacity=1, minimum width={node_cm:.2f}cm, "
                f"minimum height={node_cm*0.72:.2f}cm] "
                f"({nid}) at ({x:.3f},{y:.3f}) {{$C_{{{d}}}$}};")
        elif kind == "rxn":
            lines.append(f"  \\node[ragRxn] ({nid}) at ({x:.3f},{y:.3f}) {{}};")
        elif kind == "source":
            lines.append(f"  \\node[ragSrc] ({nid}) at ({x:.3f},{y:.3f}) "
                         f"{{$\\varnothing$}};")
        elif kind == "sink":
            lab = "$\\mu$" if nd == "DEATH" else "$\\varepsilon$"
            lines.append(f"  \\node[ragSinkSt] ({nid}) at ({x:.3f},{y:.3f}) "
                         f"{{{lab}}};")
    return lines


GROUP_SHORT = {
    "pro": "Protestant", "cat": "Catholic", "mor": "Mormon",
    "jew": "Jewish", "mus": "Muslim", "dhr": "Hindu/Buddhist",
    "oth": "Other faith", "una": "Unaffiliated",
}


def _group_backgrounds(G, pos, pad=0.35):
    """Light rounded rectangles behind each religious group's state nodes."""
    lines = ["  % group cluster backgrounds"]
    for gi, g in enumerate(R.GROUPS):
        xs = [pos[("S", d, sex, g)][0]
              for d in range(P.N_CLUSTERS) for sex in SEXES]
        ys = [pos[("S", d, sex, g)][1]
              for d in range(P.N_CLUSTERS) for sex in SEXES]
        x0, x1, y0, y1 = min(xs) - pad, max(xs) + pad, min(ys) - pad, max(ys) + pad
        lines.append(
            f"  \\draw[fill=ragG{gi}, fill opacity=0.10, draw=ragG{gi}, "
            f"draw opacity=0.35, rounded corners=4pt] "
            f"({x0:.3f},{y0:.3f}) rectangle ({x1:.3f},{y1:.3f});")
        lines.append(
            f"  \\node[font=\\tiny\\bfseries, text=ragG{gi}!70!black] "
            f"at ({(x0+x1)/2:.3f},{y1+0.30:.3f}) {{{GROUP_SHORT[g]}}};")
    return lines


def _legend(lines_spec, x, y):
    """lines_spec: list of (tikz draw opts, label). Drawn as samples."""
    out = [f"  \\node[draw=ragSink!40, fill=white, rounded corners=3pt, "
           f"font=\\scriptsize, align=left, anchor=north west] (ragLegend) "
           f"at ({x:.3f},{y:.3f}) {{%;"]
    out.append("    \\begin{tabular}{@{}l@{\\hspace{4pt}}l@{}}")
    for opts, label in lines_spec:
        out.append(f"      \\tikz\\draw[{opts}] (0,0) -- (0.7,0); & {label} \\\\")
    out.append("    \\end{tabular}};")
    return out


def write_twosex(path):
    n = 4
    G = build_twosex(n)
    pos, w, h = layout_twosex(G, n)
    node_cm = 0.62

    groupcolors = ""
    head = TIKZ_HEAD.format(
        title=f"Two-sex RAG, n={n} (all vertices and edges)",
        seed=SEED, extra_comment="", groupcolors=groupcolors)

    body = []
    body.extend(_emit_nodes_twosex(G, pos, node_cm))
    body.extend(_emit_edges(G))

    # lane + source/sink labels
    dmax_x = max(pos[("S", d, "pat")][0] for d in range(n + 1))
    body.append(
        f"  \\node[font=\\small\\itshape, text=ragSink, anchor=east] "
        f"at ({pos[('S',0,'pat')][0]-0.55:.3f},{pos[('S',0,'pat')][1]:.3f}) "
        f"{{paternal $P_{{\\mathrm{{pat}}}}$}};")
    body.append(
        f"  \\node[font=\\small\\itshape, text=ragSink, anchor=east] "
        f"at ({pos[('S',0,'mat')][0]-0.55:.3f},{pos[('S',0,'mat')][1]:.3f}) "
        f"{{maternal $P_{{\\mathrm{{mat}}}}$}};")
    body.append(
        f"  \\node[font=\\scriptsize, text=ragSink, anchor=south] "
        f"at ({pos['SRC'][0]:.3f},{pos['SRC'][1]+0.42:.3f}) {{immigration}};")
    body.append(
        f"  \\node[font=\\scriptsize, text=ragSink, anchor=south] "
        f"at ({pos['DEATH'][0]:.3f},{pos['DEATH'][1]+0.38:.3f}) {{death}};")
    body.append(
        f"  \\node[font=\\scriptsize, text=ragSink, anchor=north] "
        f"at ({pos['EMIG'][0]:.3f},{pos['EMIG'][1]-0.38:.3f}) {{emigration}};")
    # cluster-depth labels under each depth column
    for d in range(n + 1):
        xm = (pos[("S", d, "pat")][0] + pos[("S", d, "mat")][0]) / 2
        ym = min(pos[("S", d, "pat")][1], pos[("S", d, "mat")][1])
        body.append(
            f"  \\node[font=\\scriptsize, text=ragSink] at ({xm:.3f},{ym-0.62:.3f}) "
            f"{{depth ${d}$}};")

    legend = _legend([
        ("->, ragImm, very thick", "immigration $\\sigma I(t)$, $(1-\\sigma)I(t)$"),
        ("->, ragSink!70, thin", "death $\\mu$ to sink"),
        ("->, ragSink!70, thin, dotted", "emigration $\\varepsilon$ to sink"),
        ("->, ragMate!75, thin, dashed", "mating: parents catalytic (not consumed)"),
        ("->, ragShallowest, semithick", "child edge, shallowest inheritance"),
        ("->, ragDeepest, semithick", "child edge, deepest inheritance"),
        ("->, ragBoth, semithick", "child edge, both rules agree"),
    ], 0.0, -h - 0.55)
    body.extend(legend)
    # node-shape key
    body.append(
        f"  \\node[font=\\scriptsize, anchor=north west] at (8.2,{-h-0.55:.3f}) "
        f"{{\\tikz\\node[ragState, rectangle, fill=ragPatFill, minimum width=0.5cm, "
        f"minimum height=0.34cm] {{$C_d$}}; paternal $\\;\\;$"
        f"\\tikz\\node[ragState, ellipse, fill=ragMatFill, minimum width=0.5cm, "
        f"minimum height=0.34cm] {{$C_d$}}; maternal $\\;\\;$"
        f"\\tikz\\node[ragRxn] {{}}; mating-pair vertex $M_{{qr}}$}};")

    tex = head + "\n".join(body) + "\n\\end{tikzpicture}\n"
    with open(path, "w") as f:
        f.write(tex)
    return G


def write_religious(path):
    n = 2
    G = build_religious(n)
    pos, w, h = layout_religious(G, n)
    node_cm = 0.44

    groupcolors = "".join(
        f"\\definecolor{{ragG{gi}}}{{RGB}}{{{rgb[0]},{rgb[1]},{rgb[2]}}}\n"
        for gi, g in enumerate(R.GROUPS)
        for rgb in [GROUP_COLORS[g][1]])
    extra = ("% Cross-group mating: one faint representative edge per ordered\n"
             "% group pair (g1 -> g2), g1 != g2: 56 edges drawn from a fixed-seed\n"
             "% sampled (q,r); the full cross-group mating set has "
             f"{len(R.GROUPS)*(len(R.GROUPS)-1)*(n+1)**2} combos and is not drawn.\n"
             "% Disaffiliation: all 24 transfers (d,sex,g)->(d,sex,una),\n"
             "% g in {pro,cat,mor,oth}, drawn.\n")
    head = TIKZ_HEAD.format(
        title=f"Religious RAG, n={n} (all vertices; within-group mating all drawn)",
        seed=SEED, extra_comment=extra, groupcolors=groupcolors)

    body = []
    body.extend(_group_backgrounds(G, pos))
    body.extend(_emit_nodes_religious(G, pos, node_cm))
    body.extend(_emit_edges(G, imm_style="imm_thin"))

    body.append(
        f"  \\node[font=\\scriptsize, text=ragSink, anchor=south] "
        f"at ({pos['SRC'][0]:.3f},{pos['SRC'][1]+0.42:.3f}) {{immigration}};")

    legend = _legend([
        ("->, ragImm, very thick", "immigration into $C_0$ (per group, per sex)"),
        ("->, ragSink!70, thin", "death $\\mu$"),
        ("->, ragSink!70, thin, dotted", "emigration $\\varepsilon$"),
        ("->, ragMate!75, thin, dashed", "within-group mating (catalytic)"),
        ("->, ragShallowest, semithick", "child, shallowest"),
        ("->, ragDeepest, semithick", "child, deepest"),
        ("->, ragBoth, semithick", "child, rules agree"),
        ("->, ragDis, semithick, dotted", "disaffiliation $\\delta(t)\\to$ unaffiliated"),
        ("->, ragX!60, ultra thin, dashed", "cross-group mating (1 sampled edge/pair)"),
    ], 0.0, -h - 0.75)
    body.extend(legend)
    body.append(
        f"  \\node[font=\\scriptsize, anchor=north west] at (9.0,{-h-0.75:.3f}) "
        f"{{\\tikz\\node[ragState, rectangle, fill=ragPatFill, minimum width=0.42cm, "
        f"minimum height=0.3cm] {{$C_d$}}; paternal $\\;\\;$"
        f"\\tikz\\node[ragState, ellipse, fill=ragMatFill, minimum width=0.42cm, "
        f"minimum height=0.3cm] {{$C_d$}}; maternal $\\;\\;$"
        f"\\tikz\\node[ragRxn] {{}}; $M_{{qr}}$}};")

    tex = head + "\n".join(body) + "\n\\end{tikzpicture}\n"
    with open(path, "w") as f:
        f.write(tex)
    return G


def edge_report(G, name):
    from collections import Counter
    kinds = Counter()
    for u, v, dat in G.edges(data=True):
        k = dat.get("kind")
        if k == "child":
            kinds[f"child[{dat['rule']}]"] += 1
        else:
            kinds[k] += 1
    nn = Counter(dat.get("kind") for _, dat in G.nodes(data=True))
    print(f"--- {name}: {G.number_of_nodes()} vertices, "
          f"{G.number_of_edges()} directed edges")
    print("    vertices:", dict(nn))
    print("    edges:", dict(kinds))


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    p1 = os.path.join(FIGDIR, "rag_twosex.tex")
    p2 = os.path.join(FIGDIR, "rag_religious.tex")
    G1 = write_twosex(p1)
    print("wrote", p1)
    G2 = write_religious(p2)
    print("wrote", p2)
    edge_report(G1, "two-sex n=4")
    edge_report(G2, "religious n=2")


if __name__ == "__main__":
    main()
