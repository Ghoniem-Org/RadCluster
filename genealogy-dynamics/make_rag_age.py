#!/usr/bin/env python3
r"""Build the NetworkX-generated TikZ RAG for the AGE-STRUCTURED model.

Reads the REAL model structure from the genealogy package
(populations.child_cluster for shallowest/deepest, age.AGING_FLUX,
age.REPRO_BANDS, age.IMMIGRANT_AGE_PROFILE) and emits an \input-able
TikZ fragment (drive_doc/figs/rag_age.tex) in the style of rag_twosex.tex.

Representative subset (the full model has 15 clusters x 2 sexes x 18 ages
= 540 states at the n = 14 default, far too many to draw): n=2 pedigree clusters (C_0..C_2),
both sexes, age bands a=0..4 (0-4, 5-9, 10-14, 15-19, 20-24). Every drawn
edge corresponds to a real term in genealogy.age.rhs_age:

  aging:      (d,s,a) -> (d,s,a+1), rate gamma_a = 1/5/yr (AGING_FLUX)
  mate_in:    (q,pat,a),(r,mat,a) -> M_qr for a in REPRO_BANDS (here a=3,4);
              parents catalytic, weighted by beta_a (mating_weights_age)
  child:      M_qr -> (c,pat,0),(c,mat,0), c = child_cluster(q,r,rule)
              (birth_inflow_age; newborns enter band a=0), shallowest vs deepest
  imm:        SRC -> (0,s,a), I*sigma_s*imm_age_a (immigrant age profile)
  death:      (d,s,a) -> DEATH, mu_{s,a}
  emig:       (d,s,a) -> EMIG, eps

Layout: networkx.spring_layout with fixed seed; state vertices pinned on an
age (x) by cluster (y) grid, mating vertices relaxed near the reproductive
columns. A full legend names every edge type; a few representative edges
carry direct labels.
"""

import os

os.environ.setdefault("GENEALOGY_MAX_DEPTH", "12")

import numpy as np
import networkx as nx

import genealogy.populations as P
from genealogy import age as A

SEED = 20260926
FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "drive_doc", "figs")

# Representative subset -------------------------------------------------------
N_SUB = 2                    # clusters C_0..C_2  (model code run at n=2)
AGES_SUB = (0, 1, 2, 3, 4)   # 0-4, 5-9, 10-14, 15-19, 20-24
AGE_LABEL = {0: "0--4", 1: "5--9", 2: "10--14", 3: "15--19", 4: "20--24"}
SEXES = ("pat", "mat")
RULES = ("shallowest_inheritance", "deepest_inheritance")
REPRO_SUB = tuple(a for a in AGES_SUB if a in A.REPRO_BANDS)  # (3, 4)


def _set_depth(n):
    P.MAX_DEPTH = n
    P.N_CLUSTERS = n + 1


def build_age():
    """Age-structured RAG on the representative subset."""
    _set_depth(N_SUB)
    # sanity: the drawn edges must match the real model constants
    assert A.N_AGE == 18 and len(A.AGING_FLUX) == 18
    assert all(A.AGING_FLUX[a] == 1.0 / 5 for a in range(17))
    assert A.AGING_FLUX[17] == 0.0
    assert set(REPRO_SUB) <= set(A.REPRO_BANDS)

    G = nx.DiGraph()
    for d in range(N_SUB + 1):
        for sex in SEXES:
            for a in AGES_SUB:
                G.add_node(("S", d, sex, a), kind="state", d=d, sex=sex, a=a)
    G.add_node("SRC", kind="source")
    G.add_node("DEATH", kind="sink")
    G.add_node("EMIG", kind="sink")

    # aging transfers (adjacent compartments, gamma_a = 1/5/yr)
    for d in range(N_SUB + 1):
        for sex in SEXES:
            for a in AGES_SUB[:-1]:
                G.add_edge(("S", d, sex, a), ("S", d, sex, a + 1),
                           kind="aging", gamma=float(A.AGING_FLUX[a]))

    # mating-pair reaction vertices; only reproductive ages feed them
    for q in range(N_SUB + 1):
        for r in range(N_SUB + 1):
            m = ("M", q, r)
            G.add_node(m, kind="rxn", q=q, r=r)
            for a in REPRO_SUB:
                G.add_edge(("S", q, "pat", a), m, kind="mate_in")
                G.add_edge(("S", r, "mat", a), m, kind="mate_in")
            for rule in RULES:
                c = P.child_cluster(q, r, rule)  # actual model code
                G.add_edge(m, ("S", c, "pat", 0), kind="child", rule=rule)
                G.add_edge(m, ("S", c, "mat", 0), kind="child", rule=rule)

    # immigration into C_0 across the drawn ages (real age profile)
    for sex in SEXES:
        for a in AGES_SUB:
            G.add_edge("SRC", ("S", 0, sex, a), kind="imm",
                       share=float(A.IMMIGRANT_AGE_PROFILE[a]))

    # death + emigration sinks from every drawn compartment
    for d in range(N_SUB + 1):
        for sex in SEXES:
            for a in AGES_SUB:
                G.add_edge(("S", d, sex, a), "DEATH", kind="death")
                G.add_edge(("S", d, sex, a), "EMIG", kind="emig")
    return G


# ---------------------------------------------------------------- layout
def layout_age(G, width_cm=16.0):
    dx, dy, ds = 3.1, 2.1, 0.52
    init, fixed = {}, []
    for d in range(N_SUB + 1):
        for sex in SEXES:
            for a in AGES_SUB:
                k = ("S", d, sex, a)
                init[k] = (a * dx,
                           (N_SUB - d) * dy + (ds if sex == "pat" else -ds))
                fixed.append(k)
    ymid = N_SUB * dy / 2
    init["SRC"] = (-2.6, ymid); fixed.append("SRC")
    init["DEATH"] = (len(AGES_SUB) * dx - dx + 2.2, ymid + dy)
    fixed.append("DEATH")
    init["EMIG"] = (len(AGES_SUB) * dx - dx + 2.2, ymid - dy)
    fixed.append("EMIG")
    rng = np.random.default_rng(SEED)
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") == "rxn":
            q, r = dat["q"], dat["r"]
            init[nd] = (3.4 * dx + rng.normal(0, 0.3),
                        ((N_SUB - q) + (N_SUB - r)) / 2 * dy
                        + rng.normal(0, 0.4))
    clamp = {nd: (-2.6, len(AGES_SUB) * dx, -1.6, N_SUB * dy + 1.6)
             for nd, dat in G.nodes(data=True)
             if dat.get("kind") == "rxn"}
    pos = nx.spring_layout(G, pos=init, fixed=fixed, seed=SEED,
                           iterations=400)
    for k, box in clamp.items():
        if k in pos:
            x, y = pos[k]
            x0, x1, y0, y1 = box
            pos[k] = (min(max(x, x0), x1), min(max(y, y0), y1))
    # pin mating vertices' x near the reproductive columns (their child
    # edges pull them left toward the newborn targets; override so the
    # mating cloud sits over the reproductive ages it feeds from)
    mate_x = 3.55 * dx
    for nd, dat in G.nodes(data=True):
        if dat.get("kind") == "rxn":
            _, y = pos[nd]
            pos[nd] = (mate_x, y)
    xs = np.array([p[0] for p in pos.values()])
    ys = np.array([p[1] for p in pos.values()])
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    sx = width_cm / max(x1 - x0, 1e-9)
    return ({k: ((x - x0) * sx, (y - y0) * sx) for k, (x, y) in pos.items()},
            (x1 - x0) * sx, (y1 - y0) * sx)


# ---------------------------------------------------------------- tikz
def _nid(k):
    return k if isinstance(k, str) else "-".join(str(x) for x in k)


TIKZ_HEAD = """% Age-structured RAG, representative subset: n=2 clusters (C_0..C_2) x 2 sexes x ages 0-4..20-24
% Generated by make_rag_age.py from the genealogy package (NetworkX, seed {seed}).
% Vertices: representative (d,s,a) state compartments + immigration source + death/emigration
% sinks + mating-pair reaction vertices (diamonds). Child cluster from
% genealogy.populations.child_cluster (shallowest / deepest inheritance); aging rates from
% genealogy.age.AGING_FLUX (gamma_a = 1/5/yr); reproductive ages from genealogy.age.REPRO_BANDS;
% immigrant age shares from genealogy.age.IMMIGRANT_AGE_PROFILE.
% The full model has 15 clusters x 2 sexes x 18 ages = 540 states at the
% n = 14 default; only the drawn subset is shown.
% Requires: \\usepackage{{tikz}} \\usetikzlibrary{{arrows.meta,shapes}}
\\definecolor{{ragPatFill}}{{RGB}}{{219,232,251}}
\\definecolor{{ragMatFill}}{{RGB}}{{179,205,245}}
\\definecolor{{ragImm}}{{RGB}}{{0,51,153}}
\\definecolor{{ragSink}}{{RGB}}{{110,110,110}}
\\definecolor{{ragMate}}{{RGB}}{{150,90,30}}
\\definecolor{{ragShallowest}}{{RGB}}{{230,110,20}}
\\definecolor{{ragDeepest}}{{RGB}}{{0,150,130}}
\\definecolor{{ragBoth}}{{RGB}}{{60,60,60}}
\\definecolor{{ragAge}}{{RGB}}{{0,120,200}}
\\begin{{tikzpicture}}[>=Stealth,
  ragState/.style={{draw=ragSink, thick, font=\\tiny\\bfseries, inner sep=1pt}},
  ragRxn/.style={{diamond, draw=ragMate, fill=ragMate!25, inner sep=0pt, minimum size=3.2pt}},
  ragSrc/.style={{circle, draw=ragImm, fill=ragImm!12, thick, font=\\small, minimum size=0.55cm}},
  ragSinkSt/.style={{circle, draw=ragSink, fill=ragSink!12, font=\\scriptsize, minimum size=0.5cm}},
]
"""


EDGE_STYLE = {
    "imm": "->, ragImm, very thick",
    "death": "->, ragSink!70, thin",
    "emig": "->, ragSink!70, thin, dotted",
    "mate_in": "->, ragMate!75, thin, dashed",
    "aging": "->, ragAge, semithick",
    "child_shallowest": "->, ragShallowest, semithick",
    "child_deepest": "->, ragDeepest, semithick",
    "child_both": "->, ragBoth, semithick",
}


def emit(G, pos):
    lines = []
    # --- vertices ---------------------------------------------------------
    for nd, dat in G.nodes(data=True):
        x, y = pos[nd]
        kind, nid = dat.get("kind"), _nid(nd)
        if kind == "state":
            d, sex, a = dat["d"], dat["sex"], dat["a"]
            shape = "rectangle" if sex == "pat" else "ellipse"
            fill = "ragPatFill" if sex == "pat" else "ragMatFill"
            lines.append(
                f"  \\node[ragState, {shape}, fill={fill}, minimum width=0.60cm, "
                f"minimum height=0.42cm] ({nid}) at ({x:.3f},{y:.3f}) "
                f"{{$C_{{{d}}}^{{{a}}}$}};")
        elif kind == "rxn":
            lines.append(f"  \\node[ragRxn] ({nid}) at ({x:.3f},{y:.3f}) {{}};")
        elif kind == "source":
            lines.append(f"  \\node[ragSrc] ({nid}) at ({x:.3f},{y:.3f}) "
                         f"{{$\\varnothing$}};")
        elif kind == "sink":
            lab = "$\\mu$" if nd == "DEATH" else "$\\varepsilon$"
            lines.append(f"  \\node[ragSinkSt] ({nid}) at ({x:.3f},{y:.3f}) "
                         f"{{{lab}}};")
    # --- edges ------------------------------------------------------------
    # merge shallowest/deepest child edges to the same target into one dark edge
    child_rules = {}
    for u, v, dat in G.edges(data=True):
        if dat.get("kind") == "child":
            child_rules.setdefault((u, v), set()).add(dat["rule"])
    seen_both = set()
    # direct labels on one representative edge of each key type
    labeled = set()
    for u, v, dat in G.edges(data=True):
        kind = dat.get("kind")
        if kind == "child":
            rules = child_rules[(u, v)]
            if len(rules) == 2:
                if (u, v) in seen_both:
                    continue
                seen_both.add((u, v))
                style = EDGE_STYLE["child_both"]
                opts = ""
            else:
                rule = dat["rule"]
                style = (EDGE_STYLE["child_shallowest"]
                         if rule == "shallowest_inheritance"
                         else EDGE_STYLE["child_deepest"])
                opts = ""
            lines.append(f"  \\draw[{style}{opts}] ({_nid(u)}) -- ({_nid(v)});")
        elif kind == "aging":
            tag = ""
            if "aging" not in labeled:
                tag = (" node[midway, font=\\tiny, text=ragAge, above] "
                       "{$\\gamma_a{=}1/5$}")
                labeled.add("aging")
            lines.append(f"  \\draw[{EDGE_STYLE['aging']}] ({_nid(u)}) -- "
                         f"({_nid(v)}){tag};")
        elif kind == "imm":
            tag = ""
            if "imm" not in labeled:
                tag = (" node[midway, font=\\tiny, text=ragImm, above, sloped] "
                       "{immigration}")
                labeled.add("imm")
            lines.append(f"  \\draw[{EDGE_STYLE['imm']}] ({_nid(u)}) -- "
                         f"({_nid(v)}){tag};")
        elif kind == "mate_in":
            lines.append(f"  \\draw[{EDGE_STYLE['mate_in']}] ({_nid(u)}) -- "
                         f"({_nid(v)});")
        elif kind in EDGE_STYLE:
            lines.append(f"  \\draw[{EDGE_STYLE[kind]}] ({_nid(u)}) -- "
                         f"({_nid(v)});")
    return lines


def annotations(G, pos, h):
    out = []
    # age-band column labels (top)
    for a in AGES_SUB:
        xs = [pos[("S", d, sex, a)][0]
              for d in range(N_SUB + 1) for sex in SEXES]
        ys = [pos[("S", d, sex, a)][1]
              for d in range(N_SUB + 1) for sex in SEXES]
        out.append(
            f"  \\node[font=\\scriptsize\\bfseries, text=ragSink] "
            f"at ({sum(xs)/len(xs):.3f},{max(ys)+0.85:.3f}) "
            f"{{age ${AGE_LABEL[a]}$}};")
    # cluster row labels (left of the a=0 column)
    for d in range(N_SUB + 1):
        ys = [pos[("S", d, sex, 0)][1] for sex in SEXES]
        out.append(
            f"  \\node[font=\\small\\itshape, text=ragSink, anchor=east] "
            f"at ({pos[('S',d,'pat',0)][0]-0.75:.3f},{sum(ys)/2:.3f}) "
            f"{{$C_{d}$}};")
    # reproductive-age bracket
    xr = [pos[("S", 0, sex, a)][0] for a in REPRO_SUB for sex in SEXES]
    out.append(
        f"  \\node[font=\\tiny\\itshape, text=ragMate] "
        f"at ({sum(xr)/len(xr):.3f},{-0.55:.3f}) "
        f"{{reproductive ages feed mating (bands 15--49)}};")
    out.append(
        f"  \\node[font=\\scriptsize, text=ragSink, anchor=south] "
        f"at ({pos['SRC'][0]:.3f},{pos['SRC'][1]+0.42:.3f}) {{immigration}};")
    out.append(
        f"  \\node[font=\\scriptsize, text=ragSink, anchor=south] "
        f"at ({pos['DEATH'][0]:.3f},{pos['DEATH'][1]+0.38:.3f}) {{death}};")
    out.append(
        f"  \\node[font=\\scriptsize, text=ragSink, anchor=north] "
        f"at ({pos['EMIG'][0]:.3f},{pos['EMIG'][1]-0.38:.3f}) {{emigration}};")
    # legend (every edge type named)
    leg = [
        "  \\node[draw=ragSink!40, fill=white, rounded corners=3pt, "
        "font=\\scriptsize, align=left, anchor=north west] (ragLegend) "
        f"at (0.000,{-h-0.55:.3f}) {{%;",
        "    \\begin{tabular}{@{}l@{\\hspace{4pt}}l@{}}",
        "      \\tikz\\draw[->, ragAge, semithick] (0,0) -- (0.7,0); "
        "& aging transfer $\\gamma_a = 1/5\\,\\mathrm{yr}^{-1}$ (adjacent bands) \\\\",
        "      \\tikz\\draw[->, ragImm, very thick] (0,0) -- (0.7,0); "
        "& immigration $I(t)\\,\\sigma_s\\,\\mathrm{imm\\_age}_a$ into $C_0$ \\\\",
        "      \\tikz\\draw[->, ragMate!75, thin, dashed] (0,0) -- (0.7,0); "
        "& catalytic mating: reproductive-age parents (not consumed) \\\\",
        "      \\tikz\\draw[->, ragShallowest, semithick] (0,0) -- (0.7,0); "
        "& child edge, shallowest inheritance (newborns to age 0--4) \\\\",
        "      \\tikz\\draw[->, ragDeepest, semithick] (0,0) -- (0.7,0); "
        "& child edge, deepest inheritance (newborns to age 0--4) \\\\",
        "      \\tikz\\draw[->, ragBoth, semithick] (0,0) -- (0.7,0); "
        "& child edge, both rules agree \\\\",
        "      \\tikz\\draw[->, ragSink!70, thin] (0,0) -- (0.7,0); "
        "& death $\\mu_{s,a}$ to sink \\\\",
        "      \\tikz\\draw[->, ragSink!70, thin, dotted] (0,0) -- (0.7,0); "
        "& emigration $\\varepsilon$ to sink \\\\",
        "    \\end{tabular}};",
    ]
    out.extend(leg)
    out.append(
        f"  \\node[font=\\scriptsize, anchor=north west] at (8.6,{-h-0.55:.3f}) "
        f"{{\\tikz\\node[ragState, rectangle, fill=ragPatFill, minimum width=0.5cm, "
        f"minimum height=0.34cm] {{$C_d^a$}}; paternal $\\;\\;$"
        f"\\tikz\\node[ragState, ellipse, fill=ragMatFill, minimum width=0.5cm, "
        f"minimum height=0.34cm] {{$C_d^a$}}; maternal $\\;\\;$"
        f"\\tikz\\node[ragRxn] {{}}; mating-pair vertex $M_{{qr}}$}};")
    out.append(
        f"  \\node[font=\\tiny\\itshape, text=ragSink, anchor=north west] "
        f"at (0.000,{-h-2.9:.3f}) {{Representative subset of the 468-state model "
        f"(13 clusters $\\times$ 2 sexes $\\times$ 18 ages); vertex $C_d^a$ = "
        f"cluster $d$, age band $a$.}};")
    return out


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    G = build_age()
    pos, w, h = layout_age(G)
    tex = (TIKZ_HEAD.format(seed=SEED) + "\n".join(emit(G, pos))
           + "\n" + "\n".join(annotations(G, pos, h))
           + "\n\\end{tikzpicture}\n")
    path = os.path.join(FIGDIR, "rag_age.tex")
    with open(path, "w") as f:
        f.write(tex)
    print("wrote", path)

    from collections import Counter
    kinds = Counter()
    for u, v, dat in G.edges(data=True):
        k = dat.get("kind")
        kinds[f"child[{dat['rule']}]" if k == "child" else k] += 1
    nn = Counter(dat.get("kind") for _, dat in G.nodes(data=True))
    print(f"age RAG (n={N_SUB}, ages {AGES_SUB[0]}-{AGES_SUB[-1]}): "
          f"{G.number_of_nodes()} vertices, {G.number_of_edges()} edges")
    print("  vertices:", dict(nn))
    print("  edges:", dict(kinds))


if __name__ == "__main__":
    main()
