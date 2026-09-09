"""
core/to_networkx.py — explicit materialization of a RAG for diagnostics.

Host-independent.  The production path keeps the reaction admissibility
graph *implicit*: an :class:`Edge` is a size-parameterized family, and the
graph walker never instantiates a vertex.  That is the right choice for the
inner loop (see the paper, Section "Implementation"), but it means the
standard structural queries — connectivity, reachability from the monomer
pools, orphaned vertices, cycle structure — are not available for free.

This module supplies them by exporting the declared graph, truncated at a
chosen size cut-off, into an explicit ``networkx.MultiDiGraph``.  It is a
DIAGNOSTIC path only: nothing in the solver imports it, and NetworkX is an
optional dependency.

Typical use::

    from py_utils.materials import build_eurofer_rag
    from py_utils.core.to_networkx import to_networkx, structural_report

    rag, layout = build_eurofer_rag(inp, rr)
    G = to_networkx(rag, n_max=40)
    print(structural_report(rag, n_max=40))

Cite as: NetworkX, Hagberg, Schult & Swart, Proc. SciPy 2008.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple

from .cluster_identifier import Polarity, Population
from .edge_classes import EDGE_CLASS_SPEC, EdgeClass
from .rag import Edge, ReactionAdmissibilityGraph


# ─────────────────────────────────────────────────────────────────────────
# vertex naming
# ─────────────────────────────────────────────────────────────────────────
def vertex(pop: Population, n: int) -> str:
    """Canonical string key for the vertex (polarity, size, population)."""
    return f"{pop.polarity.label}:{pop.name}:{n}"


def _sizes(pop: Population, n_max: int) -> range:
    return range(pop.n_min, n_max + 1)


# ─────────────────────────────────────────────────────────────────────────
# export
# ─────────────────────────────────────────────────────────────────────────
def to_networkx(rag: ReactionAdmissibilityGraph,
                n_max: int = 40,
                include_binary: bool = True):
    """Materialize ``rag`` up to size ``n_max`` as a ``networkx.MultiDiGraph``.

    Parameters
    ----------
    rag : ReactionAdmissibilityGraph
        The host's declared graph.
    n_max : int
        Size cut-off.  Keep it small (tens): the binary edge classes are
        O(n_max**2) per family, which is exactly why the production path
        does not materialize the graph at all.
    include_binary : bool
        Expand COALESCENCE and ANNIHILATION families.  These are directed
        *hyperedges* (two precursors); NetworkX has no hyperedge type, so
        each is represented by its two tail->head arcs, both tagged
        ``hyperedge=True`` and sharing a ``hyper_id``.  Set False for a
        cleaner drawing.

    Returns
    -------
    networkx.MultiDiGraph
        Node attributes: ``polarity``, ``population``, ``size``.
        Edge attributes: ``edge_class``, ``label``, ``kernel``,
        ``hyperedge``, ``hyper_id``.
    """
    try:
        import networkx as nx
    except ImportError as exc:                       # pragma: no cover
        raise ImportError(
            "to_networkx() needs the optional dependency 'networkx' "
            "(pip install networkx). It is a diagnostic path only; the "
            "solver does not require it.") from exc

    G = nx.MultiDiGraph(name=rag.name)

    # ── vertices ────────────────────────────────────────────────────────
    for pop in rag.populations:
        for n in _sizes(pop, n_max):
            G.add_node(vertex(pop, n),
                       polarity=pop.polarity.label,
                       population=pop.name,
                       size=n)

    # Sources and sinks are represented by two sentinel vertices so that
    # reachability questions have somewhere to start and end.
    G.add_node("SOURCE", polarity="none", population="cascade", size=0)
    G.add_node("SINK", polarity="none", population="unresolved", size=0)

    hyper_id = 0
    for e in rag.edges:
        ec = e.edge_class
        pop = e.population
        attrs = dict(edge_class=ec.value, label=e.label, kernel=e.kernel,
                     hyperedge=False, hyper_id=None)

        if ec is EdgeClass.SOURCE:
            for n in _sizes(pop, n_max):
                G.add_edge("SOURCE", vertex(pop, n), **attrs)

        elif ec is EdgeClass.SINK:
            for n in _sizes(pop, n_max):
                G.add_edge(vertex(pop, n), "SINK", **attrs)

        elif ec in (EdgeClass.GROWTH,):
            for n in _sizes(pop, n_max):
                if n + 1 <= n_max:
                    G.add_edge(vertex(pop, n), vertex(pop, n + 1), **attrs)

        elif ec in (EdgeClass.SHRINKAGE, EdgeClass.DISSOCIATION,
                    EdgeClass.RECOMBINATION):
            for n in _sizes(pop, n_max):
                if n - 1 >= pop.n_min:
                    G.add_edge(vertex(pop, n), vertex(pop, n - 1), **attrs)

        elif ec is EdgeClass.INTER_POPULATION:
            dst = e.product_population
            for n in _sizes(pop, n_max):
                if dst.n_min <= n <= n_max:
                    G.add_edge(vertex(pop, n), vertex(dst, n), **attrs)

        elif ec is EdgeClass.SOLUTE_TRAPPING:
            # Composition axis is not expanded here (it would multiply the
            # node count by |C|); the trapping edge is drawn as a self-loop
            # tagged with the gas species index.
            a = dict(attrs, gas_species=e.gas_species)
            for n in _sizes(pop, n_max):
                G.add_edge(vertex(pop, n), vertex(pop, n), **a)

        elif ec is EdgeClass.COALESCENCE and include_binary:
            dst = e.product_population or pop
            partner = e.partner_population or pop
            for n in _sizes(pop, n_max):
                for m in _sizes(partner, n_max):
                    if n + m > n_max or n + m < dst.n_min:
                        continue
                    a = dict(attrs, hyperedge=True, hyper_id=hyper_id)
                    G.add_edge(vertex(pop, n), vertex(dst, n + m), **a)
                    G.add_edge(vertex(partner, m), vertex(dst, n + m), **a)
                    hyper_id += 1

        elif ec is EdgeClass.ANNIHILATION and include_binary:
            partner = e.partner_population
            for n in _sizes(pop, n_max):
                for m in _sizes(partner, n_max):
                    if n == m:
                        head = "SINK"                # mutual annihilation
                    elif n > m:
                        head = vertex(pop, n - m)
                    else:
                        head = vertex(partner, m - n)
                    if head not in G:
                        continue
                    a = dict(attrs, hyperedge=True, hyper_id=hyper_id)
                    G.add_edge(vertex(pop, n), head, **a)
                    G.add_edge(vertex(partner, m), head, **a)
                    hyper_id += 1

    return G


# ─────────────────────────────────────────────────────────────────────────
# serialization
# ─────────────────────────────────────────────────────────────────────────
def write_graphml(G, path) -> None:
    """Write ``G`` to GraphML, dropping attributes GraphML cannot represent.

    ``to_networkx`` tags every arc with ``hyper_id``, which is ``None`` on
    the unary edge classes because they are not hyperedge members, and
    ``kernel`` is ``None`` for any edge the host declared without one.
    The GraphML writer rejects ``None`` outright ("GraphML writer does not
    support <class 'NoneType'> as data values"), so exporting the graph
    unmodified always fails.  Attributes are dropped rather than coerced:
    an absent ``hyper_id`` is what "not a hyperedge" means, and writing a
    sentinel would invite a reader to treat it as an id.

    The input graph is not modified.
    """
    import networkx as nx

    H = G.copy()
    for _, _, d in H.edges(data=True):
        for k in [k for k, v in d.items() if v is None]:
            del d[k]
    for _, d in H.nodes(data=True):
        for k in [k for k, v in d.items() if v is None]:
            del d[k]
    nx.write_graphml(H, path)


# ─────────────────────────────────────────────────────────────────────────
# structural diagnostics
# ─────────────────────────────────────────────────────────────────────────
def structural_report(rag: ReactionAdmissibilityGraph,
                      n_max: int = 40) -> str:
    """Human-readable structural audit of a declared RAG.

    Reports vertex/edge counts, weak connectivity, the set of vertices not
    reachable from the monomer pools (a declaration bug: population with no
    production channel), and the set with no outgoing edge (a sink-only
    trap).  Intended to be run once per host declaration and quoted in the
    paper, not called during integration.
    """
    import networkx as nx

    G = to_networkx(rag, n_max=n_max)
    lines: List[str] = []
    lines.append(f"RAG {rag.name!r} materialized at n_max={n_max}")
    lines.append(f"  vertices : {G.number_of_nodes()}")
    lines.append(f"  arcs     : {G.number_of_edges()}")
    lines.append(f"  edge families declared: {len(rag.edges)}")
    lines.append(f"  weakly connected: {nx.is_weakly_connected(G)}")

    # reachability from the monomer pools
    roots = []
    for polarity in (Polarity.SIA, Polarity.VACANCY):
        try:
            pop = rag.monomer_population(polarity)
        except KeyError:
            continue
        roots.append(vertex(pop, pop.n_min))
    reach = set()
    for r in roots:
        if r in G:
            reach |= nx.descendants(G, r) | {r}
    orphans = [v for v in G if v not in reach and v not in ("SOURCE", "SINK")]
    lines.append(f"  monomer roots: {roots}")
    lines.append(f"  unreachable from monomers: {len(orphans)}")
    if orphans:
        lines.append(f"    e.g. {orphans[:8]}")

    dead = [v for v in G if G.out_degree(v) == 0 and v != "SINK"]
    lines.append(f"  vertices with no outgoing arc: {len(dead)}")
    if dead:
        lines.append(f"    e.g. {dead[:8]}")

    by_class: Dict[str, int] = {}
    for _, _, d in G.edges(data=True):
        by_class[d["edge_class"]] = by_class.get(d["edge_class"], 0) + 1
    lines.append("  arcs by edge class:")
    for k in sorted(by_class, key=lambda k: -by_class[k]):
        lines.append(f"    {k:<18s} {by_class[k]}")
    return "\n".join(lines)


if __name__ == "__main__":                            # pragma: no cover
    import argparse
    import sys
    from pathlib import Path

    # The relative imports above mean this file cannot be run as a script;
    # invoke it as a module from the RadCluster_2_1 directory:
    #     python -m py_utils.core.to_networkx
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from py_utils.input_data import InputData
    from py_utils.reaction_rates import ReactionRates
    from py_utils.materials import build_eurofer_rag

    ap = argparse.ArgumentParser(
        description="Materialize and audit the EUROFER-97 reaction "
                    "admissibility graph.")
    ap.add_argument("--n-max", type=int, default=30,
                    help="size cut-off; binary classes are O(n_max**2)")
    ap.add_argument("--I", type=int, default=200)
    ap.add_argument("--V", type=int, default=200)
    ap.add_argument("--equations", default="discrete",
                    choices=["discrete", "bin_moment"])
    ap.add_argument("--cascade", default="fission",
                    choices=["fission", "fusion"])
    ap.add_argument("--no-binary", action="store_true",
                    help="omit COALESCENCE/ANNIHILATION; they dominate the "
                         "arc count and make any drawing unreadable")
    ap.add_argument("--graphml", type=Path, default=None,
                    help="also export the materialized graph to this path")
    ap.add_argument("--report", type=Path, default=None,
                    help="write the structural report here as well as stdout")
    a = ap.parse_args()

    inp = InputData(I=a.I, V=a.V,
                    physics_option=f"full_CD_{a.cascade}")
    rag, _ = build_eurofer_rag(inp, ReactionRates(inp),
                               equations=a.equations, cascade=a.cascade)

    report = structural_report(rag, n_max=a.n_max)
    print(report)
    if a.report:
        a.report.write_text(report + "\n", encoding="utf-8")
        print(f"\nreport  -> {a.report}")
    if a.graphml:
        G = to_networkx(rag, n_max=a.n_max,
                        include_binary=not a.no_binary)
        write_graphml(G, a.graphml)
        print(f"graphml -> {a.graphml}  "
              f"({G.number_of_nodes()} nodes, {G.number_of_edges()} arcs)")
