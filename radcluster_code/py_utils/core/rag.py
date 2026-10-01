"""
core/rag.py — The Reaction Admissibility Graph  G = (V, E, nu, k)  (Layer 1).

Host-independent.  Implements the labelled multidigraph of
Ghoniem (2026), *A Generalized Graph-Based Cluster Dynamics Framework*,
Sections 2.2 and 3.2:

    G = (V, E, nu, k)

  V : admissible vertex set — cluster identifiers (chi, n, p, c)
  E : edge set — every admissible elementary reaction, each carrying a
      stoichiometric tuple nu and an :class:`EdgeClass`
  nu: the stoichiometric action (fixed per EdgeClass by EDGE_CLASS_SPEC)
  k : the rate-kernel library — calibration-specific, supplied by a host

The graph is *host-specific in (V, E)* and *calibration-specific in k*;
the algorithmic core (the graph walker) is host-independent.  This module
provides only the container and the admissibility logic — it declares no
populations, no edges and no kernels itself.  A host material builds a
RAG by declaring its populations, adding :class:`Edge` families, and
registering rate kernels.

Physically forbidden reactions are excluded at construction time: an edge
is admissible only if its endpoints are admissible and its stoichiometry
respects the host conservation laws.  Once built, the master equation is
the host-independent graph-walking accumulator applied edge-by-edge.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from .cluster_identifier import Polarity, Population
from .edge_classes import EDGE_CLASS_SPEC, EdgeClass


@dataclass
class Edge:
    """One admissible reaction-edge *family* in the RAG (paper Section 3.2).

    An ``Edge`` is parameterised over the size axis: it represents the
    whole family of elementary reactions of one :class:`EdgeClass` acting
    on one population (or population pair).  The graph walker applies the
    edge's kernel vectorised over n, so a single ``Edge`` stands for
    O(N) scalar reactions.

    Attributes
    ----------
    edge_class : EdgeClass
        Which of the ten elementary classes this family belongs to.
    label : str
        Human-readable identifier, unique within a RAG (diagnostics).
    population : Population
        The population the (surviving) cluster belongs to.
    kernel : str
        Key into the host rate-kernel library ``k``.  The kernel supplies
        the size-resolved rate; the abstract core never evaluates it.
    partner_population : Population, optional
        The second precursor population for binary classes
        (COALESCENCE: same polarity; ANNIHILATION: opposite polarity).
    product_population : Population, optional
        The product population for INTER_POPULATION edges (p' != p).
    gas_species : int, optional
        Index into the host gas list, for SOLUTE_TRAPPING edges.
    meta : dict
        Free-form host annotations (e.g. emission channel, sink kind).
    """

    edge_class: EdgeClass
    label: str
    population: Population
    kernel: str
    partner_population: Optional[Population] = None
    product_population: Optional[Population] = None
    gas_species: Optional[int] = None
    meta: Dict = field(default_factory=dict)

    @property
    def spec(self):
        """The host-independent stoichiometric contract for this edge."""
        return EDGE_CLASS_SPEC[self.edge_class]

    def __post_init__(self) -> None:
        sp = self.spec
        if sp.n_precursor_clusters == 2 and self.partner_population is None:
            # Binary edges default the partner to the primary population.
            self.partner_population = self.population
        if self.edge_class is EdgeClass.INTER_POPULATION:
            if self.product_population is None:
                raise ValueError(
                    f"edge {self.label!r}: INTER_POPULATION requires a "
                    "product_population")
            if self.product_population.polarity is not self.population.polarity:
                raise ValueError(
                    f"edge {self.label!r}: inter-population transfer cannot "
                    "cross the polarity layers")
        if self.edge_class is EdgeClass.ANNIHILATION:
            if (self.partner_population is not None and
                    self.partner_population.polarity is self.population.polarity):
                raise ValueError(
                    f"edge {self.label!r}: ANNIHILATION partner must be of "
                    "opposite polarity")
        if (self.edge_class is EdgeClass.COALESCENCE and
                self.product_population is not None and
                self.product_population.polarity is not self.population.polarity):
            # Optional cross-population product (e.g. <111>+<111> junction ->
            # <100>): same-polarity binary growth that deposits the product
            # into a *different* same-polarity population.  Polarity must be
            # preserved so signed-defect content q = chi*n is conserved.
            raise ValueError(
                f"edge {self.label!r}: COALESCENCE product_population must "
                "share the reactants' polarity")
        if (self.edge_class is EdgeClass.SOLUTE_TRAPPING and
                self.gas_species is None):
            raise ValueError(
                f"edge {self.label!r}: SOLUTE_TRAPPING requires gas_species")

    def __repr__(self) -> str:
        return f"<Edge {self.edge_class.value}:{self.label} @ {self.population.name}>"


class ReactionAdmissibilityGraph:
    """The host-specific RAG: vertex space, edge set, and kernel library.

    A host material populates this object by declaring populations, adding
    :class:`Edge` families, and registering rate kernels.  The abstract
    graph walker then consumes it without any further reference to the
    host parameter set.
    """

    def __init__(self, name: str, gas_species: Optional[List[str]] = None):
        self.name = name
        self.gas_species: List[str] = list(gas_species or [])
        # Population names are unique only *within* a polarity layer — the
        # paper's EUROFER instantiation has both a vacancy 'bulk' and an
        # SIA 'bulk' (Eqs. 42-43) — so the registry is keyed (polarity, name).
        self._populations: Dict[Tuple[Polarity, str], Population] = {}
        self._edges: List[Edge] = []
        self._kernels: Dict[str, Callable] = {}
        # The population whose size-1 vertex is the point-defect monomer
        # pool for each polarity (the chi-monomer / (-chi)-monomer source
        # consumed by GROWTH / SHRINKAGE / RECOMBINATION edges).
        self._monomer_pop: Dict[Polarity, Population] = {}

    # ── populations ───────────────────────────────────────────────────────
    @staticmethod
    def _pop_key(p: Population) -> "Tuple[Polarity, str]":
        return (p.polarity, p.name)

    def add_population(self, p: Population) -> Population:
        key = self._pop_key(p)
        if key in self._populations:
            raise ValueError(
                f"duplicate population {p.name!r} in the "
                f"{p.polarity.label} layer")
        self._populations[key] = p
        return p

    @property
    def populations(self) -> Tuple[Population, ...]:
        return tuple(self._populations.values())

    def populations_of(self, polarity: Polarity) -> Tuple[Population, ...]:
        return tuple(p for p in self._populations.values()
                     if p.polarity is polarity)

    def set_monomer_population(self, p: Population) -> None:
        """Declare ``p`` as the point-defect monomer pool for its polarity."""
        if self._pop_key(p) not in self._populations:
            raise ValueError(
                f"monomer population {p.name!r} must be registered first")
        self._monomer_pop[p.polarity] = p

    def monomer_population(self, polarity: Polarity) -> Population:
        """Return the population that holds the ``polarity`` monomer pool."""
        try:
            return self._monomer_pop[polarity]
        except KeyError:
            raise KeyError(
                f"RAG {self.name!r}: no monomer population set for "
                f"{polarity.label}") from None

    # ── edges ─────────────────────────────────────────────────────────────
    def add_edge(self, edge: Edge) -> Edge:
        if any(e.label == edge.label for e in self._edges):
            raise ValueError(f"duplicate edge label {edge.label!r}")
        for p in (edge.population, edge.partner_population,
                  edge.product_population):
            if p is not None and self._pop_key(p) not in self._populations:
                raise ValueError(
                    f"edge {edge.label!r} references unregistered "
                    f"{p.polarity.label} population {p.name!r}")
        self._edges.append(edge)
        return edge

    @property
    def edges(self) -> Tuple[Edge, ...]:
        return tuple(self._edges)

    def edges_of_class(self, edge_class: EdgeClass) -> Tuple[Edge, ...]:
        return tuple(e for e in self._edges if e.edge_class is edge_class)

    # ── kernel library k ──────────────────────────────────────────────────
    def register_kernel(self, name: str, kernel) -> None:
        """Bind a rate kernel to a name.

        ``kernel`` may be a precomputed array/scalar, or a zero-argument
        callable that builds and returns the kernel on first use.  Lazy
        callables let large O(N²) kernels (coalescence/annihilation, the
        junction-branching split) be registered cheaply and materialised only
        if the Python GraphWalker actually evaluates the edge; the C++ solver
        path computes them on the fly and never requests them, so production
        runs never pay the allocation.
        """
        self._kernels[name] = kernel

    def kernel(self, name: str):
        try:
            k = self._kernels[name]
        except KeyError:
            raise KeyError(
                f"RAG {self.name!r}: no kernel registered for {name!r}") from None
        # Resolve and memoise a lazy builder on first request.  Arrays and
        # scalars (the common case) are not callable and pass straight through.
        if callable(k):
            k = k()
            self._kernels[name] = k
        return k

    # ── validation ────────────────────────────────────────────────────────
    def validate(self) -> None:
        """Check structural consistency: every edge's kernel is registered.

        Signed-defect conservation is guaranteed structurally by
        EDGE_CLASS_SPEC (the per-class S_ar sign pattern), so it does not
        need to be re-verified here; only the host-supplied bindings are
        checked.
        """
        for e in self._edges:
            if e.kernel not in self._kernels:
                raise KeyError(
                    f"RAG {self.name!r}: edge {e.label!r} needs unregistered "
                    f"kernel {e.kernel!r}")

    def summary(self) -> str:
        by_class = {ec: len(self.edges_of_class(ec)) for ec in EdgeClass}
        active = {k.value: v for k, v in by_class.items() if v}
        return (f"RAG {self.name!r}: {len(self._populations)} populations, "
                f"{len(self._edges)} edges, {len(self._kernels)} kernels; "
                f"edge classes active: {active}")

    def __repr__(self) -> str:
        return (f"<ReactionAdmissibilityGraph {self.name!r} "
                f"V_pop={len(self._populations)} E={len(self._edges)}>")
