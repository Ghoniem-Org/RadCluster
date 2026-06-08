"""
core/edge_classes.py — The ten abstract reaction-edge classes (Layer 1).

Host-independent.  Implements the elementary edge-class catalogue of
Ghoniem (2026), *A Generalized Graph-Based Cluster Dynamics Framework*,
Section 2.2.

Every edge of a reaction admissibility graph (RAG) falls into exactly one
of ten classes, identified by the pattern of polarity and population at
its endpoints and by the form of its stoichiometric tuple.  This module
defines the class taxonomy and, for each class, the *abstract
stoichiometric contract* — how many precursor vertices it has, whether it
crosses the polarity layers, and how it shifts size / population /
composition.  The concrete, size-resolved rate kernels are supplied by a
host material (Layer 2); the abstract core only fixes the contract that
every host edge of a given class must obey.

The graph-walker assembles the master-equation right-hand side

    dc_a/dt = G_a + sum_r S_ar J_r(c) - D_a c_a

by iterating edges; each edge contributes its signed flux S_ar J_r to the
affected vertices.  Because the contract below fixes the sign pattern of
S_ar per class, conservation of signed lattice-defect content is a
structural property of the graph, not something re-derived per host.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class EdgeClass(Enum):
    """The ten elementary reaction-edge classes (paper Section 2.2)."""

    GROWTH = "growth"                      # 1.  chi-monomer absorption, n -> n+1
    SHRINKAGE = "shrinkage"                # 2.  (-chi)-monomer impingement, n -> n-1
    DISSOCIATION = "dissociation"          # 3.  thermal chi-monomer emission, n -> n-1
    RECOMBINATION = "recombination"        # 4.  cross-polarity monomer annihilation
    ANNIHILATION = "annihilation"          # 5.  cross-polarity binary hyperedge
    INTER_POPULATION = "inter_population"  # 6.  p -> p', size and composition fixed
    SOLUTE_TRAPPING = "solute_trapping"    # 7.  c -> c +/- e_s (trap / detrap)
    COALESCENCE = "coalescence"            # 8.  same-polarity binary growth
    SOURCE = "source"                      # 9.  cascade injection, no precursor
    SINK = "sink"                          # 10. absorption at a fixed sink

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True)
class EdgeClassSpec:
    """The abstract stoichiometric contract obeyed by every edge of a class.

    Attributes
    ----------
    n_precursor_clusters : int
        Number of distinct *cluster* vertices consumed (0 for SOURCE).
    kinetic_order : int
        Order of the reaction flux J_r in concentration: 0 for a pure
        source, 1 for a unary (first-order) transition, 2 for a binary
        collision.  SOLUTE_TRAPPING spans both directions — trapping is
        binary (cavity x gas reservoir), detrapping is unary — and is
        recorded here as order 2 for the forward (trapping) edge.
    consumes_monomer : bool
        True if the reaction consumes a same-polarity point-defect monomer
        from its pool (GROWTH) — distinct from a tracked cluster vertex.
    cross_polarity : bool
        True if the edge couples the chi = -1 and chi = +1 layers.
    size_shift : int | None
        Change in size n of the (surviving) cluster: +1, -1, 0, or None
        when it is not fixed by the class (ANNIHILATION, COALESCENCE,
        SOURCE).
    changes_population : bool
        True if the edge moves the cluster between populations.
    changes_composition : bool
        True if the edge changes the trapped-solute composition vector.
    creates_vertex : bool
        True if the edge has no precursor cluster vertex (SOURCE).
    destroys_to_sink : bool
        True if the edge removes content to an untracked fixed sink (SINK).
    description : str
        One-line statement of the elementary transformation.
    """

    edge_class: EdgeClass
    n_precursor_clusters: int
    kinetic_order: int
    consumes_monomer: bool
    cross_polarity: bool
    size_shift: "int | None"
    changes_population: bool
    changes_composition: bool
    creates_vertex: bool
    destroys_to_sink: bool
    description: str


# ── The fixed, host-independent contract for each of the ten classes ──────────
EDGE_CLASS_SPEC: "dict[EdgeClass, EdgeClassSpec]" = {
    EdgeClass.GROWTH: EdgeClassSpec(
        EdgeClass.GROWTH, 1, 2, True, False, +1, False, False, False, False,
        "(chi,n,p,c) -> (chi,n+1,p,c) by absorption of a chi-monomer"),
    EdgeClass.SHRINKAGE: EdgeClassSpec(
        EdgeClass.SHRINKAGE, 1, 2, False, True, -1, False, False, False, False,
        "(chi,n,p,c) -> (chi,n-1,p,c) by impingement of a (-chi)-monomer"),
    EdgeClass.DISSOCIATION: EdgeClassSpec(
        EdgeClass.DISSOCIATION, 1, 1, False, False, -1, False, False, False, False,
        "(chi,n,p,c) -> (chi,n-1,p,c) by thermal emission of a chi-monomer "
        "into the monomer pool; kernel involves a binding energy"),
    EdgeClass.RECOMBINATION: EdgeClassSpec(
        EdgeClass.RECOMBINATION, 1, 2, False, True, -1, False, False, False, False,
        "(chi,n,p,c) + (-chi,1) -> (chi,n-1,p,c); the (-chi)-monomer is "
        "annihilated rather than returned to any pool"),
    EdgeClass.ANNIHILATION: EdgeClassSpec(
        EdgeClass.ANNIHILATION, 2, 2, True, True, None, False, False, False, False,
        "(chi,n,p,c) + (-chi,n') -> (chi',|n-n'|,p_surv,c_surv); the "
        "opposite-polarity analogue of coalescence (binary hyperedge)"),
    EdgeClass.INTER_POPULATION: EdgeClassSpec(
        EdgeClass.INTER_POPULATION, 1, 1, False, False, 0, True, False, False, False,
        "(chi,n,p,c) -> (chi,n,p',c) with p' != p: segregation, sweeping, "
        "or capture by a moving sink; changes only the population coordinate"),
    EdgeClass.SOLUTE_TRAPPING: EdgeClassSpec(
        EdgeClass.SOLUTE_TRAPPING, 1, 2, False, False, 0, False, True, False, False,
        "(chi,n,p,c) <-> (chi,n,p,c +/- e_s): trapping (binary, cluster x gas "
        "reservoir) or detrapping (unary); changes only the composition vector"),
    EdgeClass.COALESCENCE: EdgeClassSpec(
        EdgeClass.COALESCENCE, 2, 2, False, False, None, False, False, False, False,
        "(chi,n,p,c) + (chi,n',p,c') -> (chi,n+n',p,c+c'): same-polarity "
        "binary growth (hyperedge with a quadratic rate kernel)"),
    EdgeClass.SOURCE: EdgeClassSpec(
        EdgeClass.SOURCE, 0, 0, False, False, None, False, False, True, False,
        "(empty) -> (chi,n,p,c) at a rate set by the displacement-cascade "
        "or transmutation spectrum; the only edge with no precursor vertex"),
    EdgeClass.SINK: EdgeClassSpec(
        EdgeClass.SINK, 1, 1, False, False, None, False, False, False, True,
        "(chi,n,p,c) -> (sink) at a fixed sink (free surface, grain "
        "boundary, dislocation); first-order loss, absorbs both polarities"),
}

# Sanity: the catalogue is complete and consistent.
assert set(EDGE_CLASS_SPEC) == set(EdgeClass), \
    "EDGE_CLASS_SPEC must cover every EdgeClass exactly once"
for _ec, _spec in EDGE_CLASS_SPEC.items():
    assert _spec.edge_class is _ec, f"EDGE_CLASS_SPEC key/value mismatch for {_ec}"


def spec_of(edge_class: EdgeClass) -> EdgeClassSpec:
    """Return the abstract stoichiometric contract for ``edge_class``."""
    return EDGE_CLASS_SPEC[edge_class]
