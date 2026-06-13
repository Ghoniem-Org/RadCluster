"""
RadCluster_2_1 — Layer 1: abstract, host-independent cluster-dynamics core.

This package implements the host-independent algorithmic core of the
generalized graph-based cluster-dynamics framework of

    Ghoniem, N.M. (2026), "A Generalized Graph-Based Cluster Dynamics
    Framework for Irradiated Materials"
    (docs/Formulation/Generalized_Cluster_Dynamics.pdf).

It defines the two algorithmic primitives that make a generic CD
formulation possible:

  * the polarity index  chi in {-1, +1}  that splits the cluster state
    space into vacancy-type and SIA-type sub-spaces without any
    host-specific taxonomy (``cluster_identifier``); and

  * the reaction admissibility graph (RAG)  G = (V, E, nu, k)  whose
    edges encode every elementary reaction with stoichiometric tuples
    and rate kernels (``rag``, ``edge_classes``).

The master equation reduces to the host-independent graph-walking
accumulator

    dc_a/dt = G_a + sum_r S_ar J_r(c) - D_a c_a

applied edge-by-edge over E (``graph_walker``).

Nothing in this package references a specific host material.  A host
(e.g. EUROFER-97) is introduced as an *instantiation* under
``py_utils/materials/`` — a declaration of the admissible vertex set,
population set, edge set, and rate-kernel library — never a rewrite of
this core.  This is the abstract-base-class / subclass split described
in the paper's introduction.
"""
from .cluster_identifier import (
    Polarity,
    Population,
    ClusterIdentifier,
)
from .edge_classes import (
    EdgeClass,
    EdgeClassSpec,
    EDGE_CLASS_SPEC,
    spec_of,
)
from .rag import (
    Edge,
    ReactionAdmissibilityGraph,
)
from .state_layout import (
    Block,
    StateLayout,
)
from .graph_walker import (
    GraphWalker,
)
from .reductions import (
    BinMomentReduction,
    HeReductionMode,
)

__all__ = [
    "Polarity",
    "Population",
    "ClusterIdentifier",
    "EdgeClass",
    "EdgeClassSpec",
    "EDGE_CLASS_SPEC",
    "spec_of",
    "Edge",
    "ReactionAdmissibilityGraph",
    "Block",
    "StateLayout",
    "GraphWalker",
    "BinMomentReduction",
    "HeReductionMode",
]
