"""
materials/eurofer97 — the EUROFER-97 instantiation of the CD core (Layer 2).

EUROFER-97 is the reduced-activation ferritic-martensitic steel used as the
principal worked example of Ghoniem (2026), *A Generalized Graph-Based
Cluster Dynamics Framework for Irradiated Materials*, Section 4.

This package declares the host-specific reaction admissibility graph
G_Eur — its populations, gas list, edge set, and rate-kernel library —
and the contiguous ODE-index :class:`StateLayout` that packs it.  It does
not re-implement the master-equation integrator; the abstract graph
walker (and the production C++ solver) consume the graph unchanged.

Public API
----------
build_eurofer_rag(input_data, reaction_rates, *, equations, cascade)
    Construct and return the ``(ReactionAdmissibilityGraph, StateLayout)``
    pair for an EUROFER-97 run.
"""
from .declaration import build_eurofer_rag

__all__ = ["build_eurofer_rag"]
