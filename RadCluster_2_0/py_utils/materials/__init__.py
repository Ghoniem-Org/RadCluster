"""
materials — Layer 2: host-material instantiations of the abstract CD core.

A *host material* (e.g. EUROFER-97) is introduced under this package as an
*instantiation* of the host-independent core in ``py_utils/core`` — a
declaration of the admissible vertex set, population set, edge set, and
rate-kernel library of Ghoniem (2026), *A Generalized Graph-Based Cluster
Dynamics Framework for Irradiated Materials*.

This is the "subclass" half of the paper's abstract-base-class / subclass
split: Layer 1 (``core``) fixes the algorithm; Layer 2 (``materials``)
declares *which* populations, edges, and kernels a particular alloy
admits.  A material instantiation never re-implements the integrator —
the C++ solver remains the engine — it only builds the
:class:`~py_utils.core.rag.ReactionAdmissibilityGraph` and the
:class:`~py_utils.core.state_layout.StateLayout` that describe the host.

Sub-packages
------------
eurofer97/ — the EUROFER-97 ferritic-martensitic steel instantiation
             (paper Section 4): one SIA bulk population, one vacancy
             bulk population, helium as the single resolved gas, and the
             eight EUROFER process classes P1-P8.
"""
from .eurofer97 import build_eurofer_rag

__all__ = ["build_eurofer_rag"]
