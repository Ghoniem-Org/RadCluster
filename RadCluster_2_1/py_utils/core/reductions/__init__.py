"""
core/reductions — Graph-theoretic state-space reductions (Layer 1).

The reductions in this package are stated and implemented on the abstract
reaction admissibility graph, so any host inherits them unchanged
(Ghoniem 2026, Sections 5.2-5.5).

  * ``bin_moment``   — the logarithmic size-bin moment reduction: the
                       per-size SIA / vacancy ladders are replaced by
                       O(log N) bins carrying 1-3 moments each, evolved
                       by the "reconstruct -> walk -> project" procedure
                       wrapped around any :class:`GraphWalker`.

  * ``he_reduction`` — the helium-occupancy reductions: the fast-He
                       mean-field loading mu(m) (Case 1, fusion) and the
                       decoupled scalar inventory Q_tot (Case 2, fission).

Both carry no host parameters; a host selects a reduction and supplies
its kernels.
"""
from .bin_moment import BinMomentReduction
from .he_reduction import HeReductionMode

__all__ = ["BinMomentReduction", "HeReductionMode"]
