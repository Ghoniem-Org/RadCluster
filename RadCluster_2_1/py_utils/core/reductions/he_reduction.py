"""
core/reductions/he_reduction.py — Helium-occupancy state-space reductions.

Host-independent declaration of the two helium reductions derived in
Ghoniem (2026), Sections 5.2-5.3, as systematic limits of the same
parent master equations.  Both collapse the prohibitive two-dimensional
(m, ell) vacancy-helium grid; they preserve the helium-conservation
moment exactly.

  CASE1_MEANFIELD  (fusion, fast-He equilibration, Section 5.2)
      Helium migrates far faster than cavities grow, so at each vacancy
      size m the loading equilibrates to a mean value mu(m, t).  The 2-D
      grid c_{m,ell} collapses to two 1-D arrays (c_m, Q_m) with
      mu(m) = Q_m / c_m, reducing the vacancy-He block to ~2 N_V ODEs.

  CASE2_DECOUPLED  (fission, decoupled inventory, Section 5.3)
      Helium is a weak perturbation; the per-size loading need not be
      tracked.  A single scalar trapped-He inventory Q_tot is evolved and
      distributed across cavity sizes algebraically by the capture-rate
      weight m^{1/3} c_m.  The vacancy-He block shrinks to N_V + 1 ODEs.

  DYNAMIC
      No reduction: free helium c_h is an explicit ODE variable.  (A
      quasi-steady-state variant solves dc_h/dt = 0 algebraically.)

This module only fixes the taxonomy and the block-count contract; the
size-resolved helium kinetics are supplied by the host material, which
owns the SOLUTE_TRAPPING edge class.
"""
from __future__ import annotations

from enum import Enum


class HeReductionMode(Enum):
    """The helium-occupancy reduction selected for a run."""

    DYNAMIC = "dynamic"                    # explicit free-He ODE
    CASE1_MEANFIELD = "case1_meanfield"    # fusion: mean loading mu(m) per size
    CASE2_DECOUPLED = "case2_decoupled"    # fission: scalar inventory Q_tot

    @property
    def is_fusion(self) -> bool:
        return self is HeReductionMode.CASE1_MEANFIELD

    @property
    def is_fission(self) -> bool:
        return self is HeReductionMode.CASE2_DECOUPLED

    def he_block_length(self, n_vac_sizes: int, free_he_tracked: bool) -> int:
        """ODE count contributed by the helium state for this reduction.

        Parameters
        ----------
        n_vac_sizes : int
            Number of vacancy size classes (or vacancy bins) carrying He.
        free_he_tracked : bool
            True if free helium c_h is an explicit ODE (dynamic, non-QSS).
        """
        free = 1 if free_he_tracked else 0
        if self is HeReductionMode.CASE1_MEANFIELD:
            return n_vac_sizes + free          # one Q_m per vacancy class
        if self is HeReductionMode.CASE2_DECOUPLED:
            return 1 + free                    # one scalar Q_tot
        return free                            # DYNAMIC: just c_h (if tracked)

    @classmethod
    def from_cascade(cls, cascade: str) -> "HeReductionMode":
        """Map a cascade spectrum to its canonical He reduction.

        fission -> Case 2 (decoupled);  fusion -> Case 1 (mean-field).
        """
        c = str(cascade).lower()
        if "fus" in c:
            return cls.CASE1_MEANFIELD
        if "fis" in c:
            return cls.CASE2_DECOUPLED
        raise ValueError(f"unknown cascade spectrum {cascade!r}")
