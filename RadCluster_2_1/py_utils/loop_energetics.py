"""
loop_energetics.py — ½⟨111⟩ vs ⟨100⟩ interstitial-loop free energies and the
thermodynamic driving force for the *unary* (Dudarev) ½⟨111⟩ → ⟨100⟩ conversion.

Part of the ½⟨111⟩ → ⟨100⟩ loop-conversion work
(see ``docs/design_notes/loop_111_to_100_conversion.md``).

Physics
-------
Prismatic-loop free energy (Dudarev, Bullough & Derlet, PRL 100, 135503 (2008),
Eq. 5; equivalently Marian, Wirth & Perlado, PRL 88, 255507 (2002), Eq. 3):

    E_l^X(n,T) = P_X(n) · [ F̂_X(T)·ln(4 R*_X(n) / (e·δ)) + F_δ^X + F_c^X ]

with X ∈ {111, 100}, P_X the loop perimeter, R*_X the equivalent-circle radius
of a platelet of n self-interstitials (area A = n·Ω/b_X), and δ ≈ 0.4 nm the
core cutoff.  ½⟨111⟩ loops are hexagonal {110} (6 sides, b = (√3/2)a); ⟨100⟩
loops are square {100} (4 sides, b = a).

Temperature enters only through the ⟨100⟩ prelogarithmic factor, which softens
toward the α–γ transition as the analytic ⟨100⟩[100] solution
F̂_001([100]) ∝ √(c₁₁−c₁₂) ∝ (1 − T/T_c)^{1/4} (Dudarev Eq. 4); the
½⟨111⟩ factor is treated as ~T-independent:

    F̂_100(T) = F̂_100^0 · (1 − T/T_c)^{1/4},   T_c = 1185 K.

Driving force for ½⟨111⟩ → ⟨100⟩ (favourable where positive):

    ΔF(n,T) = E_l^111(n,T) − E_l^100(n,T).

Calibration
-----------
Rather than trust digitised absolute values of the elastic constants, the
overall ⟨100⟩ prelog magnitude F̂_100^0 is *derived* from a single calibration
target — the crossover temperature T* at a reference size n_ref where
ΔF(n_ref, T*) = 0 — via :meth:`LoopEnergetics.calibrate`.  This guarantees, by
construction, that at the reference size ½⟨111⟩ is favoured below T* and ⟨100⟩
above it, matching the Dudarev Fig. 4 stability regions (sign change ≈ 350 °C,
strongly ⟨100⟩-favoured by ≈ 550 °C).  The remaining size dependence of the
crossover then *emerges* from the loop geometry.

.. note::
   The default core/prelog constants below are approximate (digitised from
   Dudarev Figs. 2–3 + best-fit table).  The *direction and magnitude* of the
   size dependence of the crossover is sensitive to them and is a physics
   calibration decision (Phase 6, against ``loop_burgers_fraction.py``); only
   the qualitative T behaviour at the reference size is guaranteed by the
   calibration.  All constants are overridable on the dataclass.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_E = float(np.e)
_C_TO_K = 273.15


@dataclass
class LoopEnergetics:
    """½⟨111⟩ / ⟨100⟩ loop energies and the unary-conversion driving force.

    Lengths are in ångström and energies in eV.  Defaults are bcc-Fe values
    consistent with :mod:`binding_energies`.
    """

    # ── lattice / geometry ───────────────────────────────────────────────────
    a: float = 2.8665          # bcc-Fe lattice parameter [Å]
    Omega: float = 11.78       # atomic volume [Å³]  (= a³/2)
    delta: float = 4.0         # dislocation core cutoff [Å]  (0.4 nm)
    Tc: float = 1185.0         # α–γ / spin-fluctuation softening temperature [K]

    # ── prelogarithmic energy factors [eV/Å] ─────────────────────────────────
    Fhat_111: float = 0.28     # ½⟨111⟩, treated ~T-independent (Dudarev Fig. 2)
    Fhat_100_0: float | None = None   # ⟨100⟩ at T=0; set by calibrate() if None

    # ── core-traction + nonlinear-core terms [eV/Å] (Dudarev best fit) ────────
    Fdelta_111: float = 0.345
    Fc_111: float = 0.46
    Fdelta_100: float = 0.387
    Fc_100: float = 0.33

    # ── default calibration target ───────────────────────────────────────────
    T_star_C: float = 450.0    # crossover temperature [°C] at n_ref
    n_ref: float = 50.0        # reference loop size [SIAs]

    def __post_init__(self) -> None:
        self.b_111 = (np.sqrt(3.0) / 2.0) * self.a
        self.b_100 = self.a
        if self.Fhat_100_0 is None:
            self.calibrate(self.T_star_C, self.n_ref)

    # ── geometry ─────────────────────────────────────────────────────────────
    def _b(self, char: int) -> float:
        return self.b_111 if char == 111 else self.b_100

    def area(self, n, char: int):
        """Platelet area A = n·Ω/b_X  [Å²]."""
        return np.asarray(n, dtype=float) * self.Omega / self._b(char)

    def R_star(self, n, char: int):
        """Equivalent-circle radius R*_X = √(A/π)  [Å]."""
        return np.sqrt(self.area(n, char) / np.pi)

    def perimeter(self, n, char: int):
        """Loop perimeter: hexagon (½⟨111⟩) or square (⟨100⟩)  [Å]."""
        A = self.area(n, char)
        if char == 111:
            return 6.0 * np.sqrt(2.0 * A / (3.0 * np.sqrt(3.0)))
        return 4.0 * np.sqrt(A)

    def _Lln(self, n, char: int):
        return np.log(4.0 * self.R_star(n, char) / (_E * self.delta))

    # ── per-character loop energies ──────────────────────────────────────────
    def E_l_111(self, n):
        """½⟨111⟩ loop free energy E_l^111(n)  [eV]  (T-independent here)."""
        P = self.perimeter(n, 111)
        return P * (self.Fhat_111 * self._Lln(n, 111)
                    + self.Fdelta_111 + self.Fc_111)

    def E_l_100_core(self, n):
        """T-independent core part of the ⟨100⟩ loop energy  [eV]."""
        return self.perimeter(n, 100) * (self.Fdelta_100 + self.Fc_100)

    def E_l_100_pre(self, n):
        """⟨100⟩ prelog part *without* the F̂_100(T) factor  [eV·Å/(eV/Å)]."""
        return self.perimeter(n, 100) * self._Lln(n, 100)

    def _softening(self, T_K):
        """(1 − T/T_c)^{1/4}, clipped to ≥ 0 above T_c (Dudarev Eq. 4)."""
        return np.maximum(1.0 - np.asarray(T_K, dtype=float) / self.Tc, 0.0) ** 0.25

    def Fhat_100(self, T_K):
        """⟨100⟩ prelogarithmic factor F̂_100(T)  [eV/Å]."""
        return self.Fhat_100_0 * self._softening(T_K)

    def E_l_100(self, n, T_K):
        """⟨100⟩ loop free energy E_l^100(n,T)  [eV]."""
        return self.E_l_100_core(n) + self.Fhat_100(T_K) * self.E_l_100_pre(n)

    # ── driving force & calibration ──────────────────────────────────────────
    def delta_F(self, n, T_K):
        """Unary-conversion driving force ΔF(n,T) [eV]; > 0 favours ⟨100⟩."""
        return self.E_l_111(n) - self.E_l_100(n, T_K)

    def calibrate(self, T_star_C: float | None = None,
                  n_ref: float | None = None) -> float:
        """Set F̂_100^0 so that ΔF(n_ref, T*) = 0; return F̂_100^0.

        Guarantees ΔF(n_ref, T) < 0 for T < T* and > 0 for T > T* (the
        softening is monotonic), i.e. ½⟨111⟩ stable below the crossover and
        ⟨100⟩ above it.  Requires E_l^111(n_ref) > E_l^100_core(n_ref) (the
        ⟨100⟩ core alone lies below the ½⟨111⟩ energy), which holds for the
        defaults.
        """
        if T_star_C is not None:
            self.T_star_C = float(T_star_C)
        if n_ref is not None:
            self.n_ref = float(n_ref)
        T_star = self.T_star_C + _C_TO_K
        soft = float(self._softening(T_star))
        num = float(self.E_l_111(self.n_ref) - self.E_l_100_core(self.n_ref))
        den = float(self.E_l_100_pre(self.n_ref)) * soft
        if not (num > 0.0 and den > 0.0):
            raise ValueError(
                "calibration ill-posed: need E_l^111(n_ref) > E_l^100_core(n_ref) "
                f"and 0 < T* < T_c (got num={num:.3g}, soft={soft:.3g}); "
                "adjust the core/prelog constants or the target.")
        self.Fhat_100_0 = num / den
        return self.Fhat_100_0

    def crossover_temperature(self, n):
        """Temperature [°C] where ΔF(n,T)=0, or NaN if none in (0, T_c).

        Solves F̂_100^0·(1−T/T_c)^{1/4}·E_pre = E_l^111 − E_l^100_core.
        """
        n = np.asarray(n, dtype=float)
        ratio = ((self.E_l_111(n) - self.E_l_100_core(n))
                 / (self.Fhat_100_0 * self.E_l_100_pre(n)))
        with np.errstate(invalid="ignore"):
            T_K = np.where((ratio > 0.0) & (ratio < 1.0),
                           self.Tc * (1.0 - ratio ** 4), np.nan)
        out = T_K - _C_TO_K
        return float(out) if out.ndim == 0 else out

    # ── interface for the Γ_uni kernel (Phase 3) ─────────────────────────────
    def conversion_mask(self, T_K, n_max: int) -> np.ndarray:
        """Boolean array over sizes 1..n_max where ΔF(n,T) > 0 (⟨100⟩ favoured).

        The unary-conversion kernel Γ_uni(n,T) is nonzero only on this support
        (and is further gated by the [1 − exp(−ΔF/k_BT)] factor in
        ``reaction_rates``).  Independent of the ⟨100⟩ population's admissible
        size floor, which is the loop-onset size, not this thermodynamic set.
        """
        n = np.arange(1, int(n_max) + 1)
        return self.delta_F(n, T_K) > 0.0

    def driving_force_array(self, T_K, n_max: int) -> np.ndarray:
        """ΔF(n,T) over sizes n = 1..n_max  [eV] (0-indexed by size−1)."""
        n = np.arange(1, int(n_max) + 1)
        return self.delta_F(n, T_K)


# ── self-test ─────────────────────────────────────────────────────────────────
def _selftest() -> int:
    """Verify the construction-guaranteed properties (not absolute eV values)."""
    le = LoopEnergetics()
    Tc = le.Tc
    Tstar = le.T_star_C + _C_TO_K
    nref = le.n_ref
    fails = 0

    def check(name, cond, *info):
        nonlocal fails
        if cond:
            print(f"  PASS  {name}   {info if info else ''}")
        else:
            fails += 1
            print(f"  FAIL  {name}   {info}")

    # 1. Calibration zeroes dF at (n_ref, T*).
    check("dF(n_ref, T*) ~ 0",
          abs(le.delta_F(nref, Tstar)) < 1e-9, le.delta_F(nref, Tstar))
    # 2. <111> favoured below T*, <100> favoured above (at n_ref).
    check("dF(n_ref, T*-100K) < 0", le.delta_F(nref, Tstar - 100.0) < 0.0,
          le.delta_F(nref, Tstar - 100.0))
    check("dF(n_ref, T*+100K) > 0", le.delta_F(nref, Tstar + 100.0) > 0.0,
          le.delta_F(nref, Tstar + 100.0))
    # 3. Monotone in T at n_ref (softening monotone => dF increasing in T).
    Ts = np.linspace(300.0, Tc - 1.0, 50)
    dF = le.delta_F(nref, Ts)
    check("dF monotone increasing in T", np.all(np.diff(dF) > -1e-12))
    # 4. Crossover temperature round-trips.
    Tc_ref = le.crossover_temperature(nref)
    check("crossover_temperature(n_ref) ~ T*", abs(Tc_ref - le.T_star_C) < 1e-6,
          Tc_ref)
    # 5. Geometry sanity: both perimeters positive, R* grows with n.
    check("perimeters > 0",
          le.perimeter(10, 111) > 0 and le.perimeter(10, 100) > 0)
    check("R* increasing with n", le.R_star(100, 111) > le.R_star(10, 111))
    # 6. Fhat_100^0 calibrated to a physical magnitude (0.1-1 eV/A).
    check("Fhat_100_0 in (0.1, 1.0) eV/A", 0.1 < le.Fhat_100_0 < 1.0,
          le.Fhat_100_0)

    # Report the emergent size-dependence direction for the user (not asserted).
    T_op = 723.15  # 450 °C operating point
    mask = le.conversion_mask(T_op, 300)
    sizes = np.arange(1, 301)[mask]
    if sizes.size:
        print(f"\n  [info] at {T_op-273.15:.0f} C, dF>0 (<100> favoured) for "
              f"sizes n in [{sizes.min()}, {sizes.max()}]  "
              f"({'small-loop' if sizes.min()==1 else 'large-loop'} biased)")
    else:
        print(f"\n  [info] at {T_op-273.15:.0f} C, no size has dF>0")
    print(f"  [info] Fhat_100_0 = {le.Fhat_100_0:.4f} eV/A "
          f"(calibrated to T*={le.T_star_C:.0f} C at n_ref={nref:.0f})")

    print(f"\n{7 - fails}/7 checks passed")
    return 1 if fails else 0


if __name__ == "__main__":
    import sys
    sys.exit(_selftest())
