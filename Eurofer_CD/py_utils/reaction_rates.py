"""
reaction_rates.py — Pre-computed rate constant arrays for Eurofer_CD.

Implements the Waite (1957) capture rate formulas and thermal emission rates
for bcc Fe / EUROFER97 cluster dynamics.

Physics reference
-----------------
Ghoniem, N.M. (2024), Sections 5.2–5.3 (Rate_Equations.pdf).

Key formulas
------------
Capture radius:
    r_n = r_0 · n^(1/3),   r_0 = (3Ω/4π)^(1/3)             [eq. 86]

3D Waite capture rate (both species 3D diffusers):
    K_AB(n) = 4π · r_AB · (D_A + D_B)                        [eq. 87]
    r_AB = r_A + r_B = r_0 · (n_A^(1/3) + n_B^(1/3))

1D-glide effective capture rate for glissile SIA clusters n ≥ n_1D:
    K_1D(n, ρ) = π · r_n · D_1D(n) / ρ_sink                 [eq. 92, simplified]
    (effective rate averaged over 1D glide segments)

Thermal emission rate:
    α_void(m) = K_vv(m) · C_v_eq_surf(m)                     [eq. 97]
    α_loop(n) = K_vi(n) · exp(−E_b_loop(n) / kT)             [eq. 98]

Array naming convention (mirrors Full_CD for easy comparison)
-------------------------------------------------------------
  KVV[m-1]  — vacancy capture rate: vacancy + vac.cluster(m) → vac.cluster(m+1)
  KVI[m-1]  — interstitial capture rate: SIA(1) + vac.cluster(m) → vac.cluster(m-1)
  KII[n-1]  — interstitial capture rate: SIA(1) + int.cluster(n) → int.cluster(n+1)
  KIV[n-1]  — vacancy capture rate: vac(1) + int.cluster(n) → int.cluster(n-1)
  GVV[m-1]  — thermal vacancy emission: vac.cluster(m) → vac.cluster(m-1) + vac(1)
  GII[n-1]  — thermal SIA emission: int.cluster(n) → int.cluster(n-1) + SIA(1)
  KHeV[m-1] — He capture rate: He(free) + vac.cluster(m) → He-vac.cluster(m,1)
"""

import numpy as np
from .binding_energies import (
    E_b_void, E_b_loop, C_v_eq_surf, capture_radius
)

_kB = 8.617333262e-5   # eV K^-1


class ReactionRates:
    """
    Pre-computed rate constant arrays for the Eurofer_CD ODE system.

    All arrays are 0-indexed; index k corresponds to cluster size k+1.

    Parameters
    ----------
    input_data : InputData
    """

    def __init__(self, input_data):
        self.inp = input_data
        self._precompute()

    # ── Pre-computation ──────────────────────────────────────────────────────

    def _precompute(self):
        p    = self.inp.material_params
        d    = self.inp.derived
        mp   = self.inp.model_params
        Nv   = self.inp.Nv
        Ni   = self.inp.Ni
        T    = float(p['T'])
        kBT  = _kB * T

        Omega   = d['Omega']
        r0      = d['r0']
        Di      = d['Di']
        Dv      = d['Dv']
        DHe     = d['DHe']
        Cv_eq   = d['Cv_eq']
        gamma_s = float(p['gamma_s'])
        E_f_v   = float(p['E_f_v'])
        E_b_2i  = float(p['E_b_2i'])
        E_b_inf = float(p['E_b_inf_loop'])
        n_trans = float(p.get('n_trans_loop', 8.0))
        n_1D    = int(float(p.get('n_1D', 4)))
        Z_i     = float(p.get('Z_i', 1.05))
        Z_v     = float(p.get('Z_v', 1.00))

        # ── Helper: Waite capture radius [m] for cluster of size n ───────────
        def r_n(n):
            return r0 * float(n)**(1.0 / 3.0)

        # ── Helper: point-defect capture rate 4π·r_AB·(D_A+D_B) ─────────────
        # For a point defect (size 1) reacting with a cluster of size n
        def K_waite_v_cluster(n):
            """Capture rate of a free vacancy (3D) by vacancy cluster(n)."""
            r_AB = r_n(1) + r_n(n)
            # For n=1 (mono + mono): both mobile → 2·Dv
            # For n≥2: cluster stationary → only the free vacancy diffuses
            D_eff = (Dv + Dv) if n == 1 else Dv
            return 4.0 * np.pi * r_AB * D_eff

        def K_waite_i_by_v_cluster(n):
            """Capture rate of a free SIA (3D) by vacancy cluster(n) [annihilation]."""
            r_AB = r_n(1) + r_n(n)
            return 4.0 * np.pi * r_AB * (Di + 0.0)  # vacancy cluster stationary

        def K_waite_i_by_i_cluster(n):
            """Capture rate of free SIA (3D) by interstitial cluster(n) → growth."""
            r_AB = r_n(1) + r_n(n)
            return 4.0 * np.pi * r_AB * (Di + 0.0)  # cluster stationary

        def K_waite_v_by_i_cluster(n):
            """Capture rate of free vacancy (3D) by interstitial cluster(n) → shrink."""
            r_AB = r_n(1) + r_n(n)
            return 4.0 * np.pi * r_AB * (Dv + 0.0)

        def K_waite_He_by_cluster(n):
            """Capture rate of free He (3D) by vacancy cluster(n)."""
            r_AB = r_n(1) + r_n(n)
            return 4.0 * np.pi * r_AB * (DHe + 0.0)

        # ── Effective capture rate for 1D-gliding SIA clusters ───────────────
        # Glissile clusters (n ≥ n_1D) diffuse along ⟨111⟩ channels.
        # Effective 3D-averaged rate (PDF Section 3, eq. 42–47):
        #   K_1D_eff(n) = K_3D(n) · (1 + K_3D(n)/(2π·r_n·D_1D(n)·ρ_sink^(1/2)))^{−1}
        # For a simplified first implementation, use K_3D with a bias correction:
        rho_d = float(p.get('rho_d', 5e14))

        # 1D diffusion coefficient for a glissile SIA cluster of size n:
        # D_1D(n) ≈ D_i / n  (empirical; faster clusters are smaller)
        def D_1D(n):
            return Di / float(n)

        def K_effective_i_cluster(n):
            """
            Effective capture rate of vacancy by a glissile SIA cluster(n).
            For n < n_1D: standard Waite (3D).
            For n ≥ n_1D: mixed 1D/3D effective rate (PDF Section 3).
            """
            K3D = K_waite_v_by_i_cluster(n)
            if n < n_1D:
                return K3D
            # 1D-3D crossover formula (Trinkaus-Singh-Foreman 1992, simplified):
            # K_eff ≈ K_3D  (conservative; full 1D expression requires sink geometry)
            # TODO: implement full eq. 42–47 when sink geometry is specified
            return K3D

        # ── Thermal emission rates ────────────────────────────────────────────

        def gamma_void(m):
            """
            Vacancy emission rate from vacancy cluster(m):
                α_void(m) = K_vv(m) · C_v_eq_surf(m)
            """
            K = K_waite_v_cluster(m)
            C_surf = C_v_eq_surf(m, E_f_v, gamma_s, Omega, T)
            return K * C_surf

        def gamma_loop(n):
            """
            SIA emission rate from interstitial cluster(n):
                α_loop(n) = K_iv(n) · exp(−E_b_loop(n) / kT)
            Only relevant for small n (thermal dissociation negligible for large loops).
            """
            K = K_waite_i_by_i_cluster(n)
            Eb = E_b_loop(n, E_b_2i, E_b_inf, n_trans)
            return K * np.exp(-Eb / kBT)

        # ── Build arrays ──────────────────────────────────────────────────────
        ms = np.arange(1, Nv + 1, dtype=float)   # vacancy cluster sizes
        ns = np.arange(1, Ni + 1, dtype=float)   # interstitial cluster sizes

        # Vacancy cluster arrays (indexed m-1) — raw Waite rates [m^3/s]
        KVV_raw   = np.array([K_waite_v_cluster(m)        for m in ms])
        KVI_raw   = np.array([K_waite_i_by_v_cluster(m)   for m in ms])
        GVV_raw   = np.array([gamma_void(m)                for m in ms])
        KHeV_raw  = np.array([K_waite_He_by_cluster(m)    for m in ms])

        # Interstitial cluster arrays (indexed n-1) — raw Waite rates [m^3/s]
        KII_raw   = np.array([K_waite_i_by_i_cluster(n)   for n in ns])
        KIV_raw   = np.array([K_effective_i_cluster(n)    for n in ns])
        GII_raw   = np.array([gamma_loop(n)                for n in ns])

        # ── Unit normalisation: divide by atomic volume Omega ─────────────────
        # Concentrations in the ODE are in atom fractions (dimensionless).
        # Waite capture rates K [m^3/s] and emission rates G [m^3/s · atom_frac]
        # must be divided by Omega [m^3] to give units of s^-1 per atom_frac
        # that are consistent with the atom-fraction rate equations:
        #
        #   dC [atom_frac/s] = K[m^3/s] / Omega[m^3] × C_A × C_B
        #
        # Dislocation sink strengths k2_disl [s^-1] and the Brinkman alpha [s^-1]
        # are already in atom-fraction-compatible units and must NOT be divided.
        #
        # Note: KVI[0]/Omega + KIV[0]/Omega ≈ K(1,-1)/Omega ≈ alpha (Brinkman),
        # so the `alpha` recombination term in rate_equations.py is REMOVED to
        # avoid double-counting the mono-SIA + mono-vacancy reaction.
        inv_Omega = 1.0 / Omega

        self.KVV  = KVV_raw  * inv_Omega   # [s^-1 per atom_frac]
        self.KVI  = KVI_raw  * inv_Omega
        self.GVV  = GVV_raw  * inv_Omega   # [s^-1] (first-order emission rate)
        self.KHeV = KHeV_raw * inv_Omega
        self.KII  = KII_raw  * inv_Omega
        self.KIV  = KIV_raw  * inv_Omega
        self.GII  = GII_raw  * inv_Omega

        # Apply Z-factors to network dislocation sinks (stored as scalars)
        self.k2_disl_v = Z_v * rho_d * Dv    # [s^-1] per unit Cv — already correct
        self.k2_disl_i = Z_i * rho_d * Di    # [s^-1] per unit Ci — already correct

        # Dislocation sink rates for SIA clusters n = 1..Ni
        # Small 3D-mobile clusters (n < n_1D): D(n) = Di / n (empirical decrease)
        # Glissile 1D clusters (n ≥ n_1D):    D_1D(n) = Di / n (same formula)
        # Effective rate: k2(n) = Z_i * rho_d * D(n)  [s^-1] — already correct
        k2_SIA = np.zeros(Ni)
        for _n in range(1, Ni + 1):
            k2_SIA[_n - 1] = Z_i * rho_d * (Di / float(_n))
        self.k2_SIA_cluster = k2_SIA   # index k → cluster size k+1

        # He dislocation sink (He is mobile; absorbed at network dislocations)
        Z_He = float(p.get('Z_He', 1.00))
        self.k2_disl_He = Z_He * rho_d * DHe  # [s^-1] per unit C_He — already correct

        # ── SIA-cluster ↔ vacancy-cluster recombination: K_IclV[n-1, m-1] ────
        # Rate for glissile SIA cluster(n≥2) sweeping into stationary vacancy
        # cluster(m).  Uses Waite formula with D_1D(n) = Di/n, then /Omega.
        # n=1 is already covered by KVI (mono-SIA + vac-cluster annihilation).
        K_IclV = np.zeros((Ni, Nv))
        for _n in range(2, Ni + 1):
            d_1d = Di / float(_n)
            for _m in range(1, Nv + 1):
                r_nm = r_n(_n) + r_n(_m)
                K_IclV[_n - 1, _m - 1] = 4.0 * np.pi * r_nm * d_1d * inv_Omega
        self.K_IclV = K_IclV   # shape (Ni, Nv); K_IclV[0,:] == 0 (n=1 excluded)

        print(f"ReactionRates: KVV[0]={self.KVV[0]:.3e}  KII[0]={self.KII[0]:.3e}"
              f"  GVV[0]={self.GVV[0]:.3e}  GVV[1]={self.GVV[1]:.3e}")
        print(f"  k2_disl_v={self.k2_disl_v:.3e}  k2_disl_i={self.k2_disl_i:.3e}"
              f"  k2_disl_He={self.k2_disl_He:.3e}")
        print(f"  K_IclV[1,0]={K_IclV[1,0]:.3e}  K_IclV[1,1]={K_IclV[1,1]:.3e}"
              f"  (n=2 SIA cluster <-> vacancy cluster)")
