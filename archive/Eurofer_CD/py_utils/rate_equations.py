"""
rate_equations.py — ODE right-hand side for Eurofer_CD cluster dynamics.

Implements the full master equations for defect cluster evolution in bcc Fe /
EUROFER97 under neutron irradiation, including helium effects via two
He-vacancy state-space reduction modes.

Physics reference
-----------------
Ghoniem, N.M. (2024), "Formulation of Cluster Dynamics Equations for
Irradiated Ferritic-Martensitic Steels," Section 5.7 (Rate_Equations.pdf).

State vector y[N]
-----------------
'decoupled' and 'fast_eq' modes (default):

    Indices  0 .. Ni-1          SIA clusters    C_{n=1..Ni}
    Indices  Ni .. Ni+Nv-1      Vacancy clusters C_{-m=1..Nv}
    Index    Ni+Nv              Free He          C_{0,1}

'full' mode (not yet implemented — falls back to 'fast_eq'):

    Explicit He-vacancy pairs (m, ℓ) for m=1..Nv, ℓ=0..L_He_max.

He-vacancy reduction modes (PDF Section 5.6.5)
-----------------------------------------------
'decoupled' (fission, low He/dpa ratio):
    He is treated as a modifier of the void binding energy only.
    Mean He loading ⟨ℓ⟩_m = (total He in m-class) / C_{-m} is tracked
    via the free-He equation and a quasi-steady-state closure.
    Effective void emission rate uses E_b_eff(m, ⟨ℓ⟩_m) < E_b_void(m).

'fast_eq' (general, higher He/dpa):
    He distribution within each void class equilibrates rapidly.
    Track marginal C_{-m}^tot = Σ_ℓ C_{-m,ℓ} (includes He-containing voids).
    Mean He loading ⟨ℓ⟩_m computed from He mass balance per void class.
    Effective emission rate uses He-pressure-modified binding energy.

Rate equation structure (Section 5.7)
--------------------------------------
For SIA clusters C_n (Class I, n = 1..Ni):
  dC_n/dt = P_n^i                               [cascade production]
           + α_loop(n+1)·C_{n+1}               [emission from n+1 → n]
           - K_II(n)·C_1·C_n                   [SIA capture → growth]
           - K_IV(n)·C_v·C_n                   [vacancy capture → shrink]
           + K_IV(n+1)·C_v·C_{n+1}             [n+1 shrinks to n]
           - K_disl_i·C_n                       [dislocation sink]

For vacancy clusters C_{-m} (Class II, m = 1..Nv):
  dC_{-m}/dt = P_m^v                            [cascade production]
              + α_void(m+1)·C_{-(m+1)}         [emission from m+1 → m]
              - K_VV(m)·C_v·C_{-m}             [vacancy capture → growth]
              + K_VV(m-1)·C_v·C_{-(m-1)}       [m-1 grows to m]
              - K_VI(m)·C_i·C_{-m}             [SIA annihilation → shrink]
              + K_VI(m-1)·C_i·C_{-(m-1)}       [m-1 loses SIA → m? No]
              - K_HeV(m)·C_He·C_{-m}           [He capture → bubble mode]

For free He C_{0,1}:
  dC_{0,1}/dt = G_He                            [transmutation production]
               - Σ_m K_HeV(m)·C_{-m}·C_{0,1}  [He captured by voids]
               - K_HeV(0)·C_v·C_{0,1}          [He trapped at monovacancy]

Note: In 'decoupled' mode, He merely shifts binding energies; the
He-vacancy explicit clusters are not tracked as separate ODE variables.
"""

import numpy as np
from .binding_energies import E_b_bubble, dE_b_He_dell, _He_fit
from .defect_production import compute_epsilon

_kB = 8.617333262e-5


class RateEquations:
    """
    ODE system for Eurofer_CD cluster dynamics.

    Dispatches to the appropriate he_mode implementation.

    Parameters
    ----------
    input_data    : InputData
    reaction_rates: ReactionRates
    """

    def __init__(self, input_data, reaction_rates):
        self.inp = input_data
        self.rr  = reaction_rates

        Nv = input_data.Nv
        Ni = input_data.Ni
        self.Nv = Nv
        self.Ni = Ni

        he_mode = input_data.derived['he_mode']
        self.he_mode = he_mode

        # State vector length:
        # Ni (SIA) + Nv (vacancy) + 1 (free He)
        self.N_SIA = Ni
        self.N_VAC = Nv
        self.N_He  = 1
        self.N     = Ni + Nv + 1

        # Index offsets
        self.i_SIA = 0        # SIA clusters:    y[0 .. Ni-1]
        self.i_VAC = Ni       # Vacancy clusters: y[Ni .. Ni+Nv-1]
        self.i_He  = Ni + Nv  # Free He:          y[Ni+Nv]

        # Human-readable names for post_process / visualization
        self.concentration_names = (
            [f'Ci{n}' for n in range(1, Ni + 1)]   +
            [f'Cv{m}' for m in range(1, Nv + 1)]   +
            ['C_He']
        )

        # Pre-build cascade production arrays
        self._build_production_arrays()

        print(f"RateEquations: N={self.N} (Ni={Ni}, Nv={Nv}, He=1)  "
              f"he_mode='{he_mode}'")

    # ── Production arrays ────────────────────────────────────────────────────

    def _build_production_arrays(self):
        """
        Compute cascade production rate arrays P_n^i and P_m^v.

        P_n^i = ε_n^(i) · G  [atom frac / s]  for n = 2 … m1
        P_m^v = ε_m^(v) · G  [atom frac / s]  for m = 2 … n1
        Also:
        P_1^i = η · (1 - f_i_cl) · G
        P_1^v = η · (1 - f_v_cl) · G

        References: PDF Section 1, eq. 1–12; defect_production.py.
        """
        p    = self.inp.material_params
        G    = float(p['G'])
        eta  = float(p.get('eta',  0.30))
        f_i  = float(p.get('f_i_cl', 0.58))
        f_v  = float(p.get('f_v_cl', 0.15))
        s_i  = float(p.get('s_i',  1.6))
        s_v  = float(p.get('s_v',  2.5))
        m1   = int(p.get('m1_spec', 20))
        n1   = int(p.get('n1_spec', 10))

        Ni = self.Ni
        Nv = self.Nv

        # Monomer production: η·(1-f_cl)·G  (surviving free point defects)
        P1i = eta * (1.0 - f_i) * G
        P1v = eta * (1.0 - f_v) * G

        # Cluster production spectrum (ε_n for n = 2 … min(m1, Ni))
        Pr_SIA = np.zeros(Ni)
        Pr_VAC = np.zeros(Nv)
        Pr_SIA[0] = P1i
        Pr_VAC[0] = P1v

        m_arr, eps_i = compute_epsilon(f_i, s_i, min(m1, Ni))
        for m_val, eps in zip(m_arr, eps_i):
            idx = int(m_val) - 1
            if 0 <= idx < Ni:
                Pr_SIA[idx] = eps * eta * G

        n_arr, eps_v = compute_epsilon(f_v, s_v, min(n1, Nv))
        for n_val, eps in zip(n_arr, eps_v):
            idx = int(n_val) - 1
            if 0 <= idx < Nv:
                Pr_VAC[idx] = eps * eta * G

        # He transmutation production rate [atom frac / s]
        G_He_per_dpa = float(p.get('G_He_per_dpa', 0.5))
        # G_He = G_He_per_dpa [appm/dpa] × G [dpa/s] × 1e-6 [appm→frac]
        self.G_He = G_He_per_dpa * G * 1.0e-6

        self.Pr_SIA = Pr_SIA
        self.Pr_VAC = Pr_VAC

    # ── Public interface ─────────────────────────────────────────────────────

    def ode_system(self, t, y):
        """Full ODE right-hand side — dispatches to he_mode implementation."""
        return self._rhs(y)

    def get_initial_conditions(self):
        """Build the initial state vector."""
        d    = self.inp.derived
        Cv_eq = d['Cv_eq']

        y0 = np.zeros(self.N)
        # Free vacancy — start near thermal equilibrium
        i0 = self.i_VAC
        y0[i0]     = Cv_eq           # C_{-1} = mono-vacancy
        y0[i0 + 1] = 1.0e-6 * Cv_eq  # C_{-2}
        y0[i0 + 2] = 1.0e-8 * Cv_eq  # C_{-3}

        # Free SIA — very low initial concentration
        y0[self.i_SIA]     = 1.0e-20  # C_{1}
        y0[self.i_SIA + 1] = 1.0e-30  # C_{2}
        y0[self.i_SIA + 2] = 1.0e-35  # C_{3}

        # Free He — essentially zero at start
        y0[self.i_He] = 1.0e-40

        return y0

    # ── ODE right-hand side ──────────────────────────────────────────────────

    def _rhs(self, y):
        """Dispatch to the appropriate he_mode implementation."""
        if self.he_mode == 'decoupled':
            return self._rhs_decoupled(y)
        elif self.he_mode == 'fast_eq':
            return self._rhs_fast_eq(y)
        else:
            # Fall back to decoupled for unknown modes
            return self._rhs_decoupled(y)

    # -- 'decoupled' mode ---------------------------------------------------- #

    def _rhs_decoupled(self, y):
        """
        'decoupled' He mode (PDF Section 5.6.5, option 2).

        He is NOT explicitly tracked as a separate cluster species.
        Instead, the mean He loading ⟨ℓ⟩_m per void class is derived from
        the He mass balance and used to shift the effective void binding energy:

            E_b_eff(m) = E_b_void(m) + ⟨ℓ⟩_m · (∂E_b_He/∂ℓ)

        This is appropriate for fission spectra where He/dpa ~ 0.5–1 appm/dpa.

        State vector:
            y[0..Ni-1]      C_i(n),  n = 1..Ni    SIA clusters
            y[Ni..Ni+Nv-1]  C_v(m),  m = 1..Nv    vacancy clusters
            y[Ni+Nv]        C_He                   free He
        """
        Ni = self.Ni
        Nv = self.Nv
        rr = self.rr
        p  = self.inp.material_params
        d  = self.inp.derived

        dydt = np.zeros(self.N)

        # Unpack state vector — no artificial floor here; clipping to ≥0 is done
        # at segment restart in simulation.py.  A 1e-100 floor creates a
        # discontinuity in the RHS that corrupts LSODA's Jacobian estimate.
        Ci_arr = np.maximum(y[self.i_SIA:self.i_SIA + Ni], 0.0)
        Cv_arr = np.maximum(y[self.i_VAC:self.i_VAC + Nv], 0.0)
        C_He   = max(y[self.i_He], 0.0)

        Ci1 = Ci_arr[0]    # mono-SIA
        Cv1 = Cv_arr[0]    # mono-vacancy
        kBT = _kB * float(p['T'])
        Omega   = d['Omega']
        gamma_s = float(p['gamma_s'])
        E_f_v   = float(p['E_f_v'])

        # ── Mean He loading per void class (quasi-steady-state closure) ──────
        # For 'decoupled' mode: QSS balance within each void class m:
        #   K_HeV(m) · C_He · C_{-m} = β_He_emit · ⟨ℓ⟩_m · C_{-m}
        #   → ⟨ℓ⟩_m = K_HeV(m) · C_He / β_He_emit
        # where β_He_emit = ν_He · exp(−(E_b_HeV + E_m_He) / kT)  [s⁻¹]
        nu_He      = float(p.get('nu_He', 6.25e12))
        E_m_He     = float(p.get('E_m_He', 0.06))
        E_b_HeV    = float(p.get('E_b_HeV', 2.60))
        beta_He    = nu_He * np.exp(-(E_b_HeV + E_m_He) / kBT)

        # Representative ell_mean using m=2 capture rate (for GVV correction)
        # Clamped to [0, L_He_max]
        L_max = float(self.inp.L_He_max)
        ell_mean_global = min(rr.KHeV[1] * C_He / max(beta_He, 1.0e-200), L_max)

        # ── Precompute GVV_eff array for all m=1..Nv (He-pressure correction) ──
        # α_eff(m) = K_vv(m) · C_v_eq_surf(m) · exp(−ΔE_He / kT)
        # ΔE_He = ⟨ℓ⟩_m · ∂E_b_He/∂ℓ   (He pressure lowers emission barrier)
        _m_arr      = np.arange(1, Nv + 1, dtype=float)
        _ell_m      = np.minimum(rr.KHeV[:Nv] * C_He / max(beta_He, 1.0e-200), L_max)
        _mask       = _ell_m >= 1.0e-10
        _ell_c      = np.maximum(_ell_m, 0.1)          # clamped for derivative
        _hd         = _He_fit['delta_He']
        _hb         = _He_fit['beta_He']
        _dE_arr     = _hd * _hb / _m_arr * (_ell_c / _m_arr)**(_hb - 1.0)
        _exp_arg    = np.clip(-_ell_m * _dE_arr / kBT, -100.0, 100.0)
        GVV_eff_arr = np.where(_mask, rr.GVV[:Nv] * np.exp(_exp_arg), rr.GVV[:Nv])

        # ── Precompute SIA-cluster ↔ vacancy-cluster recombination rates ────────
        # Vac_recom[m-1] = Σ_{n=2}^{Ni} K_IclV[n-1,m-1]·C_i(n)
        #   — effective annihilation rate density for vacancy cluster m
        #     from all glissile SIA clusters sweeping through the crystal.
        # SIA_recom[n-1] = Σ_{m=1}^{Nv} K_IclV[n-1,m-1]·C_v(m)
        #   — effective loss rate of SIA cluster n due to encounters with
        #     vacancy clusters.
        # K_IclV[0,:] == 0 (n=1 excluded; handled by KVI already).
        Vac_recom = rr.K_IclV[1:, :Nv].T @ Ci_arr[1:Ni]   # shape (Nv,)
        SIA_recom = rr.K_IclV[:Ni, :Nv] @ Cv_arr[:Nv]     # shape (Ni,)

        # ── SIA cluster equations: dCi(n)/dt, n = 1..Ni ──────────────────────
        dCi = dydt[self.i_SIA:self.i_SIA + Ni]

        # --- n = 1 (mono-SIA) ---
        # Note: mono-SIA + mono-vacancy annihilation is captured by KVI[0] (Di part)
        # and KIV[0] (Dv part), which together give K(1,-1)/Omega (Waite, eq. 102).
        # The separate `alpha` Brinkman term has been removed to avoid double-counting.
        dCi[0] = (
            self.Pr_SIA[0]                                  # cascade source
            - rr.KII[0] * Ci1**2                            # Ci + Ci → C2i (nucleation)
            + 2.0 * rr.GII[1] * Ci_arr[1]                  # C2i → 2×Ci
            - rr.KIV[0] * Cv1 * Ci1                        # Ci annihilated by mono-vac
            - rr.k2_disl_i * Ci1                            # dislocation sink
        )
        # Interaction of Ci1 with larger vacancy clusters (annihilation → cluster shrinks)
        dCi[0] -= np.dot(rr.KVI[1:Nv], Cv_arr[1:Nv]) * Ci1
        # Interaction of Ci1 with larger SIA clusters (capture → cluster grows)
        dCi[0] -= np.dot(rr.KII[1:Ni], Ci_arr[1:Ni]) * Ci1
        # Emission from n≥3 SIA clusters returns one free SIA to pool
        # (analogous to np.dot(GVV[2:Nv], Cv_arr[2:Nv]) for vacancies)
        if Ni >= 3:
            dCi[0] += np.dot(rr.GII[2:Ni], Ci_arr[2:Ni])

        # --- n = 2 (di-SIA) ---
        if Ni >= 2:
            dCi[1] = (
                self.Pr_SIA[1]                               # cascade source
                + 0.5 * rr.KII[0] * Ci1**2                  # Ci + Ci → C2i
                - (rr.KII[1] * Ci1                           # C2i captures Ci → C3i
                   + rr.KIV[1] * Cv1                         # C2i captures Cv → Ci
                   + rr.GII[1]                               # C2i emits SIA → 2Ci
                   + rr.k2_SIA_cluster[1]                    # dislocation sink
                   + SIA_recom[1]) * Ci_arr[1]               # recomb. with vac. clusters
                + rr.GII[2] * Ci_arr[2]                      # C3i emits SIA → C2i
                + rr.KIV[2] * Cv1 * Ci_arr[2]               # C3i absorbs Cv → C2i
            )

        # --- n = 3..Ni (vectorized) ---
        if Ni >= 3:
            _idx = np.arange(2, Ni)                          # n=3..Ni, idx=2..Ni-1
            dCi[_idx] = (
                self.Pr_SIA[_idx]
                + rr.KII[_idx-1] * Ci1 * Ci_arr[_idx-1]
                - rr.KII[_idx]   * Ci1 * Ci_arr[_idx]
                - rr.KIV[_idx]   * Cv1 * Ci_arr[_idx]
                - rr.GII[_idx]   * Ci_arr[_idx]
                - rr.k2_SIA_cluster[_idx] * Ci_arr[_idx]
                - SIA_recom[_idx] * Ci_arr[_idx]
            )
            # C(n+1) → Cn contributions for n=3..Ni-1
            if Ni >= 4:
                _inn = np.arange(2, Ni - 1)
                dCi[_inn] += (rr.GII[_inn+1]   * Ci_arr[_inn+1]
                              + rr.KIV[_inn+1] * Cv1 * Ci_arr[_inn+1])

        # ── Vacancy cluster equations: dCv(m)/dt, m = 1..Nv ─────────────────
        dCv = dydt[self.i_VAC:self.i_VAC + Nv]

        # --- m = 1 (mono-vacancy) ---
        # Note: KVI[0] (Di part) + KIV[0] (Dv part) together give K(1,-1)/Omega
        # for the mono-SIA + mono-vacancy reaction. No separate alpha term needed.
        dCv[0] = (
            self.Pr_VAC[0]                                   # cascade source
            - rr.KVV[0] * Cv1**2                             # Cv + Cv → C2v
            + 2.0 * rr.GVV[1] * Cv_arr[1]                   # C2v → 2×Cv
            - rr.KVI[0] * Ci1 * Cv1                         # mono-SIA annihilates Cv
            - rr.KHeV[0] * C_He * Cv1                       # He trapping at monovac.
            + rr.k2_disl_v * (d['Cv_eq'] - Cv1)             # dislocation sink/source
            - Vac_recom[0] * Cv1                             # SIA-cluster recombination
        )
        # Emission from larger vacancy clusters
        dCv[0] += np.dot(rr.GVV[2:Nv], Cv_arr[2:Nv])
        # Loss to larger vacancy clusters (Cv1 absorbed)
        dCv[0] -= np.dot(rr.KVV[2:Nv], Cv_arr[2:Nv]) * Cv1
        # Annihilation of Cv1 by larger SIA clusters
        dCv[0] -= np.dot(rr.KIV[1:Ni], Ci_arr[1:Ni]) * Cv1

        # --- m = 2 (di-vacancy) ---
        if Nv >= 2:
            dCv[1] = (
                self.Pr_VAC[1]                               # cascade source
                + 0.5 * rr.KVV[0] * Cv1**2                  # Cv + Cv → C2v
                - (rr.KVV[1] * Cv1                           # C2v captures Cv → C3v
                   + rr.KVI[1] * Ci1                         # Ci annihilates C2v → Cv
                   + GVV_eff_arr[1]                          # C2v emits Cv → 2Cv
                   + Vac_recom[1]) * Cv_arr[1]               # SIA-cluster recombination
                + GVV_eff_arr[2] * Cv_arr[2]                # C3v emits Cv → C2v
                + rr.KVI[2] * Ci1 * Cv_arr[2]               # Ci annihilates C3v → C2v
                - rr.KHeV[1] * C_He * Cv_arr[1]             # He capture by C2v
            )

        # --- m = 3..Nv (vectorized) ---
        if Nv >= 3:
            _idx = np.arange(2, Nv)                          # m=3..Nv, idx=2..Nv-1
            dCv[_idx] = (
                self.Pr_VAC[_idx]
                + rr.KVV[_idx-1] * Cv1 * Cv_arr[_idx-1]
                - rr.KVV[_idx]   * Cv1 * Cv_arr[_idx]
                - rr.KVI[_idx]   * Ci1 * Cv_arr[_idx]
                - GVV_eff_arr[_idx] * Cv_arr[_idx]
                - rr.KHeV[_idx]  * C_He * Cv_arr[_idx]
                - Vac_recom[_idx] * Cv_arr[_idx]
            )
            # C(m+1) → Cm contributions for m=3..Nv-1
            if Nv >= 4:
                _inn = np.arange(2, Nv - 1)
                dCv[_inn] += (GVV_eff_arr[_inn+1] * Cv_arr[_inn+1]
                              + rr.KVI[_inn+1] * Ci1 * Cv_arr[_inn+1])

        # ── Free He equation ─────────────────────────────────────────────────
        # dC_He/dt = G_He − Σ_m K_HeV(m)·C_He·C_v(m) − k2_disl_He·C_He
        He_capture = np.dot(rr.KHeV[:Nv], Cv_arr[:Nv]) * C_He
        dydt[self.i_He] = self.G_He - He_capture - rr.k2_disl_He * C_He

        return dydt

    # -- 'fast_eq' mode ---------------------------------------------------- #

    def _rhs_fast_eq(self, y):
        """
        'fast_eq' He mode (PDF Section 5.6.5, option 1).

        He distribution within each void class equilibrates rapidly.
        Track marginal C_{-m}^tot = C_{-m,0} + Σ_{ℓ≥1} C_{-m,ℓ}.
        Mean He loading ⟨ℓ⟩_m from He mass balance per void class.

        In this implementation, the He distribution is assumed Poisson with
        mean ⟨ℓ⟩_m, so the marginal density C_{-m}^tot = Σ_ℓ C_{-m,ℓ} is
        tracked as a single ODE variable, and ⟨ℓ⟩_m is a derived quantity.

        He mass balance per cluster class:
            d/dt (⟨ℓ⟩_m · C_m) ≈ K_HeV(m)·C_He·C_m − β_He_emit(m)·⟨ℓ⟩_m·C_m

        At quasi-steady-state:
            ⟨ℓ⟩_m = K_HeV(m)·C_He / β_He_emit(m)

        where β_He_emit(m) is the He emission rate from a void of size m.
        """
        # For now, implement as a corrected version of 'decoupled' with
        # per-class He loading estimated from local QSS
        Ni = self.Ni
        Nv = self.Nv
        rr = self.rr
        p  = self.inp.material_params
        d  = self.inp.derived

        dydt = np.zeros(self.N)

        Ci_arr = np.maximum(y[self.i_SIA:self.i_SIA + Ni], 0.0)
        Cv_arr = np.maximum(y[self.i_VAC:self.i_VAC + Nv], 0.0)
        C_He   = max(y[self.i_He], 0.0)

        Ci1 = Ci_arr[0]
        Cv1 = Cv_arr[0]
        kBT = _kB * float(p['T'])
        E_b_HeV = float(p.get('E_b_HeV', 2.60))

        # He emission rate from a void (He thermal release rate)
        nu_He  = float(p.get('nu_He', 6.25e12))
        E_m_He = float(p.get('E_m_He', 0.06))
        beta_He_emit = nu_He * np.exp(-(E_b_HeV + E_m_He) / kBT)

        # Per-class mean He loading (QSS)
        # ⟨ℓ⟩_m = K_HeV(m)·C_He / (β_He_emit + K_HeV(m)·C_He/C_m · small_correction)
        ell_mean = np.where(
            beta_He_emit > 0,
            rr.KHeV[:Nv] * C_He / (beta_He_emit + 1.0e-100),
            0.0
        )
        # Clamp to reasonable values
        ell_mean = np.minimum(ell_mean, float(self.inp.L_He_max))

        # Effective void emission rates including He-pressure correction (vectorized)
        GVV_eff = rr.GVV.copy()
        _m_arr  = np.arange(1, Nv + 1, dtype=float)
        _ell_c  = np.maximum(ell_mean[:Nv], 0.1)
        _hd     = _He_fit['delta_He']; _hb = _He_fit['beta_He']
        _dE     = _hd * _hb / _m_arr * (_ell_c / _m_arr)**(_hb - 1.0)
        _corr   = ell_mean[:Nv] * _dE
        _pos    = ell_mean[:Nv] > 0
        GVV_eff[:Nv][_pos] = rr.GVV[:Nv][_pos] * np.exp(-_corr[_pos] / kBT)

        # Reuse the same rate equation structure as 'decoupled' but with GVV_eff
        # (shared code path — call internal _rhs_with_GVV_eff)
        return self._rhs_shared(y, Ci_arr, Cv_arr, C_He, GVV_eff, ell_mean, kBT)

    # -- Shared RHS (used by both modes) ------------------------------------ #

    def _rhs_shared(self, y, Ci_arr, Cv_arr, C_He, GVV_eff, ell_mean, kBT):
        """
        Shared ODE RHS used by both 'decoupled' (after computing GVV_eff externally)
        and 'fast_eq' modes.
        """
        Ni = self.Ni
        Nv = self.Nv
        rr = self.rr
        p  = self.inp.material_params
        d  = self.inp.derived

        dydt = np.zeros(self.N)

        Ci1 = Ci_arr[0]
        Cv1 = Cv_arr[0]

        # ── Precompute SIA-cluster ↔ vacancy-cluster recombination rates ────────
        Vac_recom = rr.K_IclV[1:, :Nv].T @ Ci_arr[1:Ni]   # shape (Nv,)
        SIA_recom = rr.K_IclV[:Ni, :Nv] @ Cv_arr[:Nv]     # shape (Ni,)

        # ── SIA clusters ─────────────────────────────────────────────────────
        dCi = dydt[self.i_SIA:self.i_SIA + Ni]

        dCi[0] = (
            self.Pr_SIA[0]
            - rr.KII[0] * Ci1**2
            + 2.0 * rr.GII[1] * Ci_arr[1]
            - rr.KIV[0] * Cv1 * Ci1                         # mono-SIA annihilated by mono-vac
            - rr.k2_disl_i * Ci1
        )
        dCi[0] -= np.dot(rr.KVI[1:Nv], Cv_arr[1:Nv]) * Ci1
        dCi[0] -= np.dot(rr.KII[1:Ni], Ci_arr[1:Ni]) * Ci1
        # Emission from n≥3 SIA clusters returns one free SIA to pool
        if Ni >= 3:
            dCi[0] += np.dot(rr.GII[2:Ni], Ci_arr[2:Ni])

        if Ni >= 2:
            dCi[1] = (
                self.Pr_SIA[1]
                + 0.5 * rr.KII[0] * Ci1**2
                - (rr.KII[1] * Ci1 + rr.KIV[1] * Cv1 + rr.GII[1]
                   + rr.k2_SIA_cluster[1]
                   + SIA_recom[1]) * Ci_arr[1]               # recomb. with vac. clusters
                + rr.GII[2] * Ci_arr[2]
                + rr.KIV[2] * Cv1 * Ci_arr[2]
            )

        # n = 3..Ni (vectorized)
        if Ni >= 3:
            _idx = np.arange(2, Ni)
            dCi[_idx] = (
                self.Pr_SIA[_idx]
                + rr.KII[_idx-1] * Ci1 * Ci_arr[_idx-1]
                - rr.KII[_idx]   * Ci1 * Ci_arr[_idx]
                - rr.KIV[_idx]   * Cv1 * Ci_arr[_idx]
                - rr.GII[_idx]   * Ci_arr[_idx]
                - rr.k2_SIA_cluster[_idx] * Ci_arr[_idx]
                - SIA_recom[_idx] * Ci_arr[_idx]
            )
            if Ni >= 4:
                _inn = np.arange(2, Ni - 1)
                dCi[_inn] += (rr.GII[_inn+1] * Ci_arr[_inn+1]
                              + rr.KIV[_inn+1] * Cv1 * Ci_arr[_inn+1])

        # ── Vacancy clusters ──────────────────────────────────────────────────
        dCv = dydt[self.i_VAC:self.i_VAC + Nv]

        dCv[0] = (
            self.Pr_VAC[0]
            - rr.KVV[0] * Cv1**2
            + 2.0 * GVV_eff[1] * Cv_arr[1]
            - rr.KVI[0] * Ci1 * Cv1                          # mono-SIA annihilates Cv
            - rr.KHeV[0] * C_He * Cv1
            + rr.k2_disl_v * (d['Cv_eq'] - Cv1)
            - Vac_recom[0] * Cv1                              # SIA-cluster recombination
        )
        dCv[0] += np.dot(GVV_eff[2:Nv], Cv_arr[2:Nv])
        dCv[0] -= np.dot(rr.KVV[2:Nv], Cv_arr[2:Nv]) * Cv1
        dCv[0] -= np.dot(rr.KIV[1:Ni], Ci_arr[1:Ni]) * Cv1

        if Nv >= 2:
            dCv[1] = (
                self.Pr_VAC[1]
                + 0.5 * rr.KVV[0] * Cv1**2
                - (rr.KVV[1] * Cv1 + rr.KVI[1] * Ci1 + GVV_eff[1]
                   + Vac_recom[1]) * Cv_arr[1]               # SIA-cluster recombination
                + GVV_eff[2] * Cv_arr[2]
                + rr.KVI[2] * Ci1 * Cv_arr[2]
                - rr.KHeV[1] * C_He * Cv_arr[1]
            )

        # m = 3..Nv (vectorized)
        if Nv >= 3:
            _idx = np.arange(2, Nv)
            dCv[_idx] = (
                self.Pr_VAC[_idx]
                + rr.KVV[_idx-1] * Cv1 * Cv_arr[_idx-1]
                - rr.KVV[_idx]   * Cv1 * Cv_arr[_idx]
                - rr.KVI[_idx]   * Ci1 * Cv_arr[_idx]
                - GVV_eff[_idx]  * Cv_arr[_idx]
                - rr.KHeV[_idx]  * C_He * Cv_arr[_idx]
                - Vac_recom[_idx] * Cv_arr[_idx]
            )
            if Nv >= 4:
                _inn = np.arange(2, Nv - 1)
                dCv[_inn] += (GVV_eff[_inn+1] * Cv_arr[_inn+1]
                              + rr.KVI[_inn+1] * Ci1 * Cv_arr[_inn+1])

        # ── Free He ───────────────────────────────────────────────────────────
        # dC_He/dt = G_He − Σ_m K_HeV(m)·C_He·C_v(m) − k2_disl_He·C_He
        He_capture = np.dot(rr.KHeV[:Nv], Cv_arr[:Nv]) * C_He
        dydt[self.i_He] = self.G_He - He_capture - rr.k2_disl_He * C_He

        return dydt
