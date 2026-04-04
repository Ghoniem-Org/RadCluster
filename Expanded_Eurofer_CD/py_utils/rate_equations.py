"""
rate_equations.py — ODE right-hand side for Expanded_Eurofer_CD.

Implements the full cluster dynamics master equations (Eqs. 152, 155, 157)
and dispatches to the appropriate He-reduction mode based on physics_option.

Physics reference
-----------------
Ghoniem, N.M. (2026), Sections 5-8 (Rate_Equations.pdf):
  Eqs. 152, 155, 157 (master equations)
  Eq. 174 (He Case 1 — fusion/fast_eq)
  Eq. 175 (He Case 2 — fission/decoupled)

He options (he_options)
------------------------
'dynamic'            — free He c_h integrated as a full ODE variable (Eq. 157)
'quasi_steady_state' — c_h eliminated; computed algebraically from dc_h/dt = 0
                       at each RHS call (valid because E_m_h = 0.06 eV is small)

State vector y for full_CD modes
---------------------------------
Case 2 — fission (decoupled), dynamic:
  y[0..N-1]       : SIA clusters c_n, n=1..N
  y[N..N+M-1]     : void/bubble marginal c_m = Σ_ℓ c_{m,ℓ}, m=1..M
  y[N+M]          : total He in bubbles Q_tot = Σ_{m,ℓ} ℓ·c_{m,ℓ}
  y[N+M+1]        : free He c_h
  N_eq = N + M + 2

Case 2 — fission (decoupled), quasi_steady_state:
  y[0..N-1], y[N..N+M-1], y[N+M] as above
  c_h from QSS:  c_h = (G_He + β·Q_tot) / (Σ K_HeV·c_v + k2_He)
  N_eq = N + M + 1

Case 1 — fusion (mean-field), dynamic:
  y[0..N-1]       : SIA clusters c_n
  y[N..N+M-1]     : marginal c_m^tot
  y[N+M..N+2M-1]  : He content per class Q_m = Σ_ℓ ℓ·c_{m,ℓ}
  y[N+2M]         : free He c_h
  N_eq = N + 2M + 1

Case 1 — fusion (mean-field), quasi_steady_state:
  As above without c_h in state; c_h from QSS.
  N_eq = N + 2M

Concentration floor (C_floor)
------------------------------
Any state variable below C_floor is clamped to C_floor before evaluating
rate terms.  After computing dydt, derivatives for floor-clamped variables
are clamped to ≥ 0, preventing the integrator from driving concentrations
negative.  C_floor is set by inp.reactions['C_floor'] (default 1e-15).
"""

import numpy as np
from .defect_production import production_rates
from .binding_energies  import ell_max

_kB = 8.617333262e-5

_VALID_HE_OPTIONS = ('dynamic', 'quasi_steady_state')


class RateEquations:
    """
    ODE system for Expanded_Eurofer_CD.

    Dispatches to the appropriate He-reduction mode:
      'full_CD_fission'  → Case 2 decoupled (Eq. 175)
      'full_CD_fusion'   → Case 1 mean-field (Eq. 174)

    The bin_moment_CD modes are handled by BinMomentRateEquations (bin_moment_rates.py).

    Parameters
    ----------
    input_data     : InputData
    reaction_rates : ReactionRates
    """

    def __init__(self, input_data, reaction_rates):
        self.inp = input_data
        self.rr  = reaction_rates

        N = input_data.N
        M = input_data.M
        self.N = N
        self.M = M

        po = input_data.physics_option
        self.physics_option = po

        # Concentration floor — clamp before rate evaluation; gate negative dydt
        self.C_floor = float(input_data.reactions.get('C_floor', 1e-15))

        # Free He mode: 'dynamic' or 'quasi_steady_state'
        raw = str(input_data.reactions.get('he_options', 'dynamic')).lower()
        if raw not in _VALID_HE_OPTIONS:
            import warnings
            warnings.warn(f"Unknown he_options='{raw}'. Using 'dynamic'.")
            raw = 'dynamic'
        self.he_options = raw
        self.qss_He     = (raw == 'quasi_steady_state')

        # Pre-compute beta_He (He de-trapping attempt rate) for QSS and Q_tot RHS
        kBT_val  = input_data.derived['kBT']
        E_b_hV1  = float(input_data.energetics.get('E_b_hV_1', 2.30))
        E_m_h    = float(input_data.derived['E_m_h'])
        nu_h     = float(input_data.derived['nu_h'])
        self.beta_He = nu_h * np.exp(-(E_b_hV1 + E_m_h) / kBT_val)

        # Decide He reduction mode and build state-vector index map
        if 'fusion' in po:
            self.he_mode = 'case1'    # Eq. 174 — mean-field ℓ̄_m
            self.i_SIA   = 0
            self.i_VAC   = N
            self.i_Q     = N + M        # per-class He content Q_m
            if self.qss_He:
                self.N_eq = N + 2 * M
                self.i_He = None        # not in state vector
            else:
                self.N_eq = N + 2 * M + 1
                self.i_He = N + 2 * M  # free He
        else:
            self.he_mode = 'case2'    # Eq. 175 — decoupled scalar ℓ_tot
            self.i_SIA   = 0
            self.i_VAC   = N
            self.i_Qtot  = N + M        # scalar total He in bubbles
            if self.qss_He:
                self.N_eq = N + M + 1
                self.i_He = None
            else:
                self.N_eq = N + M + 2
                self.i_He = N + M + 1  # free He

        # Cascade production
        spectrum  = input_data.derived['spectrum']
        G         = input_data.derived['G']
        Pr_SIA, Pr_VAC, G_He = production_rates(G, spectrum, N, M)

        self.Pr_SIA = Pr_SIA[1:]   # [N] 0-indexed: index k → size k+1
        self.Pr_VAC = Pr_VAC[1:]   # [M]
        self.G_He   = G_He

        self.alpha_He = input_data.alpha_He

        print(f"RateEquations: he_mode='{self.he_mode}'  he_options='{self.he_options}'"
              f"  N_eq={self.N_eq}  N={N}  M={M}"
              f"  G_He={G_He:.3e} at.frac/s"
              f"  C_floor={self.C_floor:.1e}  beta_He={self.beta_He:.3e}")

    def get_initial_conditions(self):
        """Return y0 — all concentrations set to C_floor."""
        return np.full(self.N_eq, self.C_floor)

    # ── QSS He algebraic formula ──────────────────────────────────────────────

    def compute_c_h_qss(self, c_v, Q_tot=None, Q_m=None):
        """
        Quasi-steady-state free He (Eq. 157 with dc_h/dt = 0):

          Case 2:  c_h = (G_He + β·Q_tot) / (Σ K_HeV·c_v + k2_He)
          Case 1:  c_h = (G_He + Σ β·Q_m) / (Σ K_HeV·c_v + k2_He)

        Parameters
        ----------
        c_v   : ndarray [M] — void concentration array (clamped)
        Q_tot : float — scalar total He in voids (case2)
        Q_m   : ndarray [M] — per-class He content (case1)

        Returns
        -------
        c_h : float
        """
        rr   = self.rr
        sink = np.sum(rr.K_HeV * c_v) + rr.k2_He_scalar
        if self.he_mode == 'case2':
            source = self.G_He + self.beta_He * (Q_tot if Q_tot is not None else 0.0)
        else:
            # Case 1: use beta_He * Q_m (same scalar approximation as RHS)
            He_emit = self.beta_He * np.sum(Q_m) if Q_m is not None else 0.0
            source  = self.G_He + He_emit
        return source / max(sink, 1e-300)

    # ── ODE RHS dispatch ─────────────────────────────────────────────────────

    def ode_system(self, t, y):
        """ODE right-hand side.  Dispatches by he_mode."""
        if self.he_mode == 'case1':
            return self._rhs_case1(t, y)
        else:
            return self._rhs_case2(t, y)

    # ── Case 2 — fission (decoupled), Eq. 175 ────────────────────────────────

    def _rhs_case2(self, t, y):
        """
        Decoupled He-reduction (Case 2, fission, Eq. 175).

        he_options='dynamic':
          State: [c_1..c_N | c_{-1}..c_{-M} | Q_tot | c_h]
        he_options='quasi_steady_state':
          State: [c_1..c_N | c_{-1}..c_{-M} | Q_tot]
          c_h = (G_He + β·Q_tot) / (Σ K_HeV·c_v + k2_He)
        """
        N       = self.N
        M       = self.M
        rr      = self.rr
        C_floor = self.C_floor
        dydt    = np.zeros(self.N_eq)

        # Raw state segments (for derivative clamping guard)
        y_SIA  = y[self.i_SIA:self.i_VAC]
        y_VAC  = y[self.i_VAC:self.i_Qtot]
        y_Qtot = y[self.i_Qtot]

        # Clamp to C_floor for rate evaluation
        c_i   = np.maximum(y_SIA,  C_floor)
        c_v   = np.maximum(y_VAC,  C_floor)
        Q_tot = max(float(y_Qtot), C_floor)

        # Free He
        if self.qss_He:
            c_h = self.compute_c_h_qss(c_v, Q_tot=Q_tot)
        else:
            y_He = y[self.i_He]
            c_h  = max(float(y_He), C_floor)

        ci1 = c_i[0]
        cv1 = c_v[0]

        # Mean He loading per void (scalar, Eq. 175)
        C_vac_tot = np.sum(c_v)
        ell_bar   = Q_tot / max(C_vac_tot, 1e-200)

        # Effective thermal void emission corrected for He pressure
        G_VAC_eff = rr.G_VAC.copy()
        if ell_bar > 0.01:
            for m in range(1, min(M, 20) + 1):
                G_VAC_eff[m - 1] = rr.alpha_bubble_fn(m, ell_bar * m**(2.0/3.0))

        # ── SIA cluster equations (Eq. 152) ──────────────────────────────────
        dci = dydt[self.i_SIA:self.i_VAC]

        dci += self.Pr_SIA

        dci[:-1] += rr.G_SIA[1:] * c_i[1:]    # thermal SIA emission n+1 → n

        dci -= rr.K_SIA_grow * ci1 * c_i       # loss: n → n+1
        dci[1:] += rr.K_SIA_grow[:-1] * ci1 * c_i[:-1]  # gain: n-1 → n

        dci -= rr.G_SIA * c_i                  # thermal emission loss

        dci -= rr.K_SIA_shrink * cv1 * c_i     # vacancy annihilation loss
        dci[:-1] += rr.K_SIA_shrink[1:] * cv1 * c_i[1:]  # gain from n+1

        # 1D glide recombination (Eq. 141)
        n_max_i = self.inp.derived['n_max_i']
        L_hat   = rr.L_hat
        B_rot   = rr.B_rot
        for n in range(4, min(N, n_max_i) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            for m in range(1, M + 1):
                m_f   = float(m)
                k_eff = k_pref * m_f**(1.0/3.0) / (1.0 + B_rot * L_hat**2 * m_f**(-1.0/3.0))
                dci[n - 1] -= k_eff * c_i[n - 1] * c_v[m - 1]

        dci -= rr.k2_SIA * c_i                 # fixed sinks

        # ── Vacancy cluster equations (Eq. 155, Case 2) ──────────────────────
        # He capture does NOT change the void size class m — no K_HeV term here;
        # He balance is carried entirely by the Q_tot equation.
        dcv = dydt[self.i_VAC:self.i_Qtot]

        dcv += self.Pr_VAC

        dcv -= G_VAC_eff * c_v                 # thermal emission loss
        dcv[:-1] += G_VAC_eff[1:] * c_v[1:]   # gain from m+1 emitting

        dcv -= rr.K_VAC_grow * cv1 * c_v       # growth loss
        dcv[1:] += rr.K_VAC_grow[:-1] * cv1 * c_v[:-1]   # gain from m-1

        dcv -= rr.K_VAC_shrink * ci1 * c_v     # SIA annihilation loss
        dcv[:-1] += rr.K_VAC_shrink[1:] * ci1 * c_v[1:]

        for n in range(4, min(N, n_max_i) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            for m in range(1, M + 1):
                m_f   = float(m)
                k_eff = k_pref * m_f**(1.0/3.0) / (1.0 + B_rot * L_hat**2 * m_f**(-1.0/3.0))
                dcv[m - 1] -= k_eff * c_i[n - 1] * c_v[m - 1]

        dcv -= rr.k2_vac_scalar * c_v

        # ── Q_tot equation (total He in voids) ───────────────────────────────
        He_uptake  = np.sum(rr.K_HeV * c_h * c_v)
        He_release = self.beta_He * Q_tot
        He_sink    = rr.k2_vac_scalar * Q_tot   # He lost when voids absorbed at sinks
        dydt[self.i_Qtot] = He_uptake - He_release - He_sink

        # ── Free He equation (Eq. 157) — dynamic mode only ───────────────────
        if not self.qss_He:
            dydt[self.i_He] = (self.G_He
                               - He_uptake
                               - rr.k2_He_scalar * c_h
                               + He_release)

        # ── Derivative floor guard ────────────────────────────────────────────
        dci[(y_SIA  < C_floor) & (dci < 0.0)] = 0.0
        dcv[(y_VAC  < C_floor) & (dcv < 0.0)] = 0.0
        if float(y_Qtot) < C_floor and dydt[self.i_Qtot] < 0.0:
            dydt[self.i_Qtot] = 0.0
        if not self.qss_He and float(y[self.i_He]) < C_floor and dydt[self.i_He] < 0.0:
            dydt[self.i_He] = 0.0

        return dydt

    # ── Case 1 — fusion (mean-field), Eq. 174 ────────────────────────────────

    def _rhs_case1(self, t, y):
        """
        Mean-field He-reduction (Case 1, fusion, Eq. 174).

        he_options='dynamic':
          State: [c_1..c_N | c_{-1}..c_{-M} | Q_1..Q_M | c_h]
        he_options='quasi_steady_state':
          State: [c_1..c_N | c_{-1}..c_{-M} | Q_1..Q_M]
          c_h = (G_He + β·Σ Q_m) / (Σ K_HeV·c_v + k2_He)
        """
        N       = self.N
        M       = self.M
        rr      = self.rr
        C_floor = self.C_floor
        dydt    = np.zeros(self.N_eq)

        y_SIA = y[self.i_SIA:self.i_VAC]
        y_VAC = y[self.i_VAC:self.i_Q]
        y_Q   = y[self.i_Q:self.i_Q + M]

        c_i = np.maximum(y_SIA, C_floor)
        c_v = np.maximum(y_VAC, C_floor)
        Q_m = np.maximum(y_Q,   C_floor)

        ell_bar_m = Q_m / np.maximum(c_v, 1e-200)   # [M]

        # Free He: QSS or ODE
        if self.qss_He:
            c_h = self.compute_c_h_qss(c_v, Q_m=Q_m)
        else:
            y_He = y[self.i_He]
            c_h  = max(float(y_He), C_floor)

        ci1 = c_i[0]
        cv1 = c_v[0]

        n_max_i = self.inp.derived['n_max_i']
        L_hat   = rr.L_hat
        B_rot   = rr.B_rot

        G_VAC_eff = np.array([
            rr.alpha_bubble_fn(m + 1, ell_bar_m[m]) for m in range(M)
        ])

        # ── SIA clusters (same as Case 2) ────────────────────────────────────
        dci = dydt[self.i_SIA:self.i_VAC]
        dci += self.Pr_SIA
        dci[:-1] += rr.G_SIA[1:] * c_i[1:]
        dci -= rr.K_SIA_grow * ci1 * c_i
        dci[1:] += rr.K_SIA_grow[:-1] * ci1 * c_i[:-1]
        dci -= rr.G_SIA * c_i
        dci -= rr.K_SIA_shrink * cv1 * c_i
        dci[:-1] += rr.K_SIA_shrink[1:] * cv1 * c_i[1:]
        for n in range(4, min(N, n_max_i) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            for m in range(1, M + 1):
                m_f   = float(m)
                k_eff = k_pref * m_f**(1.0/3.0) / (1.0 + B_rot * L_hat**2 * m_f**(-1.0/3.0))
                dci[n - 1] -= k_eff * c_i[n - 1] * c_v[m - 1]
        dci -= rr.k2_SIA * c_i

        # ── Vacancy clusters — He capture does NOT change size class m ────────
        dcv = dydt[self.i_VAC:self.i_Q]
        dcv += self.Pr_VAC
        dcv -= G_VAC_eff * c_v
        dcv[:-1] += G_VAC_eff[1:] * c_v[1:]
        dcv -= rr.K_VAC_grow * cv1 * c_v
        dcv[1:] += rr.K_VAC_grow[:-1] * cv1 * c_v[:-1]
        dcv -= rr.K_VAC_shrink * ci1 * c_v
        dcv[:-1] += rr.K_VAC_shrink[1:] * ci1 * c_v[1:]
        for n in range(4, min(N, n_max_i) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            for m in range(1, M + 1):
                m_f   = float(m)
                k_eff = k_pref * m_f**(1.0/3.0) / (1.0 + B_rot * L_hat**2 * m_f**(-1.0/3.0))
                dcv[m - 1] -= k_eff * c_i[n - 1] * c_v[m - 1]
        dcv -= rr.k2_vac_scalar * c_v

        # ── Q_m equations (He content per void class, Eq. 174) ───────────────
        dQ = dydt[self.i_Q:self.i_Q + M]

        He_cap_total  = 0.0
        He_emit_total = 0.0
        for m in range(1, M + 1):
            k  = m - 1
            cm = c_v[k]
            ell = ell_bar_m[k]

            he_cap_m       = rr.K_HeV[k] * c_h * cm
            dQ[k]         += he_cap_m
            He_cap_total  += he_cap_m

            alpha_h        = rr.alpha_He_emit_fn(m, max(int(round(ell)), 0))
            he_emit_m      = alpha_h * Q_m[k]
            dQ[k]         -= he_emit_m
            He_emit_total += he_emit_m

            dQ[k] -= rr.k2_vac_scalar * Q_m[k]   # He lost with voids at sinks

        # ── Free He (Eq. 157) — dynamic mode only ────────────────────────────
        if not self.qss_He:
            dydt[self.i_He] = (self.G_He
                               - He_cap_total
                               - rr.k2_He_scalar * c_h
                               + He_emit_total)

        # ── Derivative floor guard ────────────────────────────────────────────
        dci[(y_SIA < C_floor) & (dci < 0.0)] = 0.0
        dcv[(y_VAC < C_floor) & (dcv < 0.0)] = 0.0
        dQ[ (y_Q   < C_floor) & (dQ  < 0.0)] = 0.0
        if not self.qss_He and float(y[self.i_He]) < C_floor and dydt[self.i_He] < 0.0:
            dydt[self.i_He] = 0.0

        return dydt
