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
  y[0..I-1]       : SIA clusters c_n, n=1..I
  y[I..I+V-1]     : void/bubble marginal c_m = Σ_ℓ c_{m,ℓ}, m=1..V
  y[I+V]          : total He in bubbles Q_tot = Σ_{m,ℓ} ℓ·c_{m,ℓ}
  y[I+V+1]        : free He c_h
  N_eq = I + V + 2

Case 2 — fission (decoupled), quasi_steady_state:
  y[0..I-1], y[I..I+V-1], y[I+V] as above
  c_h from QSS:  c_h = (G_He + β·Q_tot) / (Σ K_HeV·c_v + k2_He)
  N_eq = I + V + 1

Case 1 — fusion (mean-field), dynamic:
  y[0..I-1]       : SIA clusters c_n
  y[I..I+V-1]     : marginal c_m^tot
  y[I+V..I+2V-1]  : He content per class Q_m = Σ_ℓ ℓ·c_{m,ℓ}
  y[I+2V]         : free He c_h
  N_eq = I + 2V + 1

Case 1 — fusion (mean-field), quasi_steady_state:
  As above without c_h in state; c_h from QSS.
  N_eq = I + 2V

Concentration floor (C_floor)
------------------------------
Any state variable below C_floor is clamped to C_floor before evaluating
rate terms.  C_floor is set by inp.reactions['C_floor'] (default 1e-15).
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

        I = input_data.I
        V = input_data.V
        self.I = I
        self.V = V

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
            self.i_VAC   = I
            self.i_Q     = I + V        # per-class He content Q_m
            if self.qss_He:
                n_phys = I + 2 * V
                self.i_He = None
            else:
                n_phys = I + 2 * V + 1
                self.i_He = I + 2 * V
        else:
            self.he_mode = 'case2'    # Eq. 175 — decoupled scalar ℓ_tot
            self.i_SIA   = 0
            self.i_VAC   = I
            self.i_Qtot  = I + V        # scalar total He in bubbles
            if self.qss_He:
                n_phys = I + V + 1
                self.i_He = None
            else:
                n_phys = I + V + 2
                self.i_He = I + V + 1

        # Conservation accounting ODEs (5 cumulative integrals)
        self.i_J_SIA_fixed  = n_phys
        self.i_J_SIA_mutual = n_phys + 1
        self.i_J_VAC_fixed  = n_phys + 2
        self.i_J_VAC_mutual = n_phys + 3
        self.i_J_He_sink    = n_phys + 4
        self.N_eq = n_phys + 5

        # Cascade production
        spectrum  = input_data.derived['spectrum']
        G         = input_data.derived['G']
        Pr_SIA, Pr_VAC, G_He = production_rates(G, spectrum, I, V)

        self.Pr_SIA = Pr_SIA[1:]   # [I] 0-indexed: index k → size k+1
        self.Pr_VAC = Pr_VAC[1:]   # [V]
        self.G_He   = G_He

        self.alpha_He = input_data.alpha_He

        # Pre-compute geometry arrays for cavity-absorption inner loops
        self._m_arr   = np.arange(1.0, V + 1)
        self._m13     = self._m_arr ** (1.0 / 3.0)
        self._denom_m = 1.0 + reaction_rates.B_rot * reaction_rates.L_hat**2 / self._m13

        print(f"RateEquations: he_mode='{self.he_mode}'  he_options='{self.he_options}'"
              f"  N_eq={self.N_eq}  I={I}  V={V}"
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
        c_v   : ndarray [V] — void concentration array (clamped)
        Q_tot : float — scalar total He in voids (case2)
        Q_m   : ndarray [V] — per-class He content (case1)

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
          State: [c_1..c_I | c_{-1}..c_{-V} | Q_tot | c_h]
        he_options='quasi_steady_state':
          State: [c_1..c_I | c_{-1}..c_{-V} | Q_tot]
          c_h = (G_He + β·Q_tot) / (Σ K_HeV·c_v + k2_He)
        """
        I       = self.I
        V       = self.V
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
            for m in range(1, min(V, 20) + 1):
                G_VAC_eff[m - 1] = rr.alpha_bubble_fn(m, ell_bar * m**(2.0/3.0))

        # ── SIA cluster equations (Eq. ME_SIA) ─────────────────────────────
        dci = dydt[self.i_SIA:self.i_VAC]
        i_mobile = self.inp.derived['i_mobile']
        reflect  = self.inp.derived.get('boundary_flux', 'absorption') == 'reflection'
        m13     = self._m13
        denom_m = self._denom_m

        # Production
        dci += self.Pr_SIA

        # Thermal SIA emission (gain + loss)
        dci[:-1] += rr.G_SIA[1:] * c_i[1:]
        dci -= rr.G_SIA * c_i

        # i–i monomer growth (gain + loss)
        dci -= rr.K_SIA_grow * ci1 * c_i            # target loss
        dci[1:] += rr.K_SIA_grow[:-1] * ci1 * c_i[:-1]  # gain at n+1
        if reflect:  # suppress I_1 + I_I → I_{I+1} (reaction blocked at wall)
            dci[-1] += rr.K_SIA_grow[-1] * ci1 * c_i[-1]  # undo target loss
            dci[0]  += ci1 * rr.K_SIA_grow[-1] * c_i[-1]  # undo monomer depletion
        # Monomer projectile depletion: every I_n + I_1 → I_{n+1} consumes
        # one monomer as projectile.  At n=1 this gives the Becker-Döring
        # factor of 2 for self-reaction; at n≥2 it is the projectile loss.
        dci[0] -= ci1 * np.dot(rr.K_SIA_grow, c_i)

        # V–I annihilation (gain + loss)
        # n=1: P1 recombination uses K_iv (Eq. monodef_iv)
        dci[0] -= rr.K_iv * cv1 * ci1
        # n≥2: loop shrinkage by mono-vacancy (P3, gain + loss)
        dci[1:] -= rr.K_SIA_shrink[1:] * cv1 * c_i[1:]
        dci[:-1] += rr.K_SIA_shrink[1:] * cv1 * c_i[1:]

        # SIA cluster–cavity absorption: all mobile n = 1..i_mobile
        # n=1 already handled by K_iv above for m=1;
        # for m≥2 cavities absorbing the monomer SIA:
        sum_K3D_cv_m2 = rr.K_3D_cav_pref * np.dot(m13[1:], c_v[1:])
        dci[0] -= ci1 * sum_K3D_cv_m2

        # n=2,3 (3D): absorb into all cavities m=1..V
        sum_K3D_cv = rr.K_3D_cav_pref * np.dot(m13, c_v)
        for n in range(2, min(4, i_mobile + 1)):
            dci[n - 1] -= c_i[n - 1] * sum_K3D_cv

        # n=4..i_mobile (1D/3D mixed): K = K_1D_pref · m^{1/3} / denom
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            dci[n - 1] -= c_i[n - 1] * k_pref * np.dot(m13 / denom_m, c_v)

        # Fixed sinks
        dci -= rr.k2_SIA * c_i

        # ── Vacancy cluster equations (Eq. ME_vac, Case 2) ──────────────────
        dcv = dydt[self.i_VAC:self.i_Qtot]

        # Production
        dcv += self.Pr_VAC

        # Thermal vacancy emission (gain + loss)
        # V_m → V_{m-1} + V_1: loss at m, gain at m-1 (residual), gain at 1 (emitted monomer)
        dcv -= G_VAC_eff * c_v
        dcv[:-1] += G_VAC_eff[1:] * c_v[1:]
        dcv[0] += np.sum(G_VAC_eff[1:] * c_v[1:])   # emitted monomers

        # V–V monomer growth (gain + loss)
        dcv -= rr.K_VAC_grow * cv1 * c_v            # target loss
        dcv[1:] += rr.K_VAC_grow[:-1] * cv1 * c_v[:-1]  # gain at m+1
        if reflect:  # suppress V_1 + V_V → V_{V+1} (reaction blocked at wall)
            dcv[-1] += rr.K_VAC_grow[-1] * cv1 * c_v[-1]  # undo target loss
            dcv[0]  += cv1 * rr.K_VAC_grow[-1] * c_v[-1]  # undo monomer depletion
        # Vacancy monomer projectile depletion (same logic as SIA)
        dcv[0] -= cv1 * np.dot(rr.K_VAC_grow, c_v)

        # SIA-induced cavity shrinkage: all mobile n = 1..i_mobile
        #
        # m=1, n=1: P1 recombination (same K_iv as SIA side)
        dcv[0] -= rr.K_iv * ci1 * cv1
        # m=1, gain from m=2 shrunk by mono-SIA:
        if V >= 2:
            dcv[0] += rr.K_VAC_shrink[1] * ci1 * c_v[1]

        # m≥2: mono-SIA absorption (P2, gain + loss)
        dcv[1:] -= rr.K_VAC_shrink[1:] * ci1 * c_v[1:]
        if V >= 3:
            dcv[1:-1] += rr.K_VAC_shrink[2:] * ci1 * c_v[2:]

        # Vacancy monomer consumed when SIA clusters absorb it:
        # I_n + V_1 → I_{n-1} for n≥2 (n=1 already in K_iv above)
        dcv[0] -= cv1 * np.dot(rr.K_SIA_shrink[1:], c_i[1:])

        # n=2,3 (3D mobile SIA clusters): absorb into cavities
        for n in range(2, min(4, i_mobile + 1)):
            cn = c_i[n - 1]
            dcv -= rr.K_VAC_shrink * cn * c_v
            if n < V:
                dcv[:V - n] += rr.K_VAC_shrink[n:] * cn * c_v[n:]

        # n=4..i_mobile (1D/3D mixed): gain + loss
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            cn = c_i[n - 1]
            k_loss = k_pref * m13 / denom_m
            dcv -= k_loss * cn * c_v
            if n < V:
                mp    = self._m_arr[:V - n] + float(n)
                mp13  = mp ** (1.0 / 3.0)
                denom = 1.0 + rr.B_rot * rr.L_hat**2 / mp13
                dcv[:V - n] += (k_pref * mp13 / denom) * cn * c_v[n:]

        # Fixed sinks — only mobile vacancy clusters (m ≤ v_mobile) diffuse to sinks
        v_mobile = self.inp.derived['v_mobile']
        mask_mobile_v = np.zeros(V)
        mask_mobile_v[:min(v_mobile, V)] = 1.0
        dcv -= rr.k2_vac_scalar * c_v * mask_mobile_v

        # ── Q_tot equation (total He in voids) ───────────────────────────────
        He_uptake  = np.sum(rr.K_HeV * c_h * c_v)
        He_release = self.beta_He * Q_tot
        # He lost only from mobile voids reaching fixed sinks
        ms_arr = np.arange(1.0, V + 1.0)
        ell_m_arr = ell_bar * ms_arr ** (2.0 / 3.0)
        He_sink = rr.k2_vac_scalar * np.sum(ell_m_arr[:v_mobile] * c_v[:v_mobile])
        dydt[self.i_Qtot] = He_uptake - He_release - He_sink

        # ── Free He equation (Eq. 157) — dynamic mode only ───────────────────
        if not self.qss_He:
            dydt[self.i_He] = (self.G_He
                               - He_uptake
                               - rr.k2_He_scalar * c_h
                               + He_release)

        # ── Conservation accounting ODEs (cumulative integrals, exact via ODE) ─
        ns_all = np.arange(1.0, I + 1)

        # (1) SIA content to fixed sinks: Σ_n n · k2_SIA[n] · c_n
        dydt[self.i_J_SIA_fixed] = np.dot(ns_all, rr.k2_SIA * c_i)

        # (2) SIA content to mutual annihilation (recomb + cavity absorption)
        mutual = rr.K_iv * ci1 * cv1                          # V_1+I_1 recomb (1 SIA)
        mutual += ci1 * sum_K3D_cv_m2                          # I_1 → cavities m≥2 (1 SIA)
        mutual += cv1 * np.sum(rr.K_SIA_shrink[1:] * c_i[1:]) # V_1+I_n→I_{n-1} (1 SIA each)
        for n in range(2, min(4, i_mobile + 1)):
            mutual += n * c_i[n - 1] * sum_K3D_cv             # I_n → all cavities (n SIA)
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            mutual += n * c_i[n - 1] * k_pref * np.dot(m13 / denom_m, c_v)
        dydt[self.i_J_SIA_mutual] = mutual

        # (3) VAC content to fixed sinks: Σ_{m≤v_mobile} m · k2_vac · c_m
        vm = min(v_mobile, V)
        ms_mob = np.arange(1.0, vm + 1)
        dydt[self.i_J_VAC_fixed] = rr.k2_vac_scalar * np.dot(ms_mob, c_v[:vm])

        # (3b) VAC content to mutual annihilation:
        # Same as J_SIA_mutual for channel (a), but channel (b) removes m'
        # per vacancy cluster consumed (not min(m',n) as for SIA side).
        vac_mutual = mutual  # start with SIA mutual (channels a all equal)
        # Correction for channel (b): when m' > n, SIA counted n but VAC
        # should count m' (the entire vacancy cluster is consumed).
        for mp in range(2, vm + 1):
            c_mp = c_v[mp - 1]
            if c_mp < 1e-300:
                continue
            for n in range(1, min(mp, I + 1)):
                cn = c_i[n - 1]
                if cn < 1e-300:
                    continue
                # Rate constant for I_n + V_mp reaction
                if n <= 3:
                    K_s = rr.K_3D_cav_pref * float(mp) ** (1.0 / 3.0)
                else:
                    k_pref = rr.K_1D_pref[n - 1]
                    if k_pref < 1e-300:
                        continue
                    K_s = k_pref * float(mp) ** (1.0 / 3.0) / (
                        1.0 + rr.B_rot * rr.L_hat**2 / float(mp) ** (1.0 / 3.0))
                vac_mutual += (mp - n) * K_s * c_mp * cn
        dydt[self.i_J_VAC_mutual] = vac_mutual

        # (4) He to sinks: k2_He · c_h + k2_vac · Σ_{m≤v_mobile} ℓ_m · c_m
        he_sink_acc = rr.k2_He_scalar * c_h
        if C_vac_tot > 1e-300 and Q_tot > 0:
            ell_m_mob = ell_bar * ms_mob ** (2.0 / 3.0)
            he_sink_acc += rr.k2_vac_scalar * np.dot(ell_m_mob, c_v[:vm])
        dydt[self.i_J_He_sink] = he_sink_acc

        return dydt

    # ── Case 1 — fusion (mean-field), Eq. 174 ────────────────────────────────

    def _rhs_case1(self, t, y):
        """
        Mean-field He-reduction (Case 1, fusion, Eq. 174).

        he_options='dynamic':
          State: [c_1..c_I | c_{-1}..c_{-V} | Q_1..Q_V | c_h]
        he_options='quasi_steady_state':
          State: [c_1..c_I | c_{-1}..c_{-V} | Q_1..Q_V]
          c_h = (G_He + β·Σ Q_m) / (Σ K_HeV·c_v + k2_He)
        """
        I       = self.I
        V       = self.V
        rr      = self.rr
        C_floor = self.C_floor
        dydt    = np.zeros(self.N_eq)

        y_SIA = y[self.i_SIA:self.i_VAC]
        y_VAC = y[self.i_VAC:self.i_Q]
        y_Q   = y[self.i_Q:self.i_Q + V]

        c_i = np.maximum(y_SIA, C_floor)
        c_v = np.maximum(y_VAC, C_floor)
        Q_m = np.maximum(y_Q,   C_floor)

        ell_bar_m = Q_m / np.maximum(c_v, 1e-200)   # [V]

        # Free He: QSS or ODE
        if self.qss_He:
            c_h = self.compute_c_h_qss(c_v, Q_m=Q_m)
        else:
            y_He = y[self.i_He]
            c_h  = max(float(y_He), C_floor)

        ci1 = c_i[0]
        cv1 = c_v[0]

        i_mobile = self.inp.derived['i_mobile']
        L_hat   = rr.L_hat
        B_rot   = rr.B_rot

        G_VAC_eff = np.array([
            rr.alpha_bubble_fn(m + 1, ell_bar_m[m]) for m in range(V)
        ])

        # ── SIA clusters (same structure as Case 2) ────────────────────────
        dci = dydt[self.i_SIA:self.i_VAC]
        reflect = self.inp.derived.get('boundary_flux', 'absorption') == 'reflection'
        m13     = self._m13
        denom_m = self._denom_m

        dci += self.Pr_SIA
        dci[:-1] += rr.G_SIA[1:] * c_i[1:]
        dci -= rr.G_SIA * c_i
        dci -= rr.K_SIA_grow * ci1 * c_i            # target loss
        dci[1:] += rr.K_SIA_grow[:-1] * ci1 * c_i[:-1]  # gain
        if reflect:  # suppress I_1 + I_I → I_{I+1} (reaction blocked at wall)
            dci[-1] += rr.K_SIA_grow[-1] * ci1 * c_i[-1]  # undo target loss
            dci[0]  += ci1 * rr.K_SIA_grow[-1] * c_i[-1]  # undo monomer depletion
        # Monomer projectile depletion (same as Case 2)
        dci[0] -= ci1 * np.dot(rr.K_SIA_grow, c_i)

        # V–I annihilation (gain + loss)
        # n=1: P1 recombination uses K_iv (Eq. monodef_iv)
        dci[0] -= rr.K_iv * cv1 * ci1
        # n≥2: loop shrinkage by mono-vacancy (P3, gain + loss)
        dci[1:] -= rr.K_SIA_shrink[1:] * cv1 * c_i[1:]
        dci[:-1] += rr.K_SIA_shrink[1:] * cv1 * c_i[1:]

        # SIA cluster–cavity absorption (all mobile n)
        # n=1: only m≥2 cavities (m=1 already handled by K_iv)
        sum_K3D_cv_m2 = rr.K_3D_cav_pref * np.dot(m13[1:], c_v[1:])
        dci[0] -= ci1 * sum_K3D_cv_m2
        # n=2,3: all cavities m=1..V
        sum_K3D_cv = rr.K_3D_cav_pref * np.dot(m13, c_v)
        for n in range(2, min(4, i_mobile + 1)):
            dci[n - 1] -= c_i[n - 1] * sum_K3D_cv
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            dci[n - 1] -= c_i[n - 1] * k_pref * np.dot(m13 / denom_m, c_v)

        dci -= rr.k2_SIA * c_i

        # ── Vacancy clusters — He capture does NOT change size class m ────────
        dcv = dydt[self.i_VAC:self.i_Q]
        dcv += self.Pr_VAC
        dcv -= G_VAC_eff * c_v
        dcv[:-1] += G_VAC_eff[1:] * c_v[1:]
        dcv[0] += np.sum(G_VAC_eff[1:] * c_v[1:])   # emitted monomers
        dcv -= rr.K_VAC_grow * cv1 * c_v             # target loss
        dcv[1:] += rr.K_VAC_grow[:-1] * cv1 * c_v[:-1]  # gain
        if reflect:  # suppress V_1 + V_V → V_{V+1} (reaction blocked at wall)
            dcv[-1] += rr.K_VAC_grow[-1] * cv1 * c_v[-1]  # undo target loss
            dcv[0]  += cv1 * rr.K_VAC_grow[-1] * c_v[-1]  # undo monomer depletion
        # Vacancy monomer projectile depletion
        dcv[0] -= cv1 * np.dot(rr.K_VAC_grow, c_v)

        # SIA-induced cavity shrinkage (all mobile n, gain + loss)
        # m=1, n=1: P1 recombination
        dcv[0] -= rr.K_iv * ci1 * cv1
        if V >= 2:
            dcv[0] += rr.K_VAC_shrink[1] * ci1 * c_v[1]
        # m≥2: mono-SIA absorption (P2, gain + loss)
        dcv[1:] -= rr.K_VAC_shrink[1:] * ci1 * c_v[1:]
        if V >= 3:
            dcv[1:-1] += rr.K_VAC_shrink[2:] * ci1 * c_v[2:]
        # n=2,3 (3D mobile)
        for n in range(2, min(4, i_mobile + 1)):
            cn = c_i[n - 1]
            dcv -= rr.K_VAC_shrink * cn * c_v
            if n < V:
                dcv[:V - n] += rr.K_VAC_shrink[n:] * cn * c_v[n:]
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            cn = c_i[n - 1]
            k_loss = k_pref * m13 / denom_m
            dcv -= k_loss * cn * c_v
            if n < V:
                mp    = self._m_arr[:V - n] + float(n)
                mp13  = mp ** (1.0 / 3.0)
                denom = 1.0 + rr.B_rot * rr.L_hat**2 / mp13
                dcv[:V - n] += (k_pref * mp13 / denom) * cn * c_v[n:]

        # Vacancy monomer consumed by SIA loop shrinkage:
        # I_n + V_1 → I_{n-1} for n≥2 (n=1 already in K_iv)
        dcv[0] -= cv1 * np.dot(rr.K_SIA_shrink[1:], c_i[1:])

        # Fixed sinks — only mobile vacancy clusters (m ≤ v_mobile)
        v_mobile = self.inp.derived['v_mobile']
        mask_mobile_v = np.zeros(V)
        mask_mobile_v[:min(v_mobile, V)] = 1.0
        dcv -= rr.k2_vac_scalar * c_v * mask_mobile_v

        # ── Q_m equations (He content per void class, Eq. 174) ───────────────
        dQ = dydt[self.i_Q:self.i_Q + V]

        He_cap_total  = 0.0
        He_emit_total = 0.0
        for m in range(1, V + 1):
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

            # He lost only when mobile voids are absorbed at sinks
            if m <= v_mobile:
                dQ[k] -= rr.k2_vac_scalar * Q_m[k]

        # ── Free He (Eq. 157) — dynamic mode only ────────────────────────────
        if not self.qss_He:
            dydt[self.i_He] = (self.G_He
                               - He_cap_total
                               - rr.k2_He_scalar * c_h
                               + He_emit_total)

        # ── Conservation accounting ODEs (cumulative integrals, exact via ODE) ─
        ns_all = np.arange(1.0, I + 1)

        # (1) SIA content to fixed sinks
        dydt[self.i_J_SIA_fixed] = np.dot(ns_all, rr.k2_SIA * c_i)

        # (2) SIA content to mutual annihilation
        mutual = rr.K_iv * ci1 * cv1
        mutual += ci1 * sum_K3D_cv_m2
        mutual += cv1 * np.sum(rr.K_SIA_shrink[1:] * c_i[1:])
        for n in range(2, min(4, i_mobile + 1)):
            mutual += n * c_i[n - 1] * sum_K3D_cv
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            mutual += n * c_i[n - 1] * k_pref * np.dot(m13 / denom_m, c_v)
        dydt[self.i_J_SIA_mutual] = mutual

        # (3) VAC content to fixed sinks
        vm = min(v_mobile, V)
        ms_mob = np.arange(1.0, vm + 1)
        dydt[self.i_J_VAC_fixed] = rr.k2_vac_scalar * np.dot(ms_mob, c_v[:vm])

        # (3b) VAC content to mutual annihilation:
        # Same as J_SIA_mutual for channel (a), but channel (b) removes m'
        # per vacancy cluster consumed (not min(m',n) as for SIA side).
        vac_mutual = mutual  # start with SIA mutual (channels a all equal)
        for mp in range(2, vm + 1):
            c_mp = c_v[mp - 1]
            if c_mp < 1e-300:
                continue
            for n in range(1, min(mp, I + 1)):
                cn = c_i[n - 1]
                if cn < 1e-300:
                    continue
                if n <= 3:
                    K_s = rr.K_3D_cav_pref * float(mp) ** (1.0 / 3.0)
                else:
                    k_pref = rr.K_1D_pref[n - 1]
                    if k_pref < 1e-300:
                        continue
                    K_s = k_pref * float(mp) ** (1.0 / 3.0) / (
                        1.0 + rr.B_rot * rr.L_hat**2 / float(mp) ** (1.0 / 3.0))
                vac_mutual += (mp - n) * K_s * c_mp * cn
        dydt[self.i_J_VAC_mutual] = vac_mutual

        # (4) He to sinks
        C_vac_tot = np.sum(c_v)
        Q_tot_loc = np.sum(Q_m)
        he_sink_acc = rr.k2_He_scalar * c_h
        if C_vac_tot > 1e-300 and Q_tot_loc > 0:
            ell_bar_loc = Q_tot_loc / C_vac_tot
            ell_m_mob = ell_bar_loc * ms_mob ** (2.0 / 3.0)
            he_sink_acc += rr.k2_vac_scalar * np.dot(ell_m_mob, c_v[:vm])
        dydt[self.i_J_He_sink] = he_sink_acc

        return dydt
