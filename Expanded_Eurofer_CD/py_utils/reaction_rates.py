"""
reaction_rates.py — Pre-computed rate constant arrays for Expanded_Eurofer_CD.

Implements all capture, emission, trap-mutation, and re-solution rates for
the full per-size cluster dynamics system in bcc Fe / EUROFER97.

Physics reference
-----------------
Ghoniem, N.M. (2026), Sections 5-6 (Rate_Equations.pdf):
  Eqs. 109-143, Tables 25, 26, 28, 30.

Key rate formulas
-----------------
Geometric prefactors (Eq. 128):
  A_sph  = (48π²)^{1/3} ≈ 7.818
  A_loop = 8√(π/√3)    ≈ 10.78
  A_1D   = 9/(8π^{2/3}) ≈ 2.632
  B_rot  = (4/π)(8π/3)^{1/3} ≈ 2.627

Spherical 3D capture (Eq. 109, 131):
  K_sph(α, m) = A_sph · m^{1/3} · ω_α^eff  [m^3/s per unit Ω]

Loop capture of SIA by dislocation loops (Eq. 113, 132):
  K_loop(i, n) = A_loop · n^{1/2} · Z_i^loop · ω_i^eff

V–SIA recombination (Eq. 130):
  K_iv = 4√3·π · ω_v^eff  [m^3/s per unit Ω]  (= A_iv · ω_v^eff)

Mixed 1D/3D effective rate for glissile SIA clusters (Eq. 121, 141):
  K_n,m^eff = A_sph · m^{1/3} · ω_n^{1D} / (1 + B_rot · L̂² · m^{-1/3})

Thermal emission (Eq. 122, 138-140):
  α_α(m) = A_sph · (m-1)^{1/3} · ω_α^eff · exp(−E_b(m) / k_B T)

Fixed sinks (Eq. 134-137):
  D_α^d  = Z_α · ρ_d · ω_α^eff · a²
  D_α^gb = π² · D_α^eff / d_g²
  D_α^p  = Z_p · ρ_p · r_p · D_α^eff

State vector convention
-----------------------
All rate constants are dimensionless per-atom-fraction quantities (unit = s^-1)
obtained by dividing the volumetric rate [m^3/s] by Ω [m^3].

This ensures: dC [at.frac/s] = K [s^-1 per at.frac] · C_A · C_B.
"""

import numpy as np
from .binding_energies import (
    E_b_void, E_b_loop_i, E_b_loop_v, E_b_bubble, ell_max,
    Gamma_TM, Gamma_res, atomic_radius
)

_kB   = 8.617333262e-5    # eV K^-1
_J_eV = 6.241509074e18    # J → eV


class ReactionRates:
    """
    Pre-computed rate constant arrays for the Expanded_Eurofer_CD ODE system.

    Arrays are 0-indexed: index k corresponds to cluster size k+1.

    Parameters
    ----------
    input_data : InputData
    """

    def __init__(self, input_data):
        self.inp = input_data
        self._precompute()

    def _precompute(self):
        inp  = self.inp
        d    = inp.derived
        re   = inp.reactions
        ener = inp.energetics
        diff = inp.diffusion

        T      = d['T']
        kBT    = d['kBT']
        a_m    = d['a_m']
        Omega  = d['Omega']
        r0     = d['r0']
        b_111  = d['b_111']

        omega_i = d['omega_i_eff']    # effective jump frequency, Eq. 42
        omega_v = d['omega_v_eff']    # Eq. 48
        omega_h = d['omega_h_eff']

        Di_eff  = d['Di_eff']
        Dv_eff  = d['Dv_eff']
        Dh_eff  = d['Dh_eff']
        D1D     = d['D1D']            # callable D1D(n)
        s_1D    = d['s_1D']
        n_max_i = d['n_max_i']
        m_max_v = d['m_max_v']
        L_hat   = d['L_hat']
        B_rot   = d['B_rot']

        E_f_v   = d['E_f_v']
        gamma_s = d['gamma_s']
        E_s_He  = d['E_s_He']

        N = inp.N
        M = inp.M

        # Geometric prefactors (Eq. 128)
        A_sph  = d['A_sph']    # (48π²)^{1/3} ≈ 7.818
        A_loop = d['A_loop']   # 8√(π/√3) ≈ 10.78
        A_iv   = 4.0 * np.sqrt(3.0) * np.pi  # ≈ 21.77 for K_iv  (Eq. 130)

        # Dislocation sink parameters (Table 26)
        rho_d = float(re.get('rho_d', 1.0e14))
        Z_i   = float(re.get('Z_i',   1.10))
        Z_v   = float(re.get('Z_v',   1.00))
        Z_He  = float(re.get('Z_He',  1.00))
        Z_i_loop = float(re.get('Z_i', 1.10))   # loop bias factor ≈ same as Z_i

        # Grain boundary sink (Eq. 135)
        d_g   = float(re.get('d_g',   5.0e-6))

        # Precipitate sink (Eq. 136)
        rho_p  = float(re.get('rho_p', 1.0e21))
        r_p    = float(re.get('r_p',   5.0e-9))
        Z_p_i  = float(re.get('Z_p_i', 1.0))
        Z_p_v  = float(re.get('Z_p_v', 1.0))

        # Binding energy parameters
        A_111    = float(inp.dissociation.get('A_111',   0.7501))
        B_111    = float(inp.dissociation.get('B_111',   0.3873))
        A_100    = float(inp.dissociation.get('A_100',   0.7160))
        B_100    = float(inp.dissociation.get('B_100',   0.3581))
        n_tr     = float(inp.dissociation.get('n_tr',    25.0))
        sigma_tr = float(inp.dissociation.get('sigma_tr', 5.0))
        gamma_sf = float(inp.dissociation.get('gamma_sf', 0.6))
        alpha_He = inp.alpha_He
        nu0_TM   = float(inp.dissociation.get('nu0_TM', 1.0e12))

        # Re-solution parameter
        spec = d['spectrum']
        b0_key = 'b0_fission' if 'fiss' in spec else 'b0_fusion'
        b0_res = float(re.get(b0_key, 0.01 if 'fiss' in spec else 0.10))
        G      = d['G']

        # Ω^{-2/3}: factor for converting volumetric rate k [m^3/s] to K [s^-1]
        # K = k/Ω = (A_sph·m^{1/3}·D) / Ω = A_sph·m^{1/3}·D·Ω^{-2/3} / Ω^{1/3}
        # Equivalently: K = A_sph·m^{1/3}·D / Ω^{2/3}  [s^-1]  (Eq. 131)
        inv_Omega23 = Omega**(-2.0 / 3.0)   # [m^-2]

        # ── Notation helpers ─────────────────────────────────────────────────
        # K_sph_3D = A_sph · m^{1/3} · D / Ω^{2/3}  [s^-1 per at.frac]  (Eq. 131)
        # D [m^2/s] is the effective diffusivity of the mobile species.
        def K_sph(D, m):
            return A_sph * float(m)**(1.0/3.0) * D * inv_Omega23

        # K_loop = A_loop · n^{1/2} · Z_i^loop · D_i / Ω^{2/3}  (Eq. 132)
        def K_loop(n):
            return A_loop * float(n)**(1.0/2.0) * Z_i_loop * Di_eff * inv_Omega23

        # K_iv recombination = A_iv · D_v / Ω^{2/3}  (Eq. 130)
        K_iv_scalar = A_iv * Dv_eff * inv_Omega23

        # Mixed 1D/3D effective rate for SIA cluster(n) + vacancy cluster(m)
        # Eq. 141:  K_{n,m}^eff = A_sph·m^{1/3}·D_n^{1D} / (Ω^{2/3}·(1+B_rot·L̂²·m^{-1/3}))
        def K_1D_eff(n, m):
            denom = 1.0 + B_rot * L_hat**2 * float(m)**(-1.0/3.0)
            return A_sph * float(m)**(1.0/3.0) * D1D(n) * inv_Omega23 / denom

        # Thermal SIA emission from loop of size n (Eq. 138)
        def alpha_loop(n):
            if n <= 1:
                return 0.0
            Eb = E_b_loop_i(n, A_111, B_111, A_100, B_100, n_tr, sigma_tr)
            return A_sph * max(n - 1.0, 0.0)**(1.0/3.0) * Di_eff * np.exp(-Eb / kBT) * inv_Omega23

        # Thermal vacancy emission from void of size m (Eq. 139)
        def alpha_void(m):
            if m <= 1:
                return 0.0
            Eb = E_b_void(m, E_f_v, gamma_s, Omega)
            return A_sph * max(m - 1.0, 0.0)**(1.0/3.0) * Dv_eff * np.exp(-Eb / kBT) * inv_Omega23

        # Thermal vacancy emission from bubble (m, ell) (Eq. 139 modified)
        def alpha_bubble(m, ell):
            if m <= 1:
                return 0.0
            Eb = E_b_bubble(m, ell, E_f_v, gamma_s, Omega, T)
            Eb = max(Eb, 0.01)   # floor to prevent negative barriers
            return A_sph * max(m - 1.0, 0.0)**(1.0/3.0) * Dv_eff * np.exp(-Eb / kBT) * inv_Omega23

        # Thermal He emission from bubble (m, ell) (Eq. 140)
        def alpha_He_emit(m, ell):
            if ell <= 0:
                return 0.0
            from .binding_energies import E_b_He
            Eb = E_b_He(m, ell, E_s_He, Omega, T)
            Eb = max(Eb, 0.01)
            return A_sph * max(m - 1.0, 0.0)**(1.0/3.0) * Dh_eff * np.exp(-Eb / kBT) * inv_Omega23

        # ── Build arrays for SIA clusters n=1..N ────────────────────────────
        ns = np.arange(1, N + 1, dtype=float)

        # Rotational-correlation factor for 1D/3D mixed transport (Eq. 121)
        # Used in K_SIA_grow, K_SIA_loop, K_SIA_shrink, and k2_SIA below.
        rot_factor = 1.0 + B_rot * L_hat**2    # ≈ 6568 for B_rot=2.627, L_hat=50

        # SIA growth (absorbs mono-SIA): K_sph(D_eff, n)  (Eq. 131, 141)
        # For n < 4 (3D mobile): monomer diffusivity Di_eff is the relevant speed.
        # For 4 <= n <= n_max_i (1D gliders): the cluster sweeps space with its
        #   own effective 3D diffusivity D_n_3D = D1D(n) / rot_factor, which
        #   already includes loop-solute trapping (Eq. 52) and the 1D/3D rotation
        #   correction (Eq. 121). This is the appropriate rate for the cluster
        #   sweeping through the mono-SIA background.
        # For n > n_max_i (immobile large loops): only the monomer can diffuse
        #   to the fixed loop, so Di_eff is used.
        K_SIA_grow_arr = np.zeros(N)
        for ni in range(1, N + 1):
            if ni < 4:
                K_SIA_grow_arr[ni - 1] = K_sph(Di_eff, ni)
            elif ni <= n_max_i:
                D_n_3D = D1D(ni) / rot_factor   # effective 3D via rotation correction
                K_SIA_grow_arr[ni - 1] = K_sph(D_n_3D, ni)
            else:
                K_SIA_grow_arr[ni - 1] = K_sph(Di_eff, ni)
        self.K_SIA_grow = K_SIA_grow_arr

        # SIA loop-capture rate  K_loop(n)  (Eq. 132) — same mobility logic
        K_SIA_loop_arr = np.zeros(N)
        for ni in range(1, N + 1):
            if ni < 4:
                K_SIA_loop_arr[ni - 1] = K_loop(ni)
            elif ni <= n_max_i:
                D_n_3D = D1D(ni) / rot_factor
                K_SIA_loop_arr[ni - 1] = (A_loop * float(ni)**0.5
                                           * Z_i_loop * D_n_3D * inv_Omega23)
            else:
                K_SIA_loop_arr[ni - 1] = K_loop(ni)
        self.K_SIA_loop = K_SIA_loop_arr

        # SIA cluster shrinks by absorbing a vacancy  (Eq. 131)
        # For 1D gliders the cluster also sweeps through the vacancy background,
        # so the effective relative diffusivity is Dv_eff + D_n_3D.
        K_SIA_shrink_arr = np.zeros(N)
        for ni in range(1, N + 1):
            if ni < 4:
                K_SIA_shrink_arr[ni - 1] = K_sph(Dv_eff, ni)
            elif ni <= n_max_i:
                D_n_3D = D1D(ni) / rot_factor
                K_SIA_shrink_arr[ni - 1] = K_sph(Dv_eff + D_n_3D, ni)
            else:
                K_SIA_shrink_arr[ni - 1] = K_sph(Dv_eff, ni)
        self.K_SIA_shrink = K_SIA_shrink_arr

        # Thermal SIA emission from loop (Eq. 138)
        self.G_SIA = np.array([alpha_loop(n) for n in ns])

        # Dislocation sink for SIA clusters (Eq. 134)
        # For 3D-mobile n < 4: use ω_i^eff; for 1D n ≥ 4: use D1D(n)/a²
        # Effective 3D diffusivity for fixed-sink capture (Eq. 134-137).
        # For 1D-gliding clusters (n >= 4, n <= n_max_i): D1D is a 1D transport
        # coefficient; plugging it directly into a 3D spherical-capture formula
        # overestimates dislocation/GB absorption by the rotational-correlation
        # factor (1 + B_rot * L_hat^2).  rot_factor is defined above.
        k2_SIA = np.zeros(N)
        for n in range(1, N + 1):
            if n < 4:
                om = omega_i                    # 3D mobile: use ω_i^eff directly
            elif n <= n_max_i:
                # 1D glider: effective 3D diffusivity reduced by rotational factor
                om = D1D(n) / (a_m**2 * rot_factor)
            else:
                om = 0.0                        # immobile large loops
            k2_d  = Z_i * rho_d * om * a_m**2                     # disloc sink
            k2_gb = np.pi**2 * (om * a_m**2) / d_g**2             # GB sink
            k2_p  = Z_p_i * rho_p * r_p * (om * a_m**2)           # precip sink
            k2_SIA[n - 1] = (k2_d + k2_gb + k2_p)                 # [s^-1]
        self.k2_SIA = k2_SIA

        # Mixed 1D/3D cross-term coefficients for SIA cluster(n) + void(m)
        # Stored as K_1D_eff_n[n-1] — called at runtime with m argument
        # For efficiency: precompute K_1D_pref[n-1] = A_sph · D_n^{1D} / Ω^{2/3}
        K_1D_pref = np.zeros(N)
        for n in range(1, N + 1):
            if n <= n_max_i and n >= 4:
                K_1D_pref[n - 1] = A_sph * D1D(n) * inv_Omega23
        self.K_1D_pref = K_1D_pref   # multiply by m^{1/3}/(1+B_rot·L̂²·m^{-1/3})

        # ── Build arrays for vacancy clusters m=1..M ─────────────────────────
        ms = np.arange(1, M + 1, dtype=float)

        # Vacancy captured by void  K_sph(D_v, m)  (Eq. 131)
        self.K_VAC_grow = np.array([K_sph(Dv_eff, m) for m in ms])

        # SIA captured by void (annihilation)  K_sph(D_i, m)  (Eq. 131)
        self.K_VAC_shrink = np.array([K_sph(Di_eff, m) for m in ms])

        # He captured by void  K_sph(D_h, m)  (Eq. 131)
        self.K_HeV = np.array([K_sph(Dh_eff, m) for m in ms])

        # Thermal vacancy emission from pure void (Eq. 139)
        self.G_VAC = np.array([alpha_void(m) for m in ms])

        # Fixed vacancy sink (Eq. 134-137)
        k2_d_v  = Z_v * rho_d * Dv_eff                     # disloc
        k2_gb_v = np.pi**2 * Dv_eff / d_g**2               # GB
        k2_p_v  = Z_p_v * rho_p * r_p * Dv_eff             # precip
        self.k2_vac_scalar = k2_d_v + k2_gb_v + k2_p_v     # [s^-1]

        # Fixed He sink (Eq. 134-137)
        k2_d_h  = Z_He * rho_d * Dh_eff
        k2_gb_h = np.pi**2 * Dh_eff / d_g**2
        self.k2_He_scalar = k2_d_h + k2_gb_h

        # Fixed SIA (monomer) sink
        k2_d_i  = Z_i * rho_d * Di_eff
        k2_gb_i = np.pi**2 * Di_eff / d_g**2
        k2_p_i  = Z_p_i * rho_p * r_p * Di_eff
        self.k2_SIA_scalar = k2_d_i + k2_gb_i + k2_p_i

        # V–SIA recombination scalar (Eq. 130)
        self.K_iv = K_iv_scalar

        # Store callables for He-vacancy reactions (used in full 2D grid mode)
        self.alpha_bubble_fn  = alpha_bubble
        self.alpha_He_emit_fn = alpha_He_emit
        self.K_1D_eff_fn      = K_1D_eff
        self.Gamma_TM_fn      = lambda m, ell: Gamma_TM(m, ell, T, nu0_TM)
        self.Gamma_res_fn     = lambda ell: Gamma_res(ell, G, b0_res)

        # Scalar physics
        self.B_rot  = B_rot
        self.L_hat  = L_hat
        self.alpha_He = alpha_He

        print(f"ReactionRates: K_SIA_grow[0]={self.K_SIA_grow[0]:.3e}"
              f"  K_VAC_grow[0]={self.K_VAC_grow[0]:.3e}"
              f"  G_VAC[0]={self.G_VAC[0]:.3e}"
              f"  K_iv={self.K_iv:.3e}")
        print(f"  k2_SIA[0]={self.k2_SIA[0]:.3e}"
              f"  k2_vac={self.k2_vac_scalar:.3e}"
              f"  k2_He={self.k2_He_scalar:.3e}")
