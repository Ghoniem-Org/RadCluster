"""
reaction_rates.py — Pre-computed rate constant arrays for RadCluster_2_1.

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

V–SIA recombination (Eq. P1, monodef_iv):
  K_iv = 4√3·π · (ω_i^eff + ω_v^eff)  [s^-1 per (at.frac)^2]

Mixed 1D/3D effective rate for glissile SIA clusters (Eq. 121, 141):
  K_n,m^eff = A_sph · m^{1/3} · ω_n^{1D} / (1 + B_rot · L̂² · m^{-1/3})

Thermal emission (Eq. 122, 138-140):
  α_α(m) = A_sph · (m-1)^{1/3} · ω_α^eff · exp(−E_b(m) / k_B T)

Fixed sinks (Eq. 134-137):
  D_α^d  = Z_α · ρ_d · ω_α^eff · a²
  D_α^gb = π² · D_α^eff / d_g²
  D_α^p  = 4π · Z_p · ρ_p · r_p · D_α^eff        (Eq. P4p; Z_p is a pure bias)

State vector convention
-----------------------
All rate constants are dimensionless per-atom-fraction quantities (unit = s^-1)
obtained by dividing the volumetric rate [m^3/s] by Ω [m^3].

This ensures: dC [at.frac/s] = K [s^-1 per at.frac] · C_A · C_B.
"""

import numpy as np
from .binding_energies import (
    E_b_void, E_b_loop_i, E_b_loop_100, E_b_loop_v, E_b_bubble, ell_max,
    Gamma_TM, Gamma_res, atomic_radius, A_111_from_E_b_i2
)
from .loop_energetics import LoopEnergetics

_kB   = 8.617333262e-5    # eV K^-1
_J_eV = 6.241509074e18    # J → eV


def _opt_pos_float(v):
    """Workbook cell → positive float, or None.

    Returns None for an absent key, a blank cell (pandas NaN), a
    non-numeric entry, or a non-positive value.  Used for *optional override*
    keys, where "not supplied" must be distinguishable from "supplied as 0"
    and must fall back to the legacy code path rather than raise.
    """
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(f) or f <= 0.0:
        return None
    return f


_Z_WOLFER_MAX = 5.0   # clamp: bias is unphysical below r_L ~ r_c/8


def _num(v, default):
    """Workbook cell → float, falling back to ``default`` for a blank/absent
    cell.  Unlike :func:`_opt_pos_float` this accepts zero and negatives — it is
    for keys whose legitimate nominal value may be 0 (``loop_net_xi``,
    ``loop_net_K_rec``).
    """
    if v is None:
        return float(default)
    try:
        f = float(v)
    except (TypeError, ValueError):
        return float(default)
    return f if np.isfinite(f) else float(default)


def effective_n_j_min(n_j_min_junc, i_mobile, frac=0.6):
    """Junction minimum size, made consistent with the mobility cutoff.

    The Marian junction fires only when BOTH partners satisfy
    ``min(n, n') >= n_j_min_junc`` (``reaction_rates.phi_junc``,
    ``rate_kernels.cpp:conv_phi_junc``).  But coalescence partners are mobile
    clusters, so no partner can exceed ``i_mobile``.  With the shipped
    ``n_j_min_junc = 30`` and ``i_mobile < 30`` the condition is unsatisfiable
    and the junction channel is **identically zero** — measured 2026-07-31:
    ``phi_max_junc`` 0.1 and 1.0 gave bit-identical ``f_100`` at
    ``i_mobile = 10``.

    That matters because ``i_mobile`` is sampled over 5-50 in the twin, so the
    channel would switch on discontinuously at ``i_mobile = 30`` and
    ``phi_max_junc`` would screen as inert over most of the prior for a purely
    structural reason.

    Rule: ``n_j_min_eff = min(n_j_min_junc, ceil(frac * i_mobile))``, floored
    at 2 (a junction needs at least two di-interstitials).  With the default
    ``frac = 0.6`` this returns exactly ``n_j_min_junc`` at ``i_mobile = 50``
    (0.6*50 = 30), so the production configuration is **unchanged** and every
    existing result at ``i_mobile >= 50`` is bit-identical; the threshold only
    relaxes where the mobility cutoff would otherwise make it unreachable.
    """
    n_j = float(_num(n_j_min_junc, 30.0))
    i_m = float(_num(i_mobile, 50.0))
    f = float(_num(frac, 0.6))
    scaled = float(np.ceil(f * i_m)) if np.isfinite(f) and f > 0 else n_j
    return float(max(2.0, min(n_j, scaled)))


class ReactionRates:
    """
    Pre-computed rate constant arrays for the RadCluster_2_1 ODE system.

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
        s_3D    = float(diff.get('s_3D', 0.0))  # 3D cluster mobility exponent
        # TODO(Stage3): s_3D has no Excel row; document/expose it or remove.
        i_mobile = d['i_mobile']
        v_mobile = d['v_mobile']
        L_hat   = d['L_hat']
        B_rot   = d['B_rot']

        E_f_v   = d['E_f_v']
        gamma_s = d['gamma_s']
        E_s_He  = d['E_s_He']

        I = inp.I
        V = inp.V

        # Geometric prefactors (Eq. 128)
        A_sph  = d['A_sph']    # (48π²)^{1/3} ≈ 7.818
        A_loop = d['A_loop']   # 8√(π/√3) ≈ 10.78
        A_iv   = 4.0 * np.sqrt(3.0) * np.pi  # ≈ 21.77 for K_iv  (Eq. 130)
        # Spherical-precipitate sink prefactor.  Eq. P4p is
        #   D_a^p = (8*pi*r_p/a) * c_p^N * omega_a^eff
        # and with c_p^N = rho_p*Omega, omega = D/a^2, Omega ~ a^3/2 this reduces
        # to the standard spherical sink strength k^2 = 4*pi*r_p*rho_p.  Until
        # 2026-08-16 the code implemented Z_p*rho_p*r_p*D -- i.e. WITHOUT the
        # 4*pi -- so every precipitate sink was 12.57x too weak and Z_p was not
        # a pure bias factor (a "neutral" precipitate needed Z_p = 12.57).
        # Restored here so Z_p_i / Z_p_v mean what their names say.
        A_PREC = 4.0 * np.pi

        # Dislocation sink parameters (Table 26)
        rho_d = float(re.get('rho_d', 1.0e14))
        # ── Loop → network-dislocation loss (loop_network_loss.tex) ──────────
        # When the LOOP_NETWORK_LOSS channel is active, the static network
        # density rho_d is replaced by the *dynamic* network density rho_net
        # (operator-split: piecewise-constant within a segment, refreshed by
        # simulation.run_adaptive between segments via rebuild_rates).  This is
        # the saturation feedback — a growing network raises the P4 sink
        # strength.  With the flag OFF, no 'rho_net' key is present and the
        # legacy static rho_d is used verbatim (bit-identical).
        self._loop_net_on = int(_num(re.get('LOOP_NETWORK_LOSS', 0), 0)) != 0
        if self._loop_net_on:
            rho_d = float(re.get('rho_net', rho_d))   # dynamic network feeds P4 sink
        Z_i   = float(re.get('Z_i',   1.10))
        Z_v   = float(re.get('Z_v',   1.00))
        Z_He  = float(re.get('Z_He',  1.00))
        # Loop SIA bias Z_i^loop (Eq. P3_i, Table 26) — INDEPENDENT of the
        # network bias Z_i.  The pair's *ratio* sets how the SIA flux partitions
        # between loops and the dislocation network, so they must be separately
        # settable; collapsing them to one number removes that degree of freedom
        # and forces the production terms to absorb the difference.
        # Falls back to Z_i when the key is absent, so a workbook without the
        # row reproduces the old aliased behaviour bit-for-bit.
        Z_i_loop = _num(re.get('Z_i_loop', Z_i), Z_i)

        # Grain / lath boundary sink (Eq. 135)
        d_g   = float(re.get('d_g',   5.0e-6))
        # Boundary bias factors.  Until 2026-08-16 the boundary sink was
        # IDENTICAL for SIAs and vacancies, so it could carry no net bias at
        # all.  That was harmless while d_g = 1 mm made k2_gb = 9.9e6 m^-2
        # (seven orders below the network), but with d_g at the EUROFER97 lath
        # width the boundary carries 6.2e13 m^-2 and its bias is a first-order
        # term in the I/V balance.  Default 1.0 on both -> unbiased, i.e. the
        # legacy behaviour is reproduced bit-for-bit by a workbook without
        # these rows.  He is deliberately left unbiased: the bias is an elastic
        # drift of the Frenkel-pair species, and He is not part of that balance.
        Z_gb_i = float(re.get('Z_gb_i', 1.0))
        Z_gb_v = float(re.get('Z_gb_v', 1.0))

        # Precipitate sink (Eq. 136)
        rho_p  = float(re.get('rho_p', 1.0e21))
        r_p    = float(re.get('r_p',   5.0e-9))
        Z_p_i  = float(re.get('Z_p_i', 1.0))
        Z_p_v  = float(re.get('Z_p_v', 1.0))

        # Binding energy parameters
        A_111    = float(inp.dissociation.get('A_111',   0.7501))
        B_111    = float(inp.dissociation.get('B_111',   0.3873))
        # A_100 / B_100: <100> SIA-loop binding parameters.  Currently unused in
        # E_b_loop_i (only the <111> branch is active), but RESERVED for the
        # upcoming <111>/<100> SIA population split — do NOT delete.
        A_100    = float(inp.dissociation.get('A_100',   0.7160))
        B_100    = float(inp.dissociation.get('B_100',   0.3581))
        n_tr     = float(inp.dissociation.get('n_tr',    25.0))
        sigma_tr = float(inp.dissociation.get('sigma_tr', 5.0))
        gamma_sf = float(inp.dissociation.get('gamma_sf', 0.6))

        # ── Di-interstitial binding energy override (dissociation!E_b_i2) ────
        # The small-n loop branch is E_b^fit(n) = A_111·n^{+B_111}, so its value
        # at n=2 IS the di-interstitial binding energy.  Setting that
        # measurable quantity is preferable to setting the bare amplitude
        # A_111, which has no independent experimental determination.  When
        # 'E_b_i2' is present and positive, it OVERRIDES the workbook A_111:
        #
        #     A_111 <- E_b_i2 · 2^{-B_111}                (Eq. Eb_smalln_fit)
        #
        # A_100 is rescaled by the SAME factor so the <100>/<111> amplitude
        # ratio — which sets the relative stability of the two loop characters,
        # and hence f_100 — is preserved and does not become a free parameter
        # riding along with this override.
        #
        # Absent, blank or non-positive => legacy path, A_111/A_100 used
        # verbatim (bit-identical to the pre-override behaviour).
        E_b_i2 = _opt_pos_float(inp.dissociation.get('E_b_i2', None))
        self.E_b_i2 = E_b_i2
        if E_b_i2 is not None:
            A_111_new = A_111_from_E_b_i2(E_b_i2, B_111)
            if A_111 > 0.0:
                A_100 = A_100 * (A_111_new / A_111)
            A_111 = A_111_new

        # ── Void-binding blend parameters (dissociation!lambda, !A_void_0) ───
        # Read from the workbook rather than taken from the module constants in
        # binding_energies.py.  Both keys have always been present in the
        # workbook but were previously hard-coded and silently ignored — a
        # varied-but-unread parameter is indistinguishable from an inert one.
        lambda_void = float(inp.dissociation.get('lambda',    0.5756))
        A_void_0    = float(inp.dissociation.get('A_void_0',  1.2353))
        alpha_He = inp.alpha_He
        nu0_TM   = float(inp.dissociation.get('nu0_TM', 1.0e12))

        # SIA-loop continuum-elasticity parameters for E_b_loop_i (Eqs. 106-108).
        # Sourced explicitly from Excel (energetics) so the function does NOT
        # fall back to its hard-coded defaults (which carried a stale
        # E_f_i=3.77 eV; the Table-3 / Excel value is 3.64 eV).
        E_f_i    = float(ener.get('E_f_i', 3.64))        # SIA formation energy [eV]
        mu_Pa    = float(ener.get('mu', 82.0)) * 1.0e9   # shear modulus GPa → Pa
        nu_pois  = float(ener.get('nu', 0.29))           # Poisson ratio

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

        # K_iv recombination = A_iv · (D_i + D_v) / Ω^{2/3}  (Eq. P1, monodef_iv)
        # Uses mutual diffusivity: both species mobile in 3D
        K_iv_scalar = A_iv * (Di_eff + Dv_eff) * inv_Omega23

        # Mixed 1D/3D effective rate for SIA cluster(n) + vacancy cluster(m)
        # Eq. 141:  K_{n,m}^eff = A_sph·m^{1/3}·D_n^{1D} / (Ω^{2/3}·(1+B_rot·L̂²·m^{-1/3}))
        def K_1D_eff(n, m):
            denom = 1.0 + B_rot * L_hat**2 * float(m)**(-1.0/3.0)
            return A_sph * float(m)**(1.0/3.0) * D1D(n) * inv_Omega23 / denom

        # Thermal SIA emission from loop of size n (Eq. 138)
        def alpha_loop(n):
            if n <= 1:
                return 0.0
            Eb = E_b_loop_i(n, A_111, B_111, A_100, B_100, n_tr, sigma_tr,
                            E_f_i=E_f_i, G_shear=mu_Pa, b_111=b_111,
                            nu=nu_pois, gamma_sf=gamma_sf, Omega=Omega)
            return A_sph * max(n - 1.0, 0.0)**(1.0/3.0) * Di_eff * np.exp(-Eb / kBT) * inv_Omega23

        # Thermal vacancy emission from void of size m (Eq. 139)
        def alpha_void(m):
            if m <= 1:
                return 0.0
            Eb = E_b_void(m, E_f_v, gamma_s, Omega,
                          lambda_void=lambda_void, A_void_0=A_void_0)
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

        # ── Build arrays for SIA clusters n=1..I ────────────────────────────
        ns = np.arange(1, I + 1, dtype=float)

        # 3D cluster diffusivity: D_n^{3D} = Di_eff / n^{s_3D}
        # s_3D = 0 → all small clusters (n<4) diffuse at Di_eff (original)
        # s_3D > 0 → di- and tri-SIA diffuse slower than monomers
        def Di_cluster_3D(n):
            return Di_eff / float(n)**s_3D

        # Rotational-correlation factor for 1D/3D mixed transport (Eq. 121)
        # Used in K_SIA_grow, K_SIA_loop, K_SIA_shrink, and k2_SIA below.
        rot_factor = 1.0 + B_rot * L_hat**2    # ≈ 6568 for B_rot=2.627, L_hat=50

        # ══ AUTHOR-DIRECTED MODEL OPTIONS (2026-08-31) ═══════════════════════
        # Both default OFF, so every prior run reproduces bit-for-bit.
        #
        # (1) sia_D_flat -- SIZE-INDEPENDENT CLUSTER DIFFUSIVITY.
        #     Legacy: n<4 diffuse at Di_eff/n^s_3D, 4<=n<=i_mobile at
        #     D1D(n)/rot_factor (a factor 2.7e6 suppression at L_hat=1016), and
        #     n>i_mobile exactly 0.  Measured consequence: D falls 1.8e6x between
        #     n=3 and n=4, so the coalescence sum over mobile projectiles --
        #     which the kernel DOES carry, over np=1..i_mobile -- collapses to
        #     the monomer term (n=4,5 contribute ~5e-7 of it).
        #     The 1-D-glide isotropic-equivalent D1D/rot_factor is the rate a
        #     randomly placed SPHERICAL sink sees from a rotating 1-D glider;
        #     Heinisch (docs/Literature, "The effects of one-dimensional glide on
        #     the reaction kinetics of interstitial clusters") shows 1-D migration
        #     gives kinetics no scalar isotropic D reproduces, and Osetsky
        #     ("Stability and mobility of defect clusters and dislocation loops in
        #     metals") and Trinkaus, Singh & Foreman (1993) establish that
        #     1/2<111> loops in bcc Fe stay glissile at ALL sizes with migration
        #     barriers of a few tens of meV, essentially size-independent.
        #     With sia_D_flat=1 every mobile cluster n <= i_mobile diffuses at the
        #     monomer value Di_eff.  Raising i_mobile then genuinely raises the
        #     net interstitial flux instead of adding near-immobile species.
        sia_D_flat = int(_num(re.get('sia_D_flat', 0), 0))

        # (2) Z_loop_model -- SIZE-DEPENDENT LOOP BIAS (Wolfer).
        #     Legacy: Z_i^loop is a CONSTANT and Z_v^loop is hard-wired to 1.0,
        #     so the drift bracket [Z_i^loop D_i c_1 - D_v c_v] carries no n and
        #     v(n) = C n^(1/2) has no zero crossing -- no critical radius at any
        #     size.  That is why <111> mean size was a stiff invariant at n~20
        #     across six independent levers (report S17).
        #
        #     The bias originates in the elastic interaction between a point
        #     defect's relaxation volume and the loop strain field; an SIA has a
        #     much larger relaxation volume than a vacancy, hence Z_i > Z_v.  The
        #     interaction defines a DRIFT CAPTURE RADIUS r_c^j (the distance at
        #     which |E_int| = kT), and because that is a length, the capture
        #     efficiency of a loop of radius r_L necessarily depends on r_L.
        #
        #     Form used (toroidal sink, log capture efficiency):
        #         Z_j^loop(r_L) = ln(8 r_L / r_0) / ln(8 r_L / r_c^j)
        #     -> 1 as r_c^j -> r_0 (neutral sink), and -> 1 as r_L -> infinity
        #     (large loops approach the unbiased geometric limit).  r_0 is the
        #     core radius, taken as b.
        #
        #     SOURCE.  Wolfer & Ashkin, "Diffusion of vacancies and interstitials
        #     to edge dislocations", J. Appl. Phys. 47 (1976) 791; Wolfer, "The
        #     Dislocation Bias", J. Computer-Aided Mater. Des. 14 (2007) 403-417;
        #     loop/toroidal capture efficiency after Rauh & Simon, phys. stat.
        #     sol. (a) 46 (1978) 499, in the form used by Golubov, Barashev &
        #     Stoller, "Radiation Damage Theory", Comprehensive Nuclear Materials
        #     Vol. 1 Ch. 1.13 (2012).
        #
        #     CAVEAT -- READ BEFORE TRUSTING THE NUMBERS.  None of these papers is
        #     in docs/Literature.  The FUNCTIONAL FORM above is standard and I am
        #     confident in it; the capture radii r_c^i / r_c^v are exposed as
        #     parameters (in units of b) rather than hard-coded, with defaults
        #     3.0 b and 1.0 b.  Those defaults are order-of-magnitude values, NOT
        #     taken from a specific table, and should be replaced with values read
        #     from Wolfer (2007) or CNM Ch. 1.13 before any published result.
        #     r_c^v = 1.0 b reproduces the legacy Z_v^loop = 1.0 exactly.
        Z_loop_model = int(_num(re.get('Z_loop_model', 0), 0))
        r_cap_i_b    = float(_num(re.get('r_cap_i_b', 3.0), 3.0))
        r_cap_v_b    = float(_num(re.get('r_cap_v_b', 1.0), 1.0))

        def _r_loop(n):
            """Radius of a prismatic loop of n SIAs: pi r^2 = n*Omega/b."""
            return np.sqrt(max(n, 1.0) * Omega / (np.pi * b_111))

        def _Z_wolfer(n, r_cap_b):
            """Z_j^loop(r_L) = ln(8 r_L / r_0) / ln(8 r_L / r_c^j), r_0 = b."""
            r_L = _r_loop(n)
            r0  = b_111
            rc  = r_cap_b * b_111
            x0, xj = 8.0 * r_L / r0, 8.0 * r_L / rc
            # Guard the singular/inverted branch at r_L <~ r_c/8.
            if xj <= 1.05 or x0 <= 1.0:
                return _Z_WOLFER_MAX
            return float(min(np.log(x0) / np.log(xj), _Z_WOLFER_MAX))

        # SIA growth (absorbs mono-SIA): LOOP geometry for n ≥ 4 (Eq. P3_i)
        # SIA clusters of size n ≥ 4 form prismatic dislocation loops whose
        # capture cross-section scales as the circumference (∝ n^{1/2}), not
        # the surface area of an equivalent sphere (∝ n^{1/3}).
        #
        # For n < 4 (3D mobile point-defect clusters / dumbbells):
        #   Spherical geometry: K = A_sph · n^{1/3} · D_i / Ω^{2/3}
        # For n ≥ 4 (dislocation loops — both mobile and sessile):
        #   Loop geometry:  K = A_loop · n^{1/2} · Z_i^loop · D_i / Ω^{2/3}
        #   The Z_i^loop bias factor reflects preferential SIA capture by the
        #   stress field of the prismatic loop (Eq. P3_i, Table 26).
        K_SIA_grow_arr = np.zeros(I)
        for ni in range(1, I + 1):
            Zi_n = (_Z_wolfer(ni, r_cap_i_b) if Z_loop_model else Z_i_loop)
            if ni < 4 and not sia_D_flat:
                K_SIA_grow_arr[ni - 1] = K_sph(Di_eff, ni)
            elif ni <= i_mobile:
                # sia_D_flat: every mobile cluster carries the MONOMER diffusivity,
                # so the projectile term is Di_eff rather than D1D/rot_factor.
                D_n_3D = Di_eff if sia_D_flat else D1D(ni) / rot_factor
                K_SIA_grow_arr[ni - 1] = (A_loop * float(ni)**0.5
                                           * Zi_n * (Di_eff + D_n_3D) * inv_Omega23)
            else:
                K_SIA_grow_arr[ni - 1] = (A_loop * float(ni)**0.5
                                           * Zi_n * Di_eff * inv_Omega23)
        self.K_SIA_grow = K_SIA_grow_arr
        self.Z_i_loop_arr = np.array(
            [(_Z_wolfer(ni, r_cap_i_b) if Z_loop_model else Z_i_loop)
             for ni in range(1, I + 1)])

        # SIA loop-capture rate  K_loop(n)  (Eq. 132) — same mobility logic
        K_SIA_loop_arr = np.zeros(I)
        for ni in range(1, I + 1):
            if ni < 4:
                K_SIA_loop_arr[ni - 1] = K_loop(ni)
            elif ni <= i_mobile:
                D_n_3D = D1D(ni) / rot_factor
                K_SIA_loop_arr[ni - 1] = (A_loop * float(ni)**0.5
                                           * Z_i_loop * D_n_3D * inv_Omega23)
            else:
                K_SIA_loop_arr[ni - 1] = K_loop(ni)
        self.K_SIA_loop = K_SIA_loop_arr

        # SIA cluster shrinks by absorbing a vacancy  (Eq. P3_v)
        # Same loop geometry for n ≥ 4 but NO bias factor (Z_v^loop = 1.0):
        # vacancy capture by the loop is purely geometric (no elastic preference).
        K_SIA_shrink_arr = np.zeros(I)
        for ni in range(1, I + 1):
            # Legacy Z_v^loop is hard-wired to 1.0.  Under Z_loop_model the
            # vacancy side gets its OWN capture radius, so the drift bracket
            # [Z_i D_i c_1 - Z_v D_v c_v] finally carries n on BOTH sides.
            Zv_n = (_Z_wolfer(ni, r_cap_v_b) if Z_loop_model else 1.0)
            if ni < 4 and not sia_D_flat:
                K_SIA_shrink_arr[ni - 1] = K_sph(Dv_eff, ni)
            elif ni <= i_mobile:
                D_n_3D = Di_eff if sia_D_flat else D1D(ni) / rot_factor
                K_SIA_shrink_arr[ni - 1] = (A_loop * float(ni)**0.5
                                              * Zv_n * (Dv_eff + D_n_3D) * inv_Omega23)
            else:
                K_SIA_shrink_arr[ni - 1] = (A_loop * float(ni)**0.5
                                              * Zv_n * Dv_eff * inv_Omega23)
        self.K_SIA_shrink = K_SIA_shrink_arr

        # Thermal SIA emission from loop (Eq. 138)
        self.G_SIA = np.array([alpha_loop(n) for n in ns])

        # Dislocation sink for SIA clusters (Eq. 134)
        # For 3D-mobile n < 4: use ω_i^eff; for 1D n ≥ 4: use D1D(n)/a²
        # Effective 3D diffusivity for fixed-sink capture (Eq. 134-137).
        # For 1D-gliding clusters (n >= 4, n <= i_mobile): D1D is a 1D transport
        # coefficient; plugging it directly into a 3D spherical-capture formula
        # overestimates dislocation/GB absorption by the rotational-correlation
        # factor (1 + B_rot * L_hat^2).  rot_factor is defined above.
        k2_SIA = np.zeros(I)
        for n in range(1, I + 1):
            if n < 4:
                om = omega_i / float(n)**s_3D   # 3D mobile: ω_i^eff / n^{s_3D}
            elif n <= i_mobile:
                # 1D glider: effective 3D diffusivity reduced by rotational factor
                om = D1D(n) / (a_m**2 * rot_factor)
            else:
                om = 0.0                        # immobile large loops
            k2_d  = Z_i * rho_d * om * a_m**2                     # disloc sink
            k2_gb = Z_gb_i * np.pi**2 * (om * a_m**2) / d_g**2    # lath GB sink
            k2_p  = A_PREC * Z_p_i * rho_p * r_p * (om * a_m**2)  # precip sink
            k2_SIA[n - 1] = (k2_d + k2_gb + k2_p)                 # [s^-1]
        self.k2_SIA = k2_SIA

        # Mixed 1D/3D cross-term coefficients for SIA cluster(n) + void(m)
        # Stored as K_1D_eff_n[n-1] — called at runtime with m argument
        # For efficiency: precompute K_1D_pref[n-1] = A_sph · D_n^{1D} / Ω^{2/3}
        K_1D_pref = np.zeros(I)
        for n in range(1, I + 1):
            if n <= i_mobile and n >= 4:
                K_1D_pref[n - 1] = A_sph * D1D(n) * inv_Omega23
        self.K_1D_pref = K_1D_pref   # multiply by m^{1/3}/(1+B_rot·L̂²·m^{-1/3})

        # ── Build arrays for vacancy clusters m=1..V ─────────────────────────
        ms = np.arange(1, V + 1, dtype=float)

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
        k2_gb_v = Z_gb_v * np.pi**2 * Dv_eff / d_g**2      # lath GB
        k2_p_v  = A_PREC * Z_p_v * rho_p * r_p * Dv_eff    # precip
        self.k2_vac_scalar = k2_d_v + k2_gb_v + k2_p_v     # [s^-1]

        # Fixed He sink (Eq. 134-137)
        k2_d_h  = Z_He * rho_d * Dh_eff
        k2_gb_h = np.pi**2 * Dh_eff / d_g**2
        self.k2_He_scalar = k2_d_h + k2_gb_h

        # Fixed SIA (monomer) sink
        k2_d_i  = Z_i * rho_d * Di_eff
        k2_gb_i = Z_gb_i * np.pi**2 * Di_eff / d_g**2
        k2_p_i  = A_PREC * Z_p_i * rho_p * r_p * Di_eff
        self.k2_SIA_scalar = k2_d_i + k2_gb_i + k2_p_i

        # V–SIA recombination scalar (Eq. 130)
        self.K_iv = K_iv_scalar

        # Store callables for He-vacancy reactions (used in full 2D grid mode).
        # These are nested functions — excluded from pickling via __getstate__.
        self.alpha_bubble_fn  = alpha_bubble
        self.alpha_He_emit_fn = alpha_He_emit
        self.K_1D_eff_fn      = K_1D_eff

        # 3D cavity absorption prefactor: A_sph · Di_eff / Ω^{2/3}
        # Used for mobile SIA clusters n=1..3 hitting cavities (Eq. K_cav, 3D branch)
        self.K_3D_cav_pref = A_sph * Di_eff * inv_Omega23

        # ── Effective 3D diffusivities for mobile clusters (for coalescence) ──
        # D_SIA_eff[n-1] = effective 3D diffusivity of SIA cluster of size n.
        # Used by the general i–i coalescence and mobile-SIA–cavity terms.
        # ONLY mobile clusters (n ≤ i_mobile) have non-zero diffusivity.
        # For n < 4 AND n ≤ i_mobile (3D mobile): D_n = Di_eff / n^{s_3D}
        # For 4 ≤ n ≤ i_mobile (1D gliders): D_n = D1D(n) / rot_factor
        # For n > i_mobile: D_n = 0 (sessile — no coalescence as projectile)
        D_SIA_eff = np.zeros(I)
        for n in range(1, I + 1):
            if n > i_mobile:
                pass  # sessile: D = 0
            elif sia_D_flat:
                D_SIA_eff[n - 1] = Di_eff        # size-independent (author, 2026-08-31)
            elif n < 4:
                D_SIA_eff[n - 1] = Di_cluster_3D(n)
            else:
                D_SIA_eff[n - 1] = D1D(n) / rot_factor
        self.D_SIA_eff = D_SIA_eff

        # ── ½⟨111⟩ LOOP COARSENING (LOOP_COAL, default off) ───────────────────
        # D_SIA_eff above is hard-truncated to zero at i_mobile, so a loop that
        # grows past the cutoff can never take part in a coalescence again.
        # Two mobile loops of n̄ ≈ 33 merge to 66 > i_mobile = 40, and the
        # product is frozen.  That is why the ½⟨111⟩ mean size is pinned AT the
        # cutoff: measured mean_n_111/i_mobile = 1.32 ± 0.75 at i_mobile = 40 and
        # 0.84 ± 0.25 at i_mobile = 100 across 16 runs (2026-08-16), i.e. the
        # mean tracks the cutoff rather than the physics.  Real ½⟨111⟩ loops in
        # bcc Fe stay glissile to large size and coarsen by 1-D glide.
        #
        # D_loop_coal continues the SAME 1-D glide law past the cutoff instead
        # of truncating it.  It is a SEPARATE array, used ONLY by the loop–loop
        # coalescence edge: reusing D_SIA_eff would also make large loops visible
        # to the fixed-sink and cavity-absorption channels, which is a different
        # physics change and would confound the test.
        self.loop_coal = int(_num(re.get('LOOP_COAL', 0), 0))
        self.loop_coal_pref = float(_num(re.get('loop_coal_pref', 1.0), 1.0))
        D_loop_coal = np.zeros(I)
        if self.loop_coal:
            for n in range(1, I + 1):
                if n < 4:
                    D_loop_coal[n - 1] = Di_cluster_3D(n)
                else:
                    D_loop_coal[n - 1] = D1D(n) / rot_factor
            D_loop_coal *= self.loop_coal_pref
        self.D_loop_coal = D_loop_coal

        # D_VAC_eff[m-1] = effective 3D diffusivity of vacancy cluster of size m.
        # D_m = D_v / m^{s_vc} for m ≤ v_mobile; 0 for immobile clusters.
        s_vc = float(diff.get('s_vc', 1.0))
        D_VAC_eff = np.zeros(V)
        for m in range(1, V + 1):
            if m <= v_mobile:
                D_VAC_eff[m - 1] = Dv_eff / float(m) ** s_vc
            # else: 0 (sessile)
        self.D_VAC_eff = D_VAC_eff

        # Geometric prefactor for coalescence: A_sph / Ω^{2/3}
        self.A_sph_inv_O23 = A_sph * inv_Omega23
        self.A_loop_inv_O23 = A_loop * inv_Omega23

        # ── Phase-3: ½⟨111⟩ → ⟨100⟩ loop-conversion kernels ──────────────────
        # See docs/design_notes/loop_111_to_100_conversion.md.  These are the
        # 1-D ingredients owned by ReactionRates; the 2-D junction / absorption
        # coalescence kernels are assembled from D_SIA_eff + phi_junc in the
        # EUROFER declaration layer (build_eurofer_rag).
        le = LoopEnergetics(a=a_m * 1.0e10, Omega=Omega * 1.0e30,
                            T_star_C=float(re.get('T_star_conv_C', 450.0)),
                            n_ref=float(re.get('n_ref_conv', 50.0)))
        self.loop_energetics = le

        # (1) Unary (Dudarev) conversion rate Γ_uni(n) at the operating T:
        #       Γ_uni(n) = ν₀·exp(−E_a(n)/kT)·max(0, 1 − exp(−ΔF(n,T)/kT)),
        #       E_a(n)   = E_a0 + γ_a·P_111(n)/b_111   (size-dependent barrier).
        #     Nonzero only where ΔF(n,T) > 0 (small-loop biased — see note).
        n_loop_min = int(re.get('n_loop_min', 4))         # ⟨100⟩ loop-onset floor
        # Unary = DIRECT single-loop ½⟨111⟩→⟨100⟩ rotation (no junction partner to
        # enable the easy two-step path), so its barrier is Marian's direct-
        # rotation value (>2 eV) — NOT the two-step ΔH₂.  Calibrated to place the
        # f₁₁₁(T) crossover at ~350 °C at reactor dose rates (Phase 6).
        E_a0     = float(re.get('E_a0_conv',    1.8))     # eV  (direct-rotation barrier)
        gamma_a  = float(re.get('gamma_a_conv', 0.03))    # eV per perimeter-segment
        nu0_conv = float(re.get('nu0_conv',     1.0e13))  # s⁻¹ attempt frequency
        dF   = le.driving_force_array(T, I)               # eV over n = 1..I
        P111 = np.asarray(le.perimeter(ns, 111), dtype=float)   # Å
        E_a  = E_a0 + gamma_a * P111 / le.b_111           # eV
        # Driving-force gate, evaluated only on the ΔF>0 support to avoid exp
        # overflow on the large-loop (ΔF≪0) tail that np.where would discard.
        gate = np.zeros(I)
        pos  = dF > 0.0
        gate[pos] = 1.0 - np.exp(-dF[pos] / kBT)
        Gamma_uni = nu0_conv * np.exp(-E_a / kBT) * gate
        # No ⟨100⟩ loop exists below the loop-onset floor → no unary conversion.
        Gamma_uni[:max(n_loop_min - 1, 0)] = 0.0
        self.Gamma_uni     = Gamma_uni
        self.conversion_dF = dF                           # diagnostic

        # (2) Marian junction branching φ_junc[n−1, n′−1] (size comparability):
        #       φ = φ_max·exp(−(ln(n/n′))²/2σ_s²)·Θ(min(n,n′) ≥ n_j_min).
        # The dense [I, I] matrix is NOT built here: it is an O(I²) object that
        # only the Python GraphWalker consumes (the C++ solver computes φ on the
        # fly per (n,n') pair).  At production sizes (I ~ 1e5) materialising it
        # would need ~80 GB.  We store the scalar parameters now and expose the
        # matrix through the cached ``phi_junc`` property, so it is built only if
        # the Python reference RHS actually asks for it (small-I runs/tests).
        phi_max = float(re.get('phi_max_junc',  0.5))
        sigma_s = float(re.get('sigma_s_junc',  0.35))
        # Tied to the mobility cutoff: a junction partner must be mobile, so a
        # threshold above i_mobile is unsatisfiable and silently kills the
        # channel.  Unchanged at i_mobile >= 50 (see effective_n_j_min).
        n_j_min = effective_n_j_min(re.get('n_j_min_junc', 30.0), i_mobile,
                                    re.get('n_j_min_frac', 0.6))
        # Marian two-step success probability (½⟨111⟩ → ½⟨110⟩ → ⟨100⟩, Fig. 3):
        # from the metastable ½⟨110⟩ intermediate the segment either rotates
        # FORWARD to ⟨100⟩ (barrier ΔH₂) or REVERTS to ½⟨111⟩ (reverse barrier
        # ΔH_rev = peak₁ − E_⟨110⟩).  The branching is the probability that a
        # junction / absorption event yields a *stable* ⟨100⟩ rather than
        # reverting; ΔH₂ ≫ ΔH_rev ⇒ small at low T, rising with T.  This is the
        # temperature dependence the kinetic Marian channels were missing — it
        # multiplies BOTH the junction yield and the absorption rate.
        dH2_conv    = float(re.get('dH2_conv',    1.0))   # ½⟨110⟩→⟨100⟩ barrier [eV]
        dH_rev_conv = float(re.get('dH_rev_conv', 0.30))  # ½⟨110⟩→½⟨111⟩ barrier [eV]
        ef = np.exp(-dH2_conv / kBT)
        eb = np.exp(-dH_rev_conv / kBT)
        self.conv_psuccess = float(ef / (ef + eb))        # P_success(T) ∈ (0, 0.5)
        # ── Size-asymmetric success gate for ABSORPTION (2026-08-01) ─────────
        # The single P_success above gates two physically different events:
        #   (a) a JUNCTION between two loops of comparable size, where the
        #       ½⟨110⟩ intermediate must nucleate a new ⟨100⟩ character and
        #       reversion is genuinely competitive; and
        #   (b) ABSORPTION of a mobile ≤ i_mobile cluster by a ⟨100⟩ loop of
        #       hundreds of atoms, where the incorporated segment is negligible
        #       against the host and inherits — is templated by — the host's
        #       existing ⟨100⟩ character.  Reversion should be far less likely.
        # Charging (b) the same 0.70 eV as (a) throttles it by 1/P_success =
        # 4.6e5.  That is the whole ⟨100⟩ size deficit: the ½⟨111⟩ and ⟨100⟩
        # MONOMER kernels are identical (K_SIA_grow vs K_100_grow), so the only
        # asymmetry between the characters is that a sessile ½⟨111⟩ loop eats
        # mobile clusters at full rate (K_SIA_loop) while a ⟨100⟩ loop does so
        # at 2.18e-6 of it.  Hence a separate, lower forward barrier for (b).
        #   P_success^abs(T) = A_abs · [1 + exp((ΔH₂^abs − ΔH_rev)/k_BT)]^{-1}
        # Defaults reproduce legacy behaviour exactly (ΔH₂^abs = ΔH₂, A = 1).
        # Expressing it as an energy rather than a bare rate multiplier keeps
        # the Arrhenius T-dependence the digital-twin campaign needs across
        # 300–400 °C, which a constant boost factor would not supply.
        # `_num`, not float(): 'dH2_abs_conv' ships BLANK in the workbook (its
        # default is dynamic — it tracks dH2_conv), and a blank cell arrives as
        # pandas NaN rather than as an absent key, so re.get(k, default) returns
        # NaN and the default never fires.  float(NaN) would then poison the
        # gate silently.  Same trap documented for loop_net_n_inc/loop_net_w_c.
        dH2_abs = float(_num(re.get('dH2_abs_conv', dH2_conv), dH2_conv))
        A_abs   = float(_num(re.get('psucc_abs_pref', 1.0), 1.0))
        ef_abs  = np.exp(-dH2_abs / kBT)
        self.conv_psuccess_abs = float(A_abs * ef_abs / (ef_abs + eb))
        self.dH2_abs_conv      = dH2_abs
        self.psucc_abs_pref    = A_abs
        # Parameters for the lazy φ_junc property (no dense allocation here).
        self._phi_junc_params = (phi_max, sigma_s, n_j_min)
        self._phi_junc_n      = int(ns.size)
        self._phi_junc_cache  = None

        # (3) Sessile ⟨100⟩ point-defect kernels — loop geometry, immobile loop;
        #     the migrating monomer's diffusivity sets the capture rate.  ⟨100⟩
        #     loops do not migrate, so there is no fixed-sink loss (k2_100 = 0).
        K100_grow   = np.zeros(I)
        K100_shrink = np.zeros(I)
        G100        = np.zeros(I)
        # ⟨100⟩-ONLY monomer-capture scale.  `Z_i_loop` cannot serve this role:
        # it multiplies the ½⟨111⟩ loop-capture kernels too (see the K_loop
        # branches above), so raising it grows both characters equally and
        # leaves the ⟨100⟩/½⟨111⟩ ratio unchanged — which is why a 1.05→1.30
        # sweep produced no differential effect (2026-08-01).  Default 1.0 ⇒
        # bit-identical.  Reaches the C++ solver automatically because
        # cpp_bridge exports K_100_grow as an explicit per-size array.
        grow_boost_100 = float(re.get('grow_boost_100', 1.0) or 1.0)
        for n in range(1, I + 1):
            if n >= n_loop_min:
                rt = float(n) ** 0.5
                K100_grow[n - 1]   = (A_loop * rt * Z_i_loop * Di_eff
                                      * inv_Omega23 * grow_boost_100)
                K100_shrink[n - 1] = A_loop * rt * Dv_eff * inv_Omega23
            if n >= 2:
                Eb100 = E_b_loop_100(n, A_100=A_100, B_100=B_100,
                                     n_tr=n_tr, sigma_tr=sigma_tr,
                                     E_f_i=E_f_i, G_shear=mu_Pa, b_100=a_m,
                                     nu=nu_pois, gamma_sf=gamma_sf, Omega=Omega)
                G100[n - 1] = (A_loop * max(n - 1.0, 0.0) ** 0.5 * Di_eff
                               * np.exp(-Eb100 / kBT) * inv_Omega23)
        self.K_100_grow   = K100_grow
        self.K_100_shrink = K100_shrink
        self.G_100        = G100
        self.k2_100       = np.zeros(I)   # sessile: no fixed-sink loss
        self.n_loop_min   = n_loop_min

        # ── Loop → network-dislocation loss channel (loop_network_loss.tex) ──
        # Transfer of SIA loops to the pre-existing network as a *network-only*,
        # diffusivity-independent diagonal sink Λ_n^net (active even for sessile
        # loops, which dominate incorporation).  Both characters are treated;
        # the channel is folded additively into the P4 diagonal sinks
        # (k2_SIA for ½⟨111⟩, k2_100 for ⟨100⟩) per the additive-sink-strength
        # rule, so the existing GraphWalker SINK edges and the C++ k2_SIA path
        # carry it with no new edge/term, and the SIA-content ledger
        # J_SIA_fixed (which already sums k2_SIA) stays exactly conservative.
        #   Λ_n^net = ν_net · P_ℓd(n),   ν_net = v_net · ρ_net · w_c
        # See docs/Formulation/loop_network_loss.tex Eqs. (loop_diameter_from_n,
        # loop_network_spacing, elastic_interaction_zone, P_loop_dislocation,
        # Lambda_network, vnet).
        ns_f      = np.arange(1.0, I + 1.0)
        b_100_val = a_m   # ⟨100⟩ Burgers magnitude ≈ a (cf. E_b_loop_100, b_100=a_m)
        # Per-character loop diameter d_n = 2·sqrt(n·Ω/(π·b_c))  (Eq. loop_diameter_from_n)
        self.d_loop_111 = 2.0 * np.sqrt(ns_f * Omega / (np.pi * b_111))
        self.d_loop_100 = 2.0 * np.sqrt(ns_f * Omega / (np.pi * b_100_val))
        self.Lambda_net_111 = np.zeros(I)
        self.Lambda_net_100 = np.zeros(I)
        self.loop_network_loss = self._loop_net_on
        self.rho_net = rho_d          # = static rho_d when channel off
        self.K_rec   = 0.0
        if self._loop_net_on:
            # Blank workbook cells fall back to the code default: the keys with
            # a *static* default carry an explicit number in the sheet, while
            # 'loop_net_n_inc' and 'loop_net_w_c' have *dynamic* defaults
            # (i_mobile, and the per-character Burgers vector) and are shipped
            # blank.  Reading them with float()/int() directly would turn a
            # blank cell into NaN and propagate it into Λ_n^net.
            chi   = _num(re.get('loop_net_chi',  1.0), 1.0)   # geometric range (Eq. 4)
            xi    = _num(re.get('loop_net_xi',   0.0), 0.0)   # small-n floor (default off)
            n_inc = int(_num(re.get('loop_net_n_inc', i_mobile), i_mobile))  # onset
            # Recovery prefactor K_rec [m·s^-1] sets the steady-state ρ_net; it
            # is a primary CALIBRATION parameter (Phase 5) and defaults to 0
            # (no recovery → ρ_net grows monotonically) so an uncalibrated run
            # does not crash ρ_net to the floor.  See loop_network_loss.tex.
            self.K_rec = _num(re.get('loop_net_K_rec', 0.0), 0.0)
            rho_net = rho_d
            # Segment-frozen point-defect monomers set the network climb velocity
            # v_net (Eq. vnet); held constant within a segment (operator split),
            # so Λ_n^net is a pure size array — no Jacobian coupling to c_1.
            #
            # Climb velocity (Bullough–Newman): the number of defects absorbed
            # per unit dislocation line length per time is Z_α a² ω_α c_α / Ω,
            # each adding volume Ω and climbing the line by Ω/b, so
            #   v_net = (a²/b)·(Z_i ω_i^eff c_i − Z_v ω_v^eff c_v)   [m/s].
            # (The note's Eq. vnet writes Ω/b, which is dimensionally m²/s; the
            # correct prefactor is a²/b — Ω≈a³.)  |v_net| is used: the network
            # sweeps past stationary loops whichever way it climbs.
            ci1 = float(re.get('ci1_seg', 0.0))
            cv1 = float(re.get('cv1_seg', 0.0))
            v_net = (a_m ** 2 / b_111) * (Z_i * omega_i * ci1 - Z_v * omega_v * cv1)
            v_net = abs(v_net)
            L_ld  = rho_net ** -0.5 if rho_net > 0.0 else np.inf

            def _lambda_net(d_arr, b_c, w_c):
                R_int = chi * d_arr
                if xi != 0.0:
                    R_int = R_int + xi * b_c * np.sqrt(
                        mu_Pa * Omega * _J_eV * ns_f / max(kBT, 1e-30))
                s_ld = np.maximum(b_c, L_ld - 0.5 * d_arr)
                P_ld = 0.5 * (1.0 + np.tanh((R_int - s_ld) / b_c))
                Lam  = (v_net * rho_net * w_c) * P_ld
                if n_inc > 1:
                    Lam[:min(n_inc - 1, I)] = 0.0   # loops below n_inc excluded
                return Lam

            # Capture width w_c = O(b_c) per character (fixed-but-tunable).  An
            # explicit 'loop_net_w_c' overrides BOTH characters (used to amplify
            # the channel in tests); otherwise each uses its own Burgers vector.
            w_c_over = _opt_pos_float(re.get('loop_net_w_c', None))
            w_c_111 = w_c_over if w_c_over is not None else b_111
            w_c_100 = w_c_over if w_c_over is not None else b_100_val
            self.Lambda_net_111 = _lambda_net(self.d_loop_111, b_111, w_c_111)
            self.Lambda_net_100 = _lambda_net(self.d_loop_100, b_100_val, w_c_100)
            # Additive-sink-strength rule: total loss diagonal = P4 + Λ_n^net.
            self.k2_SIA = self.k2_SIA + self.Lambda_net_111
            self.k2_100 = self.k2_100 + self.Lambda_net_100
            self.rho_net = rho_net

        # ── Cavity → network-dislocation SWEEPING  (VOID_NETWORK_LOSS) ────────
        # A climbing/gliding network dislocation that intersects a cavity
        # absorbs it: the cavity's m vacancies are delivered to the line and the
        # cavity is destroyed.  Same geometry as the loop channel above
        #     Λ_m^void = |v_net| · ρ_net · w_c · P_ld(m)
        # but applied to the VACANCY block, which makes it the only cavity sink
        # in the model that does not act through c_v.  That matters: c_v sets the
        # loop drive Z_i^loop ω_i c_i − ω_v c_v, which sits at −17% of its
        # vacancy term, so every lever that moves c_v enough to shrink cavities
        # flips the loops first (measured: Z_v 1.0→1.15 gives S_I ×65).  This
        # edge removes cavity content while leaving c_v, and hence the loops,
        # untouched.
        #
        # It is also the only channel that escapes the swelling-identity floor.
        # With loops and the bias flux as the only routes, S = S_I + ΔJ^d pins
        # the cavity content to the loop content and bounds d_cavity below by
        # 4.70 nm.  Sweeping routes cavity vacancies straight to the line, so it
        # enters ΔJ^d and S may fall below S_I.
        #
        # Λ ∝ d_cav ∝ m^(1/3) removes LARGE cavities preferentially, which
        # truncates the upper tail — the size-dependent negative feedback the
        # model otherwise lacks (growth is dm/dt ∝ m^(1/3)·drive with a drive
        # that never decays, giving unbounded m ∝ t^1.5).
        #
        # Gated OFF by default: with the key absent Λ_m^void ≡ 0 and every
        # downstream expression is unchanged, so the legacy path is
        # bit-identical.
        self.void_network_loss = int(_num(re.get('VOID_NETWORK_LOSS', 0), 0)) != 0
        self.Lambda_net_void = np.zeros(V)
        if self.void_network_loss:
            v_chi = _num(re.get('void_net_chi', 1.0), 1.0)   # elastic range / d
            # Onset: mobile vacancy clusters DIFFUSE to the network and are
            # already counted in k2_disl_v.  Sweeping them too would double-count
            # that loss, so the channel starts above v_mobile.
            m_inc = int(_num(re.get('void_net_m_inc', v_mobile), v_mobile))
            rho_n = self.rho_net
            # Same Bullough–Newman climb velocity as the loop channel; recomputed
            # here so the cavity edge does not require LOOP_NETWORK_LOSS to be on.
            ci1_v = float(re.get('ci1_seg', 0.0))
            cv1_v = float(re.get('cv1_seg', 0.0))
            v_net_v = abs((a_m ** 2 / b_111)
                          * (Z_i * omega_i * ci1_v - Z_v * omega_v * cv1_v))
            self.L_ld_void = rho_n ** -0.5 if rho_n > 0.0 else np.inf  # diagnostic
            # Spherical cavity diameter d_m = 2(3mΩ/4π)^(1/3)  — the same
            # relation post_process/run_ensemble use for d_cavity_nm.
            d_cav = 2.0 * (3.0 * ms * Omega / (4.0 * np.pi)) ** (1.0 / 3.0)
            # Capture width: blank ⇒ the cavity's own diameter (bare geometric
            # cross-section).  An explicit void_net_w_c overrides it.
            w_c_v = _opt_pos_float(re.get('void_net_w_c', None))
            w_arr = d_cav if w_c_v is None else np.full(V, w_c_v)
            # NO geometric gate here, unlike the loop channel.  P_ld there asks
            # whether a loop lies INSIDE the instantaneous elastic zone of a
            # line, which is right for incorporation of a static loop.  A
            # SWEEPING dislocation traverses the volume, so v·ρ·w is already the
            # encounter rate and a gate would double-count it — with the loop
            # gate applied, spacing L_ld = 44.7 nm against a 5.5 nm cavity gives
            # P_ld = 0 and the channel vanishes at χ = 1, which is wrong.
            # χ instead widens the CAPTURE WIDTH (elastic attraction draws
            # cavities onto the passing line): w_c = χ·d_cav.
            Lam_v = (v_net_v * rho_n) * (v_chi * w_arr)
            # Zero sizes 1..m_inc INCLUSIVE: with m_inc = v_mobile those clusters
            # diffuse to the network and are already charged to k2_disl_v, so
            # sweeping them as well would count the same loss twice.
            if m_inc >= 1:
                Lam_v[:min(m_inc, V)] = 0.0
            self.Lambda_net_void = Lam_v

        # Scalar physics
        self.B_rot  = B_rot
        self.L_hat  = L_hat
        self.alpha_He = alpha_He

    @property
    def phi_junc(self):
        """Marian junction-branching matrix φ_junc[n−1, n′−1], built lazily.

        φ = P_success · φ_max · exp(−(ln(n/n'))² / 2σ_s²) · Θ(min(n,n') ≥ n_j_min).

        This is a dense [I, I] array consumed only by the Python GraphWalker
        (the reference RHS) and the loop-conversion unit tests.  The production
        C++ solver computes the same branching on the fly per (n,n') pair, so
        the matrix is built only on first access and never on a C++ run —
        avoiding an ~O(I²) allocation that is infeasible at production sizes.
        """
        if self._phi_junc_cache is not None:
            return self._phi_junc_cache
        phi_max, sigma_s, n_j_min = self._phi_junc_params
        n = np.arange(1, self._phi_junc_n + 1, dtype=float)
        ni = n[:, None]
        nj = n[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio = np.log(ni / nj)
        phi = phi_max * np.exp(-(log_ratio ** 2) / (2.0 * sigma_s ** 2))
        phi = phi * (np.minimum(ni, nj) >= n_j_min)
        self._phi_junc_cache = phi * self.conv_psuccess
        return self._phi_junc_cache

    def __getstate__(self):
        # Nested-function attributes from _precompute can't be pickled.
        # They are reconstructable but unused after unpickling (replot path
        # only reads data arrays), so drop them.  The lazily-built φ_junc
        # cache is also dropped: it is large (O(I²)) and reconstructable.
        state = self.__dict__.copy()
        for k in ('alpha_bubble_fn', 'alpha_He_emit_fn', 'K_1D_eff_fn',
                  '_phi_junc_cache'):
            state.pop(k, None)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # Restore the dropped lazy cache so the property can rebuild on demand.
        if not hasattr(self, '_phi_junc_cache'):
            self._phi_junc_cache = None

    def loop_network_drho_dt(self, c_111, c_100=None):
        """Operator-split network-density rate dρ_net/dt (loop_network_loss.tex,
        Eq. rho_net_balance):

            dρ_net/dt = (π/Ω) Σ_c Σ_n d_n^(c) Λ_n^(c) c_n^(c)  −  K_rec ρ_net^{3/2}

        ``c_111`` / ``c_100`` are atom-fraction loop distributions (size arrays);
        ``c_100`` may be None when the ⟨100⟩ population is absent.  Returns 0.0
        when the channel is off.  Called by simulation.run_adaptive between
        segments to advance ρ_net explicitly while it is held frozen inside the
        stiff cluster solve.
        """
        if not getattr(self, 'loop_network_loss', False):
            return 0.0
        Omega = float(self.inp.derived['Omega'])
        c_111 = np.asarray(c_111, dtype=float)
        gain = np.dot(self.d_loop_111[:c_111.size] * self.Lambda_net_111[:c_111.size],
                      np.maximum(c_111, 0.0))
        if c_100 is not None:
            c_100 = np.asarray(c_100, dtype=float)
            gain += np.dot(self.d_loop_100[:c_100.size] * self.Lambda_net_100[:c_100.size],
                           np.maximum(c_100, 0.0))
        gain *= np.pi / Omega
        recovery = self.K_rec * max(self.rho_net, 0.0) ** 1.5
        return float(gain - recovery)

    def format_diagnostic(self, mean_n_i=None):
        """Return key rate constants as a formatted string.

        Parameters
        ----------
        mean_n_i : float, optional
            Mean SIA cluster size <n> at the current output step.
        """
        mean_str = f"  mean_n_i={mean_n_i:.2f}" if mean_n_i is not None else ""
        lines = []
        lines.append(
            f"ReactionRates: K_SIA_grow[0]={self.K_SIA_grow[0]:.3e}"
            f"  K_VAC_grow[0]={self.K_VAC_grow[0]:.3e}"
            f"  G_VAC[0]={self.G_VAC[0]:.3e}"
            f"  K_iv={self.K_iv:.3e}"
            f"  K_3D_cav={self.K_3D_cav_pref:.3e}"
            f"{mean_str}")
        lines.append(
            f"  k2_SIA[0]={self.k2_SIA[0]:.3e}"
            f"  k2_vac={self.k2_vac_scalar:.3e}"
            f"  k2_He={self.k2_He_scalar:.3e}")
        lines.append(
            f"  D_SIA_eff: n=1: {self.D_SIA_eff[0]:.3e}"
            f"  n=2: {self.D_SIA_eff[1]:.3e}"
            f"  n=3: {self.D_SIA_eff[2]:.3e}"
            + (f"  n=4: {self.D_SIA_eff[3]:.3e}" if len(self.D_SIA_eff) > 3 else ""))
        if len(self.K_SIA_grow) >= 5:
            lines.append(
                f"  C_i5: K_grow={self.K_SIA_grow[4]:.3e}"
                f"  K_shrink={self.K_SIA_shrink[4]:.3e}"
                f"  K_loop={self.K_SIA_loop[4]:.3e}"
                f"  k2={self.k2_SIA[4]:.3e}"
                f"  G_emit={self.G_SIA[4]:.3e}")
        return '\n'.join(lines)

    def print_diagnostic(self, mean_n_i=None):
        """Print key rate constants — called at output time steps."""
        print(self.format_diagnostic(mean_n_i))
