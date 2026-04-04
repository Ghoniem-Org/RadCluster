"""
bin_moment_rates.py — Size-bin moment reduction for Expanded_Eurofer_CD.

Implements the Chapter 9 state-space reduction for the SIA cluster population
via logarithmic size bins and zeroth/first moments.

Physics reference
-----------------
Ghoniem, N.M. (2026), Chapter 9 (Rate_Equations.pdf):
  Eqs. 188-211.

Bin partition (Eq. 188)
-----------------------
Logarithmic bin edges:
  n_k = floor(n_1 · r^k),  k = 0, 1, ..., K_bins
  r > 1  (typically r = 2)

Bin B_k = {n : n_k ≤ n < n_{k+1}}  (Eq. 189)
Bin width |B_k| = n_{k+1} − n_k    (Eq. 190)

Bin moments (Eq. 192)
---------------------
  μ_k^(0) = Σ_{n ∈ B_k} c_n        (zeroth — bin density)
  μ_k^(1) = Σ_{n ∈ B_k} n · c_n    (first  — bin content)

Binned master equation (Eq. 193-197)
-------------------------------------
  dμ_k^(q)/dt = Σ_{n ∈ B_k} n^q · (dc_n/dt)_reaction

Piecewise-constant closure (Eq. 198-200)
-----------------------------------------
  c_n ≈ μ_k^(0) / |B_k|   for all n ∈ B_k

Hat-function (Galerkin) closure (Eq. 201-206)
----------------------------------------------
  c_n = φ_{k,0}(n) · μ_k^(0) + φ_{k,1}(n) · μ_k^(1)
where
  φ_{k,0}(n) = (n_{k+1} − n) / (|B_k| · (|B_k|+1)/2)   (Eq. 203)
  φ_{k,1}(n) = (n − n_k)    / (|B_k| · (|B_k|+1)/2)   (Eq. 204)

Inter-bin upwind flux (Eq. 207-208)
-------------------------------------
  F_{k→k+1} = (growth_rate at n_{k+1}) · c_{n_{k+1}}^{upwind}
  Upwind: use bin k value if growth net positive, bin k+1 if negative.

Conservation diagnostic (Eq. 211)
------------------------------------
  δ_FP^bin = |Σ_k μ_k^(1) + Σ_m m·c_m − ΔJ^d| / (G·t)
"""

import numpy as np
from .defect_production import production_rates


_kB = 8.617333262e-5


def build_bins(N_max, n1=1, r=2.0):
    """
    Build logarithmic bin partition.

    Parameters
    ----------
    N_max : int  — maximum SIA cluster size tracked
    n1    : int  — minimum cluster size (first bin starts here)
    r     : float — bin ratio n_{k+1}/n_k > 1

    Returns
    -------
    bins : list of (n_lo, n_hi) tuples
        Each tuple gives the inclusive lower and exclusive upper bounds of bin k.
        bins[k] = (n_k, n_{k+1})
    edges : ndarray of int bin edges
    """
    edges = [n1]
    while edges[-1] < N_max:
        next_edge = max(int(np.floor(edges[-1] * r)), edges[-1] + 1)
        edges.append(min(next_edge, N_max + 1))
    edges = np.array(edges, dtype=int)
    bins  = [(int(edges[k]), int(edges[k + 1])) for k in range(len(edges) - 1)]
    return bins, edges


def moments_from_distribution(c_n, bins):
    """
    Compute bin moments μ_k^(0) and μ_k^(1) from per-size concentrations.

    Parameters
    ----------
    c_n  : ndarray [N_max]  — concentrations, c_n[n-1] = c_n
    bins : list of (n_lo, n_hi)

    Returns
    -------
    mu0 : ndarray [K_bins]  — zeroth moments
    mu1 : ndarray [K_bins]  — first moments
    """
    K = len(bins)
    mu0 = np.zeros(K)
    mu1 = np.zeros(K)
    for k, (nlo, nhi) in enumerate(bins):
        ns  = np.arange(nlo, nhi)
        idx = ns - 1
        valid = (idx >= 0) & (idx < len(c_n))
        mu0[k] = np.sum(c_n[idx[valid]])
        mu1[k] = np.sum(ns[valid] * c_n[idx[valid]])
    return mu0, mu1


def distribution_from_moments_pc(mu0, mu1, bins, N_max):
    """
    Reconstruct piecewise-constant per-size concentrations from moments.

    c_n = μ_k^(0) / |B_k|   for n ∈ B_k   (Eq. 198-200)

    Parameters
    ----------
    mu0  : ndarray [K_bins]
    mu1  : ndarray [K_bins]  (unused in piecewise-constant)
    bins : list of (n_lo, n_hi)
    N_max: int

    Returns
    -------
    c_n : ndarray [N_max]
    """
    c_n = np.zeros(N_max)
    for k, (nlo, nhi) in enumerate(bins):
        bw = float(nhi - nlo)
        val = mu0[k] / bw if bw > 0 else 0.0
        ns  = np.arange(nlo, nhi) - 1
        valid = (ns >= 0) & (ns < N_max)
        c_n[ns[valid]] = val
    return c_n


def distribution_from_moments_hat(mu0, mu1, bins, N_max):
    """
    Reconstruct hat-function concentrations from moments.

    c_n = φ_{k,0}(n)·μ_k^(0) + φ_{k,1}(n)·μ_k^(1)   (Eqs. 201-206)

    Parameters
    ----------
    mu0, mu1 : ndarray [K_bins]
    bins     : list of (n_lo, n_hi)
    N_max    : int

    Returns
    -------
    c_n : ndarray [N_max]
    """
    c_n = np.zeros(N_max)
    for k, (nlo, nhi) in enumerate(bins):
        ns  = np.arange(nlo, nhi, dtype=float)
        bw  = float(nhi - nlo)
        if bw <= 0:
            continue
        norm = bw * (bw + 1.0) / 2.0
        phi0 = (float(nhi) - ns) / norm    # Eq. 203
        phi1 = (ns - float(nlo)) / norm    # Eq. 204
        vals = phi0 * mu0[k] + phi1 * mu1[k]
        idx  = (ns - 1).astype(int)
        valid = (idx >= 0) & (idx < N_max)
        c_n[idx[valid]] = np.maximum(vals[valid], 0.0)
    return c_n


class BinMomentRateEquations:
    """
    Chapter 9 bin-moment ODE system for SIA clusters.

    The vacancy cluster and He populations use the same equations as
    full_CD (either Case 1 or Case 2 He-reduction).

    State vector:
      y[0..2*K-1]        : SIA bin moments [μ_0^(0), μ_0^(1), ..., μ_{K-1}^(0), μ_{K-1}^(1)]
      y[2*K..2*K+M-1]    : void/bubble c_m or c_m^tot, m=1..M
      y[2*K+M..]         : He variables (same as Case 2 or Case 1)

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

        # Build bin partition
        n1_bin  = int(float(input_data.reactions.get('n1_bin',  1)))
        r_ratio = float(input_data.reactions.get('r_ratio', 2.0))
        self.bins, self.edges = build_bins(N, n1=n1_bin, r=r_ratio)
        self.K    = len(self.bins)
        self.n1   = n1_bin
        self.r    = r_ratio

        # He mode (same dispatch as full_CD)
        po = input_data.physics_option
        if 'fusion' in po:
            self.he_mode = 'case1'
            # Indices: [SIA moments | vac | Q_m | He]
            self.i_MOM = 0
            self.i_VAC = 2 * self.K
            self.i_Q   = 2 * self.K + M
            self.i_He  = 2 * self.K + 2 * M
            self.N_eq  = 2 * self.K + 2 * M + 1
        else:
            self.he_mode = 'case2'
            self.i_MOM  = 0
            self.i_VAC  = 2 * self.K
            self.i_Qtot = 2 * self.K + M
            self.i_He   = 2 * self.K + M + 1
            self.N_eq   = 2 * self.K + M + 2

        # Cascade production
        spectrum = input_data.derived['spectrum']
        G        = input_data.derived['G']
        Pr_SIA, Pr_VAC, G_He = production_rates(G, spectrum, N, M)
        self.Pr_SIA = Pr_SIA[1:]   # [N], index k → size k+1
        self.Pr_VAC = Pr_VAC[1:]   # [M]
        self.G_He   = G_He

        print(f"BinMomentRateEquations: K_bins={self.K}  N_eq={self.N_eq}"
              f"  he_mode='{self.he_mode}'  r={r_ratio}")
        print(f"  Bin edges (first 5): {list(self.edges[:6])}")

    def get_initial_conditions(self):
        """Return y0 — near-zero moments + concentrations."""
        C_floor = float(self.inp.reactions.get('C_floor', 1e-100))
        return np.full(self.N_eq, C_floor)

    def ode_system(self, t, y):
        """ODE RHS using piecewise-constant closure (Eq. 198-200)."""
        K   = self.K
        M   = self.M
        rr  = self.rr
        inp = self.inp
        dydt = np.zeros(self.N_eq)

        # Unpack moments: y[2k] = μ_k^(0), y[2k+1] = μ_k^(1)
        mu0 = np.maximum(y[0::2][:K], 0.0)   # [K]
        mu1 = np.maximum(y[1::2][:K], 0.0)   # [K]

        # Reconstruct per-size distribution (piecewise-constant, Eq. 198-200)
        c_n = distribution_from_moments_pc(mu0, mu1, self.bins, self.N)

        # Vacancy and He state
        if self.he_mode == 'case2':
            c_v   = np.maximum(y[self.i_VAC:self.i_Qtot], 0.0)
            Q_tot = max(y[self.i_Qtot], 0.0)
            c_h   = max(y[self.i_He], 0.0)
        else:
            c_v   = np.maximum(y[self.i_VAC:self.i_Q], 0.0)
            Q_m   = np.maximum(y[self.i_Q:self.i_He], 0.0)
            c_h   = max(y[self.i_He], 0.0)
            Q_tot = np.sum(Q_m)

        ci1 = c_n[0] if len(c_n) > 0 else 0.0
        cv1 = c_v[0]

        # Compute full dc_n/dt for each size (reuse rate terms)
        dc_n = np.zeros(self.N)

        # Cascade production (Eq. 152 source)
        dc_n += self.Pr_SIA

        # SIA emission from n+1
        dc_n[:-1] += rr.G_SIA[1:] * c_n[1:]
        # Loss by emission
        dc_n -= rr.G_SIA * c_n
        # SIA capture: growth n-1→n
        dc_n[1:] += rr.K_SIA_grow[:-1] * ci1 * c_n[:-1]
        # Loss by capture (growth n→n+1)
        dc_n -= rr.K_SIA_grow * ci1 * c_n
        # Vacancy annihilation
        dc_n -= rr.K_SIA_shrink * cv1 * c_n
        dc_n[:-1] += rr.K_SIA_shrink[1:] * cv1 * c_n[1:]
        # Fixed sinks
        dc_n -= rr.k2_SIA * c_n

        # Project dc_n onto bin-moment equations (Eqs. 193-197)
        # dμ_k^(0)/dt = Σ_{n∈B_k} dc_n
        # dμ_k^(1)/dt = Σ_{n∈B_k} n · dc_n
        for k, (nlo, nhi) in enumerate(self.bins):
            ns  = np.arange(nlo, nhi)
            idx = ns - 1
            valid = (idx >= 0) & (idx < self.N)
            dydt[2 * k]     = np.sum(dc_n[idx[valid]])
            dydt[2 * k + 1] = np.sum(ns[valid] * dc_n[idx[valid]])

        # Add inter-bin upwind flux (Eq. 207-208)
        # For each bin boundary, transfer from k to k+1 if net growth positive
        n_max_i = inp.derived['n_max_i']
        for k in range(K - 1):
            nlo2, nhi2 = self.bins[k + 1]
            n_edge = nlo2            # boundary size
            if n_edge >= self.N:
                break
            # Growth rate at boundary: K_grow · ci1
            g_rate = rr.K_SIA_grow[n_edge - 1] * ci1 if n_edge <= self.N else 0.0
            c_edge = c_n[n_edge - 1]
            # Upwind: if growing, flux leaves bin k and enters bin k+1
            flux_01 = g_rate * c_edge   # [at.frac/s] leaving bin k
            flux_10 = rr.G_SIA[n_edge - 1] * c_edge   # emission from boundary back
            dydt[2 * k]     -= flux_01
            dydt[2 * k + 1] -= float(n_edge) * flux_01
            dydt[2 * (k + 1)]     += flux_01
            dydt[2 * (k + 1) + 1] += float(n_edge) * flux_01
            # Emission flux (Eq. 208 reverse)
            dydt[2 * k]     += flux_10
            dydt[2 * k + 1] += float(n_edge) * flux_10
            dydt[2 * (k + 1)]     -= flux_10
            dydt[2 * (k + 1) + 1] -= float(n_edge) * flux_10

        # ── Vacancy cluster equations (same as full_CD Case 2) ─────────────
        C_vac_tot = np.sum(c_v)
        ell_bar   = Q_tot / max(C_vac_tot, 1e-200)
        G_VAC_eff = rr.G_VAC.copy()
        if ell_bar > 0.01:
            for m in range(1, min(M, 20) + 1):
                G_VAC_eff[m - 1] = rr.alpha_bubble_fn(m, ell_bar * m**(2.0/3.0))

        dcv = dydt[self.i_VAC:self.i_VAC + M]
        dcv += self.Pr_VAC
        dcv -= G_VAC_eff * c_v
        dcv[:-1] += G_VAC_eff[1:] * c_v[1:]
        dcv -= rr.K_VAC_grow * cv1 * c_v
        dcv[1:] += rr.K_VAC_grow[:-1] * cv1 * c_v[:-1]
        dcv -= rr.K_VAC_shrink * ci1 * c_v
        dcv[:-1] += rr.K_VAC_shrink[1:] * ci1 * c_v[1:]
        dcv -= rr.K_HeV * c_h * c_v
        dcv -= rr.k2_vac_scalar * c_v

        # ── He equations ────────────────────────────────────────────────────
        if self.he_mode == 'case2':
            He_uptake  = np.sum(rr.K_HeV * c_h * c_v)
            He_release = 0.0
            dydt[self.i_Qtot] = He_uptake - He_release
            dydt[self.i_He]   = (self.G_He
                                 - np.sum(rr.K_HeV * c_v * c_h)
                                 - rr.k2_He_scalar * c_h)
        else:
            dQ = dydt[self.i_Q:self.i_He]
            for m in range(1, M + 1):
                dQ[m - 1] += rr.K_HeV[m - 1] * c_h * c_v[m - 1]
                ell_m = Q_m[m - 1] / max(c_v[m - 1], 1e-200)
                alpha_h = rr.alpha_He_emit_fn(m, max(int(round(ell_m)), 0))
                dQ[m - 1] -= alpha_h * Q_m[m - 1]
            dydt[self.i_He] = (self.G_He
                               - np.sum(rr.K_HeV * c_v * c_h)
                               - rr.k2_He_scalar * c_h)

        return dydt

    def conservation_diagnostic(self, t, y, G_tot):
        """
        FP conservation diagnostic for bin-moment system (Eq. 211).

        δ_FP^bin = |Σ_k μ_k^(1) + Σ_m m·c_m − ΔJ^d| / (G·t)
        """
        K = self.K
        M = self.M
        mu1 = np.maximum(y[1::2][:K], 0.0)
        c_v = np.maximum(y[self.i_VAC:self.i_VAC + M], 0.0)
        ms  = np.arange(1, M + 1, dtype=float)

        SIA_content = np.sum(mu1)
        VAC_content = np.sum(ms * c_v)
        delta_J     = abs(SIA_content - VAC_content)
        denom       = max(G_tot * max(t, 1e-20), 1e-300)
        return delta_J / denom
