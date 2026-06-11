"""
bin_moment_rates.py — Size-bin moment reduction for RadCluster_2_0.

REFERENCE ONLY (ODE RHS) — the `ode_system` method of
BinMomentRateEquations is not used by the solver; the production RHS is the
C++ implementation.  May be stale.  Scheduled for replacement by the
Stage-3 graph-walker.  simulation.py integrates exclusively via the C++
solver and never calls the Python `ode_system`.  Note that the
reconstruction/projection helpers and bin-construction code in this module
ARE used (e.g. by simulation.py's adaptive domain expansion); only the
`ode_system` RHS is dead code.

Implements the Chapter 9 state-space reduction for the SIA and vacancy
cluster populations via logarithmic size bins.  Three intra-bin shape
functions (closures) are supported:

  shape_function = "constant"   — piecewise-constant (1 moment/bin: μ₀)
  shape_function = "linear"     — hat-function / dual-basis (2 moments/bin: μ₀, μ₁)
  shape_function = "lognormal"  — log-normal (3 moments/bin: μ₀, μ₁, μ₂)

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
  μ_k^(2) = Σ_{n ∈ B_k} n² · c_n   (second — used only by lognormal)

Closures
--------
Piecewise-constant (Eq. 198-200):
  c_n ≈ μ_k^(0) / |B_k|   for all n ∈ B_k

Hat-function / linear (Galerkin) (Eq. 201-206):
  c_n = φ_{k,0}(n) · μ_k^(0) + φ_{k,1}(n) · μ_k^(1)

Log-normal (Eq. 213-216):
  c_n ∝ (1/n) exp[−(ln n − m_k)² / (2σ_k²)]
  with m_k, σ_k determined from μ₀, μ₁, μ₂; normalized to preserve μ₀.
"""

import numpy as np
from .defect_production import production_rates


_kB = 8.617333262e-5

# Valid shape function names and the number of moments each requires per bin
_SHAPE_FUNCTIONS = {'constant': 1, 'linear': 2, 'lognormal': 3}

# Integer encoding for C++ bridge: constant=0, linear=1, lognormal=2
_SHAPE_FUNCTION_INT = {'constant': 0, 'linear': 1, 'lognormal': 2}


def n_moments_for_shape(shape_function):
    """Return the number of tracked moments per bin (1, 2, or 3)."""
    if shape_function not in _SHAPE_FUNCTIONS:
        raise ValueError(
            f"Unknown shape_function='{shape_function}'. "
            f"Must be one of {list(_SHAPE_FUNCTIONS.keys())}.")
    return _SHAPE_FUNCTIONS[shape_function]


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


def moments_from_distribution(c_n, bins, n_mom=2):
    """
    Compute bin moments from per-size concentrations.

    Parameters
    ----------
    c_n   : ndarray [N_max]  — concentrations, c_n[n-1] = c_n
    bins  : list of (n_lo, n_hi)
    n_mom : int — number of moments to compute (1, 2, or 3)

    Returns
    -------
    mu0 : ndarray [K_bins]  — zeroth moments
    mu1 : ndarray [K_bins] or None  — first moments (if n_mom >= 2)
    mu2 : ndarray [K_bins] or None  — second moments (if n_mom >= 3)
    """
    K = len(bins)
    mu0 = np.zeros(K)
    mu1 = np.zeros(K) if n_mom >= 2 else None
    mu2 = np.zeros(K) if n_mom >= 3 else None
    for k, (nlo, nhi) in enumerate(bins):
        ns  = np.arange(nlo, nhi)
        idx = ns - 1
        valid = (idx >= 0) & (idx < len(c_n))
        c_valid = c_n[idx[valid]]
        n_valid = ns[valid]
        mu0[k] = np.sum(c_valid)
        if mu1 is not None:
            mu1[k] = np.sum(n_valid * c_valid)
        if mu2 is not None:
            mu2[k] = np.sum(n_valid * n_valid * c_valid)
    return mu0, mu1, mu2


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


def midpoint_distribution_from_moments(mu0, mu1, bins, mu2=None,
                                       shape_function='linear'):
    """
    Evaluate concentration at bin geometric midpoints using the chosen closure.

    Parameters
    ----------
    mu0  : ndarray [K]
    mu1  : ndarray [K] or None
    bins : list of (n_lo, n_hi)
    mu2  : ndarray [K] or None  — second moments (used only by lognormal)
    shape_function : str — "constant", "linear", or "lognormal"

    Returns
    -------
    midpoints : ndarray [K]  — geometric midpoint of each bin
    conc      : ndarray [K]  — concentration at that midpoint
    """
    K = len(bins)
    midpoints = np.zeros(K)
    conc = np.zeros(K)

    for k, (nlo, nhi) in enumerate(bins):
        bw = nhi - nlo
        if bw <= 0:
            continue
        midpoints[k] = np.sqrt(float(nlo) * float(max(nhi - 1, nlo)))

        if bw == 1:
            conc[k] = max(mu0[k], 0.0)
            continue

        if shape_function == 'constant':
            conc[k] = max(mu0[k] / float(bw), 0.0)

        elif shape_function == 'lognormal' and mu2 is not None and mu0[k] > 0:
            n_bar  = mu1[k] / max(mu0[k], 1e-300)
            n2_bar = mu2[k] / max(mu0[k], 1e-300)
            ratio  = n2_bar / max(n_bar * n_bar, 1e-300)
            if ratio > 1.0 + 1e-12:
                sig2 = np.log(ratio)
                m_k  = np.log(max(n_bar, 1e-300)) - 0.5 * sig2
                n_mid = midpoints[k]
                ln_n  = np.log(n_mid)
                conc[k] = max(mu0[k] / float(bw)
                              * np.exp(-(ln_n - m_k)**2 / (2.0 * sig2)), 0.0)
            else:
                # Fallback to linear
                mu1_k = mu1[k] if mu1 is not None else mu0[k] * (nlo + nhi - 1) / 2.0
                ns = np.arange(nlo, nhi, dtype=float)
                S1, S2 = np.sum(ns), np.sum(ns * ns)
                det = bw * S2 - S1 * S1
                if abs(det) < 1e-300:
                    conc[k] = max(mu0[k] / bw, 0.0)
                else:
                    n_mid = midpoints[k]
                    phi0 = (S2 - S1 * n_mid) / det
                    phi1 = (bw * n_mid - S1) / det
                    conc[k] = max(phi0 * mu0[k] + phi1 * mu1_k, 0.0)

        else:  # 'linear' (default / hat-function)
            mu1_k = mu1[k] if mu1 is not None else mu0[k] * (nlo + nhi - 1) / 2.0
            n_mid = midpoints[k]
            ns  = np.arange(nlo, nhi, dtype=float)
            S1  = np.sum(ns)
            S2  = np.sum(ns * ns)
            det = bw * S2 - S1 * S1
            if abs(det) < 1e-300:
                conc[k] = max(mu0[k] / bw, 0.0)
            else:
                phi0 = (S2 - S1 * n_mid) / det
                phi1 = (bw * n_mid - S1) / det
                conc[k] = max(phi0 * mu0[k] + phi1 * mu1_k, 0.0)

    return midpoints, conc


def distribution_from_moments_hat(mu0, mu1, bins, N_max):
    """
    Reconstruct per-size concentrations from moments using a linear
    dual-basis that exactly preserves both μ₀ and μ₁.

    Shape functions (linear dual basis):
        φ₀(n) = (S₂ − S₁·n) / Δ
        φ₁(n) = (bw·n − S₁) / Δ
    where
        S₁ = Σ n,  S₂ = Σ n²,  Δ = bw·S₂ − S₁²

    This guarantees:
        Σ φ₀ = 1,  Σ φ₁ = 0   →  Σ cₙ = μ₀
        Σ n·φ₀ = 0, Σ n·φ₁ = 1 →  Σ n·cₙ = μ₁

    Negative values are clamped to 0 (breaks exact preservation for
    strongly skewed bins, but keeps concentrations physical).

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
        bw = float(nhi - nlo)
        if bw <= 0:
            continue
        if bw == 1:
            # Single-size bin: c_n = mu0
            idx = nlo - 1
            if 0 <= idx < N_max:
                c_n[idx] = max(mu0[k], 0.0)
            continue
        ns  = np.arange(nlo, nhi, dtype=float)
        S1  = np.sum(ns)
        S2  = np.sum(ns * ns)
        det = bw * S2 - S1 * S1
        if abs(det) < 1e-300:
            # Degenerate bin — fall back to piecewise-constant
            c_n_val = mu0[k] / bw
            idx = (ns - 1).astype(int)
            valid = (idx >= 0) & (idx < N_max)
            c_n[idx[valid]] = max(c_n_val, 0.0)
            continue
        phi0 = (S2 - S1 * ns) / det
        phi1 = (bw * ns - S1) / det
        vals = phi0 * mu0[k] + phi1 * mu1[k]
        idx  = (ns - 1).astype(int)
        valid = (idx >= 0) & (idx < N_max)
        c_n[idx[valid]] = np.maximum(vals[valid], 0.0)
    return c_n


def distribution_from_moments_lognormal(mu0, mu1, mu2, bins, N_max):
    """
    Reconstruct per-size concentrations from moments using a log-normal
    intra-bin shape function (Eq. 213-216).

    Within each bin B_k, the distribution is assumed to follow:
        c_n ∝ (1/n) exp[−(ln n − m_k)² / (2σ_k²)]

    The two log-normal parameters are determined from the moments:
        n_bar   = μ₁/μ₀
        n2_bar  = μ₂/μ₀
        σ_k²    = ln(n2_bar / n_bar²) = ln(μ₂·μ₀ / μ₁²)
        m_k     = ln(n_bar) − σ_k²/2

    The amplitude is set so that Σ_{n∈B_k} c_n = μ₀ exactly.

    Falls back to linear (hat-function) closure when σ_k² ≤ 0
    (monodisperse limit), and to piecewise-constant when the bin
    has only one size.

    Parameters
    ----------
    mu0, mu1, mu2 : ndarray [K_bins]
    bins          : list of (n_lo, n_hi)
    N_max         : int

    Returns
    -------
    c_n : ndarray [N_max]
    """
    c_n = np.zeros(N_max)
    for k, (nlo, nhi) in enumerate(bins):
        bw = nhi - nlo
        if bw <= 0 or mu0[k] <= 0.0:
            continue
        if bw == 1:
            idx = nlo - 1
            if 0 <= idx < N_max:
                c_n[idx] = max(mu0[k], 0.0)
            continue

        ns = np.arange(nlo, nhi, dtype=float)
        n_bar  = mu1[k] / max(mu0[k], 1e-300)
        n2_bar = mu2[k] / max(mu0[k], 1e-300)

        # σ² = ln(μ₂·μ₀ / μ₁²);  degenerate when ≤ 0
        ratio = n2_bar / max(n_bar * n_bar, 1e-300)
        if ratio <= 1.0 + 1e-12:
            # Monodisperse or invalid — fall back to linear (hat-function)
            S1  = np.sum(ns)
            S2  = np.sum(ns * ns)
            det = float(bw) * S2 - S1 * S1
            if abs(det) < 1e-300:
                val = mu0[k] / float(bw)
                idx = (ns - 1).astype(int)
                valid = (idx >= 0) & (idx < N_max)
                c_n[idx[valid]] = max(val, 0.0)
            else:
                phi0 = (S2 - S1 * ns) / det
                phi1 = (float(bw) * ns - S1) / det
                vals = phi0 * mu0[k] + phi1 * mu1[k]
                idx  = (ns - 1).astype(int)
                valid = (idx >= 0) & (idx < N_max)
                c_n[idx[valid]] = np.maximum(vals[valid], 0.0)
            continue

        sig2 = np.log(ratio)
        m_k  = np.log(max(n_bar, 1e-300)) - 0.5 * sig2

        # Log-normal shape: f(n) = (1/n) exp[−(ln n − m)² / (2σ²)]
        ln_ns = np.log(ns)
        log_f = -(ln_ns - m_k)**2 / (2.0 * sig2) - ln_ns
        # Subtract max for numerical stability
        log_f -= np.max(log_f)
        f = np.exp(log_f)

        # Normalize so that Σ f_n = μ₀
        f_sum = np.sum(f)
        if f_sum < 1e-300:
            val = mu0[k] / float(bw)
            idx = (ns - 1).astype(int)
            valid = (idx >= 0) & (idx < N_max)
            c_n[idx[valid]] = max(val, 0.0)
            continue

        vals = f * (mu0[k] / f_sum)
        idx  = (ns - 1).astype(int)
        valid = (idx >= 0) & (idx < N_max)
        c_n[idx[valid]] = np.maximum(vals[valid], 0.0)

    return c_n


def reconstruct_distribution(shape_function, mu0, mu1, mu2, bins, N_max,
                             *, smooth_edges=False):
    """
    Dispatch to the appropriate reconstruction function.

    Parameters
    ----------
    shape_function : str — "constant", "linear", or "lognormal"
    mu0 : ndarray [K]
    mu1 : ndarray [K] or None
    mu2 : ndarray [K] or None
    bins : list of (n_lo, n_hi)
    N_max : int
    smooth_edges : bool, default False
        If True, post-process the reconstructed distribution by averaging
        the two values straddling each internal bin boundary (geometric
        mean in log-space) and rescaling each bin so that
        Σ_{n ∈ B_k} c_n equals its pre-smoothing sum, preserving μ_k^(0)
        exactly. Use for visualization only — ODE right-hand sides and
        conservation diagnostics must pass False.

    Returns
    -------
    c_n : ndarray [N_max]
    """
    if shape_function == 'constant':
        c_n = distribution_from_moments_pc(mu0, mu1, bins, N_max)
    elif shape_function == 'linear':
        c_n = distribution_from_moments_hat(mu0, mu1, bins, N_max)
    elif shape_function == 'lognormal':
        c_n = distribution_from_moments_lognormal(mu0, mu1, mu2, bins, N_max)
    else:
        raise ValueError(f"Unknown shape_function='{shape_function}'")

    if smooth_edges:
        _smooth_bin_edges_inplace(c_n, bins)
    return c_n


def _smooth_bin_edges_inplace(c_n, bins):
    """
    Visualization post-process: kill the staircase jumps at internal bin
    boundaries while preserving each bin's μ_k^(0) exactly.

    Algorithm:
      1. Record each bin's pre-smoothing sum (= μ_k^(0) of the closure).
      2. At each internal boundary, replace both straddling values
         (last size of bin k, first size of bin k+1) with their
         geometric mean (or arithmetic mean if either is non-positive).
      3. Multiplicatively rescale each bin's interior so the sum is
         restored to step 1.

    Width-1 bins are effectively unchanged (the rescale undoes the mean).
    First/last bins keep their outer edge untouched.
    """
    K = len(bins)
    if K < 2:
        return c_n
    N_max = len(c_n)

    pre_sums = np.empty(K)
    for k, (nlo, nhi) in enumerate(bins):
        idx_lo = max(nlo - 1, 0)
        idx_hi = min(nhi - 1, N_max)
        pre_sums[k] = np.sum(c_n[idx_lo:idx_hi]) if idx_hi > idx_lo else 0.0

    for k in range(K - 1):
        n_boundary = bins[k][1]      # = bins[k+1][0]
        i_left  = n_boundary - 2     # 0-idx: last size of bin k
        i_right = n_boundary - 1     # 0-idx: first size of bin k+1
        if i_left < 0 or i_right >= N_max:
            continue
        v_l, v_r = c_n[i_left], c_n[i_right]
        if v_l > 0.0 and v_r > 0.0:
            avg = np.sqrt(v_l * v_r)
        else:
            avg = 0.5 * (v_l + v_r)
        c_n[i_left]  = avg
        c_n[i_right] = avg

    for k, (nlo, nhi) in enumerate(bins):
        idx_lo = max(nlo - 1, 0)
        idx_hi = min(nhi - 1, N_max)
        if idx_hi <= idx_lo or pre_sums[k] <= 0.0:
            continue
        new_sum = np.sum(c_n[idx_lo:idx_hi])
        if new_sum > 0.0:
            c_n[idx_lo:idx_hi] *= pre_sums[k] / new_sum
    return c_n


def distribution_from_moments_continuous(mu0, mu1, bins, N_max):
    """
    Reconstruct per-size concentrations with C⁰ continuity across bin edges.

    At each internal bin boundary, the left and right hat-function values
    are averaged.  Within each bin, linear interpolation between the
    (averaged) edge values replaces the raw hat function.  This removes
    the staircase discontinuities visible in the standard reconstruction
    while preserving the bin moments to first order.

    Parameters
    ----------
    mu0, mu1 : ndarray [K_bins]
    bins     : list of (n_lo, n_hi)
    N_max    : int

    Returns
    -------
    c_n : ndarray [N_max]
    """
    K = len(bins)
    c_n = np.zeros(N_max)

    # Step 1: compute raw hat-function value at every bin edge
    # edge_val_left[k]  = value at n_lo of bin k  (from bin k)
    # edge_val_right[k] = value at n_hi-1 of bin k (from bin k, right side)
    edge_left  = np.zeros(K)  # c(n_lo) from bin k
    edge_right = np.zeros(K)  # c(n_hi-1) from bin k (last point in bin)

    for k, (nlo, nhi) in enumerate(bins):
        bw = float(nhi - nlo)
        if bw <= 0:
            continue
        norm = bw * (bw + 1.0) / 2.0
        # Value at n = nlo (leftmost point)
        phi0_lo = float(nhi - nlo) / norm
        phi1_lo = 0.0 / norm  # (nlo - nlo) = 0
        edge_left[k] = phi0_lo * mu0[k] + phi1_lo * mu1[k]
        # Value at n = nhi - 1 (rightmost point in bin)
        phi0_hi = 1.0 / norm   # (nhi - (nhi-1)) = 1
        phi1_hi = (nhi - 1.0 - nlo) / norm
        edge_right[k] = phi0_hi * mu0[k] + phi1_hi * mu1[k]

    # Step 2: average edge values at internal boundaries
    # At boundary between bin k and k+1:
    #   edge_right[k] is bin k's value at n = n_hi[k] - 1
    #   edge_left[k+1] is bin k+1's value at n = n_lo[k+1] = n_hi[k]
    # These are adjacent sizes (n_hi-1 and n_hi), so we set:
    #   c(n_hi - 1) = average of edge_right[k] and edge_left[k+1]
    # and keep the bin's own left edge for the leftmost point.
    c_left  = np.zeros(K)  # adjusted value at bin's left edge
    c_right = np.zeros(K)  # adjusted value at bin's right edge

    for k in range(K):
        c_left[k]  = edge_left[k]
        c_right[k] = edge_right[k]

    # Smooth internal boundaries: average the right-of-k with left-of-k+1
    for k in range(K - 1):
        avg = 0.5 * (edge_right[k] + edge_left[k + 1])
        c_right[k]    = avg
        c_left[k + 1] = avg

    # Step 3: linearly interpolate within each bin between adjusted edges
    for k, (nlo, nhi) in enumerate(bins):
        bw = nhi - nlo
        if bw <= 0:
            continue
        ns = np.arange(nlo, nhi, dtype=float)
        if bw == 1:
            vals = np.array([c_left[k]])
        else:
            # Linear interpolation from c_left at nlo to c_right at nhi-1
            t = (ns - float(nlo)) / float(bw - 1)
            vals = c_left[k] * (1.0 - t) + c_right[k] * t
        idx = (ns - 1).astype(int)
        valid = (idx >= 0) & (idx < N_max)
        c_n[idx[valid]] = np.maximum(vals[valid], 0.0)

    return c_n


class BinMomentRateEquations:
    """
    Chapter 9 bin-moment ODE system for both SIA and vacancy clusters.

    Both populations are reduced to logarithmic size bins.  The intra-bin
    shape function is controlled by ``shape_function``:

      "constant"  — piecewise-constant (1 moment/bin: μ₀)
      "linear"    — hat-function / dual-basis (2 moments/bin: μ₀, μ₁)
      "lognormal" — log-normal (3 moments/bin: μ₀, μ₁, μ₂)

    State vector layout (example for Case 2, QSS He, P moments/bin):
      y[0..i_d-1]                          : discrete SIA
      y[i_d..i_d+P*I_bin-1]               : SIA bin moments
      y[i_d+P*I_bin..+v_d-1]              : discrete VAC
      y[..+v_d..+v_d+P*V_bin-1]           : VAC bin moments
      y[..]                                : He state (Q_tot or Q_k)

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
        # Backward-compat aliases
        self.N = I
        self.M = V

        # ── Shape function ────────────────────────────────────────────────
        sf_raw = str(input_data.reactions.get('shape_function', 'linear')).lower()
        if sf_raw not in _SHAPE_FUNCTIONS:
            raise ValueError(
                f"Unknown shape_function='{sf_raw}'. "
                f"Must be one of {list(_SHAPE_FUNCTIONS.keys())}.")
        self.shape_function = sf_raw
        self.n_mom = _SHAPE_FUNCTIONS[sf_raw]  # moments per bin: 1, 2, or 3

        # Discrete and binned sizes
        self.i_discrete = input_data.i_discrete
        self.v_discrete = input_data.v_discrete

        # Build SIA bin partition (covers sizes i_discrete+1..I)
        I_bin_target = input_data.I_bin
        if I_bin_target > 0 and self.i_discrete < I:
            # Compute r from target bin count
            r_i = (float(I) / max(self.i_discrete, 1)) ** (1.0 / I_bin_target)
            r_i = max(r_i, 1.01)  # ensure r > 1
            self.bins, self.edges = build_bins(
                I, n1=self.i_discrete + 1, r=r_i)
            self.I_bin = len(self.bins)
        else:
            self.bins, self.edges = [], np.array([I + 1])
            self.I_bin = 0
            r_i = 2.0
        self.r = r_i
        # Backward-compat aliases
        self.K    = self.I_bin
        self.K_i  = self.I_bin
        self.n1   = self.i_discrete

        # Build vacancy bin partition (covers sizes v_discrete+1..V)
        V_bin_target = input_data.V_bin
        if V_bin_target > 0 and self.v_discrete < V:
            r_v = (float(V) / max(self.v_discrete, 1)) ** (1.0 / V_bin_target)
            r_v = max(r_v, 1.01)
            self.vac_bins, self.vac_edges = build_bins(
                V, n1=self.v_discrete + 1, r=r_v)
            self.V_bin = len(self.vac_bins)
        else:
            self.vac_bins, self.vac_edges = [], np.array([V + 1])
            self.V_bin = 0
        # Backward-compat alias
        self.K_v = self.V_bin

        # Precompute vacancy bin midpoints for zeroth-moment-only closure
        self.vac_mid = np.array([(mlo + mhi - 1) / 2.0
                                  for mlo, mhi in self.vac_bins])

        # Free He mode: 'dynamic' or 'quasi_steady_state' (mirrors RateEquations)
        raw = str(input_data.reactions.get('he_kinetics', 'dynamic')).lower()
        self.he_kinetics = raw
        self.qss_He      = (raw == 'quasi_steady_state')

        # Pre-compute beta_He for QSS (same as RateEquations)
        kBT_val  = input_data.derived['kBT']
        E_b_hV1  = float(input_data.energetics.get('E_b_hV_1', 2.30))
        E_m_h    = float(input_data.derived['E_m_h'])
        nu_h     = float(input_data.derived['nu_h'])
        self.beta_He = nu_h * np.exp(-(E_b_hV1 + E_m_h) / kBT_val)

        P = self.n_mom  # shorthand for moments per bin

        # ── State-vector index map ────────────────────────────────────────
        # Layout: [discrete SIA | binned SIA moments | discrete VAC |
        #          binned VAC moments | He state]
        po = input_data.physics_option
        n_sia = self.i_discrete + P * self.I_bin  # SIA ODEs
        n_vac = self.v_discrete + P * self.V_bin    # VAC ODEs

        self.i_MOM  = self.i_discrete   # first binned SIA moment index
        self.i_VAC  = n_sia             # first vacancy index

        if 'fusion' in po:
            self.he_mode = 'case1'
            self.i_Q    = n_sia + n_vac           # Q_k per vac bin
            if self.qss_He:
                n_he = self.V_bin
                self.i_He = None
            else:
                n_he = self.V_bin + 1
                self.i_He = n_sia + n_vac + self.V_bin
        else:
            self.he_mode = 'case2'
            self.i_Qtot = n_sia + n_vac
            if self.qss_He:
                n_he = 1
                self.i_He = None
            else:
                n_he = 2
                self.i_He = n_sia + n_vac + 1

        # Conservation accounting ODEs (5 cumulative integrals)
        n_phys = n_sia + n_vac + n_he
        self.i_J_SIA_fixed  = n_phys      # cumulative SIA content to fixed sinks
        self.i_J_SIA_mutual = n_phys + 1   # cumulative SIA to recombination + cavity
        self.i_J_VAC_fixed  = n_phys + 2   # cumulative VAC content to fixed sinks
        self.i_J_VAC_mutual = n_phys + 3   # cumulative VAC to recombination + cavity
        self.i_J_He_sink    = n_phys + 4   # cumulative He to sinks
        self.N_eq = n_phys + 5

        # Cascade production
        # G_He_r sourced from the Excel Reactions sheet (single authoritative
        # source — see defect_production.production_rates docstring).
        spectrum = input_data.derived['spectrum']
        G        = input_data.derived['G']
        G_He_r   = input_data.derived['G_He_r']
        Pr_SIA, Pr_VAC, G_He = production_rates(G, spectrum, I, V, G_He_r)
        self.Pr_SIA = Pr_SIA[1:]   # [I], index k → size k+1
        self.Pr_VAC = Pr_VAC[1:]   # [V]
        self.G_He   = G_He

        print(f"BinMomentRateEquations: i_discrete={self.i_discrete}"
              f"  I_bin={self.I_bin}  V_bin={self.V_bin}  N_eq={self.N_eq}"
              f"  he_mode='{self.he_mode}'  r_i={r_i:.3f}"
              f"  shape_function='{self.shape_function}' ({P} mom/bin)")
        if self.I_bin > 0:
            print(f"  SIA bin edges (first 5): {list(self.edges[:6])}")
        if self.V_bin > 0:
            print(f"  VAC bin edges (first 5): {list(self.vac_edges[:6])}")

    def compute_c_h_qss(self, c_v, Q_tot=None, Q_m=None):
        """QSS free He — same formula as RateEquations.compute_c_h_qss."""
        rr   = self.rr
        sink = np.sum(rr.K_HeV * c_v) + rr.k2_He_scalar
        if self.he_mode == 'case2':
            source = self.G_He + self.beta_He * (Q_tot if Q_tot is not None else 0.0)
        else:
            He_emit = self.beta_He * np.sum(Q_m) if Q_m is not None else 0.0
            source  = self.G_He + He_emit
        return source / max(sink, 1e-300)

    def get_initial_conditions(self):
        """Return y0 — near-zero but consistent moments + concentrations.

        Bin moments are initialized so that μ₁/μ₀ = midpoint and
        μ₂/μ₀ = midpoint², ensuring the reconstruction places the
        distribution inside the bin range.  Setting all moments to
        C_floor would give μ₁/μ₀ = 1 (outside bin range), causing
        garbage reconstruction and immediate CVODE failure.

        The C_floor fallback default (1e-15) must match the value the
        C++ solver clamps to, so the initial state is consistent with
        the solver's floor.
        """
        C_floor = float(self.inp.reactions.get('C_floor', 1e-15))
        y0 = np.full(self.N_eq, C_floor)
        P = self.n_mom

        # SIA bin moments: μ₀ = C_floor, μ₁ = C_floor × midpoint, ...
        for k, (nlo, nhi) in enumerate(self.bins):
            mid = 0.5 * (nlo + nhi - 1)
            idx = self.i_discrete + P * k
            y0[idx] = C_floor                           # μ₀
            if P >= 2:
                y0[idx + 1] = C_floor * mid             # μ₁
            if P >= 3:
                y0[idx + 2] = C_floor * mid * mid       # μ₂

        # Vacancy bin moments: same consistent initialization
        iv = self.i_VAC
        for k, (mlo, mhi) in enumerate(self.vac_bins):
            mid = 0.5 * (mlo + mhi - 1)
            idx = iv + self.v_discrete + P * k
            y0[idx] = C_floor                           # μ₀
            if P >= 2:
                y0[idx + 1] = C_floor * mid             # μ₁
            if P >= 3:
                y0[idx + 2] = C_floor * mid * mid       # μ₂

        return y0

    def _unpack_moments(self, y_block, K_bin):
        """Extract mu0, mu1, mu2 from a contiguous block of P*K_bin values."""
        P = self.n_mom
        mu0 = np.maximum(y_block[0::P][:K_bin], 0.0)
        mu1 = y_block[1::P][:K_bin] if P >= 2 else None
        mu2 = y_block[2::P][:K_bin] if P >= 3 else None
        return mu0, mu1, mu2

    def _project_to_moments(self, dc, bins, N_max, dydt, offset, K_bin):
        """Project per-size rates dc onto P moments per bin in dydt."""
        P = self.n_mom
        for k, (nlo, nhi) in enumerate(bins):
            ns  = np.arange(nlo, nhi)
            idx = ns - 1
            valid = (idx >= 0) & (idx < N_max)
            dc_valid = dc[idx[valid]]
            n_valid  = ns[valid].astype(float)
            dydt[offset + P * k] = np.sum(dc_valid)                     # dμ₀/dt
            if P >= 2:
                dydt[offset + P * k + 1] = np.sum(n_valid * dc_valid)   # dμ₁/dt
            if P >= 3:
                dydt[offset + P * k + 2] = np.sum(n_valid * n_valid * dc_valid)  # dμ₂/dt

    def ode_system(self, t, y):
        """ODE RHS using hybrid discrete + bin-moment layout.

        REFERENCE ONLY — not used by the solver; the production RHS is the
        C++ implementation.  May be stale.  Scheduled for replacement by the
        Stage-3 graph-walker.
        """
        I_bin = self.I_bin
        V_bin = self.V_bin
        V   = self.V
        I   = self.I
        i_d = self.i_discrete
        v_d = self.v_discrete
        rr  = self.rr
        inp = self.inp
        P   = self.n_mom
        sf  = self.shape_function
        dydt = np.zeros(self.N_eq)

        # ── Unpack SIA: discrete + binned → full c_n[0..I-1] ─────────────
        c_n = np.zeros(I)
        # Discrete sizes 1..i_discrete
        c_n[:i_d] = np.maximum(y[:i_d], 0.0)
        # Binned sizes i_discrete+1..I (from moments)
        if I_bin > 0:
            mom_start = i_d
            sia_mu0, sia_mu1, sia_mu2 = self._unpack_moments(
                y[mom_start:mom_start + P * I_bin], I_bin)
            c_binned = reconstruct_distribution(
                sf, sia_mu0, sia_mu1, sia_mu2, self.bins, I)
            c_n[i_d:] = c_binned[i_d:]

        # ── Unpack VAC: discrete + binned → full c_v[0..V-1] ─────────────
        iv = self.i_VAC
        c_v = np.zeros(V)
        # Discrete sizes 1..v_discrete
        c_v[:v_d] = np.maximum(y[iv:iv + v_d], 0.0)
        # Binned sizes v_discrete+1..V (from tracked moments)
        if V_bin > 0:
            vac_start = iv + v_d
            vac_mu0, vac_mu1, vac_mu2 = self._unpack_moments(
                y[vac_start:vac_start + P * V_bin], V_bin)
            c_v_binned = reconstruct_distribution(
                sf, vac_mu0, vac_mu1, vac_mu2, self.vac_bins, V)
            c_v[v_d:] = c_v_binned[v_d:]

        # ── He state ──────────────────────────────────────────────────────
        if self.he_mode == 'case2':
            Q_tot = max(y[self.i_Qtot], 0.0)
            Q_k   = None
            if self.qss_He:
                c_h = self.compute_c_h_qss(c_v, Q_tot=Q_tot)
            else:
                c_h = max(y[self.i_He], 0.0)
        else:
            Q_k   = np.maximum(y[self.i_Q:self.i_Q + V_bin], 0.0)
            Q_tot = np.sum(Q_k)
            if self.qss_He:
                c_h = self.compute_c_h_qss(c_v, Q_m=Q_k)
            else:
                c_h = max(y[self.i_He], 0.0)

        ci1 = c_n[0] if len(c_n) > 0 else 0.0
        cv1 = c_v[0]
        i_mobile = inp.derived['i_mobile']
        reflect  = inp.derived.get('boundary_flux', 'absorption') == 'reflection'
        m13     = np.arange(1.0, V + 1) ** (1.0 / 3.0)
        denom_m = 1.0 + rr.B_rot * rr.L_hat**2 / m13

        # ══════════════════════════════════════════════════════════════════
        # SIA per-size dc_n/dt (full per-size, then project)
        # ══════════════════════════════════════════════════════════════════
        dc_n = np.zeros(I)

        dc_n += self.Pr_SIA
        dc_n[:-1] += rr.G_SIA[1:] * c_n[1:]
        dc_n -= rr.G_SIA * c_n
        dc_n[1:] += rr.K_SIA_grow[:-1] * ci1 * c_n[:-1]  # gain at n+1
        dc_n -= rr.K_SIA_grow * ci1 * c_n                  # target loss
        if reflect:  # suppress I_1 + I_I → I_{I+1} (reaction blocked at wall)
            dc_n[-1] += rr.K_SIA_grow[-1] * ci1 * c_n[-1]  # undo target loss
            dc_n[0]  += ci1 * rr.K_SIA_grow[-1] * c_n[-1]  # undo monomer depletion
        # Monomer projectile depletion
        dc_n[0] -= ci1 * np.dot(rr.K_SIA_grow, c_n)

        # V–I annihilation
        dc_n[0] -= rr.K_iv * cv1 * ci1
        dc_n[1:] -= rr.K_SIA_shrink[1:] * cv1 * c_n[1:]
        dc_n[:-1] += rr.K_SIA_shrink[1:] * cv1 * c_n[1:]

        # SIA cluster–cavity absorption
        sum_K3D_cv_m2 = rr.K_3D_cav_pref * np.dot(m13[1:], c_v[1:])
        dc_n[0] -= c_n[0] * sum_K3D_cv_m2
        sum_K3D_cv = rr.K_3D_cav_pref * np.dot(m13, c_v)
        for n in range(2, min(4, i_mobile + 1)):
            dc_n[n - 1] -= c_n[n - 1] * sum_K3D_cv
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            dc_n[n - 1] -= c_n[n - 1] * k_pref * np.dot(m13 / denom_m, c_v)

        dc_n -= rr.k2_SIA * c_n

        # ── Project dc_n into dydt ───────────────────────────────────────
        # Discrete sizes: direct copy
        dydt[:i_d] = dc_n[:i_d]
        # Binned sizes: project onto P moments per bin
        self._project_to_moments(dc_n, self.bins, I, dydt, i_d, I_bin)

        # ══════════════════════════════════════════════════════════════════
        # Vacancy per-size dc_v/dt → project onto vacancy bin moments
        # ══════════════════════════════════════════════════════════════════
        C_vac_tot = np.sum(c_v)
        ell_bar   = Q_tot / max(C_vac_tot, 1e-200)
        G_VAC_eff = rr.G_VAC.copy()
        if ell_bar > 0.01:
            for m in range(1, min(V, 20) + 1):
                G_VAC_eff[m - 1] = rr.alpha_bubble_fn(m, ell_bar * m**(2.0/3.0))

        dc_v = np.zeros(V)

        # Production
        dc_v += self.Pr_VAC

        # Thermal vacancy emission (gain + loss)
        # V_m → V_{m-1} + V_1: loss at m, gain at m-1 (residual), gain at 1 (emitted monomer)
        dc_v -= G_VAC_eff * c_v
        dc_v[:-1] += G_VAC_eff[1:] * c_v[1:]
        dc_v[0] += np.sum(G_VAC_eff[1:] * c_v[1:])   # emitted monomers

        # V–V monomer growth (gain + loss)
        dc_v -= rr.K_VAC_grow * cv1 * c_v             # target loss
        dc_v[1:] += rr.K_VAC_grow[:-1] * cv1 * c_v[:-1]  # gain at m+1
        if reflect:  # suppress V_1 + V_V → V_{V+1} (reaction blocked at wall)
            dc_v[-1] += rr.K_VAC_grow[-1] * cv1 * c_v[-1]  # undo target loss
            dc_v[0]  += cv1 * rr.K_VAC_grow[-1] * c_v[-1]  # undo monomer depletion
        # Vacancy monomer projectile depletion
        dc_v[0] -= cv1 * np.dot(rr.K_VAC_grow, c_v)

        # SIA-induced cavity shrinkage
        dc_v[0] -= rr.K_iv * ci1 * cv1
        if V >= 2:
            dc_v[0] += rr.K_VAC_shrink[1] * ci1 * c_v[1]
        dc_v[1:] -= rr.K_VAC_shrink[1:] * ci1 * c_v[1:]
        if V >= 3:
            dc_v[1:-1] += rr.K_VAC_shrink[2:] * ci1 * c_v[2:]
        for n in range(2, min(4, i_mobile + 1)):
            cn = c_n[n - 1]
            dc_v -= rr.K_VAC_shrink * cn * c_v
            if n < V:
                dc_v[:V - n] += rr.K_VAC_shrink[n:] * cn * c_v[n:]
        m_arr = np.arange(1.0, V + 1)
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n - 1]
            if k_pref < 1e-300:
                continue
            cn = c_n[n - 1]
            k_loss = k_pref * m13 / denom_m
            dc_v -= k_loss * cn * c_v
            if n < V:
                mp    = m_arr[:V - n] + float(n)
                mp13  = mp ** (1.0 / 3.0)
                denom = 1.0 + rr.B_rot * rr.L_hat**2 / mp13
                dc_v[:V - n] += (k_pref * mp13 / denom) * cn * c_v[n:]

        # Vacancy monomer consumed by SIA loop shrinkage:
        # I_n + V_1 → I_{n-1} for n≥2 (n=1 already in K_iv)
        dc_v[0] -= cv1 * np.dot(rr.K_SIA_shrink[1:], c_n[1:])

        # Fixed sinks — only mobile vacancy clusters (v ≤ v_mobile)
        v_mobile = self.inp.derived['v_mobile']
        mask_mobile_v = np.zeros(V)
        mask_mobile_v[:min(v_mobile, V)] = 1.0
        dc_v -= rr.k2_vac_scalar * c_v * mask_mobile_v

        # ── Project dc_v into dydt ───────────────────────────────────────
        # Discrete sizes: direct copy
        dydt[iv:iv + v_d] = dc_v[:v_d]
        # Binned sizes: project onto P moments per bin
        vac_mom_start = iv + v_d
        self._project_to_moments(dc_v, self.vac_bins, V, dydt,
                                 vac_mom_start, V_bin)

        # ══════════════════════════════════════════════════════════════════
        # He equations
        # ══════════════════════════════════════════════════════════════════
        if self.he_mode == 'case2':
            He_uptake  = np.sum(rr.K_HeV * c_h * c_v)
            He_release = rr.alpha_He * Q_tot if hasattr(rr, 'alpha_He') else 0.0
            # He lost only from mobile voids reaching fixed sinks
            ms_arr = np.arange(1.0, V + 1.0)
            ell_m_arr = ell_bar * ms_arr ** (2.0 / 3.0)
            He_sink = rr.k2_vac_scalar * np.sum(ell_m_arr[:v_mobile] * c_v[:v_mobile])
            dydt[self.i_Qtot] = He_uptake - He_release - He_sink
            if not self.qss_He:
                dydt[self.i_He] = (self.G_He
                                   - He_uptake
                                   - rr.k2_He_scalar * c_h
                                   + He_release)
        else:
            # Case 1: Q_k per vacancy bin (one per bin)
            dQ = dydt[self.i_Q:self.i_Q + V_bin]
            He_cap_total  = 0.0
            He_emit_total = 0.0
            for kv, (mlo, mhi) in enumerate(self.vac_bins):
                mu0_kv = max(vac_mu0[kv], 1e-200)
                ell_kv = Q_k[kv] / mu0_kv
                for m in range(mlo, mhi):
                    mi = m - 1
                    if mi < 0 or mi >= V:
                        continue
                    cm = c_v[mi]
                    he_cap  = rr.K_HeV[mi] * c_h * cm
                    alpha_h = rr.alpha_He_emit_fn(m, max(int(round(ell_kv)), 0))
                    q_m_approx = ell_kv * cm
                    he_emit = alpha_h * q_m_approx
                    dQ[kv] += he_cap - he_emit
                    if m <= v_mobile:   # only mobile voids
                        dQ[kv] -= rr.k2_vac_scalar * q_m_approx
                    He_cap_total  += he_cap
                    He_emit_total += he_emit
            if not self.qss_He:
                dydt[self.i_He] = (self.G_He
                                   - He_cap_total
                                   - rr.k2_He_scalar * c_h
                                   + He_emit_total)

        # ══════════════════════════════════════════════════════════════════
        # Conservation accounting ODEs (cumulative integrals, exact via ODE)
        # ══════════════════════════════════════════════════════════════════

        # (1) SIA content to fixed sinks: Σ_n n · k2_SIA[n] · c_n
        ns_all = np.arange(1.0, I + 1)
        dydt[self.i_J_SIA_fixed] = np.dot(ns_all, rr.k2_SIA * c_n)

        # (2) SIA content to mutual annihilation (recombination + cavity):
        # Channel (a): mobile SIA clusters absorbed by cavities
        # Channel (b): mobile vacancies shrinking SIA clusters
        mutual = rr.K_iv * ci1 * cv1                            # V_1+I_1 recomb (1 SIA)
        sum_K3D_cv_m2 = rr.K_3D_cav_pref * np.dot(m13[1:], c_v[1:])
        mutual += ci1 * sum_K3D_cv_m2                           # I_1 → cavities m≥2 (1 SIA)
        mutual += cv1 * np.sum(rr.K_SIA_shrink[1:] * c_n[1:])   # V_1+I_n→I_{n-1} (1 SIA each)
        for n in range(2, min(4, i_mobile + 1)):
            mutual += n * c_n[n-1] * rr.K_3D_cav_pref * np.dot(m13, c_v)
        for n in range(4, min(I, i_mobile) + 1):
            k_pref = rr.K_1D_pref[n-1]
            if k_pref < 1e-300:
                continue
            mutual += n * c_n[n-1] * k_pref * np.dot(m13 / denom_m, c_v)
        dydt[self.i_J_SIA_mutual] = mutual

        # (3) VAC content to fixed sinks: Σ_{m≤v_mobile} m · k2_vac · c_m
        vm = min(v_mobile, V)
        ms_v = np.arange(1.0, vm + 1)
        dydt[self.i_J_VAC_fixed] = rr.k2_vac_scalar * np.dot(ms_v, c_v[:vm])

        # (3b) VAC content to mutual annihilation:
        # When I_n + V_{m'} react, vacancies destroyed = min(m', n) = SIA destroyed.
        # Therefore J_VAC_mutual = J_SIA_mutual (matches C++ implementation).
        dydt[self.i_J_VAC_mutual] = mutual  # J_VAC_mutual = J_SIA_mutual (physically correct)

        # (4) He to sinks: k2_He · c_h + k2_vac · Σ_{m≤v_mobile} ℓ_m · c_m
        sum_cv = np.sum(c_v)
        he_sink = rr.k2_He_scalar * c_h
        if sum_cv > 1e-300 and Q_tot > 0:
            ell_bar_loc = Q_tot / sum_cv
            ell_m_v = ell_bar_loc * np.arange(1.0, vm + 1) ** (2.0/3.0)
            he_sink += rr.k2_vac_scalar * np.dot(ell_m_v, c_v[:vm])
        dydt[self.i_J_He_sink] = he_sink

        return dydt

    def _extract_mu1_content(self, y, offset, K_bin, bins):
        """Extract first-moment content from state vector.

        For P >= 2, μ₁ is tracked explicitly (stride P, offset 1).
        For P == 1 (constant closure), μ₁ is approximated as μ₀ × midpoint.
        """
        P = self.n_mom
        if P >= 2:
            return np.sum(y[offset + 1::P][:K_bin])
        else:
            # Approximate: μ₁ ≈ μ₀ × arithmetic midpoint
            content = 0.0
            for k, (nlo, nhi) in enumerate(bins):
                mu0_k = max(y[offset + k], 0.0)
                content += mu0_k * (nlo + nhi - 1) / 2.0
            return content

    def conservation_diagnostic(self, t, y, G_tot):
        """
        FP conservation diagnostic for hybrid discrete + bin-moment system.

        δ_FP^bin = |SIA_content − VAC_content| / (G·t)
        """
        i_d = self.i_discrete
        v_d = self.v_discrete
        iv  = self.i_VAC

        # SIA content: discrete + binned first moments
        ns_disc = np.arange(1, i_d + 1, dtype=float)
        SIA_content = np.dot(ns_disc, np.maximum(y[:i_d], 0.0))
        if self.I_bin > 0:
            SIA_content += self._extract_mu1_content(
                y, i_d, self.I_bin, self.bins)

        # VAC content: discrete + binned first moments
        ms_disc = np.arange(1, v_d + 1, dtype=float)
        VAC_content = np.dot(ms_disc, np.maximum(y[iv:iv + v_d], 0.0))
        if self.V_bin > 0:
            VAC_content += self._extract_mu1_content(
                y, iv + v_d, self.V_bin, self.vac_bins)

        delta_J = abs(SIA_content - VAC_content)
        denom   = max(G_tot * max(t, 1e-20), 1e-300)
        return delta_J / denom
