"""
core/reductions/bin_moment.py — Logarithmic size-bin moment reduction.

Host-independent.  Implements the "reconstruct -> walk -> project"
reduction of Ghoniem (2026), Section 5.4: a per-size cluster ladder of
length N is replaced by a discrete prefix plus O(log N) logarithmically
spaced bins, each carrying P in {1,2,3} moments

    mu_k^(0) = sum_{n in B_k} c_n
    mu_k^(1) = sum_{n in B_k} n  c_n
    mu_k^(2) = sum_{n in B_k} n^2 c_n

evolved by the three-step procedure executed at every RHS call:

  R  reconstruct a per-size distribution c_n from the tracked moments
     under the chosen intra-bin closure (piecewise-constant / linear
     dual-basis / log-normal);
  E  evaluate the *exact* per-size master-equation rates with a
     :class:`GraphWalker` — no bin-averaged rate constants;
  P  project the per-size rates onto the tracked moments,
     d mu_k^(p)/dt = sum_{n in B_k} n^p (dc_n/dt).

Because step E uses the exact per-size rates, every conservation
property of the master equation is inherited automatically — no
inter-bin flux has to be made upwind-consistent (Section 5.4).

The pure bin-partition and closure mathematics is reused from the
reviewed, validated ``bin_moment_rates`` module (those functions are
host-independent); this class composes them with a graph walk.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np

from ...bin_moment_rates import (  # reviewed, host-independent bin mathematics
    build_bins,
    moments_from_distribution,
    n_moments_for_shape,
    reconstruct_distribution,
)


class BinMomentReduction:
    """A logarithmic bin-moment reduction of one cluster size axis.

    Parameters
    ----------
    n_max : int
        Largest cluster size on the full per-size ladder.
    n_discrete : int
        Sizes 1..n_discrete are kept discrete (one ODE per size); sizes
        above are aggregated into geometric bins.
    n_bins : int
        Target number of geometric bins above the discrete prefix.
    shape_function : {'constant', 'linear', 'lognormal'}
        The intra-bin closure; selects P = 1, 2 or 3 moments per bin.
    """

    def __init__(self, n_max: int, n_discrete: int, n_bins: int,
                 shape_function: str = "linear"):
        self.n_max = int(n_max)
        self.n_discrete = int(n_discrete)
        self.shape_function = str(shape_function).lower()
        self.moments_per_bin = n_moments_for_shape(self.shape_function)

        if n_bins > 0 and self.n_discrete < self.n_max:
            r = (float(self.n_max) / max(self.n_discrete, 1)) ** (1.0 / n_bins)
            self.r = max(r, 1.01)
            self.bins, self.edges = build_bins(
                self.n_max, n1=self.n_discrete + 1, r=self.r)
        else:
            self.bins, self.edges = [], np.array([self.n_max + 1])
            self.r = 2.0
        self.n_bins = len(self.bins)

    # ── layout ────────────────────────────────────────────────────────────
    @property
    def length(self) -> int:
        """Number of ODE entries: discrete prefix + P moments per bin."""
        return self.n_discrete + self.moments_per_bin * self.n_bins

    # ── step R: moments -> per-size distribution ──────────────────────────
    def reconstruct(self, reduced: np.ndarray) -> np.ndarray:
        """Reconstruct the full per-size distribution c_n (length n_max)."""
        P = self.moments_per_bin
        c = np.zeros(self.n_max)
        c[:self.n_discrete] = np.maximum(reduced[:self.n_discrete], 0.0)
        if self.n_bins:
            blk = reduced[self.n_discrete:self.n_discrete + P * self.n_bins]
            mu0 = np.maximum(blk[0::P][:self.n_bins], 0.0)
            mu1 = blk[1::P][:self.n_bins] if P >= 2 else None
            mu2 = blk[2::P][:self.n_bins] if P >= 3 else None
            c_bin = reconstruct_distribution(
                self.shape_function, mu0, mu1, mu2, self.bins, self.n_max)
            c[self.n_discrete:] = c_bin[self.n_discrete:]
        return c

    # ── step P: per-size rates -> moment rates ────────────────────────────
    def project(self, dc: np.ndarray) -> np.ndarray:
        """Project per-size rates dc_n/dt onto the tracked moment rates."""
        P = self.moments_per_bin
        out = np.zeros(self.length)
        out[:self.n_discrete] = dc[:self.n_discrete]
        for k, (nlo, nhi) in enumerate(self.bins):
            ns = np.arange(nlo, nhi)
            idx = ns - 1
            valid = (idx >= 0) & (idx < self.n_max)
            dcv = dc[idx[valid]]
            nv = ns[valid].astype(float)
            base = self.n_discrete + P * k
            out[base] = dcv.sum()                       # d mu0/dt
            if P >= 2:
                out[base + 1] = float(nv @ dcv)          # d mu1/dt
            if P >= 3:
                out[base + 2] = float((nv * nv) @ dcv)   # d mu2/dt
        return out

    # ── moments of an initial per-size distribution ───────────────────────
    def moments_of(self, c: np.ndarray) -> np.ndarray:
        """Pack a per-size distribution into the reduced moment vector."""
        P = self.moments_per_bin
        reduced = np.zeros(self.length)
        reduced[:self.n_discrete] = c[:self.n_discrete]
        if self.n_bins:
            mu0, mu1, mu2 = moments_from_distribution(c, self.bins, n_mom=P)
            for k in range(self.n_bins):
                base = self.n_discrete + P * k
                reduced[base] = mu0[k]
                if P >= 2:
                    reduced[base + 1] = mu1[k]
                if P >= 3:
                    reduced[base + 2] = mu2[k]
        return reduced

    def __repr__(self) -> str:
        return (f"<BinMomentReduction n_max={self.n_max} "
                f"n_discrete={self.n_discrete} n_bins={self.n_bins} "
                f"closure={self.shape_function!r} P={self.moments_per_bin} "
                f"len={self.length}>")
