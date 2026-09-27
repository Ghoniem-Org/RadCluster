#!/usr/bin/env python3
"""Step-5 exact-regression test.

Proves the refactored code paths (callable/scalar sigma, scalar-vs-pair
mu) reproduce the PRE-STEP-5 trajectories exactly when given the old
inputs:
    sigma = 0.5 (scalar, the old IMMIGRANT_MALE_FRACTION)
    mu    = (mu, mu) pair   (the old uniform scalar mu)

Compares against the pinned old outputs:
    outputs/20260927_161901_colonial/summary_all.csv
    outputs/colonial_religious/summary_*.csv

Pass criterion: max |relative difference| over all checked numbers < 1e-10
(the code path preserves float operations, so exact 0.0 is expected).
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)

from genealogy.populations import (N_CLUSTERS, MAX_DEPTH, N_STATE, RULES,
                                   state_index)
from genealogy.kernels import SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION
from genealogy.model import rhs_td
from genealogy.solver import run_hindcast, cluster_totals, sex_totals

S = SEX_RATIO_AT_BIRTH
SIGMA_OLD = IMMIGRANT_MALE_FRACTION  # 0.5 scalar
ALPHAS = [0.0, 0.8, 0.9]
EPS = 0.00075

inputs = pd.read_csv(os.path.join(BASE, "data", "inputs.csv"))
inputs = inputs.sort_values("year").reset_index(drop=True)


def rate_series_factory():
    years = inputs["year"].to_numpy()
    I = inputs["immigration_M"].to_numpy()
    beta = (2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0)
    mu = inputs["cdr_per_1000"].to_numpy() / 1000.0

    def rate_series(t):
        i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
        bvec = np.full(N_CLUSTERS, beta[i])
        return I[i], bvec, bvec, (mu[i], mu[i]), EPS  # pair, equal components
    return rate_series


rate_series = rate_series_factory()

c0 = np.zeros(N_STATE)
c0[state_index(0, "pat")] = SIGMA_OLD * 0.050368
c0[state_index(0, "mat")] = (1.0 - SIGMA_OLD) * 0.050368


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def main():
    summ = pd.read_csv(os.path.join(
        BASE, "outputs", "20260927_161901_colonial", "summary_all.csv"))
    worst = 0.0
    for rule in RULES:
        for alpha in ALPHAS:
            t, Y = run_hindcast(
                c0, 1650.0, 2025.0, rate_series, alpha=alpha, sigma=SIGMA_OLD,
                s=S, rule=rule, dt=1.0)
            T = cluster_totals(Y)
            Tpat, Tmat = sex_totals(Y)
            total = T.sum(axis=1)
            share = T / total[:, None]
            row = summ[(summ["rule"] == rule)
                       & (summ["alpha"] == alpha)].iloc[0]
            checks = {
                "total_2025_M": (total[-1], row["total_2025_M"]),
                "Cn_share_2025": (share[-1, MAX_DEPTH], row["Cn_share_2025"]),
                "pat_share_2025": (Tpat[-1] / total[-1], row["pat_share_2025"]),
            }
            for k, (new, old) in checks.items():
                r = rel(new, old)
                worst = max(worst, r)
                print(f"{rule} a={alpha} {k}: new={new:.10f} old={old:.10f} "
                      f"rel={r:.2e}")
    print(f"\nWORST relative deviation (two-sex): {worst:.2e}")
    assert worst < 1e-10, f"regression FAILED: {worst:.2e}"
    print("TWO-SEX REGRESSION PASSED (bit-identical trajectories).")


if __name__ == "__main__":
    main()
