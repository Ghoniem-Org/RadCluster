#!/usr/bin/env python3
"""Step-5 religious-model regression: one config (shallowest, alpha=0.8) with
old inputs (sigma=0.5 scalar, mu scalar) through the refactored
rhs_religious, compared to the pinned old outputs/colonial_religious run.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)

from genealogy.populations import N_CLUSTERS, MAX_DEPTH
from genealogy.religion import (G, N_STATE_R, rstate_index, RHO_VEC,
                                immigrant_composition, disaffiliation_rate)
from genealogy.kernels import SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION
from genealogy.model import rhs_religious

S = SEX_RATIO_AT_BIRTH
SIGMA_OLD = IMMIGRANT_MALE_FRACTION
RULE, ALPHA, DELTA_MAX = "shallowest_inheritance", 0.8, 0.008
EPS = 0.00075

inputs = pd.read_csv(os.path.join(BASE, "data", "inputs.csv"))
inputs = inputs.sort_values("year").reset_index(drop=True)
years = inputs["year"].to_numpy()
I = inputs["immigration_M"].to_numpy()
beta = 2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0
mu = inputs["cdr_per_1000"].to_numpy() / 1000.0


def rate_series(t):
    i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
    bvec = np.full(N_CLUSTERS, beta[i])
    return (I[i], bvec, mu[i], EPS,  # mu scalar, as in the old driver
            immigrant_composition(t), disaffiliation_rate(t, DELTA_MAX))


c0 = np.zeros(N_STATE_R)
c0[rstate_index(0, "pat", "pro")] = SIGMA_OLD * 0.050368
c0[rstate_index(0, "mat", "pro")] = (1.0 - SIGMA_OLD) * 0.050368


def rhs_td(t, c):
    I_, beta_, mu_, eps_, iot_, dlt_ = rate_series(t)
    return rhs_religious(t, c, I=I_, beta=beta_, alpha=ALPHA, mu=mu_,
                         eps=eps_, sigma=SIGMA_OLD, s=S, rule=RULE,
                         rho=RHO_VEC, iota=iot_, delta=dlt_)


def main():
    t_eval = np.arange(1650.0, 2025.0 + 1.0, 1.0)
    sol = solve_ivp(rhs_td, (1650.0, 2025.0), c0, t_eval=t_eval,
                    method="RK45", rtol=1e-8, atol=1e-10)
    assert sol.success
    Y = sol.y.T
    Cg = Y.reshape(len(sol.t), G, 2, N_CLUSTERS)
    total = Y.sum(axis=1)
    cn = Cg[:, :, :, MAX_DEPTH].sum(axis=(1, 2)) / total
    new = {"total_2025_M": total[-1], "Cn_share_2025": cn[-1]}
    df = pd.read_csv(os.path.join(
        BASE, "outputs", "colonial_religious",
        f"results_religious_{RULE}_alpha{ALPHA}.csv"), nrows=376)
    row = df[df["year"] == 2025.0].iloc[0]
    old = {"total_2025_M": row["total_M"], "Cn_share_2025": row["Cn_share"]}
    worst = 0.0
    for k in new:
        r = abs(new[k] - old[k]) / abs(old[k])
        worst = max(worst, r)
        print(f"{k}: new={new[k]:.10f} old={old[k]:.10f} rel={r:.2e}")
    assert worst < 1e-9, f"religious regression FAILED: {worst:.2e}"
    print("RELIGIOUS REGRESSION PASSED.")


if __name__ == "__main__":
    main()
