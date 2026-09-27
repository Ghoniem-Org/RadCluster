#!/usr/bin/env python3
"""Age-uniform consistency: with age-UNIFORM rates the age-structured model
must reproduce the aggregate two-sex model to tight tolerance.

Uniform means: beta_a = const for ALL ages (no 15-49 gating), mu_{s,a} =
const, immigration spread uniformly over bands. Then:
  - aging fluxes telescope to zero on age-aggregation,
  - births aggregate to the aggregate kernel (weights w_q = beta_u * P_q),
  - immigration/death/emigration aggregate trivially,
so sum_a dc_{d,s,a}/dt == (aggregate rhs)_{d,s} EXACTLY. The two ODE
solutions must therefore coincide; we assert max |diff| < 1e-6 * total.

Also checks: the total-balance identity sum(rhs_age) ==
total_balance_age(...), and pat/mat symmetry under symmetric inputs.
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from genealogy import solver
from genealogy.age import (N_AGE, N_CLUSTERS, N_STATE_AGE, rhs_age,
                           total_balance_age, birth_inflow_age,
                           IMMIGRANT_MALE_FRACTION as SIGMA,
                           astate_index)
from genealogy.populations import N_CLUSTERS as _NC  # noqa: F401  (re-export check)
from genealogy.kernels import SEX_RATIO_AT_BIRTH as S
from genealogy.model import N_STATE as _NS  # noqa: F401


def aggregate_c0_all_in_c0(total, sigma):
    c0 = np.zeros(2 * N_CLUSTERS)
    c0[0] = sigma * total
    c0[N_CLUSTERS] = (1.0 - sigma) * total
    return c0

BETA_U = 0.030     # uniform per-capita birth rate (all ages)
MU_U = 0.012       # uniform per-capita death rate
EPS = 0.00075
I_FLOW = 1.0       # M/yr
ALPHA = 0.8
T_END = 200.0
P0 = 3.0           # initial total (M), all in C_0


def age_uniform_inputs():
    beta_a = np.full(N_AGE, BETA_U)
    mu_a = np.full(N_AGE, MU_U)
    imm = np.full(N_AGE, 1.0 / N_AGE)
    return beta_a, mu_a, mu_a, imm


def run_pair(rule):
    beta_a, mu_pat, mu_mat, imm = age_uniform_inputs()
    # aggregate IC and run
    c0 = aggregate_c0_all_in_c0(P0, SIGMA)
    _, Y_agg = solver.run_scenario(
        c0, t_end=T_END, dt=1.0, I=I_FLOW,
        beta_pat=np.full(N_CLUSTERS, BETA_U),
        beta_mat=np.full(N_CLUSTERS, BETA_U),
        mu=MU_U, eps=EPS, alpha=ALPHA, sigma=SIGMA, s=S, rule=rule)
    # age IC: same (d,s) totals, spread uniformly over bands
    c0a = np.zeros(N_STATE_AGE)
    for d in range(N_CLUSTERS):
        for si, sfrac in ((0, SIGMA), (1, 1.0 - SIGMA)):
            tot = c0[si * N_CLUSTERS + d]
            base = astate_index(d, si, 0)
            c0a[base:base + N_AGE] = tot / N_AGE
    _, Y_age = solver.run_age_scenario(
        c0a, t_end=T_END, dt=1.0, I=I_FLOW, beta_a=beta_a,
        mu_pat_a=mu_pat, mu_mat_a=mu_mat, eps=EPS, alpha=ALPHA,
        sigma=SIGMA, s=S, rule=rule, imm_age=imm)
    return Y_agg, Y_age


def main():
    ok = True
    for rule in ("shallowest_inheritance", "deepest_inheritance"):
        Y_agg, Y_age = run_pair(rule)
        T_age = solver.cluster_totals_age(Y_age)   # (steps, n)
        T_agg = Y_agg[:, :N_CLUSTERS] + Y_agg[:, N_CLUSTERS:]
        diff = np.abs(T_age - T_agg)
        scale = T_agg.sum(axis=1, keepdims=True)
        rel = (diff / np.maximum(scale, 1e-12)).max()
        absd = diff.max()
        print(f"{rule}: max abs diff = {absd:.3e} M, "
              f"max rel diff = {rel:.3e}")
        if rel > 1e-6:
            print("  FAIL: exceeds 1e-6 relative tolerance")
            ok = False
        # sex totals too
        sp_age, sm_age = solver.sex_totals_age(Y_age)
        sp_agg, sm_agg = solver.sex_totals(Y_agg)
        r2 = max(np.abs(sp_age - sp_agg).max(), np.abs(sm_age - sm_agg).max()) \
            / sp_agg.max()
        print(f"  sex totals max rel diff = {r2:.3e}")
        if r2 > 1e-6:
            print("  FAIL: sex totals")
            ok = False

    # balance identity: sum(rhs_age) == total_balance_age, random states
    rng = np.random.default_rng(0)
    beta_a, mu_pat, mu_mat, imm = age_uniform_inputs()
    for trial in range(5):
        c = rng.uniform(0, 5, N_STATE_AGE)
        dc, B = rhs_age(0.0, c, I=I_FLOW, beta_a=beta_a, mu_pat_a=mu_pat,
                        mu_mat_a=mu_mat, eps=EPS, alpha=ALPHA, sigma=SIGMA,
                        s=S, rule="shallowest_inheritance", imm_age=imm)
        C = c.reshape(2, N_CLUSTERS, N_AGE)
        bal = total_balance_age(I_FLOW, B, C, mu_pat, mu_mat, EPS)
        err = abs(dc.sum() - bal) / max(abs(bal), 1e-12)
        assert err < 1e-10, f"balance identity violated: {err}"
    print("balance identity: OK (5 random states, rel err < 1e-10)")

    # pat/mat symmetry under symmetric inputs and IC
    beta_a, mu_pat, mu_mat, imm = age_uniform_inputs()
    c0a = np.zeros(N_STATE_AGE)
    for d in range(N_CLUSTERS):
        tot = 1.0 if d == 0 else 0.0
        for si in (0, 1):
            base = astate_index(d, si, 0)
            c0a[base:base + N_AGE] = 0.5 * tot / N_AGE
    _, Y = solver.run_age_scenario(
        c0a, t_end=50.0, dt=1.0, I=I_FLOW, beta_a=beta_a, mu_pat_a=mu_pat,
        mu_mat_a=mu_mat, eps=EPS, alpha=ALPHA, sigma=0.5, s=0.5,
        rule="deepest_inheritance", imm_age=imm)
    Y3 = Y.reshape(Y.shape[0], 2, N_CLUSTERS, N_AGE)
    sym = np.abs(Y3[:, 0] - Y3[:, 1]).max() / Y3.sum(axis=(1, 2, 3)).max()
    print(f"pat/mat symmetry max rel diff = {sym:.3e}")
    assert sym < 1e-9, "symmetry broken"
    print("symmetry: OK")

    print("ALL AGE-CONSISTENCY CHECKS PASSED" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
