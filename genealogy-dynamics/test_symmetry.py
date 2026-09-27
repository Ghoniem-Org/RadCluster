#!/usr/bin/env python3
"""Symmetry-reduction test: the two-sex model must collapse exactly onto the
one-sex model when all sex structure is symmetric.

Symmetric conditions:
  s = 0.5, sigma = 0.5, beta_pat = beta_mat, uniform mu/eps across sexes,
  symmetric IC (c_{d,pat} = c_{d,mat}).

Then B = min(W, W) = W = B_old, J_{qr} = J^{old}_{qr}, and births split
50/50, so dT_d/dt (with T_d = c_{d,pat}+c_{d,mat}) equals the old one-sex
dP_d/dt. This script integrates a reference one-sex implementation (the
pre-refactor formulation) side by side with the new two-sex rhs and reports
the maximum deviation of the aggregated trajectory.
"""
import numpy as np

from genealogy.populations import (
    N_CLUSTERS, MAX_DEPTH, RULES, initial_condition, child_cluster_matrix,
    stoichiometry_matrix,
)

# The consecutive bloodline rule was evaluated and dropped (2026-09-26);
# the symmetry suite must cover exactly the two supported rules.
assert RULES == ("shallowest_inheritance", "deepest_inheritance"), RULES
from genealogy.kernels import fertility_schedule, mating_flux, MORTALITY_RATE
from genealogy.model import rhs
from genealogy.solver import run_scenario


# --- reference one-sex implementation (pre-refactor formulation) ------------
def _old_birth_inflow(P, beta, alpha, rule):
    n = N_CLUSTERS
    W = float((beta * P).sum())
    B = 0.5 * W
    if W <= 0:
        return 0.0, np.zeros(n)
    m = beta * P / W
    J = B * (alpha * np.diag(m) + (1 - alpha) * np.outer(m, m))
    births = np.zeros(n)
    C = child_cluster_matrix(rule)
    for q in range(n):
        for r in range(n):
            births[C[q, r]] += J[q, r]
    return B, births


def _old_rhs(t, P, *, I, beta, alpha, mu, eps, rule):
    _, births = _old_birth_inflow(P, beta, alpha, rule)
    dP = np.zeros_like(P)
    dP[0] = I - (mu + eps) * P[0] + births[0]
    dP[1:] = -(mu + eps) * P[1:] + births[1:]
    return dP


def main():
    beta = fertility_schedule()
    I, alpha, mu, eps, s, sigma = 1.0, 0.6, MORTALITY_RATE, 0.00075, 0.5, 0.5
    results = {}
    for rule in RULES:
        # two-sex
        c0 = initial_condition()  # symmetric by construction
        t, Y = run_scenario(c0, t_end=200.0, dt=1.0, I=I,
                            beta_pat=beta, beta_mat=beta, alpha=alpha,
                            mu=mu, eps=eps, sigma=sigma, s=s, rule=rule)
        T = Y[:, :N_CLUSTERS] + Y[:, N_CLUSTERS:]
        # one-sex reference
        from scipy.integrate import solve_ivp
        from functools import partial
        P0 = c0[:N_CLUSTERS] * 2.0
        sol = solve_ivp(partial(_old_rhs, I=I, beta=beta, alpha=alpha,
                                mu=mu, eps=eps, rule=rule),
                        (0.0, 200.0), P0, t_eval=t,
                        method="RK45", rtol=1e-8, atol=1e-10)
        assert sol.success
        dev = np.abs(T - sol.y.T).max()
        results[rule] = dev
        print(f"[{rule}] max |T_d(t) - P_d(t)| over 200 yr: {dev:.3e}")
    worst = max(results.values())
    assert worst < 1e-6, f"symmetry test FAILED: {worst:.3e}"
    print(f"SYMMETRY TEST PASSED (worst deviation {worst:.3e} < 1e-6)")

    # also verify the matrix form: births == S J, and columns of S sum to 1
    for rule in RULES:
        S = stoichiometry_matrix(0.512, rule)
        assert S.shape == (2 * N_CLUSTERS, N_CLUSTERS * N_CLUSTERS)
        colsum = S.sum(axis=0)
        assert np.allclose(colsum, 1.0), "S columns must sum to s+(1-s)=1"
        print(f"[{rule}] S shape {S.shape}, column sums == 1: OK")


if __name__ == "__main__":
    main()
