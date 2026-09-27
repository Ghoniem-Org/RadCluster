#!/usr/bin/env python3
"""Verification for the regional compartments model (Step 2).

(a) Aggregation: with inter-regional migration OFF and CONSTANT regional
    shares, the region-aggregated trajectory must reproduce the national
    two-sex model to tight tolerance (the regional ODE is exactly the
    national ODE scaled per region when rates are region-uniform).
(b) Conservation: with migration ON, the exact total-balance identity
    sum(rhs) == I + B - (mu+eps)*total holds (migration is conservative).
(c) Non-negativity on a long run (solver-enforced) and on random states.
(d) Existing suites still pass: test_symmetry.py,
    check_religious_aggregation.py (run separately).
"""
import numpy as np
import pandas as pd

from genealogy.populations import (
    N_CLUSTERS, MAX_DEPTH, RULES, N_STATE, state_index,
)
from genealogy.kernels import (
    fertility_schedule, MORTALITY_RATE,
    SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION,
)
from genealogy.model import rhs_td as rhs_td_national, total_balance
from genealogy.solver import run_hindcast, cluster_totals, sex_totals
from genealogy.regional import (
    REGIONS, N_REGIONS, rstate_index, rhs_regional_td, rhs_regional,
    total_balance_regional, region_totals, cluster_totals_regional,
    region_cluster_shares,
)
from run_colonial_regional import (
    run_regional, lam_at, M_at, EPS, S, SIGMA,
)

BASE = "/home/hatch/workspace/genealogy_dynamics"

# --- national rate series (same as run_colonial.py) ---------------------------
inputs = pd.read_csv(f"{BASE}/data/inputs.csv").sort_values("year").reset_index(drop=True)
_years = inputs["year"].to_numpy()
_I = inputs["immigration_M"].to_numpy()
_beta = (2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0)
_mu = inputs["cdr_per_1000"].to_numpy() / 1000.0


def national_rate_series(t):
    i = int(np.clip(np.searchsorted(_years, t, side="right") - 1, 0, len(_years) - 1))
    bvec = np.full(N_CLUSTERS, _beta[i])
    return _I[i], bvec, bvec, _mu[i], EPS


# --- regional rate series with migration OFF and CONSTANT shares --------------
LAM_CONST = lam_at(2023.0)          # fixed DHS-2023 regional shares
M_ZERO = np.zeros((4, 4))


def regional_nomig_constant(t):
    I_i, bpat, bmat, mu_i, eps_i = national_rate_series(t)
    return I_i, LAM_CONST, bpat, bmat, mu_i, eps_i, M_ZERO


SIGMA0 = SIGMA(1650.0) if callable(SIGMA) else SIGMA  # step 5: callable series
c0_nat = np.zeros(N_STATE)
c0_nat[state_index(0, "pat")] = SIGMA0 * 0.050368
c0_nat[state_index(0, "mat")] = (1.0 - SIGMA0) * 0.050368

c0_reg = np.zeros(4 * 2 * N_CLUSTERS)
for ri in range(N_REGIONS):
    c0_reg[rstate_index(0, "pat", ri)] = SIGMA0 * 0.050368 * LAM_CONST[ri]
    c0_reg[rstate_index(0, "mat", ri)] = (1.0 - SIGMA0) * 0.050368 * LAM_CONST[ri]


def _rk4_step(f, t, y, dt):
    k1 = f(t, y)
    k2 = f(t + dt / 2, y + dt / 2 * k1)
    k3 = f(t + dt / 2, y + dt / 2 * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def check_aggregation():
    print("== (a) migration-OFF aggregation vs national ==")
    print("  (fixed-step RK4, dt=0.25, identical grids: isolates the model,")
    print("   not the adaptive stepper)")
    from functools import partial
    worst = {}
    for rule in RULES:
        for alpha in (0.0, 0.8):
            fn = partial(rhs_td_national, rate_series=national_rate_series,
                         alpha=alpha, sigma=SIGMA, s=S, rule=rule)
            fr = partial(rhs_regional_td, rate_series=regional_nomig_constant,
                         alpha=alpha, sigma=SIGMA, s=S, rule=rule)
            yn, yr = c0_nat.copy(), c0_reg.copy()
            t = 1650.0
            dev = 0.0
            while t < 2025.0 - 1e-12:
                dt = min(0.25, 2025.0 - t)
                yn = _rk4_step(fn, t, yn, dt)
                yr = _rk4_step(fr, t, yr, dt)
                t += dt
                Tn = yn[:N_CLUSTERS] + yn[N_CLUSTERS:]
                Tr = cluster_totals_regional(yr[None, :])[0]
                dev = max(dev, np.abs(Tn - Tr).max())
            worst[(rule, alpha)] = dev
            print(f"  [{rule} a={alpha}] max|cluster dev| over 1650-2025: "
                  f"{dev:.3e} M")
            assert dev < 1e-6, f"aggregation FAILED: {dev:.3e}"
    print("  (a) PASSED")


def check_balance_migration_on():
    print("== (b) total-balance identity with migration ON ==")
    from run_colonial_regional import rate_series
    rng = np.random.default_rng(7)
    worst = 0.0
    for rule in RULES:
        for _ in range(5):
            c = rng.random(4 * 2 * N_CLUSTERS) * 3.0
            t = float(rng.uniform(1650, 2025))
            I_i, lam_i, bpat, bmat, mu_i, eps_i, M_i = rate_series(t)
            rsum = rhs_regional_td(t, c, rate_series=rate_series, alpha=0.8,
                                   sigma=SIGMA, s=S, rule=rule).sum()
            _, B_i = rhs_regional(t, c, I=I_i, lam=lam_i, beta_pat=bpat,
                                  beta_mat=bmat, alpha=0.8, mu=mu_i,
                                  eps=eps_i, M=M_i, sigma=SIGMA, s=S,
                                  rule=rule)
            pred = total_balance_regional(
                I_i, B_i, (c.reshape(4, 2, N_CLUSTERS)[:, 0, :].sum(),
                           c.reshape(4, 2, N_CLUSTERS)[:, 1, :].sum()),
                mu_i, eps_i)
            worst = max(worst, abs(rsum - pred))
            assert abs(rsum - pred) < 1e-8 * max(1.0, abs(pred))
    print(f"  worst |sum(rhs) - (I+B-(mu+eps)P)| = {worst:.3e}")
    print("  (b) PASSED")


def check_nonneg():
    print("== (c) non-negativity ==")
    # the solver raises if any state < -1e-6, so every run_regional call
    # (verification + colonial driver) is non-negativity checked.
    # Additionally: a short regional run from a random state stays >= 0.
    from run_colonial_regional import rate_series
    rng = np.random.default_rng(11)
    c0 = rng.random(4 * 2 * N_CLUSTERS) * 3.0
    t, Y = run_regional(c0, 0.8, "deepest_inheritance", dt=1.0,
                        t0=2000.0, t1=2025.0, rate_series_fn=rate_series)
    assert Y.min() >= -1e-6, Y.min()
    print(f"  min state over 2000-2025 random-IC run: {Y.min():.3e}")
    print("  (c) PASSED")


def main():
    check_aggregation()
    check_balance_migration_on()
    check_nonneg()
    print("ALL REGIONAL VERIFICATION CHECKS PASSED")


if __name__ == "__main__":
    main()
