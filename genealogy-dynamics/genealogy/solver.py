"""ODE integration driver (scipy solve_ivp) for the two-sex cluster model."""

from functools import partial

import numpy as np
from scipy.integrate import solve_ivp

from .model import rhs, rhs_td
from .populations import N_CLUSTERS


def run_scenario(c0, *, t_end=200.0, dt=1.0, I, beta_pat, beta_mat, alpha,
                 mu, eps, sigma, s, rule):
    """Integrate the two-sex ODEs with constant rates from year 0 to t_end.

    Returns (times, Y) with Y in R^{2(n+1)} (block order pat, mat).
    """
    f = partial(rhs, I=I, beta_pat=beta_pat, beta_mat=beta_mat, alpha=alpha,
                mu=mu, eps=eps, sigma=sigma, s=s, rule=rule)
    sol = solve_ivp(f, (0.0, t_end), np.asarray(c0, dtype=float),
                    t_eval=np.arange(0.0, t_end + dt, dt),
                    method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    if np.any(sol.y < -1e-6):
        raise RuntimeError("negative state detected")
    return sol.t, sol.y.T


def run_hindcast(c0, t0, t1, rate_series, alpha, sigma, s, rule, dt=0.25):
    """Integrate with time-dependent (I, beta_pat, beta_mat, mu, eps) rates."""
    f = partial(rhs_td, rate_series=rate_series, alpha=alpha,
                sigma=sigma, s=s, rule=rule)
    n_steps = int(round((t1 - t0) / dt))
    t_eval = np.linspace(t0, t1, n_steps + 1)
    sol = solve_ivp(f, (t0, t1), np.asarray(c0, dtype=float), t_eval=t_eval,
                    method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    if np.any(sol.y < -1e-6):
        raise RuntimeError("negative state detected")
    return sol.t, sol.y.T


def cluster_totals(Y):
    """Total population per cluster: T_d = c_{d,pat} + c_{d,mat}; (n+1,) rows."""
    return Y[:, :N_CLUSTERS] + Y[:, N_CLUSTERS:]


def sex_totals(Y):
    """(total_pat, total_mat): sum over clusters within each population."""
    return Y[:, :N_CLUSTERS].sum(axis=1), Y[:, N_CLUSTERS:].sum(axis=1)


# --- age-structured model -----------------------------------------------------
# State c(t) in R^{2(n+1)*18}, block order [sex: cluster: age]
# (see genealogy/age.py).

def run_age_scenario(c0, *, t_end=200.0, dt=1.0, I, beta_a, mu_pat_a,
                     mu_mat_a, eps, alpha, sigma, s, rule,
                     imm_age=None):
    """Integrate the age-structured ODEs with constant rates."""
    from .age import rhs_age_vec, N_STATE_AGE, IMMIGRANT_AGE_PROFILE
    if imm_age is None:
        imm_age = IMMIGRANT_AGE_PROFILE
    f = partial(rhs_age_vec, I=I, beta_a=np.asarray(beta_a, dtype=float),
                mu_pat_a=np.asarray(mu_pat_a, dtype=float),
                mu_mat_a=np.asarray(mu_mat_a, dtype=float), eps=eps,
                alpha=alpha, sigma=sigma, s=s, rule=rule, imm_age=imm_age)
    sol = solve_ivp(f, (0.0, t_end), np.asarray(c0, dtype=float),
                    t_eval=np.arange(0.0, t_end + dt, dt),
                    method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    if np.any(sol.y < -1e-6):
        raise RuntimeError("negative state detected")
    assert sol.y.shape[0] == N_STATE_AGE
    return sol.t, sol.y.T


def run_age_hindcast(c0, t0, t1, age_series, alpha, sigma, s, rule,
                     dt=0.25, imm_age=None, male_shift_lam=0.0):
    """Integrate the age-structured ODEs with time-dependent rates.

    age_series : callable(t) -> (I, beta_a, mu_pat_a, mu_mat_a, eps).
    male_shift_lam : male-fertility sensitivity hook (see genealogy.age).
    """
    from .age import rhs_age_td, N_STATE_AGE, IMMIGRANT_AGE_PROFILE
    if imm_age is None:
        imm_age = IMMIGRANT_AGE_PROFILE
    f = partial(rhs_age_td, age_series=age_series, alpha=alpha,
                sigma=sigma, s=s, rule=rule, imm_age=imm_age,
                male_shift_lam=male_shift_lam)
    n_steps = int(round((t1 - t0) / dt))
    t_eval = np.linspace(t0, t1, n_steps + 1)
    sol = solve_ivp(f, (t0, t1), np.asarray(c0, dtype=float), t_eval=t_eval,
                    method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    if np.any(sol.y < -1e-6):
        raise RuntimeError("negative state detected")
    assert sol.y.shape[0] == N_STATE_AGE
    return sol.t, sol.y.T


def cluster_totals_age(Y):
    """Total population per cluster: T_d = sum_{s,a} c_{d,s,a}; (n+1,) rows."""
    from .age import N_AGE
    Y3 = Y.reshape(Y.shape[0], 2, N_CLUSTERS, N_AGE)
    return Y3.sum(axis=(1, 3))


def sex_totals_age(Y):
    """(total_pat, total_mat) summed over clusters and ages."""
    from .age import N_AGE
    Y3 = Y.reshape(Y.shape[0], 2, N_CLUSTERS, N_AGE)
    return Y3[:, 0].sum(axis=(1, 2)), Y3[:, 1].sum(axis=(1, 2))
