"""Age-structured two-sex cluster model: c_{d,s,a}(t) (Leslie-type).

State: pedigree cluster d = 0..n, sex s in {pat, mat}, 5-year age band
a = 0..17 (0-4, 5-9, ..., 80-84, 85+ open-ended). Flat vector in
R^{2(n+1)*18}, block order [sex: cluster: age], i.e. flat index
(s*(n+1) + d)*18 + a.

Master equation (component form):

    dc_{d,s,a}/dt = births_{d,s} * delta_{a,0}          (newborns -> band 0)
                    + I(t)*sigma_s * pi_a * delta_{d,0} (immigration age profile)
                    + gamma_{a-1} c_{d,s,a-1} - gamma_a c_{d,s,a}   (aging)
                    - (mu_{s,a}(t) + eps) c_{d,s,a}     (mortality+emigration)

Births: B_d split s : (1-s) between the newborns' paternal/maternal
populations of the child cluster c(q,r) given by the UNCHANGED shallowest /
deepest child rules. Mating weights are fertility-weighted over
reproductive ages only:

    w^{pat}_q = sum_a beta_a(t) c_{q,pat,a},   (beta_a = 0 outside 15-49)
    w^{mat}_r = sum_a beta_a(t) c_{r,mat,a},

and J_{qr} follows the same assortative kernel as the aggregate model
(kernels.mating_flux called on the pseudo-state [wpat, wmat] with unit
beta). Fathers are thus drawn from reproductive-age males; the male
schedule is taken identical to the female ASFR schedule (ASSUMED: NCHS
does not publish a male age-specific fertility series).

Aging (documented choice): continuous transfer flux gamma_a = 1/w_a with
w_a = 5 yr for bands 0..16 and gamma_17 = 0 for the open-ended 80+ band.
Within-band age is then exponentially distributed with mean w_a -- the
standard compartmental-ODE approximation. (A discrete annual shift would
require delay-differential equations; the flux form keeps the model an
ODE system in the dc/dt = P + SJ - Dc family, with aging as an extra
linear operator A c.)

Consistency property: with age-UNIFORM rates (beta_a = const for ALL
ages, mu_{s,a} = const, immigration split across bands in fixed shares),
summing over age bands reproduces the aggregate two-sex model EXACTLY
(aging fluxes telescope; see test_age_consistency.py).
"""

import numpy as np

from .populations import N_CLUSTERS, MAX_DEPTH, RULES, child_cluster_matrix
from .kernels import (
    mating_flux, SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION,
)

# --- age grid ---------------------------------------------------------------
N_AGE = 18                     # 5-year bands 0-4 .. 80-84, plus 85+ open
BAND_W = 5.0                   # band width (years) for bands 0..16
AGE_MIDPOINTS = np.array([2.5 + 5 * a for a in range(17)] + [92.5])
# gamma_a: continuous aging flux out of band a (1/yr). Band 17 (85+) is
# open-ended: nobody ages out of it.
AGING_FLUX = np.array([1.0 / BAND_W] * 17 + [0.0])
# Reproductive bands: 15-19 .. 45-49
REPRO_BANDS = tuple(range(3, 10))

N_STATE_AGE = 2 * N_CLUSTERS * N_AGE

# --- immigrant age profile (ASSUMED) ----------------------------------------
# Fixed shares of the immigration inflow by 5-year age band. Stylized
# working-age-concentrated profile (documented assumption; no annual
# age-of-immigrant series is wired up). Future calibration source:
# DHS Yearbook of Immigration Statistics Table 8 (LPRs by age and sex).
IMMIGRANT_AGE_PROFILE = np.array([
    0.03, 0.03, 0.03,                   # 0-14
    0.08, 0.12, 0.13, 0.12, 0.10, 0.08,  # 15-44
    0.06, 0.05, 0.04,                   # 45-59
    0.03, 0.03, 0.02, 0.02, 0.015, 0.015 # 60+
])
assert abs(IMMIGRANT_AGE_PROFILE.sum() - 1.0) < 1e-12
assert len(IMMIGRANT_AGE_PROFILE) == N_AGE


def astate_index(d: int, sex, a: int) -> int:
    """Flat state index of (d, sex, a); sex in {"pat","mat"} or {0,1}."""
    s = {"pat": 0, "mat": 1}[sex] if isinstance(sex, str) else int(sex)
    return (s * N_CLUSTERS + int(d)) * N_AGE + int(a)


def mating_weights_age(C: np.ndarray, beta_a: np.ndarray,
                      beta_pat_a: np.ndarray | None = None):
    """Fertility-weighted class populations: w^{pat}_q, w^{mat}_r.

    C : (2, n, A) array; beta_a : (A,) per-woman birth rates (0 outside
    15-49). The maternal schedule beta_a is the measured ASFR series;
    the paternal schedule beta_pat_a defaults to it (ASSUMED proxy: no
    male ASFR series exists). Returns (wpat, wmat) each (n,).
    """
    beta_a = np.asarray(beta_a, dtype=float)
    beta_pat_a = beta_a if beta_pat_a is None else np.asarray(beta_pat_a,
                                                             dtype=float)
    wpat = (C[0] * beta_pat_a[None, :]).sum(axis=1)
    wmat = (C[1] * beta_a[None, :]).sum(axis=1)
    return wpat, wmat


def birth_inflow_age(C: np.ndarray, beta_a: np.ndarray, alpha: float,
                     s: float, rule: str,
                     beta_pat_a: np.ndarray | None = None):
    """Return (B, births_pat, births_mat): total births and per-cluster
    newborn vectors (deposited into age band 0 by the caller).

    The class-level mating weights come from mating_weights_age; the (q,r)
    flux then uses the identical assortative kernel as the aggregate model
    (mating_flux on the pseudo-state [wpat, wmat] with unit beta).
    """
    n = N_CLUSTERS
    wpat, wmat = mating_weights_age(C, beta_a, beta_pat_a)
    pseudo = np.concatenate([wpat, wmat])
    J, B = mating_flux(pseudo, np.ones(n), np.ones(n), alpha, rule)
    J = J.reshape(n, n)
    child = child_cluster_matrix(rule)
    bpat = np.zeros(n)
    bmat = np.zeros(n)
    np.add.at(bpat, child, s * J)
    np.add.at(bmat, child, (1.0 - s) * J)
    return B, bpat, bmat


def rhs_age(t, c, *, I, beta_a, mu_pat_a, mu_mat_a, eps,
            alpha, sigma=IMMIGRANT_MALE_FRACTION, s=SEX_RATIO_AT_BIRTH,
            rule="shallowest_inheritance", imm_age=IMMIGRANT_AGE_PROFILE,
            male_shift_lam=0.0, beta_pat_a=None):
    """Right-hand side of the age-structured master equation.

    I : immigration inflow (millions/yr). beta_a : (18,) per-woman birth
    rates. mu_pat_a / mu_mat_a : (18,) per-capita death rates by sex.
    eps : emigration rate (uniform). imm_age : (18,) immigrant age shares.
    sigma may be a scalar or a sigma(t) callable.

    Male-fertility sensitivity (step 5): the maternal schedule beta_a is
    the measured ASFR series; the paternal mating-weight schedule
    defaults to it (ASSUMED proxy: no male ASFR series exists). Pass
    beta_pat_a for an explicit male schedule, or male_shift_lam in
    (0,1] to shift fathers ~2 yr older (beta_pat_a :=
    (1-lam)*beta_a + lam*beta_a rolled one 5-yr band). Births are
    maternal; the child split (s, 1-s) and the min-form kernel are
    untouched.
    """
    n = N_CLUSTERS
    A = N_AGE
    C = np.asarray(c, dtype=float).reshape(2, n, A)
    beta_a = np.asarray(beta_a, dtype=float)
    mu_pat_a = np.asarray(mu_pat_a, dtype=float)
    mu_mat_a = np.asarray(mu_mat_a, dtype=float)
    imm_age = np.asarray(imm_age, dtype=float)
    sig = sigma(t) if callable(sigma) else sigma
    if beta_pat_a is None and male_shift_lam > 0.0:
        beta_pat_a = ((1.0 - male_shift_lam) * beta_a
                      + male_shift_lam * np.roll(beta_a, 1))

    dc = np.zeros_like(C)

    # births -> age band 0 of the child cluster (shallowest/deepest rule)
    B, bpat, bmat = birth_inflow_age(C, beta_a, alpha, s, rule,
                                     beta_pat_a=beta_pat_a)
    dc[0, :, 0] += bpat
    dc[1, :, 0] += bmat

    # immigration -> C_0 across age bands (fixed assumed profile)
    dc[0, 0, :] += sig * I * imm_age
    dc[1, 0, :] += (1.0 - sig) * I * imm_age

    # aging: continuous transfer flux gamma_a c_{d,s,a} to band a+1
    out = AGING_FLUX[None, None, :] * C
    dc -= out
    dc[:, :, 1:] += out[:, :, :-1]

    # mortality + emigration sinks (sex- and age-specific)
    dc[0] -= (mu_pat_a[None, :] + eps) * C[0]
    dc[1] -= (mu_mat_a[None, :] + eps) * C[1]

    return dc.reshape(-1), B


def rhs_age_vec(t, c, **kw):
    """Vector-only wrapper of rhs_age for solve_ivp (drops the births B)."""
    dc, _ = rhs_age(t, c, **kw)
    return dc


def rhs_age_td(t, c, *, age_series, alpha, sigma, s, rule,
               imm_age=IMMIGRANT_AGE_PROFILE, male_shift_lam=0.0):
    """Time-dependent driver. age_series(t) -> (I, beta_a, mu_pat_a,
    mu_mat_a, eps). sigma may be a scalar or a sigma(t) callable.
    male_shift_lam: male-fertility sensitivity hook (see rhs_age)."""
    I, beta_a, mu_pat_a, mu_mat_a, eps = age_series(t)
    dc, _ = rhs_age(t, c, I=I, beta_a=beta_a, mu_pat_a=mu_pat_a,
                    mu_mat_a=mu_mat_a, eps=eps, alpha=alpha,
                    sigma=sigma, s=s, rule=rule, imm_age=imm_age,
                    male_shift_lam=male_shift_lam)
    return dc


def total_balance_age(I, B, C, mu_pat_a, mu_mat_a, eps):
    """Exact total-population balance (sums rhs_age over all states):

        d/dt sum_{d,s,a} c = I + B - sum_a mu_{pat,a} C_{pat,a}
                                - sum_a mu_{mat,a} C_{mat,a} - eps*total.

    Aging fluxes telescope to zero and births conserve (each (q,r) column
    contributes s + (1-s) = 1 times its flux), so neither appears.
    """
    C = np.asarray(C, dtype=float).reshape(2, N_CLUSTERS, N_AGE)
    mu_pat_a = np.asarray(mu_pat_a, dtype=float)
    mu_mat_a = np.asarray(mu_mat_a, dtype=float)
    deaths = ((mu_pat_a[None, :] * C[0]).sum()
              + (mu_mat_a[None, :] * C[1]).sum())
    return I + B - deaths - eps * C.sum()


def stable_age_distribution(beta_a, mu_a, newborn_frac=0.5):
    """Stable age distribution (dominant eigenvector, sums to 1).

    Continuous-time 18-state projection matrix A: A[0,a] =
    ``newborn_frac * beta_a`` (sex-specific births into band 0; beta_a is
    births per *woman*, so the female matrix needs newborn_frac = 1 - s),
    A[a+1,a] = gamma_a for a = 0..16 (aging), and
    A[a,a] = -(gamma_a + mu_a) (aging-out + death). The Perron vector is
    the age pyramid of the exponentially-growing stable population; used
    for initial conditions (avoids an age-structure spin-up transient).

    NOTE (2026-09-26): an earlier version used the raw beta_a in the birth
    row, overstating single-sex births ~2x; it returned r = +1.5%/yr for
    2024 (TFR 1.62, R0 = 0.78 < 1), an impossibility. The newborn_frac
    scaling fixes this.
    """
    A = N_AGE
    beta_a = np.asarray(beta_a, dtype=float)
    mu_a = np.asarray(mu_a, dtype=float)
    M = np.zeros((A, A))
    M[0, :] = newborn_frac * beta_a
    for a in range(A - 1):
        M[a + 1, a] = AGING_FLUX[a]
    for a in range(A):
        M[a, a] -= AGING_FLUX[a] + mu_a[a]
    vals, vecs = np.linalg.eig(M)
    i = int(np.argmax(vals.real))
    v = vecs[:, i].real
    if (v < 0).all():
        v = -v
    assert (v > 0).all(), "Perron vector must be strictly positive"
    r = float(vals[i].real)
    return v / v.sum(), r


def mean_age_at_childbearing(beta_a):
    """MAC = sum_a mid_a beta_a / sum_a beta_a (years). beta_a per-woman."""
    beta_a = np.asarray(beta_a, dtype=float)
    tot = beta_a.sum()
    if tot <= 0:
        return float("nan")
    return float((AGE_MIDPOINTS * beta_a).sum() / tot)


def initial_condition_age(total_c0, beta_a, mu_pat_a, mu_mat_a,
                          sigma=IMMIGRANT_MALE_FRACTION, s=0.512):
    """IC: total_c0 (millions) in C_0, sexes split sigma/(1-sigma), each
    sex on its stable age distribution at the given rates; C_1..C_n = 0.

    The female profile uses newborn_frac = 1 - s (daughters per woman);
    the male profile uses newborn_frac = s and the female growth rate
    (male births are driven by the female population). ``s`` defaults to
    the measured NCHS sex ratio at birth.
    """
    prof_mat, r_mat = stable_age_distribution(beta_a, mu_mat_a,
                                              newborn_frac=1.0 - s)
    # Male profile conditional on the female-driven growth rate r_mat: the
    # male subpopulation obeys dM/dt = T M + b(t) with T the male
    # aging/mortality transition matrix and b(t) the newborn-son inflow
    # growing at r_mat. In stable growth M(t) = e^{r t} m, so
    # m = (r_mat I - T)^{-1} e_0 (up to scale), e_0 = newborn inflow vector.
    A = N_AGE
    T = np.zeros((A, A))
    for a in range(A - 1):
        T[a + 1, a] = AGING_FLUX[a]
    for a in range(A):
        T[a, a] -= AGING_FLUX[a] + mu_pat_a[a]
    e0 = np.zeros(A)
    e0[0] = 1.0
    v = np.linalg.solve(r_mat * np.eye(A) - T, e0)
    assert (v > 0).all(), "male conditional age profile must be positive"
    prof_pat = v / v.sum()
    c = np.zeros(N_STATE_AGE)
    c[astate_index(0, "pat", 0):astate_index(0, "pat", 0) + N_AGE] = (
        sigma * total_c0 * prof_pat)
    c[astate_index(0, "mat", 0):astate_index(0, "mat", 0) + N_AGE] = (
        (1.0 - sigma) * total_c0 * prof_mat)
    return c
