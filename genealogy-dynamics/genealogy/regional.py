"""Regional (Census-region) compartments on the two-sex cluster model.

State: c_{d,s,r}(t), d = 0..n pedigree cluster, s in {pat, mat} sex,
r in {NE, MW, S, W} Census region. At n = MAX_DEPTH = 14 this is
15 x 2 x 4 = 120 ODE states. Flat order is region-major:
index = r * 2*(n+1) + a * (n+1) + d, a = 0 (pat), 1 (mat).

Master equation (matrix form dc/dt = P + S J - D c + T c):

    dc_{d,s,r}/dt = P_{d,s,r}(t) + sum_{(q,r') pairings in r} S J_r(c_{.,.,r})
                    - (mu(t) + eps) c_{d,s,r}
                    + sum_{r' != r} m_{r'->r}(t) c_{d,s,r'}
                    - sum_{r' != r} m_{r->r'}(t) c_{d,s,r}

- P: immigration enters C_0 of each region: P_{0,pat,r} = sigma I(t) lam_r(t),
  P_{0,mat,r} = (1-sigma) I(t) lam_r(t); lam_r(t) = regional immigration
  share (MEASURED from DHS LPR state-of-residence, 2014-2023; ASSUMED before).
- Mating is WITHIN-REGION (region assortativity lambda_region = 1, ASSUMED:
  partners are drawn from the same region's paternal/maternal pools; the
  cluster assortativity alpha is applied within each region exactly as in
  the national model). Cross-region unions are future work.
- S is the same (n+1)^2-column stoichiometry matrix per region; parents are
  catalytic. Births conserve people within each region.
- mu(t), eps are the NATIONAL schedules applied uniformly across regions
  (documented choice; regional vital-rate differences are future work).
- T: inter-regional migration TRANSFER operator (linear, conservative):
  per-capita rates m_{r->r'}(t) from Census ACS state-to-state flows
  (MEASURED 2023/2024; ASSUMED before). Applies equally to all (d, s)
  (no cluster/sex-selective migration data).

Total-balance identity (exact): sum_{d,s,r} dc/dt = I + B - (mu+eps) * total,
with B = sum_r B_r. The transfer operator is exactly conservative:
sum_r T_r = 0 term by term.
"""

import numpy as np

from .populations import N_CLUSTERS, MAX_DEPTH, RULES, state_index
from .kernels import (
    mating_flux, SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION,
)

REGIONS = ("NE", "MW", "S", "W")
REGION_FULL = ("Northeast", "Midwest", "South", "West")
REGION_INDEX = {r: i for i, r in enumerate(REGIONS)}
N_REGIONS = 4

N_STATE_REGIONAL = N_REGIONS * 2 * N_CLUSTERS  # 120 at n=14

SEXES = ("pat", "mat")


def rstate_index(d: int, sex, region) -> int:
    """Flat index of state (d, sex, region); sex in {"pat","mat"} or {0,1};
    region in {"NE","MW","S","W"} or {0..3}."""
    a = {"pat": 0, "mat": 1}[sex] if isinstance(sex, str) else int(sex)
    r = REGION_INDEX[region] if isinstance(region, str) else int(region)
    return r * 2 * N_CLUSTERS + a * N_CLUSTERS + int(d)


def region_block(c: np.ndarray, region) -> np.ndarray:
    """The 2*(n+1) national-order state vector of one region (a view)."""
    r = REGION_INDEX[region] if isinstance(region, str) else int(region)
    n = N_CLUSTERS
    return c[r * 2 * n:(r + 1) * 2 * n]


def source_vector_regional(I: float, lam, sigma=IMMIGRANT_MALE_FRACTION):
    """P(t) in R^{4*2*(n+1)}: immigration into C_0 of each region by share lam_r.

    P_{(0,pat,r)} = sigma * I(t) * lam_r(t),
    P_{(0,mat,r)} = (1-sigma) * I(t) * lam_r(t); else 0.
    lam : (4,) regional shares (sums to 1).
    """
    P = np.zeros(N_STATE_REGIONAL)
    lam = np.asarray(lam, dtype=float)
    for ri in range(N_REGIONS):
        P[rstate_index(0, "pat", ri)] = sigma * I * lam[ri]
        P[rstate_index(0, "mat", ri)] = (1.0 - sigma) * I * lam[ri]
    return P


def migration_transfer_rates(M):
    """(4, 4) per-capita out-migration rate matrix m_{r->r'} (1/yr).

    M[i, j] = rate at which people in region i move to region j; M[i, i] = 0.
    The transfer operator on the flat state is
        (T c)_{d,s,r} = sum_{r'} M[r', r] c_{d,s,r'} - sum_{r'} M[r, r'] c_{d,s,r}.
    """
    return np.asarray(M, dtype=float)


def rhs_regional(t, c, *, I, lam, beta_pat, beta_mat, alpha, mu, eps, M,
                 sigma=IMMIGRANT_MALE_FRACTION, s=SEX_RATIO_AT_BIRTH,
                 rule="shallowest_inheritance"):
    """Right-hand side of the regional master equation.

    I : national immigration (M/yr); lam : (4,) regional shares.
    beta_pat/beta_mat : (n+1,) fertility schedules (national, applied to
    every region); alpha : cluster assortativity (within-region);
    mu, eps : national mortality/emigration; M : (4,4) migration-rate matrix.
    """
    from .populations import stoichiometry_matrix  # local: avoids cycle
    n = N_CLUSTERS
    C = np.asarray(c, dtype=float).reshape(N_REGIONS, 2, n)  # [r, sex, d]
    S = stoichiometry_matrix(s, rule)

    dc = np.zeros_like(C)
    B_total = 0.0
    for ri in range(N_REGIONS):
        flat = C[ri].reshape(2 * n)          # block order [pat, mat]
        J, B = mating_flux(flat, beta_pat, beta_mat, alpha, rule)
        B_total += B
        dc[ri] += (S @ J).reshape(2, n)
    # immigration source into C_0 of each region
    lam = np.asarray(lam, dtype=float)
    sig = sigma(t) if callable(sigma) else sigma
    dc[:, 0, 0] += sig * I * lam
    dc[:, 1, 0] += (1.0 - sig) * I * lam
    # sex-specific mortality + emigration sinks (national schedules,
    # uniform by region; step 5). mu may be a scalar or a
    # (mu_pat, mu_mat) pair; emigration eps stays uniform by sex.
    from .model import _as_pair  # local: avoids cycle
    mu_pat, mu_mat = _as_pair(mu, "mu")
    dc[:, 0, :] -= (mu_pat + eps) * C[:, 0, :]
    dc[:, 1, :] -= (mu_mat + eps) * C[:, 1, :]
    # inter-regional migration transfers (conservative)
    M = np.asarray(M, dtype=float)
    out_rate = M.sum(axis=1)                  # (4,) total out-rate per region
    inflow = (M.T @ C.reshape(N_REGIONS, -1)).reshape(N_REGIONS, 2, n)
    dc += inflow - out_rate[:, None, None] * C
    return dc.reshape(N_STATE_REGIONAL), B_total


def rhs_regional_td(t, c, *, rate_series, alpha, sigma, s, rule):
    """Time-dependent wrapper.

    rate_series : callable(t) -> (I, lam, beta_pat, beta_mat, mu, eps, M).
    Returns the rhs vector only (B_total discarded).
    """
    I, lam, beta_pat, beta_mat, mu, eps, M = rate_series(t)
    dcdt, _ = rhs_regional(t, c, I=I, lam=lam, beta_pat=beta_pat,
                           beta_mat=beta_mat, alpha=alpha, mu=mu, eps=eps,
                           M=M, sigma=sigma, s=s, rule=rule)
    return dcdt


def total_balance_regional(I, B, total, mu, eps):
    """Exact total-population balance for the regional model.

    With sex-specific mortality (step 5):
        d/dt sum_{d,s,r} c = I + B - (mu_pat+eps)*total_pat
                                       - (mu_mat+eps)*total_mat.
    Births conserve within each region (columns of S sum to 1) and the
    migration transfer operator is exactly conservative (sum_r T_r = 0).

    mu may be a scalar (uniform, old behavior) or a (mu_pat, mu_mat) pair;
    total may likewise be a scalar or a (total_pat, total_mat) pair.  A
    scalar total is only valid with uniform (scalar) mu, where the balance
    reduces to I + B - (mu+eps)*total; with sex-specific mu a scalar total
    raises instead of double-counting the sink.
    """
    from .model import _as_pair  # local: avoids cycle
    mu_pat, mu_mat = _as_pair(mu, "mu")
    a = np.asarray(total, dtype=float)
    if a.ndim == 0:
        if mu_pat != mu_mat:
            raise ValueError(
                "total must be a (total_pat, total_mat) pair when mu is "
                "sex-specific; a scalar total cannot be split across sexes")
        return I + B - (mu_pat + eps) * a.item()
    tpat, tmat = _as_pair(total, "total")
    return I + B - (mu_pat + eps) * tpat - (mu_mat + eps) * tmat


def region_totals(Y):
    """(steps, 4) regional population totals from a trajectory."""
    Y3 = Y.reshape(Y.shape[0], N_REGIONS, 2 * N_CLUSTERS)
    return Y3.sum(axis=2)


def cluster_totals_regional(Y):
    """(steps, n+1) cluster totals summed over sexes AND regions."""
    Y4 = Y.reshape(Y.shape[0], N_REGIONS, 2, N_CLUSTERS)
    return Y4.sum(axis=(1, 2))


def region_cluster_shares(Y):
    """(steps, 4, n+1) per-region cluster shares (cluster totals / region total)."""
    Y4 = Y.reshape(Y.shape[0], N_REGIONS, 2, N_CLUSTERS)
    T = Y4.sum(axis=(2, 3))                       # (steps, 4)
    C = Y4.sum(axis=2)                            # (steps, 4, n+1)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(T[:, :, None] > 0, C / T[:, :, None], 0.0)
    return out
