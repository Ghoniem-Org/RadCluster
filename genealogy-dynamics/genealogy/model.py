"""Two-sex master equations: dc/dt = P(t) + S J(c) - D(t) c.

Cluster states (d, alpha) with alpha in {pat, mat}; state vector c(t) in
R^{2(n+1)}, block-ordered [c_{0,pat}..c_{n,pat}, c_{0,mat}..c_{n,mat}].

Component form:

    dc_{0,pat}/dt = sigma*I - (mu_pat + eps_pat)*c_{0,pat} + s   * B_0(c)
    dc_{0,mat}/dt = (1-sigma)*I - (mu_mat + eps_mat)*c_{0,mat} + (1-s)* B_0(c)
    dc_{d,pat}/dt = -(mu_pat + eps_pat)*c_{d,pat} + s   * B_d(c),  d = 1..n
    dc_{d,mat}/dt = -(mu_mat + eps_mat)*c_{d,mat} + (1-s)* B_d(c),  d = 1..n

where B_d(c) = sum over (q,r) with child_cluster(q,r)=d of J_{qr}(c).

Each piece of the matrix form is exposed as a function so the model can be
inspected (and the report can display it) term by term:

    P(t)  = source_vector(I, sigma)            -> R^{2(n+1)}
    J(c)  = mating_flux(c, beta_pat, beta_mat, alpha, rule)  -> R^{(n+1)^2}
    S     = stoichiometry_matrix(s, rule)      -> R^{2(n+1) x (n+1)^2}
    D(t)  = loss_matrix(mu, eps)               -> R^{2(n+1) x 2(n+1)} diagonal
"""

import numpy as np

from .populations import (
    N_CLUSTERS, N_STATE, MAX_DEPTH, POPULATIONS, POP_INDEX, state_index,
    stoichiometry_matrix,
)
from .kernels import mating_flux, birth_inflow, SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION


def _as_pair(x, name):
    """Scalar -> (x, x) for both sexes; length-2 -> as given (allows
    sex-specific rates while keeping uniform the default)."""
    a = np.asarray(x, dtype=float)
    if a.ndim == 0:
        return a.item(), a.item()
    if a.shape == (2,):
        return a[0].item(), a[1].item()
    raise ValueError(f"{name} must be a scalar or a length-2 (pat, mat) array")


def source_vector(I, sigma=IMMIGRANT_MALE_FRACTION):
    """P(t) in R^{2(n+1)}: immigration splits into the two populations.

    P_{(0,pat)} = sigma * I(t),  P_{(0,mat)} = (1-sigma) * I(t),
    all other entries zero. sigma defaults to the legacy assumed value;
    the production series is data/immigrant_sex_share.csv (step 5), passed
    as a scalar or a sigma(t) callable by the drivers.
    """
    P = np.zeros(N_STATE)
    P[state_index(0, "pat")] = sigma * I
    P[state_index(0, "mat")] = (1.0 - sigma) * I
    return P


def loss_matrix(mu, eps):
    """D(t): diagonal 2(n+1) x 2(n+1) loss matrix.

    D_{(d,alpha),(d,alpha)} = mu_alpha(t) + eps_alpha (mortality + emigration).
    mu and eps may be scalars (uniform across sexes, the default) or
    length-2 (pat, mat) arrays for sex-specific sinks.
    """
    mu_pat, mu_mat = _as_pair(mu, "mu")
    eps_pat, eps_mat = _as_pair(eps, "eps")
    d = np.empty(N_STATE)
    d[:N_CLUSTERS] = mu_pat + eps_pat      # paternal block
    d[N_CLUSTERS:] = mu_mat + eps_mat      # maternal block
    return np.diag(d)


def loss_rates(mu, eps):
    """Per-state loss rates as a vector (the diagonal of D)."""
    mu_pat, mu_mat = _as_pair(mu, "mu")
    eps_pat, eps_mat = _as_pair(eps, "eps")
    d = np.empty(N_STATE)
    d[:N_CLUSTERS] = mu_pat + eps_pat
    d[N_CLUSTERS:] = mu_mat + eps_mat
    return d


def rhs(t, c, *, I, beta_pat, beta_mat, alpha, mu, eps,
        sigma=IMMIGRANT_MALE_FRACTION, s=SEX_RATIO_AT_BIRTH,
        rule="shallowest_inheritance"):
    """Right-hand side: dc/dt = P + S J(c) - D c."""
    B, births = birth_inflow(c, beta_pat, beta_mat, alpha, s, rule)
    return source_vector(I, sigma) + births - loss_rates(mu, eps) * c


def rhs_td(t, c, *, rate_series, alpha, sigma, s, rule):
    """Right-hand side for time-dependent rates (colonial/hindcast drivers).

    rate_series : callable(t) -> (I, beta_pat, beta_mat, mu, eps).
    alpha : float (constant assortativity, the default) or a callable
        alpha(t) -> float for time-varying assortativity
        (see genealogy.assortativity.alpha_schedule).
    """
    a = alpha(t) if callable(alpha) else alpha
    I, beta_pat, beta_mat, mu, eps = rate_series(t)
    sig = sigma(t) if callable(sigma) else sigma
    return rhs(t, c, I=I, beta_pat=beta_pat, beta_mat=beta_mat, alpha=a,
               mu=mu, eps=eps, sigma=sig, s=s, rule=rule)


def total_balance(I, B, total_pat, total_mat, mu, eps):
    """Checkable exact total-population balance (sums dc/dt over all states):

        d/dt sum_{d,alpha} c_{d,alpha}
            = I + B - (mu_pat+eps_pat)*total_pat - (mu_mat+eps_mat)*total_mat

    which reduces to I + B - (mu+eps)*total when sinks are uniform across
    sexes. The mating flux conserves people: sum(S J) = B exactly, because
    each column of S sums to s + (1-s) = 1.
    """
    mu_pat, mu_mat = _as_pair(mu, "mu")
    eps_pat, eps_mat = _as_pair(eps, "eps")
    return (I + B
            - (mu_pat + eps_pat) * total_pat
            - (mu_mat + eps_mat) * total_mat)


# --- religious-dimension master equation --------------------------------------
# State (d, sex, g), c(t) in R^{G*2*(n+1)}, block order [g: pat, mat].
# dc/dt = P(t) + S J(c) - D(t) c + T(t) c, where T(t) is the disaffiliation
# transfer operator (linear, internal: it conserves total people).


def source_vector_religious(I, sigma, iota):
    """P(t) in R^{G*2*(n+1)}: immigration into C_0 of each group.

    P_{(0,pat,g)} = sigma * I(t) * iota_g(t),
    P_{(0,mat,g)} = (1-sigma) * I(t) * iota_g(t); else 0.
    iota : (G,) immigrant religious composition (sums to 1).
    """
    from .religion import G, N_STATE_R, rstate_index  # local: avoids cycle
    P = np.zeros(N_STATE_R)
    iota = np.asarray(iota, dtype=float)
    for gi in range(G):
        P[rstate_index(0, "pat", gi)] = sigma * I * iota[gi]
        P[rstate_index(0, "mat", gi)] = (1.0 - sigma) * I * iota[gi]
    return P


def rhs_religious(t, c, *, I, beta, alpha, mu, eps, sigma, s, rule, rho,
                  iota, delta):
    """Right-hand side of the religious-dimension master equation.

    iota : (G,) immigrant religious composition at time t (driver-provided).
    delta : disaffiliation rate at time t (driver-provided, 1/yr).

    Births: for each pairing (g1,g2), the (q,r) flux J_{(q,g1),(r,g2)} puts
    s*J into (c(q,r),pat,g2) and (1-s)*J into (c(q,r),mat,g2) -- the child
    takes the MOTHER's religious group (ASSUMED). Disaffiliation moves
    delta*c_{d,sex,g} from g in {pro,cat,mor,oth} to (d,sex,una).
    """
    from .religion import (G, N_STATE_R, UNA_IDX, DISAFFILIATING_IDX,
                           rstate_index)
    from .kernels import mating_flux_religious
    from .populations import child_cluster_matrix
    n = N_CLUSTERS
    C = np.asarray(c, dtype=float).reshape(G, 2, n)   # [g, sex, d]
    J, B = mating_flux_religious(c, beta, alpha, rho, rule)
    child = child_cluster_matrix(rule)                # (n,n): (q,r) -> d

    dc = np.zeros_like(C)
    # births: S J accumulation (S has 2 nonzeros per (q,g1,r,g2) column)
    for g1 in range(G):
        for g2 in range(G):
            Jsub = J[g1, g2]
            if Jsub.sum() == 0.0:
                continue
            np.add.at(dc[g2, 0, :], child, s * Jsub)
            np.add.at(dc[g2, 1, :], child, (1.0 - s) * Jsub)
    # immigration source
    iota = np.asarray(iota, dtype=float)
    sig = sigma(t) if callable(sigma) else sigma
    dc[:, 0, 0] += sig * I * iota
    dc[:, 1, 0] += (1.0 - sig) * I * iota
    # sex-specific mortality + emigration sinks (step 5)
    mu_pat, mu_mat = _as_pair(mu, "mu")
    dc[:, 0, :] -= (mu_pat + eps) * C[:, 0, :]
    dc[:, 1, :] -= (mu_mat + eps) * C[:, 1, :]
    # disaffiliation transfers (internal: conserve total people)
    if delta > 0.0:
        for gi in DISAFFILIATING_IDX:
            flow = delta * C[gi]
            dc[gi] -= flow
            dc[UNA_IDX] += flow
    return dc.reshape(N_STATE_R)


def total_balance_religious(I, B, total, mu, eps):
    """Exact total-population balance for the religious model.

    With sex-specific mortality (step 5):
        d/dt sum_{d,sex,g} c = I + B - (mu_pat+eps)*total_pat
                                         - (mu_mat+eps)*total_mat.
    Births conserve (each mating column contributes s + (1-s) = 1 times its
    flux, summed to B) and disaffiliation is an internal transfer, so
    neither appears.

    mu may be a scalar (uniform, old behavior) or a (mu_pat, mu_mat) pair.
    total may likewise be a scalar or a (total_pat, total_mat) pair, with one
    restriction: a scalar total cannot be split across sexes, so it is only
    valid with uniform (scalar) mu, where the balance reduces to
    I + B - (mu+eps)*total.  Passing a scalar total with sex-specific mu
    raises instead of silently double-counting the sink.
    """
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
