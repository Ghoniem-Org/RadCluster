"""Rate kernels: fertility/mortality schedules and the two-sex mating kernel.

RAG edge classes and their kernels:

1. Immigration (source)        emptyset -> (0, alpha)      at sigma_alpha * I(t)
2. Mating/birth (binary)       (q,pat)+(r,mat) -> parents (catalytic)
                               + s*(c(q,r),pat) + (1-s)*(c(q,r),mat)
                               at J_{qr}(c)
3. Mortality (sink)            (d,alpha) -> emptyset        at mu_alpha(t) * c_{d,alpha}
4. Emigration (sink)           (d,alpha) -> emptyset        at eps_alpha * c_{d,alpha}

Births form across the two populations: the father (drawn from P_pat) and
the mother (drawn from P_mat) mate cross-population, so the mating flux is
indexed by the paternal/maternal pair (q, r). The sex ratio at birth s
splits the births between the paternal and maternal populations of the
child cluster.
"""

import numpy as np

from .populations import N_CLUSTERS, MAX_DEPTH, child_cluster_matrix

# --- stylized uniform rates (per year) ------------------------------------
MORTALITY_RATE = 0.009    # ~9 per 1,000 (illustrative; uniform over clusters)
EMIGRATION_RATE = 0.002   # stylized baseline (uniform across clusters and sexes)

_BETA_0 = 0.035           # stylized per-capita birth rate, cluster C_0
_BETA_N = 0.022           # stylized per-capita birth rate, cluster C_n
                          # (declines with pedigree depth; linear ramp between)

# --- sex parameters --------------------------------------------------------
SEX_RATIO_AT_BIRTH = 0.512     # s: male fraction among US births (MEASURED).
                               # NCHS: ~1.05 male per female at birth, i.e.
                               # 1.05/2.05 ~= 0.512 (cite NCHS where used).
IMMIGRANT_MALE_FRACTION = 0.5  # legacy default only (ASSUMED);
# step 5 supersedes it with data/immigrant_sex_share.csv, sigma(t).
# Kept for exact-regression tests and as the rhs() default.
                               # (ASSUMED; DHS sex-detail tables named as the
                               # future calibration source).


def fertility_schedule(beta0=_BETA_0, betan=_BETA_N) -> np.ndarray:
    """Per-capita birth rate of one population's cluster C_d (linear ramp,
    beta0 at C_0 down to betan at C_n). Default: same schedule for both
    populations; pass explicit beta_pat/beta_mat arrays to differentiate."""
    return np.linspace(beta0, betan, N_CLUSTERS)


def mating_flux(c: np.ndarray, beta_pat: np.ndarray, beta_mat: np.ndarray,
                alpha: float, rule: str):
    """Cross-sex mating flux J(c) and total births B(c).

    c : state vector in R^{2(n+1)}, block order [pat block, mat block].
    J is returned as a flat (n+1)^2 vector indexed by (q, r) with flat
    index q*(n+1)+r, matching the column layout of S = stoichiometry_matrix.

        J_{qr} = B * [ alpha * delta_{qr} * m_q^{pat}
                       + (1-alpha) * m_q^{pat} * m_r^{mat} ]
        B      = min( sum beta_d^{pat} c_{d,pat}, sum beta_d^{mat} c_{d,mat} )

    m_q^{pat} = beta_q^{pat} c_{q,pat} / sum_p beta_p^{pat} c_{p,pat} is the
    paternal mating weight (similarly maternal); a zero denominator returns
    zero weights rather than 0/0. The assortative term delta_{qr} m_q^{pat}
    means that with probability alpha the father's cluster equals the
    mother's; otherwise the parents are drawn independently. B is a
    marriage-function-style total birth rate: each birth needs one father
    and one mother, so the realized births are limited by the scarcer sex's
    total fertility-weighted population.
    """
    n = N_CLUSTERS
    cpat = c[:n]
    cmat = c[n:]
    wpat = beta_pat * cpat
    wmat = beta_mat * cmat
    Wpat = float(wpat.sum())
    Wmat = float(wmat.sum())
    J = np.zeros(n * n)
    if Wpat <= 0.0 or Wmat <= 0.0:
        return J, 0.0
    B = min(Wpat, Wmat)
    mpat = wpat / Wpat
    mmat = wmat / Wmat
    J = (B * (alpha * np.diag(mpat) + (1.0 - alpha) * np.outer(mpat, mmat))).ravel()
    return J, B


def birth_inflow(c: np.ndarray, beta_pat: np.ndarray, beta_mat: np.ndarray,
                 alpha: float, s: float, rule: str):
    """Return (B, births): total births B and the births vector in R^{2(n+1)}.

    births = S J(c), with S the stoichiometry matrix (splits each (q,r)
    pairing's births s : (1-s) between the paternal and maternal populations
    of the child cluster c(q,r)). Parents are catalytic and not consumed.
    """
    from .populations import stoichiometry_matrix  # local import: avoids cycle
    J, B = mating_flux(c, beta_pat, beta_mat, alpha, rule)
    births = stoichiometry_matrix(s, rule) @ J
    return B, births


def mating_flux_religious(c: np.ndarray, beta: np.ndarray, alpha: float,
                          rho: np.ndarray, rule: str):
    """Mating flux with group-specific religious assortativity.

    c : flat state vector in R^{G*2*(n+1)}, block order
        [g: pat block C_0..C_n, mat block C_0..C_n] for g in GROUPS.
    beta : (n,) per-capita fertility schedule (uniform across groups).
    rho : (G,) group endogamy = share of married people in g whose spouse
        is also in g.

    Availability: A^M_{q,g} = beta_q * c_{q,pat,g}; A^M_g = sum_q A^M_{q,g};
    A^M = sum_g A^M_g (same for F). B = min(A^M, A^F).

    Religious pairing (male-side driven):
        P(g1,g2) = pi^M_{g1} * [rho_{g1} * delta_{g1,g2}
                                + (1 - rho_{g1}) * pi^F_{g2}],
        pi^M_{g1} = A^M_{g1}/A^M,  pi^F_{g2} = A^F_{g2}/A^F.

    Class weights within group-sex: w^M_{q|g1} = A^M_{q,g1}/A^M_{g1}
    (guard 0/0 -> 0).

        J_{(q,g1),(r,g2)} = B * P(g1,g2)
            * [alpha * delta_{qr} * w^M_{q|g1}
               + (1-alpha) * w^M_{q|g1} * w^F_{r|g2}]

    Returns J as a (G, G, n, n) array and B. Normalization:
    sum_{q,r,g1,g2} J = B exactly (each bracket sums to 1), so the flux is
    conservative. Zero-denominator guards return zero flux rather than 0/0.
    """
    from .religion import G  # local import: avoids cycle
    n = N_CLUSTERS
    C = np.asarray(c, dtype=float).reshape(G, 2, n)   # [g, sex, d]
    beta = np.asarray(beta, dtype=float)
    AM = beta[:, None] * C[:, 0, :].T     # (n, G): A^M_{q,g}
    AF = beta[:, None] * C[:, 1, :].T     # (n, G): A^F_{r,g}
    AMg = AM.sum(axis=0)                  # (G,)
    AFg = AF.sum(axis=0)
    AMtot = float(AMg.sum())
    AFtot = float(AFg.sum())
    J = np.zeros((G, G, n, n))
    if AMtot <= 0.0 or AFtot <= 0.0:
        return J, 0.0
    B = min(AMtot, AFtot)
    piM = AMg / AMtot
    piF = AFg / AFtot
    rho = np.asarray(rho, dtype=float)
    for g1 in range(G):
        if AMg[g1] <= 0.0:
            continue                      # no marriageable males in g1
        wM = AM[:, g1] / AMg[g1]           # (n,)
        for g2 in range(G):
            if AFg[g2] <= 0.0:
                continue                  # no marriageable females in g2
            p = piM[g1] * (rho[g1] * (1.0 if g1 == g2 else 0.0)
                           + (1.0 - rho[g1]) * piF[g2])
            if p <= 0.0:
                continue
            wF = AF[:, g2] / AFg[g2]       # (n,)
            J[g1, g2] = (B * p * (alpha * np.diag(wM)
                                  + (1.0 - alpha) * np.outer(wM, wF)))
    return J, B
