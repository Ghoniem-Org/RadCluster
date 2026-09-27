"""Cluster / population declarations for the two-sex genealogy RAG.

Paper language (Ghoniem cluster-dynamics framework, "A Generalized
Graph-Based Cluster Dynamics Framework for Irradiated Materials"):

- Species class C: pedigree clusters, labeled by pedigree depth d = 0..n
  (the "size" analog). n = MAX_DEPTH is a free parameter.
- Populations P = {pat, mat}: P_pat (paternal/male) and P_mat
  (maternal/female), with P_pat, P_mat subset of C.
- Cluster state (d, alpha): a cluster of pedigree depth d in population
  alpha; concentration c_{d,alpha}(t) in millions of people.
- State vector c(t) in R^{2(n+1)}, block-ordered:
  [c_{0,pat}..c_{n,pat}, c_{0,mat}..c_{n,mat}].

This module declares the vertex set V = C x P of the reaction
admissibility graph (RAG). Everything else (kernels, master equations,
solver) is written against this declaration, so adding structure (e.g.
regions) only means extending it.

Three bloodline rules are supported (child cluster given paternal cluster q
and maternal cluster r, with n = MAX_DEPTH):
- "shallowest_inheritance" ("shallowest inheritance"): the child continues
  the SHALLOWER of the two parental bloodlines -- pedigree = number of
  generations with no immigrant ancestor on EITHER side. Every ancestor must
  be US-born, so both parents must be deep for the child to be deep:
  child_cluster = min(1 + min(q, r), n); only deep x deep pairings replenish
  the deep clusters and C_n stays a small, endogamous fraction.
- "deepest_inheritance" ("deepest inheritance"): the child continues the
  DEEPER of the two parental bloodlines:
  child_cluster = min(1 + max(q, r), n). A single deep parent suffices to
  keep the lineage deep, so deep mass ratchets upward and accumulates at
  the cap C_n.
"""

import os

import numpy as np

# Free parameter: pedigree depth n (clusters C_0..C_n). The user assumed
# n = 12 generations back to ~1650, but the depth should be set from the
# actual mean generation interval; see suggested_generations().
# Since Step 6 the default is n = 14 (375 yr / 14 ≈ 26.8 yr per generation),
# aligned with the religious model. Override with GENEALOGY_MAX_DEPTH to
# reproduce n = 12 (or any other depth); regression is pinned at n = 12.
MAX_DEPTH = int(os.environ.get("GENEALOGY_MAX_DEPTH", "14"))
N_CLUSTERS = MAX_DEPTH + 1          # clusters C_0 .. C_n (pedigree depth)

POPULATIONS = ("pat", "mat")        # P_pat (paternal), P_mat (maternal)
POP_INDEX = {"pat": 0, "mat": 1}
N_STATE = 2 * N_CLUSTERS            # state vector c in R^{2(n+1)}

RULES = ("shallowest_inheritance", "deepest_inheritance")
RULE_DISPLAY = {
    "shallowest_inheritance": "shallowest inheritance",
    "deepest_inheritance": "deepest inheritance",
}


def state_index(d: int, alpha) -> int:
    """Flat state index of cluster state (d, alpha); alpha in {"pat","mat"} or {0,1}."""
    a = POP_INDEX[alpha] if isinstance(alpha, str) else int(alpha)
    return a * N_CLUSTERS + int(d)


def suggested_generations(start_year, end_year, years_per_generation):
    """Pedigree depth implied by a timespan and a mean generation interval.

    years_per_generation is the free parameter, to be set from data (mean age
    at childbirth / generation-interval series); e.g. 375 yr / 27 yr ≈ 14.
    """
    return int(round((end_year - start_year) / years_per_generation))


def cluster_description(d: int, alpha: str, rule: str = "shallowest_inheritance") -> str:
    """Human-readable meaning of cluster state (d, alpha) under the rule."""
    n = MAX_DEPTH
    pop = "paternal" if alpha == "pat" else "maternal"
    if d == 0:
        return f"foreign-born (immigrants), {pop}"
    if rule == "shallowest_inheritance":
        if d == n:
            return (f"US-born, every ancestor US-born for {n}+ generations "
                    f"(both sides), {pop}")
        return (f"US-born, every ancestor US-born for {d} generation(s) "
                f"(both sides), {pop}")
    if rule == "deepest_inheritance":
        if d == n:
            return (f"US-born, deepest ancestral bloodline {n}+ generations, "
                    f"{pop}")
        return (f"US-born, deepest ancestral bloodline {d} generation(s), "
                f"{pop}")
    raise ValueError(f"unknown pedigree rule: {rule!r}")


def child_cluster(q: int, r: int, rule: str = "shallowest_inheritance") -> int:
    """Pedigree cluster of a US-born child given paternal cluster q and
    maternal cluster r (n = MAX_DEPTH).

    shallowest_inheritance: min(1 + min(q, r), n). The child continues the
        SHALLOWER parental bloodline: every ancestor must be US-born on both
        sides, so only deep x deep pairings produce deep children.
    deepest_inheritance:    min(1 + max(q, r), n). The child continues the
        DEEPER parental bloodline: a single deep parent keeps the lineage
        deep, so deep mass ratchets upward toward the cap.
    """
    if rule == "shallowest_inheritance":
        return min(1 + min(q, r), MAX_DEPTH)
    if rule == "deepest_inheritance":
        return min(1 + max(q, r), MAX_DEPTH)
    raise ValueError(f"unknown pedigree rule: {rule!r}")


def child_cluster_matrix(rule: str = "shallowest_inheritance") -> np.ndarray:
    """(N_CLUSTERS, N_CLUSTERS) matrix C[q, r] = child cluster for
    paternal cluster q and maternal cluster r."""
    q = np.arange(N_CLUSTERS)[:, None]
    r = np.arange(N_CLUSTERS)[None, :]
    if rule == "shallowest_inheritance":
        return np.minimum(1 + np.minimum(q, r), MAX_DEPTH)
    if rule == "deepest_inheritance":
        return np.minimum(1 + np.maximum(q, r), MAX_DEPTH)
    raise ValueError(f"unknown pedigree rule: {rule!r}")


_STOICH_CACHE: dict = {}


def stoichiometry_matrix(s: float, rule: str) -> np.ndarray:
    """Stoichiometry matrix S of the mating RAG edge class.

    S has shape (2(n+1), (n+1)^2): rows are cluster states (d, alpha) in
    block order [pat block, mat block], columns are ordered parental pairs
    (q, r) with flat index q*(n+1)+r. The mating edge
    (q,pat)+(r,mat) -> (q,pat)+(r,mat) + s*(c(q,r),pat) + (1-s)*(c(q,r),mat)
    is catalytic in the parents, so each column has exactly two nonzero
    entries: S[(c(q,r),pat),(q,r)] = s and S[(c(q,r),mat),(q,r)] = 1-s.
    """
    key = (round(float(s), 12), rule)
    if key not in _STOICH_CACHE:
        n = N_CLUSTERS
        child = child_cluster_matrix(rule).ravel()  # (q,r) -> flat q*n+r
        cols = np.arange(n * n)
        S = np.zeros((2 * n, n * n))
        # accumulate: several (q,r) pairs can map to the same child cluster
        np.add.at(S, (child, cols), s)
        np.add.at(S, (child + n, cols), 1.0 - s)
        _STOICH_CACHE[key] = S
    return _STOICH_CACHE[key]


def initial_condition() -> np.ndarray:
    """Stylized ~2026 US initial condition, in millions. Sums to ~330.

    C_0 = 46 (foreign-born). The remaining 284 is spread over C_1..C_n with
    weights 1..(n-1) for C_1..C_{n-1} and weight 30 for C_n (the capped
    cluster accumulates everyone at the maximum pedigree depth).
    Split evenly between the paternal and maternal populations (symmetric
    IC; the two-sex model reduces exactly to the one-sex model under
    symmetry). These numbers are illustrative, NOT calibrated to census data.

    The same vector is used for both bloodline rules (fair comparison); only
    the cluster *interpretation* changes (see cluster_description).
    """
    P = np.zeros(N_CLUSTERS)
    P[0] = 46.0
    rest = 284.0
    weights = np.arange(N_CLUSTERS, dtype=float)  # weights[d] = d
    weights[0] = 0.0
    weights[MAX_DEPTH] = 30.0                     # capped cluster gets extra mass
    P[1:] = rest * weights[1:] / weights[1:].sum()
    assert abs(P.sum() - 330.0) < 1e-9
    c = np.zeros(N_STATE)
    c[:N_CLUSTERS] = P / 2.0      # paternal block
    c[N_CLUSTERS:] = P / 2.0      # maternal block
    return c
