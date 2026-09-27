"""Religious-dimension declarations for the genealogy RAG.

Extends the two-sex model: populations are now (sex, religion),
P = {pat, mat} x G with G = 8 religious groups. A cluster state is
(d, sex, g): pedigree depth d = 0..n in the paternal/maternal population of
religious group g; concentration c_{d,sex,g}(t) in millions of people.

State vector c(t) in R^{G*2*(n+1)}, block-ordered: for g in GROUPS,
[pat block C_0..C_n, mat block C_0..C_n].

RAG edge classes added on top of the two-sex model:

5. Religious assortative mating (binary, cross-population):
   (q,pat,g1)+(r,mat,g2) -> parents (catalytic)
   + s*(c(q,r),pat,g2) + (1-s)*(c(q,r),mat,g2)
   at J_{(q,g1),(r,g2)}(c). The child takes the MOTHER's religious group
   (ASSUMED). The pairing first draws the religious pairing (g1,g2) with
   group-specific endogamy rho_g, then the pedigree classes within each
   group-sex with the usual class assortativity alpha.
6. Disaffiliation (inter-population transfer, cf. the paper's RAG reaction
   classes): (d,sex,g) -> (d,sex,una) at rate delta(t) for
   g in {pro, cat, mor, oth}. One-way transfer into the unaffiliated group.

Group endogamy rho_g = share of married people in g whose spouse is also in
g. Values are ESTIMATED from Pew surveys where available (see RHO_SOURCE);
historical constancy of rho_g is ASSUMED.
"""

import numpy as np

from .populations import N_CLUSTERS, MAX_DEPTH

# --- religious groups ---------------------------------------------------------
GROUPS = ("pro", "cat", "mor", "jew", "mus", "dhr", "oth", "una")
GROUP_INDEX = {g: i for i, g in enumerate(GROUPS)}
G = len(GROUPS)
GROUP_DISPLAY = {
    "pro": "Protestant",
    "cat": "Catholic",
    "mor": "Mormon (LDS)",
    "jew": "Jewish",
    "mus": "Muslim",
    "dhr": "Hindu/Buddhist",
    "oth": "Orthodox/other faith",
    "una": "Unaffiliated",
}

N_STATE_R = G * 2 * N_CLUSTERS      # (n+1)*2*G = 240 ODEs at n = 14


def rstate_index(d: int, sex: str, g: str) -> int:
    """Flat state index of cluster state (d, sex, g)."""
    gi = GROUP_INDEX[g] if isinstance(g, str) else int(g)
    si = 0 if sex == "pat" else 1
    return (gi * 2 + si) * N_CLUSTERS + int(d)


# --- group endogamy rho_g ------------------------------------------------------
# Share of married people in g whose spouse shares their religion.
RHO = {
    # ESTIMATED: Pew 2023-24 Religious Landscape Study, intermarriage chapter
    # ("How many married Americans have spouses of the same religion",
    # Feb 2025): 87% of married Latter-day Saints, 81% of married Protestants,
    # 75% of married Catholics, 68% of married unaffiliated have same-religion
    # spouses; 65% of married Jewish respondents (RLS) -- see jew note below.
    "pro": 0.81,
    "cat": 0.75,
    "mor": 0.87,
    "una": 0.68,
    # ESTIMATED: Pew 2020 survey of Jewish Americans: 58% of married Jews have
    # a Jewish spouse (56% in 2013); RLS 2023-24 reports 65%. Use 0.62.
    "jew": 0.62,
    # ESTIMATED: Pew 2017 survey of U.S. Muslims: 84% of married/partnered
    # Muslim Americans have a Muslim spouse/partner.
    "mus": 0.84,
    # ASSUMED: Pew 2008 RLS found 90% of married Hindus with Hindu spouses;
    # no separate Buddhist estimate exists. Combined Hindu/Buddhist group.
    "dhr": 0.85,
    # ASSUMED: no Pew intermarriage estimate; Orthodox Christians are
    # traditionally endogamous; "other faiths" mixed.
    "oth": 0.80,
}
RHO_SOURCE = {
    "pro": "ESTIMATED (Pew 2023-24 RLS)",
    "cat": "ESTIMATED (Pew 2023-24 RLS)",
    "mor": "ESTIMATED (Pew 2023-24 RLS)",
    "una": "ESTIMATED (Pew 2023-24 RLS)",
    "jew": "ESTIMATED (Pew 2020 Jewish Americans; RLS 2023-24)",
    "mus": "ESTIMATED (Pew 2017 Muslim Americans)",
    "dhr": "ASSUMED (Hindu 0.90 Pew 2008 RLS anchor)",
    "oth": "ASSUMED",
}
RHO_VEC = np.array([RHO[g] for g in GROUPS])

# --- immigrant religious composition iota_g(t) --------------------------------
# Era brackets. 1965-2025 refined from Pew 2013 "The Religious Affiliation of
# U.S. Immigrants" (legal immigrants: Christian 68% in 1992 -> 61% in 2012;
# Muslim 5%->10%; Hindu 3%->7%; Buddhist 7%->6%; unaffiliated ~14% stable;
# unauthorized immigrants 83% Christian). Earlier eras are ASSUMED brackets
# consistent with the known immigration history (colonial Protestants; Irish/
# German Catholics 1840s-60s; Ellis Island era with Jewish and Orthodox
# inflows; restriction era).
_IOTA_ERAS = [
    # (start, end, {group: share})
    (1650, 1800, {"pro": 0.93, "cat": 0.05, "oth": 0.02}),                                   # ASSUMED
    (1800, 1860, {"pro": 0.68, "cat": 0.27, "jew": 0.02, "oth": 0.03}),                      # ASSUMED
    (1860, 1920, {"pro": 0.48, "cat": 0.35, "jew": 0.10, "oth": 0.05, "una": 0.02}),         # ASSUMED
    (1920, 1965, {"pro": 0.50, "cat": 0.30, "jew": 0.10, "oth": 0.05, "una": 0.05}),         # ASSUMED
    # ESTIMATED (Pew 2013 immigrant-religion study anchors, mapped to groups)
    (1965, 2100, {"pro": 0.28, "cat": 0.34, "mus": 0.08, "dhr": 0.10,
                  "jew": 0.02, "mor": 0.01, "oth": 0.03, "una": 0.14}),
]
_IOTA_LABEL = {
    (1650, 1800): "ASSUMED", (1800, 1860): "ASSUMED", (1860, 1920): "ASSUMED",
    (1920, 1965): "ASSUMED", (1965, 2100): "ESTIMATED (Pew 2013)",
}


def immigrant_composition(t: float) -> np.ndarray:
    """iota_g(t): religious composition of immigrants in year t; (G,) vector."""
    for start, end, d in _IOTA_ERAS:
        if start <= t < end:
            v = np.zeros(G)
            for g, share in d.items():
                v[GROUP_INDEX[g]] = share
            v /= v.sum()
            return v
    return immigrant_composition(2025.0)


def iota_era_label(t: float) -> str:
    for start, end, _ in _IOTA_ERAS:
        if start <= t < end:
            return _IOTA_LABEL[(start, end)]
    return _IOTA_LABEL[(1965, 2100)]


# --- disaffiliation ------------------------------------------------------------
# Inter-population transfer (d,sex,g) -> (d,sex,una) for g in
# DISAFFILIATING_GROUPS at rate delta(t): 0 before 1970, linear ramp
# 1970-1990 to delta_max, constant after. delta_max is CALIBRATED so the
# unaffiliated share lands near the Pew 2023-24 RLS value (~29%).
# One-way transfer is a crude stand-in for cohort-based switching
# (Pew 2023-24 RLS retention table: raised Catholic 57% stay / 24% now
# unaffiliated; raised Protestant 70%/22%; raised Mormon 54%/28%).
DISAFFILIATING_GROUPS = ("pro", "cat", "mor", "oth")
DISAFFILIATING_IDX = tuple(GROUP_INDEX[g] for g in DISAFFILIATING_GROUPS)
UNA_IDX = GROUP_INDEX["una"]


def disaffiliation_rate(t: float, delta_max: float = 0.007) -> float:
    """delta(t): per-capita disaffiliation rate (1/yr)."""
    if t < 1970.0:
        return 0.0
    if t < 1990.0:
        return delta_max * (t - 1970.0) / 20.0
    return delta_max
