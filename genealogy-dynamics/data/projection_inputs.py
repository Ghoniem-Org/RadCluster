#!/usr/bin/env python3
"""CASE 2 projection inputs: 2025 -> 2050 input series.

Every projected input is labeled ASSUMED/PROJECTED -- NEVER measured.
Case 1 (1650-2025 colonial) inputs keep their own provenance labels
(data/inputs.csv, data/age_inputs_provenance.md); this module only covers
the 2025-2050 projection window.

Per-input method choices (frozen-at-2024 vs trend extrapolation):

1. beta(t) [birth rate]: TREND extrapolation. OLS linear fit on the
   measured 2015-2024 crude birth rate (World Bank series in
   data/inputs.csv; 12.4 -> 10.6 /1000, slope -0.20/1000/yr). Justification:
   the 2015-2024 decline is monotonic with no break, so the trend is the
   defensible continuation. beta(t) = 2 * CBR_proj(t)/1000 as a uniform
   (n,) vector, the same mapping the Case-1 colonial runs use (the factor
   2 is the pipeline's per-population convention). Clamped at 0 (physical
   bound only; no ad-hoc floor). PROJECTED.
   Fitted (t in years): CBR_proj(2025) = 10.26, CBR_proj(2050) = 5.26 /1000.

2. mu(t) [death rate]: TREND extrapolation on the measured 2015-2019 plus
   2023-2024 crude death rates, EXCLUDING the COVID outlier years 2020-2022
   (CDR 10.3/10.4/9.8 -- the pandemic mortality spike is not a trend).
   OLS slope +0.075/1000/yr, reflecting continued population aging.
   Justification: a full-window fit would drag the pandemic spike into the
   projection slope. mu(t) = CDR_proj(t)/1000, uniform across clusters and
   sexes (same as Case 1). PROJECTED.
   Fitted: CDR_proj(2025) = 9.20, CDR_proj(2050) = 11.09 /1000
   (2024 measured: 9.0 -- smooth handoff).

3. I(t) [immigration inflow]: FROZEN at the 2024 measured level,
   I = 1.705 M/yr (DHS LPR + unauthorized stock->flow series). Justification:
   the 2021-2024 post-COVID rebound slope (~0.2 M/yr^2) is not sustainable;
   projecting it to 2050 would give an absurd ~7 M/yr. Freezing at the
   latest measured level avoids inventing a decline/rebound path and gives
   an exact rate-continuity handoff from the Case-1 series (which holds
   2024 values through 2025). ASSUMED (frozen-at-2024).

4. eps: emigration 0.00075/yr, CALIBRATED (from the 1980-2024 hindcast),
   unchanged from Case 1. Not re-estimated for the projection window.

5. delta(t) [religious disaffiliation]: continued at the calibrated
   delta_max = 8e-3/yr via genealogy.religion.disaffiliation_rate(t, 0.008)
   (constant for t >= 1990). PROJECTED continuation. Justification: the
   calibration anchored the 2025 unaffiliated share at 28.2% vs Pew's 29%;
   holding the rate is the neutral continuation. Caveat (documented, not
   hidden): the transfer is one-way with no re-affiliation flow, so the
   unaffiliated share keeps growing over the projection; two-way switching
   remains future work.

6. iota_g(t) [immigrant religious composition]: the (1965, 2100) era
   bracket of genealogy.religion.immigrant_composition, already defined
   through 2100 and constant over the projection window. ESTIMATED
   (Pew 2013 immigrant-religion study anchors); pre-1965 brackets are
   ASSUMED (unchanged from Case 1).

7. sigma = sigma_2024 (immigrant male fraction): the 2024 blend from
   data/immigrant_sex_share.csv (step 5), held constant through 2050.
   PROJECTED (not measured).
   s = 0.512 (sex ratio at birth): MEASURED (NCHS), unchanged.

8. mu_pat(t), mu_mat(t) [sex-specific death rates]: the TREND-extrapolated
   total CDR is split by the constant 2024 measured HLD sex ratio
   r(2024) (mu_pat = mu*2r/(1+r), mu_mat = mu*2/(1+r)), preserving the
   projected total level exactly. PROJECTED + the constant-ratio
   assumption.

Rate continuity at the 2025 handoff: the Case-1 rate_series holds 2024
measured values for t >= 2024. The projection series evaluates to
CBR 10.26 (vs 10.6 measured), CDR 9.20 (vs 9.0), I 1.705 (exact) --
no material discontinuity at t = 2025.
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(BASE))  # for genealogy.sex_rates at import

# --- fitted trends (OLS on data/inputs.csv, computed once at import) ---------
_df = pd.read_csv(os.path.join(BASE, "inputs.csv"))
_recent = _df[(_df["year"] >= 2015) & (_df["year"] <= 2024)]


def _ols(x, y):
    A = np.vstack([x, np.ones(len(x))]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(slope), float(intercept)


# beta: CBR trend 2015-2024 (measured World Bank series)
_CBR_SLOPE, _CBR_INTERCEPT = _ols(
    _recent["year"].to_numpy(), _recent["cbr_per_1000"].to_numpy())

# mu: CDR trend excluding COVID outlier years 2020-2022
_cdr = _recent[~_recent["year"].isin((2020, 2021, 2022))]
_CDR_SLOPE, _CDR_INTERCEPT = _ols(
    _cdr["year"].to_numpy(), _cdr["cdr_per_1000"].to_numpy())

# I: frozen at the 2024 measured level
_I_2024 = float(_recent[_recent["year"] == 2024]["immigration_M"].iloc[0])

# eps: calibrated, same as Case 1
_EPS = 0.00075
# delta_max: calibrated continuation
_DELTA_MAX = 0.008

# sigma: the 2024 immigrant male-share blend (data/immigrant_sex_share.csv),
# held constant through the projection window. PROJECTED (not measured).
from genealogy.sex_rates import load_sigma_series
SIGMA_CASE2 = float(load_sigma_series()(2024.0))

# sex-specific mortality: the 2024 measured HLD male/female CDR ratio,
# held constant while the total CDR follows the projected trend.
_MU_R_2024 = float(pd.read_csv(os.path.join(BASE, "sex_specific_mortality.csv"))
                   .query("year == 2024")["sex_ratio_r"].iloc[0])


def projected_mu_pair(t: float):
    """(mu_pat, mu_mat) per-capita/yr: projected total CDR split by the
    constant 2024 measured HLD sex ratio (PROJECTED + constant-ratio
    assumption). The sex-average reproduces projected_cdr(t)/1000 exactly."""
    mu = projected_cdr(t) / 1000.0
    r = _MU_R_2024
    return mu * 2 * r / (1 + r), mu * 2 / (1 + r)

METHODS = {
    "beta": "TREND extrapolation (PROJECTED): OLS linear fit on measured "
            "2015-2024 CBR, beta(t)=2*CBR_proj/1000 uniform; clamped >= 0",
    "mu": "TREND extrapolation (PROJECTED): OLS linear fit on measured "
          "2015-2019 + 2023-2024 CDR (COVID outlier years 2020-2022 excluded); "
          "mu(t)=CDR_proj/1000 uniform",
    "I": f"FROZEN at 2024 measured level (ASSUMED): I={_I_2024:.3f} M/yr",
    "eps": "CALIBRATED 0.00075/yr, unchanged from Case 1",
    "delta": "PROJECTED continuation at calibrated delta_max=8e-3/yr "
             "(one-way flow caveat documented)",
    "iota": "ESTIMATED (Pew 2013) 1965-2100 era bracket, constant over window",
    "sigma": f"PROJECTED sigma_2024={SIGMA_CASE2:.4f} held, from "
             "data/immigrant_sex_share.csv",
    "s": "MEASURED 0.512 (NCHS), unchanged",
    "mu_pair": f"PROJECTED: total CDR trend split by constant 2024 measured "
               f"HLD male/female ratio r={_MU_R_2024:.4f}",
}


def projected_cbr(t: float) -> float:
    """Projected crude birth rate per 1000 (TREND, PROJECTED)."""
    return max(0.0, _CBR_SLOPE * t + _CBR_INTERCEPT)


def projected_cdr(t: float) -> float:
    """Projected crude death rate per 1000 (TREND excl. COVID, PROJECTED)."""
    return max(0.0, _CDR_SLOPE * t + _CDR_INTERCEPT)


def projected_immigration(t: float) -> float:
    """Projected immigration inflow, M/yr (FROZEN at 2024, ASSUMED)."""
    return _I_2024


def projection_rate_series(t: float):
    """(I, beta_pat, beta_mat, mu, eps) for the two-sex model, t in [2025, 2050].

    mu is a (mu_pat, mu_mat) pair (step 5): the projected total CDR split by
    the constant 2024 measured HLD sex ratio. N_CLUSTERS is read from the
    importing process's GENEALOGY_MAX_DEPTH.
    """
    from genealogy.populations import N_CLUSTERS  # local: needs env var set
    beta = 2.0 * projected_cbr(t) / 1000.0
    bvec = np.full(N_CLUSTERS, beta)
    return (projected_immigration(t), bvec, bvec, projected_mu_pair(t), _EPS)


def projection_rate_series_religious(t: float):
    """(I, beta, mu, eps, iota, delta) for the religious model, t in [2025, 2050].

    mu is a (mu_pat, mu_mat) pair (step 5)."""
    from genealogy.populations import N_CLUSTERS
    from genealogy.religion import immigrant_composition, disaffiliation_rate
    beta = 2.0 * projected_cbr(t) / 1000.0
    bvec = np.full(N_CLUSTERS, beta)
    return (projected_immigration(t), bvec, projected_mu_pair(t), _EPS,
            immigrant_composition(t), disaffiliation_rate(t, _DELTA_MAX))


def describe() -> str:
    """Provenance block for Case-2 output directories."""
    lines = ["# Case-2 (2025-2050) input methods: every projected input is",
             "# ASSUMED/PROJECTED -- none is measured.",
             f"# beta:  CBR OLS 2015-2024: slope={_CBR_SLOPE:+.4f}/1000/yr -> "
             f"CBR(2025)={projected_cbr(2025):.2f}, CBR(2050)={projected_cbr(2050):.2f}",
             f"# mu:    CDR OLS excl. 2020-22: slope={_CDR_SLOPE:+.4f}/1000/yr -> "
             f"CDR(2025)={projected_cdr(2025):.2f}, CDR(2050)={projected_cdr(2050):.2f}",
             f"# I:     frozen at 2024 measured level = {_I_2024:.3f} M/yr",
             f"# eps:   calibrated {_EPS}/yr (Case-1 value)",
             f"# delta: continued at calibrated {_DELTA_MAX}/yr (one-way flow caveat)"]
    return "\n".join(lines) + "\n"
