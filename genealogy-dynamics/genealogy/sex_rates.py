"""Sex-specific rate helpers (step 5 of the genealogy work package).

Loads:
  - data/immigrant_sex_share.csv : sigma(t), male share of the
    immigration inflow (data/build_immigrant_sex_share.py).
  - data/sex_specific_mortality.csv : mu_pat(t), mu_mat(t) (per-capita/yr)
    whose sex-average reproduces the calibrated aggregate death_rate
    (data/build_sex_specific_mortality.py).

All series use the same convention as the rate-series drivers: step
(piecewise-constant) interpolation, value of the greatest year <= t.
"""
import os

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SIGMA_CSV = os.path.join(BASE, "data", "immigrant_sex_share.csv")
_MU_CSV = os.path.join(BASE, "data", "sex_specific_mortality.csv")


def load_sigma_series(path=_SIGMA_CSV):
    """sigma(t) callable from the CSV (step interpolation)."""
    df = pd.read_csv(path).sort_values("year")
    years = df["year"].to_numpy(dtype=float)
    vals = df["sigma"].to_numpy(dtype=float)

    def sigma(t):
        i = int(np.searchsorted(years, t, side="right") - 1)
        return float(vals[max(0, min(i, len(vals) - 1))])
    return sigma


def load_mu_pair_series(path=_MU_CSV):
    """(mu_pat(t), mu_mat(t)) callable from the CSV (step interpolation)."""
    df = pd.read_csv(path).sort_values("year")
    years = df["year"].to_numpy(dtype=float)
    pat = df["mu_pat"].to_numpy(dtype=float)
    mat = df["mu_mat"].to_numpy(dtype=float)

    def mu_pair(t):
        i = int(np.searchsorted(years, t, side="right") - 1)
        i = max(0, min(i, len(pat) - 1))
        return (float(pat[i]), float(mat[i]))
    return mu_pair


def load_sex_rates(base=BASE):
    """Convenience loader -> (sigma_callable, mu_pair_callable)."""
    return (load_sigma_series(os.path.join(base, "data", "immigrant_sex_share.csv")),
            load_mu_pair_series(os.path.join(base, "data", "sex_specific_mortality.csv")))
