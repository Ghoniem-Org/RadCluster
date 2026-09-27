"""Time-varying pedigree-cluster assortativity alpha(t).

Loads data/assortativity_schedule.csv (built by
data/build_assortativity_schedule.py) and interpolates it to a callable
alpha(t) for use anywhere the model takes an `alpha` keyword:
`genealogy.model.rhs_td` accepts either a float (constant, the default
behavior, unchanged) or a callable alpha(t).

Interpolation is piecewise linear between anchor years; outside the table
range the edge values are held (flat), which is the documented neutral
continuation.
"""
import os

import numpy as np
import pandas as pd

_DEFAULT_CSV = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                            "data", "assortativity_schedule.csv")


class AssortativitySchedule:
    """Piecewise-linear alpha(t) from the schedule table."""

    def __init__(self, df):
        df = df.sort_values("year").reset_index(drop=True)
        self.years = df["year"].to_numpy(dtype=float)
        self.alphas = df["alpha"].to_numpy(dtype=float)
        self.table = df

    def __call__(self, t):
        a = np.interp(np.asarray(t, dtype=float), self.years, self.alphas)
        return a.item() if np.ndim(t) == 0 else a

    def at(self, t):
        return self(t)

    def __repr__(self):
        return (f"AssortativitySchedule({len(self.years)} anchors, "
                f"{self.years[0]:.0f}-{self.years[-1]:.0f}, "
                f"range {self.alphas.min():.3f}-{self.alphas.max():.3f})")


def load_schedule(path=None):
    """Load the CSV into a DataFrame (year, alpha, quality, basis)."""
    return pd.read_csv(_DEFAULT_CSV if path is None else path)


def alpha_schedule(path=None):
    """Return an AssortativitySchedule callable alpha(t)."""
    return AssortativitySchedule(load_schedule(path))


def constant_alpha(value):
    """A constant callable, useful for regression tests against float alpha."""
    return lambda t: float(value)
