#!/usr/bin/env python3
"""Calibrate the emigration rate eps by grid search against 1980-2024 totals.

Model totals are rule-independent (immigration, births, deaths, emigration do
not depend on the pedigree rule), so the shallowest inheritance rule is used.

Method: for each eps in a grid, run the 1980-2024 hindcast and compute the sum
of squared errors of model total population vs World Bank SP.POP.TOTL
(data/validation.csv). Report argmin(SSE), plus the endpoint-error-optimal
eps for context. Write the result to data/emigration_calibration.txt.

The calibrated eps is then hardcoded in run_hindcast.py and used for the
colonial residual in build_inputs.py (re-run it after this script).
"""
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(BASE, ".."))

from run_hindcast import load_inputs, load_ic, load_validation, ALPHA  # noqa: E402
from genealogy.solver import run_hindcast  # noqa: E402


def sse_for_eps(years, I, beta, mu, P0, actual, eps):
    _, Y = run_hindcast(
        P0, years, I_arr=I, beta_arr=beta, mu_arr=mu,
        alpha=ALPHA, eps=eps, rule="shallowest_inheritance")
    model = Y.sum(axis=0)
    return float(np.sum((model - actual) ** 2)), float(model[-1] - actual[-1])


def main():
    years, I, beta, mu = load_inputs()
    P0 = load_ic()
    val = load_validation()
    actual = val["total_pop_M"].to_numpy(float)
    assert len(actual) == len(years), (len(actual), len(years))

    grid = np.arange(0.0, 0.004001, 0.0001)
    rows = []
    print(f"{'eps':>8} {'SSE':>12} {'endpoint_err_M':>14}")
    for eps in grid:
        sse, end_err = sse_for_eps(years, I, beta, mu, P0, actual, float(eps))
        rows.append((float(eps), sse, end_err))
        print(f"{eps:8.4f} {sse:12.3f} {end_err:14.3f}")

    best_sse = min(rows, key=lambda r: r[1])
    best_end = min(rows, key=lambda r: abs(r[2]))
    print(f"\nargmin SSE:      eps={best_sse[0]:.4f}  SSE={best_sse[1]:.3f}  "
          f"endpoint err={best_sse[2]:+.3f}M")
    print(f"argmin |enderr|: eps={best_end[0]:.4f}  SSE={best_end[1]:.3f}  "
          f"endpoint err={best_end[2]:+.3f}M")

    # Refine around the SSE minimum at 5e-5 resolution.
    fine = np.arange(max(0.0, best_sse[0] - 0.0002),
                     best_sse[0] + 0.00021, 0.00005)
    fine_rows = [(float(e),) + sse_for_eps(years, I, beta, mu, P0, actual, float(e))[:2]
                 for e in fine]
    fine_best = min(fine_rows, key=lambda r: r[1])
    sse, end_err = sse_for_eps(years, I, beta, mu, P0, actual, fine_best[0])
    print(f"refined:         eps={fine_best[0]:.5f}  SSE={sse:.3f}  "
          f"endpoint err={end_err:+.3f}M")

    out = os.path.join(BASE, "emigration_calibration.txt")
    with open(out, "w") as f:
        f.write("# Emigration rate calibration (data/calibrate_emigration.py)\n")
        f.write("# Grid search on 1980-2024 total population vs World Bank "
                "SP.POP.TOTL.\n")
        f.write("# SSE = sum of squared errors over all 45 annual totals.\n")
        f.write(f"# Calibrated {__import__('datetime').date.today().isoformat()}\n")
        f.write(f"eps_calibrated={fine_best[0]:.5f}\n")
        f.write(f"sse_at_calibrated={sse:.3f}\n")
        f.write(f"endpoint_error_M={end_err:+.3f}\n")
        f.write(f"coarse_grid_argmin_sse={best_sse[0]:.4f}\n")
        f.write(f"coarse_grid_argmin_endpoint={best_end[0]:.4f}\n")
    print("wrote", out)
    print("\nNext: copy eps into run_hindcast.py (EPS) and re-run "
          "data/build_inputs.py so the colonial residual uses it.")


if __name__ == "__main__":
    main()
