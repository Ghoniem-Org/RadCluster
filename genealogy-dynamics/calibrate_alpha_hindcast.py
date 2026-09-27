#!/usr/bin/env python3
"""Hindcast alpha-grid cross-check for Step 3 (CPS calibration of alpha).

Runs the VALIDATED 1980->2024 hindcast (measured 1980 ICs, correct C0 stock)
over an alpha grid, both rules, and computes the 2013 two-parent share ratio
deepest-C1 / shallowest-C1 vs the CPS 2013 anchor (0.59, P23-214 Fig 3).

Caveat (documented): the 1980 C1 IC (20M, assumed) sits entirely in C1 in
both rules' runs; the true 1980 two-parent stock was likely < 20M and the
true 1980 second-gen stock likely > 20M, so this ratio is biased slightly
UP (toward 1) and its alpha slightly DOWN. See STEP_LOG Step 3.
"""
import os
import numpy as np
import pandas as pd

import run_hindcast as rh
from genealogy.populations import N_CLUSTERS, N_STATE
from genealogy.solver import run_hindcast, cluster_totals

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs", "20260926_212450_calibrate_alpha")

c0 = np.zeros(N_STATE)
ic = pd.read_csv(os.path.join(BASE, "data", "initial_condition_1980.csv"))
c0[:N_CLUSTERS] = ic["millions"].to_numpy() / 2.0
c0[N_CLUSTERS:] = ic["millions"].to_numpy() / 2.0
rate_series = rh.rate_series_factory()

ALPHAS = np.round(np.arange(0.55, 0.81, 0.025), 3)
rows = []
for a in ALPHAS:
    a = float(a)
    out = {}
    for rule in ("shallowest_inheritance", "deepest_inheritance"):
        t, Y = run_hindcast(c0, 1980.0, 2024.0, rate_series, alpha=a,
                            sigma=rh.SIGMA, s=rh.S, rule=rule, dt=0.25)
        T = cluster_totals(Y)
        total = T.sum(axis=1)
        share = T / total[:, None]
        i = int(np.argmin(np.abs(t - 2013.0)))
        out[rule] = (float(share[i, 0]), float(share[i, 1]), float(total[i]))
    sC1 = out["shallowest_inheritance"][1]
    lC1 = out["deepest_inheritance"][1]
    rows.append(dict(alpha=a,
                     C0_share_2013=out["shallowest_inheritance"][0],
                     strict_C1_2013=sC1, lenient_C1_2013=lC1,
                     ratio=lC1 / sC1))
    print(f"alpha={a:.3f} C0={out['shallowest_inheritance'][0]*100:.2f}% "
          f"strictC1={sC1*100:.2f}% lenientC1={lC1*100:.2f}% ratio={lC1/sC1:.3f}")
pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "fit_ratio_hindcast.csv"),
                          index=False)
print("wrote", os.path.join(OUTDIR, "fit_ratio_hindcast.csv"))
