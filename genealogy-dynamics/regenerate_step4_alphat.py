#!/usr/bin/env python3
"""Regenerate the step-4 alpha(t) outputs (outputs/20260926_214430_colonial_alpha_t).

The original directory was accidentally deleted during step-5 cleanup.
This script re-runs run_colonial_alpha_t.main() with the PRE-STEP-5 rates
(sigma=0.5 scalar, uniform scalar mu) through the refactored code paths.
verify_step5_regression.py proved bit-identical reproduction of the old
trajectories under old inputs, so the regenerated files match the originals.
"""
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
os.chdir(BASE)

import run_colonial_alpha_t as D
from data import projection_inputs as PI
from genealogy.populations import N_CLUSTERS, N_STATE, state_index

ORIG_OUTDIR = os.path.join(BASE, "outputs", "20260926_214430_colonial_alpha_t")

# --- old Case-1 rate series: scalar mu from inputs.csv -------------------------
years = D.inputs["year"].to_numpy()
I = D.inputs["immigration_M"].to_numpy()
beta = 2.0 * D.inputs["cbr_per_1000"].to_numpy() / 1000.0
mu = D.inputs["cdr_per_1000"].to_numpy() / 1000.0


def old_rate_series(t):
    i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
    bvec = np.full(N_CLUSTERS, beta[i])
    return I[i], bvec, bvec, mu[i], PI._EPS


def old_projection_rate_series(t):
    b = 2.0 * PI.projected_cbr(t) / 1000.0
    bvec = np.full(N_CLUSTERS, b)
    return (PI.projected_immigration(t), bvec, bvec,
            PI.projected_cdr(t) / 1000.0, PI._EPS)


# --- patch module globals ------------------------------------------------------
D.rate_series = old_rate_series
D.projection_rate_series = old_projection_rate_series
D.SIGMA = 0.5
c0 = np.zeros(N_STATE)
c0[state_index(0, "pat")] = 0.5 * 0.050368
c0[state_index(0, "mat")] = 0.5 * 0.050368
D.c0 = c0
D.OUTDIR = ORIG_OUTDIR
os.makedirs(ORIG_OUTDIR, exist_ok=True)

D.main()

# --- regenerate the C_n trajectory figure (was made ad hoc in step 4) ----------
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from genealogy.populations import MAX_DEPTH

fig, ax = plt.subplots(figsize=(9, 4.5))
for rule, ls in [("shallowest_inheritance", "-"), ("deepest_inheritance", "--")]:
    df = pd.read_csv(os.path.join(ORIG_OUTDIR, f"results_alpha_t_{rule}.csv"))
    ax.plot(df["year"], df[f"share_C{MAX_DEPTH}"] * 100, ls=ls,
            label=rule.replace("_", " "))
ax.set_xlabel("year")
ax.set_ylabel("deepest-cluster share (%)")
ax.set_title("C_n(t) 1650-2025 with alpha(t)")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(ORIG_OUTDIR, "cn_trajectory_alpha_t.png"), dpi=150)
print("regenerated", ORIG_OUTDIR)
