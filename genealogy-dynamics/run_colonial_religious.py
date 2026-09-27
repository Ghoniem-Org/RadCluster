#!/usr/bin/env python3
"""Colonial run with the religious dimension: 1650 -> 2025, one (rule, alpha).

Usage: GENEALOGY_MAX_DEPTH=14 python3 run_colonial_religious.py <rule> <alpha> <outdir> [delta_max]

- Rates: data/inputs.csv (1650-2024; 2025 holds 2024 values); emigration
  eps calibrated (data/emigration_calibration.txt).
- IC 1650: 50,368 people, all Protestant, all C_0, split sigma/(1-sigma)
  pat/mat (ASSUMED spin-up IC).
- Immigration religious composition iota_g(t): era brackets
  (genealogy/religion.py; pre-1965 ASSUMED, 1965+ ESTIMATED from Pew 2013).
- Disaffiliation (d,sex,g)->(d,sex,una) at delta(t), ramp 1970-1990 to
  delta_max (default 0.007/yr, calibrated for ~29% unaffiliated by 2025).
- Writes results_religious_<rule>_alpha<a>.csv + snapshots + heatmaps into
  <outdir>.
"""
import os
import sys

# MAX_DEPTH must be fixed before genealogy modules are imported.
os.environ.setdefault("GENEALOGY_MAX_DEPTH", "14")

import numpy as np
import pandas as pd

from genealogy.populations import N_CLUSTERS, MAX_DEPTH, RULE_DISPLAY
from genealogy.religion import (
    G, GROUPS, GROUP_INDEX, GROUP_DISPLAY, N_STATE_R, rstate_index,
    RHO_VEC, immigrant_composition, disaffiliation_rate,
)
from genealogy.kernels import (
    fertility_schedule, SEX_RATIO_AT_BIRTH,
    mating_flux_religious,
)
from genealogy.model import rhs_religious, total_balance_religious
from genealogy.visualization import plot_pedigree_heatmap
from scipy.integrate import solve_ivp
from functools import partial

BASE = os.path.dirname(os.path.abspath(__file__))

from genealogy.sex_rates import load_sigma_series, load_mu_pair_series

S = SEX_RATIO_AT_BIRTH
SIGMA = load_sigma_series()  # step 5: sigma(t)
MU_PAIR = load_mu_pair_series()  # step 5: (mu_pat(t), mu_mat(t))

inputs = pd.read_csv(os.path.join(BASE, "data", "inputs.csv"))
inputs = inputs.sort_values("year").reset_index(drop=True)

calib = {}
with open(os.path.join(BASE, "data", "emigration_calibration.txt")) as f:
    for line in f:
        line = line.strip()
        if line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        try:
            calib[k.strip()] = float(v)
        except ValueError:
            pass
EPS = calib.get("eps_calibrated", 0.00075)

RULE = sys.argv[1]
ALPHA = float(sys.argv[2])
OUTDIR = sys.argv[3]
DELTA_MAX = float(sys.argv[4]) if len(sys.argv) > 4 else 0.008
os.makedirs(OUTDIR, exist_ok=True)


def rate_series_factory():
    years = inputs["year"].to_numpy()
    I = inputs["immigration_M"].to_numpy()
    beta = 2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0
    def rate_series(t):
        i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
        bvec = np.full(N_CLUSTERS, beta[i])
        return (I[i], bvec, MU_PAIR(t), EPS,
                immigrant_composition(t), disaffiliation_rate(t, DELTA_MAX))
    return rate_series


rate_series = rate_series_factory()

# --- IC 1650: all Protestant, C_0 --------------------------------------------
c0 = np.zeros(N_STATE_R)
c0[rstate_index(0, "pat", "pro")] = SIGMA(1650.0) * 0.050368
c0[rstate_index(0, "mat", "pro")] = (1.0 - SIGMA(1650.0)) * 0.050368


def rhs_td(t, c):
    I, beta, mu, eps, iot, dlt = rate_series(t)
    return rhs_religious(t, c, I=I, beta=beta, alpha=ALPHA, mu=mu, eps=eps,
                         sigma=SIGMA, s=S, rule=RULE, rho=RHO_VEC,
                         iota=iot, delta=dlt)


t_eval = np.arange(1650.0, 2025.0 + 1.0, 1.0)
sol = solve_ivp(rhs_td, (1650.0, 2025.0), c0, t_eval=t_eval,
                method="RK45", rtol=1e-8, atol=1e-10)
if not sol.success:
    raise RuntimeError(f"solve_ivp failed: {sol.message}")
if np.any(sol.y < -1e-6):
    raise RuntimeError("negative state detected")
t, Y = sol.t, sol.y.T
print(f"[{RULE} a={ALPHA}] integrated {len(t)} steps, total_2025={Y[-1].sum():.1f}M",
      flush=True)

# --- balance check ------------------------------------------------------------
dev = 0.0
for i in range(0, len(t), 50):
    I_i, beta_i, mu_i, eps_i, iot_i, dlt_i = rate_series(t[i])
    rsum = rhs_religious(t[i], Y[i], I=I_i, beta=beta_i, alpha=ALPHA,
                         mu=mu_i, eps=eps_i, sigma=SIGMA, s=S, rule=RULE,
                         rho=RHO_VEC, iota=iot_i, delta=dlt_i).sum()
    _, B_i = mating_flux_religious(Y[i], beta_i, ALPHA, RHO_VEC, RULE)
    Ysex = Y[i].reshape(G, 2, N_CLUSTERS)
    pred = total_balance_religious(I_i, B_i,
                                   (Ysex[:, 0, :].sum(), Ysex[:, 1, :].sum()),
                                   mu_i, eps_i)
    dev = max(dev, abs(rsum - pred))
print(f"[{RULE} a={ALPHA}] balance max dev = {dev:.2e}", flush=True)

# --- outputs -------------------------------------------------------------------
Cg = Y.reshape(len(t), G, 2, N_CLUSTERS)          # [t, g, sex, d]
total = Y.sum(axis=1)
group_tot = Cg.sum(axis=(2, 3))                    # [t, g]
group_share = group_tot / total[:, None]
Cn = Cg[:, :, :, MAX_DEPTH].sum(axis=(1, 2)) / total   # deepest-cluster total share
Cn_by_g = Cg[:, :, :, MAX_DEPTH].sum(axis=2) / group_tot  # [t, g]

SUM = pd.DataFrame({"year": t, "total_M": total, "Cn_share": Cn})
for gi, g in enumerate(GROUPS):
    SUM[f"group_{g}_M"] = group_tot[:, gi]
    SUM[f"group_{g}_share"] = group_share[:, gi]
    SUM[f"Cn_share_{g}"] = Cn_by_g[:, gi]
tag = f"religious_{RULE}_alpha{ALPHA}"
np.savez(os.path.join(OUTDIR, f"finalstate_{tag}.npz"), t=t, Y=Y)
SUM.to_csv(os.path.join(OUTDIR, f"results_{tag}.csv"), index=False)

snap_years = list(range(1650, 2026, 25))
idx = [int(np.argmin(np.abs(t - yy))) for yy in snap_years]
snap = pd.DataFrame({"year": t[idx]})
for gi, g in enumerate(GROUPS):
    snap[f"group_{g}_share"] = group_share[idx, gi]
snap.to_csv(os.path.join(OUTDIR, f"snapshots_{tag}.csv"), index=False)

# heatmaps per group: cluster-share distribution (n, snaps) per group
for gi, g in enumerate(GROUPS):
    Zg = (Cg[idx][:, gi].sum(axis=1) / total[idx][:, None]).T  # (n, snaps)
    out = os.path.join(OUTDIR, f"heatmap_{tag}_{g}.png")
    plot_pedigree_heatmap(
        np.array(snap_years), Zg, out,
        title=(f"Cluster distribution 1650-2025 [{RULE_DISPLAY[RULE]} rule] "
               f"-- {GROUP_DISPLAY[g]}"),
        subtitle=(f"eps={EPS}/yr, α={ALPHA}, s={S}, religious dimension"),
        vmin=0.0, vmax=0.40, population=g)
    print("wrote", out, flush=True)

row = dict(rule=RULE, alpha=ALPHA, total_2025_M=float(total[-1]),
           Cn_share_2025=float(Cn[-1]), balance_max_dev=float(dev),
           delta_max=DELTA_MAX)
for gi, g in enumerate(GROUPS):
    row[f"share2025_{g}"] = float(group_share[-1, gi])
    row[f"Cn2025_{g}"] = float(Cn_by_g[-1, gi])
pd.DataFrame([row]).to_csv(os.path.join(OUTDIR, f"summary_{tag}.csv"), index=False)
print("done", tag, flush=True)
