#!/usr/bin/env python3
"""Figure for Step 3: CPS parental-nativity calibration of alpha.

drive_doc/figs/cps_calibration.png — two panels:
  left:  CPS second-generation share 2005-2024 (measured, MOE whiskers) vs
         model shallowest-C1 share at alpha = 0.0 / 0.66 (calibrated) / 0.8 / 0.9.
  right: 2013 two-foreign-born-parent share among the second generation vs
         alpha: model ratio deepest-C1/shallowest-C1 (colonial Case 1) with the
         CPS 2013 anchor (0.59, P23-214 Fig 3); verticals at calibrated 0.66
         and central 0.8.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
CAL = os.path.join(BASE, "outputs", "20260927_162519_calibrate_alpha")
FIG = os.path.join(BASE, "drive_doc", "figs", "cps_calibration.png")
os.makedirs(os.path.dirname(FIG), exist_ok=True)

cps = pd.read_csv(os.path.join(BASE, "data", "cps_secondgen_series.csv"))
grid = pd.read_csv(os.path.join(CAL, "grid_C1_series.csv"))
ratio = pd.read_csv(os.path.join(CAL, "fit_ratio_fine.csv"))

# CPS MOEs (percentage points) from Table 1: 0.2 (0.1 in 2007)
moe = np.where(cps.year == 2007, 0.1, 0.2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# --- left: second-gen share 2005-2024 ---------------------------------------
ax1.errorbar(cps.year, cps.second_gen_share_pct, yerr=moe, fmt="ko",
             ms=4, capsize=3, elinewidth=1.2,
             label="CPS ASEC, 2nd-gen share (measured ± MOE)")
for a, ls, col, lab in ((0.00, "--", "0.6", "model shallowest C1, α=0.0 (null)"),
                        (0.66, "-", "C3", "model shallowest C1, α=0.66 (CPS-calibrated)"),
                        (0.80, "-", "C0", "model shallowest C1, α=0.8 (central)"),
                        (0.90, ":", "C1", "model shallowest C1, α=0.9 (bound)")):
    coln = f"shallowest_C1_a{a:.2f}"
    yy = cps.year.to_numpy()
    if coln not in grid.columns:  # interpolate between grid points
        lo, hi = f"shallowest_C1_a{int(a*20)/20:.2f}", f"shallowest_C1_a{(int(a*20)+1)/20:.2f}"
        w = a * 20 - int(a * 20)
        vals = np.interp(yy, grid.year.to_numpy(),
                         ((1 - w) * grid[lo] + w * grid[hi]).to_numpy()) * 100.0
    else:
        vals = np.interp(yy, grid.year.to_numpy(), grid[coln].to_numpy()) * 100.0
    ax1.plot(yy, vals, ls=ls, color=col, lw=1.8, label=lab)
ax1.set_xlim(2004.5, 2024.5)
ax1.set_xlabel("year")
ax1.set_ylabel("second-generation share of population (%)")
ax1.set_title("Second-generation share: CPS vs colonial model (shallowest C1)\n"
              "CPS = US-born with ≥1 foreign-born parent (measured)")
ax1.legend(fontsize=8, loc="upper left")
ax1.grid(alpha=0.3)

# --- right: 2013 two-parent share vs alpha -----------------------------------
ax2.plot(ratio.alpha, ratio.ratio * 100.0, "C0-", lw=2,
         label="model deepest-C1 / shallowest-C1, 2013 (colonial)")
ax2.axhline(59.0, color="k", ls="--", lw=1.5,
            label="CPS 2013: 59% of 2nd gen has two FB parents\n(P23-214, Fig 3; measured)")
ax2.axvline(0.66, color="C3", ls="-", lw=2, label="calibrated α=0.66")
ax2.axvline(0.80, color="C0", ls=":", lw=2, label="central α=0.8 (kept)")
ax2.set_xlim(0.54, 0.81)
ax2.set_ylim(42, 72)
ax2.set_xlabel("class assortativity α")
ax2.set_ylabel("two-foreign-born-parent share of 2nd generation (%)")
ax2.set_title("2013 parental-nativity split vs α\n"
              "deepest C1 = US-born children of two immigrants")
ax2.legend(fontsize=8, loc="upper left")
ax2.grid(alpha=0.3)

fig.suptitle("Step 3 — CPS parental-nativity calibration of α  "
             "(two-sex colonial Case 1, deepest rule)",
             fontsize=12, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG, dpi=150)
print("wrote", FIG)

# --- CSV of plotted series ----------------------------------------------------
out = pd.DataFrame({
    "year": cps.year,
    "cps_second_gen_share_pct": cps.second_gen_share_pct,
    "cps_moe_pp": moe,
    "model_shallowest_C1_alpha0.66_pct": np.interp(
        cps.year, grid.year,
        grid["shallowest_C1_a0.65"] * 0.8 + grid["shallowest_C1_a0.70"] * 0.2) * 100.0,
    "model_shallowest_C1_alpha0.80_pct": np.interp(
        cps.year, grid.year, grid["shallowest_C1_a0.80"]) * 100.0,
})
out.to_csv(os.path.join(CAL, "cps_vs_model_series.csv"), index=False)
r = ratio[["alpha", "shallowest_C1_2013", "deepest_C1_2013", "ratio"]].copy()
r.to_csv(os.path.join(CAL, "ratio_vs_alpha_2013.csv"), index=False)
print("wrote CSVs")
