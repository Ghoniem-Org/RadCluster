#!/usr/bin/env python3
"""Step 3: CPS parental-nativity calibration of the class-assortativity alpha.

Targets (all MEASURED, CPS ASEC published estimates):
  T1: second-generation share of total population, 2005-2024, annual
      (Census 2024 CPS table package, Table 1 "Population by Generation:
      2005 to 2024", released Sept 2025; data/data/raw/...xlsx,
      data/cps_secondgen_series.csv). Model analog: SHALLOWEST C1 share
      (US-born with >=1 foreign-born parent; child(q,0)=C1 for all q).
  T2: two-foreign-born-parent share of total population in 2013
      = 0.117 (second-gen share, Table 1) x 0.59 (two-FB-parent share of
      the second generation, Census report P23-214 "Characteristics of the
      U.S. Population by Generational Status: 2013", Trevelyan 2016,
      Figure 3: 59% two FB parents / 20% FB mother only / 21% FB father
      only) = 6.90%. Model analog: DEEPEST C1 share (US-born children of
      two immigrants; child(0,0)=C1 only).

Method: two-sex colonial model (Case 1, 1650->2025), both rules, alpha grid
0.00..1.00 step 0.05 (21 points). SSE(alpha) computed separately on T1
(shallowest C1 vs annual series -> alpha_A) and T2 (deepest C1 vs 2013 anchor
-> alpha_B); joint minimizer reported as the calibrated alpha.

Outputs: outputs/<stamp>_calibrate_alpha/ with per-alpha CSV, fit summary,
and data for the figure (figure itself made by plot_cps_calibration.py).
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd

import run_colonial as rc
from genealogy.populations import N_CLUSTERS, state_index  # noqa: F401

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_calibrate_alpha")
os.makedirs(OUTDIR, exist_ok=True)

GRID = np.round(np.arange(0.0, 1.001, 0.05), 2)

# --- CPS anchors (measured) ---------------------------------------------------
cps = pd.read_csv(os.path.join(BASE, "data", "cps_secondgen_series.csv"))
cps_years = cps["year"].to_numpy()
cps_2nd = cps["second_gen_share"].to_numpy()
# 2013 two-FB-parent share of second generation (P23-214, Fig 3), measured:
SPLIT_2013 = 0.59
CPS_2ND_2013 = cps.loc[cps.year == 2013, "second_gen_share"].iloc[0]
CPS_TWOFB_2013 = CPS_2ND_2013 * SPLIT_2013   # 6.90% of total population


def series_for(rule, alpha):
    t, Y, T, Tpat, Tmat, total, share = rc.run_rule(rule, alpha)
    return t, share


def main():
    rows = []
    shallow_cache, deep_cache = {}, {}
    for alpha in GRID:
        t_s, sh_s = series_for("shallowest_inheritance", float(alpha))
        t_l, sh_l = series_for("deepest_inheritance", float(alpha))
        shallow_cache[float(alpha)] = (t_s, sh_s)
        deep_cache[float(alpha)] = (t_l, sh_l)
        # T1: shallowest C1 vs annual CPS second-gen series 2005-2024
        idx = np.array([int(np.argmin(np.abs(t_s - y))) for y in cps_years])
        resid = sh_s[idx, 1] - cps_2nd
        sse1 = float((resid ** 2).sum())
        rmse1 = float(np.sqrt((resid ** 2).mean()))
        # T2: deepest C1 (2013) vs CPS two-FB-parent absolute share
        i2013 = int(np.argmin(np.abs(t_l - 2013.0)))
        c1_2013 = float(sh_l[i2013, 1])
        sse2 = (c1_2013 - CPS_TWOFB_2013) ** 2
        rows.append(dict(alpha=float(alpha), sse_shallow=sse1, rmse_shallow=rmse1,
                         deep_C1_2013=c1_2013, sse_deep_2013=sse2))
        print(f"alpha={alpha:4.2f}  shallowest RMSE(05-24)={rmse1*100:5.3f}pp  "
              f"deepest C1(2013)={c1_2013*100:5.2f}% (CPS {CPS_TWOFB_2013*100:.2f}%)")
    fit = pd.DataFrame(rows)
    fit.to_csv(os.path.join(OUTDIR, "fit_grid.csv"), index=False)

    aA = float(fit.loc[fit.sse_shallow.idxmin(), "alpha"])
    aB = float(fit.loc[fit.sse_deep_2013.idxmin(), "alpha"])
    # joint: normalized SSEs (each series scaled to unit weight)
    s1 = fit.sse_shallow / fit.sse_shallow.max()
    s2 = fit.sse_deep_2013 / fit.sse_deep_2013.max()
    joint = s1 + s2
    aJ = float(fit.loc[joint.idxmin(), "alpha"])
    rA = fit.loc[fit.alpha == aA].iloc[0]
    rB = fit.loc[fit.alpha == aB].iloc[0]
    rJ = fit.loc[fit.alpha == aJ].iloc[0]

    # series at joint alpha for the figure
    t_s, sh_s = shallow_cache[aJ]
    t_l, sh_l = deep_cache[aJ]
    fig = pd.DataFrame({"year": t_s,
                        "shallowest_C1": sh_s[:, 1],
                        "deepest_C1": sh_l[:, 1]})
    fig.to_csv(os.path.join(OUTDIR, "modeled_series_joint.csv"), index=False)
    # also the full grid of C1 series (for sensitivity bands in figure)
    g = pd.DataFrame({"year": t_s})
    for a in GRID:
        g[f"shallowest_C1_a{a:.2f}"] = shallow_cache[float(a)][1][:, 1]
        g[f"deepest_C1_a{a:.2f}"] = deep_cache[float(a)][1][:, 1]
    g.to_csv(os.path.join(OUTDIR, "grid_C1_series.csv"), index=False)

    summary = {
        "CPS_twoFB_2013": CPS_TWOFB_2013,
        "alpha_A_shallowest_series": aA,
        "rmse_A_pp": rA["rmse_shallow"] * 100,
        "alpha_B_deepest_2013": aB,
        "deep_C1_2013_at_B_pp": rB["deep_C1_2013"] * 100,
        "alpha_joint": aJ,
        "rmse_joint_pp": rJ["rmse_shallow"] * 100,
        "deep_C1_2013_at_joint_pp": rJ["deep_C1_2013"] * 100,
    }
    pd.DataFrame([summary]).to_csv(os.path.join(OUTDIR, "fit_summary.csv"),
                                   index=False)
    print("\n=== calibration summary ===")
    print(f"CPS 2013 two-FB-parent share (measured): {CPS_TWOFB_2013*100:.2f}%")
    print(f"alpha_A (shallowest C1 vs 2005-2024 series): {aA:.2f}, RMSE {rA['rmse_shallow']*100:.3f} pp")
    print(f"alpha_B (deepest C1 vs 2013 anchor):     {aB:.2f}, "
          f"model C1(2013)={rB['deep_C1_2013']*100:.2f}%")
    print(f"alpha_joint:                            {aJ:.2f}")
    print("outdir:", OUTDIR)


if __name__ == "__main__":
    main()
