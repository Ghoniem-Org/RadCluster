#!/usr/bin/env python3
"""CASE 2 regional: forward projection 2025 -> 2050, regional two-sex model.

- IC: the 2025 state vector from the Case-1 regional colonial runs, obtained
  by re-running 1650-2025 in memory (identical code path, so continuity is
  exact by construction).
- Inputs: national projected 2025-2050 series from data/projection_inputs
  (projected_cbr/cdr, frozen immigration -- every projected input labeled
  assumed/projected, never measured); inter-regional migration held at the
  2024 ACS matrix (ASSUMED); immigration regional allocation held at the
  2023 DHS LPR pattern (ASSUMED).
- alpha = 0.0 (random-mating null) / 0.8 (enforced central endogamy) /
  0.9 (strong-endogamy bound); shallowest + deepest inheritance rules.

Writes to outputs/<stamp>_case2_regional/:
  results_case2_regional_<rule>_alpha<a>.csv, summary_case2_regional.csv.
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd

os.environ.setdefault("GENEALOGY_MAX_DEPTH", "14")

from genealogy.populations import N_CLUSTERS, MAX_DEPTH, RULE_DISPLAY
from genealogy.kernels import SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION
from genealogy.regional import (
    REGIONS, N_REGIONS, rstate_index, total_balance_regional, rhs_regional,
    cluster_totals_regional, region_totals,
)
from genealogy.solver import run_hindcast
from data.projection_inputs import projection_rate_series, SIGMA_CASE2

import run_colonial_regional as C1
from run_colonial_regional import run_regional

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_case2_regional")
os.makedirs(OUTDIR, exist_ok=True)

S = SEX_RATIO_AT_BIRTH
SIGMA = SIGMA_CASE2
ALPHAS = [0.0, 0.8, 0.9]
T0, T1 = 2025.0, 2050.0

# Held (ASSUMED) regional inputs for 2025-2050: lam_at/M_at already hold
# their last benchmark constant outside the measured window.
LAM_HOLD = C1.lam_at(2025.0)   # 2023 DHS LPR regional shares
M_HOLD = C1.M_at(2025.0)       # 2024 ACS inter-regional migration matrix


def projection_rate_series_regional(t):
    I, bpat, bmat, mu, eps = projection_rate_series(t)
    return I, LAM_HOLD, bpat, bmat, mu, eps, M_HOLD


def run_case2_regional(rule, alpha):
    # Case-1 re-run in memory -> exact 2025 handoff state.
    _, Y1 = run_regional(C1.c0, alpha, rule, dt=1.0, t0=1650.0, t1=2025.0,
                         rate_series_fn=C1.rate_series)
    c0 = Y1[-1]
    t, Y = run_regional(c0, alpha, rule, dt=0.5, t0=T0, t1=T1,
                        rate_series_fn=projection_rate_series_regional)
    # exact total-balance identity on the projection window
    from genealogy.regional import rhs_regional_td
    dev = 0.0
    for i in range(0, len(t), 5):
        I_i, lam_i, bpat, bmat, mu_i, eps_i, M_i = \
            projection_rate_series_regional(t[i])
        rvec, B_i = rhs_regional(
            t[i], Y[i], I=I_i, lam=lam_i, beta_pat=bpat, beta_mat=bmat,
            alpha=alpha, mu=mu_i, eps=eps_i, M=M_i, sigma=SIGMA, s=S,
            rule=rule)
        rsum = rvec.sum()
        rsum_td = rhs_regional_td(t[i], Y[i],
                                  rate_series=projection_rate_series_regional,
                                  alpha=alpha, sigma=SIGMA, s=S,
                                  rule=rule).sum()
        assert abs(rsum - rsum_td) < 1e-9
        pred = total_balance_regional(
            I_i, B_i,
            (Y[i].reshape(N_REGIONS, 2, N_CLUSTERS)[:, 0, :].sum(),
             Y[i].reshape(N_REGIONS, 2, N_CLUSTERS)[:, 1, :].sum()),
            mu_i, eps_i)
        dev = max(dev, abs(rsum - pred))
    assert dev < 1e-8, f"balance check failed: {dev:.2e}"

    T = cluster_totals_regional(Y)
    total = T.sum(axis=1)
    share = T / total[:, None]
    RT = region_totals(Y)  # (steps, 4); region order = REGIONS
    SUM = pd.DataFrame({"year": t, "total_M": total,
                        "Cn_share": share[:, MAX_DEPTH]})
    for ri, r in enumerate(REGIONS):
        rtot = RT[:, ri]
        SUM[f"total_{r}_M"] = rtot
        for d in range(N_CLUSTERS):
            col = (Y[:, rstate_index(d, "pat", r)]
                   + Y[:, rstate_index(d, "mat", r)]) / total
            SUM[f"share_C{d}_{r}"] = col
    SUM.to_csv(os.path.join(
        OUTDIR, f"results_case2_regional_{rule}_alpha{alpha}.csv"), index=False)

    reg = {}
    for ri, r in enumerate(REGIONS):
        rtot = RT[:, ri]
        c12 = (Y[:, rstate_index(MAX_DEPTH, "pat", r)]
               + Y[:, rstate_index(MAX_DEPTH, "mat", r)])
        reg[r] = (float(rtot[-1]), float(c12[-1] / rtot[-1]))
    print(f"[{rule} a={alpha}] 2050 total={total[-1]:.1f}M "
          f"Cn={share[-1, MAX_DEPTH]*100:.2f}% balance={dev:.1e} | " +
          " ".join(f"{r}={reg[r][1]*100:.2f}%" for r in REGIONS), flush=True)
    return dict(rule=rule, alpha=alpha, total_2050_M=float(total[-1]),
                Cn_share_2050=float(share[-1, MAX_DEPTH]),
                **{f"Cn_share_2050_{r}": reg[r][1] for r in REGIONS},
                **{f"total_2050_{r}_M": reg[r][0] for r in REGIONS})


def main():
    rows = []
    for rule in ("shallowest_inheritance", "deepest_inheritance"):
        for alpha in ALPHAS:
            rows.append(run_case2_regional(rule, alpha))
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "summary_case2_regional.csv"),
                              index=False)
    with open(os.path.join(OUTDIR, "case2_regional_inputs.md"), "w") as f:
        f.write(
            "# Regional Case-2 inputs (2025-2050)\n\n"
            "- National I(t), b(t), mu(t), eps(t): data/projection_inputs.py "
            "(projected CBR/CDR trends, immigration frozen at 2024: "
            "ASSUMED/projected).\n"
            "- Inter-regional migration matrix: held at 2024 ACS 1-year "
            "state-to-state flows aggregated to regions (ASSUMED).\n"
            "- Immigration regional allocation: held at 2023 DHS Yearbook "
            "Table 4 LPR state-of-residence pattern (ASSUMED; LPR pattern "
            "applied to all of I(t)).\n"
            "- Mating within region (region assortativity 1, ASSUMED); "
            "national vital rates applied uniformly (regional vital "
            "differences: future work).\n")
    print("outdir:", OUTDIR)


if __name__ == "__main__":
    main()
