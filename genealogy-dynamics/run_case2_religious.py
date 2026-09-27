#!/usr/bin/env python3
"""CASE 2: forward projection 2025 -> 2050 with the religious dimension.

- IC: the 2025 state vector of the Case-1 religious colonial runs
  (outputs/colonial_religious/finalstate_religious_<rule>_alpha<a>.npz),
  used AS-IS (the old 2225 driver rescaled to a 342M anchor; Case 2
  deliberately continues the colonial state with no rescaling). A continuity
  check asserts the handoff state matches the Case-1 saved 2025 aggregates
  to < 1e-9 relative.
- Inputs: PROJECTED 2025-2050 series (data/projection_inputs.py) --
  disaffiliation continues at the calibrated delta_max = 8e-3/yr.
- alpha = 0.8 (enforced central) / 0.0 (null) / 0.9 (bound); n = 14
  (GENEALOGY_MAX_DEPTH=14), 8 groups, 240 ODEs.

Writes to outputs/<stamp>_case2_religious/: CSVs, summary_case2_religious.csv,
case2_inputs.md, verification.md.

Usage: python3 run_case2_religious.py
"""
import os

# MAX_DEPTH must be fixed before genealogy modules are imported.
os.environ.setdefault("GENEALOGY_MAX_DEPTH", "14")

import sys
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from genealogy.populations import N_CLUSTERS, MAX_DEPTH, RULES, RULE_DISPLAY
from genealogy.religion import (
    G, GROUPS, N_STATE_R, RHO_VEC, GROUP_DISPLAY,
)
from genealogy.kernels import (
    mating_flux_religious, SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION,
)
from genealogy.model import rhs_religious, total_balance_religious
from data.projection_inputs import (
    projection_rate_series_religious, describe as describe_inputs,
    SIGMA_CASE2,
)

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_case2_religious")
os.makedirs(OUTDIR, exist_ok=True)

COLONIAL_REL = os.environ.get("GENEALOGY_CASE1_RELIGIOUS_DIR",
                              os.path.join(BASE, "outputs", "colonial_religious"))

S = SEX_RATIO_AT_BIRTH
SIGMA = SIGMA_CASE2
ALPHAS = [0.0, 0.8, 0.9]
T0, T1 = 2025.0, 2050.0


def handoff_state(rule, alpha):
    """2025 state vector from the Case-1 religious colonial final state.

    Continuity check: the loaded state must reproduce the Case-1 saved
    2025 totals and group shares to < 1e-9 relative.
    """
    tag = f"religious_{rule}_alpha{alpha}"
    src = np.load(os.path.join(COLONIAL_REL, f"finalstate_{tag}.npz"))
    t_c, Yc = src["t"], src["Y"]
    assert abs(t_c[-1] - 2025.0) < 1e-9
    c0 = Yc[-1]
    df = pd.read_csv(os.path.join(COLONIAL_REL, f"results_{tag}.csv"))
    last = df.iloc[-1]
    rel = abs(c0.sum() - last["total_M"]) / last["total_M"]
    assert rel < 1e-9, f"continuity total failed for {rule}/a={alpha}: {rel:.2e}"
    Cg = c0.reshape(G, 2, N_CLUSTERS)
    gs = Cg.sum(axis=(1, 2)) / c0.sum()
    worst = 0.0
    for gi, g in enumerate(GROUPS):
        want = float(last[f"group_{g}_share"])
        worst = max(worst, abs(gs[gi] - want) / want)
    assert worst < 1e-9, f"continuity groups failed for {rule}/a={alpha}: {worst:.2e}"
    return c0, max(rel, worst)


def rhs_td(t, c, alpha, rule):
    I, beta, mu, eps, iota, delta = projection_rate_series_religious(t)
    return rhs_religious(t, c, I=I, beta=beta, alpha=alpha, mu=mu, eps=eps,
                         sigma=SIGMA, s=S, rule=rule, rho=RHO_VEC,
                         iota=iota, delta=delta)


def run_case2(rule, alpha):
    c0, cont = handoff_state(rule, alpha)
    from functools import partial
    f = partial(rhs_td, alpha=alpha, rule=rule)
    t_eval = np.arange(T0, T1 + 1.0, 1.0)
    sol = solve_ivp(f, (T0, T1), c0, t_eval=t_eval,
                    method="RK45", rtol=1e-8, atol=1e-10)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    if np.any(sol.y < -1e-6):
        raise RuntimeError("negative state detected")
    t, Y = sol.t, sol.y.T

    # exact total-balance identity at sampled points
    dev = 0.0
    for i in range(0, len(t), 5):
        I_i, beta_i, mu_i, eps_i, iota_i, delta_i = \
            projection_rate_series_religious(t[i])
        rsum = rhs_religious(t[i], Y[i], I=I_i, beta=beta_i, alpha=alpha,
                             mu=mu_i, eps=eps_i, sigma=SIGMA, s=S, rule=rule,
                             rho=RHO_VEC, iota=iota_i, delta=delta_i).sum()
        _, B_i = mating_flux_religious(Y[i], beta_i, alpha, RHO_VEC, rule)
        Ysex = Y[i].reshape(G, 2, N_CLUSTERS)
        pred = total_balance_religious(
            I_i, B_i, (Ysex[:, 0, :].sum(), Ysex[:, 1, :].sum()), mu_i, eps_i)
        dev = max(dev, abs(rsum - pred))
    assert dev < 1e-8, f"balance failed for {rule}/a={alpha}: {dev:.2e}"

    Cg = Y.reshape(len(t), G, 2, N_CLUSTERS)
    total = Y.sum(axis=1)
    group_tot = Cg.sum(axis=(2, 3))
    group_share = group_tot / total[:, None]
    Cn = Cg[:, :, :, MAX_DEPTH].sum(axis=(1, 2)) / total
    Cn_by_g = Cg[:, :, :, MAX_DEPTH].sum(axis=2) / group_tot
    SUM = pd.DataFrame({"year": t, "total_M": total, "Cn_share": Cn})
    for gi, g in enumerate(GROUPS):
        SUM[f"group_{g}_M"] = group_tot[:, gi]
        SUM[f"group_{g}_share"] = group_share[:, gi]
        SUM[f"Cn_share_{g}"] = Cn_by_g[:, gi]
    SUM.to_csv(os.path.join(OUTDIR,
                            f"results_case2_religious_{rule}_alpha{alpha}.csv"),
               index=False)
    print(f"[{rule} a={alpha}] 2050 total={total[-1]:.1f}M "
          f"Cn={Cn[-1]*100:.2f}% una={group_share[-1, GROUPS.index('una')]*100:.1f}% "
          f"continuity={cont:.1e} balance={dev:.1e}", flush=True)
    row = dict(rule=rule, alpha=alpha, total_2050_M=float(total[-1]),
               Cn_share_2050=float(Cn[-1]), Cn_share_2025=float(Cn[0]),
               continuity_rel_err=float(cont), balance_max_dev=float(dev))
    for gi, g in enumerate(GROUPS):
        row[f"share2050_{g}"] = float(group_share[-1, gi])
        row[f"share2025_{g}"] = float(group_share[0, gi])
        row[f"Cn2050_{g}"] = float(Cn_by_g[-1, gi])
    return row


def main():
    rows = [run_case2(rule, alpha) for rule in RULES for alpha in ALPHAS]
    pd.DataFrame(rows).to_csv(
        os.path.join(OUTDIR, "summary_case2_religious.csv"), index=False)
    with open(os.path.join(OUTDIR, "case2_inputs.md"), "w") as f:
        f.write(describe_inputs())
        f.write("\nIC: 2025 state from the Case-1 religious colonial runs "
                "(outputs/colonial_religious/finalstate_*.npz), used as-is "
                "(no 342M-anchor rescaling).\n")
        f.write(f"alpha: {ALPHAS} (0.8 central/enforced, 0.0 null, 0.9 bound); "
                f"rho_g: Pew endogamy estimates (ESTIMATED) / assumed dhr, oth; "
                f"child takes mother's group (ASSUMED); sigma={SIGMA} (assumed); "
                f"s={S} (measured NCHS); n=14.\n")
    with open(os.path.join(OUTDIR, "verification.md"), "w") as f:
        f.write("# Case-2 (religious) verification\n\n")
        for r in rows:
            f.write(f"- {r['rule']} alpha={r['alpha']}: continuity rel err "
                    f"{r['continuity_rel_err']:.2e}, total-balance max dev "
                    f"{r['balance_max_dev']:.2e}, non-negativity enforced by "
                    "solver (raises if any state < -1e-6)\n")
    print("outdir:", OUTDIR)


if __name__ == "__main__":
    main()
