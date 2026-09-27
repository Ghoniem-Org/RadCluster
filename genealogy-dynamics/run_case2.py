#!/usr/bin/env python3
"""CASE 2: forward projection 2025 -> 2050, two-sex model, both bloodline rules.

- IC: the 2025 state vector of the Case-1 colonial runs
  (outputs/20260927_161901_colonial), reconstructed from the saved
  per-cluster paternal/maternal shares x total population. A continuity
  check asserts the handoff state matches the Case-1 2025 totals to
  < 1e-9 relative.
- Inputs: PROJECTED 2025-2050 series (data/projection_inputs.py) -- every
  projected input is labeled assumed/projected, never measured.
- alpha = 0.8 (enforced central endogamy) / 0.0 (random-mating null) /
  0.9 (strong-endogamy bound).
- sigma = 0.5 (assumed); s = 0.512 (measured NCHS).

Writes to outputs/<stamp>_case2/: results_case2_<rule>_alpha<a>.csv,
summary_case2.csv, case2_inputs.md, verification.md.

The old 2025->2225 200-year projection (run_religious_scenario.py) and the
stylized 200-year scenarios (run.py) were deleted; this driver replaces
them with the 2025->2050 Case-2 study.
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd

from genealogy.populations import (
    N_CLUSTERS, MAX_DEPTH, N_STATE, RULES, RULE_DISPLAY, state_index,
)
from genealogy.kernels import (
    mating_flux, SEX_RATIO_AT_BIRTH, IMMIGRANT_MALE_FRACTION,
)
from genealogy.model import rhs_td, total_balance
from genealogy.solver import run_hindcast, cluster_totals, sex_totals
from data.projection_inputs import (
    projection_rate_series, describe as describe_inputs, SIGMA_CASE2,
)

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_case2")
os.makedirs(OUTDIR, exist_ok=True)

# Case-1 two-sex colonial run (pinned: this is the featured historical case).
# Case-1 two-sex colonial run (pinned; GENEALOGY_CASE1_DIR overrides).
COLONIAL_DIR = os.environ.get("GENEALOGY_CASE1_DIR",
                              os.path.join(BASE, "outputs", "20260927_161901_colonial"))
SUMMARY = pd.read_csv(os.path.join(COLONIAL_DIR, "summary_all.csv"))

S = SEX_RATIO_AT_BIRTH
SIGMA = SIGMA_CASE2
ALPHAS = [0.0, 0.8, 0.9]
T0, T1 = 2025.0, 2050.0


def handoff_state(rule, alpha):
    """2025 state vector from the Case-1 colonial outputs.

    Continuity check: the reconstructed total must match the Case-1
    summary 2025 total to < 1e-9 relative (the shares round-trip through
    float64 CSV columns, so this is a tight check).
    """
    df = pd.read_csv(
        os.path.join(COLONIAL_DIR, f"results_{rule}_alpha{alpha}.csv"))
    row = df[df["year"] == 2025.0].iloc[0]
    total = float(row["total_M"])
    c0 = np.zeros(N_STATE)
    for d in range(N_CLUSTERS):
        c0[state_index(d, "pat")] = float(row[f"share_C{d}_pat"]) * total
        c0[state_index(d, "mat")] = float(row[f"share_C{d}_mat"]) * total
    summ = SUMMARY[(SUMMARY["rule"] == rule)
                   & (SUMMARY["alpha"] == alpha)].iloc[0]
    rel = abs(c0.sum() - summ["total_2025_M"]) / summ["total_2025_M"]
    assert rel < 1e-9, f"continuity check failed for {rule}/a={alpha}: {rel:.2e}"
    return c0, rel


def run_case2(rule, alpha):
    c0, cont = handoff_state(rule, alpha)
    t, Y = run_hindcast(c0, T0, T1, projection_rate_series, alpha=alpha,
                        sigma=SIGMA, s=S, rule=rule, dt=0.25)
    # exact total-balance identity: sum(rhs) == I + B - (mu+eps)*total
    dev = 0.0
    for i in range(0, len(t), 10):
        I_i, bpat, bmat, mu_i, eps_i = projection_rate_series(t[i])
        rsum = rhs_td(t[i], Y[i], rate_series=projection_rate_series,
                      alpha=alpha, sigma=SIGMA, s=S, rule=rule).sum()
        _, B_i = mating_flux(Y[i], bpat, bmat, alpha, rule)
        pred = total_balance(I_i, B_i, Y[i, :N_CLUSTERS].sum(),
                             Y[i, N_CLUSTERS:].sum(), mu_i, eps_i)
        dev = max(dev, abs(rsum - pred))
    assert dev < 1e-8, f"balance check failed for {rule}/a={alpha}: {dev:.2e}"

    T = cluster_totals(Y)
    Tpat, Tmat = sex_totals(Y)
    total = T.sum(axis=1)
    share = T / total[:, None]
    SUM = pd.DataFrame({"year": t, "total_M": total,
                        "total_pat_M": Tpat, "total_mat_M": Tmat,
                        "Cn_share": share[:, MAX_DEPTH]})
    for d in range(N_CLUSTERS):
        SUM[f"share_C{d}_pat"] = Y[:, state_index(d, "pat")] / total
        SUM[f"share_C{d}_mat"] = Y[:, state_index(d, "mat")] / total
        SUM[f"share_C{d}"] = share[:, d]
    SUM.to_csv(os.path.join(OUTDIR,
                            f"results_case2_{rule}_alpha{alpha}.csv"),
               index=False)
    print(f"[{rule} a={alpha}] 2050 total={total[-1]:.1f}M "
          f"Cn={share[-1, MAX_DEPTH]*100:.2f}% continuity={cont:.1e} "
          f"balance={dev:.1e}", flush=True)
    return dict(rule=rule, alpha=alpha, total_2050_M=float(total[-1]),
                Cn_share_2050=float(share[-1, MAX_DEPTH]),
                Cn_share_2025=float(share[0, MAX_DEPTH]),
                continuity_rel_err=float(cont), balance_max_dev=float(dev))


def main():
    rows = []
    for rule in RULES:
        for alpha in ALPHAS:
            rows.append(run_case2(rule, alpha))
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "summary_case2.csv"),
                              index=False)
    with open(os.path.join(OUTDIR, "case2_inputs.md"), "w") as f:
        f.write(describe_inputs())
        f.write("\nIC: 2025 state from the Case-1 colonial runs "
                f"({os.path.basename(COLONIAL_DIR)}), reconstructed from saved "
                "per-cluster paternal/maternal shares x total.\n")
        f.write(f"alpha: {ALPHAS} (0.8 central/enforced, 0.0 null, 0.9 bound); "
                f"sigma={SIGMA} (assumed); s={S} (measured NCHS).\n")
    with open(os.path.join(OUTDIR, "verification.md"), "w") as f:
        f.write("# Case-2 (two-sex) verification\n\n")
        for r in rows:
            f.write(f"- {r['rule']} alpha={r['alpha']}: continuity rel err "
                    f"{r['continuity_rel_err']:.2e}, total-balance max dev "
                    f"{r['balance_max_dev']:.2e}, non-negativity enforced by "
                    "solver (raises if any state < -1e-6)\n")
    print("outdir:", OUTDIR)


if __name__ == "__main__":
    main()
