#!/usr/bin/env python3
"""Step 4: colonial Case 1 (1650->2025) and Case 2 (2025->2050) with
time-varying assortativity alpha(t), both bloodline rules.

- alpha(t): genealogy.assortativity.alpha_schedule(), loaded from
  data/assortativity_schedule.csv (built by
  data/build_assortativity_schedule.py). Era anchors documented in the
  builder; the schedule is OPT-IN (constant-alpha drivers are untouched).
- Regression test (exit nonzero if it fails): alpha = constant callable
  0.8 must reproduce the float-alpha=0.8 run exactly (max |diff| < 1e-9).
- Case 2: hands off from the alpha(t) Case-1 2025 state in THIS output
  dir; projection rates from data/projection_inputs.py; alpha(t) is flat
  at alpha(2025)=0.653 over 2025-2050 (schedule's neutral continuation).
- Baseline for comparison: the constant-alpha Case-1 run in
  outputs/20260927_161901_colonial/summary_all.csv (alpha=0.8 central).

Writes to outputs/<stamp>_colonial_alpha_t/:
  results_alpha_t_<rule>.csv   per-cluster shares + alpha_t column (Case 1)
  summary_alpha_t.csv          2025 totals, C_n shares, vs constant-0.8 baseline
  results_case2_alpha_t_<rule>.csv   Case-2 trajectories
  summary_case2_alpha_t.csv
  verification.md, colonial_alpha_t.md
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd

from genealogy.populations import (
    N_CLUSTERS, MAX_DEPTH, N_STATE, RULES, RULE_DISPLAY, state_index,
)
from genealogy.kernels import (
    SEX_RATIO_AT_BIRTH, mating_flux,
)
from genealogy.assortativity import alpha_schedule, constant_alpha
from genealogy.model import rhs_td, total_balance
from genealogy.solver import run_hindcast, cluster_totals, sex_totals
from data.projection_inputs import projection_rate_series

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_colonial_alpha_t")
os.makedirs(OUTDIR, exist_ok=True)

# Pinned constant-alpha baseline (the featured Case-1 run, alpha=0.8 central).
BASELINE = pd.read_csv(os.path.join(
    BASE, "outputs", "20260927_161901_colonial", "summary_all.csv"))

S = SEX_RATIO_AT_BIRTH
from genealogy.sex_rates import load_sigma_series, load_mu_pair_series

SIGMA = load_sigma_series()  # step 5: sigma(t)
MU_PAIR = load_mu_pair_series()

# --- inputs (identical to run_colonial.py) ------------------------------------
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


def rate_series_factory():
    years = inputs["year"].to_numpy()
    I = inputs["immigration_M"].to_numpy()
    beta = (2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0)
    def rate_series(t):
        i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
        bvec = np.full(N_CLUSTERS, beta[i])
        return I[i], bvec, bvec, MU_PAIR(t), EPS
    return rate_series


rate_series = rate_series_factory()
ALPHA_T = alpha_schedule()

c0 = np.zeros(N_STATE)
c0[state_index(0, "pat")] = SIGMA(1650.0) * 0.050368
c0[state_index(0, "mat")] = (1.0 - SIGMA(1650.0)) * 0.050368


def balance_check(t, Y, rs, rule, alpha_of_t, step):
    """Max |sum(rhs) - (I + B - (mu+eps)*total)| over sampled points."""
    dev = 0.0
    for i in range(0, len(t), step):
        a = alpha_of_t(t[i]) if callable(alpha_of_t) else alpha_of_t
        I_i, bpat, bmat, mu_i, eps_i = rs(t[i])
        rsum = rhs_td(t[i], Y[i], rate_series=rs, alpha=a,
                      sigma=SIGMA, s=S, rule=rule).sum()
        _, B_i = mating_flux(Y[i], bpat, bmat, a, rule)
        pred = total_balance(I_i, B_i, Y[i, :N_CLUSTERS].sum(),
                             Y[i, N_CLUSTERS:].sum(), mu_i, eps_i)
        dev = max(dev, abs(rsum - pred))
    return dev


def run_case1(rule, alpha, dt=1.0):
    t, Y = run_hindcast(c0, 1650.0, 2025.0, rate_series, alpha=alpha,
                        sigma=SIGMA, s=S, rule=rule, dt=dt)
    return t, Y


def summarize(t, Y):
    T = cluster_totals(Y)
    Tpat, Tmat = sex_totals(Y)
    total = T.sum(axis=1)
    share = T / total[:, None]
    return T, Tpat, Tmat, total, share


def write_case1(rule, alpha, tag, dt=1.0):
    t, Y = run_case1(rule, alpha, dt=dt)
    T, Tpat, Tmat, total, share = summarize(t, Y)
    alpha_col = (np.array([alpha(tt) for tt in t]) if callable(alpha)
                 else np.full_like(t, float(alpha)))
    SUM = pd.DataFrame({"year": t, "total_M": total,
                        "total_pat_M": Tpat, "total_mat_M": Tmat,
                        "alpha_t": alpha_col})
    for d in range(N_CLUSTERS):
        SUM[f"share_C{d}_pat"] = Y[:, state_index(d, "pat")] / total
        SUM[f"share_C{d}_mat"] = Y[:, state_index(d, "mat")] / total
        SUM[f"share_C{d}"] = share[:, d]
    SUM.to_csv(os.path.join(OUTDIR, f"results_{tag}_{rule}.csv"), index=False)
    dev = balance_check(t, Y, rate_series, rule, alpha, step=25)
    return t, Y, total, share, dev


def main():
    ver = []

    # --- regression: constant callable 0.8 vs float 0.8 -----------------------
    tA, YA = run_case1("shallowest_inheritance", 0.8)
    tB, YB = run_case1("shallowest_inheritance", constant_alpha(0.8))
    reg_diff = float(np.abs(YA - YB).max())
    reg_ok = reg_diff < 1e-9
    ver.append(f"regression (shallowest, callable-0.8 vs float-0.8): "
               f"max|diff|={reg_diff:.3e} -> {'PASS' if reg_ok else 'FAIL'}")
    print(ver[-1], flush=True)
    assert reg_ok, "regression FAILED"

    # --- Case 1 with alpha(t) ---------------------------------------------------
    rows = []
    for rule in RULES:
        t, Y, total, share, dev = write_case1(rule, ALPHA_T, "alpha_t")
        base = BASELINE[(BASELINE["rule"] == rule)
                        & (BASELINE["alpha"] == 0.8)].iloc[0]
        rows.append(dict(
            rule=rule,
            total_2025_M=float(total[-1]),
            Cn_share_2025=float(share[-1, MAX_DEPTH]),
            Cn_share_2025_baseline_alpha08=float(base["Cn_share_2025"]),
            Cn_share_delta_pp=(float(share[-1, MAX_DEPTH])
                               - float(base["Cn_share_2025"])) * 100,
            alpha_2025=float(ALPHA_T(2025.0)),
            balance_max_dev=float(dev)))
        print(f"[alpha(t) {rule}] 2025 total={total[-1]:.1f}M "
              f"Cn={share[-1, MAX_DEPTH]*100:.3f}% "
              f"(baseline a=0.8: {base['Cn_share_2025']*100:.3f}%) "
              f"balance={dev:.1e}", flush=True)
        ver.append(f"Case1 alpha(t) {rule}: balance max dev = {dev:.3e} "
                   f"(assert < 1e-8: {'PASS' if dev < 1e-8 else 'FAIL'})")
        assert dev < 1e-8
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "summary_alpha_t.csv"),
                              index=False)

    # --- Case 2 with alpha(t): handoff from the alpha(t) Case-1 2025 state ----
    c2rows = []
    for rule in RULES:
        df = pd.read_csv(os.path.join(OUTDIR, f"results_alpha_t_{rule}.csv"))
        row = df[df["year"] == 2025.0].iloc[0]
        total25 = float(row["total_M"])
        c0b = np.zeros(N_STATE)
        for d in range(N_CLUSTERS):
            c0b[state_index(d, "pat")] = float(row[f"share_C{d}_pat"]) * total25
            c0b[state_index(d, "mat")] = float(row[f"share_C{d}_mat"]) * total25
        rel = abs(c0b.sum() - total25) / total25
        assert rel < 1e-12, f"handoff reconstruction failed: {rel:.2e}"
        t2, Y2 = run_hindcast(c0b, 2025.0, 2050.0, projection_rate_series,
                              alpha=ALPHA_T, sigma=SIGMA, s=S, rule=rule,
                              dt=0.25)
        T, Tpat, Tmat, total, share = summarize(t2, Y2)
        dev = balance_check(t2, Y2, projection_rate_series, rule, ALPHA_T,
                            step=10)
        assert dev < 1e-8, f"Case-2 balance failed: {dev:.2e}"
        SUM = pd.DataFrame({"year": t2, "total_M": total,
                            "total_pat_M": Tpat, "total_mat_M": Tmat,
                            "alpha_t": [ALPHA_T(tt) for tt in t2],
                            "Cn_share": share[:, MAX_DEPTH]})
        for d in range(N_CLUSTERS):
            SUM[f"share_C{d}_pat"] = Y2[:, state_index(d, "pat")] / total
            SUM[f"share_C{d}_mat"] = Y2[:, state_index(d, "mat")] / total
            SUM[f"share_C{d}"] = share[:, d]
        SUM.to_csv(os.path.join(OUTDIR,
                                f"results_case2_alpha_t_{rule}.csv"),
                   index=False)
        c2rows.append(dict(rule=rule, total_2050_M=float(total[-1]),
                           Cn_share_2050=float(share[-1, MAX_DEPTH]),
                           Cn_share_2025=float(share[0, MAX_DEPTH]),
                           alpha_2025_2050=float(ALPHA_T(2025.0)),
                           continuity_rel_err=float(rel),
                           balance_max_dev=float(dev)))
        print(f"[case2 alpha(t) {rule}] 2050 total={total[-1]:.1f}M "
              f"Cn={share[-1, MAX_DEPTH]*100:.3f}% "
              f"(2025: {share[0, MAX_DEPTH]*100:.3f}%) balance={dev:.1e}",
              flush=True)
        ver.append(f"Case2 alpha(t) {rule}: continuity {rel:.1e}, "
                   f"balance {dev:.1e} -> PASS")
    pd.DataFrame(c2rows).to_csv(os.path.join(OUTDIR, "summary_case2_alpha_t.csv"),
                                index=False)

    with open(os.path.join(OUTDIR, "verification.md"), "w") as f:
        f.write("# Verification -- Step 4 alpha(t) runs\n\n")
        f.write("Solver enforces non-negativity (raises if any state < -1e-6); "
                "no violations in any run.\n\n")
        f.write("\n".join(f"- {line}" for line in ver) + "\n")

    with open(os.path.join(OUTDIR, "colonial_alpha_t.md"), "w") as f:
        f.write("# Colonial 1650-2025 + Case-2 2025-2050 with alpha(t)\n\n")
        f.write("alpha(t) from data/assortativity_schedule.csv "
                "(genealogy.assortativity.alpha_schedule); era anchors and "
                "quality labels documented in "
                "data/build_assortativity_schedule.py.\n\n")
        f.write("results_alpha_t_<rule>.csv: Case-1 per-cluster shares + "
                "alpha_t column (the schedule value at each year).\n")
        f.write("results_case2_alpha_t_<rule>.csv: Case-2 trajectories; "
                "alpha(t) held at alpha(2025)=0.653 over 2025-2050.\n")
        f.write("summary_alpha_t.csv: 2025 C_n shares vs the constant "
                "alpha=0.8 baseline (outputs/20260927_161901_colonial).\n")
    print("outdir:", OUTDIR)


if __name__ == "__main__":
    main()
