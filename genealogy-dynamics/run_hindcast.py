#!/usr/bin/env python3
"""Hindcast driver: 1980 -> 2024 with real historical rates, both bloodline rules.

Reads data/inputs.csv (immigration, birth/death rates), data/validation.csv
(actual totals, migrant-stock anchors, decennial foreign-born), and
data/initial_condition_1980.csv (cluster populations, split evenly M/F).

Writes outputs/<YYYYMMDD_HHMMSS>_hindcast_<rule>/ per rule:
  results.csv, stacked_share.png, provenance.md
plus outputs/<YYYYMMDD_HHMMSS>_hindcast_validation/:
  hindcast_validation.png (rule-independent totals + C_0 vs migrant stock)
and prints a validation table.
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from genealogy.populations import (
    N_CLUSTERS, MAX_DEPTH, N_STATE, RULES, RULE_DISPLAY,
)
from genealogy.kernels import (
    mating_flux, SEX_RATIO_AT_BIRTH,
)
from genealogy.model import total_balance
from genealogy.solver import run_hindcast, cluster_totals, sex_totals
from genealogy.visualization import plot_stacked_area

BASE = os.path.dirname(os.path.abspath(__file__))
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUT = os.path.join(BASE, "outputs")

S = SEX_RATIO_AT_BIRTH
from genealogy.sex_rates import load_sigma_series, load_mu_pair_series

SIGMA = load_sigma_series()  # step 5: sigma(t)
MU_PAIR = load_mu_pair_series()
ALPHA = 0.6  # central working assumption

# --- inputs ------------------------------------------------------------------
inputs = pd.read_csv(os.path.join(BASE, "data", "inputs.csv"))
inputs = inputs[(inputs["year"] >= 1980) & (inputs["year"] <= 2024)].reset_index(drop=True)
valid = pd.read_csv(os.path.join(BASE, "data", "validation.csv"))
ic = pd.read_csv(os.path.join(BASE, "data", "initial_condition_1980.csv"))
if len(ic) != N_CLUSTERS:
    raise RuntimeError(f"initial_condition_1980.csv has {len(ic)} rows, need {N_CLUSTERS}")

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

# 1980 IC: cluster totals from the CSV, split evenly between the sexes.
c0 = np.zeros(N_STATE)
c0[:N_CLUSTERS] = ic["millions"].to_numpy() / 2.0   # paternal block
c0[N_CLUSTERS:] = ic["millions"].to_numpy() / 2.0   # maternal block


def rate_series_factory():
    years = inputs["year"].to_numpy()
    I = inputs["immigration_M"].to_numpy()
    beta = (2.0 * inputs["cbr_per_1000"].to_numpy() / 1000.0)
    mu = inputs["cdr_per_1000"].to_numpy() / 1000.0

    def rate_series(t):
        i = int(np.clip(np.searchsorted(years, t, side="right") - 1, 0, len(years) - 1))
        bvec = np.full(N_CLUSTERS, beta[i])
        return I[i], bvec, bvec, MU_PAIR(t), EPS
    return rate_series


rate_series = rate_series_factory()


def check_balance_td(t, Y, alpha, rule, label):
    """Exact identity sum(rhs) == I + B - (mu+eps)*total_pat/mat at samples."""
    from genealogy.model import rhs_td
    worst = 0.0
    for k in range(0, len(t), 40):
        I_i, bpat, bmat, mu_i, eps_i = rate_series(t[k])
        rsum = rhs_td(t[k], Y[k], rate_series=rate_series, alpha=alpha,
                      sigma=SIGMA, s=S, rule=rule).sum()
        _, B = mating_flux(Y[k], bpat, bmat, alpha, rule)
        Tpat = Y[k, :N_CLUSTERS].sum(); Tmat = Y[k, N_CLUSTERS:].sum()
        pred = total_balance(I_i, B, Tpat, Tmat, mu_i, eps_i)
        worst = max(worst, abs(rsum - pred))
    print(f"[{label}] time-dependent total-balance max deviation: {worst:.2e}")
    assert worst < 1e-8, f"[{label}] balance check failed: {worst}"


def run_rule(rule):
    d = os.path.join(OUT, f"{STAMP}_hindcast_{rule}")
    os.makedirs(d, exist_ok=True)
    t, Y = run_hindcast(c0, 1980.0, 2024.0, rate_series, alpha=ALPHA,
                        sigma=SIGMA, s=S, rule=rule, dt=0.25)
    check_balance_td(t, Y, ALPHA, rule, rule)
    T = cluster_totals(Y)
    Tpat, Tmat = sex_totals(Y)
    total = T.sum(axis=1)
    share = T / total[:, None]
    cols = {"year": t, "total_M": total,
            "total_pat_M": Tpat, "total_mat_M": Tmat}
    for dd in range(N_CLUSTERS):
        cols[f"C{dd}_pat"] = Y[:, dd]
        cols[f"C{dd}_mat"] = Y[:, N_CLUSTERS + dd]
        cols[f"C{dd}_total"] = T[:, dd]
    pd.DataFrame(cols).to_csv(os.path.join(d, "results.csv"), index=False)
    plot_stacked_area(t, share, os.path.join(d, "stacked_share.png"),
                      title=f"Hindcast 1980-2024 [{RULE_DISPLAY[rule]} rule] "
                            f"alpha={ALPHA}",
                      ylabel="Share of population")
    with open(os.path.join(d, "provenance.md"), "w") as f:
        f.write(f"# Hindcast 1980-2024 [{rule}]\n\n")
        f.write(f"alpha={ALPHA}, eps={EPS}/yr (calibrated), s={S} (measured NCHS), "
                f"sigma(t) from data/immigrant_sex_share.csv (step 5), "
                f"mu=(mu_pat,mu_mat) from data/sex_specific_mortality.csv (step 5)\n")
        f.write(f"IC: 1980 cluster totals from initial_condition_1980.csv, "
                f"split evenly M/F ({c0.sum():.1f} M)\n")
        f.write(f"final: total={total[-1]:.2f} M, "
                f"C_0={T[-1, 0]:.2f} M, pat={Tpat[-1]:.2f} M / mat={Tmat[-1]:.2f} M\n")
    print("wrote", d)
    return t, T, total


def main():
    outs = {rule: run_rule(rule) for rule in RULES}

    # --- validation (rule-independent totals; C_0 = foreign-born) ------------
    vd = os.path.join(OUT, f"{STAMP}_hindcast_validation")
    os.makedirs(vd, exist_ok=True)
    rule0 = RULES[0]
    t, T, total = outs[rule0]
    act = valid[["year", "total_pop_M", "migrant_stock_M",
                 "fb_census_M"]].dropna(subset=["total_pop_M"])

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    a1.plot(t, total, lw=2.5, label="model total (two-sex)")
    a1.scatter(act["year"], act["total_pop_M"], s=22, color="black",
               zorder=5, label="Census/WDI actual")
    a1.set_ylabel("Millions")
    a1.set_title(f"Hindcast 1980-2024: total population "
                 f"(rule-independent; alpha={ALPHA}, eps={EPS}/yr)")
    a1.legend(); a1.grid(True, alpha=0.3)
    a2.plot(t, T[:, 0], lw=2.5, label="model C_0 (foreign-born)")
    ms = valid.dropna(subset=["migrant_stock_M"])
    a2.scatter(ms["year"], ms["migrant_stock_M"], s=22, color="black",
               zorder=5, label="DHS migrant stock")
    fb = valid.dropna(subset=["fb_census_M"])
    a2.scatter(fb["year"], fb["fb_census_M"], s=30, marker="x",
               color="darkred", zorder=5, label="decennial foreign-born")
    a2.set_xlabel("Year"); a2.set_ylabel("Millions")
    a2.set_title("C_0 (foreign-born cluster) vs migrant-stock anchors")
    a2.legend(); a2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(vd, "hindcast_validation.png"), dpi=150)
    plt.close(fig)

    # totals are nearly rule-independent (same I, beta, mu, sinks); the tiny
    # residual (~1e-3 M) comes from B depending on the cluster distribution
    # through the beta ramp and the min() marriage function
    for rule in RULES[1:]:
        dmax = np.abs(outs[rule][2] - total).max()
        print(f"total(rule={rule}) vs total(rule={RULES[0]}): max dev {dmax:.2e} M")
        assert dmax < 1e-2, "totals should be essentially rule-independent"

    yearly = pd.DataFrame({"year": t, "total_M": total})
    yearly = yearly[np.isclose(yearly["year"] % 1, 0)].copy()
    yearly["year"] = yearly["year"].astype(int)
    merged = pd.merge(yearly, act.assign(year=act["year"].astype(int)),
                      on="year", how="inner")
    merged["err_M"] = merged["total_M"] - merged["total_pop_M"]
    merged["rel_err_pct"] = 100 * merged["err_M"] / merged["total_pop_M"]
    merged.to_csv(os.path.join(vd, "hindcast_errors.csv"), index=False)
    worst = merged.loc[merged["err_M"].abs().idxmax()]
    print(f"max total error: {worst['err_M']:+.2f} M ({worst['rel_err_pct']:+.2f}%) "
          f"in {int(worst['year'])}")
    print(f"2024 error: {merged[merged.year == 2024]['err_M'].values[0]:+.2f} M "
          f"({merged[merged.year == 2024]['rel_err_pct'].values[0]:+.2f}%)")
    print("validation dir:", vd)


if __name__ == "__main__":
    main()
