#!/usr/bin/env python3
"""Recalibrate pre-1900 mortality in data/age_inputs.csv.

Problem (found 2026-09-26): the original backcast scaled the 1900-02 HLD
mortality pattern by lam = CDR_target / 17.2 so that the STATIONARY-population
CDR hit the historical CDR target. But the colonial population is growing, so
its pyramid is younger than stationary and its ENDOGENOUS stable CDR came out
2-4/1000 BELOW target (e.g. 1650: 24.2 vs 28.0). Compounded over 250 years,
that 3/1000/yr mortality shortfall inflated the 2025 population to ~440M.

Fix: for each year 1650-1899, solve for lam(t) such that the STABLE
population at (beta(t), lam(t) * mu_1900_02) has endogenous CDR exactly equal
to the historical CDR target from data/inputs.csv. The level is therefore
*calibrated* to Haines' CDR series; the age PATTERN remains the assumed
1900-02 schedule. 1900-2023 measured HLD schedules are untouched.

Outputs: overwrites data/age_inputs.csv (mu_* columns for 1650-1899 only)
and appends to data/age_inputs_provenance.md.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from genealogy.age import (stable_age_distribution, AGING_FLUX, N_AGE,
                           REPRO_BANDS)
from genealogy.kernels import SEX_RATIO_AT_BIRTH as S

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BANDS = ["0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
         "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69",
         "70_74", "75_79", "80_84", "85p"]


def stable_cdr(beta_a, mu_pat_a, mu_mat_a):
    """Endogenous CDR (per 1000) of the stable two-sex population."""
    p, r = stable_age_distribution(beta_a, mu_mat_a, newborn_frac=1 - S)
    T = np.zeros((N_AGE, N_AGE))
    for a in range(N_AGE - 1):
        T[a + 1, a] = AGING_FLUX[a]
    for a in range(N_AGE):
        T[a, a] -= AGING_FLUX[a] + mu_pat_a[a]
    e0 = np.zeros(N_AGE)
    e0[0] = 1.0
    vm = np.linalg.solve(r * np.eye(N_AGE) - T, e0)
    vm /= vm.sum()
    F = (1 - S) * p
    M = S * vm
    tot = (F + M).sum()
    return 1000.0 * ((mu_mat_a * F).sum() + (mu_pat_a * M).sum()) / tot


def main():
    age = pd.read_csv(os.path.join(BASE, "data", "age_inputs.csv"))
    inp = pd.read_csv(os.path.join(BASE, "data", "inputs.csv")).set_index("year")

    ref = age[(age["year"] >= 1900) & (age["year"] <= 1902)]
    mu_pat_ref = ref[[f"mu_pat_{b}" for b in BANDS]].to_numpy().mean(axis=0)
    mu_mat_ref = ref[[f"mu_mat_{b}" for b in BANDS]].to_numpy().mean(axis=0)

    out = []
    for idx, row in age.iterrows():
        y = int(row["year"])
        if 1650 <= y <= 1899:
            beta_a = row[[f"beta_{b}" for b in BANDS]].to_numpy(dtype=float)
            target = float(inp.loc[y, "cdr_per_1000"])

            def f(lam):
                return stable_cdr(beta_a, lam * mu_pat_ref,
                                  lam * mu_mat_ref) - target

            # stable CDR is increasing in lam; bracket and solve
            lo, hi = 0.5, 4.0
            assert f(lo) < 0 < f(hi), f"bracket failed for {y}"
            lam = brentq(f, lo, hi, xtol=1e-6)
            for b, vp, vm_ in zip(BANDS, mu_pat_ref, mu_mat_ref):
                age.loc[idx, f"mu_pat_{b}"] = lam * vp
                age.loc[idx, f"mu_mat_{b}"] = lam * vm_
            out.append((y, lam, target))
        # else: keep measured 1900-2024 rows unchanged
    age.to_csv(os.path.join(BASE, "data", "age_inputs.csv"), index=False)

    lam1650 = out[0][1]
    lam1899 = out[-1][1]
    print(f"lam(1650)={lam1650:.3f} lam(1800)={out[150][1]:.3f} "
          f"lam(1850)={out[200][1]:.3f} lam(1899)={lam1899:.3f}")
    # verify
    for y in (1650, 1800, 1850, 1899):
        r = age[age["year"] == y].iloc[0]
        beta_a = r[[f"beta_{b}" for b in BANDS]].to_numpy(dtype=float)
        mp = r[[f"mu_pat_{b}" for b in BANDS]].to_numpy(dtype=float)
        mm = r[[f"mu_mat_{b}" for b in BANDS]].to_numpy(dtype=float)
        print(f"  {y}: stable CDR={stable_cdr(beta_a, mp, mm):.2f} "
              f"target={inp.loc[y,'cdr_per_1000']:.1f}")

    prov = os.path.join(BASE, "data", "age_inputs_provenance.md")
    with open(prov) as f:
        existing = f.read()
    if "2026-09-26 recalibration" not in existing:
        with open(prov, "a") as f:
            f.write("\n## 2026-09-26 recalibration (pre-1900 mortality level)\n\n")
            f.write("The original backcast set lam = CDR_target/17.2 to match the "
                    "*stationary*-population CDR. Because the colonial population "
                    "is growing (young pyramid), its endogenous stable CDR came "
                    "out 2-4/1000 below target, inflating 2025 population to "
                    "~440M. Recalibrated: for each year 1650-1899, lam(t) solves "
                    "stable_CDR(beta(t), lam*mu_1900_02) = CDR_target(t) "
                    "(brentq). lam(1650)=%.3f, lam(1899)=%.3f. Level is now "
                    "*calibrated* to Haines CDR via inputs.csv; age pattern "
                    "remains *assumed* (1900-02 schedule).\n" % (lam1650, lam1899))
    print("rewrote data/age_inputs.csv (1650-1899 mu columns)")


if __name__ == "__main__":
    main()
