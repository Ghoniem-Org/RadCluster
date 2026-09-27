#!/usr/bin/env python3
"""Build per-sex aggregate mortality mu_pat(t), mu_mat(t), 1650-2025.

Method:
  - 1900-2023: the sex RATIO r(t) = mu_pat(t)/mu_mat(t) is MEASURED from the
    HLD national life tables (Region=0, Ethnicity=0; later-publication
    table preferred, same table choice as the age model,
    data/build_age_inputs.py): r = sum m_male*L_male / sum L_male divided
    by the female equivalent, i.e. the exposure-weighted (population-
    weighted) aggregate of exactly the mu_{s,a} the age model uses.
    IMPORTANT: the life-table stationary population is older than the
    actual growing US population, so its absolute CDR runs 2-4/1000 above
    the calibrated aggregate mu(t) (e.g. 2000: HLD 13.1 vs World Bank
    8.5). The measured quantity carried into the aggregate model is the
    sex ratio; the LEVEL stays the calibrated mu(t) from data/inputs.csv
    (World Bank/Haines), so totals remain exactly consistent with the
    calibrated baseline and the age model's calibrated totals:
        mu_pat = mu * 2*r/(1+r),  mu_mat = mu * 2/(1+r).
    2024 and 2025 hold r(2023) (interpolated).
  - 1650-1899: ESTIMATED. Same split with the mean measured 1900-1919
    ratio r_bar (no sex-specific Haines series exists).

Output: data/sex_specific_mortality.csv: year, mu_pat, mu_mat (per yr),
sex_ratio_r, quality, quality_note.
"""
import csv
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(BASE, "data", "raw", "hld_usa", "USA.csv")
INP = os.path.join(BASE, "data", "inputs.csv")
OUT = os.path.join(BASE, "data", "sex_specific_mortality.csv")


def load_hld_sex_cdr():
    """(year, sex) -> exposure-weighted mean m(x): the per-sex crude
    death rate implied by the HLD life table (1=male, 2=female)."""
    # pass 1: choose the national table per year (same rule as build_age_inputs)
    cands = {}
    with open(RAW) as f:
        for r in csv.DictReader(f):
            if r["Year1"] != r["Year2"]:
                continue
            y = int(r["Year1"])
            if not (1900 <= y <= 2023):
                continue
            if r["Region"] != "0" or r["Ethnicity"] != "0":
                continue
            ref = r["Ref-ID"]
            cands.setdefault(y, set()).add((float(ref.split(".")[0]), ref))
    chosen = {y: max(s)[1] for y, s in cands.items()}
    # pass 2: accumulate m*L and L by (year, sex)
    num, den = {}, {}
    with open(RAW) as f:
        for r in csv.DictReader(f):
            if r["Year1"] != r["Year2"]:
                continue
            y = int(r["Year1"])
            if not (1900 <= y <= 2023) or r["Ref-ID"] != chosen[y]:
                continue
            if r["Region"] != "0" or r["Ethnicity"] != "0":
                continue
            k = (y, int(r["Sex"]))
            num[k] = num.get(k, 0.0) + float(r["m(x)"]) * float(r["L(x)"])
            den[k] = den.get(k, 0.0) + float(r["L(x)"])
    return {k: num[k] / den[k] for k in num}


def main():
    cdr = load_hld_sex_cdr()
    years = sorted({y for (y, s) in cdr})
    mu_pat = {y: cdr[(y, 1)] for y in years}
    mu_mat = {y: cdr[(y, 2)] for y in years}
    r_bar = float(np.mean([mu_pat[y] / mu_mat[y] for y in range(1900, 1920)]))
    print(f"mean 1900-1919 male/female CDR ratio r_bar = {r_bar:.4f}")
    print(f"  mu_pat(1900)={mu_pat[1900]:.6f} mu_mat(1900)={mu_mat[1900]:.6f}; "
          f"mu_pat(2023)={mu_pat[2023]:.6f} mu_mat(2023)={mu_mat[2023]:.6f}")

    inp = pd.read_csv(INP).set_index("year")
    r_bar = float(np.mean([mu_pat[y] / mu_mat[y] for y in range(1900, 1920)]))
    rows = []
    for y in range(1650, 2026):
        mu = float(inp.loc[min(y, 2024), "death_rate"])   # 2025 holds 2024
        if y <= 1899:
            r, q, note = r_bar, "estimated", \
                ("calibrated aggregate mu(t) split by mean 1900-1919 "
                 "measured HLD male/female CDR ratio")
        elif y <= 2023:
            r, q, note = mu_pat[y] / mu_mat[y], "measured", \
                ("calibrated aggregate mu(t) split by measured HLD "
                 "male/female CDR ratio (exposure-weighted life tables)")
        else:
            r, q, note = mu_pat[2023] / mu_mat[2023], "interpolated", \
                "2023 measured HLD sex ratio held"
        mp, mm = mu * 2 * r / (1 + r), mu * 2 / (1 + r)
        rows.append((y, mp, mm, r, q, note))
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["year", "mu_pat", "mu_mat", "sex_ratio_r", "quality",
                    "quality_note"])
        for y, mp, mm, r, q, note in rows:
            w.writerow([y, f"{mp:.8f}", f"{mm:.8f}", f"{r:.6f}", q, note])
    print(f"wrote {OUT} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
