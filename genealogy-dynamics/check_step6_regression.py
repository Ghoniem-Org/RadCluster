#!/usr/bin/env python3
"""Step-6 regression: GENEALOGY_MAX_DEPTH=12 re-runs must reproduce the pinned
step-5 outputs.

Two-sex: compare new n=12 run's summary_all.csv against
  outputs/20260926_215055_colonial/summary_all.csv
Regional: compare new n=12 run's summary against
  outputs/20260926_220052_colonial_regional/summary_*.csv
Age: outputs/age_step1/colonial_age (n=12 re-run, overwrote the dir)
  vs outputs/age_step1/colonial_age_step5_n12 (pinned copy).

Pass: max |relative difference| < 1e-10 on totals and Cn shares.
Usage: python3 check_step6_regression.py <n12_twosex_dir> <n12_regional_dir>
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-10


def rel(a, b):
    return abs(a - b) / max(1e-300, abs(b))


def check_two_sex(newdir):
    pin = os.path.join(BASE, "outputs", "20260926_215055_colonial",
                       "summary_all.csv")
    pinned = pd.read_csv(pin)
    new = pd.read_csv(os.path.join(newdir, "summary_all.csv"))
    worst = 0.0
    cols = ["total_2025_M", "pat_2025_M", "mat_2025_M", "pat_share_2025",
            "Cn_share_2025", "balance_max_dev"]
    for _, pr in pinned.iterrows():
        nr = new[(new["rule"] == pr["rule"])
                 & (new["alpha"] == pr["alpha"])].iloc[0]
        for c in cols:
            if c == "balance_max_dev":
                worst = max(worst, abs(nr[c] - pr[c]))
                continue
            worst = max(worst, rel(nr[c], pr[c]))
    print(f"two-sex n=12 regression worst rel dev: {worst:.3e} "
          f"({'PASS' if worst < TOL else 'FAIL'})")
    return worst < TOL


def check_regional(newdir):
    worst = 0.0
    pin_dir = os.path.join(BASE, "outputs", "20260926_220052_colonial_regional")
    pin_files = [f for f in os.listdir(pin_dir) if f.startswith("summary_")]
    for pf in pin_files:
        pinned = pd.read_csv(os.path.join(pin_dir, pf))
        new = pd.read_csv(os.path.join(newdir, pf))
        for c in pinned.columns:
            if pinned[c].dtype == object or "dev" in c:
                continue
            for i in range(len(pinned)):
                worst = max(worst, rel(float(new[c].iloc[i]),
                                       float(pinned[c].iloc[i])))
    print(f"regional n=12 regression worst rel dev: {worst:.3e} "
          f"({'PASS' if worst < TOL else 'FAIL'})")
    return worst < TOL


def check_age():
    import filecmp
    pin = os.path.join(BASE, "outputs", "age_step1", "colonial_age_step5_n12")
    new = os.path.join(BASE, "outputs", "age_step1", "colonial_age")
    csvs = [f for f in os.listdir(pin) if f.endswith(".csv")]
    worst = 0.0
    ok_files = True
    for f in csvs:
        a = pd.read_csv(os.path.join(pin, f))
        b = pd.read_csv(os.path.join(new, f))
        if list(a.columns) != list(b.columns) or len(a) != len(b):
            print(f"  age: {f} shape mismatch -> FAIL")
            ok_files = False
            continue
        num = a.select_dtypes("number")
        for c in num.columns:
            d = np.abs(b[c].to_numpy() - a[c].to_numpy())
            denom = np.maximum(np.abs(a[c].to_numpy()), 1e-300)
            worst = max(worst, float(np.max(d / denom)))
    print(f"age n=12 regression worst rel dev over {len(csvs)} CSVs: "
          f"{worst:.3e} ({'PASS' if ok_files and worst < TOL else 'FAIL'})")
    return ok_files and worst < TOL


def main():
    n12_twosex = sys.argv[1] if len(sys.argv) > 1 else None
    n12_regional = sys.argv[2] if len(sys.argv) > 2 else None
    ok = True
    if n12_twosex:
        ok &= check_two_sex(n12_twosex)
    if n12_regional:
        ok &= check_regional(n12_regional)
    ok &= check_age()
    print("REGRESSION:", "ALL PASS" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
