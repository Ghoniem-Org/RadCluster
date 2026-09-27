#!/usr/bin/env python3
"""Step-5 male-fertility sensitivity (age model).

Births are maternal: the child sex split (s, 1-s) and the min-form mating
kernel are untouched. The question is only the paternal mating-weight
schedule: the DEFAULT weights fathers with the same ASFR shape as mothers
(beta_pat_a = beta_a, ASSUMED proxy -- no male ASFR series exists). Here
we compare, on the identical Case-1 1650-2025 age-model trajectory:

    lam = 0.0 : equal-schedule proxy (the production default)
    lam = 0.4 : fathers ~2 yr older -- beta_pat_a = (1-lam)*beta_a +
                lam*beta_a rolled one 5-yr band (genealogy/age.py
                male_shift_lam hook)

Report: the difference in 2025 deepest-cluster (C_n) total shares. If
negligible, the equal-schedule proxy is retained with its label.
Outputs: outputs/<ts>_male_fertility/{summary.csv, male_fertility.md}.
"""
import os
from datetime import datetime

import numpy as np
import pandas as pd

from genealogy.populations import N_CLUSTERS, RULES, MAX_DEPTH
from genealogy.kernels import SEX_RATIO_AT_BIRTH
from genealogy import solver
from genealogy.age import N_AGE, initial_condition_age
from genealogy.sex_rates import load_sigma_series
from run_colonial_age import AgeSeries  # local age rate series (step 1)

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs",
                      f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_male_fertility")
os.makedirs(OUTDIR, exist_ok=True)

S = SEX_RATIO_AT_BIRTH
SIGMA = load_sigma_series()
ALPHAS = [0.0, 0.8, 0.9]
T0, T1 = 1650.0, 2025.0
EPS = 0.00075


def run_variant(rule, alpha, lam):
    """Run with the male-fertility shift applied (lam) or not."""
    series = AgeSeries()
    _, beta0, mu_p0, mu_m0, _ = series.at_year(1650)
    c0 = initial_condition_age(0.050368, beta0, mu_p0, mu_m0,
                               sigma=SIGMA(1650.0))
    t, Y = solver.run_age_hindcast(
        c0, T0, T1, series, alpha, SIGMA, S, rule, dt=1.0,
        male_shift_lam=lam)
    T = solver.cluster_totals_age(Y)
    total = T.sum(axis=1)
    return total[-1], (T[-1, MAX_DEPTH] / total[-1])


def main():
    rows = []
    for rule in RULES:
        for alpha in ALPHAS:
            tot0, cn0 = run_variant(rule, alpha, 0.0)
            tot1, cn1 = run_variant(rule, alpha, 0.4)
            rows.append(dict(rule=rule, alpha=alpha,
                             total_proxy_M=tot0, Cn_proxy=cn0,
                             total_shift_M=tot1, Cn_shift=cn1,
                             d_Cn_pp=(cn1 - cn0) * 100,
                             d_total_pct=(tot1 - tot0) / tot0 * 100))
            r = rows[-1]
            print(f"[{rule} a={alpha}] Cn proxy={cn0*100:.3f}% "
                  f"shifted={cn1*100:.3f}% d={r['d_Cn_pp']:+.4f}pp "
                  f"total {r['d_total_pct']:+.4f}%", flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "summary.csv"), index=False)
    with open(os.path.join(OUTDIR, "male_fertility.md"), "w") as f:
        f.write("# Male-fertility sensitivity (step 5, age model)\n\n")
        f.write("Equal-schedule proxy (lam=0.0) vs fathers ~2 yr older "
                "(lam=0.4). Births maternal, min-form kernel kept.\n\n")
        f.write("| rule | alpha | Cn proxy | Cn shifted | d (pp) |\n"
                "|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['rule']} | {r['alpha']} | {r['Cn_proxy']*100:.3f}% | "
                    f"{r['Cn_shift']*100:.3f}% | {r['d_Cn_pp']:+.4f} |\n")
    print("wrote", OUTDIR)


if __name__ == "__main__":
    main()
