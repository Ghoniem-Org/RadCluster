#!/usr/bin/env python3
"""Colonial run with the age-structured two-sex model: 1650 -> 2025.

- Age-specific rates: data/age_inputs.csv (1650-2024; 2025 holds 2024).
  Immigration inflow: data/inputs.csv (same residual series as the
  aggregate colonial run). Emigration eps = 0.00075/yr (calibrated, same).
- IC 1650: 50,368 people in C_0, sexes split sigma/(1-sigma), each sex on
  its stable age distribution at 1650 rates (avoids an age-structure
  spin-up transient); C_1..C_n = 0.
- alpha in {0.0, 0.8, 0.9} (0.8 central); rules shallowest/deepest.
- Immigrant age profile: the assumed fixed shares
  (genealogy/age.py:IMMIGRANT_AGE_PROFILE).

Writes outputs/age_step1/colonial_age/: per-rule/per-alpha cluster-share
CSVs, 25-year snapshots, heatmaps, anchor validation, colonial_age.md.
"""
import os

import numpy as np
import pandas as pd

from genealogy.populations import N_CLUSTERS, RULES
from genealogy.kernels import SEX_RATIO_AT_BIRTH
from genealogy.sex_rates import load_sigma_series

SIGMA = load_sigma_series()  # step 5: sigma(t) immigrant male share
from genealogy import solver
from genealogy.age import (N_AGE, initial_condition_age,
                           IMMIGRANT_AGE_PROFILE)

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs", "age_step1", "colonial_age")
os.makedirs(OUTDIR, exist_ok=True)

S = SEX_RATIO_AT_BIRTH
EPS = 0.00075  # calibrated (data/emigration_calibration.txt)
T0, T1 = 1650, 2025
ALPHAS = (0.0, 0.8, 0.9)
BANDS = ["0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
         "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69",
         "70_74", "75_79", "80_84", "85p"]


class AgeSeries:
    """Callable(t) -> (I, beta_a, mu_pat_a, mu_mat_a, eps); linear in time
    between annual rows of age_inputs.csv; immigration from inputs.csv."""

    def __init__(self):
        age = pd.read_csv(os.path.join(BASE, "data", "age_inputs.csv"))
        age = age.sort_values("year").reset_index(drop=True)
        self.years = age["year"].to_numpy(dtype=float)
        self.beta = age[[f"beta_{b}" for b in BANDS]].to_numpy(dtype=float)
        self.mu_p = age[[f"mu_pat_{b}" for b in BANDS]].to_numpy(dtype=float)
        self.mu_m = age[[f"mu_mat_{b}" for b in BANDS]].to_numpy(dtype=float)
        inp = pd.read_csv(os.path.join(BASE, "data", "inputs.csv"))
        inp = inp.sort_values("year").reset_index(drop=True)
        self.iy = inp["year"].to_numpy(dtype=float)
        self.imm = inp["immigration_M"].to_numpy(dtype=float)

    def _interp(self, xs, vs, t):
        return np.array([np.interp(t, xs, vs[:, j])
                         for j in range(vs.shape[1])])

    def __call__(self, t):
        tc = min(max(t, self.years[0]), self.years[-1])
        beta_a = self._interp(self.years, self.beta, tc)
        mu_p = self._interp(self.years, self.mu_p, tc)
        mu_m = self._interp(self.years, self.mu_m, tc)
        I = float(np.interp(tc, self.iy, self.imm))
        return I, beta_a, mu_p, mu_m, EPS

    def at_year(self, y):
        return self(min(max(y, T0), 2024))


def main():
    series = AgeSeries()
    anchors = pd.read_csv(os.path.join(BASE, "data", "validation_long.csv"))
    anchors = anchors[["year", "pop_M"]]

    # IC: stable age pyramid at 1650 rates
    _, beta0, mu_p0, mu_m0, _ = series.at_year(1650)
    c0 = initial_condition_age(0.050368, beta0, mu_p0, mu_m0,
                               sigma=SIGMA(1650.0))

    summary = []
    for rule in RULES:
        for alpha in ALPHAS:
            tag = f"{rule}_alpha{alpha}"
            print(f"running {tag} ...", flush=True)
            t, Y = solver.run_age_hindcast(
                c0, T0, T1, series, alpha, SIGMA, S, rule, dt=1.0)
            T = solver.cluster_totals_age(Y)          # (steps, n)
            total = T.sum(axis=1)
            share = T / total[:, None]
            years = np.round(t).astype(int)

            df = pd.DataFrame({"year": years, "total_M": total})
            for d in range(N_CLUSTERS):
                df[f"share_C{d}"] = share[:, d]
            df.to_csv(os.path.join(OUTDIR, f"results_{tag}.csv"), index=False)

            snap = df[df["year"] % 25 == 0].copy()
            snap.to_csv(os.path.join(OUTDIR, f"snapshots_{tag}.csv"),
                        index=False)

            # anchor validation
            m = pd.merge(anchors, df[["year", "total_M"]], on="year",
                         how="left")
            m["rel_err"] = (m["total_M"] - m["pop_M"]) / m["pop_M"]
            m.to_csv(os.path.join(OUTDIR, f"anchor_validation_{tag}.csv"),
                     index=False)
            worst = m.loc[m["rel_err"].abs().idxmax()]

            # heatmap (cluster share vs time)
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(9, 5))
            im = ax.imshow(share.T, aspect="auto", origin="lower",
                           extent=[T0, T1, -0.5, N_CLUSTERS - 0.5],
                           cmap="magma", norm=matplotlib.colors.LogNorm(
                               vmin=1e-4, vmax=1.0))
            ax.set_xlabel("year")
            ax.set_ylabel("pedigree cluster d")
            ax.set_title(f"Age-structured colonial run: {rule}, "
                         f"alpha={alpha} (cluster shares, log scale)")
            fig.colorbar(im, ax=ax, label="share")
            fig.tight_layout()
            fig.savefig(os.path.join(OUTDIR, f"heatmap_{tag}.png"), dpi=100)
            plt.close(fig)

            fin = share[-1]
            summary.append((rule, alpha, total[-1], fin[-1],
                            worst["year"], worst["rel_err"]))
            print(f"  2025 total {total[-1]:.1f}M, C_n share {fin[-1]*100:.2f}%, "
                  f"worst anchor {worst['rel_err']*100:+.1f}% ({worst['year']:.0f})")

    with open(os.path.join(OUTDIR, "colonial_age.md"), "w") as f:
        f.write("# Colonial run 1650-2025 (age-structured two-sex model)\n\n")
        f.write(f"{2 * N_CLUSTERS * 18} states ({N_CLUSTERS} clusters x 2 sexes "
                f"x 18 age bands). Rates: "
                "data/age_inputs.csv; immigration: data/inputs.csv; "
                f"eps={EPS}/yr (calibrated); s={S} (measured NCHS); "
                "sigma(t) from data/immigrant_sex_share.csv (step 5). "
                "Immigrant age profile: assumed fixed "
                "shares (see genealogy/age.py).\n\n")
        f.write("## 2025 deepest-cluster (C_n) total share\n\n")
        f.write("| rule | alpha=0.0 | **0.8 (central)** | 0.9 |\n")
        f.write("|---|---|---|---|\n")
        for rule in RULES:
            row = [s for s in summary if s[0] == rule]
            cells = " | ".join(f"{s[3]*100:.2f}%" for s in row)
            f.write(f"| {rule} | {cells} |\n")
        f.write("\nAggregate-model baseline (same scenario, no age structure):\n\n")
        f.write("| rule | alpha=0.0 | **0.8** | 0.9 |\n")
        f.write("|---|---|---|---|\n")
        f.write("| shallowest_inheritance | 0.00% | 1.33% | 3.31% |\n")
        f.write("| deepest_inheritance | 0.24% | 3.46% | 5.01% |\n")
        f.write("\n## 2025 total population (M) and anchor validation\n\n")
        for rule, alpha, tot, cn, wy, we in summary:
            f.write(f"- {rule} alpha={alpha}: 2025 total {tot:.1f}M; worst "
                    f"anchor rel err {we*100:+.1f}% ({wy:.0f})\n")
    print("wrote", OUTDIR)


if __name__ == "__main__":
    main()
