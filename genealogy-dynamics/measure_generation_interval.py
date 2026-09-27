#!/usr/bin/env python3
"""Quantify the aggregate model's missing generation-delay artifact.

Experiment: closed population (I = 0, eps = 0), constant vital rates
(1650 colonial rates and 2024 modern rates), everyone in C_0 at t = 0
with a stable age pyramid (age model) vs the same total in C_0
(aggregate model). The aggregate model's per-capita rates are calibrated
so total births and deaths match the age model at t = 0 -- the ONLY
difference is age structure. We then time the pedigree wavefront:
t_d = first year the share of cluster C_d crosses a threshold, and the
per-rung advance Delta_d = t_{d+1} - t_d over rungs 2..10.

Expected: the age model's front advances ~one rung per human generation
(bounded below by the minimum parental age, 15 yr; the leading edge
tracks the youngest parents so it runs slightly faster than the mean,
MAC). The aggregate model has no such floor: its rung rate scales with
the birth rate, so under high colonial fertility it climbs several times
faster than one rung per generation -- the known artifact.

Outputs: outputs/age_step1/generation_interval.md + figures.
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from genealogy import solver
from genealogy.age import (N_AGE, N_CLUSTERS, N_STATE_AGE, astate_index,
                           stable_age_distribution, mean_age_at_childbearing,
                           birth_inflow_age)
from genealogy.kernels import (SEX_RATIO_AT_BIRTH as S,
                               IMMIGRANT_MALE_FRACTION as SIGMA)

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "outputs", "age_step1")
os.makedirs(OUTDIR, exist_ok=True)

BANDS = ["0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
         "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69",
         "70_74", "75_79", "80_84", "85p"]

P0 = 330.0      # initial total, millions, all in C_0
T_END = 600.0
THRESHOLDS = (1e-3, 1e-2)


def load_rates(year):
    import csv
    with open("data/age_inputs.csv") as f:
        for r in csv.DictReader(f):
            if int(r["year"]) == year:
                beta = np.array([float(r[f"beta_{b}"]) for b in BANDS])
                mu_p = np.array([float(r[f"mu_pat_{b}"]) for b in BANDS])
                mu_m = np.array([float(r[f"mu_mat_{b}"]) for b in BANDS])
                return beta, mu_p, mu_m
    raise RuntimeError(f"{year} row not found")


def wavefront_times(T, thr):
    """T: (steps, n) cluster totals -> [t_1..t_n] first passage (years)."""
    share = T / T.sum(axis=1, keepdims=True)
    ts = []
    for d in range(1, N_CLUSTERS):
        hit = np.nonzero(share[:, d] >= thr)[0]
        ts.append(float(hit[0]) if len(hit) else np.nan)
    return np.array(ts)


def run_front(beta_a, mu_pat_a, mu_mat_a, alpha):
    """Return (mac, tfr, {thr: (age_rung, agg_rung)}, beta_agg)."""
    mac = mean_age_at_childbearing(beta_a)
    tfr = 5.0 * beta_a[3:10].sum()
    prof_p, _ = stable_age_distribution(beta_a, mu_pat_a)
    prof_m, _ = stable_age_distribution(beta_a, mu_mat_a)
    c0a = np.zeros(N_STATE_AGE)
    b = astate_index(0, 0, 0)
    c0a[b:b + N_AGE] = SIGMA * P0 * prof_p
    b = astate_index(0, 1, 0)
    c0a[b:b + N_AGE] = (1 - SIGMA) * P0 * prof_m

    C0 = c0a.reshape(2, N_CLUSTERS, N_AGE)
    B0, _, _ = birth_inflow_age(C0, beta_a, 1.0, S, "shallowest_inheritance")
    D0 = (mu_pat_a[None, :] * C0[0]).sum() + (mu_mat_a[None, :] * C0[1]).sum()
    beta_agg = B0 / C0[0].sum()
    mu_agg = D0 / C0.sum()

    _, Y_age = solver.run_age_scenario(
        c0a, t_end=T_END, dt=1.0, I=0.0, beta_a=beta_a,
        mu_pat_a=mu_pat_a, mu_mat_a=mu_mat_a, eps=0.0, alpha=alpha,
        sigma=SIGMA, s=S, rule="shallowest_inheritance")
    T_age = solver.cluster_totals_age(Y_age)
    c0g = np.zeros(2 * N_CLUSTERS)
    c0g[0] = SIGMA * P0
    c0g[N_CLUSTERS] = (1 - SIGMA) * P0
    _, Y_agg = solver.run_scenario(
        c0g, t_end=T_END, dt=1.0, I=0.0,
        beta_pat=np.full(N_CLUSTERS, beta_agg),
        beta_mat=np.full(N_CLUSTERS, beta_agg),
        mu=mu_agg, eps=0.0, alpha=alpha, sigma=SIGMA, s=S,
        rule="shallowest_inheritance")
    T_agg = Y_agg[:, :N_CLUSTERS] + Y_agg[:, N_CLUSTERS:]
    out = {}
    for thr in THRESHOLDS:
        ta = wavefront_times(T_age, thr)
        tg = wavefront_times(T_agg, thr)
        out[thr] = (float(np.nanmean(np.diff(ta[1:10]))),
                    float(np.nanmean(np.diff(tg[1:10]))))
    return mac, tfr, out, beta_agg, (T_age, T_agg)


def main():
    all_res = {}
    trajs = {}
    for year in (1650, 2024):
        beta_a, mu_pat_a, mu_mat_a = load_rates(year)
        for alpha in (1.0, 0.8):
            mac, tfr, out, b_agg, TT = run_front(beta_a, mu_pat_a,
                                                 mu_mat_a, alpha)
            all_res[(year, alpha)] = (mac, tfr, out, b_agg)
            trajs[(year, alpha)] = TT
            era = "1650" if year == 1650 else "2024"
            print(f"{era} a={alpha}: MAC={mac:.1f} TFR={tfr:.2f} "
                  f"beta_agg={b_agg:.4f}")
            for thr in THRESHOLDS:
                da, dg = out[thr]
                print(f"    thr={thr}: age {da:5.1f} yr/rung | "
                      f"agg {dg:5.1f} yr/rung | age/agg {da/dg:.2f}x")

    with open(os.path.join(OUTDIR, "generation_interval.md"), "w") as f:
        f.write("# Generation-interval artifact: age-structured vs aggregate\n\n")
        f.write("Closed population (I = 0, eps = 0), constant vital rates, "
                "shallowest rule, P0 = 330M all in C_0 (age model: stable age "
                "pyramid). The aggregate model's per-capita rates are "
                "calibrated so total births/deaths match the age model at "
                "t = 0 -- the ONLY difference is age structure.\n\n")
        f.write("`rung` = mean per-rung advance (years) of the C_d share "
                "wavefront over rungs 2..10.\n\n")
        f.write("| era | MAC | TFR | alpha | thr | age rung | agg rung | "
                "age/agg |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        for (year, alpha), (mac, tfr, out, b_agg) in all_res.items():
            era = "1650" if year == 1650 else "2024"
            for thr in THRESHOLDS:
                da, dg = out[thr]
                f.write(f"| {era} | {mac:.1f} | {tfr:.2f} | {alpha} | {thr} "
                        f"| {da:.1f} | {dg:.1f} | {da/dg:.2f}x |\n")
        f.write("\nReading: the age-structured front advances ~one rung per "
                "human generation (bounded below by the minimum parental "
                "age of 15 yr; the leading edge tracks the youngest "
                "parents, so it runs slightly faster than the mean MAC). "
                "The aggregate model has no such floor -- its rung rate "
                "scales with the birth rate. Under high colonial fertility "
                "the aggregate ladder climbs several times faster than one "
                "rung per generation: the known artifact of the ageless "
                "model. Under low modern fertility the aggregate front is "
                "instead mass-accumulation-limited and can run slower; the "
                "artifact is specifically the missing generation-delay "
                "floor.\n")
    print("wrote", os.path.join(OUTDIR, "generation_interval.md"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, year in zip(axes, (1650, 2024)):
        T_age, T_agg = trajs[(year, 0.8)]
        ta = wavefront_times(T_age, 1e-3)
        tg = wavefront_times(T_agg, 1e-3)
        d = np.arange(1, N_CLUSTERS)
        ax.plot(d, ta, "o-", label="age-structured", ms=4)
        ax.plot(d, tg, "s-", label="aggregate", ms=4)
        mac = all_res[(year, 0.8)][0]
        ax.set_xlabel("pedigree cluster d")
        ax.set_title(f"{year} vital rates, $\\alpha=0.8$, shallowest\n"
                     f"(MAC = {mac:.1f} yr; C_d share first $\\geq 10^{{-3}}$)")
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("years since t=0")
    axes[0].legend()
    fig.suptitle("Pedigree wavefront: years to reach 0.1% share by cluster")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "wavefront_comparison.png"), dpi=110)
    print("wrote wavefront_comparison.png")


if __name__ == "__main__":
    main()
