#!/usr/bin/env python3
"""Figures for Step 1 (age structure): input schedules and their history.

- asfr_schedules.png : ASFR by single year of age, 1800 (backcast),
  1950, 1980, 2023 (WPP estimates).
- mac_tfr_history.png : mean age at childbearing and TFR, 1650-2024.
- mortality_schedules.png : sex-specific death rates by age band,
  1900 / 1950 / 2023 (HLD measured), log scale.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from genealogy.age import (AGE_MIDPOINTS, REPRO_BANDS,
                           mean_age_at_childbearing)

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "outputs", "age_step1")
os.makedirs(OUTDIR, exist_ok=True)
BANDS = ["0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
         "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69",
         "70_74", "75_79", "80_84", "85p"]

age = pd.read_csv(os.path.join(BASE, "data", "age_inputs.csv"))
age = age.sort_values("year").reset_index(drop=True)


def row(year):
    return age[age["year"] == year].iloc[0]


def beta_of(r):
    return np.array([r[f"beta_{b}"] for b in BANDS])


def mu_of(r, sex):
    return np.array([r[f"mu_{sex}_{b}"] for b in BANDS])


# --- ASFR schedules ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4.8))
for y, style, lab in ((1800, "--", "1800 (backcast: TFR 7.04 x 1950 shape)"),
                      (1950, "-", "1950 (WPP estimate)"),
                      (1980, "-", "1980 (WPP estimate)"),
                      (2023, "-", "2023 (WPP estimate)")):
    b = beta_of(row(y))
    x = AGE_MIDPOINTS[list(REPRO_BANDS)]
    ax.plot(x, 1000 * b[list(REPRO_BANDS)], style, label=lab)
ax.set_xlabel("age of mother (years)")
ax.set_ylabel("births per 1,000 women per year")
ax.set_title("Age-specific fertility rate by mother's age")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "asfr_schedules.png"), dpi=110)

# --- MAC and TFR history ------------------------------------------------------
years = age["year"].to_numpy()
mac = np.array([mean_age_at_childbearing(beta_of(age.iloc[i]))
                for i in range(len(age))])
tfr = np.array([5.0 * beta_of(age.iloc[i])[list(REPRO_BANDS)].sum()
                for i in range(len(age))])
fig, ax1 = plt.subplots(figsize=(9, 4.5))
ax1.plot(years, mac, "b-", lw=1.2, label="MAC")
ax1.set_xlabel("year")
ax1.set_ylabel("mean age at childbearing (years)", color="b")
ax1.tick_params(axis="y", labelcolor="b")
ax1.grid(True, alpha=0.3)
ax2 = ax1.twinx()
ax2.plot(years, tfr, "r-", lw=1.2, label="TFR")
ax2.set_ylabel("total fertility rate", color="r")
ax2.tick_params(axis="y", labelcolor="r")
ax2.axvspan(1650, 1800, color="gray", alpha=0.12, label="assumed")
ax2.axvspan(1800, 1950, color="orange", alpha=0.10, label="interpolated")
ax2.legend(fontsize=8, loc="upper right")
fig.suptitle("Mean age at childbearing and TFR, 1650-2024 "
             "(shading = input quality)")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "mac_tfr_history.png"), dpi=110)

# --- mortality schedules ------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
for ax, sex, title in zip(axes, ("pat", "mat"), ("Male", "Female")):
    for y, style in ((1900, "--"), (1950, "-"), (2023, "-")):
        mu = mu_of(row(y), sex)
        ax.plot(AGE_MIDPOINTS, mu, style, label=f"{y}")
    ax.set_yscale("log")
    ax.set_xlabel("age (years)")
    ax.set_title(title)
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=8)
axes[0].set_ylabel("death rate per person-year (log scale)")
fig.suptitle("Age-specific mortality (HLD life tables, measured 1900-2023)")
fig.tight_layout()
fig.savefig(os.path.join(OUTDIR, "mortality_schedules.png"), dpi=110)

print("figures written to", OUTDIR)
print("MAC 1800/1900/1950/1980/2023:",
      " ".join(f"{m:.1f}" for m in
               [mean_age_at_childbearing(beta_of(row(y)))
                for y in (1800, 1900, 1950, 1980, 2023)]))
