#!/usr/bin/env python3
"""Step 8 — Independent validation of class shares against Census data.

Checks (all numbers read from archived outputs / public downloads; no new
model runs):
  A. REGIONAL ANCESTRY (primary): 2023 ACS 1-year B04006 "American" ancestry %
     by Census region vs regional model 2025 deepest-cluster (C_14) share by
     region, shallowest/deepest inheritance, alpha=0.8.
  C. FOREIGN-BORN CONSISTENCY (fallback): 2023 ACS 1-year B05002 foreign-born %
     by region vs regional model 2025 C_0 share by region.  Labelled a
     consistency check (immigration data fed the model), NOT validation.
  B. SURNAME CHECK: genuine attempt documented in docs/STEP_LOG.md; no
     defensible operationalization with the public 2010 surname list alone
     (national-only, no surname dating, model predicts shares not surname
     frequencies) -> falls back to A + C per the work-package spec.

Writes:
- drive_doc/figs/validation_american_ancestry.png
- drive_doc/figs/validation_foreignborn.png
- outputs/<stamp>_step8_validation/results_step8.csv
"""

import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "data", "raw", "validation_step8")
FIGDIR = os.path.join(HERE, "drive_doc", "figs")
OUT = os.path.join(HERE, "outputs", "20260926_224500_step8_validation")
os.makedirs(OUT, exist_ok=True)

REGIONS = ["Northeast", "Midwest", "South", "West"]
GIDS = {"Northeast": "0200000US1", "Midwest": "0200000US2",
        "South": "0200000US3", "West": "0200000US4"}


def load_dat(fn):
    with open(os.path.join(RAW, fn), encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="|"))


def ratio_pct_moe(num, num_moe, den, den_moe):
    """ACS handbook MOE for a proportion (numerator subset of denominator).

    Negative MOE = Census sentinel for a controlled total (MOE undefined);
    controlled denominators are fixed, so their sampling variance is 0.
    """
    if den_moe < 0:
        den_moe = 0.0
    if num_moe < 0:
        num_moe = 0.0
    p = num / den
    var = max(num_moe ** 2 - (p ** 2) * (den_moe ** 2), 0.0)
    return 100.0 * p, 100.0 * np.sqrt(var) / den


# ---- observed: ACS 2023 1-year -------------------------------------------
b04 = {r["GEO_ID"]: r for r in load_dat("acsdt1y2023-b04006.dat")}
b05 = {r["GEO_ID"]: r for r in load_dat("acsdt1y2023-b05002.dat")}

obs_ame, obs_fb = {}, {}
for r in REGIONS + ["US"]:
    gid = GIDS[r] if r != "US" else "0100000US"
    a = b04[gid]; f = b05[gid]
    obs_ame[r] = ratio_pct_moe(float(a["B04006_E005"]), float(a["B04006_M005"]),
                              float(a["B04006_E001"]), float(a["B04006_M001"]))
    obs_fb[r] = ratio_pct_moe(float(f["B05002_E013"]), float(f["B05002_M013"]),
                              float(f["B05002_E001"]), float(f["B05002_M001"]))

# ---- modeled: regional colonial run 2025, alpha=0.8 ------------------------
REGIONAL_DIR = os.path.join(HERE, "outputs", "20260927_161913_colonial_regional")


def regional_summary(rule):
    fn = os.path.join(REGIONAL_DIR, "summary_all_regional.csv")
    with open(fn) as fh:
        for row in csv.DictReader(fh):
            if row["rule"] == rule and row["alpha"] == "0.8":
                return row
    raise RuntimeError(rule)

SHALLOW, DEEP = "shallowest_inheritance", "deepest_inheritance"
mod_c14, mod_c0 = {}, {}
for rule in (SHALLOW, DEEP):
    s = regional_summary(rule)
    mod_c14[rule] = {r: 100.0 * float(s[f"Cn_share_{r}_2025"]) for r in REGIONS}
    mod_c14[rule]["US"] = 100.0 * float(s["Cn_share_2025"])
fn = os.path.join(REGIONAL_DIR,
                  "results_regional_shallowest_inheritance_alpha0.8.csv")
with open(fn) as fh:
    rows = list(csv.DictReader(fh))
last = rows[-1]
for r in REGIONS:
    mod_c0[r] = 100.0 * float(last[f"share_C0_{r}"])
mod_c0["US"] = 100.0 * float(last["share_C0"])
# national cumulative depth shares (shallowest) for proxy-breadth context
cum = {}
for d0 in (8, 10, 12):
    cum[d0] = 100.0 * (sum(float(last[f"share_C{d}"]) for d in range(d0, 13))
                       + float(last["Cn_share"]))
cum[14] = 100.0 * float(last["Cn_share"])


def spearman(a, b):
    """Spearman rho between two dicts over the same keys (rank by value)."""
    keys = list(a)
    ra = {k: i for i, k in enumerate(sorted(keys, key=lambda k: a[k]))}
    rb = {k: i for i, k in enumerate(sorted(keys, key=lambda k: b[k]))}
    d2 = sum((ra[k] - rb[k]) ** 2 for k in keys)
    n = len(keys)
    return 1.0 - 6.0 * d2 / (n * (n ** 2 - 1))


cats = REGIONS  # rank checks use regions only
rho_ame_shallow = spearman({r: obs_ame[r][0] for r in cats},
                            {r: mod_c14[SHALLOW][r] for r in cats})
rho_ame_deep = spearman({r: obs_ame[r][0] for r in cats},
                        {r: mod_c14[DEEP][r] for r in cats})
rho_fb = spearman({r: obs_fb[r][0] for r in cats},
                  {r: mod_c0[r] for r in cats})

# ---- results table ---------------------------------------------------------
with open(os.path.join(OUT, "results_step8.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["region", "obs_american_pct", "obs_american_moe",
                "model_c14_shallowest_pct", "model_c14_deepest_pct",
                "obs_foreignborn_pct", "obs_foreignborn_moe", "model_c0_pct"])
    for r in REGIONS + ["US"]:
        w.writerow([r, round(obs_ame[r][0], 2), round(obs_ame[r][1], 2),
                    round(mod_c14[SHALLOW][r], 3), round(mod_c14[DEEP][r], 3),
                    round(obs_fb[r][0], 2), round(obs_fb[r][1], 2),
                    round(mod_c0[r], 2)])

print("American-ancestry % (obs) vs C_14 % (model), by region:")
for r in REGIONS + ["US"]:
    print(f"  {r:10s} obs {obs_ame[r][0]:5.2f}±{obs_ame[r][1]:.2f} | "
          f"shallowest {mod_c14[SHALLOW][r]:6.3f} | deepest {mod_c14[DEEP][r]:6.3f}")
print(f"  Spearman rank rho (regions): shallowest {rho_ame_shallow:.2f}, "
      f"deepest {rho_ame_deep:.2f}")
print("  shallowest cumulative national: " +
      ", ".join(f"C{d}+={cum[d]:.2f}%" for d in (8, 10, 12, 14)))
print("Foreign-born % (obs) vs C_0 % (model), by region:")
for r in REGIONS + ["US"]:
    d = mod_c0[r] - obs_fb[r][0]
    print(f"  {r:10s} obs {obs_fb[r][0]:5.2f}±{obs_fb[r][1]:.2f} | "
          f"model {mod_c0[r]:6.2f} | gap {d:+5.2f}pp "
          f"({100*d/obs_fb[r][0]:+.0f}% rel)")
print(f"  Spearman rank rho (regions): {rho_fb:.2f}")

# ---- figures ----------------------------------------------------------------
x = np.arange(len(REGIONS) + 1)
labels = ["NE", "MW", "S", "W", "US"]
keys = REGIONS + ["US"]

fig, ax = plt.subplots(figsize=(10, 5.6))
w = 0.24
o = np.array([obs_ame[k][0] for k in keys]); oe = np.array([obs_ame[k][1] for k in keys])
s = np.array([mod_c14[SHALLOW][k] for k in keys])
ln = np.array([mod_c14[DEEP][k] for k in keys])
ax.bar(x - w, o, w, yerr=oe, capsize=3, color="#2ca02c", edgecolor="black",
       linewidth=0.5, label='Observed "American" ancestry (ACS 2023 1-yr B04006)')
ax.bar(x, s, w, color="#1f77b4", edgecolor="black", linewidth=0.5,
       label="Model C$_{14}$ share, shallowest (α=0.8)")
ax.bar(x + w, ln, w, color="#ff7f0e", edgecolor="black", linewidth=0.5,
       label="Model C$_{14}$ share, deepest (α=0.8)")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("percent of population")
ax.set_title('Check A — "American" ancestry (observed) vs deepest-cluster share (modeled)\n'
             f"by Census region; rank-order Spearman ρ={rho_ame_shallow:.2f} (shallowest)")
ax.legend(fontsize=8.5); ax.grid(True, axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "validation_american_ancestry.png"), dpi=150)
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 5.6))
o = np.array([obs_fb[k][0] for k in keys]); oe = np.array([obs_fb[k][1] for k in keys])
m = np.array([mod_c0[k] for k in keys])
ax.bar(x - w / 2, o, w, yerr=oe, capsize=3, color="#2ca02c", edgecolor="black",
       linewidth=0.5, label="Observed foreign-born (ACS 2023 1-yr B05002)")
ax.bar(x + w / 2, m, w, color="#9467bd", edgecolor="black", linewidth=0.5,
       label="Model C$_{0}$ share (α=0.8, rule-independent)")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("percent of population")
ax.set_title("Check C — foreign-born (observed) vs C$_{0}$ (modeled), by Census region\n"
             f"consistency check only; rank-order Spearman ρ={rho_fb:.2f}")
ax.legend(fontsize=9); ax.grid(True, axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "validation_foreignborn.png"), dpi=150)
plt.close(fig)

with open(os.path.join(OUT, "provenance.md"), "w") as fh:
    fh.write(
        "# Step 8 validation — provenance\n\n"
        "Observed data (all measured survey estimates, downloaded 2026-09-26):\n"
        "- 2023 ACS 1-year Detailed Table B04006 (People Reporting Ancestry), "
        "variable B04006_005E = \"American\" ancestry, B04006_001E = total; "
        "table-based summary file acsdt1y2023-b04006.dat from "
        "https://www2.census.gov/programs-surveys/acs/summary_file/2023/table-based-SF/data/1YRData/ "
        "(Census API data endpoint required a key as of 2026-09; summary file used instead).\n"
        "- 2023 ACS 1-year Detailed Table B05002 (Place of Birth by Nativity), "
        "B05002_013E = foreign-born, B05002_001E = total; acsdt1y2023-b05002.dat, same source.\n"
        "- Variable labels from https://api.census.gov/data/2023/acs/acs1/variables.json "
        "(metadata endpoint, no key).\n"
        "- Region rows read directly (GEO_ID 0200000US1-4); 6 small states "
        "(AK, ND, SD, VT, WV, WY) absent at state level in the B04006 .dat but "
        "included in the published region rows.\n"
        "- MOEs via the ACS handbook ratio formula (numerator subset of denominator).\n"
        "- Census 2010 surname list Names_2010Census.csv (151,671 surnames, >=100 bearers) "
        "from https://www2.census.gov/topics/genealogy/2010surnames/names.zip — "
        "downloaded and inspected; no defensible surname check operationalized (see STEP_LOG step 8).\n"
        "Modeled data: outputs/20260927_161913_colonial_regional "
        "(summary_all_regional.csv alpha=0.8 rows; "
        "results_regional_shallowest_inheritance_alpha0.8.csv 2025 row). No new model runs.\n")
print("wrote", OUT)
