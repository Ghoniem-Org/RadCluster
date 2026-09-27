#!/usr/bin/env python3
"""Build data/assortativity_schedule.csv: historically grounded alpha(t).

Pedigree-cluster assortativity alpha is not directly observed historically,
so the schedule is built from the intermarriage literature via the
approximate proxy

    endogamy ~= 1 - k * (intermarriage share among newlyweds).

Construction (era by era; every value labeled assumed / estimated /
measured-series / interpolated -- an assumed anchor is NEVER presented
as an estimate):

1. PRE-1960 anchors (no continuous intermarriage series exists in
   comparable units; values estimated from qualitative literature):
   - 1650: 0.95 ASSUMED -- colonial: small, homogeneous settler
     communities; no marriage data exist.
   - 1776: 0.94 ASSUMED -- independence; colonial endogamy assumed
     persistent.
   - 1850: 0.93 ESTIMATED -- 19th-century European ethnic endogamy strong
     (Alba & Golden 1986).
   - 1880: 0.93 ESTIMATED -- 60-70% of immigrants married endogamously
     (Spoerlein et al. 2014); strong third-generation endogamy
     (Logan & Shin 2012).
   - 1910: 0.94 ESTIMATED -- endogamy "castelike" for new ethnics at the
     turn of the century (Pagnini & Morgan 1990); Draschler (1921, NYC
     1908-12 licenses): exogamy rose first->second generation but the
     first generation stayed strongly endogamous.
   - 1940: 0.945 ESTIMATED -- pre-WWII; intermarriage rose "dramatically
     after World War II" (Alba 1983), so pre-WWII endogamy was still high.
   Interpolated linearly between anchors.

2. 1960-2015 (MEASURED series + CALIBRATED mapping):
   alpha(t) = 1 - k * x(t), where x(t) is the Pew newlywed
   interracial/interethnic share:
     1960: 2.4, 1965: 2.9, 1970: 4.0, 1975: 5.0, 1980: 6.7, 1985: 7.3,
     1990: 8.3, 1995: 9.3, 2000: 11.2, 2005: 13.1, 2008: 14.5,
     2009: 15.4, 2010: 15.1, 2011: 15.5 (Pew 2023, "Interracial Marriage
     Grows"); 1967: 3.0 (Pew 2017, Loving v. Virginia); 2015: 17.0
     (Pew 2017).
   The slope k is CALIBRATED: normalized so the schedule reproduces the
   CPS-calibrated cluster assortativity alpha_CPS = 0.70 at 2013 (Step 3).
   x(2013) = 16.25% (linear between 2011 and 2015) -> k = 0.30/0.1625
   = 1.8462. Interpretation: 1 pp of racial/ethnic intermarriage maps to
   ~1.85 pp of lost same-cluster mating -- the extra ~0.85 captures
   same-race/ethnicity couples who still differ in pedigree cluster
   (different immigration waves), consistent with alpha_CPS (0.70)
   sitting below 1 - x(2013) = 0.84.
   Consistency property (checked by this script): alpha(2013) = 0.700.

3. 2015-2025 (ESTIMATED, moderated trend): 2015: 17.0% -> 2025: 18.8%,
   slope 0.18 pp/yr = half the 2008-2015 rate (0.36 pp/yr). Pew (2017)
   notes the rise was slowing (Asian intermarriage ticked down,
   Hispanic stable). alpha(2025) = 0.653.

4. 2025-2050 (ASSUMED): held flat at alpha(2025) = 0.653 -- the neutral
   continuation for Case 2 (avoids inventing a decline path).

Writes data/assortativity_schedule.csv with columns:
year, alpha, quality, basis.
"""
import os
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "assortativity_schedule.csv")

# --- Pew newlywed interracial/interethnic share (% / 100) ---------------------
# Pew Research Center, "Interracial Marriage Grows," June 15, 2023
# (https://www.pewresearch.org/short-reads/2023/06/15/... chart table),
# plus Pew 2017 "Intermarriage in the U.S. 50 Years After Loving v. Virginia"
# (1967: 3%, 2015: 17%).
PEW_X = [
    (1960, 0.024), (1965, 0.029), (1967, 0.030), (1970, 0.040),
    (1975, 0.050), (1980, 0.067), (1985, 0.073), (1990, 0.083),
    (1995, 0.093), (2000, 0.112), (2005, 0.131), (2008, 0.145),
    (2009, 0.154), (2010, 0.151), (2011, 0.155), (2015, 0.170),
]

# Slope calibrated to the CPS-calibrated anchor alpha_CPS = 0.70 at 2013
# (Step 3, docs/STEP_LOG.md). x(2013) = 0.1625 by linear interpolation
# between 2011 and 2015.
K_SLOPE = (1.0 - 0.70) / 0.1625   # = 1.846153...

# 2015-2025 moderated trend for x (ASSUMED moderation): 0.18 pp/yr,
# half the 2008-2015 rate (Pew 2017: rise slowing).
X2020 = 0.170 + 5 * 0.0018
X2025 = 0.170 + 10 * 0.0018

# --- pre-1960 anchors: (year, alpha, quality, basis) --------------------------
PRE1960 = [
    (1650, 0.950, "assumed",
     "colonial: small homogeneous settler communities; no marriage data exist"),
    (1700, 0.950, "assumed", "colonial; held at 1650 value"),
    (1776, 0.940, "assumed", "independence; colonial endogamy assumed persistent"),
    (1800, 0.940, "assumed", "interpolated between 1776 and 1850 anchors"),
    (1850, 0.930, "estimated",
     "19th-c. European ethnic endogamy strong (Alba & Golden 1986)"),
    (1880, 0.930, "estimated",
     "60-70% of immigrants endogamous (Spoerlein et al. 2014); strong "
     "third-generation endogamy (Logan & Shin 2012)"),
    (1900, 0.935, "interpolated", "between the 1880 and 1910 anchors"),
    (1910, 0.940, "estimated",
     "endogamy 'castelike' for new ethnics ~1910 (Pagnini & Morgan 1990); "
     "first-generation exogamy low (Draschler 1921, NYC 1908-12)"),
    (1920, 0.940, "interpolated", "between the 1910 and 1940 anchors"),
    (1930, 0.940, "interpolated", "between the 1910 and 1940 anchors"),
    (1940, 0.945, "estimated",
     "pre-WWII; intermarriage rose 'dramatically after World War II' "
     "(Alba 1983), so pre-WWII endogamy still high"),
    (1950, 0.950, "interpolated",
     "between the 1940 anchor and the 1960 measured-series value"),
]


def main():
    rows = list(PRE1960)

    pew_years = np.array([y for y, _ in PEW_X])
    pew_x = np.array([x for _, x in PEW_X])

    def alpha_from_x(year, x, note):
        return (year, round(1.0 - K_SLOPE * x, 4),
                "measured-series",
                f"Pew newlywed intermarriage x={x*100:.1f}% -> "
                f"alpha = 1 - {K_SLOPE:.4f}*x ({note})")

    for y, x in PEW_X:
        note = ("Pew 2023 'Interracial Marriage Grows'" if y != 1967 and y != 2015
                else "Pew 2017 'Intermarriage ... 50 Years After Loving v. Virginia'")
        if y == 2015:
            note = ("Pew 2017 'Intermarriage ... 50 Years After Loving v. Virginia'; "
                    "mapping normalized so alpha(2013)=0.70 = alpha_CPS (Step 3)")
        rows.append(alpha_from_x(y, x, note))

    # consistency: the calibrated mapping reproduces the CPS anchor
    x2013 = float(np.interp(2013, pew_years, pew_x))
    a2013 = 1.0 - K_SLOPE * x2013
    rows.append((2013, round(a2013, 4), "measured-series",
                 f"x(2013)={x2013*100:.2f}% interpolated 2011-2015; "
                 f"alpha(2013)={a2013:.4f} = alpha_CPS (0.70) by construction"))

    for y, x in ((2020, X2020), (2025, X2025)):
        rows.append((y, round(1.0 - K_SLOPE * x, 4), "estimated",
                     f"x={x*100:.1f}%: moderated trend (0.18 pp/yr = half the "
                     f"2008-2015 rate; Pew 2017 notes the rise slowing)"))

    # Case-2 neutral continuation: held flat at the 2025 value (assumed)
    a2025 = 1.0 - K_SLOPE * X2025
    for y in (2030, 2040, 2050):
        rows.append((y, round(a2025, 4), "assumed",
                     f"Case-2 neutral continuation: held at alpha(2025)={a2025:.4f}"))

    df = pd.DataFrame(rows, columns=["year", "alpha", "quality", "basis"])
    df = df.sort_values("year").reset_index(drop=True)
    df.to_csv(OUT, index=False)

    # checks
    assert abs(a2013 - 0.70) < 1e-3, f"calibration consistency failed: {a2013}"
    assert df["alpha"].between(0, 1).all()
    print(f"k_slope = {K_SLOPE:.6f}")
    print(f"alpha(2013) = {a2013:.4f} (target 0.70)")
    print(f"alpha(1960) = {df.loc[df.year == 1960, 'alpha'].iloc[0]:.4f}")
    print(f"alpha(2015) = {df.loc[df.year == 2015, 'alpha'].iloc[0]:.4f}")
    print(f"alpha(2025) = {df.loc[df.year == 2025, 'alpha'].iloc[0]:.4f}")
    print("wrote", OUT, f"({len(df)} rows)")


if __name__ == "__main__":
    main()
