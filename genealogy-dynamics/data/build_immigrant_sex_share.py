#!/usr/bin/env python3
"""Build the immigrant male-share series sigma(t), 1650-2025.

sigma(t) = male fraction of the immigration inflow I(t). It replaces the
old sigma = 0.5 ASSUMED constant (genealogy/kernels.py:IMMIGRANT_MALE_FRACTION
kept as a legacy default only).

Construction:
  - sigma_LPR(t): male share of legal/permanent immigration inflow.
      * 1650-1819 : 0.65 ASSUMED. No annual series exists; colonial-era
        migration was heavily male-skewed indentured labor plus forced
        migration (the demographic-residual I(t) captures both). The
        1820-1917 measured 0.629 anchors the plausible range; the colonial
        value is set at 0.65, inside the 0.6-0.7 documented male skew.
      * 1820-1917 : 0.629 MEASURED. Historical Statistics of the United
        States (Millennial Ed., Carter et al., 2006): 62.9% of recorded
        1820-1917 immigrants were male.
      * 1918-1949 : linear 0.629 -> 0.50 ESTIMATED/INTERPOLATED. The quota
        era narrowed the male margin; HSUS documents the foreign-born
        sex ratio reaching ~1 by 1950.
      * 1950-2014 : linear 0.50 -> 0.458 ESTIMATED/INTERPOLATED, landing on
        the measured 2015 anchor.
      * FY2015 : 0.4581 MEASURED (DHS Yearbook of Immigration Statistics,
        Table 9: 481,485 male / 1,051,031).
      * FY2020 : 0.4615 MEASURED (Table 9: 326,414 / 707,362).
      * FY2021 : 0.4547 MEASURED (OHSS "U.S. Lawful Permanent Residents:
        2023": 336,450 / 740,002).
      * FY2022 : 0.4657 MEASURED (Table 9: 474,242 / 1,018,349).
      * FY2023 : 0.4523 MEASURED (Table 9: 530,550 / 1,172,910).
      * FY2024-2025 : 0.4523 INTERPOLATED (FY2024 Table 9 not extracted;
        2023 value held).
  - sigma_U(t): male share of the unauthorized inflow (1980+ only).
      * 2018 : 0.5195, 2019 : 0.5203, 2020 : 0.5243 MEASURED (OHSS
        "Estimates of the Unauthorized Immigrant Population ... Jan 2018
        to Jan 2022", Table 4).
      * 2021 : 0.5333 INTERPOLATED (2020->2022).
      * 2022 : 0.5423 MEASURED (OHSS Table 4); 2023-2025 held
        (INTERPOLATED).
      * 1980-2017 : 0.5215 ESTIMATED (mean of the measured 2018-2020
        values; the unauthorized stream is consistently male-skewed
        relative to legal immigration).
  - sigma(t) : flow-weighted blend
        sigma = (L*sigma_LPR + U*sigma_U) / (L + U),
    with L(t) = DHS LPR inflow (M/yr) and U(t) = unauthorized inflow
    (M/yr) from data/build_inputs.py (U = 0 before 1980, so sigma =
    sigma_LPR there; pre-1820 the residual I(t) takes sigma_LPR = 0.65).

Output: data/immigrant_sex_share.csv with per-year quality = weakest link
across the components actually entering the blend.
"""
import csv
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "data"))

from build_inputs import _load_lpr, _unauth_inflow, YEARS  # noqa: E402

OUT = os.path.join(BASE, "data", "immigrant_sex_share.csv")

# --- measured anchors -------------------------------------------------------
LPR_MEASURED = {  # fiscal year -> male share of LPR grants (DHS Yearbook Table 9)
    2015: 481485 / 1051031,   # 0.4581
    2020: 326414 / 707362,    # 0.4615
    2021: 336450 / 740002,    # 0.4547 (OHSS "U.S. Lawful Permanent Residents: 2023")
    2022: 474242 / 1018349,   # 0.4657
    2023: 530550 / 1172910,   # 0.4523
}
UNAUTH_MEASURED = {  # year -> male share (OHSS unauthorized estimates, Table 4)
    2018: 6010000 / 11570000,   # 0.5195
    2019: 5780000 / 11110000,   # 0.5203
    2020: 5510000 / 10510000,   # 0.5243
    2022: 5960000 / 10990000,   # 0.5423
}
SIGMA_COLONIAL = 0.65          # 1650-1819 ASSUMED (male-skewed indentured labor)
SIGMA_1820_1917 = 0.629        # MEASURED (HSUS Millennial Ed.: 62.9% male, 1820-1917)
SIGMA_1950 = 0.50              # ESTIMATED (HSUS: foreign-born sex ratio ~1 by 1950)
SIGMA_UNAUTH_EARLY = sum(UNAUTH_MEASURED[y] for y in (2018, 2019, 2020)) / 3.0


def _lin(a, b, t, t0, t1):
    return a + (b - a) * (t - t0) / (t1 - t0)


def sigma_lpr(y):
    """Male share of legal immigration inflow + quality label."""
    if y <= 1819:
        return SIGMA_COLONIAL, "assumed", \
            "no annual series; male-skewed indentured labor migration"
    if y <= 1917:
        return SIGMA_1820_1917, "measured", \
            "HSUS Millennial Ed.: 62.9% of 1820-1917 recorded immigrants male"
    if y <= 1949:
        return _lin(SIGMA_1820_1917, SIGMA_1950, y, 1917, 1950), \
            "estimated", "interpolated quota-era convergence to ~parity"
    if y <= 2014:
        if y in LPR_MEASURED:
            return LPR_MEASURED[y], "measured", "DHS Yearbook Table 9 (LPR by sex)"
        if y <= 1950:
            return SIGMA_1950, "estimated", \
                "HSUS: foreign-born sex ratio ~1 by 1950 (quota-era anchor)"
        v = _lin(SIGMA_1950, LPR_MEASURED[2015], y, 1950, 2015)
        q = "interpolated" if y >= 2016 else "estimated"
        return v, q, "between 1950 parity anchor and measured 2015 DHS anchor"
    if 2016 <= y <= 2019:
        return _lin(LPR_MEASURED[2015], LPR_MEASURED[2020], y, 2015, 2020), \
            "interpolated", "between measured 2015 and 2020 DHS Table 9 anchors"
    if y in LPR_MEASURED:
        return LPR_MEASURED[y], "measured", "DHS Yearbook Table 9 (LPR by sex)"
    return LPR_MEASURED[2023], "interpolated", "2023 measured value held"


def sigma_unauth(y):
    """Male share of unauthorized inflow + quality label (None before 1980)."""
    if y < 1980:
        return None, None, None
    if y in UNAUTH_MEASURED:
        return UNAUTH_MEASURED[y], "measured", \
            "OHSS unauthorized estimates 2018-2022, Table 4"
    if y == 2021:
        return _lin(UNAUTH_MEASURED[2020], UNAUTH_MEASURED[2022], 2021, 2020, 2022), \
            "interpolated", "2020->2022 OHSS Table 4"
    if y > 2022:
        return UNAUTH_MEASURED[2022], "interpolated", "2022 measured value held"
    return SIGMA_UNAUTH_EARLY, "estimated", \
        "mean of measured 2018-2020 OHSS values"


def main():
    lpr = _load_lpr()                      # fiscal year -> count
    # unauthorized inflow uses the same per-capita death rates as build_inputs
    import pandas as pd
    _inp = pd.read_csv(os.path.join(BASE, "data", "inputs.csv"))
    mu = {int(r.year): float(r.death_rate) for r in _inp.itertuples()}
    mu.update({y: mu[2024] for y in YEARS + [2025] if y not in mu})
    unauth = _unauth_inflow(mu)
    unauth[2025] = unauth[2024]            # 2025 holds 2024 (Case-1 convention)

    order = {"assumed": 0, "estimated": 1, "interpolated": 2, "measured": 3}
    rows = []
    for y in YEARS + [2025]:
        s_lpr, q_lpr, n_lpr = sigma_lpr(y)
        L = lpr.get(y, 0.0) / 1e6
        s_u, q_u, n_u = sigma_unauth(y)
        U = float(unauth.get(y, 0.0))
        if y == 2025:
            # LPR series ends 2024: hold the 2024 blend (INTERPOLATED)
            sig, q, note = rows[-1][1], "interpolated", "2024 blend held"
        elif y < 1980 or U <= 0.0:
            sig, q, note = s_lpr, q_lpr, f"LPR-component: {n_lpr} [{q_lpr}]"
        else:
            sig = (L * s_lpr + U * s_u) / (L + U)
            q = q_lpr if order[q_lpr] < order[q_u] else q_u
            note = (f"LPR {L:.3f}M x {s_lpr:.4f} [{q_lpr}] + unauth "
                    f"{U:.3f}M x {s_u:.4f} [{q_u}]")
        rows.append((y, sig, s_lpr, s_u if s_u else 0.0, L, U, q, note))
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["year", "sigma", "sigma_lpr", "sigma_unauth",
                    "lpr_M", "unauth_inflow_M", "quality", "quality_note"])
        for y, sig, s_lpr, s_u, L, U, q, note in rows:
            w.writerow([y, f"{sig:.6f}", f"{s_lpr:.6f}", f"{s_u:.6f}",
                        f"{L:.6f}", f"{U:.6f}", q, note])
    print(f"wrote {OUT} ({len(rows)} rows)")
    print(f"  sigma(1650)={rows[0][1]:.4f}, sigma(1900)={rows[250][1]:.4f}, "
          f"sigma(1950)={rows[300][1]:.4f}, sigma(2000)={rows[350][1]:.4f}, "
          f"sigma(2024)={rows[-2][1]:.4f}")


if __name__ == "__main__":
    main()
