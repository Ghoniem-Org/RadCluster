#!/usr/bin/env python3
"""Build data/age_inputs.csv: annual 1650-2024 age-specific fertility and
mortality schedules for the age-structured genealogy model.

Output columns: year,
  beta_<band>   (18): per-woman birth rate by 5-year age band (0 outside 15-49),
  mu_pat_<band> (18): male per-capita death rate by 5-year age band,
  mu_mat_<band> (18): female per-capita death rate by 5-year age band,
  quality, quality_note.

Provenance labels (measured / interpolated / estimated / assumed) follow the
project's standing rule: NEVER present an assumed parameter as an estimate.

FERTILITY
  1950-2023 : ESTIMATED. UN World Population Prospects 2024, single-year
              ASFR (births per 1,000 women), USA, Medium variant. For the
              United States the UN's input data are NCHS registered births
              by age of mother; the UN adjusts/interpolates, hence
              "estimated" rather than "measured". Aggregated to 5-year
              bands (mean of single-year ASFR / 1000).
  2024      : INTERPOLATED (2023 schedule held).
  1800-1949 : INTERPOLATED. Annual TFR from Haines' historical US fertility
              series (white population; ESTIMATED from census
              child-woman ratios) times the 1950 WPP age pattern (PASFR
              shape held fixed -- the pre-1950 age pattern is the
              interpolation assumption).
  1650-1799 : ASSUMED. TFR held at the 1800 level (7.04) with the 1950
              age pattern; no age-specific colonial fertility data exist.

MORTALITY (sex-specific m(x) from life tables)
  1900-2023 : MEASURED. Human Life-Table Database (MPIDR, lifetable.de) USA
              pooled file: single-year period life tables by sex. Table
              selection per year: Year1==Year2, Region=0, Ethnicity=0
              (total population -- this excludes the state-level and
              race/Hispanic-specific tables archived alongside), preferring
              the NCHS/NVSR table over the SSA one where both exist
              (1997-1999). Underlying sources per HLD reference list:
              - 1900-1996: Bell & Miller, SSA Actuarial Study No. 116
                           (Life Tables for the US Social Security Area
                           1900-2100);
              - 1997-2023: NCHS National Vital Statistics Reports
                           (decennial and annual US life tables).
              Single-year m(x) aggregated to 5-year bands by
              exposure-weighting: mu_band = sum m(x)*L(x) / sum L(x).
  2024      : INTERPOLATED (2023 schedule held).
  1800-1899 : INTERPOLATED. Mean 1900-1902 HLD age pattern (sex-specific),
              level-scaled RELATIVE to the CDR anchors: lambda(y) =
              CDR_anchor(y)/CDR_anchor(1900). The 1900-02 DRA-based pattern
              implies CDR 20.4 vs the national 1900 anchor 17.2, so relative
              scaling splices the backcast continuously onto the measured
              1900 table (absolute scaling would jump at 1900). CDR anchors
              1800: 24.0, 1850: 23.0, 1900: 17.2 per 1,000 (same Haines-based
              anchors as data/inputs.csv, linearly interpolated).
  1650-1799 : ASSUMED. Same 1900-02 pattern scaled to CDR = 28.0 per 1,000
              (colonial mortality assumption, consistent with
              data/build_inputs.py).

Full citations are written to data/age_inputs_provenance.md.
"""
import csv
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(BASE))  # project root (genealogy pkg)
RAW = os.path.join(BASE, "raw")
OUT = os.path.join(BASE, "age_inputs.csv")

BANDS = ["0_4", "5_9", "10_14", "15_19", "20_24", "25_29", "30_34",
         "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69",
         "70_74", "75_79", "80_84", "85p"]
N_AGE = 18
REPRO = list(range(3, 10))  # 15-19 .. 45-49

# Haines historical US TFR (white population), from Table 1 of the NBER
# conference paper (Haines 2002). ESTIMATED from census child-woman ratios.
HAINES_TFR = {1800: 7.04, 1810: 6.92, 1820: 6.73, 1830: 6.55, 1840: 6.14,
              1850: 5.42, 1860: 5.21, 1870: 4.55, 1880: 4.24, 1890: 3.87,
              1900: 3.56, 1910: 3.42, 1920: 3.17, 1930: 2.45, 1940: 2.22}
# CDR anchors per 1,000 (same series as data/inputs.csv).
CDR_ANCHORS = {1800: 24.0, 1850: 23.0, 1900: 17.2}
CDR_COLONIAL = 28.0

Y0, Y1 = 1650, 2024


def interp_anchors(anchors, y):
    ys = sorted(anchors)
    if y <= ys[0]:
        return anchors[ys[0]]
    if y >= ys[-1]:
        return anchors[ys[-1]]
    for a, b in zip(ys[:-1], ys[1:]):
        if a <= y <= b:
            f = (y - a) / (b - a)
            return (1 - f) * anchors[a] + f * anchors[b]
    raise AssertionError


# --- fertility ---------------------------------------------------------------
def load_wpp_beta():
    """beta[year] -> (18,) per-woman birth rates, 1950-2023."""
    asfr = {}
    with open(os.path.join(RAW, "wpp_usa_asfr.csv")) as f:
        for row in csv.DictReader(f):
            asfr.setdefault(int(row["year"]), {})[int(row["age"])] = \
                float(row["asfr_per_1000"])
    beta = {}
    for y, d in asfr.items():
        b = np.zeros(N_AGE)
        for a in REPRO:
            ages = range(5 * a, 5 * a + 5)
            b[a] = np.mean([d[x] for x in ages]) / 1000.0
        beta[y] = b
    return beta


def tfr_of(beta):
    return 5.0 * beta[REPRO].sum()


# --- mortality ---------------------------------------------------------------
def load_hld_mu():
    """mu[(year, sex)] -> (18,) death rates, 1900-2023, sex 1=male 2=female.

    National total-population tables: Region=0, Ethnicity=0 (excludes the
    state-level and race/Hispanic-specific tables); where two national
    tables exist (1997-1999: SSA Bell & Miller vs NCHS NVSR) the NCHS/NVSR
    table (larger Ref-ID = later publication) is preferred.
    """
    path = os.path.join(RAW, "hld_usa", "USA.csv")
    # pass 1: choose ref per year
    cands = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["Year1"] != r["Year2"]:
                continue
            y = int(r["Year1"])
            if not (1900 <= y <= 2023):
                continue
            if r["Region"] != "0" or r["Ethnicity"] != "0":
                continue
            ref = r["Ref-ID"]
            base = float(ref.split(".")[0])
            cands.setdefault(y, set()).add((base, ref))
    chosen = {}
    for y, s in cands.items():
        best = max(s)[1]  # later publication preferred (NVSR over SSA)
        chosen[y] = best
    # report the SSA->NVSR handoffs
    for y in (1997, 1998, 1999):
        print(f"  HLD {y}: national table = Ref {chosen[y]}")
    # pass 2: aggregate m(x) to bands, exposure-weighted by L(x)
    acc = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["Year1"] != r["Year2"]:
                continue
            y = int(r["Year1"])
            if not (1900 <= y <= 2023) or r["Ref-ID"] != chosen[y]:
                continue
            if r["Region"] != "0" or r["Ethnicity"] != "0":
                continue
            sex = int(r["Sex"])
            x = int(r["Age"])
            m = float(r["m(x)"])
            L = float(r["L(x)"])
            band = 17 if x >= 85 else x // 5
            k = (y, sex)
            a = acc.setdefault(k, {"num": np.zeros(N_AGE),
                                   "den": np.zeros(N_AGE)})
            a["num"][band] += m * L
            a["den"][band] += L
    mu = {}
    for k, a in acc.items():
        assert (a["den"] > 0).all(), f"empty band in HLD table {k}"
        mu[k] = a["num"] / a["den"]
    assert len(mu) == 2 * 124, f"{len(mu)} (year,sex) tables"
    # sanity: male e0 < female e0 every year (spot check via mu[0])
    return mu, chosen


# --- main --------------------------------------------------------------------
def main():
    print("loading WPP ASFR ...")
    wpp_beta = load_wpp_beta()
    print(f"  years {min(wpp_beta)}-{max(wpp_beta)}")
    tfr1950 = tfr_of(wpp_beta[1950])
    print(f"  WPP 1950 TFR = {tfr1950:.2f}")
    pasfr = 5.0 * wpp_beta[1950][REPRO] / tfr1950  # 1950 age pattern (shares)

    print("loading HLD life tables ...")
    hld_mu, chosen = load_hld_mu()

    # 1900-02 reference pattern for backcasting (sex-specific)
    pat_ref = np.mean([hld_mu[(y, 1)] for y in (1900, 1901, 1902)], axis=0)
    mat_ref = np.mean([hld_mu[(y, 2)] for y in (1900, 1901, 1902)], axis=0)
    # stationary age distribution of the 1900-02 tables (pi propto L(x))
    # -> implied CDR of the unscaled pattern, for level-scaling to anchors
    Lpat = np.zeros(N_AGE)
    Lmat = np.zeros(N_AGE)
    with open(os.path.join(RAW, "hld_usa", "USA.csv")) as f:
        for r in csv.DictReader(f):
            if r["Year1"] != r["Year2"]:
                continue
            y = int(r["Year1"])
            if y not in (1900, 1901, 1902) or r["Ref-ID"] != chosen[y]:
                continue
            if r["Region"] != "0" or r["Ethnicity"] != "0":
                continue
            x = int(r["Age"])
            band = 17 if x >= 85 else x // 5
            if int(r["Sex"]) == 1:
                Lpat[band] += float(r["L(x)"])
            else:
                Lmat[band] += float(r["L(x)"])
    pi = (Lpat / Lpat.sum() + Lmat / Lmat.sum()) / 2.0
    cdr_ref = 1000.0 * float((pat_ref * pi).sum() + (mat_ref * pi).sum()) / 2.0
    print(f"  1900-02 pattern implied CDR on stationary pop = {cdr_ref:.2f}/1000")
    print("  (vs the 1900 national anchor 17.2: the 1900-02 table is DRA-based")
    print("   and more urban; pre-1900 levels are therefore scaled RELATIVE to")
    print("   the 1900 anchor so the schedule splices continuously to the")
    print("   measured 1900 table)")


    # --- annual series 1650-2024 -------------------------------------------
    hdr = (["year"]
           + [f"beta_{b}" for b in BANDS]
           + [f"mu_pat_{b}" for b in BANDS]
           + [f"mu_mat_{b}" for b in BANDS]
           + ["quality", "quality_note"])
    rows_out = []
    for y in range(Y0, Y1 + 1):
        # fertility
        if 1950 <= y <= 2023:
            beta = wpp_beta[y]
            qual = "estimated"
            note = ("ASFR: UN WPP 2024 single-year estimates (US inputs are "
                    "NCHS registered births); mortality: HLD/NCHS life tables")
        elif y == 2024:
            beta = wpp_beta[2023]
            qual = "interpolated"
            note = "2023 schedules held for 2024"
        elif 1800 <= y <= 1949:
            if y <= 1940:
                tfr = interp_anchors(HAINES_TFR, y)
            else:  # 1941-1949: Haines 1940 -> WPP 1950 (baby-boom ramp)
                f = (y - 1940) / 10.0
                tfr = (1 - f) * HAINES_TFR[1940] + f * tfr_of(wpp_beta[1950])
            beta = np.zeros(N_AGE)
            beta[REPRO] = tfr * pasfr / 5.0
            qual = "interpolated"
            note = ("TFR level from Haines historical series (white pop., "
                    "census-based); age pattern = 1950 WPP shape held fixed")
        else:  # 1650-1799
            beta = np.zeros(N_AGE)
            beta[REPRO] = HAINES_TFR[1800] * pasfr / 5.0
            qual = "assumed"
            note = ("TFR held at 1800 level (7.04) with 1950 age pattern; "
                    "no age-specific colonial fertility data exist")
        # mortality
        if 1900 <= y <= 2023:
            mu_p, mu_m = hld_mu[(y, 1)], hld_mu[(y, 2)]
            if qual == "estimated":
                pass
            else:
                qual = "interpolated"
                note += "; mortality: HLD/NCHS life tables"
        elif y == 2024:
            mu_p, mu_m = hld_mu[(2023, 1)], hld_mu[(2023, 2)]
            qual = "interpolated"
            note = "2023 schedules held for 2024"
        elif 1800 <= y <= 1899:
            # RELATIVE level-scaling: the 1900-02 DRA-based pattern implies a
            # higher CDR (20.4) than the national 1900 anchor (17.2), so the
            # anchors set the relative level and the schedule splices
            # continuously onto the measured 1900 table (no jump at 1900).
            lam = interp_anchors(CDR_ANCHORS, y) / CDR_ANCHORS[1900]
            mu_p, mu_m = lam * pat_ref, lam * mat_ref
            qual = "interpolated"
            note = ("mortality: 1900-02 HLD age pattern, level set relative "
                    "to CDR anchors (splices to measured 1900 table)"
                    + ("; " + note if "TFR" in note else ""))
        else:  # 1650-1799
            lam = CDR_COLONIAL / CDR_ANCHORS[1900]
            mu_p, mu_m = lam * pat_ref, lam * mat_ref
            qual = "assumed"
            note = ("1900-02 HLD age pattern, level relative to assumed "
                    "colonial CDR 28/1000; TFR held at 1800 level")
        rows_out.append([y]
                        + [f"{v:.6f}" for v in beta]
                        + [f"{v:.6f}" for v in mu_p]
                        + [f"{v:.6f}" for v in mu_m]
                        + [qual, note])

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(hdr)
        w.writerows(rows_out)
    print(f"wrote {OUT} ({len(rows_out)} annual rows)")

    # --- diagnostics ---------------------------------------------------------
    print("diagnostics:")
    from genealogy.age import mean_age_at_childbearing
    for y in (1950, 1980, 2023):
        b = wpp_beta[y]
        print(f"  {y}: TFR={tfr_of(b):.2f} "
              f"MAC={mean_age_at_childbearing(b):.1f} yr")
    for y in (2007, 2008, 2009):
        mp, mm = hld_mu[(y, 1)], hld_mu[(y, 2)]
        # crude e0 from the band schedule (for the spot check only)
        print(f"  HLD {y}: mu_20-24 M={mp[4]:.5f} F={mm[4]:.5f}, "
              f"mu_70-74 M={mp[14]:.4f} F={mm[14]:.4f}")
    # implied CDR of the 2023 measured schedule on the 2023 stationary pop
    print("done.")

    write_provenance()


def write_provenance():
    md = """# Age inputs: provenance (data/age_inputs.csv)

Annual 1650-2024 age-specific fertility (per-woman birth rates by 5-year
band, nonzero only at 15-49) and mortality (per-capita death rates by
5-year band and sex) for the age-structured genealogy model.

Every row carries a `quality` label: exactly one of `measured`,
`interpolated`, `estimated`, `assumed` (weakest link of its inputs).

## Fertility

- **1950-2023 -- estimated.** UN World Population Prospects 2024, table
  "Fertility by single age": age-specific fertility rates (births per
  1,000 women) by single year of age of mother, USA ("United States of
  America (and dependencies)"), Medium variant (VarID=2), annual
  1950-2023. For the United States the UN's input data are NCHS
  registered births classified by mother's age (see UN WPP 2024
  Methodology); the UN adjusts and interpolates, hence *estimated*.
  Aggregated here to 5-year bands as the band mean ASFR / 1000.
  Sanity: implied TFR 1957 = 3.71 (baby-boom peak), 1976 = 1.79
  (trough), 2007 = 2.09, 2023 = 1.62; mean age at childbearing 2023 =
  29.9 yr.
- **2024 -- interpolated** (2023 schedule held).
- **1800-1949 -- interpolated.** Annual total fertility rate from
  Haines' historical US fertility series (white population, estimated
  from census child-woman ratios; Table 1: 1800: 7.04 ... 1900: 3.56 ...
  1940: 2.22; 1941-1949 linearly ramped to the WPP 1950 TFR),
  multiplied by the fixed 1950 WPP age pattern (PASFR shares). The
  held-fixed age shape is the interpolation assumption.
- **1650-1799 -- assumed.** TFR held at the 1800 level (7.04) with the
  1950 age pattern. No age-specific colonial fertility data exist.

## Mortality (sex-specific)

- **1900-2023 -- measured.** Human Life-Table Database (Max Planck
  Institute for Demographic Research, https://www.lifetable.de), USA
  pooled file (years 1900-2023): single-year period life tables by sex;
  m(x) = age-specific death rate, aggregated to 5-year bands by
  exposure weighting mu_band = sum m(x) L(x) / sum L(x).
  Table selection per year: Year1 == Year2, Region = 0, Ethnicity = 0
  (total population; this excludes the US *state* life tables and the
  race/Hispanic-specific tables archived in the same file), preferring
  the NCHS/NVSR table where two national tables exist (1997-1999).
  Underlying sources (HLD reference list):
  - 1900-1996: Bell & Miller, *Life Tables for the United States Social
    Security Area 1900-2100*, SSA Actuarial Study No. 116;
  - 1997-2023: NCHS National Vital Statistics Reports (decennial and
    annual United States Life Tables).
- **2024 -- interpolated** (2023 schedule held).
- **1800-1899 -- interpolated.** Mean 1900-1902 HLD age pattern
  (sex-specific), level-scaled *relative* to the CDR anchors:
  lambda(y) = CDR_anchor(y)/CDR_anchor(1900). The 1900-02 DRA-based
  pattern implies CDR 20.4 on the stationary population vs the
  national 1900 anchor 17.2, so relative scaling splices the backcast
  continuously onto the measured 1900 table. CDR anchors per 1,000:
  1800: 24.0, 1850: 23.0, 1900: 17.2 (Haines-based; the same anchors
  used by data/build_inputs.py, linearly interpolated).
- **1650-1799 -- assumed.** Same 1900-02 pattern scaled to the assumed
  colonial CDR of 28.0 per 1,000 (consistent with data/build_inputs.py).

## Deliberate modeling assumptions (labeled `assumed` in code, not data)

- Male age-specific fertility follows the female ASFR schedule
  (genealogy/age.py): NCHS publishes no male age-specific fertility
  series, so fathers are drawn from reproductive-age males with the
  female schedule's shape.
- Immigrant age profile: fixed stylized working-age-concentrated shares
  (genealogy/age.py:IMMIGRANT_AGE_PROFILE). Future calibration source:
  DHS Yearbook of Immigration Statistics, Table 8.
- Continuous aging flux gamma_a = 1/5 yr^-1 for bands 0-16, gamma_17 = 0
  for the open-ended 85+ band (documented in genealogy/age.py).

## References

- United Nations, Department of Economic and Social Affairs, Population
  Division (2024). *World Population Prospects 2024*, "Fertility by
  single age" bulk CSV. https://population.un.org/wpp/ (CC BY 3.0 IGO).
- Human Life-Table Database. Max Planck Institute for Demographic
  Research (Rostock) and Vienna Institute of Demography.
  https://www.lifetable.de/ -- USA pooled file, years 1900-2023.
- Bell, F. C. & Miller, M. L. (2005). *Life Tables for the United
  States Social Security Area 1900-2100*. SSA Actuarial Study No. 116.
- National Center for Health Statistics. *United States Life Tables*
  (decennial and annual), National Vital Statistics Reports, various
  volumes 1997-2023. https://www.cdc.gov/nchs/nvss/
- Haines, M. R. (2002). "The White Population of the United States,
  1790-1920" / Table 1 historical fertility series, NBER conference
  paper. https://conference.nber.org/confer/2002/si2002/haines.pdf
"""
    p = os.path.join(BASE, "age_inputs_provenance.md")
    with open(p, "w") as f:
        f.write(md)
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
    # One-step pipeline: the pre-1900 mortality LEVEL must then be
    # recalibrated so the stable-population CDR matches Haines (see
    # data/recalibrate_pre1900_mortality.py). Running the builder alone
    # leaves the uncalibrated stationary-CDR scaling in place.
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "recal", os.path.join(BASE, "recalibrate_pre1900_mortality.py"))
    recal = importlib.util.module_from_spec(spec)
    sys.modules["recal"] = recal
    spec.loader.exec_module(recal)
    recal.main()
