#!/usr/bin/env python3
"""Build clean model inputs from data/raw/ — full span 1650-2024.

Outputs (in data/):
- inputs.csv: year, immigration_M, birth_rate, death_rate, tfr,
  cbr_per_1000, cdr_per_1000, quality, quality_note
- validation.csv: (unchanged) year, total_pop_M (World Bank 1980-2024),
  migrant_stock_M, fb_census_M, source
- validation_long.csv: year, pop_M, source (population anchors 1650-2024:
  colonial estimates, decennial census, World Bank annual)
- initial_condition_1980.csv: (unchanged)

Rate conventions (documented):
- immigration_M (millions/yr), by era:
  * 1650-1819: DEMOGRAPHIC RESIDUAL. For each decade between population
    anchors, migration = anchor(t1) - anchor(t0)*exp(int(b-mu-eps)).
    This captures all net migration — including enslaved imports, which
    the DHS series misses by construction. quality=interpolated.
  * 1820-1979: DHS fiscal-year LPR count / 1e6, mapped to calendar year t.
    (Fiscal years run Oct-Sep; treating FY t as year t is a <=6-month shift.)
    quality=measured.
  * 1980-2024: DHS LPR + unauthorized inflow from DHS OHSS stock estimates:
        I_unauth(t) = max(0, S(t)-S(t-1)) + mu(t)*S(t-1),
    i.e. annualized positive net stock growth plus deaths of the
    unauthorized stock (mu = death_rate). This is NET stock growth turned
    into an approximate GROSS inflow; departures/legalizations of the
    unauthorized stock are netted, not added back. quality=interpolated
    (the stock->flow conversion and interpolation dominate the error).
    CAVEAT: DHS LPR counts include adjustments of status, so some
    legalized unauthorized immigrants are counted in BOTH series; the
    2000-2024 inflow is a modest overcount of true gross entries.
- birth_rate (per-capita /yr, uniform across classes):
      beta(t) = 2 * CBR(t)/1000
  The model's births are B = 0.5 * sum_p beta_p * P_p (the 1/2 corrects for
  two parents per birth), so a uniform beta gives B = CBR/1000 * total,
  i.e. total births exactly track the crude birth rate by construction.
- death_rate (per-capita /yr, uniform): mu(t) = CDR(t)/1000.
- Vital-rate anchors, pre-1960 (per 1000):
  * CBR: Haines white crude birth rates, EH.net Table 1 (1800: 55.0,
    1810: 54.3, 1820: 52.8, 1830: 51.4, 1840: 48.3, 1850: 43.3,
    1860: 41.4, 1870: 38.3, 1880: 35.2, 1890: 31.5, 1900: 30.1,
    1910: 29.2, 1920: 26.9, 1930: 20.6, 1940: 18.6, 1950: 23.0);
    black CBRs ran higher but the white series is the standard continuous
    one — see README.
  * CDR: 1800 ~24 (Blodget/Haines national range 24-26); 1850 ~23;
    1900 17.2 (Death Registration Area); 1930 11.3 (DRA). Linear between.
  * 1650-1799: CBR=52, CDR=28 ASSUMED (documented; colonial demography).
    quality=assumed for 1650-1799 rows.
- Emigration is NOT in inputs.csv: kept as a stylized constant in the
  hindcast driver (calibrated, see data/emigration_calibration.txt), since
  no reliable annual US emigration series exists. The colonial residual
  immigration uses the same calibrated eps for consistency.
- quality per row = weakest link across (immigration, vitals):
  1650-1799 assumed; 1800-1819 interpolated; 1820-1959 interpolated;
  1960-1979 measured; 1980-2024 interpolated (unauthorized series).
"""
import argparse
import csv
import json
import math
import os

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "raw")

YEARS = list(range(1650, 2025))

# --- Population anchors (people) -------------------------------------------
# Colonial estimates: U.S. Census Bureau, Historical Statistics of the United
# States, Colonial Times to 1970, Series Z 1-19.
COLONIAL_POP = {
    1650: 50368, 1660: 75058, 1670: 111935, 1680: 151507, 1690: 210372,
    1700: 250888, 1710: 331711, 1720: 466185, 1730: 629445, 1740: 905563,
    1750: 1170760, 1760: 1593625, 1770: 2148076, 1780: 2780369,
}
# Decennial census totals: U.S. Census Bureau (1790-2020).
DECENNIAL_POP = {
    1790: 3929214, 1800: 5308483, 1810: 7239881, 1820: 9638453,
    1830: 12866020, 1840: 17069453, 1850: 23191876, 1860: 31443321,
    1870: 38558371, 1880: 50189209, 1890: 62979766, 1900: 76212168,
    1910: 92228496, 1920: 106021537, 1930: 123202624, 1940: 132164569,
    1950: 151325798, 1960: 179323175, 1970: 203211926, 1980: 226656805,
    1990: 248709873, 2000: 281421906, 2010: 308745538, 2020: 331449281,
}

# --- Vital-rate anchors (per 1000 per year) --------------------------------
CBR_ANCHORS = {1800: 55.0, 1810: 54.3, 1820: 52.8, 1830: 51.4, 1840: 48.3,
               1850: 43.3, 1860: 41.4, 1870: 38.3, 1880: 35.2, 1890: 31.5,
               1900: 30.1, 1910: 29.2, 1920: 26.9, 1930: 20.6, 1940: 18.6,
               1950: 23.0}
CBR_COLONIAL = 52.0   # assumed, 1650-1799
CDR_ANCHORS = {1800: 24.0, 1850: 23.0, 1900: 17.2, 1930: 11.3}
CDR_COLONIAL = 28.0   # assumed, 1650-1799

# --- Unauthorized immigrant stock, Jan-1 (millions) -------------------------
# DHS Office of Homeland Security Statistics lineage. PDFs cached in
# data/raw/dhs_unauth/ (see fetch_data.py for URLs).
UNAUTH_STOCK = {
    1980: (2.00, "assumed (INS-era; DHS 1990-2000 report starts at 1990=3.5M)"),
    1990: (3.50, "DHS 1990-2000 report"),
    1996: (5.00, "DHS 1990-2000 report (Oct 1996)"),
    2000: (8.46, "Baker & Warren 2024, Table A2-1"),
    2005: (10.49, "Baker & Warren 2024, Table A2-1"),
    2006: (11.31, "Baker & Warren 2024, Table A2-1"),
    2007: (11.78, "Baker & Warren 2024, Table A2-1"),
    2008: (11.60, "Baker & Warren 2024, Table A2-1"),
    2009: (10.75, "Baker & Warren 2024, Table A2-1"),
    2010: (11.59, "Baker & Warren 2024, Table A2-1 (2010* revised)"),
    2011: (11.51, "Baker & Warren 2024, Table A2-1"),
    2012: (11.43, "Baker & Warren 2024, Table A2-1"),
    2013: (11.21, "Baker & Warren 2024, Table A2-1"),
    2014: (11.46, "Baker & Warren 2024, Table A2-1"),
    2015: (11.44, "Baker & Warren 2024, Table A2-1 (2015** revised)"),
    2016: (11.75, "Baker & Warren 2024, Table A2-1"),
    2017: (11.41, "Baker & Warren 2024, Table A2-1"),
    2018: (11.57, "Baker & Warren 2024, Table A2-1 (2018** revised)"),
    2019: (11.11, "Baker & Warren 2024, Table A2-1"),
    2020: (10.51, "Baker & Warren 2024, Table A2-1"),
    2022: (10.99, "Baker & Warren 2024, Table A2-1"),
    2023: (11.23, "assumed: linear extrapolation of 2020-2022 DHS trend"),
    2024: (11.47, "assumed: linear extrapolation of 2020-2022 DHS trend"),
}

# Decennial census foreign-born counts (millions), published figures
# (U.S. Census Bureau, Historical Census Statistics on the Foreign-Born
# Population of the United States, 1850-2000, Working Paper 81).
FB_CENSUS = {1980: 14.08, 1990: 19.77, 2000: 31.11}

CALIB_FILE = os.path.join(BASE, "emigration_calibration.txt")


def _load_lpr():
    d = {}
    with open(os.path.join(RAW, "dhs_lpr_fy.csv")) as f:
        for row in csv.DictReader(f):
            d[int(row["fiscal_year"])] = int(row["lpr"])
    return d


def _load_wb(indicator):
    d = {}
    payload = json.load(open(os.path.join(RAW, f"wb_{indicator}.json")))
    for r in payload[1]:
        if r["value"] is not None:
            d[int(r["date"])] = float(r["value"])
    return d


def _interp_anchors(anchors, years):
    """Linear interpolation of an anchor dict over a year list."""
    out = {}
    keys = sorted(anchors)
    for y in years:
        if y <= keys[0]:
            out[y] = anchors[keys[0]]
        elif y >= keys[-1]:
            out[y] = anchors[keys[-1]]
        else:
            lo = max(k for k in keys if k <= y)
            hi = min(k for k in keys if k >= y)
            if lo == hi:
                out[y] = anchors[lo]
            else:
                f = (y - lo) / (hi - lo)
                out[y] = anchors[lo] + f * (anchors[hi] - anchors[lo])
    return out


def _unauth_stock_series():
    """Annual Jan-1 unauthorized stock (millions), interpolated."""
    years = list(range(1980, 2025))
    anchors = {y: v for y, (v, _) in UNAUTH_STOCK.items()}
    return _interp_anchors(anchors, years)


def _unauth_inflow(mu):
    """Annual unauthorized inflow (millions/yr), 1980-2024.

    I(t) = max(0, S(t)-S(t-1)) + mu(t)*S(t-1); S(1979) := S(1980).
    mu: dict year -> per-capita death rate.
    """
    S = _unauth_stock_series()
    out = {}
    for t in range(1980, 2025):
        prev = S[t - 1] if t - 1 in S else S[1980]
        out[t] = max(0.0, S[t] - prev) + mu[t] * prev
    return out


def _colonial_residual_immigration(cbr, cdr, eps):
    """Demographic-residual immigration, 1650-1819 (millions/yr).

    For each decade [t0, t1) between population anchors, solve for the
    constant annual inflow I such that the continuous model
        dP/dt = I + r(t)*P,  r(t) = b(t)-d(t)-eps
    carries P(t0)=anchor(t0) to P(t1)=anchor(t1):

        I = max(anchor(t1) - anchor(t0)*exp(G), 0) / S,

    where G is the decade's cumulative natural growth and S integrates the
    within-decade growth of immigrants arriving uniformly through each year
    (r taken piecewise-constant at the Jan-1 annual values, matching the
    model's time-dependent integration convention).  Dividing by S (rather
    than by the decade length) corrects the first-order error of ignoring
    immigrants' own natural increase inside the decade.
    """
    anchors = {y: v / 1e6 for y, v in COLONIAL_POP.items()}
    anchors.update({y: v / 1e6 for y, v in DECENNIAL_POP.items() if y <= 1820})
    out = {}
    for t0, t1 in zip(sorted(anchors), sorted(anchors)[1:]):
        n = t1 - t0
        r = [cbr[t] / 1000.0 - cdr[t] / 1000.0 - eps for t in range(t0, t1)]
        G = sum(r)
        # S = sum over arrival years k of (growth from year k to t1) *
        #     (within-year arrival averaging factor).
        S = 0.0
        for k in range(n):
            grow = math.exp(sum(r[k:]))
            avg = (1.0 - math.exp(-r[k])) / r[k] if abs(r[k]) > 1e-12 else 1.0
            S += grow * avg
        M = anchors[t1] - anchors[t0] * math.exp(G)
        annual = max(M, 0.0) / S
        for t in range(t0, t1):
            out[t] = annual
    return out


def _resolve_eps(cli_eps):
    if cli_eps is not None:
        return cli_eps, "CLI --eps"
    if os.path.exists(CALIB_FILE):
        with open(CALIB_FILE) as f:
            for line in f:
                if line.startswith("eps_calibrated"):
                    v = float(line.split("=")[1])
                    return v, f"calibration file ({CALIB_FILE})"
    print("NOTE: no calibration file yet; using eps=0.0002 placeholder "
          "for the colonial residual. Re-run after calibrate_emigration.py.")
    return 0.0002, "placeholder"


def build_inputs(eps):
    lpr = _load_lpr()
    wb_tfr = _load_wb("SP.DYN.TFRT.IN")
    wb_cbr = _load_wb("SP.DYN.CBRT.IN")
    wb_cdr = _load_wb("SP.DYN.CDRT.IN")

    # Vital rates: colonial assumed -> Haines anchors -> World Bank.
    cbr = {}
    cdr = {}
    for y in YEARS:
        if y < 1800:
            cbr[y] = CBR_COLONIAL
            cdr[y] = CDR_COLONIAL
        elif y < 1960:
            cbr[y] = _interp_anchors({**CBR_ANCHORS, 1960: wb_cbr[1960]}, [y])[y]
            cdr[y] = _interp_anchors({**CDR_ANCHORS, 1960: wb_cdr[1960]}, [y])[y]
        else:
            cbr[y] = wb_cbr[y]
            cdr[y] = wb_cdr[y]
    mu = {y: cdr[y] / 1000.0 for y in YEARS}

    # Immigration by era.
    imm = {}
    imm_note = {}
    residual = _colonial_residual_immigration(cbr, cdr, eps)
    unauth = _unauth_inflow(mu)
    for y in YEARS:
        if y <= 1819:
            imm[y] = residual[y]
            imm_note[y] = "demographic residual vs population anchors"
        elif y <= 1979:
            imm[y] = lpr[y] / 1e6
            imm_note[y] = "DHS LPR (FY as calendar year)"
        else:
            imm[y] = lpr[y] / 1e6 + unauth[y]
            imm_note[y] = "DHS LPR + unauthorized stock->flow"

    def quality(y):
        if y <= 1799:
            imm_q = "interpolated"
            vit_q = "assumed"
        elif y <= 1819:
            imm_q = "interpolated"
            vit_q = "interpolated"
        elif y <= 1959:
            imm_q = "measured"
            vit_q = "interpolated"
        elif y <= 1979:
            imm_q = "measured"
            vit_q = "measured"
        else:
            imm_q = "interpolated"
            vit_q = "measured"
        order = {"assumed": 0, "interpolated": 1, "measured": 2}
        q = imm_q if order[imm_q] < order[vit_q] else vit_q
        note = (f"immigration: {imm_note[y]} [{imm_q}]; "
                f"vitals: {'CBR=%g/CDR=%g assumed colonial' % (cbr[y], cdr[y]) if y < 1800 else ('Haines/anchored' if y < 1960 else 'World Bank')} [{vit_q}]")
        return q, note

    path = os.path.join(BASE, "inputs.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["year", "immigration_M", "birth_rate", "death_rate",
                    "tfr", "cbr_per_1000", "cdr_per_1000",
                    "quality", "quality_note"])
        for y in YEARS:
            beta = 2.0 * cbr[y] / 1000.0
            q, qnote = quality(y)
            w.writerow([y, f"{imm[y]:.6f}", f"{beta:.6f}", f"{mu[y]:.6f}",
                        f"{wb_tfr.get(y, '')}", f"{cbr[y]:.2f}",
                        f"{cdr[y]:.2f}", q, qnote])
    print("wrote", path, f"({len(YEARS)} rows, eps={eps})")
    return imm, unauth


def build_validation():
    pop = _load_wb("SP.POP.TOTL")
    mig = _load_wb("SM.POP.TOTL")
    path = os.path.join(BASE, "validation.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["year", "total_pop_M", "migrant_stock_M", "fb_census_M",
                    "source"])
        for y in range(1980, 2025):
            w.writerow([
                y,
                f"{pop[y] / 1e6:.3f}" if y in pop else "",
                f"{mig[y] / 1e6:.3f}" if y in mig else "",
                f"{FB_CENSUS[y]:.2f}" if y in FB_CENSUS else "",
                "WB SP.POP.TOTL; WB SM.POP.TOTL (5-yr); Census decennial",
            ])
    print("wrote", path)


def build_validation_long():
    """Population anchors 1650-2024 for the full-span validation plot."""
    pop = _load_wb("SP.POP.TOTL")
    rows = []
    for y, v in sorted(COLONIAL_POP.items()):
        rows.append((y, v / 1e6, "Census HSUS Series Z 1-19 (colonial est.)"))
    for y, v in sorted(DECENNIAL_POP.items()):
        rows.append((y, v / 1e6, "U.S. Census decennial"))
    for y in sorted(pop):
        if y >= 1980:
            rows.append((y, pop[y] / 1e6, "World Bank SP.POP.TOTL"))
    path = os.path.join(BASE, "validation_long.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["year", "pop_M", "source"])
        for y, v, s in rows:
            w.writerow([y, f"{v:.3f}", s])
    print("wrote", path, f"({len(rows)} anchors)")


def build_initial_condition():
    # 1980 IC (millions), n-generic. P0 is data (1980 Census foreign-born).
    # The P1..Pn split is a documented parametric assumption: no dataset
    # measures pedigree depth, and CPS parental-nativity data only begin
    # in 1994. P1 ~= second generation (rough estimate); P2..P9 thin
    # (immigration was low 1930-1970, so few mid-depth lineages);
    # P10..P(n-1) = old stock; Pn is the residual so the total matches
    # the 1980 population. For n = 12 this reproduces the documented
    # (10, 20), (11, 30) old-stock split exactly.
    n = int(os.environ.get("GENEALOGY_MAX_DEPTH", "14"))
    total_1980 = 227.225  # World Bank SP.POP.TOTL 1980 (Census: 226.54M)
    rows = [
        (0, 14.08, "1980 Census foreign-born (data)"),
        (1, 20.0, "assumption: ~2nd generation (no 1980 measurement exists)"),
        (2, 4.0, "assumption: thin mid-depth (low immigration 1930-1970)"),
        (3, 4.0, "assumption: thin mid-depth"),
        (4, 4.0, "assumption: thin mid-depth"),
        (5, 5.0, "assumption: thin mid-depth"),
        (6, 6.0, "assumption: thin mid-depth"),
        (7, 7.0, "assumption: thin mid-depth"),
        (8, 8.0, "assumption: thin mid-depth"),
        (9, 9.0, "assumption: thin mid-depth"),
    ]
    if n == 12:
        rows += [(10, 20.0, "assumption: old stock"),
                 (11, 30.0, "assumption: old stock")]
    else:
        top = list(range(10, n))  # old-stock classes below the cap
        wsum = sum(range(1, len(top) + 1))
        for k, c in enumerate(top, start=1):
            rows.append((c, round(50.0 * k / wsum, 3),
                         "assumption: old stock"))
    assigned = sum(v for _, v, _ in rows)
    pn = round(total_1980 - assigned, 3)
    rows.append((n, pn, f"residual: 1980 total minus classes 0-{n - 1}"))
    assert abs(sum(v for _, v, _ in rows) - total_1980) < 1e-6
    path = os.path.join(BASE, "initial_condition_1980.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["class", "millions", "basis"])
        w.writerows(rows)
    print("wrote", path, f"(total {sum(v for _, v, _ in rows):.3f}M)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps", type=float, default=None,
                    help="emigration rate for the colonial residual "
                         "(default: read data/emigration_calibration.txt)")
    a = ap.parse_args()
    eps, eps_src = _resolve_eps(a.eps)
    print(f"colonial residual uses eps={eps} ({eps_src})")
    imm, unauth = build_inputs(eps)
    # Report the unauthorized addition for the record.
    tot_unauth = sum(unauth.values())
    print(f"unauthorized inflow 1980-2024 total: {tot_unauth:.2f}M "
          f"(gross, approx)")
    build_validation()
    build_validation_long()
    build_initial_condition()
    print("done ->", BASE)


if __name__ == "__main__":
    main()
