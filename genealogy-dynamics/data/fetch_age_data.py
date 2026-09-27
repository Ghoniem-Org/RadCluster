#!/usr/bin/env python3
"""Download age-structured demographic inputs into data/raw/ (idempotent).

Sources (all public, no credentials):
- UN World Population Prospects 2024, "Fertility by single age" bulk CSV:
  age-specific fertility rates (births per 1,000 women) by single year of
  age of mother, USA, annual 1950-2023 (estimates; VarID=2/Medium).
  For the United States the UN's input data are NCHS registered births
  classified by age of mother (see WPP2024 Data Sources / Methodology).
  -> raw/WPP2024_Fertility_by_Age1.csv.gz, then USA rows extracted to
     raw/wpp_usa_asfr.csv by extract_usa_asfr().
- Human Life-Table Database (Max Planck Institute for Demographic Research,
  https://www.lifetable.de), USA pooled file (years 1900-2023): period life
  tables with m(x) = age-specific death rate by single year of age and sex.
  Underlying sources per HLD reference list: NCHS National Vital Statistics
  Reports Vol. 57 No. 1 (decennial life tables), SSA Actuarial Study No. 116
  (Bell & Miller 1900-2100), Keyfitz & Flieger.
  -> raw/HLD_USA.zip -> raw/hld_usa/USA.csv
"""
import argparse
import csv
import gzip
import os
import urllib.request
import zipfile

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "raw")

UA = "genealogy-dynamics-prototype/1.0 (research; contact via repo)"

WPP_FERT_URL = ("https://population.un.org/wpp/assets/Excel%20Files/"
                "1_Indicator%20(Standard)/CSV_FILES/"
                "WPP2024_Fertility_by_Age1.csv.gz")
HLD_USA_URL = "https://www.lifetable.de/File/GetDocument/data%5CUSA%5CUSA.zip"

USA_LOCATION = "United States of America (and dependencies)"


def _get(url: str, path: str, refresh: bool = False) -> None:
    if os.path.exists(path) and os.path.getsize(path) > 0 and not refresh:
        print("cached:", path)
        return
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r, open(path, "wb") as f:
        f.write(r.read())
    print("fetched:", url, "->", path, f"({os.path.getsize(path)} bytes)")


def fetch_wpp_fertility(refresh: bool = False) -> str:
    path = os.path.join(RAW, "WPP2024_Fertility_by_Age1.csv.gz")
    _get(WPP_FERT_URL, path, refresh)
    return path


def fetch_hld_usa(refresh: bool = False) -> str:
    zpath = os.path.join(RAW, "HLD_USA.zip")
    _get(HLD_USA_URL, zpath, refresh)
    dest = os.path.join(RAW, "hld_usa")
    os.makedirs(dest, exist_ok=True)
    if refresh or not os.path.exists(os.path.join(dest, "USA.csv")):
        with zipfile.ZipFile(zpath) as z:
            z.extractall(dest)
        print("extracted ->", dest)
    else:
        print("cached:", os.path.join(dest, "USA.csv"))
    return os.path.join(dest, "USA.csv")


def extract_usa_asfr(refresh: bool = False) -> str:
    """USA, VarID=2 (Medium = estimates for 1950-2023), ages 10-54, annual.

    Writes raw/wpp_usa_asfr.csv: year, age, asfr_per_1000.
    """
    out = os.path.join(RAW, "wpp_usa_asfr.csv")
    if os.path.exists(out) and os.path.getsize(out) > 0 and not refresh:
        print("cached:", out)
        return out
    src = fetch_wpp_fertility(refresh)
    n_in = n_out = 0
    with gzip.open(src, "rt", encoding="utf-8-sig") as f, \
            open(out, "w", newline="") as g:
        r = csv.DictReader(f)
        w = csv.writer(g)
        w.writerow(["year", "age", "asfr_per_1000"])
        for row in r:
            n_in += 1
            if (row["Location"] == USA_LOCATION and row["VarID"] == "2"
                    and row["Variant"] == "Medium"):
                y = int(row["Time"])
                if 1950 <= y <= 2023:
                    w.writerow([y, int(row["AgeGrpStart"]),
                                f"{float(row['ASFR']):.3f}"])
                    n_out += 1
    print(f"extracted {n_out} USA ASFR rows -> {out} (scanned {n_in})")
    return out


def sanity_check():
    """Print implied TFR and mean age at childbearing for spot years."""
    path = os.path.join(RAW, "wpp_usa_asfr.csv")
    by_year = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            by_year.setdefault(int(row["year"]), []).append(
                (int(row["age"]), float(row["asfr_per_1000"])))
    for y in (1950, 1957, 1976, 2007, 2023):
        rows = sorted(by_year[y])
        tfr = sum(v for _, v in rows) / 1000.0
        tot = sum(v for _, v in rows)
        mac = sum((a + 0.5) * v for a, v in rows) / tot
        print(f"{y}: TFR={tfr:.2f}  MAC={mac:.1f} yr  "
              f"(ASFR 20-24={dict(rows).get(20, float('nan')):.1f})")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    fetch_wpp_fertility(a.refresh)
    fetch_hld_usa(a.refresh)
    extract_usa_asfr(a.refresh)
    sanity_check()
    print("done ->", RAW)


if __name__ == "__main__":
    main()
