#!/usr/bin/env python3
"""Download public demographic data into data/raw/ (idempotent; use --refresh).

Sources (all public, no credentials):
- DHS Office of Homeland Security Statistics, Yearbook of Immigration
  Statistics Table 1 (persons obtaining lawful permanent resident status,
  fiscal years 1820-2024) — parsed from the published HTML table.
- World Bank World Development Indicators API (no key required):
  SP.DYN.TFRT.IN (total fertility rate), SP.DYN.CBRT.IN (crude birth rate),
  SP.DYN.CDRT.IN (crude death rate), SP.POP.TOTL (total population),
  SM.POP.TOTL (international migrant stock).
- DHS OHSS, Estimates of the Unauthorized Immigrant Population Residing in
  the United States (January 2005-January 2022 reports; PDFs). Headline
  Jan-1 totals are transcribed in build_inputs.py from these reports
  (Baker 2021 Table A2-1 for 2000-2018; Baker & Warren 2024 Table A2-1
  for the full updated 2000-2022 series). NOTE: the January 2018-January
  2022 report URL carries a double-encoded en dash (%25E2%2580%2593);
  the plain-encoded form 404s.
"""
import argparse
import csv
import json
import os
import re
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "raw")

DHS_TABLE1_URL = "https://ohss.dhs.gov/topics/immigration/yearbook/2024/table1"
WB_INDICATORS = [
    "SP.DYN.TFRT.IN",   # fertility rate, total (births per woman)
    "SP.DYN.CBRT.IN",   # birth rate, crude (per 1,000 people)
    "SP.DYN.CDRT.IN",   # death rate, crude (per 1,000 people)
    "SP.POP.TOTL",      # population, total
    "SM.POP.TOTL",      # international migrant stock, total
]
UA = "genealogy-dynamics-prototype/1.0 (research; contact via repo)"


def _get(url: str, path: str, refresh: bool = False) -> None:
    if os.path.exists(path) and os.path.getsize(path) > 0 and not refresh:
        print("cached:", path)
        return
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=90) as r, open(path, "wb") as f:
        f.write(r.read())
    print("fetched:", url, "->", path, f"({os.path.getsize(path)} bytes)")


def fetch_dhs_lpr(refresh: bool = False) -> None:
    """DHS OHSS Yearbook Table 1 -> raw/dhs_lpr_fy.csv (fiscal_year, lpr)."""
    html_path = os.path.join(RAW, "dhs_ohss_table1_2024.html")
    _get(DHS_TABLE1_URL, html_path, refresh)
    txt = open(html_path, encoding="utf-8", errors="replace").read()
    # Strip footnote anchors inside year cells (e.g. 1976 carries endnote
    # [1]: the federal fiscal year end changed, so FY1976 includes the
    # transition quarter; DHS reports 499,090 for it).
    txt = re.sub(r'<a href="#endnote[^>]*>.*?</a>', "", txt)
    rows = re.findall(
        r"<tr[^>]*>\s*<td[^>]*>\s*(\d{4})\s*</td>\s*<td[^>]*>\s*([\d,]+)\s*</td>",
        txt,
    )
    if len(rows) < 100:
        raise RuntimeError(f"DHS Table 1 parse failed ({len(rows)} rows)")
    out = os.path.join(RAW, "dhs_lpr_fy.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["fiscal_year", "lpr"])
        for y, n in sorted(rows, key=lambda r: int(r[0])):
            w.writerow([y, n.replace(",", "")])
    print("wrote", out, f"({len(rows)} rows)")


def fetch_worldbank(refresh: bool = False) -> None:
    for ind in WB_INDICATORS:
        url = (f"https://api.worldbank.org/v2/country/USA/indicator/{ind}"
               f"?format=json&date=1960:2024&per_page=120")
        _get(url, os.path.join(RAW, f"wb_{ind}.json"), refresh)


# --- DHS OHSS unauthorized-immigrant-population reports (PDFs) -----------
# Only the headline Jan-1 totals are used downstream (see build_inputs.py).
DHS_UNAUTH_REPORTS = [
    ("1990_to_2000",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-1990-january-2000.pdf"),
    ("2005",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2005.pdf"),
    ("2006",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2006.pdf"),
    ("2007",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2007.pdf"),
    ("2008",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2008.pdf"),
    ("2009",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2009.pdf"),
    ("2010",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2010.pdf"),
    ("2011",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2011.pdf"),
    ("2012",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2012.pdf"),
    ("2014",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2014.pdf"),
    ("2015_2018",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2015-january-2018.pdf"),
    # The January 2018-January 2022 report URL carries a double-encoded
    # en dash (%25E2%2580%2593); the plain-encoded form 404s.
    ("2018_2022",
     "https://ohss.dhs.gov/sites/default/files/2024-06/2024_0418_ohss_estimates-of-the-unauthorized-immigrant-population-residing-in-the-united-states-january-2018%25E2%2580%2593january-2022.pdf"),
]


def fetch_dhs_unauth(refresh: bool = False) -> None:
    dest = os.path.join(RAW, "dhs_unauth")
    os.makedirs(dest, exist_ok=True)
    for name, url in DHS_UNAUTH_REPORTS:
        path = os.path.join(dest, name + ".pdf")
        try:
            _get(url, path, refresh)
        except Exception as e:  # keep going; note the gap
            print(f"WARNING: could not download {name}: {e}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true",
                    help="re-download even if cached files exist")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    fetch_dhs_lpr(a.refresh)
    fetch_worldbank(a.refresh)
    fetch_dhs_unauth(a.refresh)
    print("done ->", RAW)


if __name__ == "__main__":
    main()
