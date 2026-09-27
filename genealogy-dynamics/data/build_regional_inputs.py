#!/usr/bin/env python3
"""Build regional (Census-region) inputs for the regional compartments model.

Two input sets:

1. Inter-regional migration rates (4x4 per-capita out-migration rates, 1/yr).
   Source: U.S. Census Bureau, American Community Survey 1-year
   State-to-State Migration Flows tables (2023, 2024), aggregated from the
   50 states + DC to the 4 Census regions (Northeast, Midwest, South, West).
   Benchmark years 2023/2024 are MEASURED; pre-2023 values are held at the
   2023 rates (ASSUMED); 2025 holds 2024 (INTERPOLATED/held).

2. Immigration regional allocation shares lambda_r(t).
   Source: DHS Office of Homeland Security Statistics, Yearbook of
   Immigration Statistics, Table 4 "Persons Obtaining Lawful Permanent
   Resident Status by State or Territory of Residence", FY 2014-2023,
   aggregated to Census regions. 2014-2023 MEASURED; pre-2014 held at 2014
   (ASSUMED); 2024-2025 held at 2023 (INTERPOLATED/held).

Caveats (documented, not hidden):
- ACS flows cover the population 1 year and over (newborns excluded) and
  come from survey estimates with margins of error (estimates used as-is).
- DHS LPR state-of-residence shares are used as the allocation pattern for
  ALL immigration I(t); the unauthorized-immigration component's regional
  geography is not identified in the DHS yearbook and is ASSUMED to follow
  the LPR pattern.
- DHS cells are rounded to the nearest 10 (disclosure protection).
- Puerto Rico / island areas / foreign country are excluded from the
  inter-regional system (the model covers the 50 states + DC).

Outputs:
- data/regional_migration_rates.csv : year, from_region, to_region,
  rate_per_yr, flow_M, from_pop1yr_M, quality
- data/regional_immigration_shares.csv : year, NE, MW, S, W, quality
- data/regional_inputs_provenance.md
"""

import os

import numpy as np
import openpyxl
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "raw", "regional")

# --- Census regions (state -> region) -----------------------------------------
REGION_STATES = {
    "NE": ["Connecticut", "Maine", "Massachusetts", "New Hampshire",
           "Rhode Island", "Vermont", "New Jersey", "New York", "Pennsylvania"],
    "MW": ["Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin",
           "Iowa", "Kansas", "Minnesota", "Missouri", "Nebraska",
           "North Dakota", "South Dakota"],
    "S": ["Delaware", "Florida", "Georgia", "Maryland", "North Carolina",
          "South Carolina", "Virginia", "District of Columbia",
          "West Virginia", "Alabama", "Kentucky", "Mississippi", "Tennessee",
          "Arkansas", "Louisiana", "Oklahoma", "Texas"],
    "W": ["Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico",
          "Utah", "Wyoming", "Alaska", "California", "Hawaii", "Oregon",
          "Washington"],
}
REGIONS = ["NE", "MW", "S", "W"]
STATE2REGION = {s.strip(): r for r, ss in REGION_STATES.items() for s in ss}
assert len(STATE2REGION) == 51  # 50 states + DC

REGION_FULL = {"NE": "Northeast", "MW": "Midwest", "S": "South", "W": "West"}


# --- 1. ACS state-to-state flows ----------------------------------------------
def parse_acs_table_wide(path):
    """Parse the 2023-style WIDE ACS table.

    Returns (flow, pop1yr): flow is a (51 x 51) DataFrame of movers
    (index = origin, columns = destination) in people/yr; pop1yr is a
    Series of 'Population 1 year and over' by state.
    """
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb["Table"]

    # origin states: row 7, estimate columns 10, 12, 14, ... (MOE at +1)
    hdr = next(ws.iter_rows(min_row=7, max_row=7, values_only=True))
    origins, ocols = [], []
    c = 10
    while c < len(hdr):
        name = hdr[c - 1]
        if name is None:
            break
        name = str(name).strip()
        if name in ("Total", "Puerto Rico", "U.S. Island Area3",
                    "U.S. Island Area", "Foreign Country"):
            break
        origins.append(name)
        ocols.append(c)
        c += 2
    assert len(origins) == 51, f"expected 51 origin cols, got {len(origins)}"

    # destination rows: column A; data rows 12..62 (Alabama..Wyoming)
    rows = list(ws.iter_rows(min_row=12, max_row=62, max_col=max(ocols) + 1,
                             values_only=True))
    dests = [str(r[0]).strip() for r in rows]
    assert len(dests) == 51 and set(dests) <= set(STATE2REGION), dests[:5]

    flow = pd.DataFrame(0.0, index=origins, columns=dests)  # origin x dest
    for r, dst in zip(rows, dests):
        for org, cc in zip(origins, ocols):
            if dst == org:
                continue  # diagonal is N/A ("different state" table)
            v = r[cc - 1]
            if v is None or (isinstance(v, str)
                             and v.strip() in ("N/A", "**", "***", "X")):
                continue
            flow.loc[org, dst] = float(v)
    pop1yr = pd.Series({dst: float(r[1]) for r, dst in zip(rows, dests)})
    return flow, pop1yr


def parse_acs_table_long(path):
    """Parse the 2024-style LONG ACS table (+ supplemental origin populations).

    Returns (flow, pop1yr): flow (origin x destination) movers in people/yr;
    pop1yr 'Population 1 year and over' by state (current-residence sheet).
    """
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)

    ws = wb["Table"]
    flow = pd.DataFrame(0.0, index=list(STATE2REGION), columns=list(STATE2REGION))
    for r in ws.iter_rows(min_row=9, values_only=True):
        dst, org = r[0], r[1]
        if dst is None:
            continue
        dst, org = str(dst).strip(), str(org).strip()
        if dst.startswith("Footnotes") or dst.startswith("Source"):
            break
        if dst not in STATE2REGION or org not in STATE2REGION or dst == org:
            continue
        v = r[2]
        if v is None or isinstance(v, str):
            continue  # 'X' (same state), 'N' (suppressed), etc.
        flow.loc[org, dst] = float(v)

    ws2 = wb["Supplemental - Current Res"]
    pop1yr = {}
    for r in ws2.iter_rows(min_row=9, values_only=True):
        st = r[0]
        if st is None:
            continue
        st = str(st).strip()
        if st.startswith("Footnotes") or st.startswith("Source"):
            break
        if st in STATE2REGION:
            pop1yr[st] = float(r[1])
    assert len(pop1yr) == 51, f"expected 51 pop rows, got {len(pop1yr)}"
    return flow, pd.Series(pop1yr)


def build_migration_rates():
    tables = {2023: ("State_to_State_Migration_Table_2023_T13.xlsx", parse_acs_table_wide),
              2024: ("State_to_State_Migration_Table_2024_T13.xlsx", parse_acs_table_long)}
    frames = []
    bench = {}  # year -> (flow44 in M, pop4 in M)
    for year, (fname, parser) in tables.items():
        flow, pop1yr = parser(os.path.join(RAW, fname))
        # aggregate to 4 regions (flow is origin x destination)
        fi = flow.copy()
        fi.index = fi.index.map(STATE2REGION)      # origin -> region
        fi.columns = fi.columns.map(STATE2REGION)  # destination -> region
        f44 = fi.groupby(level=0).sum().T.groupby(level=0).sum().T
        f44 = f44.reindex(index=REGIONS, columns=REGIONS).fillna(0.0)
        p4 = pop1yr.groupby(pop1yr.index.map(STATE2REGION)).sum().reindex(REGIONS)
        bench[year] = (f44, p4)
        for fr in REGIONS:
            for to in REGIONS:
                if fr == to:
                    continue
                fl = float(f44.loc[fr, to]) / 1e6      # M/yr, origin->dest
                frames.append(dict(year=year, from_region=REGION_FULL[fr],
                                   to_region=REGION_FULL[to],
                                   rate_per_yr=fl * 1e6 / float(p4.loc[fr]),
                                   flow_M=fl,
                                   from_pop1yr_M=float(p4.loc[fr]) / 1e6,
                                   quality="measured"))
    # benchmark-year matrices for the model
    matrices = {}
    for year in tables:
        f44, p4 = bench[year]
        M = np.zeros((4, 4))
        for i, fr in enumerate(REGIONS):
            for j, to in enumerate(REGIONS):
                if i == j:
                    continue
                M[i, j] = (float(f44.loc[fr, to]) / float(p4.loc[fr]))
        matrices[year] = M
    df = pd.DataFrame(frames)
    df.to_csv(os.path.join(BASE, "regional_migration_rates.csv"), index=False)
    return matrices, df


# --- 2. DHS LPR-by-state -> regional immigration shares ----------------------
def build_immigration_shares():
    url = "https://ohss.dhs.gov/topics/immigration/yearbook/2023/table4"
    tabs = pd.read_html(url)
    t = tabs[0]
    # first column header may vary; normalize
    t = t.rename(columns={t.columns[0]: "state"})
    t["state"] = t["state"].str.strip()
    t = t[t["state"].isin(STATE2REGION)].copy()
    assert len(t) == 51, f"expected 51 states, got {len(t)}"
    year_cols = [c for c in t.columns if str(c).isdigit()]
    t[year_cols] = t[year_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    t["region"] = t["state"].map(STATE2REGION)
    g = t.groupby("region")[year_cols].sum()
    shares = g.div(g.sum(axis=0), axis=1)  # region x year -> shares
    rows = []
    for y in year_cols:
        row = {"year": int(y), "quality": "measured"}
        for r in REGIONS:
            row[r] = float(shares.loc[r, y])
        rows.append(row)
    df = pd.DataFrame(rows).sort_values("year")
    df.to_csv(os.path.join(BASE, "regional_immigration_shares.csv"), index=False)
    return df


def main():
    print("parsing ACS state-to-state tables ...")
    matrices, mdf = build_migration_rates()
    print("fetching DHS LPR-by-state table ...")
    sdf = build_immigration_shares()

    # --- summary diagnostics -------------------------------------------------
    print("\n== inter-regional migration rates (per-capita, 1/yr) ==")
    for year in sorted(matrices):
        M = matrices[year]
        print(f"{year}: total inter-regional flow = "
              f"{(M * 1).sum():.6f}/yr (rate-weighted)")
    # print human-readable matrix for 2023 (rates per 1000)
    M = matrices[2023]
    print("\n2023 per-capita out-migration rates per 1000 (from -> to):")
    print(pd.DataFrame(M * 1000, index=REGIONS, columns=REGIONS)
          .round(2).to_string())

    print("\n== immigration regional shares (DHS LPR) ==")
    print(sdf.round(4).to_string(index=False))

    # --- provenance -----------------------------------------------------------
    with open(os.path.join(BASE, "regional_inputs_provenance.md"), "w") as f:
        f.write(PROVENANCE)
    print("\nwrote data/regional_migration_rates.csv, "
          "data/regional_immigration_shares.csv, "
          "data/regional_inputs_provenance.md")


PROVENANCE = """# Regional inputs provenance

## Inter-regional migration rates (`regional_migration_rates.csv`)

- **Source:** U.S. Census Bureau, American Community Survey 1-year
  State-to-State Migration Flows tables:
  - 2023: `State_to_State_Migration_Table_2023_T13.xlsx`
    (https://www2.census.gov/programs-surveys/demo/tables/geographic-mobility/2023/state-to-state-migration/State_to_State_Migration_Table_2023_T13.xlsx)
  - 2024: `State_to_State_Migration_Table_2024_T13.xlsx`
    (https://www2.census.gov/programs-surveys/demo/tables/geographic-mobility/2024/state-to-state-migration/State_to_State_Migration_Table_2024_T13.xlsx)
- **Aggregation:** the 50-state + DC origin x destination mover matrix is
  summed to the 4 Census regions (Northeast, Midwest, South, West) using
  the Census region definitions. Puerto Rico, island areas, and movers from
  a foreign country are excluded (the model covers the 50 states + DC).
- **Rates:** per-capita out-migration rate from region r to r' =
  flow(r->r') / (population 1 year and over in r), with the denominator
  taken from the same table's "Population 1 year and over" row
  (MEASURED, same source).
- **Quality labels:** 2023/2024 benchmark values = `measured`; values held
  at the 2023 benchmark before 2023 = `assumed`; 2025 held at 2024 =
  `interpolated` (held).
- **Caveats:** ACS flows cover the population 1 year and over only
  (newborns excluded); estimates are survey-based with published margins
  of error (point estimates used as-is).

## Immigration regional allocation shares (`regional_immigration_shares.csv`)

- **Source:** Department of Homeland Security, Office of Homeland Security
  Statistics, *Yearbook of Immigration Statistics*, Table 4 "Persons
  Obtaining Lawful Permanent Resident Status by State or Territory of
  Residence: Fiscal Years 2014 to 2023"
  (https://ohss.dhs.gov/topics/immigration/yearbook/2023/table4).
  50 states + DC aggregated to Census regions; shares sum to 1.
- **Quality labels:** FY 2014-2023 = `measured`; pre-2014 held at 2014 =
  `assumed`; 2024-2025 held at 2023 = `interpolated` (held).
- **Caveats:** DHS cells are rounded to the nearest 10 (disclosure
  protection). The national immigration inflow I(t) includes
  unauthorized immigration (DHS estimates folded into `data/inputs.csv`);
  its regional geography is NOT separately identified in the yearbook, so
  the LPR pattern is ASSUMED to apply to all of I(t).
"""

if __name__ == "__main__":
    main()
