# Regional inputs provenance

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
  of error (point estimates used as-is). Cells marked `N` in the workbook
  (e.g. Connecticut<->Alabama 2024) mean "the estimate or margin of error
  cannot be displayed because there were an insufficient number of sample
  cases" (workbook footnote, MEASURED documentation); these negligible
  flows are treated as zero (ASSUMED).

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
