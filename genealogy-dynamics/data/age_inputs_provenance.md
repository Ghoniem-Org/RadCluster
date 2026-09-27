# Age inputs: provenance (data/age_inputs.csv)

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

## 2026-09-26 recalibration (pre-1900 mortality level)

The original backcast set lam = CDR_target/17.2 to match the *stationary*-population CDR. Because the colonial population is growing (young pyramid), its endogenous stable CDR came out 2-4/1000 below target, inflating 2025 population to ~440M. Recalibrated: for each year 1650-1899, lam(t) solves stable_CDR(beta(t), lam*mu_1900_02) = CDR_target(t) (brentq). lam(1650)=1.917, lam(1899)=0.929. Level is now *calibrated* to Haines CDR via inputs.csv; age pattern remains *assumed* (1900-02 schedule).

## 2026-09-26 recalibration (pre-1900 mortality level)

The original backcast set lam = CDR_target/17.2 to match the *stationary*-population CDR. Because the colonial population is growing (young pyramid), its endogenous stable CDR came out 2-4/1000 below target, inflating 2025 population to ~440M. Recalibrated: for each year 1650-1899, lam(t) solves stable_CDR(beta(t), lam*mu_1900_02) = CDR_target(t) (brentq). lam(1650)=1.917, lam(1899)=0.929. Level is now *calibrated* to Haines CDR via inputs.csv; age pattern remains *assumed* (1900-02 schedule).
