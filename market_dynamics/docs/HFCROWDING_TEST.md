# Hedge-Fund Crowding → Returns: 13F-based test

**Target paper:** Brown, Howard & Lundblad (2022), "Crowded Trades and Tail Risk",
*Review of Financial Studies* 35(7):3231–3271.

**Status: DESIGN FROZEN BEFORE ANY PREDICTIVE REGRESSION** (this section written
2026-09-28, before predictive tests ran; results appended after).

## 1. What BHL did and why we cannot replicate it exactly

BHL build **four crowdedness measures from a proprietary hedge-fund holdings
database** (~6,000 stocks, 2004–2016): holdings value, number of funds,
percent of shares outstanding, and HF holdings / average daily volume.
Their headline: high-minus-low crowdedness earns **+2.8%/yr value-weighted
(t=6.2)** and +9.1%/yr equal-weighted (t=8.0) — crowding is a **compensated
risk factor**, not a mispricing to fade — but crowdedness explains downside
tail risk in distress episodes.

We have no access to their database. The closest **public-data** analogue is
SEC EDGAR **13F-HR** filings (quarterly, 45-day filing lag, long US equity
positions only — shorts, non-US holdings, and sub-$100M managers are
invisible). This test is therefore a **13F-based approximation**, not a
replication.

## 2. Our construction (documented, assumed choices)

- **HF subset (assumed choice):** 18 hedge-fund managers' 13F-HR filings,
  CIKs verified via EDGAR 2026-09-28: Bridgewater (0001350694), Millennium
  (0001273087), AQR (0001167557), Elliott (0001048445), Lone Pine (0001061165),
  Tiger Global (0001167483), Coatue (0001135730), Renaissance (0001037389),
  Third Point (0001040273), Pershing Square (0001336528), Greenlight (0001079114),
  Soros (0001029160), Farallon (0000909661), Davidson Kempner (0000937617),
  Sculptor (0001054587), King Street (0001218199), Two Sigma (0001179392),
  SAC Capital (0001018103, covers the 1999–2013 era before becoming Point72).
  Dropped for insufficient coverage: D.E. Shaw, Viking, Baupost, Appaloosa,
  Moore, Angelo Gordon, Citadel (partial), Point72-era, Anchorage, Och-Ziff.
- **Parsing:** `infotable.xml` preferred (lxml, namespace-aware, SH shares
  only); best-effort CUSIP-anchored regex fallback for pre-XML paper-format
  filings; amendments (13F-HR/A) excluded; per (manager, report-quarter) the
  latest filingDate is kept. Code: `src/hf13f_build.py`.
- **CUSIP→ticker:** OpenFIGI batch mapping (vintage 2026-09-28); kept only
  tickers that are point-in-time S&P 500 members in our panel.
- **Universe:** panel members at formation (`data/panel_v2.csv`, clean
  post-unit-fix, 1999-12–2022-12). Large-cap only — BHL's effect may live in
  smaller stocks; documented limitation.
- **Crowding measures** per (ticker, holdings-quarter q):
  1. `n_hf` — number of HF managers holding (≈ BHL "# funds")
  2. `hf_own` — HF shares / shares outstanding (≈ BHL "% shares outstanding";
     shares out ≈ mcap/adjclose at quarter-end)
  3. `hf_adv_ratio` — HF shares / trailing-3-month average daily volume
     (≈ BHL "HF holdings / ADV")
  4. `hf_hhi` — Herfindahl of HF holder share fractions (extra concentration
     measure, no BHL analogue)

## 3. Timing (45-day 13F lag, honest)

- Holdings as of quarter-end q (e.g. 2000-03-31).
- Filing deadline is +45 days → public by ~mid-May.
- **Formation:** close of month q+2 (May 31 for Q1 holdings).
- **Holding periods:** 1 month and 3 months from formation.
- Forward returns from panel month-end `adjclose` (months are last trading
  days, verified in the dispersion test).

## 4. Pre-registered tests

**Sign hypothesis (from BHL's published result):** high-minus-low crowdedness
earns **positive** average returns. We test **two-sided** (crowding could also
be zero or negative in our large-cap universe), but a "replication" requires
the positive sign.

- **T1 — Quintile spreads:** each formation date, sort stocks into quintiles
  by each measure; Q5−Q1 long-short forward return, value- and equal-weighted.
  Report mean (annualized), Newey-West t-stat (3 lags quarterly / 6 monthly).
- **T2 — Fama-MacBeth:** cross-sectional regression of forward returns on
  crowdedness rank scaled to [0,1]; report mean slope, t-stat. Robustness:
  with log-size control.
- **Sample:** formations from holdings quarters 1999Q4–2021Q2 (last formation
  Aug-2021, 3-mo hold ends Nov-2021, inside panel).
- **Split:** estimate 2000–2007 formations; **walk-forward 2008–2021**.
- **Pass criteria (all required):** (a) in-sample |NW t| > 2.0 with the
  positive sign; (b) walk-forward mean spread > 0 with positive sign;
  (c) walk-forward hit rate > 0.5.
- **Benchmarks:** zero (naive); BHL's published +2.8%/yr VW (t=6.2).

**Explicit non-goals:** no trading mechanism is built unless all pass criteria
are met; a failure is reported as an honest negative (non-replication, not
refutation — different DB, universe, and period).

## 5. Results

*(appended 2026-09-28 after predictive tests ran; design above frozen beforehand)*

**Verdict: HONEST NON-REPLICATION — all three pass criteria fail.**
No forecasting mechanism or paper portfolio is built.

### 5.1 Q5–Q1 crowdedness spreads (value-weighted, annualized)

| measure | hor | full ann. | NW t | IS (2000–07) ann. | IS t | WF (2008–21) ann. | WF hit |
|---|---|---|---|---|---|---|---|
| n_hf | 1m | −6.5% | −1.10 | −7.8% | −0.56 | −5.7% | 0.39 |
| n_hf | 3m | −3.6% | −0.86 | −8.6% | −0.91 | −0.2% | 0.56 |
| hf_own | 1m | +7.8% | +0.94 | +19.5% | +1.08 | −0.4% | 0.44 |
| hf_own | 3m | +8.2% | +1.70 | +16.9% | +1.80 | +2.1% | 0.44 |
| hf_adv_ratio | 1m | +14.2% | +1.48 | +30.2% | +1.48 | +3.5% | 0.47 |
| hf_adv_ratio | 3m | +4.4% | +0.89 | +16.9% | +1.69 | −4.0% | 0.42 |
| hf_hhi | 1m | +8.9% | +1.15 | +18.9% | +1.16 | +1.9% | 0.58 |
| hf_hhi | 3m | +4.7% | +1.10 | +11.7% | +1.39 | −0.2% | 0.33 |

61 quarterly formations (1999Q4–2021Q3 holdings → formations Feb-2000…Nov-2021);
16,415 formation-stock observations. Against BHL's published +2.8%/yr VW
(t=6.2): our closest full-sample analogue (hf_own 3m, +8.2%, t=+1.70) is
nowhere near significance, and its walk-forward leg is +2.1% with a 0.44 hit
rate — indistinguishable from noise.

### 5.2 Fama–MacBeth slopes (rank[0,1] crowdedness → forward return, annualized)

| measure | hor | slope | NW t | + log-size control slope | t |
|---|---|---|---|---|---|
| n_hf | 1m | −0.110 | −0.96 | −0.063 | −0.77 |
| hf_own | 1m | +0.147 | +1.59 | +0.044 | +1.49 |
| hf_adv_ratio | 1m | +0.130 | +1.43 | +0.060 | +1.36 |
| hf_hhi | 1m | +0.118 | +1.07 | +0.058 | +0.69 |
| n_hf | 3m | −0.035 | −0.46 | −0.009 | −0.17 |
| hf_own | 3m | +0.090 | +2.13 | +0.028 | +1.22 |
| hf_adv_ratio | 3m | +0.030 | +0.59 | −0.022 | −0.79 |
| hf_hhi | 3m | +0.045 | +0.60 | +0.015 | +0.31 |

The single marginal signal (hf_own 3m, t=+2.13) collapses to t=+1.22 once
log market-cap is controlled — the crowdedness proxy was loading on size.

### 5.3 Pass criteria audit

- (a) IS |t|>2.0 positive: **FAIL** (max IS t = +1.80, hf_own 3m spread).
- (b) WF mean spread > 0: **FAIL** (mixed signs, largest |WF| = −5.7% n_hf 1m).
- (c) WF hit rate > 0.5: **FAIL** (5 of 8 cells below 0.5; best 0.58).

### 5.4 Why this is a non-replication, not a refutation

1. BHL's four measures come from a **proprietary hedge-fund holdings database**;
   ours come from public 13F-HR long positions of an **assumed** 18-manager
   subset — a different, much narrower lens (no shorts, no non-US, 45-day lag).
2. CUSIP→ticker mapping covered only 25.6% of distinct CUSIPs (Yahoo search;
   OpenFIGI unreachable), so n_hf and the ownership ratios are understated,
   especially for early-period delisted names.
3. Different universe (point-in-time S&P 500 vs BHL's broader cross-section)
   and different crowding construction throughout.

The sign pattern is weakly BHL-consistent for hf_own / hf_adv_ratio / hf_hhi
(positive but insignificant, not walk-forward robust) and outright negative
for n_hf. Nothing here clears the bar for a forecasting mechanism.

**Deviation note:** holdings window ran 1999Q4–2021Q3 (pre-reg said –2021Q2);
the extra 2021Q3 quarter adds one formation (2021-11-30) inside the
walk-forward window. Pre-1999Q4 filings (1999Q1–Q3, real but outside the panel
start) were dropped.

## References

- Brown, Howard & Lundblad (2022), *Rev. Financ. Stud.* 35:3231–3271 —
  crowdedness measures, +2.8%/yr VW premium, compensated-risk interpretation.
- SEC EDGAR 13F-HR filings (point of use: §2 holdings build; accessed
  2026-09-28 via sec.gov EDGAR full-text/archive endpoints).
- Yahoo Finance search API (point of use: §2 CUSIP→ticker map; OpenFIGI
  unreachable from this network during the build).
- Project panel `data/panel_v2.csv` (point of use: §2–3 point-in-time S&P 500
  membership, prices, shares outstanding, volume).
- Results: `data/hf_crowding_quarterly.csv`,
  `data/hf_crowding_test_results.json`; code `src/hf13f_{build,map,agg,test}.py`.
