# Comomentum → future momentum-lane returns: pre-registered test

**Date (design):** 2026-09-28 · **Code:** `src/comomentum.py` (to be built AFTER this design is frozen) · **Data:** `data/comomentum_monthly.csv`, `data/comomentum_test_results.json` (to be produced)
**Universe:** true historical S&P 500 membership, `data/panel_v2.csv` (clean post-unit-fix panel).
**Reference:** Lou & Polk (2022), "Comomentum: Inferring Arbitrage Activity from Return Correlations," *Review of Financial Studies* 35(7):3272–3302. Full literature note: `~/workspace/research_notes/crowding-dispersion-return-predictabilit-20260928-1652/notes/source-02-lou-polk-comomentum.md`.

## Question

The rev-5 concentration test found no reliable concentration→return link; the dispersion test found no dispersion→return link. Comomentum — abnormal return correlation among momentum cohorts — is the one crowding→return candidate the literature says survives (replicated, not overturned). Does it work in our panel, mapped onto our architecture?

## Measure definition (frozen before any predictive regression is run)

For each month-end *t*:

1. **Cohorts = the model's own momentum lanes.** Winners = lane `mom=2` members at *t* (`mom ≥ 0.30` trailing 12m return); losers = lane `mom=0` (`mom < 0`). `mom` is the panel's measured 12–1 month return on adjclose (`src/build_panel.py:100`). Rationale: the model IS the lane system; the target below is the lane's own future return, so the cohort definition is architecturally faithful.
2. **Residuals:** for each cohort stock, daily returns over the trailing 52 weeks (252 trading days ending at *t*; require ≥200 non-missing daily obs, else excluded — exclusion rate reported). Regress daily excess return on the Fama–French 3 factors (daily, Ken French data library — URL and download vintage documented in §Data provenance). Residuals *e_it*.
3. **Comomentum:** `CoMOM_t = mean pairwise correlation of e_it` across all stock pairs in the cohort. Primary series = average of winner and loser cohort values; the two components reported separately (pre-registered).
4. **Coverage floor:** months where either cohort has <15 qualifying stocks are skipped and listed (expected in deep-crash months when the winner lane empties).

## Pre-registered deviations from Lou & Polk (with reasons — fixed BEFORE running)

| # | Lou & Polk | This test | Reason |
|---|---|---|---|
| 1 | Winner/loser = top/bottom return **deciles** | Model's **absolute momentum lanes** (mom=2 / mom=0) | Architectural mapping: cohorts are the model's own measured states; target is the lane's future return |
| 2 | Formation = past **12–2** month returns | Panel `mom` = **12–1** | Uses the model's measured momentum state as-is; 1-month difference immaterial for a 52-week correlation window |
| 3 | NYSE/AMEX/NASDAQ universe | **S&P 500 true historical membership** | Our panel's universe; large-cap only (stated limitation) |
| 4 | CRSP delisting returns | Stocks missing ≥53 trading days in the holding window are **excluded**; exclusion rate reported | Yahoo daily has no delisting returns (stated survivorship caveat) |

**Robustness (pre-registered, not tuned after):** cross-sectional **tercile** cohorts (top/bottom third of `mom` at *t*) as an alternative cohort definition, reported alongside — not selected ex post.

## Targets (frozen)

- **Primary:** 12-month-ahead equal-weighted buy-and-hold return of the winner cohort formed at *t*, from daily data (t+1 day … t+12 months).
- **Secondary:** winner-minus-loser (WML) spread, same construction.
- Holding window uses the same ≥200-trading-day presence rule; exclusions reported.

## Sample split (frozen)

- **Series start:** first month-end with 252 trailing trading days of daily data. Daily files start 2000-01-03 → first *t* ≈ **2001-01-31** (verified in build; exact date recorded).
- **In-sample estimation:** 2001-01 – 2007-12.
- **Walk-forward validation:** 2008-01 – 2021-12 (expanding estimation window, zero lookahead; last formation month with a full 12-month target in the daily data).

## Estimation (frozen)

- `y_t = a + b · CoMOM_t + e_t`, OLS. Overlapping 12-month targets → **Newey-West t-stats, 11 lags**. Report *b*, NW *t*, R². Run for primary and secondary targets; report winner/loser components of CoMOM separately.
- **Walk-forward:** expanding window from 2008-01; forecast *y*; report corr(forecast, realized), OOS R² vs expanding historical mean, and sign hit rate (sign of forecast deviation from expanding mean vs sign of realized deviation).
- **Horse race:** `y_t = a + b1·CoMOM_t + b2·TopShare_t + b3·HHI_t + e_t` where TopShare/HHI are the rev-5 concentration measures at *t* (`src/concentration_test.py`). Question: does comomentum add anything beyond concentration?

## Validation criteria (ALL must hold → build the mechanism; ANY failure → honest negative)

1. **In-sample sign + significance:** *b* < 0 with Newey-West |t| > 2.0 (Lou & Polk sign: high comomentum → lower subsequent momentum returns).
2. **Component consistency:** winner and loser CoMOM components both show *b* < 0 (guards against a one-cohort fluke).
3. **Walk-forward:** OOS R² > 0 **and** corr(forecast, realized) > 0.15.
4. **Horse race:** CoMOM retains negative sign with |t| > 1.5 alongside the concentration controls.

## If validated → build

A comomentum-modulated expected-return overlay on the momentum lane: equations in REPORT.md, a comomentum time-series figure with high/low regimes marked, overlay added to the paper-trading protocol from 2008, REPORT.md + LaTeX report/slides updated, PDFs recompiled, Drive updated in place, commit on Stock-Market-Dynamics, never push.

## If failed → honest negative

No mechanism, no overlay, no paper track. This document gets a Results section with the statistics. A rejected mechanism is a result.

## Results (2026-09-28 — design above was frozen before these were computed)

### Series summary

`data/comomentum_monthly.csv`: **234 months, 2001-01-31 – 2022-12-30** (first month-end with a full 252-trading-day window; daily files start 2000-01-03). Mean 0.036, std 0.023, range 0.002–0.140. Winner cohort 0.039 / loser cohort 0.032 on average; ~83 winners / ~99 losers per month qualify (≥200 daily obs). Magnitudes sit in the Lou & Polk range — the measurement is sane; it just doesn't predict.

**43 months skipped** (<15 qualifying stocks in a cohort): 1999-12–2000-12 (no 252-day window yet); 2002-08–2004-03 (bear-market winner-lane depletion); **2008-09–2009-08** (crash: winner lane emptied); scattered months (2010-02, 2010-04, 2016-02/03, 2019-01, 2020-04). Stated limitation: the lane-cohort definition cannot measure comomentum exactly in the high-stress episodes where a crowding signal would matter most. Comomentum is essentially orthogonal to the rev-5 concentration measures (corr with top-lane share 0.04, with HHI −0.06) — a different, but non-predictive, quantity.

### Pre-registered tests

| Test | Result | Criterion | Verdict |
|---|---|---|---|
| In-sample win12 (2001–2007, n=72) | b=−0.46, **NW-t=−0.27**, R²=0.0006 | b<0, \|t\|>2.0 | **FAIL** |
| In-sample wml12 | b=−1.62, **NW-t=−1.22**, R²=0.023 | b<0, \|t\|>2.0 | **FAIL** |
| Component consistency (win12) | winner comp. b=**+0.63** (t=0.50, wrong sign); loser b=−4.66 (t=−1.05) | both b<0 | **FAIL** |
| Tercile robustness (win12) | b=−2.79, NW-t=−1.35 | — (descriptive) | no signal |
| Walk-forward (152 forecasts) | corr=**−0.095**, **OOS R²=−0.082**, hit=0.487 | OOS R²>0 and corr>0.15 | **FAIL** |
| Horse race, full sample (n=224) | comom b=−2.64, t=−2.34; C t=0.98; hhi t=1.06 | comom keeps sign, \|t\|>1.5 | passes arithmetically, **but on the full sample, not the pre-registered window** |

### Verdict: HONEST NEGATIVE

Criteria 1, 2, and 3 all fail: no in-sample significance (t=−0.27), the winner component has the wrong sign, and walk-forward is actively negative (correlation −0.095, OOS R² −0.08, hit rate below a coin flip). The full-sample horse-race t=−2.34 is reported for completeness but does not meet the pre-registered decision rule — and a predictor that cannot forecast walk-forward is not a mechanism. **No comomentum overlay is built; no paper track is started; the deployment verdict is unchanged.**

### Interpretation

The series itself behaves plausibly — comomentum spikes in 2008 (max 0.14) and 2020 (era mean 0.063 vs 0.021 in 2004–2007): crowded trading does rise in stress, exactly as Lou & Polk's mechanism says. But in this S&P 500-only, 2001–2022 sample, elevated comomentum does not reliably precede 12-month momentum-lane reversals. Differences from Lou & Polk that may matter: (a) large-cap-only universe vs all US stocks; (b) 2001–2007 in-sample (~72 overlapping obs) vs 1965–2015; (c) absolute-lane cohorts vs return deciles; (d) the 2008–09 crash months — arguably the most informative episode — are unmeasurable under the lane definition because the winner lane empties. This is a non-replication in our universe/sample, not a refutation of Lou & Polk.

## Data provenance (actual)

| Input | Label |
|---|---|
| Daily prices (Yahoo, `data/raw/*.csv`, 661 tickers) | measured |
| Panel membership + `mom` (`data/panel_v2.csv`) | measured |
| FF3 daily factors, Ken French data library, downloaded 2026-09-28 17:01 UTC from `https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_daily_CSV.zip` (file vintage through 2026-08-31) | measured (third-party) |
| 252-day window, ≥200-day rule, 15-stock floor, NW(11) | assumed (pre-registered) |
| Lane cohorts; 12–1 mom; early-coverage thinning (2001: ~57% of members have daily files) | assumed/limitation (pre-registered) |
| Test statistics, walk-forward metrics | mechanical |

## Known limitations (stated upfront)

- Early-sample daily coverage is thin (2001-01: 494 members, only ~280 with raw daily files) — early CoMOM is computed on ~57% of members; coverage series reported.
- S&P 500 large-cap universe only; Lou & Polk's result spans all US stocks.
- In-sample window (2001–2007, ~84 overlapping monthly obs) is short vs Lou & Polk's 1965–2015.
- No delisting returns — mild survivorship tilt in buy-and-hold targets, exclusion rate reported.
- Overlapping 12-month targets → effective sample far smaller than nominal *n*; Newey-West mitigates but does not cure.
