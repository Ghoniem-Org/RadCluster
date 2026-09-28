# Short-Interest Index (SII) Return-Forecast Test — PRE-REGISTERED DESIGN

**Status: DESIGN FROZEN BEFORE ANY PREDICTIVE REGRESSION WAS RUN (2026-09-28).**
Results section to be appended after.

## 1. Hypothesis (Rapach, Ringgenberg & Zhou 2016, JFE 121:46–65)

Aggregate short-interest index (SII) — log of the equal-weighted mean firm-level
short-interest ratio, detrended — **negatively** predicts aggregate stock market
returns at horizons h = 1, 3, 6, 12 months. Paper: annual R² ≈ 13% in-sample /
≈ 11% out-of-sample; beats 14 Goyal–Welch predictors; mechanism is the cash-flow
channel (informed shorts). Full note: `research_notes/.../notes/source-03-rpz-short-interest.md`.

## 2. Index construction (pre-registered)

- **Source:** the authors' published replication data (`Returns_short_interest_data.xlsx`,
  "Short interest" sheet: EWSI = equal-weighted mean short-interest ratio, monthly,
  1973:01–2014:12), downloaded 2026-09-28 from the data page's Drive link
  (matthewringgenberg.com/data). This is the PUBLISHED index input, not our own build.
  The separate "SII update to 2021:12" link could not be retrieved (fetch failure);
  the test window is therefore limited by the published series end (2014:12).
- **SII_t = log(EWSI_t) − trend_t**, where trend_t is a **quadratic time trend fit by
  OLS on data strictly through month t** (recursive/expanding from 1973:01 — no
  lookahead anywhere in the index). Linear-trend version as pre-registered robustness.
- **Deviation from paper (documented):** the paper's in-sample tests use full-sample
  detrending; we use recursive detrending throughout (stricter, real-time).
- **Timing alignment:** SII observed at end of month t predicts log excess returns
  over months t+1 … t+h (standard predictive-regression timing; the published monthly
  series is already aligned this way).

## 3. Data

- **Dependent variable:** S&P 500 log EXCESS return (paper's own "Returns" sheet,
  S&P 500 VW, minus Rfree from the "GW variables" sheet), 1973:01–2014:12.
- **Test window:** estimation 2000:01–2007:12 (n=96); walk-forward 2008:01–2014:12
  (first origin 2007:12; for h=12 the last origin is 2013:12). Chosen to mirror the
  rev-5/dispersion/comomentum test windows (pre-2008 estimation, post-2008 validation).

## 4. Estimators and inference (pre-registered)

- OLS predictive regression: r_{t+1:t+h} = α + β·SII_t + ε, h ∈ {1,3,6,12}.
- Overlapping observations → **Newey–West t-statistics with lag = h** (h=1: lag 1).
- **In-sample:** slope sign, NW t-stat, R² (2000–2007).
- **Walk-forward:** at each origin τ (expanding from 2000:01–2007:12), re-estimate
  (α,β) on data through τ, forecast r_{τ+1:τ+h}; benchmark = expanding historical
  mean excess return through τ. Metrics: OOS R² = 1 − MSE(model)/MSE(mean)
  (Campbell–Thompson), OOS forecast–realized correlation, directional hit rate.

## 5. Pass / fail criteria (pre-registered)

**PASS** (all required):
1. In-sample (2000–2007): β < 0 with NW |t| > 2.0 at h = 12 (the paper's headline horizon).
2. Sign consistency: β < 0 at h = 1, 3, 6, with at least one additional |t| > 1.65.
3. Walk-forward: OOS R² > 0 at h = 12 (beats the historical mean out-of-sample).
4. No sign flip: walk-forward realized slope (regression of realized on predicted) > 0.

**Otherwise: HONEST NEGATIVE** — no mechanism built, no overlay, deployment verdict unchanged.

## 6. Known limitations (pre-registered)

- Published series ends 2014:12 → walk-forward is 2008–2014 (7 years), not 2008–2021.
- SII is built on the full CRSP universe (all NYSE/AMEX/NASDAQ), not S&P 500 members;
  this tests the published aggregate index, which is the point.
- Pre-2003 short interest is exchange-reported mid-month; monthly alignment is the
  paper's own convention.

## 7. Results (2026-09-28)

### 7.1 In-sample, 2000:01–2007:12 (real-time recursive-quadratic SII)

| h | β (SII) | NW t (lag=h) | R² | n |
|---|---------|--------------|----|---|
| 1 | +0.036 | +0.97 | 0.011 | 95 |
| 3 | +0.071 | +0.98 | 0.015 | 93 |
| 6 | +0.143 | +0.82 | 0.028 | 90 |
| 12 | +0.424 | +1.31 | 0.094 | 84 |

**The sign is POSITIVE at every horizon — the opposite of the paper's negative
predictor.** Criterion 1 FAILS (β>0, |t|=1.31<2.0 at h=12). Criterion 2 FAILS
(all β>0).

### 7.2 Walk-forward, origins 2007:12 → 2014:12−h (expanding from 2000–2007)

| h | n | OOS R² | corr(fc, realized) | hit rate | realized-on-predicted slope |
|---|---|--------|--------------------|----------|-----------------------------|
| 1 | 84 | +0.014 | +0.069 | 0.607 | +0.42 |
| 3 | 82 | +0.088 | +0.197 | 0.695 | +0.90 |
| 6 | 79 | +0.175 | +0.379 | 0.759 | +1.36 |
| 12 | 73 | +0.165 | +0.071 | 0.685 | +0.15 |

OOS R² > 0 at every horizon (criterion 3 passes arithmetically) — but with the
**positive** sign, i.e. the index predicts returns in the *opposite direction*
to the published theory. Criterion 4 passes (+0.15): forecasts are directionally
consistent in-sample → out-of-sample.

### 7.3 Diagnostics

- **Paper's own recipe** (full-sample quadratic detrend, their return series):
  1973–2014, h=12: β=−0.25, NW t=−1.87 (negative — directionally replicates the
  paper, marginally below |t|=2 with this exact specification).
  1990–2014: β=−0.24, t=−1.58.
- **Linear-detrended SII, 2000–2007, h=12:** β=+0.51, t=+3.01 — the positive sign
  in the modern window is significant under linear detrending, insignificant
  (t=+1.31) under quadratic. The *reversal* is not a detrending artifact; its
  *significance* is detrending-sensitive — a fragility flag against overclaiming it.

### 7.4 Verdict: HONEST NEGATIVE on the published hypothesis

Pre-registered criteria 1 and 2 fail: in the 2000–2014 window the published
**negative** short-interest effect does not hold — it reverses sign. No mechanism
built, no overlay, deployment verdict unchanged.

**Interpretation (post-hoc, not pre-registered, not traded on):** the paper's
effect is real directionally in the long 1973–2014 sample but unstable across
subsamples — consistent with a structural break around 2000 (secular growth of
hedge-fund shorting as hedging/market-making rather than informed pessimism;
persistently elevated short interest through the 2009–2014 "wall of worry" rally).
The positive-sign predictability in the modern window (OOS R² up to +0.17) is
intriguing but was not pre-registered and is not robust to detrending choice;
it is reported as a diagnostic, not a signal.

### 7.5 Limitations realized

- Published series ends 2014:12 → walk-forward covers 2008–2014 only (7 years);
  the authors' SII update to 2021:12 could not be retrieved (drive fetch failure).
  No FINRA-based splice was attempted (splice risk onto a CRSP-universe index).
- One bug found and fixed mid-build: walk-forward origins initially ran past
  2014:12−h, giving empty (zero) realized returns at the tail; fixed before
  final results.

## Files

- `data/sii_monthly.csv` — real-time SII (recursive quadratic + linear), 1978–2014
- `data/shortinterest_test_results.json` — in-sample + walk-forward statistics
- `/tmp/sii_test.py`, `/tmp/rrz_sii/` — reproducible script + authors' replication
  package (EWSI 1973–2014, paper's returns, GW predictors)
