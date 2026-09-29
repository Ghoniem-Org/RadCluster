# Systematic Buffett Screen — design, results, limitations

**Status:** paper portfolio backtest, 2009-07 → 2022-06 (156 months). Separate sleeve;
does not touch the regime model's results or verdicts.

## Academic grounding

Frazzini, Kabiller & Pedersen (2018), "Buffett's Alpha": Berkshire Hathaway's
returns decompose into exposures to **value (HML)**, **quality (QMJ:
profitability, growth, safety, payout)** and **betting-against-beta (BAB)**,
levered ~1.6× through Berkshire's insurance float. This build implements the
**unlevered stock-selection leg only**. The float leverage is noted, not applied.

## Design

**Universe:** true historical S&P 500 members each June (panel_v2.csv, clean
post-unit-fix), 2009–2021 rebalances.

**Scores** (cross-sectional z-scores, winsorized at ±3, point-in-time only):

| Component | Formula | z |
|---|---|---|
| Quality | ROE = NetIncomeLoss / StockholdersEquity (FY annual); ROA = NetIncomeLoss / Assets secondary | mean(z(ROE), z(ROA)) |
| Value | B/M = Equity / mcap(June); E/P = NetIncome / mcap(June) | mean(z(B/M), z(E/P)) |
| Low beta | trailing 252-trading-day market beta vs SPY (≥126 daily obs) | −z(beta) |

**Composite** = mean of the three component z-scores. **Portfolio:** top 30 by
composite, equal-weighted, rebalanced annually at end of June (Buffett holds
for years; monthly churn would be unfaithful). Negative book equity excluded.
Delisted holdings exit at the last available price with proceeds redistributed
(assumed).

**Timing (point-in-time):** June Y rebalance uses the latest fiscal year ending
≥120 days before June 30 Y (10-K filed by then in practice), end within 18
months — the Fama-French end-date convention. Note: SEC companyfacts
`filed` dates are unreliable for pre-2010 years (restated facts carry the
restating filing's date — verified on AAPL FY2007); for 2010+ the end-date
rule coincides with filed-date PIT.

## Data provenance

- Prices/returns/membership/market cap: **measured** (panel_v2.csv; monthly
  total returns from dividend-adjusted prices; true historical membership —
  no survivorship bias by construction).
- Fundamentals (NetIncomeLoss, StockholdersEquity, Assets): **measured** (SEC
  XBRL companyfacts, 10-K facts only, USD). 584/585 tickers fetched; 19
  values rejected by absolute-anchor bounds (e.g. >$1T net income). Tag
  fallbacks used for equity (incl. NCI) and revenues.
- Beta: **measured** (daily prices vs SPY).
- Coverage: XBRL mandate began 2009 — FY2006: 56 tickers, FY2007: 293,
  FY2008+: 424–569. June 2008 rebalance had only 3 eligible stocks after
  filters → **skipped**; backtest starts June 2009 (152 eligible).
  11 tickers (e.g. XOM, PNC, BSX) have 10-Q-only facts in companyfacts and
  are excluded — documented gap, not filled.
- FF factors (Mkt-Rf, SMB, HML, RF): **measured** (Ken French Data Library,
  202608 CRSP vintage). AQR QMJ/BAB downloads are bot-walled (verified
  2026-09-28) — not used; no proxy substituted.

## Results (2009-07 → 2022-06)

| | Buffett screen | SPY |
|---|---|---|
| CAGR | **14.8%** | 13.6% |
| Sharpe (excess, ann.) | **1.27** | 0.94 |
| Ann. vol | **11.1%** | 14.2% |
| Max drawdown | **−19.6%** | −20.0% |
| Growth of $1 | **$6.01** | $5.26 |
| Calendar-year hit rate | 7/14 | — |
| One-way annual turnover | 53% | — |

FF3 regression (excess returns, HC1 SE, R²=0.66, n=156):
alpha **+0.46%/mo (t=2.71)**, +5.5% annualized; Mkt-Rf **0.64** (t=15.3);
SMB **−0.31** (t=−4.5, large-cap); HML **+0.00** (t=0.06).

Calendar years (screen vs SPY): 2009 +16.3/+22.4, 2010 +15.0/+15.1,
2011 **+16.8/+1.9**, 2012 +8.5/+16.0, 2013 +29.3/+32.3, 2014 **+23.5/+13.5**,
2015 +3.2/+1.2, 2016 +14.5/+12.0, 2017 +21.2/+21.7, 2018 **+1.9/−4.6**,
2019 +26.4/+31.2, 2020 −0.6/+18.3, 2021 +31.0/+28.7, 2022 **−7.5/−20.0**
(2009 and 2022 partial).

## Reading

The screen **beats SPY on a risk-adjusted basis in this backtest** (Sharpe
1.27 vs 0.94, alpha t=2.71). The mechanism is the low-beta tilt: it wins
down years (2011, 2018, 2022) and lags raging bulls (2020: −0.6% vs +18.3%,
missing mega-cap tech — as Buffett himself did). Notably the HML loading is
zero: the quality component (high ROE) neutralizes the value tilt, so this
behaves as **quality + low-beta**, not deep value. Holdings pass the smell
test (MO, CL, PEP, BRK-B, WMT, MCD, KMB across years).

## Honest limitations

1. **Backtest, not pre-registered.** Design choices (top 30, equal weight,
   June rebalance, composite formula) are standard-academic but were mine,
   not frozen before seeing data — weaker evidence than the pre-registered
   negative tests. Researcher degrees of freedom apply.
2. **No float leverage.** Berkshire's ~1.6× leverage via insurance float is
   the other half of "Buffett's Alpha" and is not replicable here.
3. **No crisis deal-flow.** Buffett's preferred-stock-plus-warrants rescue
   deals are unavailable to a screen.
4. **No QMJ/BAB attribution** (AQR data inaccessible); factor story rests on
   FF3 + the 0.64 market beta.
5. **Coverage gaps:** pre-2009 XBRL thin (backtest starts 2009); 11 tickers
   excluded for 10-Q-only facts; no transaction costs modeled (53% annual
   one-way turnover would cost roughly 0.1–0.3%/yr at institutional rates —
   assumed, not modeled).
6. One 14-year window, dominated by the post-2009 bull market; the 2020
   underperformance shows the strategy's known failure mode.

## Files

- `src/fetch_fundamentals.py` — SEC XBRL fetch (10-K facts, absolute bounds)
- `src/build_fundamentals_panel.py` — annual panel builder
- `src/buffett_screen.py` — scoring + backtest + attribution
- `data/fundamentals/` — 584 per-ticker companyfacts extracts (measured)
- `data/fundamentals_annual.csv` — annual panel
- `data/buffett_eligibility.csv`, `data/buffett_holdings.json`,
  `data/buffett_portfolio_monthly.csv`, `data/buffett_results.json`
- `figs/buffett_screen_equity.png` — equity curve vs SPY
