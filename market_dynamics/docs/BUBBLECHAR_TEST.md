# Bubble characteristics test — PRE-REGISTERED DESIGN (frozen 2026-09-28, before any predictive test)

## Objective
Test whether Greenwood–Shleifer–You (2019, JFE 133(2):291–308) bubble characteristics
predict crashes / future returns for **single S&P 500 stocks** in our panel (1999–2022),
as an additional return-forecasting mechanism. Same honest-negative discipline as rev-5.

## Source definitions (from GSY 2017 WP, verified against the paper text 2026-09-28)
- Episode: industry with >100% raw AND net-of-market return over trailing 2 years. 40 US episodes, 1926–2014.
- Crash: 40% drawdown within 2 years of identification.
- Acceleration: R_{t-24→t} − R_{t-24→t-12} (2-yr return minus first-year return, simple returns, natural units).
  "Measures the convexity of the price path ... how much of the price appreciation has occurred most recently."
- Change in volatility: 1-year change in the cross-sectional percentile rank of (daily-return) volatility.
  Crash episodes: Δ=+0.09 vs non-crash Δ=−0.03 (diff 0.11, p<5%).
- Age tilt: equal-weighted industry return minus age-weighted industry return
  (run-up concentrated in young firms → crash). Firm age = years since first CRSP/Compustat appearance.
- Issuance: % of firms with split-adjusted share-count increase ≥5% in past year.
- Their crash-prob curve: 50%→100% 2-yr net run-up ⇒ crash prob 20%→53%; 150% run-up ⇒ 80%.
- Return regression: R_{i,t→t+24} = a + b·Char_{i,t} + u; Δvol, age tilt, issuance, acceleration
  predict raw/excess returns; only Δvol and age tilt predict net-of-market returns (US).

## Adaptation to our panel (every deviation documented)
| # | GSY | This test | Reason |
|---|-----|-----------|--------|
| 1 | FF49 industries | Single S&P 500 member stocks | Our model universe |
| 2 | All CRSP 1926–2014 | panel_v2.csv (clean, post-unit-fix), 1999-12–2022-12 | Data availability |
| 3 | Daily-return volatility | Trailing-12m std of **monthly** log returns, cross-sectional percentile rank each month; Δ = rank_t − rank_{t-12} | Single data source; monthly is adequate for a 1-yr change |
| 4 | Age = first CRSP/Compustat appearance | Age = months since first appearance in panel_v2 (**left-censored at 1999-12** — documented proxy) | No CRSP/Compustat access |
| 5 | Age tilt = EW − age-weighted industry return | Stock's own cross-sectional **age percentile rank** (higher = older); secondary: young-firm dummy (age < 60m) | Within-industry young/old split impossible at single-stock level |
| 6 | Issuance from Compustat | **OMITTED** | No Compustat access. A share-count proxy (split-adjusted shares_use) was rejected: dividend drift contaminates the split adjustment and the panel has known share-unit history. Documented omission, not faked. |
| 7 | Turnover, CAPE, B/M, sales growth | **OMITTED** | Turnover: GSY null result. CAPE: market-level, constant across stocks at a formation date. B/M, sales: need Compustat. |
| 8 | 40 episodes | ~180 usable episodes (see power note) | Single stocks run up more often than industries |
| 9 | SUR clustering by calendar time | HC1 robust SEs; no calendar clustering | N too small for clustering |

**Central adaptation caveat:** a 40% drawdown is routine for a single stock but rare for an industry.
Base crash rates will be far higher than GSY's 53% — the comparison is directional, not level.

## Episode definition (locked)
- Universe: panel_v2 rows with has_price==1 and non-null adjclose.
- Trailing 24m raw simple return R24 = adjclose_t/adjclose_{t-24} − 1 > 1.00
- Trailing 24m net-of-market return (1+R24)/(1+Rmkt24) − 1 > 1.00, Rmkt from data/spx_monthly.csv log returns.
- Requires 25 consecutive valid months (t−24…t).
- **First-crossing**: first month a stock crosses both thresholds; 24-month embargo before the same stock may form another episode.
- Formation window: 2001-12 … 2020-12 (24m history back to 1999-12; forward window to 2022-12).
- Forward-data rule: ≥12 months of forward prices required, else episode dropped (document N dropped).
  Crash is evaluated over the available forward window capped at 24m. 24m forward *return*
  requires the full 24 months (smaller N, documented).
- **Attrition bias note (pre-registered):** delistings are likely informative (distress ⇒ crash
  understated; acquisitions ⇒ overstated). Net direction unknown — documented, not adjusted.

## Characteristics (locked formulas, measured with data through month t only)
1. **Acceleration** = R_{t-24→t} − R_{t-24→t-12} (simple returns). Expected sign on crash: **+**.
2. **Δvol_rank** = vol_pct_rank_t − vol_pct_rank_{t-12}, where vol = trailing-12m std of monthly
   log returns, ranked cross-sectionally each month (0–1). Expected sign on crash: **+**.
3. **age_rank** = cross-sectional percentile rank of firm age in months (higher = older).
   Expected sign on crash: **−** (younger ⇒ crash). Secondary: young = 1{age < 60m}.

## Crash & return definitions (locked)
- **Crash** = 1{min_{k=1..min(24,avail)} adjclose_{t+k}/adjclose_t − 1 ≤ −0.40}.
- **R24_fwd_raw** = adjclose_{t+24}/adjclose_t − 1 (requires 24 fwd months).
- **R24_fwd_net** = (1+R24_fwd_raw)/(1+Rmkt_fwd24) − 1.

## Tests (locked)
- **(a) Crash logit:** crash_i = Λ(α + β1·accel + β2·Δvol + β3·age_rank). Report coef, z, p,
  pseudo-R². Plus univariate crash/non-crash mean comparisons with t-tests (GSY Table 4 style).
- **(b) Return OLS:** R24_fwd_raw and R24_fwd_net on the three characteristics, HC1 SEs.
- **(c) Calibration:** overall crash rate vs GSY ~53% at 100% run-up; split episodes into
  100–150% vs >150% net run-up bins vs their 53%/80%.
- **Pass criteria (pre-registered):** (a) logit pseudo-R² > 0.05 AND ≥2 of 3 characteristics
  significant at 10% with expected signs; (b) ≥1 characteristic significant at 10% with
  expected sign for forward returns. Otherwise: HONEST NEGATIVE, no mechanism built.

## Power statement (written before seeing any predictive result)
~180 usable episodes (247 first-crossings; 67 lack ≥12m forward data). GSY had 40 US episodes
and still found significance, so power is comparable-to-better *if* the effect transfers —
but the transfer is the question: single-stock 40% drawdowns are common (high base rate),
which compresses signal. A null here is **weak evidence** (small N, large-cap-only universe,
1999–2022 window covering only dot-com/2008/2020, issuance omitted) — not a refutation of GSY.

## What was NOT examined before freezing this design
No crash/return relationship, no characteristic distribution by outcome, no regression
was run. Only episode counts (for the power statement) were computed.

---

# RESULTS (added 2026-09-28, after frozen design; no tuning after seeing outcomes)

## Sample realized
- 247 first-crossing episodes (183 distinct stocks); all 247 have ≥12m forward data
  (0 dropped at the 12m rule); 180 have full 24m forward returns.
- Formation years cluster: 2002 (46), 2010 (35), 2011 (26) — post-bust rebounds, not bubble peaks.
- **Key limitation:** the dot-com peak run-ups (1999–2000) are OUT OF SAMPLE — the panel
  starts 1999-12, so the earliest possible formation is 2001-12 (first realized: 2002-01).
  The canonical US bubble episode is missing.

## Crash calibration vs GSY
- Our crash rate: **0.186** (46/247) vs GSY ~0.53 at 100% run-up.
- By run-up size: 100–150% net → 30/161 = **0.186**; >150% net → 16/86 = **0.186**.
  No gradient with run-up size (GSY: 53% → 80%). The identical rates are a verified
  coincidence, not a bug (30/161=0.1863, 16/86=0.1860).
- Median forward drawdown: −56.8% (crash=1) vs −11.6% (crash=0) — the crash
  definition discriminates; the base rate is just low. Single-stock S&P 500 run-ups in
  1999–2022 rarely ended in 40% drawdowns within 2 years (many are post-bust rebounds:
  LMT/DLX in 2002, NVDA +444% forward in 2016, AMD +192% in 2019).

## (a) Crash logit (IRLS; HC not needed for logit z)
| char | coef | z | p | expected sign | verdict |
|------|------|---|---|---------------|---------|
| accel | −0.390 | −1.48 | 0.138 | + | WRONG SIGN |
| dvol | −0.048 | −0.08 | 0.937 | + | null |
| age_rank | +0.551 | +0.70 | 0.484 | − | WRONG SIGN |
- Pseudo-R² = **0.0134**. Significant-with-expected-sign: **0/3**.
- Univariate (GSY Table 4 style): accel crash-mean 0.880 vs non-crash 1.088 (t=−1.74,
  p=0.085, wrong sign); dvol −0.054 vs −0.041 (t=−0.35, p=0.73); age_rank 0.535 vs 0.512
  (t=+0.66, p=0.51). The only near-significant difference goes the WRONG way:
  non-crashing run-ups were MORE abrupt.

## (b) Forward-return OLS (HC1)
- 24m raw return (n=180): accel +0.088 (t=+1.35, p=0.18); dvol +0.048 (t=+0.40, p=0.69);
  age_rank −0.071 (t=−0.38, p=0.70). R²=0.012.
- 24m net-of-market (n=180): accel −0.007 (t=−0.12, p=0.90); dvol −0.016 (t=−0.10, p=0.92);
  age_rank +0.091 (t=+0.53, p=0.60). R²=0.001.

## Verdict vs pre-registered criteria
- Pass (a) required pseudo-R²>0.05 AND ≥2/3 significant with expected signs → **FAIL**
  (0.0134; 0/3).
- Pass (b) required ≥1 characteristic significant at 10% with expected sign → **FAIL**.
- **HONEST NEGATIVE: no mechanism built.** Nothing in GSY's price-based characteristic
  set (acceleration, Δvolatility, age/new-old tilt) predicts crashes or 2-yr returns for
  single S&P 500 stocks, 1999–2022.

## Interpretation (weak evidence, not refutation — per the pre-registered power statement)
1. Wrong universe for the question: GSY's industries diversify away idiosyncratic noise;
   single-stock run-ups are dominated by post-bust rebounds (2002, 2010–11), where high
   acceleration rationally reflects recovery, not bubble abruptness.
2. Missing the canonical episode: dot-com peak out of sample.
3. Issuance omitted (no Compustat) — one of GSY's strongest characteristics untested.
4. Large-cap-only, 1999–2022: three stress episodes, N=247 but only 46 crashes.

## Files
- `src/bubblechar.py` — implementation of the frozen design (IRLS logit + HC1 OLS in
  numpy/scipy; statsmodels unavailable in this environment).
- `data/bubblechar_episodes.csv` — 247 episodes with characteristics and outcomes.
- `data/bubblechar_results.json` — full statistics.
- Uncommitted, per task instructions.
