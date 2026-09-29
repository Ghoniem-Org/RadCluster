# Dispersion → next-month market return: empirical test

**Date:** 2026-09-28 · **Code:** `src/dispersion_test.py` · **Data:** `data/dispersion_monthly.csv`, `data/dispersion_test_results.json`
**Universe:** true historical S&P 500 membership, `data/panel_v2.csv` (clean post-unit-fix panel; max mcap $2.91T — sanity-checked, no absurd values).

**Question (from Nasr):** the rev-5 concentration test found no reliable link between *concentration* (top-bin share) and next-month returns. Does *cross-sectional dispersion* — the width of the return distribution — predict aggregate market returns? This is the Goyal & Santa-Clara (2003) hypothesis, tested on our own panel.

## Pre-registered design (and honest deviations)

- **Intended:** estimate 1996–2007, walk-forward validate 2008–2026.
- **Actual:** the panel runs 1999-12–2022-12 and SPX monthly through 2021-12, so: **in-sample 2000-03–2007-12 (n=94)**, **walk-forward 2008-01–2021-11 (167 forecasts)**. (2000-01/02 dropped: membership sparse before 2000-03.)
- Regress next-month S&P 500 simple return (from `data/spx_monthly.csv` log-returns) on each dispersion measure separately (univariate, as pre-registered).
- Walk-forward uses a strictly expanding estimation window — zero lookahead.
- Asymmetric version: logit for P(down month) on each measure, same split.

## Measures (month *t*, ~341 member stocks/month avg)

| Symbol | Definition |
|---|---|
| `D_ew` | Equal-weighted cross-sectional std of monthly stock returns (mean 0.080) |
| `D_vw` | Cap-weighted cross-sectional std, weights = mcap at *t*−1 renormalized over available (mean 0.063; avg mcap coverage 82%) |
| `V_gs` | Goyal–Santa-Clara average variance: per stock, Σ daily² returns within month *t* (Yahoo daily `adjclose`, ≥15 trading days required), averaged across stocks (mean 0.0126, avg 343 stocks) |

Correlations among measures: 0.60–0.78 (they track the same stress episodes).

## Results

### OLS: R(m,t+1) on dispersion(t) — in-sample 2000-03–2007-12

| Measure | slope | t (Newey-West, 3 lags) | R² |
|---|---|---|---|
| D_ew | −0.077 | −0.40 | 0.005 |
| D_vw | −0.301 | −2.43 | 0.041 |
| V_gs | −0.456 | −1.05 | 0.011 |

Nothing predicts with the G&S **positive** sign. D_vw is nominally significant but **negative** — and fragile: split the in-sample window and it vanishes (2000–2003: t=−1.32; 2004–2007: t=−0.69). With three measures × two specs, one |t|>2 of the wrong sign is what chance produces.

### Walk-forward 2008–2021 (expanding window, 167 forecasts)

| Measure | corr(forecast, realized) | hit rate | OOS R² vs hist. mean |
|---|---|---|---|
| D_ew | +0.03 | 0.64 | −0.004 |
| D_vw | +0.04 | 0.61 | −0.002 |
| V_gs | +0.09 | 0.64 | −0.055 |

Every OOS R² is **negative** — all three models lose to the historical mean. Hit rates sit *below* the naive always-long benchmark (up-month frequency = 0.66). There is no usable signal here.

### Logit: P(R(m,t+1) < 0) on dispersion(t)

| Measure | in-sample coef (LR p) | McFadden R² | OOS AUC | hit@0.5 (base down-rate 0.34) |
|---|---|---|---|---|
| D_ew | +5.24 (0.39) | 0.006 | 0.52 | 0.65 |
| D_vw | +15.65 (0.054) | 0.029 | 0.51 | 0.65 |
| V_gs | +22.18 (0.34) | 0.007 | 0.52 | 0.65 |

OOS AUC ≈ 0.51–0.52 is a coin flip; the 0.65 "hit rate" just reflects predicting "up" most of the time. Dispersion does not predict down months either.

## Relation to the Goyal & Santa-Clara debate

Goyal & Santa-Clara (2003, JFE) reported that *average stock variance* — essentially our `V_gs` — positively and significantly predicts aggregate market returns (1963–1999, all CRSP stocks, monthly). The subsequent literature substantially qualified this: Bali, Cakici, Yan & Zhang (2005) showed the result is driven by small/NASDAQ stocks and the late-1990s episode — value-weighting the average or excluding small caps kills it; Guo & Savickas and others found it weakens in extended samples; Welch & Goyal (2008) found it fails out-of-sample like most return predictors.

**Our test does not replicate G&S, and that is exactly what the critics would predict.** Our universe is S&P 500 members only — large caps, the segment where Bali et al. showed the effect is absent. Our `V_gs` has the *wrong* (negative) sign in-sample (t=−1.05) and zero walk-forward correlation (+0.09, OOS R² −0.06). The sample (2000–2021) also sits entirely after G&S's original window, consistent with the "did not survive its own out-of-sample" reading. This is a clean non-replication in the precise universe where the skeptical literature said to expect one.

## Verdict

**Honest negative.** Neither cross-sectional dispersion (equal- or cap-weighted) nor Goyal–Santa-Clara average variance reliably predicts next-month S&P 500 returns in our panel — not in-sample with the hypothesized sign, not walk-forward, not for down-month probability. Combined with the rev-5 concentration test, the panel supports neither a concentration→return nor a dispersion→return link. No mechanism built; nothing here changes the deployment verdict.
