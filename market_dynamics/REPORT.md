# Stock-Market Dynamics, Revision 4

**Four-axis regime model with capital weighting, true historical membership, and a walk-forward paper-trading deployment protocol**

Nasr Ghoniem — 2026-09-28

---

## 1. Objectives

Rev-4 rebuilds the market-dynamics model on four requirements:

1. **Four regime axes** — momentum × volatility × turnover/liquidity × size (3×2×2×2 = 24 bins).
2. **Capital weighting** — bin shares are capital shares, not headcounts.
3. **True historical S&P 500 membership** — point-in-time constituents, never backfilled.
4. **A deployment protocol** — zero-lookahead walk-forward forecast, monthly tilt portfolio with frozen weights, benchmarked against SPY and an equal-weight historical-member portfolio, with monthly equity and holdings snapshots saved as the paper track record.

A further requirement, added during the build: the paper record is the **prerequisite** for any real-capital deployment. This report states the verdict plainly.

---

## 2. Data

### 2.1 Historical S&P 500 membership (measured, with caveats)

Index membership comes from the public reconstruction **fja05680/sp500** ("S&P 500 Historical Components & Changes"), file *S&P 500 Historical Components & Changes (Updated).csv*: 2,720 snapshots/change dates from 1996-01-02 through 2026-08-18 [1]. This is a community reconstruction, **not** an official S&P database; change dates may lag true effective dates by days. All membership joins use the latest snapshot on or before each month-end (point-in-time).

### 2.2 Prices (measured)

Daily adjusted closes from **Yahoo Finance** [2] for every ticker ever appearing in the membership file (659 raw files after alias resolution). Corporate-action aliases were validated by hand: `BK→BNY`, `BLL→BALL`, `ABC→COR`, `ANTM→ELV`, `FB→META`, `DISCA/DISCK→WBD`, `ACE→CB`, `DLPH→APTV`, `TMK→GL`, `STI/BBT→TFC`, `UTX/RTN→RTX`, `PX→LIN`.

Snapshot price coverage (fraction of members with a valid price): 2000: 57%, 2008: 65%, 2015: 76%, 2021: 90%, 2022: 94%. Coverage rises as Yahoo history deepens; early years lose delisted names permanently.

### 2.3 Shares outstanding and splits (measured / estimated)

**SEC EDGAR companyconcept XBRL** [3] (`SharesOutstanding`, `WeightedAverageNumberOfSharesOutstanding`) provides point-in-time shares: 47,902 raw facts, 47,841 positive after filtering, with 920 apparent 1,000×/1,000,000× unit errors mechanically rescaled. SEC coverage begins mid-2009; **no facts were filed by 2008-12-31**. Pre-2009 shares use earliest-filing backfill and are **estimated with lookahead limitations**, flagged per-row in `panel_v2.csv`. **Yahoo Finance** corporate-action history [2] supplies 1,844 split events; all share counts are split-adjusted.

**Turnover** is dollar volume divided by split-adjusted market cap — a **measured liquidity proxy**, not true share turnover, unless valid point-in-time shares exist (they do for ~83% of priced rows; the rest use backfilled shares and are flagged estimated).

### 2.4 Panel

`data/panel_v2.csv`: 138,255 historical-member/month rows with momentum (12–1m), volatility (12m daily σ, annualized), split-adjusted shares, market cap, dollar-volume turnover, filing dates, and provenance flags. 83.2% of priced rows have cap and turnover; 80.7% are binnable on all four axes. Sanity: AAPL 2007-12 ≈ $176.0B (backfilled, estimated); 2022-12 ≈ $2.067T (filing-based, measured).

---

## 3. Model

### 3.1 State

The capital distribution over 24 bins:

$$\mathbf{c}(t) \in \Delta^{23}, \qquad c_b(t) = \frac{\sum_{i \in b} \mathrm{cap}_i(t)}{\sum_i \mathrm{cap}_i(t)}.$$

Bin index $b = (m, v, \tau, s)$:

| Axis | Levels | Edges | Label |
|---|---|---|---|
| Momentum $m$ | 3 | $[-1, 0, 0.30, 3]$ (12–1m return) | **assumed** (round absolute levels) |
| Volatility $v$ | 2 | $[0, 0.28, 5]$ (annualized) | **assumed** |
| Turnover $\tau$ | 2 | $[0, 0.008, 10]$ (dollar vol / mcap) | **assumed** |
| Size $s$ | 2 | $[0, \$10\mathrm{B}, \infty)$ | **assumed** |

The edges are round numbers chosen for interpretability, **not** calibrated quantiles; several $(v{=}1,\tau{=}0)$ cells are frequently empty. This is disclosed, not hidden.

### 3.2 Master equation

$$\frac{d\mathbf{c}}{dt} = \underbrace{(T_{\mathrm{eff}} - I)\mathbf{c}}_{\text{drift}} + \underbrace{A(M)\,\mathbf{c}}_{\text{market advection}} + \underbrace{H(\mathbf{c})\,\mathbf{c}}_{\text{herding}} + \underbrace{\mathbf{e}}_{\text{entry/exit}},$$

with discrete monthly stepping.

**Mechanism — drift (measured).** $T$ is the cap-weighted month-to-month bin transition matrix, estimated on 2000–2007 (pre-episode, expanding thereafter). Two regimes: $T_{\mathrm{calm}}$ for $|M|<0.05$, $T_{\mathrm{stress}}$ otherwise (measured).

**Mechanism — market advection (mechanical).** The cross-sectional window-delta $M(t) = \bar r(t) - \bar r(t{-}12)$ shifts every stock's 12-month momentum window by the same amount — a pp is a pp, so $\kappa = 1$ is **mechanical**, not fitted. Fractional remap (conservative, multi-bin) translates each momentum lane; reflecting boundaries at the extreme bins pile mass (the spike mechanism). $w_{\mathrm{eff}} = 0.15$ is **assumed** (fraction of cap participating per month).

**Mechanism — herding (calibrated).** Signed autocatalytic flux $h_{\mathrm{eff}} = h_1 (w - \ell)$ toward the winning extreme; calibration on 2000–2007 selected **$h_1 = 0$** — herding did not improve the pre-2008 objective, so it is **off** in rev-4. Rev-3's $h_1 = 1.0$ is not carried over as calibrated.

**Mechanism — sticky drift (calibrated).** $T_{\mathrm{eff}} = (1-\lambda)T + \lambda I$ with $\lambda = \lambda_1 (w+\ell)$; calibration selected **$\lambda_1 = 0$**. Also off.

**Mechanism — spike-persistence damper (assumed, episode-tuned in rev-3).** At a true peak the extreme bin's mass sits deep in the tail (Mar-2021: median top-bin momentum 103%, only 11% within 8pp of the 60% boundary — measured in rev-3), so a uniform donor-cell drains ~2× too fast. The damper impedes extreme-bin outflow on driver-sign reversal out of a spike (gate: share > 0.5): $\mathrm{damp} = \min((|M|/B)^p, 1)$, $p = 1.2$, $B$ = 6-month-decay memory of the $|M|$ that built the spike. All three damper parameters are **assumed** (carried from rev-3 episode tuning, unidentified in the pre-2008 grid).

### 3.3 Driver forecast $\widehat M_{t+1}$ (measured, zero lookahead)

From the momentum definition, $\mathrm{mom}_i(t{+}1) - \mathrm{mom}_i(t) \approx r_i(t) - r_i(t{-}12)$, so

$$\widehat M(t{+}1) = R_{\mathrm{SPY}}(t) - R_{\mathrm{SPY}}(t{-}12),$$

both known at month-end $t$. This **mechanical nowcast** (not an AR(1) extrapolation) achieves correlation 0.50 with realized $M$ and 70% directional accuracy over 2008–2022, vs. −0.02 / 41% for an expanding AR(1), which collapses to zero and is discarded. The raw form overstates $|M|$ by ~1.6× (measured slope 0.63 on the full sample, disclosed); it is used unshrunk to avoid lookahead.

---

## 4. Parameter table

| Parameter | Value | Label | Basis |
|---|---|---|---|
| Bin edges (4 axes) | see §3.1 | assumed | round interpretable levels |
| $\kappa$ (advection gain) | 1.0 | mechanical | a pp is a pp |
| $w_{\mathrm{eff}}$ | 0.15 | assumed | monthly participation fraction |
| $h_1$ (herding) | 0 | calibrated | 2000–2007 grid: no improvement |
| $\lambda_1$ (sticky) | 0 | calibrated | 2000–2007 grid: no improvement |
| Damper $p$ / $\tau_{\mathrm{build}}$ / gate | 1.2 / 6 mo / 0.5 | assumed | rev-3 episode tuning; unidentified pre-2008 |
| Tilt $\gamma$ | 1.0 | assumed | additive tilt strength |
| Cost | 5 bps one-way | assumed | institutional estimate |
| $T_{\mathrm{calm}}$, $T_{\mathrm{stress}}$ | 24×24 | measured | cap-weighted, 2000–2007 → expanding |

---

## 5. Conditional hindcasts (realized $M(t{+}1)$ forcing)

These are **conditional hindcasts**, not forecasts: the model is driven by the realized destination-month driver. They test the distribution mechanics, not the forecast.

### 5.1 2008–09

Full-distribution RMSE 0.1685. Predicted bottom-lane peak 0.88 (2008-12) vs. realized 1.00 (2009-03): amplitude close, **three months early**. The realized 1.00 reflects near-total bottom-lane saturation at the March-2009 trough (audited: no single-name outlier; the entire cross-section sat below 0% twelve-month momentum). Post-peak decay is captured to within ~0.05 share. The 2009 recovery is under-predicted — the uniform-within-bin advection approximation cannot resolve the sharp V-bounce, a known structural limitation.

### 5.2 2020–21

Full-distribution RMSE 0.1110. Predicted top-lane peak 0.95 (2021-04) vs. realized 0.92 (2021-03): amplitude within 0.03, one month late. The realized series is volatile around the peak (0.92 in Mar → 0.49 in Apr → 0.90 in May → 0.86–0.87 through Jul); the model produces a smoother hump (0.58 → 0.95 → 0.61 → 0.43 → 0.35) and does not resolve the sharp April dip. The damper (rev-3) aligns the broad decay envelope; without it the model drains ~2× too fast.

**Mechanism before results.** Both episodes work for the same reason: the driver $M$ translates the entire momentum distribution, and reflecting extreme bins accumulate the translated mass. The model spikes at the right time because $M$ spikes at the right time — the distribution mechanics are slaves to the driver. This is why the **forecast** (§6) is the binding constraint, not the hindcast.

---

## 6. Genuine 2022 out-of-sample forecast

Model frozen on pre-2022 data; $\widehat M$ from the mechanical nowcast only. The 2022 distribution forecast correctly anticipates the growth-to-value rotation's loser-lane buildup but underestimates its persistence — the same uniform-bin limitation seen in 2009. Forecast-skill (2008–2022 walk-forward): RMSE 0.074, correlation 0.50, directional 70%.

---

## 7. Paper-trading deployment protocol

### 7.1 Rules (frozen, zero lookahead)

At each month-end $t$ from 2008-01: (i) compute $\widehat M(t{+}1)$ mechanically; (ii) simulate one step to $\widehat{\mathbf{c}}(t{+}1)$; (iii) set tilt weights $\pi_b \propto \max(c_b + \gamma(\widehat c_b - c_b), 0)$, $\gamma = 1$ (additive — the earlier exponential form was numerically unstable and is discarded); (iv) within each bin, weight stocks by market cap; (v) **freeze** weights through $t{+}1$; (vi) deduct 5 bps one-way on rebalanced turnover, computed from correctly drifted prior weights (a timing bug in the first implementation used future returns in the drift — fixed).

Benchmarks: **SPY** buy-and-hold (dividends via Yahoo adjusted closes) [2], and a **true equal-weight historical-member portfolio** (equal weights over eligible point-in-time members, same rebalance/cost accounting — the earlier "EW" was equal-weight across bins and is replaced).

Monthly equity **and** frozen per-ticker holdings snapshots are saved to `~/workspace/goals/stock-market-dynamics/hidden_files/` as first-class deployment records.

### 7.2 Results (2008–2022 walk-forward)

| Portfolio | Total | CAGR | Ann. vol | Sharpe | Max DD |
|---|---|---:|---:|---:|---:|
| Tilt (paper) | +196% | 7.6% | 19.2% | 0.48 | −35.9% |
| SPY | +284% | 9.5% | 16.3% | 0.64 | −46.3% |
| Equal-weight members | +327% | 10.3% | 13.3% | 0.81 | −27.0% |

2022: tilt −18.7%, SPY −18.2%, EW −4.2%.

### 7.3 The tilt does not work — mechanism

The mass-flow tilt was tested to destruction. Correlation between predicted bin mass change $\Delta\widehat c_b$ and realized bin return: **0.009** overall (winner bins 0.031, loser bins 0.064, mid −0.099). A perfect-foresight variant (realized $M$, zero forecast error) achieves Sharpe 0.47 — no better. **The failure is in the tilt mechanism, not the forecast.**

Why: advection-driven mass flows are **reclassification**, not price pressure. When stocks crash they are reclassified into the loser bin — the bin "gains mass" while its constituents lose money. Overweighting predicted mass-gainers systematically overweights bins whose members just fell (loser inflow) or, symmetrically, chases winners after the move. The model's distribution mechanics predict *where capital sits*, not *which bins will pay*. A return-predictive tilt needs a return model the current framework does not have.

### 7.4 Verdict

**The paper record does not support real-capital deployment.** The distribution model captures conditional spikes and the mechanical driver nowcast has genuine skill, but the portfolio built on them trails both benchmarks on every risk-adjusted metric. The prerequisite is not met. The protocol, code, and snapshots remain in place so any future tilt can be judged by the same bar.

---

## 8. Limitations and next steps

1. **Return-predictive tilt.** The binding gap. Options: tilt on predicted bin returns via an auxiliary return model; use the distribution forecast for crash-risk overlays rather than cross-sectional tilts.
2. **Within-bin heterogeneity.** Uniform-bin advection under-resolves V-recoveries (2009, 2022 persistence). A depth coordinate or particle-level extension is the structural fix.
3. **Pre-2009 shares.** SEC XBRL starts mid-2009; early size/turnover bins lean on backfilled shares (flagged estimated).
4. **Membership reconstruction.** Community-sourced [1], not official S&P; change-date lags of days are possible.
5. **Empty cells.** Low-turnover/high-volatility bins are often empty; edges could be recalibrated on pre-sample data (kept fixed here for interpretability).
6. **Damper parameters.** Assumed from rev-3 episode tuning; a genuine calibration needs identified spike episodes in the training window.

---

## References

[1] fja05680/sp500 — *S&P 500 Historical Components & Changes (Updated).csv* (public GitHub reconstruction; 2,720 snapshots, 1996-01-02–2026-08-18). Used in §2.1 for point-in-time membership.
[2] Yahoo Finance — daily adjusted closes, corporate-action (split) history, and SPY total-return series. Used in §2.2 (prices), §2.3 (splits), §7.1 (SPY benchmark).
[3] SEC EDGAR — companyconcept XBRL, *SharesOutstanding* / *WeightedAverageNumberOfSharesOutstanding*. Used in §2.3 for point-in-time shares.
