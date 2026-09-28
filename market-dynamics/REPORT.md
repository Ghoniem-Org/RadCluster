# Market Cluster Dynamics — S&P 500 Prototype (2026-09-28, rev. 2)

## The problem with rev. 1

The stationary-drift model (rev. 1, 2026-09-27) reproduced calm-market
distributions (RMSE 0.045) but could not generate the observed spikes:
March 2021 saw 55% of stocks pile into the top momentum bin (>60%
twelve-month return); late 2008 saw 98% pile into the bottom bin (<0%).
A monthly Markov transition matrix mean-reverts too fast to build or
hold such concentrations. Herding (bilinear flux toward
above-average-momentum bins) was a weak amplifier at best — it can only
redistribute mass the drift already moves.

## New reaction kernels (the mechanism first)

Three kernels were added. Each has a physical mechanism, not just a
fitted parameter.

### 1. Market advection (the spike generator)

**Mechanism.** When the market's trailing twelve-month return shifts by
M in a month, every stock's momentum shifts by ~M in return space — a
rigid translation of the whole distribution along the momentum axis.
This happens mechanically at window rollover: the month entering the
12-month window replaces the month dropping out. The driver is

    M(t) = r(t) - r(t-12)

the *synchronized change in twelve-month momentum*, measured
cross-sectionally (equal-weighted over all index stocks, `M_cs`).
March 2021: M_cs = +0.256 (the -12% March-2020 crash month exited the
window as the +15% March-2021 month entered). October 2008:
M_cs = -0.256.

**Implementation.** Donor-cell advection within each volatility lane:
a shift M moves fraction f = kappa*|M|/w_eff of each bin to the next,
kappa = 1.0 (the *mechanical* value, not fitted), w_eff = 0.15 assumed
(typical momentum-bin width). Boundaries reflect: the extreme bins
(<0%, >60%) are unbounded, so mass that reaches them accumulates —
this is what generates spikes. M_cs is **measured** from Yahoo Finance
daily closes.

**Timing note.** The transition t0 -> t1 uses the *destination-month*
driver M(t1), because the window rollover that defines M(t1) is what
moves stocks between bins during that month. This makes the hindcast
*conditional* on the realized market move — an exogenous-forcing run,
not a lookahead-free forecast. A genuine forecast needs a forecast (or
scenario) for M(t1).

### 2. Autocatalytic herding (the amplifier)

**Mechanism.** Winners attract flows; flows push prices; prices make
winners. The herding rate should grow with the winner concentration
itself — autocatalysis, not a constant rate.

**Implementation.** h_eff = h0 + h1*(w - l), where w, l are the
top/bottom momentum shares. h0 = 0, h1 = 1.0 **calibrated** (one-step
hindcast error on pre-episode data: 2000-2019 for the 2020-21 case,
2000-2007 for the 2008-09 case). The signed form is deliberate: in a
crash (l >> w) it drives mass toward losers (panic selling); in a
melt-up (w >> l) toward winners (FOMO). It is symmetric by
construction — the data decide the sign.

### 3. Sticky drift (the persistence kernel)

**Mechanism.** A concentrated distribution decays slower: when 98% of
stocks are losers, the "average" transition matrix (estimated over all
regimes) overstates the escape rate, because it mixes crash months with
normal months.

**Implementation.** T_eff = (1-lambda)*T + lambda*I, with
lambda = min(lam1*max(w,l), 0.9), lam1 = 1.0 **calibrated** (same
pre-episode protocol). When concentration is low, ordinary drift;
when the distribution piles into an extreme bin, the drift freezes
toward identity and the spike persists.

## Master equation

    n(t+1) = T_eff' n(t) + advect(n(t), M_cs) + h_eff * H(n(t)) + s - d*n(t)

Drift (measured T, sticky), advection (measured M_cs), herding
(calibrated), entry/exit (measured).

## Calibration (pre-episode only — no peeking)

| episode | T window | h1, lam1 fit on |
|---|---|---|
| 2020-21 | 2000-2019 (measured) | one-step error, 2015-2019 |
| 2008-09 | 2000-2007 (measured) | one-step error, 2005-2007 |

kappa = 1.0 is mechanical (a 15pp market move shifts momentum 15pp =
one bin width); w_eff = 0.15 assumed. No parameter was fit on the
spike months themselves.

## Results — the spikes are reproduced

### 2020-21 (COVID crash + recovery)

- Peak: **actual 0.546 @ Mar-2021 → predicted 0.579 @ Apr-2021**
  (amplitude error +6%, one month late).
- Full-distribution RMSE: **0.0915** (rev. 1 best: 0.0821 — slightly
  worse on RMSE, but rev. 1 never produced any spike at all: its peak
  was 0.09 vs the actual 0.55).
- Decay after peak: actual 0.41/0.33/0.36 @ +1/+2/+3 mo;
  predicted 0.26/0.17/0.16. **Too fast** — see limitations.

### 2008-09 (financial crisis)

- Peak: **actual 0.978 @ Mar-2009 → predicted 0.983 @ Nov-2008**
  (amplitude error +0.5%; the actual is a plateau from late 2008
  through Mar-2009 — the model captures its onset and height).
- Full-distribution RMSE: **0.1182** (rev. 1 best: 0.1110).
- Decay after peak: actual 0.97/0.94/0.88; predicted 0.98/0.61/0.72.
  +1 mo is exact; then too fast — see limitations.

### What changed vs rev. 1

| | rev. 1 (regime T, h=0.2) | rev. 2 (new kernels) |
|---|---|---|
| 2021 peak | 0.09 (no spike) | **0.58** (actual 0.55) |
| 2008 peak | 0.30 (no spike) | **0.98** (actual 0.98) |

The advection kernel is doing the work: the March-2021 spike is the
+25.6pp window-delta mechanically translating the distribution into
the unbounded top bin; the 2008 spike is the -25.6pp October delta
translating it into the bottom bin. Herding amplifies (+0.16 on the
2008 peak at h1=1 vs h1=0) and sticky drift holds the plateaus.

## Limitations (stated plainly)

1. **Decay too fast.** After each peak the model drains the extreme bin
   ~2x faster than observed. Mechanism: the donor-cell moves a uniform
   fraction f of the bin, but at a real peak the bin's mass sits deep
   in the tail (Mar-2021: median top-bin momentum was 103%, only 11%
   within 8pp of the 60% boundary — measured). A 15-bin distribution
   cannot resolve within-bin depth without memory of how the mass got
   there. Fixing this needs either sub-binning the extremes or an
   age-structured extreme bin.
2. **Conditional, not lookahead-free.** The hindcast uses the realized
   M(t1). It answers "given the market moved M, where does the
   distribution go?" — the forcing is exogenous. A live forecast needs
   a market-move forecast alongside.
3. **2021 peak one month late.** The model puts the top-bin maximum in
   Apr-2021 vs the actual Mar-2021 — the mechanical translation lags
   the true (faster-than-monthly) repricing by one step.
4. **2008 plateau too short.** The model holds 0.98 for ~2 months;
   the actual plateau lasted ~5. Same within-bin-depth cause as (1).

## Data provenance

- Prices: Yahoo Finance daily closes, S&P 500 constituents, 2000-2021
  (**measured**; survivorship bias noted — current constituents
  backfilled).
- Index additions/removals: entry/exit rates **measured** per episode.
- M_cs: cross-sectional window-delta **measured** from the same prices.
- T: monthly transition matrices **measured** per calibration window.
- kappa=1.0: mechanical identity (a pp is a pp). w_eff=0.15 **assumed**.
- h1=1.0, lam1=1.0: **calibrated** on pre-episode one-step errors.
- Spike-threshold/tail experiments (concentration-dependent extreme-bin
  outflow): tested, did not improve both episodes jointly — dropped.

## Next steps

1. Sub-bin or age-structure the extreme bins to fix the decay rate.
2. Forecast M(t+1) (VIX term structure? options-implied?) to make the
   hindcast a genuine forecast.
3. Test on a third episode (2000-02 dot-com, 2022 rate shock).
4. Push to GitHub (branch Stock-Market-Dynamics) once approved.

## References

- Yahoo Finance (finance.yahoo.com) — daily OHLCV, S&P 500
  constituents; used for momentum/volatility bins, transition
  matrices, and the M_cs driver.
- CBOE VIX, FRED series VIXCLS — regime labels in rev. 1 experiments.
- Project code: market_dynamics/src/model.py (master equation),
  src/hindcast_final.py (frozen-config runs).
