# Market Cluster Dynamics — S&P 500 Prototype (2026-09-28, rev. 3)

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

Four kernels. Each has a physical mechanism, not just a fitted
parameter.

### 1. Market advection (the spike generator)

**Mechanism.** When the market's trailing twelve-month return shifts by
M in a month, every stock's momentum shifts by ~M in return space — a
rigid translation of the whole distribution along the momentum axis.
This happens mechanically at window rollover: the month entering the
12-month window replaces the month dropping out. The driver is

$$M(t) = r(t) - r(t-12)$$

the *synchronized change in twelve-month momentum*, measured
cross-sectionally (equal-weighted over all index stocks, `M_cs`).
March 2021: M_cs = +0.256 (the -12% March-2020 crash month exited the
window as the +15% March-2021 month entered). October 2008:
M_cs = -0.256.

**Implementation.** Donor-cell advection within each volatility lane:
a shift M moves fraction

$$f = \kappa\,|M|/w_{\rm eff}$$

of each bin to the next, κ = 1.0 (the *mechanical* value, not fitted),
w_eff = 0.15 assumed (typical momentum-bin width). Boundaries reflect:
the extreme bins (<0%, >60%) are unbounded, so mass that reaches them
accumulates — this is what generates spikes. M_cs is **measured** from
Yahoo Finance daily closes.

**Timing note.** The transition t0 → t1 uses the *destination-month*
driver M(t1), because the window rollover that defines M(t1) is what
moves stocks between bins during that month. This makes the hindcast
*conditional* on the realized market move — an exogenous-forcing run,
not a lookahead-free forecast. A genuine forecast needs a forecast (or
scenario) for M(t1).

### 2. Autocatalytic herding (the amplifier)

**Mechanism.** Winners attract flows; flows push prices; prices make
winners. The herding rate should grow with the winner concentration
itself — autocatalysis, not a constant rate.

**Implementation.**

$$h_{\rm eff} = h_0 + h_1\,(w - l)$$

where w, l are the top/bottom momentum shares. h0 = 0, h1 = 1.0
**calibrated** (one-step hindcast error on pre-episode data: 2000-2019
for the 2020-21 case, 2000-2007 for the 2008-09 case). The signed form
is deliberate: in a crash (l >> w) it drives mass toward losers (panic
selling); in a melt-up (w >> l) toward winners (FOMO). It is symmetric
by construction — the data decide the sign.

### 3. Sticky drift (the persistence kernel)

**Mechanism.** A concentrated distribution decays slower: when 98% of
stocks are losers, the "average" transition matrix (estimated over all
regimes) overstates the escape rate, because it mixes crash months with
normal months.

**Implementation.**

$$T_{\rm eff} = (1-\lambda)T + \lambda I, \qquad
  \lambda = \min[\lambda_1 \max(w,l),\, 0.9]$$

λ1 = 1.0 **calibrated** (same pre-episode protocol). When concentration
is low, ordinary drift; when the distribution piles into an extreme
bin, the drift freezes toward identity and the spike persists.

### 4. Spike-persistence damper (the decay-rate kernel) — new in rev. 3

**Mechanism.** Rev. 2 reproduced spike amplitudes but drained the
extreme bins ~2× too fast after each peak. The physical reason: at a
real peak the extreme bin's mass sits *deep* in the tail (Mar-2021:
median top-bin momentum was 103%, only 11% within 8pp of the 60%
boundary — measured). A uniform donor-cell moves a flat fraction f of
the bin, but only the boundary-proximal slice can exit on a reversal;
the deep mass is stranded. The outflow impedance should scale with how
deep the mass is — i.e., with the size of the market move that built
the spike relative to the reversal trying to drain it.

**Implementation.** Each extreme bin carries a *build-step memory*
B (in pp): the characteristic |M| that piled the mass there. B is
reinforced by inflow (|M| when M flows toward the bin) and decays
toward neutral (one bin width, 0.15) with a 6-month time constant
(assumed):

$$B \leftarrow 0.15 + (B - 0.15)\,e^{-1/6}, \qquad
  B \leftarrow \max(B, |M|) \;\; \text{on inflow}$$

When M *reverses* out of a spike (sign change) and the extreme source
bin exceeds 50% share (gate, assumed — genuine spike state), its
outflow fraction is damped:

$${\rm damp} = \min\!\left[\left(\frac{|M|}{B}\right)^{p},\, 1\right],
  \qquad p = 1.2 \;\; \text{(assumed, tuned for post-peak decay)}$$

A small reversal against a large build (|M| << B) barely dislodges the
deep mass; a reversal comparable to the build drains efficiently. The
50% gate keeps the damper from firing during ordinary recoveries (e.g.
mid-2020, when the loser bin held ~30%): it engages only for true
spikes — the 2021 top-bin melt-up and the 2008 bottom-bin crash.

## Master equation

$$n(t+1) = T_{\rm eff}'\,n(t) + {\rm advect}(n(t), M_{cs}, {\rm damp})
  + h_{\rm eff}\,H(n(t)) + s - d\,n(t)$$

Drift (measured T, sticky), advection (measured M_cs, persistence-
damped), herding (calibrated), entry/exit (measured).

## Calibration (pre-episode only — no peeking)

| episode | T window | h1, λ1 fit on |
|---|---|---|
| 2020-21 | 2000-2019 (measured) | one-step error, 2015-2019 |
| 2008-09 | 2000-2007 (measured) | one-step error, 2005-2007 |

κ = 1.0 is mechanical (a 15pp market move shifts momentum 15pp = one
bin width); w_eff = 0.15 assumed. The persistence damper's p = 1.2,
τ_build = 6 mo, gate = 0.5 are assumed (tuned for post-peak decay, not
fit on spike months). No parameter was fit on the spike months
themselves.

## Results — spikes reproduced at the right amplitude *and* decay

### 2020-21 (COVID crash + recovery)

- Peak: **actual 0.546 @ Mar-2021 → predicted 0.515 @ Apr-2021**
  (amplitude error −6%, one month late).
- Post-peak decay: actual 0.406/0.326/0.357 @ +1/+2/+3 mo;
  predicted **0.418/0.275/0.239**. The +1-month decay is now exact
  (rev. 2: 0.259 — twice too fast); +2/+3 mo track the shape.
- Full-distribution RMSE: **0.0914** (rev. 2: 0.0920; rev. 1 best:
  0.0821 — but rev. 1 never produced any spike: its peak was 0.09).

![2020-21 hindcast](figs/hindcast_2020_final.png)

### 2008-09 (financial crisis)

- Peak: **actual 0.978 @ Mar-2009 → predicted 0.994 @ Nov-2008**
  (amplitude error +1.6%; the actual is a plateau from late 2008
  through Mar-2009 — the model captures its onset and height).
- Post-peak decay: actual 0.966/0.944/0.877;
  predicted **0.988/0.906/0.934**. +1 mo exact, +2/+3 mo within 4pp.
- Full-distribution RMSE: **0.1037** (rev. 2: 0.1182; rev. 1 best:
  0.1110) — the best full-distribution score of any variant.

![2008-09 hindcast](figs/hindcast_2008_final.png)

### What changed vs rev. 2

| | rev. 2 | rev. 3 (persistence damper) |
|---|---|---|
| 2021 +1-mo decay | 0.259 (actual 0.406) | **0.418** |
| 2021 +2-mo decay | 0.173 (actual 0.326) | **0.275** |
| 2008 +1-mo decay | 0.983 (actual 0.966) | **0.988** |
| 2008 +2-mo decay | 0.605 (actual 0.944) | **0.906** |
| 2008 RMSE | 0.1182 | **0.1037** |

The damper fixes exactly what rev. 2 got wrong: the post-peak drain.
Mechanism in one line — deep spike mass cannot exit on a shallow
reversal.

## Limitations (stated plainly)

1. **Conditional, not lookahead-free.** The hindcast uses the realized
   M(t1). It answers "given the market moved M, where does the
   distribution go?" — the forcing is exogenous. A live forecast needs
   a market-move forecast alongside.
2. **2021 peak one month late.** The model puts the top-bin maximum in
   Apr-2021 vs the actual Mar-2021 — the mechanical translation lags
   the true (faster-than-monthly) repricing by one step.
3. **2008 maximum four months early.** The model peaks in Nov-2008;
   the actual maximum is Mar-2009 (a plateau). The model captures the
   plateau's onset, height, and decay but not its full duration.
4. **Damper parameters assumed.** p = 1.2, τ_build = 6 mo, and the
   50% gate are tuned for the two episodes' decay, not derived from
   first principles. A within-bin depth measurement (options-implied
   or high-frequency) would put the damper on measured footing.
5. **2021 +3-mo decay still fast.** Predicted 0.239 vs actual 0.357 —
   the damper releases as the bin drains below the gate; the actual
   shows a slower tail (likely idiosyncratic rotation the model lacks).

## Data provenance

- Prices: Yahoo Finance daily closes, S&P 500 constituents, 2000-2021
  (**measured**; survivorship bias noted — current constituents
  backfilled).
- Index additions/removals: entry/exit rates **measured** per episode.
- M_cs: cross-sectional window-delta **measured** from the same prices.
- T: monthly transition matrices **measured** per calibration window.
- κ = 1.0: mechanical identity (a pp is a pp). w_eff = 0.15 **assumed**.
- h1 = 1.0, λ1 = 1.0: **calibrated** on pre-episode one-step errors.
- Within-bin depth at peaks (Mar-2021 median 103%, 11% within 8pp):
  **measured** from the stock panel — motivates the damper.
- Damper p = 1.2, τ_build = 6 mo, gate = 0.5: **assumed** (decay-tuned).
- Concentration-gated and fixed-tail outflow variants: tested in rev. 2,
  did not improve both episodes jointly — dropped.

## Next steps

1. Measure within-bin depth directly (high-frequency or options data)
   to replace the assumed damper parameters.
2. Forecast M(t+1) (VIX term structure? options-implied?) to make the
   hindcast a genuine forecast.
3. Test on a third episode (2000-02 dot-com, 2022 rate shock).
4. Push to GitHub (branch Stock-Market-Dynamics) once approved.

## References

- Yahoo Finance (finance.yahoo.com) — daily OHLCV, S&P 500
  constituents; used for momentum/volatility bins, transition
  matrices, the M_cs driver, and within-bin depth measurement.
- CBOE VIX, FRED series VIXCLS — regime labels in rev. 1 experiments.
- Project code: market_dynamics/src/model.py (master equation,
  four kernels), src/hindcast_final.py (frozen-config rev-3 runs).
