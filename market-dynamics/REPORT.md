# Market Cluster Dynamics — S&P 500 Prototype (2026-09-27)

## Setup
- **Vertices (15):** momentum × volatility regime bins. Momentum = 12–1 month
  return on adjclose (absolute bins: <0%, 0–15%, 15–30%, 30–60%, >60%);
  volatility = trailing-60d annualized stdev (absolute terciles: <20%, 20–32%, >32%).
- **State:** n_k(t) = fraction of index stocks in bin k.
- **Edges:** (1) drift — empirical monthly Markov transitions T (measured,
  2015–2019); (2) herding flux H_k = n_k(q_k − q̄)/4, q_k = momentum quintile
  of bin k (zero-sum, the bilinear "mating" analog); (3) index entry/exit
  source/sink (measured: 0.0008/mo in, 0.0000/mo out).
- **Master equation:** n(t+1) = T′n(t) + h·H(n(t)) + s − d·n(t).
- **Data:** 460 S&P 500 tickers, daily 2015–2021, Yahoo Finance (measured).
  Survivorship bias noted: current constituents backfilled.

## Calibration (2016–2019, calm market)
| h | 0.0 | 0.2 | 0.5 | 1.0 | 2.0 | 3.0 |
|---|-----|-----|-----|-----|-----|-----|
| RMSE | **0.0451** | 0.0509 | 0.0689 | 0.0974 | 0.1293 | 0.1444 |

Best h = 0. Drift alone explains the calm-period distribution; herding adds nothing.

## Stress test (2020–2021, COVID regime shift)
| h | 0.0 | 0.2 | 0.5 | 1.0 |
|---|-----|-----|-----|-----|
| RMSE | 0.0934 | **0.0917** | 0.0920 | 0.0984 |

Best h = 0.2 — herding helps slightly, but the stationary drift model cannot
generate the March-2021 concentration (55% of stocks >60% twelve-month winners
off the crash bottom). Dec-2021 extreme-winner share: actual 8.8%, model 5.1%.

## Mechanism reading
The herding term is a ratchet *amplifier*, not a regime generator: it can only
redistribute mass the drift already moves. A stationary T estimated in calm
markets has no knowledge of crash-recovery dynamics, so no h repairs the
hindcast. This mirrors the genealogy finding that the deepest-inheritance
ratchet needs the right *rule*, not just a stronger parameter.

## Regime-dependent drift (VIX-split T)
Tested 2026-09-27: T_calm (VIX≤20, 20,596 transitions) vs T_stress (VIX>20,
1,858 transitions), both measured 2015–2019; hindcast switches T by actual
monthly VIX (perfect-foresight regime label).
- Stationary T, h=0: RMSE 0.0934 | Regime T, h=0: RMSE 0.1077 (worse)
- 3-regime (calm/crash/recovery): RMSE 0.1074 (worse)
- Mechanism: T_stress's stationary distribution puts **0%** in the top momentum
  bin. The 2015–19 "stress" months were all *drawdowns* (flow into losers);
  2021 was a *recovery* (flow into winners off the crash bottom). VIX level
  conflates two opposite flows — crash ≠ recovery. Only 5 crash + 3 recovery
  months exist in 2015–19, too few for a clean stress operator.
- Next step: extend history to 2000/2008 for genuine stress calibration, or a
  parametric T(VIX) instead of hard regime splits.

## Files
- `src/fetch_data.py`, `src/build_regimes.py`, `src/model.py`, `src/hindcast.py`
- `data/regimes.csv` (35k rows), `outputs/T_drift.npy`, `outputs/calibration_*.csv`
- `figs/rag.png`, `figs/calibration.png`, `figs/hindcast_topq.png`, `figs/dist_*.png`
