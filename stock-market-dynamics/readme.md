# Stock Market Cluster Dynamics (prototype)

S&P 500 market-regime cluster dynamics, ported from the RadCluster/genealogy-dynamics
master-equation framework: dc/dt = P + SJ − Dc.

- **Vertices:** 15 regime bins = momentum (12–1m return, absolute bins) × volatility (trailing-60d, terciles)
- **State:** fraction of index stocks per bin
- **Edges:** measured drift (Markov transitions), herding flux (bilinear winner-chasing term, rate h), index entry/exit
- **Validation:** 2020–21 COVID concentration hindcast; regime-dependent (VIX-split) drift tested

See REPORT.md for equations, calibration results, and mechanism findings.
Data: Yahoo Finance daily OHLCV (measured), 460 tickers 2015–2021. Raw cache excluded (90MB);
refetch with `python3 src/fetch_data.py`. Note survivorship bias: current constituents backfilled.

Pipeline: `src/fetch_data.py` → `src/build_regimes.py` → `src/hindcast.py` (+ `src/hindcast_regime.py`)
