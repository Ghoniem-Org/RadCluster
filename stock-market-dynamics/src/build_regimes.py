"""Build monthly regime panel: momentum quintile x volatility tercile -> 15 bins.

For each ticker and each month-end (2015-12 through 2024-12):
  momentum = 12-1 month total return on adjclose (measured)
  volatility = annualized stdev of daily log returns, trailing 60 trading days (measured)
Bins are cross-sectional: momentum quintiles x vol terciles, recomputed each month.
Output: data/regimes.csv  (date, ticker, mom, vol, mom_q, vol_t, bin)
"""
import os, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RAW = os.path.join(DATA, 'raw')

def load_prices():
    frames = []
    for fn in glob.glob(os.path.join(RAW, '*.csv')):
        t = os.path.basename(fn)[:-4]
        try:
            df = pd.read_csv(fn, parse_dates=['date'])
        except Exception:
            continue
        df = df.dropna(subset=['adjclose']).sort_values('date')
        if len(df) < 300:
            continue
        df['ticker'] = t
        frames.append(df[['date', 'ticker', 'adjclose']])
    px = pd.concat(frames, ignore_index=True)
    px = px.pivot(index='date', columns='ticker', values='adjclose').sort_index()
    return px

def main():
    px = load_prices()
    print(f"price panel: {px.shape}", flush=True)
    dates = px.index
    # month-ends present in data
    me = px.resample('M').last().dropna(how='all').index
    me = [d for d in me if d >= pd.Timestamp('2000-12-31')]
    rows = []
    logp = np.log(px)
    for d in me:
        # momentum: 12-1 month return
        try:
            p_now_1m = px.loc[:d].iloc[-22]      # ~1 month ago
            p_12m = px.loc[:d].iloc[-252]        # ~12 months ago
        except IndexError:
            continue
        mom = p_now_1m / p_12m - 1.0
        # volatility: trailing 60 trading days of daily log returns
        rets = logp.loc[:d].diff().iloc[-60:]
        vol = rets.std() * np.sqrt(252)
        valid = mom.dropna().index.intersection(vol.dropna().index)
        if len(valid) < 100:
            continue
        m = mom[valid]; v = vol[valid]
        # absolute bins (cluster-dynamics style: fixed states, not ranks)
        mom_bins = [-9, 0.0, 0.15, 0.30, 0.60, 9]
        vol_bins = [0.0, 0.20, 0.32, 9.0]
        mom_q = pd.cut(m, mom_bins, labels=False)
        vol_t = pd.cut(v, vol_bins, labels=False)
        ok = mom_q.notna() & vol_t.notna()
        for t in valid[ok]:
            mq, vt = int(mom_q[t]), int(vol_t[t])
            rows.append((d.date().isoformat(), t, float(m[t]), float(v[t]),
                         mq, vt, mq * 3 + vt))
    reg = pd.DataFrame(rows, columns=['date', 'ticker', 'mom', 'vol', 'mom_q', 'vol_t', 'bin'])
    reg.to_csv(os.path.join(DATA, 'regimes.csv'), index=False)
    print(f"regimes: {len(reg)} rows, {reg.date.nunique()} months, bins {sorted(reg.bin.unique())}", flush=True)

if __name__ == '__main__':
    main()
