"""Generate monthly frozen per-ticker holdings snapshots for the paper track record."""
import os, json
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from model4 import N_BIN, assign_bins, cap_state, simulate, market_driver, mech_forecast
from paper_trade import tilt_weights

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.expanduser('~/workspace/goals/stock-market-dynamics/hidden_files/holdings')

os.makedirs(OUT, exist_ok=True)
p = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'), parse_dates=['month'], low_memory=False)
p['bin'] = assign_bins(p)
P = json.load(open(os.path.join(DATA, 'params_rev4.json')))
Tcalm = np.array(P['_Tcalm']); emix = np.array(P['_emix'])
M_hist = market_driver(p)
spy = pd.read_csv(os.path.join(DATA, 'raw', 'SPY.csv'), parse_dates=['date']).sort_values('date')
spy = spy.set_index('date')['adjclose']
months = sorted(p[p['month'] >= '2008-01-01']['month'].unique())
spy_ret = {}
for m in months:
    d = spy[spy.index <= m]
    if len(d) >= 22:
        spy_ret[m.strftime('%Y-%m')] = float(d.iloc[-1] / d.iloc[-22] - 1)

for k in range(len(months) - 1):
    m0 = months[k]
    ym0 = m0.strftime('%Y-%m')
    d0 = p[p['month'] == m0]
    c, _ = cap_state(d0)
    Mhat, _ = mech_forecast(spy_ret, ym0)
    ch = simulate(c, Tcalm, [Mhat], P, 1, erate=P['_erate'], emix=emix,
                  p_damp=P['p_damp'], tau_build=P['tau_build'], gate=P['gate'])[1]
    pi = tilt_weights(c, ch, 1.0)
    d0 = d0[d0['bin'] >= 0].copy()
    bin_mcap = d0.groupby('bin')['mcap'].sum()
    d0['w'] = d0.apply(lambda r: pi[int(r['bin'])] * r['mcap'] / bin_mcap[int(r['bin'])], axis=1)
    d0['w'] /= d0['w'].sum()
    out = d0[['ticker', 'bin', 'w', 'mcap']].copy()
    out['origin'] = ym0
    out.to_csv(os.path.join(OUT, f'holdings_{ym0}.csv'), index=False)
    if k % 40 == 0:
        print(f'wrote {ym0}', flush=True)
print('done', len(months) - 1, 'holdings files')
