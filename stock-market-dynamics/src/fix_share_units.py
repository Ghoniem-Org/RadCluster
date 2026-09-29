"""REV-5 data fix: unit errors in shares via segment-level absolute anchoring.

clean_shares.py only fixes 1000x/1e6x JUMPS vs the previous filing, and only
within narrow bands (900-1100x, 9e5-1.1e6x). Residual problems found 2026-09-28:
  (A) tickers whose filings are CONSISTENTLY in the wrong unit (SRE backfill
      1e6x -> $5,000T mcap, dominated all cap-weighted states);
  (B) unit changes WITHIN the filed era missed because a genuine share change
      rides on top (HAL: 9e8 -> 880; C: 2.3e7 thousands -> 2.8e10 ones,
      ratio 1246x, outside the 900-1100 band);
  (C) garbage placeholders (shares = 1.0 / 100.0: FOX/FOXA, BKR).

Method, per ticker:
  1. sh_adj = shares_use * S_f; split into segments at month-to-month jumps
     > 20x or < 1/20x (split-adjusted shares cannot move like that).
  2. Per segment (>=3 months): med_mcap and med daily turnover.
     Sane = mcap in [$200M, $4T] AND turnover in [0.02%, 100%].
     - f=1 both sane -> keep.
     - else if some f in {1e3,1e-3,1e6,1e-6} makes both sane -> apply the one
       closest to $30B mcap (this catches C seg0: $248M mcap but 286% turnover).
     - else if f=1 mcap-sane -> keep (turnover can be genuinely extreme).
     - else if some f makes mcap sane -> apply closest-to-$30B.
     - else drop the segment (NaN).
     Segments < 3 months -> dropped (transient spikes).
  3. Recompute mcap exactly; turnover rescaled inversely (volume/denom).
Backs up original to data/panel_v2_preunitfix.csv. Idempotent.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

LO, HI = 2e8, 4e12          # sane S&P-500 mcap bounds ($)
TO_LO, TO_HI = 0.0002, 1.0  # sane average-daily turnover bounds
TARGET = 3e10
FACTORS = [1e3, 1e-3, 1e6, 1e-6]
JUMP = 20.0


def pick_factor(med_mcap, med_to):
    def mcap_ok(f):
        return np.isfinite(med_mcap) and LO <= med_mcap * f <= HI
    def to_ok(f):
        return (not np.isfinite(med_to)) or (TO_LO <= med_to / f <= TO_HI)
    if mcap_ok(1.0) and to_ok(1.0):
        return 1.0, 'keep'
    both = [f for f in FACTORS if mcap_ok(f) and to_ok(f)]
    if both:
        best = min(both, key=lambda f: abs(np.log10(med_mcap * f / TARGET)))
        return best, 'both_sane'
    if mcap_ok(1.0):
        return 1.0, 'keep_mcap_ok'
    mcaponly = [f for f in FACTORS if mcap_ok(f)]
    if mcaponly:
        best = min(mcaponly, key=lambda f: abs(np.log10(med_mcap * f / TARGET)))
        return best, 'mcap_only'
    return None, 'drop'


def main():
    pin = os.path.join(DATA, 'panel_v2.csv')
    bkp = os.path.join(DATA, 'panel_v2_preunitfix.csv')
    if not os.path.exists(bkp):
        pd.read_csv(pin, low_memory=False).to_csv(bkp, index=False)
        print('backed up original panel', flush=True)
    p = pd.read_csv(pin, parse_dates=['month'], low_memory=False)
    p = p.sort_values(['ticker', 'month']).reset_index(drop=True)
    p['sh_adj'] = p['shares_use'] * p['S_f']
    p['mcap0'] = p['close_raw'] * p['sh_adj']

    # segment at large jumps
    p['prev'] = p.groupby('ticker')['sh_adj'].shift(1)
    ratio = p['sh_adj'] / p['prev']
    p['seg'] = ((ratio > JUMP) | (ratio < 1 / JUMP)).groupby(p['ticker']).cumsum()

    stats = {'keep': 0, 'both_sane': 0, 'keep_mcap_ok': 0, 'mcap_only': 0, 'drop': 0,
             'short_drop': 0}
    detail = []
    for (t, s), g in p.groupby(['ticker', 'seg'], sort=False):
        idx = g.index
        if len(g) < 3:
            p.loc[idx, 'shares_use'] = np.nan
            stats['short_drop'] += 1
            continue
        med_mcap = float(g['mcap0'].median())
        med_to = float(g['turnover'].median()) if g['turnover'].notna().any() else np.nan
        f, why = pick_factor(med_mcap, med_to)
        stats[why] += 1
        if f is None:
            p.loc[idx, 'shares_use'] = np.nan
            detail.append(f'DROP {t} seg{s} n={len(g)} mcap={med_mcap:.2e} to={med_to:.2e}')
        elif f != 1.0:
            p.loc[idx, 'shares_use'] = p.loc[idx, 'shares_use'] * f
            detail.append(f'{t} seg{s} n={len(g)} x{f:g} ({why}) '
                          f'mcap {med_mcap:.2e}->{med_mcap*f:.2e}')
    print('segment outcomes:', stats, flush=True)
    for d in detail:
        print('  ', d, flush=True)

    # recompute
    p['sh_adj'] = p['shares_use'] * p['S_f']
    p['mcap'] = p['close_raw'] * p['sh_adj']
    b = pd.read_csv(bkp, parse_dates=['month'], low_memory=False)
    b = b.sort_values(['ticker', 'month']).reset_index(drop=True)
    b['den_old'] = b['shares_use'] * b['S_f']
    mrg = p[['ticker', 'month']].merge(
        b[['ticker', 'month', 'den_old', 'turnover']].rename(columns={'turnover': 'to_old'}),
        on=['ticker', 'month'], how='left')
    f = p['sh_adj'].values / mrg['den_old'].values
    ok = np.isfinite(f) & (f > 0) & np.isfinite(mrg['to_old'].values)
    p['turnover'] = np.nan
    p.loc[ok, 'turnover'] = (mrg.loc[ok, 'to_old'].values / f[ok])

    p['prov'] = p['prov'].astype(str) + '+unitfix2'

    # verify
    from model4 import assign_bins
    p['bin'] = assign_bins(p)
    q = p[p['mcap'].notna()]
    print('post-fix: max mcap %.2e; #>2T pre-2018: %d' % (
        q['mcap'].max(), int(((q['mcap'] > 2e12) & (q['month'] < '2018-01-01')).sum())),
        flush=True)
    print('tickers w/ median mcap outside [$200M,$4T]:',
          int((((q.groupby('ticker')['mcap'].median() < LO) |
                (q.groupby('ticker')['mcap'].median() > HI))).sum()), flush=True)
    a = p[(p['ticker'] == 'AAPL') & (p['month'] == '2007-12-31')]
    if len(a):
        print('AAPL 2007-12 mcap=$%.1fB (expect ~176)' % (a['mcap'].iloc[0] / 1e9), flush=True)
    s_ = p[(p['ticker'] == 'SRE') & (p['month'] == '2001-01-31')]
    if len(s_):
        print('SRE 2001-01 mcap=$%.1fB (expect ~5)' % (s_['mcap'].iloc[0] / 1e9), flush=True)
    for ym, ref in [('2007-12-31', 13e12), ('2021-12-31', 40e12)]:
        tot = p[(p['month'] == ym) & (p['bin'] >= 0)]['mcap'].sum()
        print(f'total binnable cap {ym}: ${tot/1e12:.1f}T (ref ~${ref/1e12:.0f}T)', flush=True)

    p = p.drop(columns=['sh_adj', 'mcap0', 'prev', 'seg', 'bin'])
    p.to_csv(pin, index=False)
    print('wrote fixed panel', flush=True)


if __name__ == '__main__':
    main()
