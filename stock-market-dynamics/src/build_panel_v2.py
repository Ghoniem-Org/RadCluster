"""Rebuild panel shares/mcap/turnover with correct split adjustment.
mcap(t) = close_splitadj(t) * shares_filed * S(f),  S(f) = cum splits filing->present
turnover(t) = volume_yahoo(t) / (shares_filed * S(f))   [volume is split-adj]
Pre-2009: backfill from earliest filing (estimated, buyback error not corrected).
Writes data/panel_v2.csv
"""
import os
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

def main():
    p = pd.read_csv(os.path.join(DATA, 'panel.csv'), parse_dates=['month'], low_memory=False)
    sh = pd.read_csv(os.path.join(DATA, 'shares_clean.csv'), parse_dates=['end', 'filed'])
    sp = pd.read_csv(os.path.join(DATA, 'splits_yahoo.csv'), parse_dates=['date'])
    print('panel', len(p), 'shares', len(sh), 'splits', len(sp), flush=True)

    # S(d): cumulative split factor from date d to present, per ticker
    split_by_t = {}
    for t, g in sp.groupby('ticker'):
        g = g.sort_values('date')
        # for filing date f: S(f) = prod(ratio for split_date > f)
        dates = g['date'].values
        ratios = g['ratio'].values
        # suffix products
        suf = np.cumprod(ratios[::-1])[::-1]
        split_by_t[t] = (dates, suf, ratios)

    def S_of(ticker, fdate):
        """Cumulative splits with split_date > fdate."""
        if pd.isna(fdate):
            return 1.0
        if ticker not in split_by_t:
            return 1.0
        dates, suf, ratios = split_by_t[ticker]
        fdate = np.datetime64(pd.Timestamp(fdate).date())
        i = np.searchsorted(dates, fdate, side='right')
        return float(suf[i]) if i < len(suf) else 1.0

    sh = sh.sort_values(['ticker', 'filed'])
    sh_by_t = {t: g for t, g in sh.groupby('ticker')}

    # vectorized per ticker-month via merge_asof per ticker
    p = p.sort_values(['px_ticker', 'month']).reset_index(drop=True)
    p['shares_use'] = np.nan
    p['filed_use'] = pd.NaT
    p['backfilled'] = False
    for t, g in p.groupby('px_ticker', sort=False):
        sg = sh_by_t.get(t)
        idx = g.index.values
        months = g['month'].values
        if sg is None or len(sg) == 0:
            continue
        filed = sg['filed'].values
        svals = sg['shares'].values.astype(float)
        # PIT: latest filed <= month
        j = np.searchsorted(filed, months, side='right') - 1
        has = j >= 0
        p.loc[idx[has], 'shares_use'] = svals[j[has]]
        p.loc[idx[has], 'filed_use'] = filed[j[has]]
        # backfill pre-first-filing with earliest
        nb = ~has
        if nb.any():
            p.loc[idx[nb], 'shares_use'] = svals[0]
            p.loc[idx[nb], 'filed_use'] = filed[0]
            p.loc[idx[nb], 'backfilled'] = True

    # S(f) per row
    p['S_f'] = [S_of(t, f) for t, f in zip(p['px_ticker'], p['filed_use'])]
    denom = p['shares_use'] * p['S_f']
    p['mcap'] = p['close_raw'] * denom
    # turnover: 21d mean volume / denom ; volume col is month-end daily volume,
    # need trailing mean -> approximate with month-end volume / denom? No:
    # recompute properly below per ticker.
    p['turnover'] = np.nan

    # proper trailing-21d mean volume turnover per ticker
    px_cache = {}
    import glob
    for t, g in p.groupby('px_ticker', sort=False):
        idx = g.index.values
        dnom = (g['shares_use'] * g['S_f']).values
        f = os.path.join(DATA, 'raw', f'{t}.csv')
        try:
            df = pd.read_csv(f, parse_dates=['date']).sort_values('date')
        except Exception:
            continue
        df = df.set_index('date')
        for k, (ii, m, dn) in enumerate(zip(idx, g['month'], dnom)):
            if not np.isfinite(dn) or dn <= 0:
                continue
            mend = pd.Timestamp(m)
            d = df[df.index <= mend]
            if len(d) < 5:
                continue
            vv = d['volume'].iloc[-21:].values
            vv = vv[np.isfinite(vv)]
            if len(vv):
                p.loc[ii, 'turnover'] = np.mean(vv) / dn

    p['prov'] = np.where(p['backfilled'], 'shares_backfill_estimated',
                 np.where(p['shares_use'].notna(), 'shares_pit_measured', 'shares_missing'))
    p.to_csv(os.path.join(DATA, 'panel_v2.csv'), index=False)
    print('wrote panel_v2', flush=True)
    ok = p[p['has_price'] == 1]
    print('mcap coverage:', ok['mcap'].notna().mean().round(3), flush=True)
    print('turnover coverage:', ok['turnover'].notna().mean().round(3), flush=True)
    # sanity: AAPL mcap Dec 2007 and Dec 2022
    for m in ['2007-12-31', '2022-12-30']:
        r = ok[(ok['ticker'] == 'AAPL') & (ok['month'] == m)]
        if len(r):
            print(m, 'AAPL mcap=$%.1fB backfilled=%s' % (r['mcap'].iloc[0]/1e9, r['backfilled'].iloc[0]), flush=True)

if __name__ == '__main__':
    main()
