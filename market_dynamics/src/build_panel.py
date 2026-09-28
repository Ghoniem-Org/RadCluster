"""Build monthly member panel 2000-2022 for rev-4.
Output: data/panel.csv with per (month-end, member):
  features (mom, vol, turnover, mcap), bin coords, provenance flags.
Labels: measured / estimated per field (see PROV).
Point-in-time discipline: only shares with filed <= month-end are used.
"""
import csv, glob, os
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

ALIAS = {'BK': 'BNY', 'BLL': 'BALL', 'ABC': 'COR', 'ANTM': 'ELV', 'FB': 'META',
         'DISCA': 'WBD', 'DISCK': 'WBD', 'ACE': 'CB', 'DLPH': 'APTV', 'TMK': 'GL',
         'STI': 'TFC', 'BBT': 'TFC', 'UTX': 'RTX', 'RTN': 'RTX', 'PX': 'LIN'}

def load_prices():
    px = {}
    for f in glob.glob(os.path.join(DATA, 'raw', '*.csv')):
        t = os.path.basename(f)[:-4]
        try:
            df = pd.read_csv(f, parse_dates=['date'])
            df = df.sort_values('date').drop_duplicates('date')
            df['ticker'] = t
            px[t] = df
        except Exception:
            pass
    return px

def load_membership():
    snaps = []
    with open('/tmp/sp500_changes.csv') as f:
        for r in csv.DictReader(f):
            if '1999-01-01' <= r['date'] <= '2022-12-31':
                snaps.append((r['date'], set(x.replace('.', '-') for x in r['tickers'].split(','))))
    snaps.sort()
    return snaps

def load_shares():
    sh = pd.read_csv(os.path.join(DATA, 'shares_pit.csv'), parse_dates=['end', 'filed'])
    sh = sh.sort_values(['ticker', 'filed'])
    return sh

def month_ends(px):
    alld = set()
    for df in px.values():
        alld.update(df['date'].dt.strftime('%Y-%m-%d'))
    alld = sorted(alld)
    me = []
    for d in alld:
        if '1999-12-01' <= d <= '2022-12-31':
            ym = d[:7]
            # last trading day of that calendar month
            cands = [x for x in alld if x[:7] == ym]
            if d == max(cands):
                me.append(d)
    return me

def main():
    px = load_prices()
    print('price series:', len(px), flush=True)
    snaps = load_membership()
    sh = load_shares()
    print('shares rows:', len(sh), flush=True)
    me = month_ends(px)
    print('month-ends:', len(me), me[0], me[-1], flush=True)

    # membership per month-end: latest snapshot <= month-end
    mem_by_me = {}
    si = 0
    for m in me:
        while si + 1 < len(snaps) and snaps[si + 1][0] <= m:
            si += 1
        mem_by_me[m] = snaps[si][1] if snaps[si][0] <= m else set()

    # shares lookup: per ticker sorted by filed
    sh_by_t = {t: g for t, g in sh.groupby('ticker')}

    out_rows = []
    for mi, m in enumerate(me):
        mend = pd.Timestamp(m)
        members = mem_by_me[m]
        for t in members:
            pt = ALIAS.get(t, t)
            df = px.get(pt)
            if df is None:
                out_rows.append((m, t, pt, 0, 1, 0, *[np.nan]*9, 'no_price'))
                continue
            d = df[df['date'] <= mend]
            if len(d) < 30:
                out_rows.append((m, t, pt, 0, 1, 0, *[np.nan]*9, 'short_history'))
                continue
            last = d.iloc[-1]
            # trading-day offsets
            def val_at(k, col):
                return d[col].iloc[-1-k] if len(d) > k else np.nan
            adj_now = last['adjclose']
            adj_21 = val_at(21, 'adjclose'); adj_252 = val_at(252, 'adjclose')
            mom = adj_21/adj_252 - 1 if adj_252 and adj_252 > 0 else np.nan
            rets = np.log(d['adjclose'].iloc[-60:].values[1:] / d['adjclose'].iloc[-60:].values[:-1])
            rets = rets[np.isfinite(rets)]
            vol = np.std(rets, ddof=1)*np.sqrt(252) if len(rets) > 20 else np.nan
            # point-in-time shares
            g = sh_by_t.get(pt)
            shares = np.nan; filed = ''
            if g is not None:
                gg = g[g['filed'] <= mend]
                if len(gg):
                    shares = gg.iloc[-1]['shares']; filed = str(gg.iloc[-1]['filed'].date())
            vols = d['volume'].iloc[-21:].values
            turnover = np.nanmean(vols/shares) if shares and shares > 0 else np.nan
            mcap = last['close']*shares if shares and shares > 0 else np.nan
            prov = []
            prov.append('px_measured')
            prov.append('shares_pit_measured' if np.isfinite(shares) else 'shares_missing')
            out_rows.append((m, t, pt, 1, 0, 0, adj_now, last['close'], last['volume'],
                             mom, vol, turnover, shares, mcap, filed, '|'.join(prov)))
        if (mi+1) % 24 == 0:
            print(f'{mi+1}/{len(me)} {m}', flush=True)

    cols = ['month', 'ticker', 'px_ticker', 'has_price', 'no_price', 'short_hist',
            'adjclose', 'close_raw', 'volume', 'mom', 'vol', 'turnover',
            'shares_pit', 'mcap', 'shares_filed', 'prov']
    panel = pd.DataFrame(out_rows, columns=cols)
    panel.to_csv(os.path.join(DATA, 'panel.csv'), index=False)
    print('panel rows:', len(panel), flush=True)
    # coverage summary
    cov = panel.groupby('month').agg(n_mem=('ticker', 'size'),
                                     n_px=('has_price', 'sum'),
                                     n_mcap=('mcap', lambda s: s.notna().sum()))
    cov.to_csv(os.path.join(DATA, 'coverage.csv'), index=True)
    print(cov.iloc[::36].to_string(), flush=True)

if __name__ == '__main__':
    main()
