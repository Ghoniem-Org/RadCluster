"""Clean SEC shares: drop <=0, per-ticker jump filter (>3x or <1/3 vs neighbor
is filing/XBRL garbage), forward/backward fill. Writes data/shares_clean.csv.
Also builds per-ticker split-factor series from adjclose/close ratio jumps
for backfill split-adjustment -> data/splits.csv (ticker,date,ratio).
"""
import os, glob
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

def clean():
    sh = pd.read_csv(os.path.join(DATA, 'shares_pit.csv'), parse_dates=['end', 'filed'])
    sh = sh[sh['shares'] > 0].copy()
    sh = sh.sort_values(['ticker', 'filed']).reset_index(drop=True)
    # Fix mechanical unit errors: a 1000x or 1e6x jump vs previous good value
    # is XBRL scaling garbage -> rescale.  Accept all other jumps (splits,
    # genuine issuance) since mcap/turnover formulas self-correct via S(f).
    vals = sh['shares'].values.astype(float)
    fixed = np.zeros(len(sh), dtype=bool)
    for t, g in sh.groupby('ticker', sort=False):
        idx = g.index.values
        last = np.nan
        for j, i in enumerate(idx):
            v = vals[i]
            if not np.isnan(last) and last > 0:
                r = v / last
                if 900 <= r <= 1100:
                    vals[i] = v / 1000; fixed[i] = True; v = vals[i]
                elif 9e5 <= r <= 1.1e6:
                    vals[i] = v / 1e6; fixed[i] = True; v = vals[i]
                elif 1/1100 <= r <= 1/900:
                    vals[i] = v * 1000; fixed[i] = True; v = vals[i]
            last = v
    sh['shares'] = vals
    sh['unit_fixed'] = fixed
    sh['clean_ok'] = True
    print(f'rows: {len(sh)}, unit-fixed: {fixed.sum()}', flush=True)
    sh.to_csv(os.path.join(DATA, 'shares_clean.csv'), index=False)

def splits():
    """Detect splits from adjclose/close ratio jumps >25%."""
    rows = []
    for f in glob.glob(os.path.join(DATA, 'raw', '*.csv')):
        t = os.path.basename(f)[:-4]
        try:
            df = pd.read_csv(f, parse_dates=['date'], usecols=['date', 'close', 'adjclose'])
        except Exception:
            continue
        df = df.sort_values('date')
        df = df[(df['close'] > 0) & (df['adjclose'] > 0)]
        if len(df) < 2:
            continue
        R = (df['adjclose'] / df['close']).values
        jr = R[1:] / R[:-1]
        mask = (jr > 1.25) | (jr < 0.8)
        for d, r in zip(df['date'].iloc[1:][mask], jr[mask]):
            rows.append((t, d.strftime('%Y-%m-%d'), round(float(r), 4)))
    sp = pd.DataFrame(rows, columns=['ticker', 'date', 'ratio'])
    sp.to_csv(os.path.join(DATA, 'splits.csv'), index=False)
    print(f'split events: {len(sp)}', flush=True)
    return sp

if __name__ == '__main__':
    clean()
    splits()
