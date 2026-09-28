"""Build annual point-in-time fundamentals panel from per-ticker companyfacts JSONs.

For each ticker: annual records (end, filed, NI, Equity, Assets, Revenues) from
10-K facts. Writes data/fundamentals_annual.csv with one row per ticker/fy-end.
Labels: measured (SEC XBRL companyfacts, 10-K).
"""
import json, os, glob
import pandas as pd
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
FDIR = os.path.join(DATA, 'fundamentals')

def pick(records, asof):
    """Latest end date with filed <= asof and end within 18 months before asof."""
    asof_d = date.fromisoformat(asof)
    best = None
    for end, filed, v in records:
        try:
            e = date.fromisoformat(end[:10])
        except Exception:
            continue
        f = filed[:10] if filed else '0000-00-00'
        if f > asof:
            continue
        if (asof_d - e).days > 548 or (asof_d - e).days < 0:
            continue
        if best is None or e > best[0]:
            best = (e, f, v, end[:10])
    return best  # (end_date, filed, value, end_str)

def main():
    rows = []
    files = glob.glob(os.path.join(FDIR, '*.json'))
    for fp in files:
        ticker = os.path.basename(fp)[:-5]
        try:
            d = json.load(open(fp))
        except Exception:
            continue
        # collect all distinct end dates across tags
        ends = set()
        for key in ('NetIncomeLoss', 'Equity', 'Assets', 'Revenues'):
            for end, filed, v in d.get(key, []):
                ends.add(end[:10])
        for end in sorted(ends):
            rec = {'ticker': ticker, 'fy_end': end}
            ok = True
            for key in ('NetIncomeLoss', 'Equity', 'Assets', 'Revenues'):
                match = [(e, f, v) for (e, f, v) in d.get(key, []) if e[:10] == end]
                if not match:
                    rec[key] = None
                else:
                    # earliest filed for this end
                    match.sort(key=lambda x: x[1])
                    rec[key] = match[0][2]
                    rec[key + '_filed'] = match[0][1][:10]
            rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, 'fundamentals_annual.csv'), index=False)
    print('rows:', len(df), 'tickers:', df['ticker'].nunique())
    # coverage by fiscal year
    df['yr'] = df['fy_end'].str[:4]
    cov = df.groupby('yr')['ticker'].nunique()
    print(cov.to_string())
    # completeness of the 4 tags
    for k in ('NetIncomeLoss', 'Equity', 'Assets', 'Revenues'):
        print(k, 'non-null:', df[k].notna().sum())

if __name__ == '__main__':
    main()
