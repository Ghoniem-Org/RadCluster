"""Fetch/extend Yahoo daily data: extend existing to 2022-12-31, fetch missing historical members."""
import csv, json, time, urllib.request, os, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RAW = os.path.join(DATA, 'raw')
os.makedirs(RAW, exist_ok=True)

P1_NEW, P2_NEW = 946684800, 1672531200  # 2000-01-01 to 2023-01-01
UA = {"User-Agent": "Mozilla/5.0"}

def yahoo_bars(t, p1, p2):
    url = (f"https://query2.finance.yahoo.com/v8/finance/chart/{t}"
           f"?interval=1d&period1={p1}&period2={p2}&events=div%2Csplit")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=30) as r:
        d = json.load(r)
    res = d['chart']['result']
    if not res:
        return None
    r0 = res[0]
    ts = r0.get('timestamp') or []
    q = r0['indicators']['quote'][0]
    adj = r0['indicators'].get('adjclose', [{}])[0].get('adjclose')
    rows = []
    for k in range(len(ts)):
        if q['close'][k] is None:
            continue
        rows.append([
            time.strftime('%Y-%m-%d', time.gmtime(ts[k])),
            q['open'][k], q['high'][k], q['low'][k], q['close'][k],
            adj[k] if adj and adj[k] else q['close'][k], q['volume'][k],
        ])
    return rows

def fetch_full(t):
    """Full history fetch for a missing ticker; keep if >=200 bars."""
    try:
        rows = yahoo_bars(t, P1_NEW, P2_NEW)
        if not rows or len(rows) < 200:
            return (t, 'short/none')
        with open(os.path.join(RAW, f"{t}.csv"), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['date', 'open', 'high', 'low', 'close', 'adjclose', 'volume'])
            w.writerows(rows)
        return (t, f'ok {len(rows)}')
    except Exception as e:
        return (t, f'err {str(e)[:60]}')

def extend_one(t):
    """Append 2022 data to an existing CSV."""
    fn = os.path.join(RAW, f"{t}.csv")
    try:
        with open(fn) as f:
            lines = f.readlines()
        last = lines[-1].split(',')[0]
        if last >= '2022-12-01':
            return (t, 'already-ok')
        p1 = int(time.mktime(time.strptime('2021-12-15', '%Y-%m-%d')))
        rows = yahoo_bars(t, p1, P2_NEW)
        if not rows:
            return (t, 'noresult')
        have = {ln.split(',')[0] for ln in lines[1:]}
        new = [r for r in rows if r[0] not in have and r[0] > last]
        with open(fn, 'a', newline='') as f:
            w = csv.writer(f)
            w.writerows(new)
        return (t, f'extended +{len(new)}')
    except Exception as e:
        return (t, f'err {str(e)[:60]}')

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'extend'
    if mode == 'extend':
        todo = [fn[:-4] for fn in os.listdir(RAW) if fn.endswith('.csv')]
        fn_ = extend_one
    else:  # fetch-missing: tickers from file
        tickers = [l.strip() for l in open(sys.argv[2]) if l.strip()]
        have = {fn[:-4] for fn in os.listdir(RAW) if fn.endswith('.csv')}
        todo = [t for t in tickers if t not in have]
        fn_ = fetch_full
    print(f"{len(todo)} tickers, mode={mode}", flush=True)
    n_ok = 0
    with ThreadPoolExecutor(max_workers=8) as ex:
        for i, (t, st) in enumerate(ex.map(fn_, todo)):
            if st.startswith('ok') or st.startswith('extended') or st == 'already-ok':
                n_ok += 1
            else:
                print(f"{t}: {st}", flush=True)
            if (i + 1) % 60 == 0:
                print(f"  {i+1}/{len(todo)} ok={n_ok}", flush=True)
    print(f"DONE ok={n_ok}", flush=True)

if __name__ == '__main__':
    main()
