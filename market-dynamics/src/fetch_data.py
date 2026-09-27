"""Parallel fetch of daily OHLCV for S&P 500 tickers from Yahoo Finance API."""
import csv, json, time, urllib.request, os
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RAW = os.path.join(DATA, 'raw')
os.makedirs(RAW, exist_ok=True)

P1, P2 = 1420070400, 1640995200  # 2015-01-01 to 2022-01-01

def fetch_one(t):
    url = (f"https://query2.finance.yahoo.com/v8/finance/chart/{t}"
           f"?interval=1d&period1={P1}&period2={P2}&events=div%2Csplit")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
        res = d['chart']['result']
        if not res:
            return (t, 'noresult')
        r0 = res[0]
        ts = r0.get('timestamp') or []
        if len(ts) < 300:
            return (t, 'short')
        q = r0['indicators']['quote'][0]
        adj = r0['indicators'].get('adjclose', [{}])[0].get('adjclose')
        with open(os.path.join(RAW, f"{t}.csv"), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['date', 'open', 'high', 'low', 'close', 'adjclose', 'volume'])
            for k in range(len(ts)):
                w.writerow([
                    time.strftime('%Y-%m-%d', time.gmtime(ts[k])),
                    q['open'][k], q['high'][k], q['low'][k], q['close'][k],
                    adj[k] if adj else q['close'][k], q['volume'][k],
                ])
        return (t, 'ok')
    except Exception as e:
        return (t, f'err {e}')

def main():
    with open(os.path.join(DATA, 'constituents.csv')) as f:
        tickers = [row['Symbol'].replace('.', '-') for row in csv.DictReader(f)]
    have = {fn[:-4] for fn in os.listdir(RAW) if fn.endswith('.csv')}
    todo = [t for t in tickers if t not in have]
    print(f"{len(tickers)} tickers, {len(todo)} to fetch", flush=True)
    n_ok = 0
    with ThreadPoolExecutor(max_workers=8) as ex:
        for i, (t, st) in enumerate(ex.map(fetch_one, todo)):
            if st == 'ok':
                n_ok += 1
            else:
                print(f"{t}: {st}", flush=True)
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(todo)} ok={n_ok}", flush=True)
    print(f"DONE ok={n_ok}", flush=True)

if __name__ == '__main__':
    main()
