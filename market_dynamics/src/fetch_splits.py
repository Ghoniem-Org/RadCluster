"""Fetch split events from Yahoo v8 for all raw tickers -> data/splits_yahoo.csv."""
import csv, glob, json, os, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"}

def get_splits(t):
    url = f"https://query1.finance.yahoo.com/v8/finance/chart/{t}?interval=1d&period1=0&period2=9999999999&events=split"
    try:
        req = urllib.request.Request(url, headers=UA)
        d = json.loads(urllib.request.urlopen(req, timeout=30).read())
        ev = d['chart']['result'][0].get('events', {}).get('splits', {})
        out = []
        for k, v in ev.items():
            ts = int(float(k))
            out.append((t, time.strftime('%Y-%m-%d', time.gmtime(ts)),
                        v['numerator']/v['denominator']))
        return out
    except Exception as e:
        return None

def main():
    tickers = [os.path.basename(f)[:-4] for f in glob.glob(os.path.join(DATA, 'raw', '*.csv'))]
    rows = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for i, res in enumerate(ex.map(get_splits, tickers)):
            if res:
                rows.extend(res)
            if (i+1) % 100 == 0:
                print(f'{i+1}/{len(tickers)} events={len(rows)}', flush=True)
    with open(os.path.join(DATA, 'splits_yahoo.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['ticker', 'date', 'ratio'])
        w.writerows(sorted(rows))
    print(f'DONE tickers={len(tickers)} split_events={len(rows)}', flush=True)

if __name__ == '__main__':
    main()
