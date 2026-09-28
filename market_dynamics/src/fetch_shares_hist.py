"""Fetch EDGAR shares for historical (non-current) tickers. Appends to data/shares_out.csv."""
import csv, json, time, urllib.request, os, gzip
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
UA = {"User-Agent": "MuseResearch/1.0 contact@example.com", "Accept-Encoding": "gzip"}

def fetch_cik(args):
    cik, ticker = args
    cikp = str(cik).zfill(10)
    url = f"https://data.sec.gov/api/xbrl/companyfacts/CIK{cikp}.json"
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=40) as r:
            raw = r.read()
            if r.headers.get('Content-Encoding') == 'gzip':
                raw = gzip.decompress(raw)
            d = json.loads(raw)
        facts = d.get('facts', {}).get('us-gaap', {})
        tag = facts.get('CommonStockSharesOutstanding')
        if not tag:
            return (ticker, 'notag', [])
        pts = tag['units'].get('shares', [])
        best = {}
        for p in pts:
            if p.get('form') not in ('10-K', '10-Q'):
                continue
            e = p['end']
            if e not in best or p['filed'] > best[e][0]:
                best[e] = (p['filed'], p['val'])
        rows = sorted((e, v[1]) for e, v in best.items())
        return (ticker, 'ok', rows)
    except Exception as e:
        time.sleep(0.3)
        return (ticker, f'err {str(e)[:50]}', [])

def main():
    pairs = []
    seen_cik = set()
    # skip CIKs already fetched
    with open(os.path.join(DATA, 'shares_out.csv')) as f:
        for row in csv.DictReader(f):
            seen_cik.add(row['cik'].lstrip('0') or '0')
    with open('/tmp/hist_cik_pairs.csv') as f:
        for line in f:
            cik, t = line.strip().split(',')
            if cik.lstrip('0') not in seen_cik:
                pairs.append((cik, t))
    print(f"{len(pairs)} historical CIKs to fetch", flush=True)
    out_rows = []
    n_ok = 0
    with ThreadPoolExecutor(max_workers=5) as ex:
        for i, (ticker, st, rows) in enumerate(ex.map(fetch_cik, pairs)):
            if st == 'ok':
                n_ok += 1
                for e, v in rows:
                    out_rows.append((pairs[i][0], ticker, e, v))
            if (i + 1) % 40 == 0:
                print(f"  {i+1}/{len(pairs)} ok={n_ok}", flush=True)
            time.sleep(0.08)
    out_rows.sort(key=lambda r: (r[1], r[2]))
    with open(os.path.join(DATA, 'shares_out.csv'), 'a', newline='') as f:
        w = csv.writer(f)
        w.writerows(out_rows)
    print(f"DONE hist_ok={n_ok} new_points={len(out_rows)}", flush=True)

if __name__ == '__main__':
    main()
