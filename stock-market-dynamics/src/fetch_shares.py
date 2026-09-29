"""Download EDGAR XBRL shares-outstanding history per CIK.
Writes data/shares_out.csv: cik,ticker,date,shares (as-reported, unadjusted).
Label: measured (SEC filings, point-in-time) from first XBRL filing (~2009);
pre-2009 values are backfilled from earliest filing -> estimated.
"""
import csv, json, time, urllib.request, os, sys
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
                import gzip as _gz
                raw = _gz.decompress(raw)
            d = json.loads(raw)
        facts = d.get('facts', {}).get('us-gaap', {})
        tag = facts.get('CommonStockSharesOutstanding')
        if not tag:
            return (ticker, 'notag', 0)
        pts = tag['units'].get('shares', [])
        # keep 10-K/10-Q, latest per end date
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
        time.sleep(0.2)
        return (ticker, f'err {str(e)[:50]}', 0)

def main():
    pairs = []
    seen = set()
    # current constituents with CIK
    with open(os.path.join(DATA, 'constituents.csv')) as f:
        for row in csv.DictReader(f):
            cik, t = row['CIK'].strip(), row['Symbol'].replace('.', '-').strip()
            if cik and cik not in seen:
                seen.add(cik)
                pairs.append((cik, t))
    # historical tickers via SEC ticker map
    try:
        req = urllib.request.Request("https://www.sec.gov/files/company_tickers.json", headers=UA)
        with urllib.request.urlopen(req, timeout=40) as r:
            tmap = json.load(r)
        for k, v in tmap.items():
            cik, t = str(v['cik_str']), v['ticker'].replace('.', '-')
            if cik not in seen:
                seen.add(cik)
                pairs.append((cik, t))
        print('ticker map entries:', len(tmap), flush=True)
    except Exception as e:
        print('ticker map failed:', str(e)[:80], flush=True)
    print(f"{len(pairs)} CIKs to fetch", flush=True)
    out_rows = []
    n_ok = n_pts = 0
    with ThreadPoolExecutor(max_workers=6) as ex:
        for i, (ticker, st, rows) in enumerate(ex.map(fetch_cik, pairs)):
            if st == 'ok':
                n_ok += 1
                n_pts += len(rows)
                for e, v in rows:
                    out_rows.append((pairs[i][0], ticker, e, v))
            elif not st.startswith('err'):
                print(f"{ticker}: {st}", flush=True)
            if (i + 1) % 60 == 0:
                print(f"  {i+1}/{len(pairs)} ok={n_ok} pts={n_pts}", flush=True)
            time.sleep(0.05)
    out_rows.sort(key=lambda r: (r[1], r[2]))
    with open(os.path.join(DATA, 'shares_out.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['cik', 'ticker', 'date', 'shares'])
        w.writerows(out_rows)
    print(f"DONE cik_ok={n_ok} points={n_pts}", flush=True)

if __name__ == '__main__':
    main()
