"""Download EDGAR XBRL shares-outstanding history per CIK, point-in-time safe.
Writes data/shares_out.csv: cik,ticker,end,filed,form,shares (as-reported, unadjusted).
Labels: measured (SEC filings) for end>=2009 filing-based; pre-2009 estimated.
Key: retains FILED date; panel must only use rows with filed <= asof date.
Tag fallback: us-gaap:CommonStockSharesOutstanding -> dei:EntityCommonStockSharesOutstanding.
"""
import csv, json, time, urllib.request, os
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
        facts = d.get('facts', {})
        pts = []
        g = facts.get('us-gaap', {}).get('CommonStockSharesOutstanding')
        if g: pts += g['units'].get('shares', [])
        dei = facts.get('dei', {}).get('EntityCommonStockSharesOutstanding')
        if dei: pts += dei['units'].get('shares', [])
        if not pts:
            return (ticker, 'notag', [])
        # keep 10-K/10-Q, earliest filing per (end,form) to avoid restatement lookahead
        best = {}
        for p in pts:
            if p.get('form') not in ('10-K', '10-Q'):
                continue
            key = (p['end'], p['form'])
            if key not in best or p['filed'] < best[key][0]:
                best[key] = (p['filed'], p['val'])
        rows = sorted((e, f, fm, v) for (e, fm), (f, v) in best.items())
        return (ticker, 'ok', rows)
    except Exception as e:
        time.sleep(0.3)
        return (ticker, f'err {str(e)[:60]}', [])

def main():
    pairs, seen = [], set()
    with open(os.path.join(DATA, 'constituents.csv')) as f:
        for row in csv.DictReader(f):
            cik, t = row['CIK'].strip(), row['Symbol'].replace('.', '-').strip()
            if cik and t not in seen:
                seen.add(t); pairs.append((cik, t))
    # historical CIKs
    try:
        r = urllib.request.Request('https://www.sec.gov/files/company_tickers.json', headers=UA)
        with urllib.request.urlopen(r, timeout=60) as h:
            sec = {v['ticker'].upper(): str(v['cik_str']).zfill(10)
                   for v in json.loads(h.read()).values()}
    except Exception as e:
        print('sec map fail', e); sec = {}
    import glob
    need = set()
    for f in glob.glob(os.path.join(DATA, 'raw', '*.csv')):
        t = os.path.basename(f)[:-4]
        if t not in seen and t in sec:
            need.add(t)
    for t in sorted(need):
        pairs.append((sec[t], t))
    print('total pairs:', len(pairs), flush=True)
    out = os.path.join(DATA, 'shares_pit.csv')
    have = set()  # fresh point-in-time file; always refetch
    todo = pairs
    print('todo:', len(todo), flush=True)
    ok = err = 0
    with open(out, 'a', newline='') as f:
        w = csv.writer(f)
        if not have:
            w.writerow(['cik', 'ticker', 'end', 'filed', 'form', 'shares'])
        with ThreadPoolExecutor(max_workers=6) as ex:
            for i, (ticker, status, rows) in enumerate(ex.map(fetch_cik, todo)):
                if status == 'ok':
                    ok += 1
                    for (e, filed, form, v) in rows:
                        w.writerow([todo[i][0], ticker, e, filed, form, int(v)])
                else:
                    err += 1
                    if err < 15: print(ticker, status, flush=True)
                if (i+1) % 50 == 0:
                    print(f'{i+1}/{len(todo)} ok={ok} err={err}', flush=True)
    print(f'DONE ok={ok} err={err}', flush=True)

if __name__ == '__main__':
    main()
