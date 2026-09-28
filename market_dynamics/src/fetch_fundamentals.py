"""Fetch annual fundamentals from SEC XBRL companyfacts per CIK (point-in-time safe).

Tags: NetIncomeLoss, StockholdersEquity (+fallback incl. NCI), Assets,
      Revenues (+fallbacks SalesRevenueNet, RevenueFromContractWithCustomerExcludingAssessedTax).
Keeps 10-K facts only; duration tags require ~1yr span; earliest filed per
(end,form) to avoid restatement lookahead. Applies absolute-anchor bounds
(the $5,000T lesson): values outside plausible ranges are dropped and counted.
Writes data/fundamentals/<TICKER>.json (compact) + data/fundamentals_summary.csv.
Labels: measured (SEC filings).
"""
import csv, json, os, pickle, time, urllib.request, gzip
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
FDIR = os.path.join(DATA, 'fundamentals')
os.makedirs(FDIR, exist_ok=True)
UA = {"User-Agent": "MuseResearch/1.0 contact@example.com", "Accept-Encoding": "gzip"}

TAGS = {
    'NetIncomeLoss': ('us-gaap', ['NetIncomeLoss'], 'duration'),
    'Equity': ('us-gaap', ['StockholdersEquity',
                            'StockholdersEquityIncludingPortionAttributableToNoncontrollingInterest'], 'instant'),
    'Assets': ('us-gaap', ['Assets'], 'instant'),
    'Revenues': ('us-gaap', ['Revenues', 'SalesRevenueNet',
                              'RevenueFromContractWithCustomerExcludingAssessedTax'], 'duration'),
}
# absolute anchors (USD); outside -> drop. Documented, not tuned.
BOUNDS = {
    'NetIncomeLoss': (-1e12, 1e12),
    'Equity': (-1e12, 1e13),
    'Assets': (0.0, 2e13),
    'Revenues': (0.0, 1e13),
}

def get_json(url, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        raw = r.read()
        if r.headers.get('Content-Encoding') == 'gzip':
            raw = gzip.decompress(raw)
    return json.loads(raw)

def extract(cik, ticker):
    cikp = str(cik).zfill(10)
    d = get_json(f"https://data.sec.gov/api/xbrl/companyfacts/CIK{cikp}.json")
    facts = d.get('facts', {})
    out, dropped = {}, 0
    for key, (ns, names, kind) in TAGS.items():
        pts = []
        for nm in names:
            g = facts.get(ns, {}).get(nm)
            if not g:
                continue
            for unit, arr in g.get('units', {}).items():
                if unit != 'USD':
                    continue
                for p in arr:
                    if p.get('form') != '10-K':
                        continue
                    try:
                        v = float(p['val'])
                    except (TypeError, ValueError):
                        continue
                    lo, hi = BOUNDS[key]
                    if not (lo <= v <= hi):
                        dropped += 1
                        continue
                    if kind == 'duration':
                        try:
                            from datetime import date
                            s = date.fromisoformat(p['start'][:10]); e = date.fromisoformat(p['end'][:10])
                            if not (340 <= (e - s).days <= 400):
                                continue
                        except Exception:
                            continue
                    pts.append(p)
            if pts:
                break  # first tag name with data wins
        # earliest filed per end date
        best = {}
        for p in pts:
            e = p['end'][:10]
            if e not in best or p.get('filed', '') < best[e][0]:
                best[e] = (p.get('filed', ''), float(p['val']))
        out[key] = sorted((e, f, v) for e, (f, v) in best.items())
    return ticker, out, dropped

def main():
    cikmap = pickle.load(open('/tmp/cikmap.pkl', 'rb'))
    todo = sorted(cikmap.items())
    print('todo:', len(todo), flush=True)
    ok = err = drop_tot = 0
    with ThreadPoolExecutor(max_workers=8) as ex:
        for i, res in enumerate(ex.map(lambda kv: safe(kv), todo)):
            ticker, payload, status, dropped = res
            drop_tot += dropped
            if status == 'ok':
                ok += 1
                with open(os.path.join(FDIR, f'{ticker}.json'), 'w') as f:
                    json.dump(payload, f)
            else:
                err += 1
                if err <= 10:
                    print(ticker, status, flush=True)
            if (i + 1) % 50 == 0:
                print(f'{i+1}/{len(todo)} ok={ok} err={err} dropped_vals={drop_tot}', flush=True)
    print(f'DONE ok={ok} err={err} dropped_vals={drop_tot}', flush=True)

def safe(kv):
    ticker, cik = kv
    try:
        t, payload, dropped = extract(cik, ticker)
        return t, payload, 'ok', dropped
    except Exception as e:
        time.sleep(0.2)
        return ticker, {}, f'err {str(e)[:70]}', 0

if __name__ == '__main__':
    main()
