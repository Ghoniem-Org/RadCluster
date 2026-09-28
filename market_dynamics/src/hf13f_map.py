"""Map 13F CUSIPs -> tickers via Yahoo search API (OpenFIGI unreachable).
Incremental cache: data/cusip_yahoo_map.csv
"""
import pandas as pd, requests, os, time, threading
from concurrent.futures import ThreadPoolExecutor

BASE = os.path.expanduser("~/workspace/market_dynamics")
CACHE = os.path.join(BASE, "data", "cusip_yahoo_map.csv")
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64)"}

def lookup(cusip):
    for a in range(3):
        try:
            r = requests.get("https://query2.finance.yahoo.com/v1/finance/search",
                             params={"q": cusip, "quotesCount": 6}, headers=UA, timeout=20)
            if r.status_code == 429:
                time.sleep(10); continue
            r.raise_for_status()
            for q in r.json().get("quotes", []):
                if q.get("quoteType") == "EQUITY" and q.get("symbol"):
                    return q["symbol"], q.get("longname") or q.get("shortname") or ""
            return "", ""
        except Exception:
            time.sleep(2 * (a + 1))
    return "", ""

def main():
    h = pd.read_csv(os.path.join(BASE, "data", "hf_holdings_raw.csv"), dtype={"cusip": str})
    cusips = sorted(h["cusip"].unique())
    done = {}
    if os.path.exists(CACHE):
        d = pd.read_csv(CACHE, dtype=str).fillna("")
        done = dict(zip(d["cusip"], zip(d["ticker"], d["name"])))
    todo = [c for c in cusips if c not in done]
    print(f"{len(cusips)} cusips, {len(todo)} to map")
    lock = threading.Lock()
    n = [0]
    def work(c):
        t, nm = lookup(c)
        with lock:
            done[c] = (t, nm)
            n[0] += 1
            if n[0] % 200 == 0:
                pd.DataFrame([(k, v[0], v[1]) for k, v in done.items()],
                             columns=["cusip", "ticker", "name"]).to_csv(CACHE, index=False)
                print(f"  {n[0]}/{len(todo)} mapped", flush=True)
    with ThreadPoolExecutor(max_workers=6) as ex:
        list(ex.map(work, todo))
    pd.DataFrame([(k, v[0], v[1]) for k, v in done.items()],
                 columns=["cusip", "ticker", "name"]).to_csv(CACHE, index=False)
    hit = sum(1 for v in done.values() if v[0])
    print(f"done: {hit}/{len(done)} mapped to a ticker")

if __name__ == "__main__":
    main()
