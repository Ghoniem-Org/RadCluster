"""Build 13F-based hedge-fund crowding panel for S&P 500 members, 1999-2021.

Closest PUBLIC-data analogue of Brown-Howard-Lundblad (2022) crowdedness:
their 4 measures use a proprietary HF holdings DB (cannot replicate exactly).
We use 13F-HR filings of a documented 18-manager hedge-fund subset (assumed
choice, see MANAGERS), parsed from SEC EDGAR. 13F covers long US equity
positions only (>$100M managers), filed <=45d after quarter-end; shorts and
non-US holdings are invisible. Amendments (13F-HR/A) excluded; per
(manager, report-quarter) the latest filingDate is kept.

Outputs:
  data/hf_crowding_quarterly.csv  (ticker, yq, n_hf, hf_shares, hf_hhi,
                                   hf_own, hf_adv_ratio, n_mgr_filed)
  data/hf_cusip_map.csv           (cusip -> ticker via OpenFIGI)
  data/raw13f_parsed/             per-filing parsed positions (cache)
"""
import requests, re, json, time, os, io, sys
import pandas as pd
from lxml import etree

H = {"User-Agent": "Muse research contact@example.com", "Accept-Encoding": "gzip"}
BASE = os.path.expanduser("~/workspace/market_dynamics")
RAW = os.path.join(BASE, "data", "raw13f_parsed")
os.makedirs(RAW, exist_ok=True)

# (name, CIK) -- assumed HF subset; CIKs verified via EDGAR 2026-09-28.
# SAC covers the 1999-2013 era (became Point72); others cover 2000s-2021.
MANAGERS = [
    ("Bridgewater Associates LP", "0001350694"),
    ("Millennium Management LLC", "0001273087"),
    ("AQR Capital Management LLC", "0001167557"),
    ("Elliott Management Corp", "0001048445"),
    ("Lone Pine Capital LLC", "0001061165"),
    ("Tiger Global Management LLC", "0001167483"),
    ("Coatue Management LLC", "0001135730"),
    ("Renaissance Technologies LLC", "0001037389"),
    ("Third Point LLC", "0001040273"),
    ("Pershing Square", "0001336528"),
    ("Greenlight Capital Inc", "0001079114"),
    ("Soros Fund Management LLC", "0001029160"),
    ("Farallon Capital Management LLC", "0000909661"),
    ("Davidson Kempner", "0000937617"),
    ("Sculptor Capital LP", "0001054587"),
    ("King Street Capital Management LP", "0001218199"),
    ("Two Sigma Investments LP", "0001179392"),
    ("SAC Capital Advisors LLC", "0001018103"),
]

ITNS = "http://www.sec.gov/edgar/document/thirteenf/informationtable"
NS = {"it": ITNS}

def get(url, timeout=40, retries=5):
    last = None
    for a in range(retries):
        try:
            r = requests.get(url, headers=H, timeout=timeout)
            r.raise_for_status()
            return r
        except Exception as e:
            last = e
            time.sleep(4 * (a + 1))
    raise last

def filing_list(cik):
    """All 13F-HR (not /A) accessions 1999-2022 for a CIK."""
    out = []
    d = get(f"https://data.sec.gov/submissions/CIK{cik}.json").json()
    buckets = [d["filings"]["recent"]]
    for fl in d["filings"].get("files", []):
        buckets.append(get(f"https://data.sec.gov/submissions/{fl['name']}").json())
    for b in buckets:
        for form, dt, acc in zip(b["form"], b["filingDate"], b["accessionNumber"]):
            if form == "13F-HR" and "1999" <= dt <= "2022":
                out.append((dt, acc))
    return out

def parse_infotable(xml_bytes):
    root = etree.fromstring(xml_bytes)
    rows = []
    for t in root.findall("it:infoTable", NS):
        try:
            cusip = t.find("it:cusip", NS).text.strip()
            amt = t.find("it:shrsOrPrnAmt", NS)
            typ = amt.find("it:sshPrnamtType", NS).text.strip()
            if typ != "SH":
                continue
            sh = int(amt.find("it:sshPrnamt", NS).text.strip())
            if len(cusip) == 9 and sh > 0:
                rows.append((cusip, sh))
        except Exception:
            continue
    return rows

def parse_old_txt(txt):
    """Best-effort pre-XML 13F info-table parse via CUSIP-anchored regex.
    Old paper-format tables print one row per holding with a 9-char CUSIP,
    value ($000), share count, and SH/PRN flag. We anchor on the CUSIP."""
    rows = []
    pat = re.compile(r"\b([0-9A-Z]{9})\b\s+([\d,]+)\s+([\d,]+)\s+(SH|PRN)\b")
    for ln in txt.splitlines():
        m = pat.search(ln)
        if m:
            cusip, _val, shs, typ = m.groups()
            if typ == "SH":
                rows.append((cusip, int(shs.replace(",", ""))))
    # dedupe (table may repeat in summary)
    seen, out = set(), []
    for c, s in rows:
        if (c, s) not in seen:
            seen.add((c, s))
            out.append((c, s))
    return out

def parse_filing(cik, dt, acc):
    """Returns (period_yq, [(cusip, shares)]). period_yq like '2001Q1'."""
    accn = acc.replace("-", "")
    idx = get(f"https://www.sec.gov/Archives/edgar/data/{int(cik)}/{accn}/index.json").json()
    docs = {it["name"]: it for it in idx["directory"]["item"]}
    # period of report: prefer primary_doc.xml, else SGML header of .txt
    yq = None
    names = list(docs.keys())
    pdoc = next((n for n in names if n == "primary_doc.xml"), None)
    if pdoc:
        xml = get(f"https://www.sec.gov/Archives/edgar/data/{int(cik)}/{accn}/{pdoc}").content
        m = re.search(rb"<periodOfReport>([^<]+)", xml)
        if m:
            p = m.group(1).decode()
            mm, dd, yyyy = p.split("-")
            yq = f"{yyyy}Q{(int(mm)-1)//3+1}"
    info = next((n for n in names if "infotable" in n.lower() and n.lower().endswith(".xml")), None)
    cand = []
    if info:
        cand = [info]
    else:
        # some filers use bespoke names (e.g. MMLLC.xml): scan all xml docs
        cand = [n for n in names if n.lower().endswith(".xml") and n != "primary_doc.xml"
                and "index" not in n.lower()]
    rows = []
    for doc in cand:
        try:
            xml = get(f"https://www.sec.gov/Archives/edgar/data/{int(cik)}/{accn}/{doc}").content
            rws = parse_infotable(xml)
            if rws:
                rows = rws
                break
        except Exception:
            continue
    if rows:
        if yq is None:
            # no primary_doc.xml: get period from SGML header of the .txt
            txtn = next((n for n in names if n.endswith(".txt") and "index" not in n), None)
            if txtn:
                try:
                    txt = get(f"https://www.sec.gov/Archives/edgar/data/{int(cik)}/{accn}/{txtn}").text
                    m = re.search(r"CONFORMED PERIOD OF REPORT:\s*(\d{8})", txt)
                    if m:
                        p = m.group(1)
                        yq = f"{p[:4]}Q{(int(p[4:6])-1)//3+1}"
                except Exception:
                    pass
        return yq, rows
    # old-format fallback: the .txt primary doc
    txtn = next((n for n in names if n.endswith(".txt") and "index" not in n), None)
    if txtn:
        txt = get(f"https://www.sec.gov/Archives/edgar/data/{int(cik)}/{accn}/{txtn}").text
        if yq is None:
            m = re.search(r"CONFORMED PERIOD OF REPORT:\s*(\d{8})", txt)
            if m:
                p = m.group(1)
                yq = f"{p[:4]}Q{(int(p[4:6])-1)//3+1}"
        return yq, parse_old_txt(txt)
    return yq, []

def main():
    log = open("/tmp/hf13f_build.log", "a")
    def say(s):
        print(s, flush=True); log.write(s + "\n"); log.flush()
    allp = []
    for mi, (name, cik) in enumerate(MANAGERS):
        fl = filing_list(cik)
        say(f"[{mi+1}/{len(MANAGERS)}] {name}: {len(fl)} 13F-HR filings")
        time.sleep(0.3)
        for dt, acc in sorted(fl):
            cache = os.path.join(RAW, f"{cik}_{acc.replace('-','')}.csv")
            if os.path.exists(cache):
                df = pd.read_csv(cache, dtype={"cusip": str})
                if len(df):
                    allp.append((name, cik, dt, df["yq"].iloc[0], df))
                continue
            try:
                yq, rows = parse_filing(cik, dt, acc)
                df = pd.DataFrame(rows, columns=["cusip", "shares"])
                df["yq"] = yq
                df.to_csv(cache, index=False)
                if yq and len(df):
                    allp.append((name, cik, dt, yq, df))
                say(f"    {dt} {acc} -> {yq} {len(df)} positions")
            except Exception as e:
                say(f"    {dt} {acc} ERR {str(e)[:100]}")
            time.sleep(0.25)
    # dedupe per (manager, yq): keep latest filingDate
    best = {}
    for name, cik, dt, yq, df in allp:
        k = (name, yq)
        if k not in best or dt > best[k][0]:
            best[k] = (dt, df)
    say(f"filing-quarters kept: {len(best)}")
    # aggregate per (yq, cusip)
    recs = []
    for (name, yq), (dt, df) in best.items():
        g = df.groupby("cusip")["shares"].sum()
        for cusip, sh in g.items():
            recs.append((yq, name, cusip, int(sh)))
    h = pd.DataFrame(recs, columns=["yq", "manager", "cusip", "shares"])
    h["cusip"] = h["cusip"].astype(str).str.zfill(9)
    h.to_csv(os.path.join(BASE, "data", "hf_holdings_raw.csv"), index=False)
    say(f"raw holdings rows: {len(h)}, unique cusips: {h['cusip'].nunique()}")
    # CUSIP -> ticker via OpenFIGI (batch)
    cusips = sorted(h["cusip"].unique())
    say(f"mapping {len(cusips)} CUSIPs via OpenFIGI...")
    cmap = {}
    for i in range(0, len(cusips), 100):
        batch = [{"idType": "ID_CUSIP", "idValue": c} for c in cusips[i:i+100]]
        for attempt in range(3):
            try:
                r = requests.post("https://api.openfigi.com/v3/mapping", headers={**H, "Content-Type": "application/json"},
                                  data=json.dumps(batch), timeout=60)
                res = r.json()
                for c, jr in zip(cusips[i:i+100], res):
                    if isinstance(jr, dict) and jr.get("data"):
                        d0 = jr["data"][0]
                        if d0.get("securityType") in ("Common Stock", "ADR", "REIT", "Units"):
                            cmap[c] = (d0.get("ticker"), d0.get("name"))
                break
            except Exception as e:
                say(f"    figi batch {i} attempt {attempt} ERR {str(e)[:80]}")
                time.sleep(5)
        time.sleep(2.5)
    pd.DataFrame([(c, v[0], v[1]) for c, v in cmap.items()], columns=["cusip", "ticker", "name"])\
      .to_csv(os.path.join(BASE, "data", "hf_cusip_map.csv"), index=False)
    say(f"mapped {len(cmap)}/{len(cusips)} CUSIPs")
    h["ticker"] = h["cusip"].map({c: v[0] for c, v in cmap.items()})
    h = h.dropna(subset=["ticker"])
    # aggregate per (yq, ticker)
    agg = h.groupby(["yq", "ticker"]).agg(
        n_hf=("manager", "nunique"),
        hf_shares=("shares", "sum"),
        hf_hhi=("shares", lambda s: float(((s / s.sum()) ** 2).sum())),
    ).reset_index()
    # panel merge: shares outstanding + ADV at quarter-end
    p = pd.read_csv(os.path.join(BASE, "data", "panel_v2.csv"),
                    usecols=["month", "ticker", "adjclose", "mcap", "volume"])
    p["month"] = pd.to_datetime(p["month"])
    p["yq"] = p["month"].dt.to_period("Q").astype(str).str.replace("Q", "Q")
    # quarter-end month per yq
    qend = p.sort_values("month").groupby(["ticker", "yq"]).tail(1)[["ticker", "yq", "adjclose", "mcap", "volume"]]
    qend = qend.rename(columns={"volume": "vol_qend"})
    # ADV: mean volume over quarter
    adv = p.groupby(["ticker", "yq"])["volume"].mean().reset_index().rename(columns={"volume": "adv"})
    m = agg.merge(qend, on=["ticker", "yq"], how="left").merge(adv, on=["ticker", "yq"], how="left")
    m["shrout"] = m["mcap"] / m["adjclose"]
    m["hf_own"] = m["hf_shares"] / m["shrout"]
    m["hf_adv_ratio"] = m["hf_shares"] / m["adv"]
    # n managers filing that quarter (for coverage diagnostics)
    cov = pd.DataFrame([(yq, name) for (name, yq) in best]).drop_duplicates()
    cov = cov.groupby(0).size().reset_index().rename(columns={0: "yq", 1: "n_mgr_filed"})
    m = m.merge(cov, on="yq", how="left")
    m.to_csv(os.path.join(BASE, "data", "hf_crowding_quarterly.csv"), index=False)
    say(f"crowding panel: {len(m)} rows, {m['yq'].nunique()} quarters, tickers {m['ticker'].nunique()}")
    say(f"quarters: {m['yq'].min()} .. {m['yq'].max()}")
    say(json.dumps(m.groupby("yq")["n_mgr_filed"].first().to_dict()))
    log.close()

if __name__ == "__main__":
    main()
