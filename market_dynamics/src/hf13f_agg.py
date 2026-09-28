"""Aggregate 13F holdings -> quarterly crowding panel for S&P 500 members.

Measures per (holdings quarter, ticker):
  n_hf         # sampled managers holding
  hf_shares    total sampled-HF shares
  hf_own       hf_shares / shares outstanding (point-in-time, shares_use)
  hf_adv_ratio hf_shares / ADV (ADV = mean monthly volume over quarter / 21)
  hf_hhi       HHI of managers' shares of sampled-HF holdings
"""
import pandas as pd, numpy as np, os

BASE = os.path.expanduser("~/workspace/market_dynamics")

def main():
    h = pd.read_csv(f"{BASE}/data/hf_holdings_raw.csv", dtype={"cusip": str})
    mp = pd.read_csv(f"{BASE}/data/cusip_yahoo_map.csv", dtype=str).fillna("")
    mp = mp[mp["ticker"] != ""]
    print(f"CUSIP->ticker mapping rate: {len(mp)}/{h['cusip'].nunique()} "
          f"= {len(mp)/h['cusip'].nunique():.1%}")
    h = h.merge(mp[["cusip", "ticker"]], on="cusip", how="inner")
    panel = pd.read_csv(f"{BASE}/data/panel_v2.csv",
                        usecols=["month", "ticker", "volume", "shares_use", "mcap"])
    panel["month"] = pd.to_datetime(panel["month"])
    univ = set(panel["ticker"].unique())
    print(f"panel tickers: {len(univ)}")
    h = h[h["ticker"].isin(univ)].copy()
    # point-in-time membership: keep only (yq, ticker) actually in the panel that quarter
    panel["yq"] = panel["month"].dt.year.astype(str) + "Q" + \
        ((panel["month"].dt.month - 1) // 3 + 1).astype(str)
    members = set(zip(panel["yq"], panel["ticker"]))
    h = h[[(y, t) in members for y, t in zip(h["yq"], h["ticker"])]].copy()
    # pre-registered window 1999Q4-2021Q2, extended one quarter to 2021Q3 (real filings;
    # formation 2021-11-30 still inside the 2008-2021 walk-forward window) -- noted in doc
    h = h[(h["yq"] >= "1999Q4") & (h["yq"] <= "2021Q3")].copy()
    print(f"holdings rows in PIT panel universe, window: {len(h)}")
    g = h.groupby(["yq", "ticker"])
    agg = g.agg(n_hf=("manager", "nunique"), hf_shares=("shares", "sum")).reset_index()
    # HHI of manager shareholdings
    hh = h.copy()
    hh["tot"] = hh.groupby(["yq", "ticker"])["shares"].transform("sum")
    hh["s"] = hh["shares"] / hh["tot"]
    hhi = hh.groupby(["yq", "ticker"])["s"].apply(lambda s: (s ** 2).sum()).reset_index(name="hf_hhi")
    agg = agg.merge(hhi, on=["yq", "ticker"])
    # quarter-end month fundamentals
    qend = panel.sort_values("month").groupby(["yq", "ticker"]).tail(1) \
        [["yq", "ticker", "shares_use", "mcap"]].rename(columns={"shares_use": "shrout"})
    # ADV: trailing 3-month mean monthly volume / 21 trading days
    panel = panel.sort_values(["ticker", "month"])
    panel["adv"] = panel.groupby("ticker")["volume"].transform(
        lambda s: s.rolling(3, min_periods=1).mean() / 21)
    adv = panel.sort_values("month").groupby(["yq", "ticker"]).tail(1)[["yq", "ticker", "adv"]]
    agg = agg.merge(qend, on=["yq", "ticker"], how="left").merge(adv, on=["yq", "ticker"], how="left")
    agg["hf_own"] = agg["hf_shares"] / agg["shrout"]
    agg["hf_adv_ratio"] = agg["hf_shares"] / agg["adv"]
    agg = agg.sort_values(["yq", "ticker"]).reset_index(drop=True)
    agg.to_csv(f"{BASE}/data/hf_crowding_quarterly.csv", index=False)
    print(f"crowding panel: {len(agg)} rows, {agg['yq'].nunique()} quarters "
          f"({agg['yq'].min()}->{agg['yq'].max()})")
    print("median n_hf:", agg["n_hf"].median(), "| median hf_own:",
          round(agg["hf_own"].median(), 4), "| median hf_adv_ratio:",
          round(agg["hf_adv_ratio"].median(), 2))
    print("quarters with <5 managers:", (agg.groupby("yq")["n_hf"].max() < 5).sum())
    print("hf_own>0.5 frac:", (agg["hf_own"] > 0.5).mean(), "| hf_own>1 frac:", (agg["hf_own"] > 1).mean())
    print("missing shrout:", agg["shrout"].isna().mean(), "| missing adv:", agg["adv"].isna().mean())

if __name__ == "__main__":
    main()
