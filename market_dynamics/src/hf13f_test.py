"""Pre-registered tests: 13F HF crowding -> forward returns.

Design frozen in docs/HFCROWDING_TEST.md BEFORE this ran.
Universe: point-in-time S&P 500 members (panel_v2.csv).
Timing: holdings quarter q -> formation at close of month q+2 (45d 13F lag) ->
  1m and 3m forward buy-and-hold returns.
Split: estimate formations 2000-2007; walk-forward 2008-2021.
Pass: (a) IS |NW t|>2.0 positive sign; (b) WF mean>0 positive; (c) WF hit>0.5.
"""
import pandas as pd, numpy as np, json, os
from scipy import stats

BASE = os.path.expanduser("~/workspace/market_dynamics")

def q_to_form(q):
    y, qq = q.split("Q"); qq = int(qq)
    m = qq * 3 + 2  # q_end month + 2
    y = int(y) + (m - 1) // 12; m = (m - 1) % 12 + 1
    return pd.Timestamp(y, m, 1) + pd.offsets.MonthEnd(0)

def nw_t(x, lags):
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 10:
        return np.nan, np.nan
    mu = x.mean()
    xc = x - mu
    v = (xc ** 2).sum() / n
    for L in range(1, lags + 1):
        w = 1 - L / (lags + 1)
        v += 2 * w * (xc[L:] * xc[:-L]).sum() / n
    se = np.sqrt(v / n)
    return mu, mu / se if se > 0 else np.nan

def main():
    crowd = pd.read_csv(f"{BASE}/data/hf_crowding_quarterly.csv")
    crowd = crowd[["yq", "ticker", "n_hf", "hf_shares", "hf_own", "hf_adv_ratio", "hf_hhi"]]
    p = pd.read_csv(f"{BASE}/data/panel_v2.csv", usecols=["month", "ticker", "adjclose", "mcap"])
    p["month"] = pd.to_datetime(p["month"])
    p = p.sort_values(["ticker", "month"])
    p["ret1"] = p.groupby("ticker")["adjclose"].pct_change().shift(-1)
    p["ret3"] = p.groupby("ticker")["adjclose"].pct_change(3).shift(-3)
    # formation month per holdings quarter
    crowd["form"] = crowd["yq"].map(q_to_form)
    form_months = p[["ticker", "month"]].rename(columns={"month": "form"})
    # forward returns at formation
    fwd = p[["ticker", "month", "ret1", "ret3", "mcap"]].rename(columns={"month": "form"})
    d = crowd.merge(fwd, on=["ticker", "form"], how="inner")
    d = d.dropna(subset=["ret1"])
    print(f"formation-stock obs: {len(d)}, formations: {d['form'].nunique()}")
    measures = ["n_hf", "hf_own", "hf_adv_ratio", "hf_hhi"]
    res = {}
    for meas in measures:
        for h, rcol in [("1m", "ret1"), ("3m", "ret3")]:
            dd = d.dropna(subset=[rcol, meas]).copy()
            # quintile by measure within formation
            dd["q"] = dd.groupby("form")[meas].transform(
                lambda s: pd.qcut(s.rank(method="first"), 5, labels=False, duplicates="drop"))
            spr = []
            for fm, g in dd.groupby("form"):
                q5 = g[g.q == 4]; q1 = g[g.q == 0]
                if len(q5) < 3 or len(q1) < 3:
                    continue
                r5 = (q5[rcol] * q5["mcap"] / q5["mcap"].sum()).sum()
                r1 = (q1[rcol] * q1["mcap"] / q1["mcap"].sum()).sum()
                spr.append((fm, r5 - r1,
                            q5[rcol].mean() - q1[rcol].mean()))
            s = pd.DataFrame(spr, columns=["form", "vw", "ew"]).sort_values("form")
            s["yr"] = s["form"].dt.year
            out = {}
            for w in ["vw", "ew"]:
                x = s[w].dropna().values
                mu, t = nw_t(x, 6 if h == "1m" else 3)
                per = 12 if h == "1m" else 4
                isx = s[s.yr <= 2007][w].dropna().values
                wfx = s[s.yr >= 2008][w].dropna().values
                mu_is, t_is = nw_t(isx, 6 if h == "1m" else 3)
                mu_wf, t_wf = nw_t(wfx, 6 if h == "1m" else 3)
                out[w] = {
                    "n": len(x), "ann_mean": float(mu * per),
                    "t_nw": float(t),
                    "is_ann": float(mu_is * per), "is_t": float(t_is), "is_n": len(isx),
                    "wf_ann": float(mu_wf * per), "wf_t": float(t_wf), "wf_n": len(wfx),
                    "wf_hit": float((wfx > 0).mean()) if len(wfx) else None,
                }
            res[f"{meas}_{h}"] = out
            print(f"{meas} {h}: VW ann {out['vw']['ann_mean']:+.3f} t={out['vw']['t_nw']:+.2f} | "
                  f"IS {out['vw']['is_ann']:+.3f} t={out['vw']['is_t']:+.2f} | "
                  f"WF {out['vw']['wf_ann']:+.3f} hit={out['vw']['wf_hit']}")
    # Fama-MacBeth: forward 3m return on rank[0,1]
    fm = []
    for h, rcol in [("1m", "ret1"), ("3m", "ret3")]:
        for meas in measures:
            dd = d.dropna(subset=[rcol, meas]).copy()
            dd["rk"] = dd.groupby("form")[meas].rank(pct=True)
            betas = []
            for fm_, g in dd.groupby("form"):
                if len(g) < 20:
                    continue
                X = np.column_stack([np.ones(len(g)), g["rk"].values])
                y = g[rcol].values
                b, *_ = np.linalg.lstsq(X, y, rcond=None)
                betas.append(b[1])
            betas = np.array(betas)
            mu, t = nw_t(betas, 6 if h == "1m" else 3)
            per = 12 if h == "1m" else 4
            fm.append({"h": h, "meas": meas, "ann_slope": float(mu * per),
                       "t": float(t), "n": len(betas)})
            print(f"FM {meas} {h}: slope ann {mu*per:+.4f} t={t:+.2f} n={len(betas)}")
    # Fama-MacBeth with log-size control (pre-registered robustness)
    fmc = []
    for h, rcol in [("1m", "ret1"), ("3m", "ret3")]:
        for meas in measures:
            dd = d.dropna(subset=[rcol, meas, "mcap"]).copy()
            dd = dd[dd["mcap"] > 0]
            dd["rk"] = dd.groupby("form")[meas].rank(pct=True)
            dd["lsz"] = np.log(dd["mcap"])
            betas = []
            for fm_, g in dd.groupby("form"):
                if len(g) < 30:
                    continue
                X = np.column_stack([np.ones(len(g)), g["rk"].values, g["lsz"].values])
                y = g[rcol].values
                b, *_ = np.linalg.lstsq(X, y, rcond=None)
                betas.append(b[1])
            betas = np.array(betas)
            mu, t = nw_t(betas, 6 if h == "1m" else 3)
            per = 12 if h == "1m" else 4
            fmc.append({"h": h, "meas": meas, "ann_slope": float(mu * per),
                        "t": float(t), "n": len(betas)})
            print(f"FM+size {meas} {h}: slope ann {mu*per:+.4f} t={t:+.2f} n={len(betas)}")
    json.dump({"spreads": res, "fm": fm, "fm_size_control": fmc},
              open(f"{BASE}/data/hf_crowding_test_results.json", "w"), indent=1)
    print("wrote data/hf_crowding_test_results.json")

if __name__ == "__main__":
    main()
