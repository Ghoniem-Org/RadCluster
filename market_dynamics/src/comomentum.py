"""COMOMENTUM (Lou & Polk 2022, RFS): build the series and run the pre-registered test.

Design frozen in docs/COMOMENTUM_TEST.md BEFORE this script ran.
Deviations from Lou & Polk are documented there (lane cohorts, 12-1 mom,
S&P 500 universe, no delisting returns).

Inputs (labels):
  data/panel_v2.csv      measured: membership + mom (12-1m return on adjclose)
  data/raw/*.csv         measured: Yahoo daily adjclose
  data/ff3_daily.zip     measured (third-party): Ken French FF3 daily factors,
                         downloaded 2026-09-28 17:01 UTC from
                         https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_daily_CSV.zip
Assumed: 252-day window, >=200-day rule, 15-stock floor, lane cohorts.
"""
import os
import json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

import sys
sys.path.insert(0, HERE)
from model4 import N_BIN, assign_bins, cap_state, mom_of

TOP_BINS = [b for b in range(N_BIN) if mom_of(b) == 2]
WIN_EDGE, LOSE_EDGE = 0.30, 0.0   # model's absolute momentum lane edges
MIN_DAYS = 200                     # assumed: min daily obs in 252-day window
MIN_STOCKS = 15                    # assumed: min qualifying stocks per cohort
NW_LAGS = 11                       # Newey-West lags for 12-mo overlapping targets


# ---------------- data ----------------
def load_daily():
    """ticker -> DataFrame(date, adjclose) sorted."""
    import glob
    out = {}
    for f in glob.glob(os.path.join(DATA, 'raw', '*.csv')):
        t = os.path.basename(f)[:-4]
        try:
            d = pd.read_csv(f, usecols=['date', 'adjclose'], parse_dates=['date'])
        except Exception:
            continue
        d = d.dropna().sort_values('date')
        if len(d) > 50:
            out[t] = d.reset_index(drop=True)
    return out


def load_ff3():
    z = os.path.join(DATA, 'ff3_daily.zip')
    df = pd.read_csv(z, skiprows=3)
    df = df.rename(columns={'Unnamed: 0': 'yyyymmdd'})
    df = df[pd.to_numeric(df['yyyymmdd'], errors='coerce').notna()].copy()
    df['date'] = pd.to_datetime(df['yyyymmdd'].astype(int).astype(str), format='%Y%m%d')
    for c in ['Mkt-RF', 'SMB', 'HML', 'RF']:
        df[c] = df[c].astype(float) / 100.0
    return df.set_index('date').sort_index()[['Mkt-RF', 'SMB', 'HML', 'RF']]


def pairwise_corr(M):
    """Pairwise-complete Pearson correlations of columns of M (rows may have NaN). Exact, vectorized."""
    mu = np.nanmean(M, axis=0)
    sd = np.nanstd(M, axis=0)
    sd[sd == 0] = np.nan
    Z = (M - mu) / sd
    Z[~np.isfinite(Z)] = 0.0
    V = np.isfinite(M).astype(float)
    cnt = V.T @ V
    C = (Z.T @ Z) / np.maximum(cnt, 1)
    np.fill_diagonal(C, 1.0)
    iu = np.triu_indices(C.shape[0], 1)
    return float(np.mean(C[iu])) if len(iu[0]) else np.nan, int(C.shape[0])


def resid_corr(tickers, daily, ff3, T):
    """Mean pairwise FF3-residual correlation for a cohort at month-end T."""
    cal = ff3.index[ff3.index <= T][-252:]
    if len(cal) < MIN_DAYS:
        return np.nan, 0
    f = ff3.loc[cal]
    X = np.column_stack([np.ones(len(cal)), f['Mkt-RF'].values,
                         f['SMB'].values, f['HML'].values])
    cols = []
    for t in tickers:
        d = daily.get(t)
        if d is None:
            continue
        s = d[d['date'] <= T].tail(252)
        if len(s) < MIN_DAYS:
            continue
        r = s.set_index('date')['adjclose'].pct_change()
        r = r.reindex(cal)
        if r.notna().sum() < MIN_DAYS:
            continue
        y = (r - f['RF']).values
        ok = np.isfinite(y)
        if ok.sum() < MIN_DAYS:
            continue
        b, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
        e = np.full(len(cal), np.nan)
        e[ok] = y[ok] - X[ok] @ b
        cols.append(e)
    if len(cols) < MIN_STOCKS:
        return np.nan, len(cols)
    M = np.column_stack(cols)
    return pairwise_corr(M)


def buy_hold(tickers, daily, T):
    """Equal-weighted 12-month buy-and-hold from T; per-stock needs >=200 trading days."""
    E = T + pd.DateOffset(months=12)
    rets = []
    for t in tickers:
        d = daily.get(t)
        if d is None:
            continue
        pre = d[d['date'] <= T]
        if len(pre) == 0:
            continue
        p0 = float(pre.iloc[-1]['adjclose'])
        win = d[(d['date'] > T) & (d['date'] <= E)]
        if len(win) < MIN_DAYS or p0 <= 0:
            continue
        p1 = float(win.iloc[-1]['adjclose'])
        if p1 <= 0:
            continue
        rets.append(p1 / p0 - 1)
    if len(rets) < MIN_STOCKS:
        return np.nan, len(rets)
    return float(np.mean(rets)), len(rets)


# ---------------- Newey-West ----------------
def nw_tstats(X, y, L):
    X = np.asarray(X, float); y = np.asarray(y, float)
    n, k = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)
    e = y - X @ beta
    S = np.zeros((k, k))
    for l in range(L + 1):
        w = 1.0 - l / (L + 1)
        a = (X[l:] * e[l:, None]).T @ (X[:n - l] * e[:n - l, None])
        S += w * (a + a.T) if l > 0 else a
    V = XtX_inv @ S @ XtX_inv
    se = np.sqrt(np.diag(V))
    r2 = 1 - e @ e / ((y - y.mean()) @ (y - y.mean()))
    return beta, se, beta / se, r2


def main():
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'), parse_dates=['month'])
    months = sorted(panel['month'].unique())
    print('panel months:', months[0].date(), '->', months[-1].date(), flush=True)
    daily = load_daily()
    print('daily tickers:', len(daily), flush=True)
    ff3 = load_ff3()
    print('ff3:', ff3.index.min().date(), '->', ff3.index.max().date(), flush=True)

    rows = []
    skipped = []
    for m in months:
        T = pd.Timestamp(m)
        if (ff3.index <= T).sum() < 252:
            continue  # no full 252-day window yet
        g = panel[panel['month'] == m]
        gm = g[g['mom'].notna()]
        win_t = gm[gm['mom'] >= WIN_EDGE]['ticker'].tolist()
        lose_t = gm[gm['mom'] < LOSE_EDGE]['ticker'].tolist()
        cw, nw_ = resid_corr(win_t, daily, ff3, T)
        cl, nl = resid_corr(lose_t, daily, ff3, T)
        if np.isnan(cw) or np.isnan(cl):
            skipped.append(str(m.date()))
            continue
        # robustness cohort: cross-sectional terciles (pre-registered)
        q = gm['mom'].quantile([1 / 3, 2 / 3])
        tw = gm[gm['mom'] >= q.iloc[1]]['ticker'].tolist()
        tl = gm[gm['mom'] <= q.iloc[0]]['ticker'].tolist()
        ctw, _ = resid_corr(tw, daily, ff3, T)
        ctl, _ = resid_corr(tl, daily, ff3, T)
        # targets: 12m buy-and-hold
        rw, qnw = buy_hold(win_t, daily, T)
        rl, qnl = buy_hold(lose_t, daily, T)
        # rev-5 concentration controls
        d0 = g.copy()
        d0['bin'] = assign_bins(d0)
        c, _ = cap_state(d0[d0['bin'] >= 0])
        C = float(c[TOP_BINS].sum()); hhi = float((c ** 2).sum())
        rows.append(dict(month=m, comom_win=cw, comom_lose=cl,
                         comom=(cw + cl) / 2, n_win=nw_, n_lose=nl,
                         comom_terc_win=ctw, comom_terc_lose=ctl,
                         ret_win_12m=rw, ret_lose_12m=rl,
                         wml_12m=(rw - rl if np.isfinite(rw) and np.isfinite(rl) else np.nan),
                         C=C, hhi=hhi, n_qw=qnw, n_ql=qnl))
    s = pd.DataFrame(rows)
    s.to_csv(os.path.join(DATA, 'comomentum_monthly.csv'), index=False)
    print('comomentum months:', len(s), s['month'].min(), '->', s['month'].max(), flush=True)
    print('skipped months (<15/cohort):', len(skipped), skipped[:8], flush=True)
    print(s[['comom', 'comom_win', 'comom_lose', 'n_win', 'n_lose']].describe().to_string(), flush=True)

    # ---- pre-registered tests ----
    d = s.dropna(subset=['comom', 'ret_win_12m']).copy()
    est = d[d['month'] < '2008-01-01'].copy()
    print('\nin-sample n =', len(est), '| full n =', len(d), flush=True)
    X = np.column_stack([np.ones(len(est)), est['comom'].values])

    res = {}
    for target, name in [('ret_win_12m', 'win12'), ('wml_12m', 'wml12')]:
        dd = d.dropna(subset=[target])
        ee = dd[dd['month'] < '2008-01-01']
        Xe = np.column_stack([np.ones(len(ee)), ee['comom'].values])
        b, se, t, r2 = nw_tstats(Xe, ee[target].values, NW_LAGS)
        # component consistency
        comp = {}
        for cc in ['comom_win', 'comom_lose']:
            Xc = np.column_stack([np.ones(len(ee)), ee[cc].values])
            bc, sec, tc, _ = nw_tstats(Xc, ee[target].values, NW_LAGS)
            comp[cc] = dict(b=round(float(bc[1]), 4), t=round(float(tc[1]), 2))
        res[name] = dict(n=len(ee), b=round(float(b[1]), 4), t_nw=round(float(t[1]), 2),
                         r2=round(float(r2), 4), components=comp)
        print('%s: b=%.4f  NW-t=%.2f  R2=%.4f  comps=%s' % (name, b[1], t[1], r2, comp), flush=True)

    # tercile robustness (in-sample, descriptive per pre-reg)
    ee = d[d['month'] < '2008-01-01'].dropna(subset=['comom_terc_win'])
    ee['comom_terc'] = (ee['comom_terc_win'] + ee['comom_terc_lose']) / 2
    Xt = np.column_stack([np.ones(len(ee)), ee['comom_terc'].values])
    bt, set_, tt, r2t = nw_tstats(Xt, ee['ret_win_12m'].values, NW_LAGS)
    res['tercile_robust'] = dict(n=len(ee), b=round(float(bt[1]), 4), t_nw=round(float(tt[1]), 2))
    print('tercile cohort: b=%.4f NW-t=%.2f (n=%d)' % (bt[1], tt[1], len(ee)), flush=True)

    # ---- walk-forward (expanding, zero lookahead) ----
    wf_rows = []
    allm = sorted(d['month'].unique())
    for m in allm:
        if m < pd.Timestamp('2008-01-01'):
            continue
        cut = m - pd.DateOffset(months=12)
        tr = d[(d['month'] < m) & (d['month'] <= cut)]
        if len(tr) < 24:
            continue
        Xtr = np.column_stack([np.ones(len(tr)), tr['comom'].values])
        bh = np.linalg.lstsq(Xtr, tr['ret_win_12m'].values, rcond=None)[0]
        fc = bh[0] + bh[1] * d.loc[d['month'] == m, 'comom'].values[0]
        rl = d.loc[d['month'] == m, 'ret_win_12m'].values[0]
        bm = tr['ret_win_12m'].mean()
        wf_rows.append(dict(month=m, fc=float(fc), realized=float(rl), bench=float(bm)))
    wf = pd.DataFrame(wf_rows)
    corr = float(wf['fc'].corr(wf['realized']))
    oos_r2 = 1 - ((wf['realized'] - wf['fc']) ** 2).sum() / ((wf['realized'] - wf['bench']) ** 2).sum()
    hit = float((((wf['fc'] - wf['bench']) > 0) == ((wf['realized'] - wf['bench']) > 0)).mean())
    res['walkforward'] = dict(n=len(wf), corr=round(corr, 3), oos_r2=round(float(oos_r2), 4),
                              hit=round(hit, 3))
    print('walk-forward n=%d: corr=%.3f OOS_R2=%.4f hit=%.3f' % (len(wf), corr, oos_r2, hit), flush=True)

    # ---- horse race vs rev-5 concentration ----
    dd = d.dropna(subset=['ret_win_12m', 'C', 'hhi'])
    Xh = np.column_stack([np.ones(len(dd)), dd['comom'].values, dd['C'].values, dd['hhi'].values])
    bh, seh, th, r2h = nw_tstats(Xh, dd['ret_win_12m'].values, NW_LAGS)
    res['horserace'] = dict(n=len(dd),
                            comom_b=round(float(bh[1]), 4), comom_t=round(float(th[1]), 2),
                            C_b=round(float(bh[2]), 4), C_t=round(float(th[2]), 2),
                            hhi_b=round(float(bh[3]), 4), hhi_t=round(float(th[3]), 2))
    print('horse race (n=%d): comom b=%.4f t=%.2f | C b=%.4f t=%.2f | hhi b=%.4f t=%.2f'
          % (len(dd), bh[1], th[1], bh[2], th[2], bh[3], th[3]), flush=True)

    with open(os.path.join(DATA, 'comomentum_test_results.json'), 'w') as f:
        json.dump(res, f, indent=2)
    print('\nwrote data/comomentum_monthly.csv + data/comomentum_test_results.json', flush=True)


if __name__ == '__main__':
    main()
