"""Systematic Buffett screen: quality + value + low-beta, annual June rebalance.

Academic grounding: Frazzini, Kabiller & Pedersen (2018), "Buffett's Alpha" —
Berkshire ~= value (HML) + quality (QMJ) + betting-against-beta (BAB), levered
~1.6x via insurance float. This implements the UNLEVERED stock-selection leg.
No float leverage applied (documented limitation).

Point-in-time: June Y rebalance uses 10-K facts with filed <= Jun 30 Y
(standard 6-mo lag on annuals). Scores are cross-sectional z, winsorized ±3.
Top 30 equal-weighted. Monthly total returns from panel adjclose.
"""
import os, json, glob
import numpy as np
import pandas as pd
from datetime import date, timedelta

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RAW = os.path.join(DATA, 'raw')
FIG = os.path.join(HERE, '..', 'figs')
os.makedirs(FIG, exist_ok=True)

N_PICKS = 30
REB_YEARS = list(range(2008, 2022))  # rebalance end-June Y, held Jul Y - Jun Y+1


def zscore_winsor(s):
    s = s.astype(float)
    mu, sd = s.mean(), s.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return pd.Series(0.0, index=s.index)
    z = (s - mu) / sd
    return z.clip(-3, 3)


def load_panel():
    p = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'),
                    usecols=['month', 'ticker', 'mcap', 'adjclose'])
    p['month'] = pd.to_datetime(p['month'])
    return p


def load_fund():
    f = pd.read_csv(os.path.join(DATA, 'fundamentals_annual.csv'))
    f['fy_end'] = pd.to_datetime(f['fy_end'])
    for k in ('NetIncomeLoss', 'Equity', 'Assets'):
        f[k + '_filed'] = pd.to_datetime(f[k + '_filed'])
    return f


def pit_fundamentals(f, ticker, asof):
    """End-date timing convention (Fama-French standard): latest fiscal year
    ending >=120 days before asof (10-K filed by then in practice), end within
    18 months. Filed-date PIT is unreliable in companyfacts for pre-2010 years
    (restated facts carry the restating filing's date); verified that for
    2010+ this selection coincides with filed<=asof PIT."""
    asof = pd.Timestamp(asof)
    cutoff = asof - timedelta(days=120)
    sub = f[f['ticker'] == ticker].copy()
    sub = sub[sub['NetIncomeLoss'].notna() & sub['Equity'].notna() & sub['Assets'].notna()]
    sub = sub[(sub['fy_end'] <= cutoff) & ((asof - sub['fy_end']).dt.days <= 548)]
    if len(sub) == 0:
        return None
    return sub.sort_values('fy_end').iloc[-1]


# ---- daily data cache for betas ----
_daily_cache = {}
_spy_daily = None

def get_spy_daily():
    global _spy_daily
    if _spy_daily is None:
        s = pd.read_csv(os.path.join(RAW, 'SPY.csv'), parse_dates=['date'])
        s = s.sort_values('date')
        s['lr'] = np.log(s['adjclose'] / s['adjclose'].shift(1))
        _spy_daily = s[['date', 'lr']].rename(columns={'lr': 'spy_lr'})
    return _spy_daily


def beta_252(ticker, asof):
    """Trailing 252-trading-day market beta vs SPY. Needs >=126 overlapping days."""
    asof = pd.Timestamp(asof)
    if ticker not in _daily_cache:
        fp = os.path.join(RAW, f'{ticker}.csv')
        if not os.path.exists(fp):
            _daily_cache[ticker] = None
        else:
            d = pd.read_csv(fp, parse_dates=['date']).sort_values('date')
            d['lr'] = np.log(d['adjclose'] / d['adjclose'].shift(1))
            _daily_cache[ticker] = d[['date', 'lr']]
    d = _daily_cache[ticker]
    if d is None:
        return np.nan
    win = d[(d['date'] > asof - timedelta(days=370)) & (d['date'] <= asof)]
    m = win.merge(get_spy_daily(), on='date', how='inner').dropna()
    if len(m) < 126:
        return np.nan
    cov = np.cov(m['lr'], m['spy_lr'], ddof=0)[0, 1]
    var = m['spy_lr'].var(ddof=0)
    return cov / var if var > 0 else np.nan


def build_scores(p, f):
    """Score table per rebalance year. Returns dict year -> DataFrame of picks."""
    months = sorted(p['month'].unique())
    june_months = [m for m in months if m.month == 6 and m.year in REB_YEARS]
    picks, diag = {}, []
    for jm in june_months:
        Y = jm.year
        asof = jm  # last trading day of June
        mem = p[p['month'] == jm][['ticker', 'mcap', 'adjclose']].copy()
        mem = mem[mem['mcap'] > 0]
        rows = []
        for _, r in mem.iterrows():
            t = r['ticker']
            pf = pit_fundamentals(f, t, asof)
            if pf is None or pf['Equity'] <= 0:
                continue
            b = beta_252(t, asof)
            if np.isnan(b):
                continue
            rows.append(dict(ticker=t,
                             ROE=pf['NetIncomeLoss'] / pf['Equity'],
                             ROA=pf['NetIncomeLoss'] / pf['Assets'],
                             BM=pf['Equity'] / r['mcap'],
                             EP=pf['NetIncomeLoss'] / r['mcap'],
                             beta=b))
        sc = pd.DataFrame(rows).set_index('ticker')
        diag.append(dict(year=Y, june_members=len(mem), eligible=len(sc)))
        if len(sc) < N_PICKS:
            picks[Y] = []
            print(f'{Y}: members={len(mem)} eligible={len(sc)} -- SKIPPED (<{N_PICKS})', flush=True)
            continue
        q_z = (zscore_winsor(sc['ROE']) + zscore_winsor(sc['ROA'])) / 2
        v_z = (zscore_winsor(sc['BM']) + zscore_winsor(sc['EP'])) / 2
        b_z = -zscore_winsor(sc['beta'])
        comp = (q_z + v_z + b_z) / 3
        picks[Y] = comp.sort_values(ascending=False).head(N_PICKS).index.tolist()
        print(f'{Y}: members={len(mem)} eligible={len(sc)}', flush=True)
    pd.DataFrame(diag).to_csv(os.path.join(DATA, 'buffett_eligibility.csv'), index=False)
    return picks


def monthly_returns(p):
    piv = p.pivot(index='month', columns='ticker', values='adjclose').sort_index()
    return piv.pct_change()


def backtest(p, picks):
    rets = monthly_returns(p)
    months = rets.index
    port_rets, holdings_hist, turnovers = [], {}, []
    prev_w = None
    for Y in REB_YEARS:
        held = picks[Y]
        if not held:
            continue
        # holding window: Jul Y .. Jun Y+1
        start = next(m for m in months if m.year == Y and m.month == 7)
        end = next(m for m in months if m.year == Y + 1 and m.month == 6)
        window = months[(months >= start) & (months <= end)]
        w = pd.Series(1.0 / len(held), index=held)
        if prev_w is not None:
            overlap = w.reindex(prev_w.index).fillna(0).add(
                prev_w.reindex(w.index).fillna(0), fill_value=0)
            # one-way turnover = share of new portfolio not overlapping old (post-drift w)
            inter = set(w.index) & set(prev_w.index)
            turnovers.append(1 - prev_w.reindex(list(inter)).fillna(0).sum())
        for m in window:
            r = rets.loc[m, w.index]
            alive = r.dropna()
            if len(alive) == 0:
                port_rets.append((m, 0.0))
                w = pd.Series(dtype=float)
                continue
            w = w.reindex(alive.index)
            w = w / w.sum()  # delisted: exit at last price, redistribute
            port_rets.append((m, float((w * alive).sum())))
            # drift weights
            w = w * (1 + alive)
            w = w / w.sum()
        prev_w = w
        holdings_hist[Y] = held
    pr = pd.DataFrame(port_rets, columns=['month', 'ret']).set_index('month').sort_index()
    json.dump({str(k): v for k, v in holdings_hist.items()},
              open(os.path.join(DATA, 'buffett_holdings.json'), 'w'), indent=1)
    return pr, float(np.mean(turnovers)) if turnovers else np.nan


def spy_monthly():
    s = pd.read_csv(os.path.join(RAW, 'SPY.csv'), parse_dates=['date']).sort_values('date')
    m = s.set_index('date')['adjclose'].resample('M').last().pct_change()
    m.index = m.index.to_period('M')  # period index: aligns with panel months
    return m


def to_period_index(df):
    df = df.copy()
    if not isinstance(df.index, pd.PeriodIndex):
        df.index = df.index.to_period('M')
    return df


def stats(pr, spy, rf):
    pr = to_period_index(pr)
    r = pr['ret'].dropna()
    s = spy.reindex(r.index).dropna()
    r = r.reindex(s.index)
    exc = r - rf.reindex(r.index).fillna(0)
    eq = (1 + r).cumprod()
    dd = (eq / eq.cummax() - 1).min()
    n = len(r)
    cagr = (1 + r).prod() ** (12 / n) - 1
    sharpe = exc.mean() / exc.std() * np.sqrt(12) if exc.std() > 0 else np.nan
    # calendar-year hit rate vs SPY
    yr = pd.DataFrame({'p': r, 's': s})
    yh = (yr.groupby(yr.index.year)['p'].apply(lambda x: (1 + x).prod() - 1) >
          yr.groupby(yr.index.year)['s'].apply(lambda x: (1 + x).prod() - 1))
    return dict(n_months=n, CAGR=float(cagr), Sharpe=float(sharpe),
                MaxDD=float(dd), ann_vol=float(r.std() * np.sqrt(12)),
                hit_rate=float(yh.mean()), years=int(yh.sum()), years_n=int(len(yh)),
                final_eq=float(eq.iloc[-1]))


def load_ff():
    ff = pd.read_csv(os.path.join(DATA, 'ff_factors_monthly.csv'), skiprows=4)
    ff.columns = ['ym', 'MktRf', 'SMB', 'HML', 'RF']
    ff['ym'] = ff['ym'].astype(str).str.strip()
    ff = ff[ff['ym'].str.match(r'^\d{6}$')]
    ff['date'] = pd.to_datetime(ff['ym'], format='%Y%m') + pd.offsets.MonthEnd(0)
    ff = ff.set_index('date')[['MktRf', 'SMB', 'HML', 'RF']].astype(float) / 100
    ff.index = ff.index.to_period('M')  # period index: aligns with panel months
    return ff


def ff_regression(pr, spy):
    ff = load_ff()
    prp = to_period_index(pr)
    r = prp['ret'].reindex(ff.index).dropna()
    X = ff.reindex(r.index)[['MktRf', 'SMB', 'HML']]
    y = r - ff.reindex(r.index)['RF']
    Xc = pd.concat([pd.Series(1.0, index=X.index, name='alpha'), X], axis=1)
    beta, *_ = np.linalg.lstsq(Xc.values, y.values, rcond=None)
    resid = y.values - Xc.values @ beta
    # HC1 standard errors
    n, k = Xc.shape
    XtXinv = np.linalg.inv(Xc.values.T @ Xc.values)
    meat = (Xc.values * resid[:, None]).T @ (Xc.values * resid[:, None])
    V = XtXinv @ meat @ XtXinv * (n / (n - k))
    se = np.sqrt(np.diag(V))
    out = {}
    for i, name in enumerate(['alpha', 'MktRf', 'SMB', 'HML']):
        out[name] = dict(coef=float(beta[i]), se=float(se[i]),
                         t=float(beta[i] / se[i]) if se[i] > 0 else np.nan)
    out['R2'] = float(1 - (resid ** 2).sum() / ((y - y.mean()) ** 2).sum())
    out['alpha_ann'] = float(beta[0] * 12)
    out['n'] = n
    return out


def main():
    print('loading...', flush=True)
    p = load_panel()
    f = load_fund()
    print('scoring...', flush=True)
    picks = build_scores(p, f)
    print('backtesting...', flush=True)
    pr, turnover = backtest(p, picks)
    spy = spy_monthly()  # period-indexed
    rf = load_ff()['RF']  # period-indexed
    prp = to_period_index(pr)
    sp = spy.reindex(prp.index)
    res = dict(
        portfolio=stats(prp, sp, rf),
        spy=stats(pd.DataFrame({'ret': sp}), sp, rf),
        turnover_oneway=float(turnover),
        factor_regression=ff_regression(prp, spy),
    )
    pr.to_csv(os.path.join(DATA, 'buffett_portfolio_monthly.csv'))
    json.dump(res, open(os.path.join(DATA, 'buffett_results.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1)[:2000])
    # figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    eq_p = (1 + prp['ret']).cumprod()
    eq_s = (1 + sp).cumprod()
    fig, ax = plt.subplots(figsize=(10, 5))
    x = [pd.Timestamp(p.year, p.month, 1) + pd.offsets.MonthEnd(0) for p in eq_p.index]
    ax.plot(x, eq_p.values, label='Buffett screen (top-30, ann. rebalance)')
    ax.plot(x, eq_s.values, label='SPY total return')
    ax.set_yscale('log')
    ax.set_title('Systematic Buffett screen vs SPY (paper portfolio, 2009-2022)')
    ax.set_ylabel('Growth of $1 (log scale)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'buffett_screen_equity.png'), dpi=120)
    print('figure saved')

if __name__ == '__main__':
    main()
