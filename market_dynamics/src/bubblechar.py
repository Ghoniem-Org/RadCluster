"""Bubble characteristics test (GSY 2019 adaptation) — implements the frozen design in
docs/BUBBLECHAR_TEST.md. Run AFTER the design doc was written; do not tune after seeing results."""
import json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

BASE = '/home/hatch/workspace/market_dynamics'

# ---------------- load ----------------
p = pd.read_csv(f'{BASE}/data/panel_v2.csv',
                usecols=['month', 'ticker', 'adjclose', 'has_price'])
p['month'] = pd.to_datetime(p['month'])
p = p.sort_values(['ticker', 'month']).reset_index(drop=True)

s = pd.read_csv(f'{BASE}/data/spx_monthly.csv')
s['ym'] = pd.to_datetime(s['ym'])
s = s.sort_values('ym').reset_index(drop=True)
s['Rmkt24'] = np.exp(s['lr'].rolling(24).sum()) - 1          # trailing 24m simple
s['Rmkt_fwd24'] = np.exp(s['lr'].shift(-1).rolling(24).sum()) - 1  # next 24m simple
mkt24 = dict(zip(s['ym'].dt.strftime('%Y-%m'), s['Rmkt24']))
mktfwd24 = dict(zip(s['ym'].dt.strftime('%Y-%m'), s['Rmkt_fwd24']))

# ---------------- per-ticker series ----------------
p['valid'] = (p['has_price'] == 1) & p['adjclose'].notna()
p['cvalid'] = p.groupby('ticker')['valid'].transform(
    lambda x: x.rolling(25, min_periods=25).sum())
p['lret'] = np.log(p['adjclose']).diff()  # no ffill artifacts (log of NA stays NA)
p.loc[p.groupby('ticker')['adjclose'].shift(1).isna(), 'lret'] = np.nan
p['px24'] = p.groupby('ticker')['adjclose'].shift(24)
p['px12'] = p.groupby('ticker')['adjclose'].shift(12)
p['R24'] = p['adjclose'] / p['px24'] - 1
p['R_yr1'] = p['px12'] / p['px24'] - 1          # first year of the 2-yr window
p['accel'] = p['R24'] - p['R_yr1']              # GSY acceleration (locked)
# volatility: trailing-12m std of monthly log returns -> cross-sectional pct rank
p['vol12'] = p.groupby('ticker')['lret'].transform(lambda x: x.rolling(12, min_periods=12).std())
p['vol_rank'] = p.groupby('month')['vol12'].rank(pct=True)
p['vol_rank12'] = p.groupby('ticker')['vol_rank'].shift(12)
p['dvol'] = p['vol_rank'] - p['vol_rank12']     # 1-yr change in vol rank (locked)
# age: months since first appearance in panel (left-censored 1999-12; documented proxy)
first = p.groupby('ticker')['month'].min().rename('first_month')
p = p.join(first, on='ticker')
p['age_m'] = ((p['month'] - p['first_month']).dt.days / 30.44).round().astype('Int64')
p['age_rank'] = p.groupby('month')['age_m'].rank(pct=True)
p['young'] = (p['age_m'] < 60).astype(int)
# market mapping
p['ym'] = p['month'].dt.strftime('%Y-%m')
p['Rmkt24'] = p['ym'].map(mkt24)
p['Rmkt_fwd24'] = p['ym'].map(mktfwd24)
p['Rnet'] = (1 + p['R24']) / (1 + p['Rmkt24']) - 1

# ---------------- episodes: first-crossing with 24m embargo ----------------
cand = p[(p['cvalid'] >= 25) & (p['R24'] > 1.0) & (p['Rnet'] > 1.0)
         & (p['month'] >= '2001-12') & (p['month'] <= '2020-12')
         & p['accel'].notna() & p['dvol'].notna() & p['age_rank'].notna()].copy()
cand = cand.sort_values(['ticker', 'month'])
eps, last = [], {}
for _, r in cand.iterrows():
    tk, m = r['ticker'], r['month']
    if tk not in last or (m - last[tk]).days >= 720:
        eps.append(r); last[tk] = m
e = pd.DataFrame(eps).reset_index(drop=True)
print(f'first-crossing episodes: {len(e)}')

# ---------------- forward outcomes ----------------
px = p.set_index(['ticker', 'month'])['adjclose']
hp = p.set_index(['ticker', 'month'])['has_price']
records = []
for _, r in e.iterrows():
    tk, m = r['ticker'], r['month']
    fwd_end = m + pd.offsets.DateOffset(months=24)
    sub = p[(p['ticker'] == tk) & (p['month'] > m) & (p['month'] <= fwd_end)
            & (p['has_price'] == 1) & p['adjclose'].notna()].sort_values('month')
    n_fwd = len(sub)
    crash, dd = np.nan, np.nan
    if n_fwd >= 12:
        path = sub['adjclose'].values / r['adjclose'] - 1
        dd = path.min()
        crash = int(dd <= -0.40)
    ret24_raw, ret24_net = np.nan, np.nan
    if n_fwd >= 24:
        ret24_raw = sub['adjclose'].iloc[23] / r['adjclose'] - 1
        rm = r['Rmkt_fwd24']
        if pd.notna(rm):
            ret24_net = (1 + ret24_raw) / (1 + rm) - 1
    records.append(dict(ticker=tk, month=m.strftime('%Y-%m-%d'),
                        R24=round(float(r['R24']), 4), Rnet=round(float(r['Rnet']), 4),
                        accel=round(float(r['accel']), 4), dvol=round(float(r['dvol']), 4),
                        age_rank=round(float(r['age_rank']), 4), age_m=int(r['age_m']),
                        young=int(r['young']), n_fwd=int(n_fwd),
                        crash=(None if pd.isna(crash) else int(crash)),
                        drawdown=None if pd.isna(dd) else round(float(dd), 4),
                        ret24_raw=None if pd.isna(ret24_raw) else round(float(ret24_raw), 4),
                        ret24_net=None if pd.isna(ret24_net) else round(float(ret24_net), 4)))
ep = pd.DataFrame(records)
ep.to_csv(f'{BASE}/data/bubblechar_episodes.csv', index=False)
print(f'episodes with >=12m forward: {ep["crash"].notna().sum()} '
      f'(dropped {(ep["crash"].isna()).sum()})')
print(f'episodes with 24m forward return: {ep["ret24_raw"].notna().sum()}')

d = ep[ep['crash'].notna()].copy()
print(f'\nN_crash_test = {len(d)}, crash rate = {d["crash"].mean():.3f} '
      f'(GSY ~0.53 at 100% run-up)')

# ---------------- (c) calibration bins ----------------
for lo, hi, lbl in [(1.0, 1.5, '100-150% net'), (1.5, 99, '>150% net')]:
    b = d[(d['Rnet'] >= lo) & (d['Rnet'] < hi)]
    if len(b):
        print(f'  bin {lbl}: n={len(b)}, crash rate={b["crash"].mean():.3f} '
              f'(GSY: 0.53 / 0.80)')

# ---------------- univariate crash vs non-crash (GSY Table 4 style) ----------------
uni = {}
for col, exp in [('accel', '+'), ('dvol', '+'), ('age_rank', '-')]:
    a = d.loc[d['crash'] == 1, col]; b = d.loc[d['crash'] == 0, col]
    t, pv = stats.ttest_ind(a, b, equal_var=False)
    uni[col] = dict(mean_crash=round(a.mean(), 4), mean_nocrash=round(b.mean(), 4),
                    diff=round(a.mean() - b.mean(), 4), t=round(t, 3), p=round(pv, 4),
                    expected_sign=exp,
                    sign_ok=bool((a.mean() - b.mean() > 0) == (exp == '+')))
    print(f'{col}: crash mean={a.mean():.4f} vs non-crash {b.mean():.4f}, '
          f't={t:.2f}, p={pv:.3f}, sign_ok={uni[col]["sign_ok"]}')

# ---------------- (a) crash logit (IRLS, no statsmodels available) ----------------
def fit_logit(X, y):
    X = np.column_stack([np.ones(len(X)), X])
    def nll(b):
        z = X @ b
        # stable log-likelihood
        ll = np.sum(y * z - np.logaddexp(0, z))
        return -ll
    res = minimize(nll, np.zeros(X.shape[1]), method='BFGS')
    b = res.x
    p = 1 / (1 + np.exp(-(X @ b)))
    W = p * (1 - p)
    XtWX = (X * W[:, None]).T @ X
    cov = np.linalg.inv(XtWX)
    se = np.sqrt(np.diag(cov))
    llf = -res.fun
    pbar = y.mean()
    llnull = len(y) * (pbar * np.log(pbar) + (1 - pbar) * np.log(1 - pbar))
    return b, se, llf, llnull

Xc = d[['accel', 'dvol', 'age_rank']].values
yc = d['crash'].values
b, se, llf, llnull = fit_logit(Xc, yc)
names = ['const', 'accel', 'dvol', 'age_rank']
print('--- crash logit (IRLS) ---')
for n_, bi, sei in zip(names, b, se):
    z = bi / sei
    pv = 2 * (1 - stats.norm.cdf(abs(z)))
    print(f'  {n_:10s} coef={bi:+.4f} se={sei:.4f} z={z:+.2f} p={pv:.4f}')
pseudo_r2 = float(1 - llf / llnull)
lres = {c: dict(coef=round(float(b[i + 1]), 4),
                z=round(float(b[i + 1] / se[i + 1]), 3),
                p=round(float(2 * (1 - stats.norm.cdf(abs(b[i + 1] / se[i + 1])))), 4))
        for i, c in enumerate(['accel', 'dvol', 'age_rank'])}
lres['pseudo_R2'] = round(pseudo_r2, 4)
n_sig = sum(1 for c, e_ in [('accel', '+'), ('dvol', '+'), ('age_rank', '-')]
            if lres[c]['p'] < 0.10 and ((lres[c]['coef'] > 0) == (e_ == '+')))
print(f'pseudo-R2={pseudo_r2:.4f}, significant-with-expected-sign: {n_sig}/3')

# ---------------- (b) forward-return OLS with HC1 robust SEs ----------------
def fit_ols_hc1(X, y):
    X = np.column_stack([np.ones(len(X)), X])
    XtX_inv = np.linalg.inv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    n, k = X.shape
    meat = (X * (e ** 2)[:, None]).T @ X
    cov = XtX_inv @ meat @ XtX_inv * (n / (n - k))  # HC1
    se = np.sqrt(np.diag(cov))
    r2 = 1 - (e @ e) / ((y - y.mean()) @ (y - y.mean()))
    return b, se, r2

ores = {}
for dep in ['ret24_raw', 'ret24_net']:
    dd2 = d[d[dep].notna()]
    Xa = dd2[['accel', 'dvol', 'age_rank']].values
    ya = dd2[dep].values
    b_, se_, r2_ = fit_ols_hc1(Xa, ya)
    names2 = ['const', 'accel', 'dvol', 'age_rank']
    print(f'\n--- {dep} OLS-HC1 (n={len(dd2)}) ---')
    coefs = {}
    for n_, bi, sei in zip(names2, b_, se_):
        t = bi / sei
        pv = 2 * (1 - stats.t.cdf(abs(t), len(ya) - 4))
        print(f'  {n_:10s} coef={bi:+.4f} se={sei:.4f} t={t:+.2f} p={pv:.4f}')
        if n_ != 'const':
            coefs[n_] = dict(coef=round(float(bi), 4), t=round(float(t), 3),
                             p=round(float(pv), 4))
    print(f'  R2={r2_:.4f}')
    ores[dep] = dict(n=len(dd2), coefs=coefs, r2=round(float(r2_), 4))

# ---------------- verdict ----------------
pass_a = (pseudo_r2 > 0.05) and (n_sig >= 2)
sig_b = {dep: sum(1 for c, e_ in [('accel', '+'), ('dvol', '+'), ('age_rank', '-')]
                  if ores[dep]['coefs'][c]['p'] < 0.10
                  and ((ores[dep]['coefs'][c]['coef'] > 0) == (e_ == '+')))
         for dep in ores}
pass_b = any(v >= 1 for v in sig_b.values())
verdict = 'PASS - mechanism warranted' if (pass_a or pass_b) else 'HONEST NEGATIVE - no mechanism'
print(f'\nVERDICT: {verdict}')

results = dict(
    n_episodes_total=int(len(ep)), n_crash_test=int(len(d)),
    n_dropped_forward=int(ep['crash'].isna().sum()),
    n_return_test=int(d['ret24_raw'].notna().sum()),
    crash_rate=round(float(d['crash'].mean()), 4),
    calibration_bins={lbl: dict(n=int(len(b)), crash_rate=round(float(b['crash'].mean()), 4))
                      for (lo, hi, lbl), b in
                      [((1.0, 1.5, '100-150pct'), d[(d['Rnet'] >= 1.0) & (d['Rnet'] < 1.5)]),
                       ((1.5, 99, 'gt150pct'), d[d['Rnet'] >= 1.5])] if len(b)},
    univariate=uni, logit=lres,
    ols=ores, n_sig_logit=n_sig, sig_b_return=sig_b,
    pass_a=bool(pass_a), pass_b=bool(pass_b), verdict=verdict)
with open(f'{BASE}/data/bubblechar_results.json', 'w') as f:
    json.dump(results, f, indent=2)
print('saved data/bubblechar_episodes.csv + data/bubblechar_results.json')
