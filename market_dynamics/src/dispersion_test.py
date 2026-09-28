"""Dispersion -> next-month market return test (pre-registered design).

UNIVERSE: true historical S&P 500 membership, panel_v2.csv (clean, post-unit-fix).
SAMPLE: dispersion months 2000-01..2021-11 (panel starts 1999-12-31; spx ends 2021-12).
  NOTE: task brief said 1996-2026; data supports 2000-2021. Design adapted honestly:
  in-sample 2000-01..2007-12 (96 obs), walk-forward 2008-01..2021-11 (167 forecasts).
MEASURES (month t):
  D_ew : equal-weighted cross-sectional std of monthly stock returns
  D_vw : cap-weighted cross-sectional std (weights = mcap at t-1, renormalized over available)
  V_gs : Goyal-Santa-Clara average variance = mean over stocks of sum of squared daily returns in month t
TARGET: R_{m,t+1}, S&P 500 simple monthly return from data/spx_monthly.csv (lr -> exp-1).
No trading mechanism built here -- measurement only.
"""
import numpy as np, pandas as pd, os, json, warnings
from scipy.optimize import minimize
warnings.filterwarnings('ignore')

BASE = os.path.expanduser('~/workspace/market_dynamics')
OUT = {}

# ---------------------------------------------------------------- panel
panel = pd.read_csv(BASE + '/data/panel_v2.csv',
                    usecols=['month', 'px_ticker', 'has_price', 'adjclose', 'mcap'])
panel['month'] = pd.to_datetime(panel['month'])
px   = panel.pivot_table(index='month', columns='px_ticker', values='adjclose')
mem  = panel.pivot_table(index='month', columns='px_ticker', values='has_price', aggfunc='max').fillna(0) > 0
mcap = panel.pivot_table(index='month', columns='px_ticker', values='mcap')
rets = px.pct_change(fill_method=None)
months = px.index

# sanity: clean panel?
mx = np.nanmax(mcap.values)
print('max mcap %.3e  (clean iff < ~4e12)' % mx)
assert mx < 4e12, 'panel looks corrupted'

# ------------------------------------------------- dispersion (a): monthly xs std
rec = []
for t in months[1:]:
    r = rets.loc[t]; m = mem.loc[t]
    rr = r[m & r.notna()]
    n = len(rr)
    ew = rr.std(ddof=1) if n > 10 else np.nan
    w = mcap.shift(1).loc[t]
    ww = w[m & r.notna() & w.notna()]
    if len(ww) > 10:
        wv = ww / ww.sum()
        mu = (wv * r[ww.index]).sum()
        vw = np.sqrt((wv * (r[ww.index] - mu) ** 2).sum() * len(ww) / (len(ww) - 1))
        cov = len(ww) / n
    else:
        vw, cov = np.nan, 0.0
    rec.append((t, ew, vw, n, cov))
disp = pd.DataFrame(rec, columns=['month', 'D_ew', 'D_vw', 'N', 'vw_cov']).set_index('month')
print('avg N/month %.0f, avg vw mcap coverage %.2f' % (disp['N'].mean(), disp['vw_cov'].mean()))

# ------------------------------------------- dispersion (b): Goyal-Santa-Clara
# NOTE: panel months are LAST TRADING DAYS (e.g. 2000-04-28), not calendar
# month-ends. Key daily aggregates by (year, month) period to avoid label mismatch.
rawdir = BASE + '/data/raw'
ssq_frames = []
for f in sorted(os.listdir(rawdir)):
    if not f.endswith('.csv'):
        continue
    tkr = f[:-4]
    d = pd.read_csv(os.path.join(rawdir, f), usecols=['date', 'adjclose'], parse_dates=['date'])
    d = d.sort_values('date').dropna(subset=['adjclose'])
    d['r'] = d['adjclose'].pct_change()
    d['ym'] = d['date'].dt.to_period('M')
    g = d.dropna(subset=['r']).groupby('ym')['r'].agg(['count', lambda s: float((s ** 2).sum())])
    g.columns = ['nday', 'ssq']
    g['ticker'] = tkr
    ssq_frames.append(g.reset_index())
ssq = pd.concat(ssq_frames, ignore_index=True)
ssq = ssq[ssq['nday'] >= 15]
print('ssq rows (>=15d):', len(ssq))

vrec = []
for t in months[1:]:
    ym = t.to_period('M')
    s = ssq[ssq['ym'] == ym]
    members_t = set(mem.columns[mem.loc[t]])
    s = s[s['ticker'].isin(members_t)]
    vrec.append((t, s['ssq'].mean() if len(s) > 10 else np.nan, len(s)))
vgs = pd.DataFrame(vrec, columns=['month', 'V_gs', 'N_gs']).set_index('month')
print('V_gs months:', vgs['V_gs'].notna().sum(), 'avg stocks %.0f' % vgs['N_gs'].mean())

# ---------------------------------------------------------------- market
# NOTE: spx ym are calendar months; panel months are last trading days.
# Map via (year, month) period to avoid label mismatch.
spx = pd.read_csv(BASE + '/data/spx_monthly.csv')
spx['ym_p'] = pd.to_datetime(spx['ym'] + '-01').dt.to_period('M')
per2label = {t.to_period('M'): t for t in months}
spx['month'] = spx['ym_p'].map(per2label)
spx['Rm'] = np.exp(spx['lr']) - 1
spx = spx.dropna(subset=['month', 'Rm']).set_index('month')[['Rm']].sort_index()
print('spx months mapped:', len(spx), spx.index.min().date(), '->', spx.index.max().date())

# ---------------------------------------------------------------- align
df = disp[['D_ew', 'D_vw']].join(vgs[['V_gs']]).join(spx)
df['Rm_next'] = df['Rm'].shift(-1)
df['down_next'] = (df['Rm_next'] < 0).astype(float)
df = df.dropna(subset=['D_ew', 'D_vw', 'V_gs', 'Rm_next'])
print('analysis months:', len(df), df.index.min().date(), '->', df.index.max().date())
df.to_csv(BASE + '/data/dispersion_monthly.csv')

IS = df.loc['2000-01':'2007-12']          # in-sample 96 obs
print('in-sample n =', len(IS))

# ---------------------------------------------------------------- helpers
def nw_t(y, x, L=3):
    X = np.column_stack([np.ones(len(x)), x])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ b
    n, k = X.shape
    Xe = X * e[:, None]
    S = Xe.T @ Xe
    for l in range(1, L + 1):
        wgt = 1 - l / (L + 1)
        G = Xe[l:].T @ Xe[:-l]
        S += wgt * (G + G.T)
    V = np.linalg.inv(X.T @ X) @ S @ np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(V))
    r2 = 1 - (e ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return b, se, b / se, r2

def logit_fit(y, x):
    def nll(b):
        z = b[0] + b[1] * x
        p = 1 / (1 + np.exp(-z))
        p = np.clip(p, 1e-9, 1 - 1e-9)
        return -(y * np.log(p) + (1 - y) * np.log(1 - p)).sum()
    r = minimize(nll, [0., 0.], method='BFGS')
    b = r.x
    # LR test vs intercept-only
    b0 = np.log(y.mean() / (1 - y.mean()))
    ll0 = (y * b0 - np.log1p(np.exp(b0))).sum()
    lr = 2 * (r.fun * -1 - ll0)
    from scipy.stats import chi2
    pval = float(chi2.sf(lr, 1))
    z = b[0] + b[1] * x; p = 1 / (1 + np.exp(-z))
    mcf = 1 - (-r.fun) / ll0
    return b, pval, mcf, p

def auc_manual(y, p):
    o = np.argsort(p); r = np.argsort(o) + 1
    n1, n0 = y.sum(), (1 - y).sum()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

res = {}
# ------------------------------------------------- in-sample OLS + walk-forward
for name in ['D_ew', 'D_vw', 'V_gs']:
    x, y = IS[name].values, IS['Rm_next'].values
    b, se, tstat, r2 = nw_t(y, x)
    # walk-forward expanding window
    fcs, acts = [], []
    idx = df.index
    oos = df.loc['2008-01':]
    for tt in oos.index:
        tr = df.loc[df.index < tt]          # pairs (D_s, R_{s+1}) fully observed at end of tt
        X = np.column_stack([np.ones(len(tr)), tr[name].values])
        bb, *_ = np.linalg.lstsq(X, tr['Rm_next'].values, rcond=None)
        fcs.append(bb[0] + bb[1] * df.loc[tt, name])
        acts.append(df.loc[tt, 'Rm_next'])
    fcs, acts = np.array(fcs), np.array(acts)
    # expanding historical-mean benchmark
    hmean = np.array([df.loc[df.index < tt]['Rm_next'].mean() for tt in oos.index])
    oos_r2 = 1 - ((acts - fcs) ** 2).sum() / ((acts - hmean) ** 2).sum()
    corr = np.corrcoef(fcs, acts)[0, 1]
    hit = np.mean(np.sign(fcs) == np.sign(acts))
    res[name] = dict(n_in=len(IS), slope=float(b[1]), se_nw=float(se[1]),
                     t_nw=float(tstat[1]), r2_in=float(r2),
                     n_oos=len(fcs), corr_oos=float(corr),
                     hit_rate=float(hit), oos_r2=float(oos_r2))
    print('%s: slope=%+.4f t_nw=%+.2f R2_in=%.4f | OOS corr=%+.3f hit=%.3f oosR2=%+.3f'
          % (name, b[1], tstat[1], r2, corr, hit, oos_r2))

# ------------------------------------------------- logit: down-month probability
    yl = IS['down_next'].values; xl = IS[name].values
    bl, pval, mcf, _ = logit_fit(yl, xl)
    # walk-forward predicted probs
    probs, yacts = [], []
    for tt in oos.index:
        tr = df.loc[df.index < tt]
        bb, _, _, _ = logit_fit(tr['down_next'].values, tr[name].values)
        z = bb[0] + bb[1] * df.loc[tt, name]
        probs.append(1 / (1 + np.exp(-z))); yacts.append(df.loc[tt, 'down_next'])
    probs, yacts = np.array(probs), np.array(yacts)
    auc = auc_manual(yacts, probs)
    hit05 = np.mean((probs > 0.5) == (yacts == 1))
    base = yacts.mean()
    res[name + '_logit'] = dict(coef=float(bl[1]), lr_pval=float(pval),
                                mcfadden_r2=float(mcf), auc_oos=float(auc),
                                hit05=float(hit05), base_rate=float(base))
    print('%s logit: coef=%+.3f LR-p=%.3f McFaddenR2=%.4f | OOS AUC=%.3f hit@.5=%.3f base=%.3f'
          % (name, bl[1], pval, mcf, auc, hit05, base))

json.dump(res, open(BASE + '/data/dispersion_test_results.json', 'w'), indent=1)
print('saved results + data/dispersion_monthly.csv')
