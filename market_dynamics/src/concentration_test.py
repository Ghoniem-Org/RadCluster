"""REV-5 TASK 1: Does concentration predict subsequent returns? MEASURE FIRST.

Pre-registered design (written before estimation; mirrored in REPORT.md rev-5):
  CONCENTRATION C(t): top-momentum-lane (m=2) capital share at month-end t. Measured.
  DEPENDENT r_top^e(t+1): cap-weighted total return of top-lane constituents
      t->t+1 (weights = cap weights at t) minus SPY total return t->t+1. Measured.
  ESTIMATION: origins t in [2001-01-31, 2007-12-31] (first month with >100
      binnable rows through last pre-2008 origin). Walk-forward: 2008-01..2022-11.
  FORM 1 (threshold-linear, pre-specified): r^e = a + b*max(0, C-0.5) + e.
      C*=0.5 fixed (= damper gate; not fitted). Pass bars: b<0, |t|>2.0,
      walk-forward corr(pred, realized) > 0.10 with frozen (a,b).
  FORM 2 (quintile non-parametric, pre-specified): quintile cutoffs of C from
      the estimation window only; mean r^e by quintile. Pass bars:
      mean(Q5)-mean(Q1) <= -0.25 %/mo in-sample AND same sign walk-forward.
  If neither passes -> HONEST NEGATIVE, no return mechanism is built.

Secondary/descriptive only (not pass/fail): HHI version, bottom-lane bounce,
full-sample correlations.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'outputs')

import sys
sys.path.insert(0, HERE)
from model4 import N_BIN, assign_bins, cap_state, mom_of

TOP_BINS = [b for b in range(N_BIN) if mom_of(b) == 2]
BOT_BINS = [b for b in range(N_BIN) if mom_of(b) == 0]


def build_monthly(panel):
    """Monthly series: C(t), HHI(t), bot share, next-month lane returns, SPY."""
    months = sorted(panel['month'].unique())
    # per-ticker next-month returns
    d = panel[(panel['bin'] >= 0) & (panel['mcap'] > 0)].copy()
    d = d.sort_values(['ticker', 'month'])
    d['an'] = d.groupby('ticker')['adjclose'].shift(-1)
    d['mn'] = d.groupby('ticker')['month'].shift(-1)
    d['dm'] = (d['mn'].dt.year * 12 + d['mn'].dt.month
               - (d['month'].dt.year * 12 + d['month'].dt.month))
    d = d[(d['dm'] == 1) & d['an'].notna() & (d['adjclose'] > 0)].copy()
    d['ret'] = d['an'] / d['adjclose'] - 1
    # SPY monthly returns, month-end to month-end to match stock timing
    # (stock ret = adjclose[next month-end]/adjclose[this month-end] - 1)
    spy = pd.read_csv(os.path.join(DATA, 'raw', 'SPY.csv'), parse_dates=['date'])
    spy = spy.sort_values('date')
    spy_close = spy.set_index('date')['adjclose']

    def spy_ret(m0, m1):
        c0 = spy_close[spy_close.index <= m0]
        c1 = spy_close[spy_close.index <= m1]
        if len(c0) == 0 or len(c1) == 0:
            return np.nan
        return float(c1.iloc[-1] / c0.iloc[-1] - 1)
    rows = []
    for m0 in months[:-1]:
        m1 = months[months.index(m0) + 1]
        d0 = panel[panel['month'] == m0]
        if (d0['bin'] >= 0).sum() < 100:
            continue
        c, _ = cap_state(d0)
        C = float(c[TOP_BINS].sum())
        hhi = float((c ** 2).sum())
        bot = float(c[BOT_BINS].sum())
        g = d[d['month'] == m0]
        if len(g) == 0:
            continue
        tot = g['mcap'].sum()
        # lane returns: cap-weighted over lane members at t
        def lane_ret(bins):
            gg = g[g['bin'].isin(bins)]
            if len(gg) == 0 or gg['mcap'].sum() == 0:
                return np.nan
            w = gg['mcap'].values / gg['mcap'].sum()
            return float(np.dot(w, gg['ret'].values))
        r_top = lane_ret(TOP_BINS)
        r_bot = lane_ret(BOT_BINS)
        r_mkt = float(np.dot(g['mcap'].values / tot, g['ret'].values))
        r_spy = spy_ret(m0, m1)
        rows.append(dict(month=m0, C=C, hhi=hhi, bot=bot,
                         r_top=r_top, r_bot=r_bot, r_mkt=r_mkt, r_spy=r_spy,
                         n_bin=int((d0['bin'] >= 0).sum())))
    s = pd.DataFrame(rows)
    s['re_top'] = s['r_top'] - s['r_spy']   # excess over SPY
    s['re_bot'] = s['r_bot'] - s['r_spy']
    s['sprd'] = s['r_top'] - s['r_bot']     # winner-minus-loser
    return s


def ols_tstats(X, y):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    XtX = X.T @ X
    beta = np.linalg.solve(XtX, X.T @ y)
    resid = y - X @ beta
    n, k = X.shape
    s2 = resid @ resid / (n - k)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(XtX)))
    return beta, se, beta / se, 1 - resid @ resid / ((y - y.mean()) @ (y - y.mean()))


def main():
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'),
                        parse_dates=['month'], low_memory=False)
    panel['bin'] = assign_bins(panel)
    s = build_monthly(panel)
    s.to_csv(os.path.join(OUT, 'conc_ret_monthly.csv'), index=False)
    print('months:', len(s), s['month'].min(), '->', s['month'].max(), flush=True)

    est = s[(s['month'] >= '2001-01-01') & (s['month'] < '2008-01-01')].copy()
    wf = s[s['month'] >= '2008-01-01'].copy()
    print(f'estimation n={len(est)}, walk-forward n={len(wf)}', flush=True)
    print('C(t) estimation: min %.3f  p50 %.3f  p90 %.3f  max %.3f  #>0.5: %d' % (
        est['C'].min(), est['C'].median(),
        est['C'].quantile(0.9), est['C'].max(), int((est['C'] > 0.5).sum())), flush=True)

    y = est['re_top'].values
    pen = np.maximum(est['C'].values - 0.5, 0.0)
    X = np.column_stack([np.ones_like(y), pen])
    beta, se, tstat, r2 = ols_tstats(X, y)
    print('FORM1 in-sample: a=%.4f (t=%.2f)  b=%.4f (t=%.2f)  R2=%.3f  n_pen=%d' % (
        beta[0], tstat[0], beta[1], tstat[1], r2, int((pen > 0).sum())), flush=True)
    # walk-forward with frozen params
    pen_wf = np.maximum(wf['C'].values - 0.5, 0.0)
    pred_wf = beta[0] + beta[1] * pen_wf
    rw = wf['re_top'].values
    okm = np.isfinite(pred_wf) & np.isfinite(rw)
    corr_wf = float(np.corrcoef(pred_wf[okm], rw[okm])[0, 1]) if okm.sum() > 5 else np.nan
    print('FORM1 walk-forward corr(pred,realized) = %.3f  (n=%d)' % (corr_wf, okm.sum()), flush=True)
    f1_pass = (beta[1] < 0) and (abs(tstat[1]) > 2.0) and (corr_wf > 0.10)
    print('FORM1 PASS' if f1_pass else 'FORM1 FAIL', flush=True)

    # FORM 2: quintiles from estimation window
    qs = est['C'].quantile([0.2, 0.4, 0.6, 0.8]).values
    est['q'] = np.clip(np.digitize(est['C'].values, qs), 0, 4)
    wf['q'] = np.clip(np.digitize(wf['C'].values, qs), 0, 4)
    m_est = est.groupby('q')['re_top'].agg(['mean', 'count'])
    m_wf = wf.groupby('q')['re_top'].agg(['mean', 'count'])
    print('FORM2 in-sample quintile means (%/mo) [n]:', flush=True)
    for qi in sorted(m_est.index):
        print('   Q%d: %+.3f%%  n=%d' % (qi + 1, 100 * m_est.loc[qi, 'mean'],
                                        int(m_est.loc[qi, 'count'])), flush=True)
    spread_is = float(m_est.loc[4, 'mean'] - m_est.loc[0, 'mean'])
    spread_wf = float(m_wf.loc[4, 'mean'] - m_wf.loc[0, 'mean'])
    print('FORM2 spread Q5-Q1: in-sample %.4f%%/mo, walk-forward %.4f%%/mo' % (
        100 * spread_is, 100 * spread_wf), flush=True)
    print('FORM2 walk-forward quintile means (%/mo) [n]:', flush=True)
    for qi in sorted(m_wf.index):
        print('   Q%d: %+.3f%%  n=%d' % (qi + 1, 100 * m_wf.loc[qi, 'mean'],
                                        int(m_wf.loc[qi, 'count'])), flush=True)
    f2_pass = (spread_is <= -0.0025) and (spread_wf < 0)
    print('FORM2 PASS' if f2_pass else 'FORM2 FAIL', flush=True)

    # descriptive secondaries
    print('--- descriptive ---', flush=True)
    print('corr(C, re_top) est=%.3f wf=%.3f' % (
        est['C'].corr(est['re_top']), wf['C'].corr(wf['re_top'])), flush=True)
    print('corr(hhi, re_top) est=%.3f wf=%.3f' % (
        est['hhi'].corr(est['re_top']), wf['hhi'].corr(wf['re_top'])), flush=True)
    print('corr(bot_share, re_bot) est=%.3f wf=%.3f' % (
        est['bot'].corr(est['re_bot']), wf['bot'].corr(wf['re_bot'])), flush=True)
    print('corr(C, sprd next-mo) est=%.3f wf=%.3f' % (
        est['C'].corr(est['sprd']), wf['C'].corr(wf['sprd'])), flush=True)
    print('VERDICT:', 'BUILD return mechanism' if (f1_pass or f2_pass)
          else 'HONEST NEGATIVE - stop', flush=True)


if __name__ == '__main__':
    main()
