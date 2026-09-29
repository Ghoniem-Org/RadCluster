"""PAPER-TRADING deployment protocol (rev-4 addendum).

At each month-end t (forecast origin):
  1. Observe measured state c(t) (cap-weighted 24-bin, PIT data only).
  2. Forecast M(t+1) via expanding AR(1) on M <= t. Zero lookahead.
  3. One model step -> predicted bin shares chat(t+1).
  4. Tilt weights pi_b(t) propto c_b(t)*exp(gamma*dchat_b/max(c_b,eps)).
     Within bin: constituent cap weights at t. WEIGHTS FROZEN at t.
  5. Hold to t+1 at realized prices. No intra-month trading.
  6. At t+1: rebalance to new frozen weights; pay 5bps one-way on turnover.
Benchmarks with identical accounting: SPY buy-and-hold (total return),
equal-weight bins (1/24, monthly rebalance, same costs).
Monthly snapshots -> goal workspace (paper_track_record.csv).
This is the prerequisite track record before any real capital.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
GOAL = os.path.expanduser('~/workspace/goals/stock-market-dynamics')
SNAP = os.path.join(GOAL, 'hidden_files', 'paper_track_record.csv')

from model4 import (N_BIN, assign_bins, cap_state, transition_matrix, entry_stats,
                    simulate, market_driver, mech_forecast, top_share, bot_share,
                    rmse, PARAM_LABELS)

EPS = 1e-6

def tilt_weights(c, c_hat, gamma=1.0):
    """Additive tilt: overweight bins predicted to gain mass, underweight
    predicted losers, in proportion to predicted flow. Stable when bins are
    small (no exponential blowup). pi_b propto max(c_b + gamma*dchat_b, 0)."""
    d = c_hat - c
    pi = np.maximum(c + gamma * d, 0.0)
    s = pi.sum()
    return pi / s if s > 0 else c

def bin_returns(panel, months):
    """Cap-weighted total bin returns between consecutive months.
    Returns dict (m0, m1) -> R[24]."""
    d = panel[(panel['bin'] >= 0) & (panel['mcap'] > 0)].copy()
    d = d.sort_values(['ticker', 'month'])
    d['adj_next'] = d.groupby('ticker')['adjclose'].shift(-1)
    d['month_next'] = d.groupby('ticker')['month'].shift(-1)
    d['dm'] = (d['month_next'].dt.year * 12 + d['month_next'].dt.month
               - (d['month'].dt.year * 12 + d['month'].dt.month))
    d = d[(d['dm'] == 1) & d['adj_next'].notna() & (d['adjclose'] > 0)]
    d['ret'] = d['adj_next'] / d['adjclose'] - 1
    out = {}
    for (m0, m1), g in d.groupby(['month', 'month_next']):
        R = np.zeros(N_BIN)
        tot = g['mcap'].sum()
        for b, gg in g.groupby('bin'):
            wb = gg['mcap'].values
            wb = wb / wb.sum()
            R[int(b)] = float(np.dot(wb, gg['ret'].values))
        out[(m0.strftime('%Y-%m-%d'), m1.strftime('%Y-%m-%d'))] = R
    return out

def stock_returns(panel, m0, m1):
    """Per-ticker realized returns m0->m1 for stocks present at m0 with prices."""
    a = panel[(panel['month'] == m0) & (panel['has_price'] == 1)
              & (panel['bin'] >= 0) & (panel['mcap'] > 0)]
    b = panel[(panel['month'] == m1) & (panel['has_price'] == 1)][['ticker', 'adjclose']]
    m = a.merge(b, on='ticker', suffixes=('', '_n'))
    m = m[m['adjclose_n'].notna() & (m['adjclose'] > 0)]
    m['ret'] = m['adjclose_n'] / m['adjclose'] - 1
    return m[['ticker', 'bin', 'mcap', 'ret']]

def run(panel, P, months, gamma=1.0, tcost=0.0005, vix=None):
    """Walk-forward paper trade. months: list of month-end Timestamps (origins).
    Returns snapshots DataFrame."""
    M_hist = market_driver(panel)
    # regime T's fit on pre-start data are passed in via P['_Tcalm'], P['_Tstress']
    Tcalm, Tstress = np.array(P['_Tcalm']), np.array(P['_Tstress'])
    erate, emix = P['_erate'], np.array(P['_emix'])
    # VIX regime per origin month
    regime_stress = {}
    if vix is not None:
        for m in months:
            regime_stress[m] = vix.get(m.strftime('%Y-%m'), 0) > 30.0

    # SPY total-return series + monthly returns for mechanical M nowcast
    spy = pd.read_csv(os.path.join(DATA, 'raw', 'SPY.csv'), parse_dates=['date'])
    spy = spy.sort_values('date').set_index('date')['adjclose']
    spy_ret = {}
    for m in months:
        d = spy[spy.index <= m]
        if len(d) >= 22:
            spy_ret[m.strftime('%Y-%m')] = float(d.iloc[-1] / d.iloc[-22] - 1)

    snaps = []
    val_paper, val_spy, val_ew = 1.0, 1.0, 1.0
    w_prev = None  # paper: drifted weights before rebalance (ticker, w)
    w_prev_ew = None
    for k, m0 in enumerate(months[:-1]):
        m1 = months[k + 1]
        ym0 = m0.strftime('%Y-%m')
        ym1 = m1.strftime('%Y-%m')
        # 1. measured state
        d0 = panel[panel['month'] == m0]
        c, cov = cap_state(d0)
        # 2. mechanical M nowcast (zero lookahead)
        Mhat, ar = mech_forecast(spy_ret, ym0)
        Mreal = M_hist.get(ym1, np.nan)
        # 3. one-step prediction with regime T
        T = Tstress if regime_stress.get(m0, False) else Tcalm
        traj = simulate(c, T, [Mhat], P, 1, erate=erate, emix=emix,
                        p_damp=P['p_damp'], tau_build=P['tau_build'], gate=P['gate'])
        c_hat = traj[1]
        # 4. frozen tilt weights (bin level) -> stock level
        pi = tilt_weights(c, c_hat, gamma)
        sret = stock_returns(panel, m0.strftime('%Y-%m-%d'), m1.strftime('%Y-%m-%d'))
        bin_mcap = sret.groupby('bin')['mcap'].sum()
        sret['w'] = sret.apply(lambda r: pi[int(r['bin'])] * r['mcap'] / bin_mcap[int(r['bin'])], axis=1)
        # 5. realized portfolio return (frozen weights)
        rp = float((sret['w'] * sret['ret']).sum())
        # rebalance cost: turnover vs drifted previous weights
        cost = 0.0
        if w_prev is not None:
            # drift w_prev by realized returns to m1 then compare
            drift = w_prev.merge(sret[['ticker', 'ret']], on='ticker', how='left')
            drift['ret'] = drift['ret'].fillna(0.0)
            drift['w_d'] = drift['w'] * (1 + drift['ret'])
            drift['w_d'] /= drift['w_d'].sum()
            new = sret[['ticker', 'w']].merge(drift[['ticker', 'w_d']], on='ticker', how='outer').fillna(0)
            turnover = float(np.abs(new['w'] - new['w_d']).sum()) / 2
            cost = tcost * turnover
        val_paper *= (1 + rp - cost)
        # EW benchmark: 1/24 bins, same mechanics and cost accounting
        pi_ew = np.full(N_BIN, 1.0 / N_BIN)
        sret['w_ew'] = sret.apply(lambda r: pi_ew[int(r['bin'])] * r['mcap'] / bin_mcap[int(r['bin'])], axis=1)
        r_ew = float((sret['w_ew'] * sret['ret']).sum())
        cost_ew = 0.0
        if w_prev_ew is not None:
            drift = w_prev_ew.merge(sret[['ticker', 'ret']], on='ticker', how='left')
            drift['ret'] = drift['ret'].fillna(0.0)
            drift['w_d'] = drift['w_ew'] * (1 + drift['ret'])
            drift['w_d'] /= drift['w_d'].sum()
            new = sret[['ticker', 'w_ew']].merge(drift[['ticker', 'w_d']], on='ticker', how='outer').fillna(0)
            turnover_ew = float(np.abs(new['w_ew'] - new['w_d']).sum()) / 2
            cost_ew = tcost * turnover_ew
        val_ew *= (1 + r_ew - cost_ew)
        # SPY buy-and-hold
        sp0 = spy[spy.index <= m0].iloc[-1] if (spy.index <= m0).any() else np.nan
        sp1 = spy[spy.index <= m1].iloc[-1] if (spy.index <= m1).any() else np.nan
        r_spy = float(sp1 / sp0 - 1) if np.isfinite(sp0) and np.isfinite(sp1) else 0.0
        val_spy *= (1 + r_spy)
        w_prev = sret[['ticker', 'w']]
        w_prev_ew = sret[['ticker', 'w_ew']]
        d1 = panel[panel['month'] == m1]
        c1, _ = cap_state(d1)
        snaps.append(dict(date=m1.strftime('%Y-%m-%d'), origin=ym0,
                          M_hat=round(Mhat, 4), M_real=round(Mreal, 4) if np.isfinite(Mreal) else None,
                          M_method='mech',
                          paper=round(val_paper, 4), spy=round(val_spy, 4), ew=round(val_ew, 4),
                          r_paper=round(rp, 4), cost=round(cost, 5),
                          top_pred=round(float(c_hat[[b for b in range(N_BIN) if b // 8 == 2]].sum()), 4),
                          top_real=round(float(c1[[b for b in range(N_BIN) if b // 8 == 2]].sum()), 4),
                          n_hold=int(len(sret))))
    return pd.DataFrame(snaps)

def summarize(snaps):
    s = snaps.copy()
    out = {}
    for col in ['paper', 'spy', 'ew']:
        v = s[col].values
        rets = v[1:] / v[:-1] - 1
        n = len(rets) / 12
        out[col] = dict(
            total_return=round(v[-1] / v[0] - 1, 4),
            cagr=round((v[-1] / v[0]) ** (1 / n) - 1, 4) if n > 0 else None,
            ann_vol=round(float(np.std(rets, ddof=1) * np.sqrt(12)), 4),
            sharpe=round(float(np.mean(rets) / np.std(rets, ddof=1) * np.sqrt(12)), 3) if np.std(rets) > 0 else None,
            max_dd=round(float(np.min(v / np.maximum.accumulate(v) - 1)), 4),
        )
    return out

if __name__ == '__main__':
    import json
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'), parse_dates=['month'], low_memory=False)
    panel['bin'] = assign_bins(panel)
    P = json.load(open(os.path.join(DATA, 'params_rev4.json')))
    months = sorted(panel[panel['month'] >= '2008-01-01']['month'].unique())
    print('origins:', len(months), flush=True)
    snaps = run(panel, P, months)
    os.makedirs(os.path.dirname(SNAP), exist_ok=True)
    snaps.to_csv(SNAP, index=False)
    print(snaps.tail(3).to_string(), flush=True)
    print(json.dumps(summarize(snaps), indent=1), flush=True)
