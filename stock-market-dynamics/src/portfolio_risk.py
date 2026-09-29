"""Rev-5 risk metric 3: portfolio risk diagnostics.

(a) Predicted volatility of the rev-5 tilt portfolio from the bin-return
    covariance matrix Sigma (24x24, measured on 2001-2022 monthly cap-weighted
    bin total returns). Tilt = rev-4-style mass-flow tilt recomputed with
    params_rev5 at 2022-12, weights frozen (zero lookahead; M via mech_forecast).
    Reported vs SPY realized vol. This is a diagnostic, not a recommendation.
(b) Reversal probability: P(top momentum lane underperforms the cap-weighted
    market next month | C(t) in top quintile of its history) vs bottom quintile.
    Measured empirically on the fixed panel (2001-2022). Zero lookahead
    (quintiles from expanding history).
Saves outputs/portfolio_risk.json.
"""
import os, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'outputs')
from model4 import (N_BIN, assign_bins, cap_state, simulate, market_driver,
                    mech_forecast, top_share)
from paper_trade import tilt_weights, bin_returns


def spy_monthly():
    spy = pd.read_csv(os.path.join(DATA, 'raw', 'SPY.csv'), parse_dates=['date'])
    spy = spy.sort_values('date')
    sc = spy.set_index('date')['adjclose']
    months = pd.date_range('2001-01-31', '2022-12-31', freq='M')
    rets = []
    for m in months:
        c1 = sc[sc.index <= m]
        pm = m - pd.offsets.MonthEnd(1)
        c0 = sc[sc.index <= pm]
        if len(c1) and len(c0):
            rets.append((m, float(c1.iloc[-1] / c0.iloc[-1] - 1)))
    return pd.DataFrame(rets, columns=['month', 'r_spy'])


def main():
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'),
                        parse_dates=['month'], low_memory=False)
    panel['bin'] = assign_bins(panel)
    P = json.load(open(os.path.join(DATA, 'params_rev5.json')))
    months = sorted(panel['month'].unique())

    # ---- (a) predicted vol ----
    # monthly bin returns (cap-weighted, measured)
    Rdict = bin_returns(panel, months)  # (m0,m1) -> R[24]
    keys = sorted(Rdict.keys())
    R = pd.DataFrame([Rdict[k] for k in keys],
                     columns=[f'b{b}' for b in range(N_BIN)])
    R = R[(R != 0).any(axis=1)]
    mu = R.mean(axis=0).values
    Sigma = np.cov(R.values, rowvar=False)
    # tilt at 2022-12 with zero-lookahead M
    m_last = months[-1]
    c, _ = cap_state(panel[panel['month'] == m_last])
    M_hist = market_driver(panel)
    spy = spy_monthly()
    spy_ret = dict(zip(spy['month'].dt.strftime('%Y-%m'), spy['r_spy']))
    Mhat, _ = mech_forecast(spy_ret, m_last.strftime('%Y-%m'))
    T = np.array(P['_Tcalm'])
    traj = simulate(c, T, [Mhat], P, 1, erate=P['_erate'],
                    emix=np.array(P['_emix']), p_damp=P['p_damp'],
                    tau_build=P['tau_build'], gate=P['gate'])
    c_hat = traj[1]
    pi = tilt_weights(c, c_hat, gamma=P.get('gamma', 2.0))
    # predicted monthly vol and annualized
    port_var = float(pi @ Sigma @ pi)
    port_vol_m = float(np.sqrt(port_var))
    port_vol_a = port_vol_m * np.sqrt(12)
    # SPY realized vol for scale
    spy_vol_a = float(spy['r_spy'].std() * np.sqrt(12))
    # market (cap-weighted bin) predicted vol
    mkt_var = float(c @ Sigma @ c)
    print('predicted vol (annualized): tilt=%.1f%%  market=%.1f%%  SPY realized=%.1f%%' % (
        100 * port_vol_a, 100 * np.sqrt(mkt_var * 12), 100 * spy_vol_a), flush=True)
    # concentration of the tilt itself
    print('tilt: max bin weight=%.3f  n_bins>0=%d' % (pi.max(), int((pi > 1e-6).sum())), flush=True)

    # ---- (b) reversal probability ----
    # monthly top-lane excess return and C(t), expanding quintiles
    s = pd.read_csv(os.path.join(OUT, 'conc_ret_monthly.csv'), parse_dates=['month'])
    s = s.sort_values('month').reset_index(drop=True)
    # note: conc_ret_monthly uses SPY month-end timing (fixed); re_top - r_mkt for unwind
    s['unwind'] = (s['r_top'] < s['r_mkt']).astype(int)  # top lane trails market
    s['C_q_exp'] = np.nan
    for i in range(len(s)):
        hist = s['C'].iloc[:i + 1]
        s.loc[s.index[i], 'C_q_exp'] = float((hist <= s['C'].iloc[i]).mean())
    s['C_q'] = pd.qcut(s['C_q_exp'], 5, labels=False, duplicates='drop')
    rev = s.groupby('C_q')['unwind'].agg(['mean', 'count'])
    print('P(top lane unwinds next month | C quintile):', flush=True)
    for q in sorted(rev.index):
        print('  Q%d (crowding %s): %.1f%%  n=%d' % (
            q + 1, 'low->high'[0] if False else '', 100 * rev.loc[q, 'mean'],
            int(rev.loc[q, 'count'])), flush=True)
    p_high = float(rev.loc[4, 'mean']) if 4 in rev.index else float('nan')
    p_low = float(rev.loc[0, 'mean']) if 0 in rev.index else float('nan')

    out = {
        'tilt_predicted_vol_annual': port_vol_a,
        'market_predicted_vol_annual': float(np.sqrt(mkt_var * 12)),
        'spy_realized_vol_annual': spy_vol_a,
        'tilt_max_bin_weight': float(pi.max()),
        'reversal_prob_C_Q5': p_high,
        'reversal_prob_C_Q1': p_low,
        'reversal_n_Q5': int(rev.loc[4, 'count']) if 4 in rev.index else 0,
        'reversal_n_Q1': int(rev.loc[0, 'count']) if 0 in rev.index else 0,
        'note': ('Predicted vol from 24-bin return covariance 2001-2022 (measured). '
                 'Tilt = mass-flow tilt, params_rev5, frozen at 2022-12 (diagnostic). '
                 'Reversal = top lane trails market next month; quintiles expanding.'),
    }
    json.dump(out, open(os.path.join(OUT, 'portfolio_risk.json'), 'w'), indent=1)
    print('wrote outputs/portfolio_risk.json', flush=True)


if __name__ == '__main__':
    main()
