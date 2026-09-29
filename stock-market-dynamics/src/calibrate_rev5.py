"""Calibrate rev-5 params on 2000-2007 (pre-2008, no lookahead into test episodes).
Grid search: h1, lam1, p_damp on conditional hindcast (realized M),
minimizing full-distribution RMSE. T_calm/T_stress fit on same window.
Writes data/params_rev5.json with provenance labels.
"""
import os, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
from model4 import (N_BIN, assign_bins, cap_state, transition_matrix, entry_stats,
                    simulate, market_driver, rmse, PARAM_LABELS)

def monthly_vix():
    v = pd.read_csv(os.path.join(DATA, 'raw', 'VIX.csv'), parse_dates=['date'])
    v['ym'] = v['date'].dt.strftime('%Y-%m')
    return v.groupby('ym')['close'].mean().to_dict()

def conditional_hindcast(panel, months, Tcalm, Tstress, vix, M_hist, P, erate, emix):
    """Realized-M hindcast. Returns (pred_top, real_top, full_rmse)."""
    preds, reals = [], []
    c, _ = cap_state(panel[panel['month'] == months[0]])
    for k, m0 in enumerate(months[:-1]):
        m1 = months[k + 1]
        T = Tstress if vix.get(m0.strftime('%Y-%m'), 0) > 30.0 else Tcalm
        M = M_hist.get(m0.strftime('%Y-%m'), 0.0)
        traj = simulate(c, T, [M], P, 1, erate=erate, emix=emix,
                        p_damp=P['p_damp'], tau_build=P['tau_build'], gate=P['gate'])
        c = traj[1]
        preds.append(c)
        cr, _ = cap_state(panel[panel['month'] == m1])
        reals.append(cr)
        c = cr  # re-anchor monthly (conditional hindcast)
    preds = np.array(preds); reals = np.array(reals)
    return preds, reals, rmse(preds, reals)

def main():
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'), parse_dates=['month'], low_memory=False)
    panel['bin'] = assign_bins(panel)
    vix = monthly_vix()
    M_hist = market_driver(panel)
    d0, d1 = '2000-01-01', '2007-12-31'
    Tcalm, n_calm, m_calm = transition_matrix(panel, d0, d1, vix, regime='calm')
    Tstress, n_stress, m_stress = transition_matrix(panel, d0, d1, vix, regime='stress')
    print(f'T fit 2000-07: calm cap-transitions=${n_calm:.2e} ({m_calm}mo), '
          f'stress=${n_stress:.2e} ({m_stress}mo)', flush=True)
    erate, emix = entry_stats(panel, d0, d1)
    print(f'entry rate={erate:.4f}/mo', flush=True)
    months = sorted(panel[(panel['month'] >= '2001-01-01') & (panel['month'] <= '2007-12-31')]['month'].unique())
    print('hindcast months:', len(months), flush=True)

    best = None
    results = []
    for h1 in [0.0, 0.5, 1.0, 1.5, 2.0]:
        for lam1 in [0.0, 0.5, 1.0]:
            for p_damp in [0.8, 1.2, 1.6]:
                P = dict(kappa=1.0, h0=0.0, h1=h1, lam1=lam1,
                         p_damp=p_damp, tau_build=6.0, gate=0.5)
                preds, reals, e = conditional_hindcast(panel, months, Tcalm, Tstress, vix, M_hist, P, erate, emix)
                # also top-bin spike amplitude error (dot-com bust has a loser spike)
                bot_pred = preds[:, [b for b in range(N_BIN) if b // 8 == 0]].sum(axis=1)
                bot_real = reals[:, [b for b in range(N_BIN) if b // 8 == 0]].sum(axis=1)
                e_bot = rmse(bot_pred, bot_real)
                score = e + 0.5 * e_bot
                results.append((score, e, e_bot, h1, lam1, p_damp))
                print(f'h1={h1} lam1={lam1} p={p_damp}: rmse={e:.4f} bot={e_bot:.4f} score={score:.4f}', flush=True)
    results.sort()
    score, e, e_bot, h1, lam1, p_damp = results[0]
    print(f'BEST h1={h1} lam1={lam1} p_damp={p_damp} rmse={e:.4f}', flush=True)
    P = dict(kappa=1.0, h0=0.0, h1=h1, lam1=lam1, p_damp=p_damp,
             tau_build=6.0, gate=0.5, gamma=2.0, tcost=0.0005,
             _erate=erate, _emix=emix.tolist(),
             _Tcalm=Tcalm.tolist(), _Tstress=Tstress.tolist(),
             _cal_window='2000-2007', _labels=PARAM_LABELS,
             _note='T regime-split by VIX>30 at origin month; calibrated on conditional hindcast (realized M), full-distribution RMSE + 0.5*bottom-bin RMSE')
    json.dump(P, open(os.path.join(DATA, 'params_rev5.json'), 'w'), indent=1)
    # save T separately as npy for convenience
    np.save(os.path.join(DATA, 'Tcalm_rev5.npy'), Tcalm)
    np.save(os.path.join(DATA, 'Tstress_rev5.npy'), Tstress)
    print('wrote params_rev5.json', flush=True)

if __name__ == '__main__':
    main()
