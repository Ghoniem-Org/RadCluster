"""Conditional hindcasts (diagnostic, realized M) for 2008-09 and 2020-21.
Free-running from pre-episode state; M(t) realized (exogenous forcing).
NOT a forecast -- the paper-trading walk-forward is the forecast track record.
"""
import os, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
from model4 import (N_BIN, assign_bins, cap_state, simulate, market_driver,
                    top_share, bot_share, rmse)

def run_episode(panel, P, start_ym, end_ym, label):
    M_hist = market_driver(panel)
    vix = pd.read_csv(os.path.join(DATA, 'raw', 'VIX.csv'), parse_dates=['date'])
    vix['ym'] = v['date'].dt.strftime('%Y-%m') if False else vix['date'].dt.strftime('%Y-%m')
    vx = vix.groupby('ym')['close'].mean().to_dict()
    Tcalm = np.array(P['_Tcalm']); Tstress = np.array(P['_Tstress'])
    months = sorted(panel[(panel['month'].dt.strftime('%Y-%m') >= start_ym)
                          & (panel['month'].dt.strftime('%Y-%m') <= end_ym)]['month'].unique())
    c, _ = cap_state(panel[panel['month'] == months[0]])
    Ms, Ts = [], []
    for m in months[:-1]:
        ym = m.strftime('%Y-%m')
        Ms.append(M_hist.get(ym, 0.0))
        Ts.append(Tstress if vx.get(ym, 0) > 30.0 else Tcalm)
    # free-running with per-month regime T: simulate manually
    traj = [c.copy()]
    cc = c.copy()
    bdec = np.exp(-1.0 / P['tau_build'])
    build_lo, build_hi = 0.15, 0.15
    M_prev = Ms[0]
    from model4 import step_full
    for t in range(len(Ms)):
        M = Ms[t]; F = abs(M)
        w, l = top_share(cc), bot_share(cc)
        build_lo = 0.15 + (build_lo - 0.15) * bdec
        build_hi = 0.15 + (build_hi - 0.15) * bdec
        if M < 0: build_lo = max(build_lo, F)
        elif M > 0: build_hi = max(build_hi, F)
        dlo, dhi = 1.0, 1.0
        if M > 0 and M_prev < 0 and l > P['gate']: dlo = min((F / build_lo) ** P['p_damp'], 1.0)
        elif M < 0 and M_prev > 0 and w > P['gate']: dhi = min((F / build_hi) ** P['p_damp'], 1.0)
        cc = step_full(cc, Ts[t], M, P, damp_lo=dlo, damp_hi=dhi,
                       erate=P['_erate'], emix=np.array(P['_emix']))
        M_prev = M
        traj.append(cc.copy())
    traj = np.array(traj)
    real = []
    for m in months:
        cr, _ = cap_state(panel[panel['month'] == m])
        real.append(cr)
    real = np.array(real)
    top_p = traj[:, [b for b in range(N_BIN) if b // 8 == 2]].sum(axis=1)
    top_r = real[:, [b for b in range(N_BIN) if b // 8 == 2]].sum(axis=1)
    bot_p = traj[:, [b for b in range(N_BIN) if b // 8 == 0]].sum(axis=1)
    bot_r = real[:, [b for b in range(N_BIN) if b // 8 == 0]].sum(axis=1)
    out = dict(label=label, months=[m.strftime('%Y-%m') for m in months],
               top_pred=top_p.tolist(), top_real=top_r.tolist(),
               bot_pred=bot_p.tolist(), bot_real=bot_r.tolist(),
               full_rmse=rmse(traj, real),
               top_rmse=rmse(top_p, top_r), bot_rmse=rmse(bot_p, bot_r))
    # peak stats
    for nm, pr, rr in [('top', top_p, top_r), ('bot', bot_p, bot_r)]:
        ip, ir = int(np.argmax(pr)), int(np.argmax(rr))
        out[f'{nm}_peak'] = dict(pred=float(pr[ip]), pred_ym=out['months'][ip],
                                 real=float(rr[ir]), real_ym=out['months'][ir])
    json.dump(out, open(os.path.join(DATA, f'hindcast_rev5_{label}.json'), 'w'), indent=1)
    print(label, 'full_rmse=%.4f top_peak pred=%.2f@%s real=%.2f@%s' % (
        out['full_rmse'], out['top_peak']['pred'], out['top_peak']['pred_ym'],
        out['top_peak']['real'], out['top_peak']['real_ym']), flush=True)
    print(label, 'bot_peak pred=%.2f@%s real=%.2f@%s' % (
        out['bot_peak']['pred'], out['bot_peak']['pred_ym'],
        out['bot_peak']['real'], out['bot_peak']['real_ym']), flush=True)
    return out

if __name__ == '__main__':
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'), parse_dates=['month'], low_memory=False)
    panel['bin'] = assign_bins(panel)
    P = json.load(open(os.path.join(DATA, 'params_rev5.json')))
    run_episode(panel, P, '2007-12', '2009-12', '2008-09')
    run_episode(panel, P, '2019-12', '2021-12', '2020-21')
