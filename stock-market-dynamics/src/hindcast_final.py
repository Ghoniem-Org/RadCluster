"""Rev-3 hindcasts with spike-persistence damper. Config (frozen):
  drift: stationary T (measured, per-episode window)
  advection: M_cs(t) cross-sectional window-delta (measured), kappa=1.0
             (mechanical), w_eff=0.15 (assumed)
  herding: h_eff = h1*(w-l), h1=1.0, h0=0 (calibrated 1-step pre-episode)
  stickiness: lam = min(1.0*max(w,l), 0.9) (lam1=1.0)
  spike persistence: damp=(|M|/build)^1.2 on reversal out of spike (gate 0.5),
             build memory tau=6mo (assumed)
  source/sink: measured
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model import (load_regimes, monthly_distributions, transition_matrix,
                   load_market_driver, simulate_full, entry_exit_rates,
                   rmse, N_BIN, bin_mom_q)

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', 'figs')
TOP_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 4]
BOT_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 0]
K, H0, H1, LAM1 = 1.0, 0.0, 1.0, 1.0


def run(d0, d1, t0, t1, focus, fig, ylabel, title):
    s, d = entry_exit_rates(reg, d0, d1)
    T = transition_matrix(reg, d0, d1)
    idx = [t for t in dist.index if pd.Timestamp(t0) <= t <= pd.Timestamp(t1)]
    actual = dist.loc[idx].values
    steps = [Mcs.get(t.strftime('%Y-%m'), 0.0) for t in idx[:-1]]
    traj = simulate_full(actual[0], T, steps, K, K, H0, H1, s, d,
                         len(idx) - 1, LAM1)
    e = rmse(traj, actual)
    fa = actual[:, focus].sum(1)
    fp = traj[:, focus].sum(1)
    ia, ip = int(np.argmax(fa)), int(np.argmax(fp))
    print(f"RMSE={e:.4f}", flush=True)
    print(f"peak actual {fa.max():.3f} @ {idx[ia].date()}, "
          f"predicted {fp.max():.3f} @ {idx[ip].date()}", flush=True)
    for k in [1, 2, 3]:
        if ia + k < len(idx) and ip + k < len(idx):
            print(f"  +{k}mo: actual {fa[ia+k]:.3f}, predicted {fp[ip+k]:.3f}",
                  flush=True)
    # full mom-marginal snapshot at peak for the report
    fig2, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(idx))
    ax.plot(x, fa, 'k-', lw=2, label='actual')
    ax.plot(x, fp, 'r--', lw=1.6,
            label='model: advection + auto-herding + sticky drift + spike persistence')
    ax.set_xticks(x[::4])
    ax.set_xticklabels([str(d_.date())[:7] for d_ in idx[::4]], rotation=30)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, fig), dpi=120)
    return e


if __name__ == '__main__':
    reg = load_regimes()
    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    Mcs = load_market_driver(source='panel')

    print("=== 2020-21 ===", flush=True)
    e1 = run('2000-12-01', '2019-12-31', '2020-01-01', '2021-12-31',
             TOP_BINS, 'hindcast_2020_final.png',
             'fraction in top momentum bin (>60%)',
             '2020-21 hindcast (rev 3): advection + autocatalytic herding + sticky drift + spike persistence')
    print("=== 2008-09 ===", flush=True)
    e2 = run('2000-12-01', '2007-12-31', '2008-01-01', '2009-12-31',
             BOT_BINS, 'hindcast_2008_final.png',
             'fraction in bottom momentum bin (<0%)',
             '2008-09 hindcast (rev 3): advection + autocatalytic herding + sticky drift + spike persistence')
    print(f"SUMMARY: 2020-21 RMSE={e1:.4f}, 2008-09 RMSE={e2:.4f}", flush=True)
