"""Calibrate (kappa, h0, h1) 1-step-ahead on pre-episode data, then free-run
hindcasts for 2020-21 and 2008-09 with the new kernels:
  - market advection: donor-cell shift along momentum driven by M(t)
  - autocatalytic herding: h_eff = h0 + h1*(winner_share - loser_share)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model import (load_regimes, monthly_distributions, transition_matrix,
                   load_market_driver, step_full, simulate_full,
                   entry_exit_rates, rmse, N_BIN, bin_mom_q)

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', 'figs')
TOP_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 4]
BOT_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 0]


def one_step_rmse(dist, Ms, T, s, d, k, h0, h1, d0, d1):
    """1-step-ahead RMSE over [d0, d1]: reset to actual each month."""
    idx = [t for t in dist.index if pd.Timestamp(d0) <= t <= pd.Timestamp(d1)]
    errs = []
    s_vec = np.full(N_BIN, s / N_BIN)
    for t0, t1 in zip(idx[:-1], idx[1:]):
        pred = step_full(dist.loc[t0].values, T, Ms.get(t0.strftime('%Y-%m'), 0.0),
                         k, k, h0, h1, s_vec, d)
        errs.append(rmse(pred, dist.loc[t1].values))
    return float(np.mean(errs))


def calibrate(dist, Ms, T, s, d, d0, d1, label):
    best = (1e9, None)
    for k in [0, 1, 2, 3, 4, 5]:
        for h0 in [0.0, 0.1]:
            for h1 in [0.0, 0.5, 1.0, 2.0]:
                e = one_step_rmse(dist, Ms, T, s, d, k, h0, h1, d0, d1)
                if e < best[0]:
                    best = (e, (k, h0, h1))
    print(f"[{label}] best 1-step: kappa={best[1][0]}, h0={best[1][1]}, "
          f"h1={best[1][2]}, RMSE={best[0]:.4f}", flush=True)
    return best[1]


def hindcast(dist, Ms, T, s, d, params, h0, h1, label, focus_bins, focus_name):
    k, hh0, hh1 = params
    idx = [t for t in dist.index if pd.Timestamp(h0) <= t <= pd.Timestamp(h1)]
    actual = dist.loc[idx].values
    steps = [Ms.get(t.strftime('%Y-%m'), 0.0) for t in idx[:-1]]
    traj = simulate_full(actual[0], T, steps, k, k, hh0, hh1, s, d, len(idx) - 1)
    e = rmse(traj, actual)
    pk_a = actual[:, focus_bins].sum(axis=1).max()
    pk_p = traj[:, focus_bins].sum(axis=1).max()
    t_a = idx[int(np.argmax(actual[:, focus_bins].sum(axis=1)))]
    t_p = idx[int(np.argmax(traj[:, focus_bins].sum(axis=1)))]
    print(f"[{label}] free-run RMSE={e:.4f} | peak actual {pk_a:.2f} @ {t_a.date()}, "
          f"predicted {pk_p:.2f} @ {t_p.date()}", flush=True)
    return idx, actual, traj, e


def main():
    reg = load_regimes()
    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    Ms = load_market_driver()

    # ---- Case 1: 2020-21, calibrate 2000-2019 ----
    s, d = entry_exit_rates(reg, '2000-12-01', '2019-12-31')
    T = transition_matrix(reg, '2000-12-01', '2019-12-31')
    p1 = calibrate(dist, Ms, T, s, d, '2000-12-01', '2019-12-31', 'cal 2000-19')
    i1, a1, t1, e1 = hindcast(dist, Ms, T, s, d, p1, '2020-01-01', '2021-12-31',
                              '2020-21', TOP_BINS, 'top')
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(i1))
    ax.plot(x, a1[:, TOP_BINS].sum(axis=1), 'k-', lw=2, label='actual')
    ax.plot(x, t1[:, TOP_BINS].sum(axis=1), 'r--', lw=1.6, label='model (new kernels)')
    ax.set_xticks(x[::4])
    ax.set_xticklabels([str(d.date())[:7] for d in i1[::4]], rotation=30)
    ax.set_ylabel('fraction in top momentum bin')
    ax.set_title('2020-21 hindcast with advection + autocatalytic herding')
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'hindcast_2020_newkernels.png'), dpi=120)

    # ---- Case 2: 2008-09, calibrate 2000-2007 ----
    s2, d2 = entry_exit_rates(reg, '2000-12-01', '2007-12-31')
    T2 = transition_matrix(reg, '2000-12-01', '2007-12-31')
    p2 = calibrate(dist, Ms, T2, s2, d2, '2000-12-01', '2007-12-31', 'cal 2000-07')
    i2, a2, t2, e2 = hindcast(dist, Ms, T2, s2, d2, p2, '2008-01-01', '2009-12-31',
                              '2008-09', BOT_BINS, 'bottom')
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(i2))
    ax.plot(x, a2[:, BOT_BINS].sum(axis=1), 'k-', lw=2, label='actual')
    ax.plot(x, t2[:, BOT_BINS].sum(axis=1), 'r--', lw=1.6, label='model (new kernels)')
    ax.set_xticks(x[::4])
    ax.set_xticklabels([str(d.date())[:7] for d in i2[::4]], rotation=30)
    ax.set_ylabel('fraction in bottom momentum bin (losers)')
    ax.set_title('2008-09 hindcast with advection + autocatalytic herding')
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'hindcast_2008_newkernels.png'), dpi=120)
    print("DONE", flush=True)


if __name__ == '__main__':
    main()
