"""V2: idiosyncratic drift (M-neutral months) + mechanical advection
f = kappa*M/w_eff + autocatalytic herding h_eff = h0 + h1*(w - l).
Calibrate 1-step-ahead, then free-run hindcasts.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model import (load_regimes, monthly_distributions, transition_matrix,
                   transition_matrix_mneutral, load_market_driver, step_full,
                   simulate_full, entry_exit_rates, rmse, N_BIN, bin_mom_q)

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', 'figs')
TOP_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 4]
BOT_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 0]


def one_step_rmse(dist, Ms, T, s, d, k, h0, h1, d0, d1):
    idx = [t for t in dist.index if pd.Timestamp(d0) <= t <= pd.Timestamp(d1)]
    errs = []
    s_vec = np.full(N_BIN, s / N_BIN)
    for t0, t1 in zip(idx[:-1], idx[1:]):
        pred = step_full(dist.loc[t0].values, T, Ms.get(t0.strftime('%Y-%m'), 0.0),
                         k, k, h0, h1, s_vec, d)
        errs.append(rmse(pred, dist.loc[t1].values))
    return float(np.mean(errs))


def calibrate(dist, Ms, T, s, d, d0, d1, label):
    results = []
    for k in [0.5, 1.0, 1.5, 2.0]:
        for h0 in [0.0, 0.1]:
            for h1 in [0.0, 0.5, 1.0, 2.0]:
                e = one_step_rmse(dist, Ms, T, s, d, k, h0, h1, d0, d1)
                results.append((e, k, h0, h1))
    results.sort()
    for e, k, h0, h1 in results[:5]:
        print(f"[{label}] kappa={k}, h0={h0}, h1={h1}: 1-step RMSE={e:.4f}", flush=True)
    return results[0][1:]


def hindcast(dist, Ms, T, s, d, params, h0, h1, label, focus_bins, fig_name,
             ylabel, title):
    k, hh0, hh1 = params
    idx = [t for t in dist.index if pd.Timestamp(h0) <= t <= pd.Timestamp(h1)]
    actual = dist.loc[idx].values
    steps = [Ms.get(t.strftime('%Y-%m'), 0.0) for t in idx[:-1]]
    traj = simulate_full(actual[0], T, steps, k, k, hh0, hh1, s, d, len(idx) - 1)
    e = rmse(traj, actual)
    fa = actual[:, focus_bins].sum(axis=1)
    fp = traj[:, focus_bins].sum(axis=1)
    t_a = idx[int(np.argmax(fa))]; t_p = idx[int(np.argmax(fp))]
    print(f"[{label}] kappa={k},h0={hh0},h1={hh1}: free-run RMSE={e:.4f} | "
          f"peak actual {fa.max():.2f} @ {t_a.date()}, predicted {fp.max():.2f} @ {t_p.date()}",
          flush=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(idx))
    ax.plot(x, fa, 'k-', lw=2, label='actual')
    ax.plot(x, fp, 'r--', lw=1.6, label='model (advection + auto-herding)')
    ax.set_xticks(x[::4])
    ax.set_xticklabels([str(d.date())[:7] for d in idx[::4]], rotation=30)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, fig_name), dpi=120)
    return e


def main():
    reg = load_regimes()
    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    Ms = load_market_driver()

    # ---- Case 1: 2020-21 ----
    s, d = entry_exit_rates(reg, '2000-12-01', '2019-12-31')
    T_full = transition_matrix(reg, '2000-12-01', '2019-12-31')
    T_idio, ntr, nmo = transition_matrix_mneutral(reg, Ms, '2000-12-01', '2019-12-31')
    print(f"T_idio 2000-19: {ntr} transitions, {nmo} neutral months", flush=True)
    print("== calibrate T_full ==", flush=True)
    p_full = calibrate(dist, Ms, T_full, s, d, '2000-12-01', '2019-12-31', 'T_full')
    print("== calibrate T_idio ==", flush=True)
    p_idio = calibrate(dist, Ms, T_idio, s, d, '2000-12-01', '2019-12-31', 'T_idio')
    e_full = hindcast(dist, Ms, T_full, s, d, p_full, '2020-01-01', '2021-12-31',
                      '2020-21 T_full', TOP_BINS, 'hindcast_2020_v2_full.png',
                      'fraction in top momentum bin',
                      '2020-21 hindcast: full drift + advection + auto-herding')
    e_idio = hindcast(dist, Ms, T_idio, s, d, p_idio, '2020-01-01', '2021-12-31',
                      '2020-21 T_idio', TOP_BINS, 'hindcast_2020_v2_idio.png',
                      'fraction in top momentum bin',
                      '2020-21 hindcast: idio drift + advection + auto-herding')

    # ---- Case 2: 2008-09 ----
    s2, d2 = entry_exit_rates(reg, '2000-12-01', '2007-12-31')
    T2_full = transition_matrix(reg, '2000-12-01', '2007-12-31')
    T2_idio, ntr2, nmo2 = transition_matrix_mneutral(reg, Ms, '2000-12-01', '2007-12-31')
    print(f"T_idio 2000-07: {ntr2} transitions, {nmo2} neutral months", flush=True)
    print("== calibrate T2_full ==", flush=True)
    p2_full = calibrate(dist, Ms, T2_full, s2, d2, '2000-12-01', '2007-12-31', 'T2_full')
    print("== calibrate T2_idio ==", flush=True)
    p2_idio = calibrate(dist, Ms, T2_idio, s2, d2, '2000-12-01', '2007-12-31', 'T2_idio')
    hindcast(dist, Ms, T2_full, s2, d2, p2_full, '2008-01-01', '2009-12-31',
             '2008-09 T2_full', BOT_BINS, 'hindcast_2008_v2_full.png',
             'fraction in bottom momentum bin (losers)',
             '2008-09 hindcast: full drift + advection + auto-herding')
    hindcast(dist, Ms, T2_idio, s2, d2, p2_idio, '2008-01-01', '2009-12-31',
             '2008-09 T2_idio', BOT_BINS, 'hindcast_2008_v2_idio.png',
             'fraction in bottom momentum bin (losers)',
             '2008-09 hindcast: idio drift + advection + auto-herding')
    print("DONE", flush=True)


if __name__ == '__main__':
    main()
