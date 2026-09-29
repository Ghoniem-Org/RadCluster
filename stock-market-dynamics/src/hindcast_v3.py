"""V3: cross-sectional advection driver M_cs + 2-stage calibration.
Stage 1: kappa on high-|M| months 1-step (the kernel's job).
Stage 2: (h0, h1) on all months 1-step.
Then free-run hindcasts.
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


def one_step(dist, Ms, T, s, d, k, h0, h1, d0, d1, m_min=0.0):
    idx = [t for t in dist.index if pd.Timestamp(d0) <= t <= pd.Timestamp(d1)]
    errs = []
    s_vec = np.full(N_BIN, s / N_BIN)
    for t0, t1 in zip(idx[:-1], idx[1:]):
        M = Ms.get(t0.strftime('%Y-%m'), 0.0)
        if abs(M) < m_min:
            continue
        pred = step_full(dist.loc[t0].values, T, M, k, k, h0, h1, s_vec, d)
        errs.append(rmse(pred, dist.loc[t1].values))
    return float(np.mean(errs)) if errs else float('nan')


def calibrate(dist, Ms, T, s, d, d0, d1, label):
    print(f"== {label}: stage 1 (kappa on |M|>0.08) ==", flush=True)
    r1 = [(one_step(dist, Ms, T, s, d, k, 0, 0, d0, d1, 0.08), k)
          for k in [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]]
    for e, k in sorted(r1):
        print(f"  kappa={k}: high-M 1-step RMSE={e:.4f}", flush=True)
    kbest = sorted(r1)[0][1]
    print(f"== {label}: stage 2 (h0,h1 on all months, kappa={kbest}) ==", flush=True)
    r2 = []
    for h0 in [0.0, 0.1]:
        for h1 in [0.0, 0.5, 1.0, 2.0]:
            e = one_step(dist, Ms, T, s, d, kbest, h0, h1, d0, d1, 0.0)
            r2.append((e, h0, h1))
    for e, h0, h1 in sorted(r2)[:4]:
        print(f"  h0={h0}, h1={h1}: all-month 1-step RMSE={e:.4f}", flush=True)
    ebest, h0b, h1b = sorted(r2)[0]
    return kbest, h0b, h1b


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
    print(f"[{label}] k={k},h0={hh0},h1={hh1}: free-run RMSE={e:.4f} | "
          f"peak actual {fa.max():.2f} @ {t_a.date()}, predicted {fp.max():.2f} @ {t_p.date()}",
          flush=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(idx))
    ax.plot(x, fa, 'k-', lw=2, label='actual')
    ax.plot(x, fp, 'r--', lw=1.6, label='model (M_cs advection + auto-herding)')
    ax.set_xticks(x[::4])
    ax.set_xticklabels([str(d.date())[:7] for d in idx[::4]], rotation=30)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, fig_name), dpi=120)


def main():
    reg = load_regimes()
    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    Ms = load_market_driver(source='panel')

    s, d = entry_exit_rates(reg, '2000-12-01', '2019-12-31')
    T = transition_matrix(reg, '2000-12-01', '2019-12-31')
    p1 = calibrate(dist, Ms, T, s, d, '2000-12-01', '2019-12-31', 'cal 2000-19')
    hindcast(dist, Ms, T, s, d, p1, '2020-01-01', '2021-12-31', '2020-21',
             TOP_BINS, 'hindcast_2020_v3.png', 'fraction in top momentum bin',
             '2020-21 hindcast: M_cs advection + autocatalytic herding')

    s2, d2 = entry_exit_rates(reg, '2000-12-01', '2007-12-31')
    T2 = transition_matrix(reg, '2000-12-01', '2007-12-31')
    p2 = calibrate(dist, Ms, T2, s2, d2, '2000-12-01', '2007-12-31', 'cal 2000-07')
    hindcast(dist, Ms, T2, s2, d2, p2, '2008-01-01', '2009-12-31', '2008-09',
             BOT_BINS, 'hindcast_2008_v3.png',
             'fraction in bottom momentum bin (losers)',
             '2008-09 hindcast: M_cs advection + autocatalytic herding')
    print("DONE", flush=True)


if __name__ == '__main__':
    main()
