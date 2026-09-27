"""Regime-switching hindcast: T depends on VIX regime (calm vs stress)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model import (load_regimes, monthly_distributions, transition_matrix,
                   transition_matrix_by_regime, entry_exit_rates, simulate,
                   simulate_regime, monthly_vix, rmse, N_BIN, bin_mom_q)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
FIG = os.path.join(HERE, '..', 'figs')

TOP_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 4]

def main():
    reg = load_regimes()
    mvix = monthly_vix()
    s_rate, d_rate = entry_exit_rates(reg, '2015-12-01', '2019-12-31')

    T_stat = transition_matrix(reg, '2015-12-01', '2019-12-31')
    T_calm, T_stress, n_calm, n_stress = transition_matrix_by_regime(
        reg, mvix, '2015-12-01', '2019-12-31', thresh=20.0)
    print(f"T_calm from {n_calm} transitions, T_stress from {n_stress} transitions", flush=True)
    np.save(os.path.join(OUT, 'T_calm.npy'), T_calm)
    np.save(os.path.join(OUT, 'T_stress.npy'), T_stress)

    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    hd = [d for d in dist.index if pd.Timestamp('2020-01-01') <= d <= pd.Timestamp('2021-12-31')]
    actual = dist.loc[hd].values
    n0 = actual[0]
    # regime label per step: stress if origin month VIX > 20 (perfect-foresight label)
    regimes = [mvix.get(d.strftime('%Y-%m'), 0) > 20.0 for d in hd[:-1]]
    print("stress months in hindcast:", sum(regimes), "of", len(regimes), flush=True)

    results = {}
    for name, fn in [
        ("stationary h=0", lambda h: simulate(n0, T_stat, h, s_rate, d_rate, len(hd) - 1)),
        ("regime h=0", lambda h: simulate_regime(n0, T_calm, T_stress, regimes, h, s_rate, d_rate)),
    ]:
        for h in [0.0, 0.2, 0.5, 1.0]:
            traj = fn(h)
            e = rmse(traj, actual)
            top = traj[:, TOP_BINS].sum(axis=1)[-1]
            results[(name, h)] = (e, top)
            print(f"{name} h={h}: RMSE={e:.4f} Dec21-top={top:.3f}", flush=True)
    print(f"actual Dec21 top-quintile share: {actual[:, TOP_BINS].sum(axis=1)[-1]:.3f}", flush=True)

    # plot: best of each family
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(hd))
    ax.plot(x, actual[:, TOP_BINS].sum(axis=1), 'k-', lw=2, label='actual')
    ax.plot(x, simulate(n0, T_stat, 0, s_rate, d_rate, len(hd)-1)[:, TOP_BINS].sum(axis=1),
            'b:', lw=1.6, label='stationary T, h=0')
    ax.plot(x, simulate_regime(n0, T_calm, T_stress, regimes, 0.5, s_rate, d_rate)[:, TOP_BINS].sum(axis=1),
            'r--', lw=1.6, label='regime T, h=0.5')
    ax.set_xticks(x[::4]); ax.set_xticklabels([str(d.date())[:7] for d in hd[::4]], rotation=30)
    ax.set_ylabel('fraction of stocks in top momentum bin')
    ax.set_title('Hindcast: stationary vs regime-switching drift')
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hindcast_regime.png'), dpi=120)
    print("DONE", flush=True)

if __name__ == '__main__':
    main()
