"""Long-history analysis: calibrate T_calm/T_stress on 2000-2019, hindcast 2020-21.
Also: hindcast 2008-09 (calibrated 2000-2007) as second validation."""
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

def run_hindcast(reg, mvix, s_rate, d_rate, T_stat, T_calm, T_stress,
                 h0, h1, label):
    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    hd = [d for d in dist.index if pd.Timestamp(h0) <= d <= pd.Timestamp(h1)]
    actual = dist.loc[hd].values
    n0 = actual[0]
    regimes = [mvix.get(d.strftime('%Y-%m'), 0) > 20.0 for d in hd[:-1]]
    print(f"[{label}] {len(hd)} months, {sum(regimes)} stress", flush=True)
    out = {}
    for name, h, use_reg in [
        ("stationary", 0.0, False), ("stationary+h", 0.2, False),
        ("regime", 0.0, True), ("regime+h", 0.2, True),
    ]:
        traj = (simulate_regime(n0, T_calm, T_stress, regimes, h, s_rate, d_rate)
                if use_reg else simulate(n0, T_stat, h, s_rate, d_rate, len(hd) - 1))
        e = rmse(traj, actual)
        out[name] = (e, traj)
        print(f"  {name}: RMSE={e:.4f}", flush=True)
    return hd, actual, out

def main():
    reg = load_regimes()
    mvix = monthly_vix()
    print(f"regimes: {len(reg)} rows", flush=True)

    # --- Test A: 2020-21 hindcast, T from 2000-2019 ---
    s, d = entry_exit_rates(reg, '2000-12-01', '2019-12-31')
    T_stat = transition_matrix(reg, '2000-12-01', '2019-12-31')
    T_calm, T_stress, nc, ns = transition_matrix_by_regime(
        reg, mvix, '2000-12-01', '2019-12-31', thresh=20.0)
    print(f"T_calm {nc} transitions, T_stress {ns} transitions", flush=True)
    # stationary dist of T_stress: does it now know about winners?
    w, v = np.linalg.eig(T_stress.T)
    pi = np.real(v[:, np.argmax(np.real(w))]); pi /= pi.sum()
    print("T_stress stationary top-bin share: %.3f" % pi[TOP_BINS].sum(), flush=True)
    hd, actual, out = run_hindcast(reg, mvix, s, d, T_stat, T_calm, T_stress,
                                   '2020-01-01', '2021-12-31', '2020-21')

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(hd))
    ax.plot(x, actual[:, TOP_BINS].sum(axis=1), 'k-', lw=2, label='actual')
    for name, ls, c in [("stationary", 'b:', None), ("regime+h", 'r--', None)]:
        _, traj = out[name]
        ax.plot(x, traj[:, TOP_BINS].sum(axis=1), ls, lw=1.6, label=name)
    ax.set_xticks(x[::4]); ax.set_xticklabels([str(d.date())[:7] for d in hd[::4]], rotation=30)
    ax.set_ylabel('fraction in top momentum bin')
    ax.set_title('2020-21 hindcast, T calibrated 2000-2019')
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hindcast_2020_longT.png'), dpi=120)

    # --- Test B: 2008-09 hindcast, T from 2000-2007 ---
    s2, d2 = entry_exit_rates(reg, '2000-12-01', '2007-12-31')
    T_stat2 = transition_matrix(reg, '2000-12-01', '2007-12-31')
    T_calm2, T_stress2, nc2, ns2 = transition_matrix_by_regime(
        reg, mvix, '2000-12-01', '2007-12-31', thresh=20.0)
    print(f"2000-07: T_calm {nc2}, T_stress {ns2}", flush=True)
    hd2, actual2, out2 = run_hindcast(reg, mvix, s2, d2, T_stat2, T_calm2, T_stress2,
                                      '2008-01-01', '2009-12-31', '2008-09')
    BOT_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 0]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(hd2))
    ax.plot(x, actual2[:, BOT_BINS].sum(axis=1), 'k-', lw=2, label='actual')
    for name, ls in [("stationary", 'b:'), ("regime+h", 'r--')]:
        _, traj = out2[name]
        ax.plot(x, traj[:, BOT_BINS].sum(axis=1), ls, lw=1.6, label=name)
    ax.set_xticks(x[::4]); ax.set_xticklabels([str(d.date())[:7] for d in hd2[::4]], rotation=30)
    ax.set_ylabel('fraction in bottom momentum bin (losers)')
    ax.set_title('2008-09 hindcast, T calibrated 2000-2007')
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hindcast_2008.png'), dpi=120)
    print("DONE", flush=True)

if __name__ == '__main__':
    main()
