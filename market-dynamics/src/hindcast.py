"""Calibrate herding h (2016-2019), then hindcast 2020-2021 concentration episode."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from model import (load_regimes, monthly_distributions, transition_matrix,
                   entry_exit_rates, simulate, calibrate_h, rmse, N_BIN, bin_mom_q)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'outputs')
FIG = os.path.join(HERE, '..', 'figs')
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)

TOP_BINS = [b for b in range(N_BIN) if bin_mom_q(b) == 4]  # top momentum quintile

def main():
    reg = load_regimes()
    print(f"regimes: {len(reg)} rows", flush=True)

    # --- drift + source/sink measured on 2015-12..2019-12 ---
    T = transition_matrix(reg, '2015-12-01', '2019-12-31')
    s_rate, d_rate = entry_exit_rates(reg, '2015-12-01', '2019-12-31')
    print(f"entry {s_rate:.4f}/mo  exit {d_rate:.4f}/mo", flush=True)
    np.save(os.path.join(OUT, 'T_drift.npy'), T)

    # --- calibrate h on 2016-2019 (calm) AND on 2020-2021 (stress) ---
    cal_calm, _, _ = calibrate_h(reg, T, s_rate, d_rate, '2016-01-01', '2019-12-31')
    cal_stress, _, _ = calibrate_h(reg, T, s_rate, d_rate, '2020-01-01', '2021-12-31')
    print("calm 2016-2019:", " ".join(f"h={h:.1f}:{e:.4f}" for h, e in cal_calm), flush=True)
    print("stress 2020-21:", " ".join(f"h={h:.1f}:{e:.4f}" for h, e in cal_stress), flush=True)
    best_h = min(cal_calm, key=lambda x: x[1])[0]
    best_h_stress = min(cal_stress, key=lambda x: x[1])[0]
    print(f"best h calm = {best_h}, best h stress = {best_h_stress}", flush=True)
    cal = cal_calm
    pd.DataFrame(cal_calm, columns=['h', 'rmse']).to_csv(os.path.join(OUT, 'calibration_calm.csv'), index=False)
    pd.DataFrame(cal_stress, columns=['h', 'rmse']).to_csv(os.path.join(OUT, 'calibration_stress.csv'), index=False)

    fig, ax = plt.subplots(figsize=(7, 4))
    hs, es = zip(*cal)
    ax.plot(hs, es, 'o-')
    ax.axvline(best_h, ls='--', c='r', label=f'best h={best_h}')
    ax.set_xlabel('herding h'); ax.set_ylabel('RMSE (bin shares, 2016-2019)')
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(FIG, 'calibration.png'), dpi=120)

    # --- hindcast 2020-01 -> 2021-12 ---
    dist = monthly_distributions(reg)
    dist.index = pd.to_datetime(dist.index)
    hd = [d for d in dist.index if pd.Timestamp('2020-01-01') <= d <= pd.Timestamp('2021-12-31')]
    actual = dist.loc[hd].values
    n0 = actual[0]
    traj_best = simulate(n0, T, best_h, s_rate, d_rate, len(hd) - 1)
    traj_null = simulate(n0, T, 0.0, s_rate, d_rate, len(hd) - 1)

    top_actual = actual[:, TOP_BINS].sum(axis=1)
    top_best = traj_best[:, TOP_BINS].sum(axis=1)
    top_null = traj_null[:, TOP_BINS].sum(axis=1)

    print(f"hindcast RMSE  h={best_h}: {rmse(traj_best, actual):.4f}", flush=True)
    print(f"hindcast RMSE  h=0 (null): {rmse(traj_null, actual):.4f}", flush=True)
    print(f"top-quintile share Dec-2021: actual {top_actual[-1]:.3f}  "
          f"h={best_h} {top_best[-1]:.3f}  null {top_null[-1]:.3f}", flush=True)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(hd))
    ax.plot(x, top_actual, 'k-', lw=2, label='actual')
    ax.plot(x, top_best, 'r--', lw=1.6, label=f'model h={best_h}')
    ax.plot(x, top_null, 'b:', lw=1.6, label='null h=0 (drift only)')
    ax.set_xticks(x[::4]); ax.set_xticklabels([str(d.date())[:7] for d in hd[::4]], rotation=30)
    ax.set_ylabel('fraction of stocks in top momentum quintile')
    ax.set_title('Hindcast: 2020-21 momentum concentration')
    ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(FIG, 'hindcast_topq.png'), dpi=120)

    # distribution snapshots: actual vs model at 2020-06, 2021-06, 2021-12
    for j, lab in [(5, '2020-06'), (17, '2021-06'), (23, '2021-12')]:
        if j >= len(hd): continue
        fig, ax = plt.subplots(figsize=(9, 3.5))
        w = 0.25; xs = np.arange(N_BIN)
        ax.bar(xs - w, actual[j], w, label='actual')
        ax.bar(xs, traj_best[j], w, label=f'model h={best_h}')
        ax.bar(xs + w, traj_null[j], w, label='null h=0')
        ax.set_xticks(xs); ax.set_xticklabels([f'M{m}V{v}' for m in range(5) for v in range(3)], rotation=60, fontsize=7)
        ax.set_ylabel('fraction of stocks'); ax.set_title(f'regime distribution {lab}')
        ax.legend(fontsize=8); fig.tight_layout()
        fig.savefig(os.path.join(FIG, f'dist_{lab}.png'), dpi=120)

    # RAG schematic
    fig, ax = plt.subplots(figsize=(8, 6))
    rng = np.random.default_rng(0)
    pos = {}
    for b in range(N_BIN):
        mq, vt = bin_mom_q(b), b % 3
        pos[b] = (mq + 0.15 * rng.standard_normal(), vt + 0.15 * rng.standard_normal())
    for i in range(N_BIN):
        for j in range(N_BIN):
            if i != j and T[i, j] > 0.08:
                xi, yi = pos[i]; xj, yj = pos[j]
                ax.annotate('', xy=(xj, yj), xytext=(xi, yi),
                            arrowprops=dict(arrowstyle='->', color='steelblue',
                                            alpha=min(0.9, T[i, j] * 4), lw=0.8))
    for b in range(N_BIN):
        x, y = pos[b]
        c = plt.cm.Reds(0.25 + 0.6 * bin_mom_q(b) / 4)
        ax.scatter([x], [y], s=520, c=[c], ec='k', zorder=3)
        ax.text(x, y, f'{b}', ha='center', va='center', fontsize=8, zorder=4)
    # herding annotation: red drift toward high-momentum column
    ax.annotate('', xy=(4.3, 1), xytext=(2.5, 1),
                arrowprops=dict(arrowstyle='->', color='crimson', lw=2.5))
    ax.text(3.4, 1.25, 'herding h: flow follows winners', color='crimson', fontsize=9, ha='center')
    ax.set_xlim(-0.6, 4.9); ax.set_ylim(-0.6, 2.9)
    ax.set_xlabel('momentum quintile ->'); ax.set_ylabel('volatility tercile')
    ax.set_title('Regime RAG: 15 clusters, drift (blue) + herding ratchet (red)')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'rag.png'), dpi=120)
    print("FIGURES DONE", flush=True)

if __name__ == '__main__':
    main()
