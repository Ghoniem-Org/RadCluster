"""Rev-5 figures: ensemble forecast intervals, fragility index,
concentration-test result, data-fix before/after."""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'outputs')
FIG = os.path.join(HERE, '..', 'doc', 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 11, 'axes.grid': True, 'grid.alpha': 0.3})


def fig_ensemble():
    d = json.load(open(os.path.join(OUT, 'ensemble_forecast.json')))
    h = pd.to_datetime(d['horizons'])
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for ax, key, lab in [(axes[0], 'top_lane', 'top momentum lane share'),
                         (axes[1], 'hhi', 'HHI (24-bin concentration)')]:
        arr = np.array(d[key])  # (H, 3): p5, p50, p95
        ax.fill_between(h, arr[:, 0], arr[:, 2], alpha=0.25, label='90% ensemble interval')
        ax.plot(h, arr[:, 1], 'r-', lw=2, label='ensemble median')
        ax.set_ylabel(lab)
        ax.legend(fontsize=10)
    axes[0].set_title('Rev-5 ensemble distribution forecast from 2022-12 (200 scenarios)\n'
                      'Dirichlet-T + AR(1)-M + VIX-regime Markov; intervals on every forecast')
    fig.autofmt_xdate(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'ensemble_forecast.png'), dpi=150)
    plt.close(fig)
    print('wrote ensemble_forecast', flush=True)


def fig_fragility():
    f = pd.read_csv(os.path.join(OUT, 'fragility.csv'), parse_dates=['month'])
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax = axes[0]
    ax.plot(f['month'], f['C_top'], 'k-', lw=1.5, label='top-lane cap share C(t)')
    ax.axhline(f['C_top'].quantile(0.9), color='r', ls='--', lw=1, label='90th pctile (history)')
    ax.axhline(f['C_top'].quantile(0.5), color='b', ls=':', lw=1, label='median')
    cur = f.iloc[-1]
    ax.plot(cur['month'], cur['C_top'], 'ro', ms=8, label='2022-12: %.3f (%dth pctile)' % (
        cur['C_top'], round(100 * cur['C_pctile_hist'])))
    ax.set_ylabel('C(t)'); ax.legend(fontsize=10)
    ax.set_title('Crowding / fragility index: top-lane share vs its own history')
    ax = axes[1]
    ax.plot(f['month'], f['hhi'], 'k-', lw=1.5, label='HHI(t)')
    ax.axhline(f['hhi'].quantile(0.9), color='r', ls='--', lw=1, label='90th pctile')
    ax.plot(cur['month'], cur['hhi'], 'ro', ms=8, label='2022-12: %.4f (%dth pctile)' % (
        cur['hhi'], round(100 * cur['hhi_pctile_hist'])))
    ax.set_ylabel('HHI'); ax.legend(fontsize=10)
    fig.autofmt_xdate(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'fragility.png'), dpi=150)
    plt.close(fig)
    print('wrote fragility', flush=True)


def fig_concentration_test():
    s = pd.read_csv(os.path.join(OUT, 'conc_ret_monthly.csv'), parse_dates=['month'])
    est = s[(s['month'] >= '2001-01-01') & (s['month'] < '2008-01-01')].copy()
    wf = s[s['month'] >= '2008-01-01'].copy()
    qs = est['C'].quantile([0.2, 0.4, 0.6, 0.8]).values
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    for ax, d, title in [(axes[0], est, 'Estimation 2001-2007 (n=84)'),
                         (axes[1], wf, 'Walk-forward 2008-2022 (n=179)')]:
        d['q'] = np.clip(np.digitize(d['C'].values, qs), 0, 4)
        means, ses, ns = [], [], []
        for q in range(5):
            v = d[d['q'] == q]['re_top'].dropna().values
            means.append(v.mean()); ses.append(v.std(ddof=1) / np.sqrt(len(v))); ns.append(len(v))
        means, ses = np.array(means), np.array(ses)
        x = np.arange(5)
        ax.bar(x, 100 * means, yerr=100 * ses, capsize=4, alpha=0.7)
        ax.set_xticks(x); ax.set_xticklabels([f'Q{q+1}\n(n={n})' for q, n in enumerate(ns)])
        ax.axhline(0, color='k', lw=0.8)
        ax.set_title(title)
        ax.set_ylabel('next-month top-lane excess return (%/mo)\n±1 SE')
    axes[0].set_xlabel('concentration quintile (from estimation)')
    axes[1].set_xlabel('concentration quintile (from estimation)')
    fig.suptitle('Concentration -> return test: no monotone crowding penalty (honest negative)\n'
                 'Q5-Q1 spread: in-sample -0.59% (t=-0.40), walk-forward -0.52% (t=-0.52)')
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'concentration_test.png'), dpi=150)
    plt.close(fig)
    print('wrote concentration_test', flush=True)


def fig_datafix():
    b = pd.read_csv(os.path.join(DATA, 'panel_v2_preunitfix.csv'),
                    parse_dates=['month'], low_memory=False)
    p = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'),
                    parse_dates=['month'], low_memory=False)
    g1 = b.groupby('month')['mcap'].max()
    g2 = p.groupby('month')['mcap'].max()
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.semilogy(g1.index, g1.values / 1e12, 'r-', lw=1.5, label='before fix (max single-stock mcap)')
    ax.semilogy(g2.index, g2.values / 1e12, 'g-', lw=1.5, label='after fix')
    ax.axhline(4, color='k', ls='--', lw=1, label='$4T sanity bound')
    ax.set_ylabel('max single-stock mcap ($T, log scale)')
    ax.legend(fontsize=10)
    ax.set_title('Data integrity fix: consistent unit errors in filed share counts\n'
                 'SRE backfill read $5,000T; segment-level rescaling to absolute anchors')
    fig.autofmt_xdate(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'data_fix.png'), dpi=150)
    plt.close(fig)
    print('wrote data_fix', flush=True)


if __name__ == '__main__':
    fig_ensemble()
    fig_fragility()
    fig_concentration_test()
    fig_datafix()
