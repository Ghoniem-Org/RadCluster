"""Rev-4 figures: hindcasts, 2022 forecast, paper-trading equity curves, coverage."""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
FIG = os.path.join(HERE, '..', 'doc', 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 11, 'axes.grid': True, 'grid.alpha': 0.3})

def fig_hindcast(label, title):
    d = json.load(open(os.path.join(DATA, f'hindcast_{label}.json')))
    months = pd.to_datetime(d['months'])
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax = axes[0]
    ax.plot(months, d['top_pred'], 'r-', lw=2, label='predicted (conditional, realized M)')
    ax.plot(months, d['top_real'], 'k--', lw=2, label='realized')
    ax.set_ylabel('top-lane cap share (mom>30%)')
    ax.legend(fontsize=10); ax.set_title(title + ' — winner lane')
    ax = axes[1]
    ax.plot(months, d['bot_pred'], 'b-', lw=2, label='predicted (conditional, realized M)')
    ax.plot(months, d['bot_real'], 'k--', lw=2, label='realized')
    ax.set_ylabel('bottom-lane cap share (mom<0%)')
    ax.legend(fontsize=10); ax.set_title(title + ' — loser lane')
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f'hindcast_{label}.png'), dpi=150)
    plt.close(fig)
    print('wrote', label, 'rmse=%.4f' % d['full_rmse'], flush=True)

def fig_paper():
    s = pd.read_csv(os.path.expanduser('~/workspace/goals/stock-market-dynamics/hidden_files/paper_track_record.csv'),
                    parse_dates=['date'])
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    ax = axes[0]
    for col, lab, st in [('paper', 'paper tilt (frozen weights, 5bps)', '-'),
                         ('spy', 'SPY buy-and-hold', '--'),
                         ('ew', 'equal-weight bins', ':')]:
        ax.plot(s['date'], s[col], st, lw=2, label=lab)
    ax.set_ylabel('equity (start=1)')
    ax.legend(fontsize=10); ax.set_title('Paper-trading deployment track record (walk-forward, zero lookahead)')
    ax = axes[1]
    ax.plot(s['date'], s['M_hat'], 'r-', lw=1.5, label='M nowcast')
    ax.plot(s['date'], s['M_real'], 'k-', lw=1, alpha=0.6, label='M realized')
    ax.set_ylabel('market driver M')
    ax.legend(fontsize=10); ax.set_title('Driver nowcast vs realized')
    fig.autofmt_xdate(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'paper_equity.png'), dpi=150)
    plt.close(fig)
    # 2022 zoom
    s22 = s[s['date'] >= '2022-01-01']
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for col, lab, st in [('paper', 'paper tilt', '-'), ('spy', 'SPY', '--'), ('ew', 'equal-weight', ':')]:
        v = s22[col].values / s22[col].values[0]
        ax.plot(s22['date'], v, st, lw=2, label=lab)
    ax.set_ylabel('equity (2022-01=1)'); ax.legend(fontsize=10)
    ax.set_title('2022 genuine out-of-sample: tilt vs benchmarks')
    fig.autofmt_xdate(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'forecast_2022.png'), dpi=150)
    plt.close(fig)
    print('wrote paper figures', flush=True)

def fig_coverage():
    import sys; sys.path.insert(0, os.path.join(HERE, 'src'))
    from model4 import assign_bins
    p = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'), parse_dates=['month'], low_memory=False)
    p['bin'] = assign_bins(p)
    g = p.groupby(p['month'].dt.strftime('%Y-%m')).agg(
        n_mem=('ticker', 'size'), n_px=('has_price', 'sum'),
        n_bin=('bin', lambda s: (s >= 0).sum()))
    g.index = pd.to_datetime(g.index)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(g.index, g['n_mem'], 'k-', lw=1.5, label='index members')
    ax.plot(g.index, g['n_px'], 'b-', lw=1.5, label='with prices')
    ax.plot(g.index, g['n_bin'], 'g-', lw=1.5, label='binnable (all 4 axes)')
    ax.set_ylabel('stocks'); ax.legend(fontsize=10)
    ax.set_title('Historical coverage: true S&P 500 membership vs available data')
    fig.autofmt_xdate(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, 'coverage.png'), dpi=150)
    plt.close(fig)
    print('wrote coverage', flush=True)

if __name__ == '__main__':
    fig_hindcast('2008-09', '2008–09 conditional hindcast')
    fig_hindcast('2020-21', '2020–21 conditional hindcast')
    fig_paper()
    fig_coverage()
