"""Rev-5 risk metric 2: crowding / fragility index.

C(t) = top-momentum-lane cap share; HHI(t) = sum of squared bin shares.
Both measured monthly on the fixed panel (cap-weighted, PIT).
Fragility index = percentile of the current reading against its own history.
Reported two ways:
  - full-history percentile (diagnostic, uses all 2001-2022),
  - expanding-window percentile (what was knowable at each month-end t,
    zero lookahead).
Saves outputs/fragility.csv.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'outputs')
from model4 import assign_bins, cap_state


def main():
    panel = pd.read_csv(os.path.join(DATA, 'panel_v2.csv'),
                        parse_dates=['month'], low_memory=False)
    panel['bin'] = assign_bins(panel)
    months = sorted(panel['month'].unique())
    rows = []
    for m in months:
        c, info = cap_state(panel[panel['month'] == m])
        if info['n_bin'] < 100:
            continue
        top = float(c[[b for b in range(24) if b // 8 == 2]].sum())
        hhi = float((c ** 2).sum())
        rows.append((m, top, hhi, info['n_bin']))
    f = pd.DataFrame(rows, columns=['month', 'C_top', 'hhi', 'n_bin'])
    # full-history percentiles (diagnostic)
    f['C_pctile_hist'] = f['C_top'].rank(pct=True)
    f['hhi_pctile_hist'] = f['hhi'].rank(pct=True)
    # expanding-window percentiles (zero lookahead)
    f['C_pctile_exp'] = [float((f['C_top'].iloc[:i + 1] <= v).mean())
                         for i, v in enumerate(f['C_top'].values)]
    f['hhi_pctile_exp'] = [float((f['hhi'].iloc[:i + 1] <= v).mean())
                           for i, v in enumerate(f['hhi'].values)]
    f.to_csv(os.path.join(OUT, 'fragility.csv'), index=False)
    cur = f.iloc[-1]
    print('fragility index (%d months):' % len(f), flush=True)
    print('  current %s: C_top=%.3f (hist %dth pctile, expanding %dth)' % (
        cur['month'].date(), cur['C_top'],
        round(100 * cur['C_pctile_hist']), round(100 * cur['C_pctile_exp'])), flush=True)
    print('  current %s: HHI=%.4f (hist %dth pctile, expanding %dth)' % (
        cur['month'].date(), cur['hhi'],
        round(100 * cur['hhi_pctile_hist']), round(100 * cur['hhi_pctile_exp'])), flush=True)
    # most fragile episodes (top-5 by expanding C percentile)
    top5 = f.nlargest(5, 'C_pctile_exp')[['month', 'C_top', 'C_pctile_exp']]
    print('  most crowded month-ends (expanding):', flush=True)
    for _, r in top5.iterrows():
        print('    %s C=%.3f pctile=%d' % (
            r['month'].date(), r['C_top'], round(100 * r['C_pctile_exp'])), flush=True)
    print('wrote outputs/fragility.csv', flush=True)


if __name__ == '__main__':
    main()
