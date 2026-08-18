"""Score a results .jsonl against the six-observable target set.

Usage:  python score_six_targets.py results/S8_15dpa.jsonl [design/S8_labels.json]

The six are N_100, d_100, N_111, d_111, N_voids, d_cavity at EUROFER97
330 C / 15 dpa (targets_330C_15dpa.json).  f_100 is NOT scored as an
independent observable -- it is N_100/(N_100+N_111), so it is reported as a
derived check only.

A row "passes" an observable when the model value lies inside [lo, hi].  The
headline number is how many of the six are in range; ratio-to-target is shown
so a near-miss is distinguishable from a three-order-of-magnitude miss.

Rows that did not reach the target dose, or that piled at the grid ceiling,
are marked -- their size readouts are not measurements and must not be scored
as if they were.
"""
import sys, json, math
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

HERE = Path(__file__).resolve().parent
TGT = json.load(open(HERE / 'targets_330C_15dpa.json'))
OBS = TGT['observables']
ORDER = ['N_loops_100', 'd_100_nm', 'N_loops_111', 'd_111_nm', 'N_voids', 'd_cavity_nm']
SHORT = {'N_loops_100': 'N_100', 'd_100_nm': 'd_100', 'N_loops_111': 'N_111',
         'd_111_nm': 'd_111', 'N_voids': 'N_void', 'd_cavity_nm': 'd_void'}

res = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / 'results/S8_15dpa.jsonl'
labf = Path(sys.argv[2]) if len(sys.argv) > 2 else res.parent.parent / 'design/S8_labels.json'
labels = json.load(open(labf)) if labf.exists() else {}
rows = sorted((json.loads(l) for l in open(res)), key=lambda r: r['row_id'])

print('target set: %s %s C / %s dpa  (%s)' % (
    TGT['condition']['material'], TGT['condition']['T_C'],
    TGT['condition']['dose_dpa'], TGT['condition']['irradiation']))
print()
hdr = '%-20s %6s' % ('row', 'dose') + ''.join('%12s' % SHORT[k] for k in ORDER) + '%7s %8s %6s' % ('in-rng', 'f_100', 'flags')
print(hdr); print('-' * len(hdr))

for r in rows:
    lab = labels.get(str(r['row_id']), str(r['row_id']))
    cells, n_ok = [], 0
    for k in ORDER:
        v = r.get(k)
        o = OBS[k]
        if v is None:
            cells.append('%12s' % '-'); continue
        ok = o['lo'] <= v <= o['hi']
        n_ok += ok
        rat = v / o['target'] if o['target'] else float('nan')
        s = ('%.2e' % v) if v >= 1e4 else ('%.2f' % v)
        cells.append('%12s' % (s + ('*' if ok else '')))
    N1, N0 = r.get('N_loops_111'), r.get('N_loops_100')
    f100 = N0 / (N0 + N1) if (N0 and N1) else float('nan')
    flags = []
    if r.get('grid_limited'):
        flags.append('GRID')
    if abs(r.get('dose_reached', 0) - TGT['condition']['dose_dpa']) > 0.02 * TGT['condition']['dose_dpa']:
        flags.append('DOSE')
    if r.get('starved'):
        flags.append('STARV')
    print('%-20s %6.3g' % (lab, r.get('dose_reached', 0)) + ''.join(cells)
          + '%7s %8.3f %6s' % ('%d/6' % n_ok, f100, ','.join(flags) or '-'))

print()
print('  * = inside the experimental range.  GRID = piled at the grid ceiling')
print('  (size readout invalid).  DOSE = did not reach 15 dpa (nothing is')
print('  comparable).  Ratios to target:')
print()
print('%-20s' % 'row' + ''.join('%12s' % SHORT[k] for k in ORDER))
for r in rows:
    lab = labels.get(str(r['row_id']), str(r['row_id']))
    out = []
    for k in ORDER:
        v, t = r.get(k), OBS[k]['target']
        out.append('%12s' % ('-' if v is None else '%.3gx' % (v / t)))
    print('%-20s' % lab + ''.join(out))
