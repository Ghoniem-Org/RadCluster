"""Four-axis RAG PNG: 2x2 panels (turnover x size), each 3 mom x 2 vol."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'doc', 'figures')

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle('Reaction-Admissibility Graph: 24 regime bins (3 momentum x 2 volatility x 2 turnover x 2 size)',
             fontsize=13, fontweight='bold')
panels = [('turnover LOW', 'size SMALL'), ('turnover LOW', 'size LARGE'),
          ('turnover HIGH', 'size SMALL'), ('turnover HIGH', 'size LARGE')]
mom_labels = ['m0: r<0%', 'm1: 0-30%', 'm2: r>30%']
vol_labels = ['v0: vol<28%', 'v1: vol>28%']
for ax, (t_lab, s_lab) in zip(axes.flat, panels):
    ax.set_xlim(-0.6, 2.6); ax.set_ylim(-0.6, 3.6)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f'{t_lab} / {s_lab}', fontsize=11, fontweight='bold')
    # 3 rows (momentum) x 2 cols (volatility)
    pos = {}
    for m in range(3):
        for v in range(2):
            x, y = v * 1.6 + 0.3, (2 - m) * 1.1 + 0.3
            pos[(m, v)] = (x, y)
            ax.add_patch(plt.Circle((x, y), 0.28, fc='lightblue', ec='navy', lw=1.5))
            ax.text(x, y, f'{m}{v}', ha='center', va='center', fontsize=9, fontweight='bold')
    # momentum advection edges (vertical, labeled)
    for v in range(2):
        for m in [0, 1]:
            x0, y0 = pos[(m, v)]; x1, y1 = pos[(m + 1, v)]
            ax.add_patch(FancyArrowPatch((x0, y0 - 0.28), (x1, y1 + 0.28),
                                         arrowstyle='->', mutation_scale=14, color='red', lw=1.5))
            ax.add_patch(FancyArrowPatch((x1, y1 + 0.28), (x0, y0 - 0.28),
                                         arrowstyle='->', mutation_scale=14, color='blue', lw=1.5))
    # vol drift edges (horizontal, gray)
    for m in range(3):
        x0, y0 = pos[(m, 0)]; x1, y1 = pos[(m, 1)]
        ax.add_patch(FancyArrowPatch((x0 + 0.28, y0), (x1 - 0.28, y1),
                                     arrowstyle='<->', mutation_scale=10, color='gray', lw=1, ls='--'))
    # axis labels
    for m in range(3):
        ax.text(-0.45, pos[(m, 0)][1], mom_labels[m], ha='right', va='center', fontsize=8, color='darkred')
    for v in range(2):
        ax.text(pos[(0, v)][0], 3.35, vol_labels[v], ha='center', va='center', fontsize=8, color='darkgreen')
# legend
fig.text(0.5, 0.02,
         'Red/blue arrows: market advection along momentum (kappa*M, fractional remap) | '
         'Gray dashed: volatility drift (measured T) | Drift T is dense 24x24 (not all shown)',
         ha='center', fontsize=9, style='italic')
fig.tight_layout(rect=[0, 0.04, 1, 0.95])
fig.savefig(os.path.join(FIG, 'rag_4axis.png'), dpi=150)
print('wrote rag_4axis.png')
