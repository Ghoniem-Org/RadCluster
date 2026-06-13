import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

# ── Figure setup ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(28, 21))
ax.set_xlim(0, 28)
ax.set_ylim(0, 21)
ax.axis('off')
fig.patch.set_facecolor('white')

# ── Font sizes ────────────────────────────────────────────────────────────────
FS_TITLE   = 30
FS_SUBT    = 24
FS_NODE    = 20
FS_BOX     = 20
FS_LABEL   = 20
FS_EQ      = 20
FS_KEY_HDR = 20
FS_KEY     = 20

# ── Colours ───────────────────────────────────────────────────────────────────
C_SIA   = '#C62828'
C_CAV   = '#1565C0'
C_HE    = '#1B5E20'
C_SRC   = '#004D40'
C_SINK  = '#37474F'
C_P3    = '#C62828'
C_P5    = '#E65100'
C_P1    = '#6A1B9A'
C_P2V   = '#1565C0'
C_P2H   = '#2E7D32'
C_P2I   = '#00695C'
C_P7    = '#558B2F'
C_P8    = '#546E7A'
C_P4    = '#37474F'
C_SRC_ARROW = '#90A4AE'   # source-box connector arrows (was white, invisible on white bg)


# ─────────────────────────────────────────────────────────────────────────────
# helpers
# ─────────────────────────────────────────────────────────────────────────────
def circle(ax, cx, cy, r, color, zorder=3):
    c = plt.Circle((cx, cy), r, color=color, zorder=zorder, linewidth=0)
    ax.add_patch(c)

def node_text(ax, cx, cy, txt, fs=FS_NODE, color='white', zorder=4):
    ax.text(cx, cy, txt, ha='center', va='center', fontsize=fs,
            color=color, fontweight='bold', zorder=zorder)

def fancy_box(ax, x, y, w, h, fc, ec='white', lw=2.0, radius=0.15, zorder=2):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad={radius}",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
    ax.add_patch(p)

def arr(ax, x1, y1, x2, y2, color, lw=2.2, ls='-', rad=0.0, zorder=2):
    cs = f'arc3,rad={rad}' if rad != 0 else 'arc3,rad=0'
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                linestyle=ls,
                                connectionstyle=cs),
                zorder=zorder)

def labeled_box(ax, x, y, w, h, title, sub, fc=C_SRC):
    fancy_box(ax, x, y, w, h, fc=fc, ec='white', lw=2.0)
    ax.text(x+w/2, y+h*0.65, title, ha='center', va='center',
            fontsize=FS_BOX, color='white', fontweight='bold')
    ax.text(x+w/2, y+h*0.25, sub,  ha='center', va='center',
            fontsize=FS_LABEL, color='white')


# ─────────────────────────────────────────────────────────────────────────────
# TITLE
# ─────────────────────────────────────────────────────────────────────────────
ax.text(14, 20.5,
        'EUROFER97 Reaction Network — Graph-Walker Representation',
        ha='center', va='center', fontsize=FS_TITLE, fontweight='bold')
ax.text(14, 19.85,
        r'$c_n$: SIA clusters,    $c_{m,\ell}$: cavities,    $c_h$: free He',
        ha='center', va='center', fontsize=FS_SUBT)


# ─────────────────────────────────────────────────────────────────────────────
# SOURCE BOXES
# ─────────────────────────────────────────────────────────────────────────────
# Displacement cascade  (left)
labeled_box(ax, 0.4, 18.0, 6.0, 1.5,
            'Displacement cascade', r'$G_n,\ G_{m,\ell}$')

# He production  (right)
labeled_box(ax, 21.5, 18.0, 6.0, 1.5,
            'He production / transmutation', r'$G_{\rm He}$')

# Graph-Walker equation  (centre)
fancy_box(ax, 8.8, 17.8, 10.3, 1.7, fc='#FFF9C4', ec='#F57F17', lw=2.5)
ax.text(13.95, 19.1,
        'Graph Walker (abstract master equation)',
        ha='center', va='center', fontsize=FS_LABEL, fontweight='bold',
        color='#E65100')
ax.text(13.95, 18.3,
        r'$\dfrac{dc_a}{dt} = G_a + \sum_{r\in\{P1\ldots P8\}}'
        r'S_{ad}\,r_r(\mathbf{c}) - D_a c_a$',
        ha='center', va='center', fontsize=FS_EQ)


# ─────────────────────────────────────────────────────────────────────────────
# SIA CLUSTER CHAIN   c_1  c_2  c_3  ...  c_n
# ─────────────────────────────────────────────────────────────────────────────
SIA_Y  = 14.5
NR     = 0.60           # node radius
SIA_X  = [1.4, 3.2, 5.0, 6.8]
SIA_L  = [r'$c_1$', r'$c_2$', r'$c_3$', r'$c_n$']

ax.text(0.3, 16.2, r'$c_I\;(\sigma=+1)$',
        ha='left', va='center', fontsize=FS_LABEL+1,
        color=C_SIA, fontweight='bold')

for x, lbl in zip(SIA_X, SIA_L):
    circle(ax, x, SIA_Y, NR, C_SIA)
    node_text(ax, x, SIA_Y, lbl)

# dots between c_3 and c_n
ax.text(5.9, SIA_Y, r'$\cdots$', ha='center', va='center',
        fontsize=22, color=C_SIA, zorder=3)

# P3 growth arrows  (forward, top half of node gap)
for i in range(len(SIA_X)-1):
    x1 = SIA_X[i] + NR + 0.05
    x2 = SIA_X[i+1] - NR - 0.05
    if i == 2:
        x1 += 0.25; x2 -= 0.1
    arr(ax, x1, SIA_Y+0.18, x2, SIA_Y+0.18, C_P3, lw=2.5)

ax.text(3.2, SIA_Y+1.15,
        r'$P3:\ \kappa_{\ell\ell}^{loop}c_n c_1$',
        ha='center', va='center', fontsize=FS_LABEL,
        color=C_P3, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=C_P3, lw=1.8))

# P5i emission arrows (backward, bottom half)
for i in range(1, len(SIA_X)):
    x1 = SIA_X[i] - NR - 0.05
    x2 = SIA_X[i-1] + NR + 0.05
    if i == 3:
        x1 -= 0.1; x2 += 0.25
    arr(ax, x1, SIA_Y-0.20, x2, SIA_Y-0.20, C_P5, lw=2.0, ls='dashed')

ax.text(4.0, SIA_Y-1.1,
        r'$P5i:\ \varepsilon_I(q)c_n$',
        ha='center', va='center', fontsize=FS_LABEL,
        color=C_P5, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=C_P5, lw=1.8))


# ─────────────────────────────────────────────────────────────────────────────
# FREE He NODE
# ─────────────────────────────────────────────────────────────────────────────
HE_X, HE_Y = 10.5, 14.5
fancy_box(ax, HE_X-1.2, HE_Y-0.82, 2.4, 1.65,
          fc=C_HE, ec='white', lw=2.5, zorder=3)
ax.text(HE_X, HE_Y+0.22, r'$c_h$',
        ha='center', va='center', fontsize=FS_NODE+4,
        color='white', fontweight='bold', zorder=5)
ax.text(HE_X, HE_Y-0.40, '(free He)',
        ha='center', va='center', fontsize=FS_LABEL-1,
        color='white', zorder=5)


# ─────────────────────────────────────────────────────────────────────────────
# CAVITY GRID   c_{m,l}
# ─────────────────────────────────────────────────────────────────────────────
CR       = 0.58          # cavity node radius
CAV_MX   = [13.5, 17.0, 22.5]   # m=1, m=2, m columns
CAV_LY   = [14.5, 12.0,  9.5]   # l=0, l=1, l=2 rows
M_LABELS = ['m=1', 'm=2', 'm']
L_LABELS = [r'$\ell=0$', r'$\ell=1$', r'$\ell=2$']

ax.text(12.5, 16.3, r'$c_V\;(\sigma=-1)$',
        ha='left', va='center', fontsize=FS_LABEL+1,
        color=C_CAV, fontweight='bold')

for j, (mx, ml) in enumerate(zip(CAV_MX, M_LABELS)):
    ax.text(mx, 15.55, ml, ha='center', va='center',
            fontsize=FS_LABEL, color=C_CAV, fontweight='bold')

for i, (ly, ll) in enumerate(zip(CAV_LY, L_LABELS)):
    ax.text(12.0, ly, ll, ha='center', va='center',
            fontsize=FS_LABEL, color=C_CAV)

# draw cavity nodes
for i, ly in enumerate(CAV_LY):
    for j, mx in enumerate(CAV_MX):
        m_str = str(j+1) if j < 2 else 'm'
        lbl   = rf'$c_{{{m_str},{i}}}$'
        circle(ax, mx, ly, CR, C_CAV)
        node_text(ax, mx, ly, lbl, fs=FS_NODE-1)

# dots between col m=2 and col m
for ly in CAV_LY:
    ax.text(19.8, ly, r'$\cdots$', ha='center', va='center',
            fontsize=22, color=C_CAV, zorder=3)
# dotted extension line on l=0 row
ax.plot([CAV_MX[2]+CR, CAV_MX[2]+1.4], [CAV_LY[0], CAV_LY[0]],
        ':', color=C_CAV, lw=2.2)

# P2v cavity growth arrows  (horizontal, l=0 row)
arr(ax, CAV_MX[0]+CR+0.05, CAV_LY[0]+0.18,
       CAV_MX[1]-CR-0.05, CAV_LY[0]+0.18, C_P2V, lw=2.5)
arr(ax, CAV_MX[1]+CR+0.05, CAV_LY[0]+0.18,
       CAV_MX[2]-CR-0.9,  CAV_LY[0]+0.18, C_P2V, lw=2.5)

ax.text(19.5, 16.6,
        r'$P2v:\ \kappa_{V,m}^{cav}c_{1,0}\,\tilde{c}_m$',
        ha='center', va='center', fontsize=FS_LABEL, color=C_P2V,
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=C_P2V, lw=1.8))

# P2h He-trapping arrows  (vertical, downward, left side of node)
for mx in CAV_MX:
    for i in range(len(CAV_LY)-1):
        y1 = CAV_LY[i]   - CR - 0.05
        y2 = CAV_LY[i+1] + CR + 0.05
        arr(ax, mx-0.12, y1, mx-0.12, y2, C_P2H, lw=2.5)

ax.text(24.5, 13.2,
        r'$P2h:\ \kappa_{h,m}^{cav}c_h\,c_{m,\ell}$',
        ha='left', va='center', fontsize=FS_LABEL, color=C_P2H,
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=C_P2H, lw=1.8))
ax.text(24.5, 12.55,
        'P2h: He trap\nP5h+P8: emit/resolute',
        ha='left', va='center', fontsize=FS_LABEL-1, color=C_P2H)

# P5v emission arrows  (vertical, upward, right side of node)
for mx in CAV_MX:
    for i in range(1, len(CAV_LY)):
        y1 = CAV_LY[i]   + CR + 0.05
        y2 = CAV_LY[i-1] - CR - 0.05
        arr(ax, mx+0.12, y1, mx+0.12, y2, C_P5, lw=2.0, ls='dashed')

ax.text(15.25, 13.25,
        r'$P5v:\ \varepsilon_v(m,\ell)c_{m,\ell}$',
        ha='center', va='center', fontsize=FS_LABEL, color=C_P5,
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=C_P5, lw=1.8))

# P2i/P6  SIA absorption by cavity  (He → cavity l=0)
arr(ax, HE_X+1.2, HE_Y+0.25, CAV_MX[0]-CR-0.05, CAV_LY[0]+0.20,
    C_P2I, lw=2.5)
ax.text(10.5, 13.15,
        r'$P2i/P6:\ \kappa_{i,m}^{cav}c_n\,c_m$',
        ha='center', va='center', fontsize=FS_LABEL, color=C_P2I,
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.25', fc='#E8F5E9', ec=C_P2I, lw=1.8))

# P1  V-I recombination  (c_1 → c_{1,0})
arr(ax, SIA_X[0]+NR*0.7, SIA_Y-NR*0.7,
       CAV_MX[0]-CR*0.7, CAV_LY[0]-CR*0.7,
    C_P1, lw=2.5, rad=0.15)
ax.text(7.8, 12.1,
        r'$P1:\ \mathcal{K}_{iV}c_1 c_{1,0}$',
        ha='center', va='center', fontsize=FS_LABEL, color=C_P1,
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=C_P1, lw=1.8))

# P7  trap mutation  V_{m,l} → V_{m+1,l} + I_1
# arrow from cavity to SIA side (escaping interstitial)
ax.annotate('', xy=(SIA_X[0]-NR*0.5, SIA_Y-NR*0.85),
            xytext=(CAV_MX[0]+CR*0.5, CAV_LY[1]-CR*0.5),
            arrowprops=dict(arrowstyle='->', color=C_P7, lw=2.2,
                            connectionstyle='arc3,rad=-0.35'),
            zorder=2)
ax.text(8.5, 11.0,
        r'$P7:\ \Gamma_{TM}c_{m,\ell}\!\to\!I_1$',
        ha='center', va='center', fontsize=FS_LABEL, color=C_P7,
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=C_P7, lw=1.8))


# ─────────────────────────────────────────────────────────────────────────────
# SOURCE → NETWORK arrows
# ─────────────────────────────────────────────────────────────────────────────
# Displacement cascade → c_2 (SIA)
ax.annotate('', xy=(SIA_X[1], SIA_Y+NR),
            xytext=(3.2, 18.0),
            arrowprops=dict(arrowstyle='->', color=C_SRC_ARROW, lw=1.8,
                            linestyle='dashed',
                            connectionstyle='arc3,rad=0.0'), zorder=5)
# Displacement cascade → c_{1,0}
ax.annotate('', xy=(CAV_MX[0], CAV_LY[0]+CR),
            xytext=(6.0, 18.0),
            arrowprops=dict(arrowstyle='->', color=C_SRC_ARROW, lw=1.8,
                            linestyle='dashed',
                            connectionstyle='arc3,rad=-0.15'), zorder=5)
# He source → c_h
ax.annotate('', xy=(HE_X, HE_Y+0.83),
            xytext=(22.0, 18.0),
            arrowprops=dict(arrowstyle='->', color=C_SRC_ARROW, lw=1.8,
                            linestyle='dashed',
                            connectionstyle='arc3,rad=0.2'), zorder=5)


# ─────────────────────────────────────────────────────────────────────────────
# SINK BOXES   — placed lower to avoid overlap with process key
# ─────────────────────────────────────────────────────────────────────────────
SINK_Y = 5.6      # top of sink boxes
SINK_H = 1.4
SINK_W = 5.8

labeled_box(ax,  0.8, SINK_Y-SINK_H, SINK_W, SINK_H,
            'Dislocation network', r'$D_\alpha^d$',   fc=C_SINK)
labeled_box(ax, 10.9, SINK_Y-SINK_H, SINK_W, SINK_H,
            'Grain boundaries',    r'$D_\alpha^{gb}$', fc=C_SINK)
labeled_box(ax, 21.2, SINK_Y-SINK_H, SINK_W, SINK_H,
            'Precipitates',        r'$D_\alpha^p$',   fc=C_SINK)

# P4 label
ax.text(14.0, 7.2,
        r'$P4:\ D_\alpha c_\alpha$',
        ha='center', va='center', fontsize=FS_LABEL+1, color=C_P4,
        fontweight='bold')

# Dashed P4 arrows from all cluster nodes down to sinks
sink_centres = [3.7, 13.8, 24.1]
source_nodes = (
    [(x, SIA_Y-NR)   for x in SIA_X] +
    [(HE_X, HE_Y-0.83)] +
    [(mx, CAV_LY[-1]-CR) for mx in CAV_MX]
)
sink_targets = (
    [(3.7,  SINK_Y-SINK_H)] * 4 +
    [(13.8, SINK_Y-SINK_H)] +
    [(24.1, SINK_Y-SINK_H)] * 3
)
for (x1,y1),(x2,y2) in zip(source_nodes, sink_targets):
    arr(ax, x1, y1, x2, y2, C_P4, lw=1.5, ls='dashed', zorder=1)


# ─────────────────────────────────────────────────────────────────────────────
# PROCESS KEY  — placed BELOW all sink boxes
# ─────────────────────────────────────────────────────────────────────────────
KY = 3.9    # top of key area (well below SINK_Y-SINK_H = 4.2)

# Three columns of three keys, each aligned under its sink box
key_columns = [
    # under Dislocation network
    [
        (C_P3,  '-',  r'P3: loop growth  $I_{n-1}+I_1\to I_n$'),
        (C_P5,  '--', r'P5: SIA emission (inverse of P3)'),
        (C_P1,  '-',  r'P1: V-I recombination  $I_1+V_1\to\emptyset$'),
    ],
    # under Grain boundaries
    [
        (C_P2V, '-',  r'P2v: cavity growth  $v+V_m\to V_{m+1}$'),
        (C_P2H, '-',  r'P2h: He trapping  $c_h+V_{m,\ell}\to V_{m,\ell+1}$'),
        (C_P2I, '-',  r'P2i/P6: SIA absorption by cavity'),
    ],
    # under Precipitates
    [
        (C_P7,  '-',  r'P7: trap mutation  $V_{m,\ell}\to V_{m+1,\ell}+I_1$'),
        (C_P8,  '--', r'P8: radiation re-solution'),
        (C_P4,  '--', r'P4: absorption at fixed sinks  $D_\alpha c_\alpha$'),
    ],
]
col_x = [0.5, 10.6, 20.9]  # left edges, aligned with sink box columns

ax.text(col_x[0], KY, 'Process key',
        fontsize=FS_KEY_HDR, fontweight='bold', va='top')

ROW_H = 0.45
for col_idx, entries in enumerate(key_columns):
    KX = col_x[col_idx]
    for k, (color, ls, label) in enumerate(entries):
        y_row = KY - 0.55 - k * ROW_H
        lx = KX + 0.15
        if ls == '--':
            ax.plot([lx, lx+0.9], [y_row, y_row],
                    color=color, lw=2.2, linestyle='dashed')
        else:
            ax.annotate('', xy=(lx+0.9, y_row), xytext=(lx, y_row),
                        arrowprops=dict(arrowstyle='->', color=color, lw=2.2))
        ax.text(lx+1.15, y_row, label,
                va='center', fontsize=FS_KEY, color='black')


# ─────────────────────────────────────────────────────────────────────────────
# SAVE
# ─────────────────────────────────────────────────────────────────────────────
import os, shutil

module_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out_dir     = os.path.join(module_root, 'output', 'schematics')
os.makedirs(out_dir, exist_ok=True)

out_path = os.path.join(out_dir, 'eurofer_reaction_network.png')
old_path = os.path.join(out_dir, 'eurofer_reaction_network_old.png')

if os.path.exists(out_path):
    shutil.move(out_path, old_path)

plt.savefig(out_path, dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("PNG saved to", out_path)
