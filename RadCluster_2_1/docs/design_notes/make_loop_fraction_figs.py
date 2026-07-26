"""Generate the summary figure for loop_fraction_two_mechanisms.tex.

Run:  python docs/design_notes/make_loop_fraction_figs.py
Produces: docs/design_notes/loop_fraction_results.pdf
All data are measured values from the calibration campaign (see the .tex).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito colourblind-safe categorical order (validated: CVD dE 11.0 worst adjacent)
C_BLUE, C_VERM, C_GREEN, C_ORANGE, C_PURPLE = (
    "#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7")
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "figure.dpi": 150,
})

OUT = Path(__file__).parent

# ---------------------------------------------------------------- panel (a)
# f_100(T). filled marker = converged (reached full 3 dpa); open = dose-starved.
T18 = [350, 400, 450, 500, 550]
F18 = [0.0820, 0.1332, 0.1620, 0.3249, 0.7762]
K18 = [True, False, False, True, True]          # converged?
T16 = [350, 400, 450, 500]
F16 = [0.1579, 0.2767, 0.2562, 0.9483]
K16 = [True, False, False, True]

fig, axes = plt.subplots(1, 3, figsize=(11.6, 3.5))
ax = axes[0]
# experimental target band: significant by 400 C, ~complete by 500 C
ax.axvspan(400, 500, color=C_GREEN, alpha=0.07, lw=0)
ax.annotate("experimental\ntarget window", xy=(450, 0.055), ha="center",
            color=C_GREEN, fontsize=7.5)

for T, F, K, col, lab in ((T18, F18, K18, C_BLUE, r"$E_a^0=1.8$ eV"),
                          (T16, F16, K16, C_VERM, r"$E_a^0=1.6$ eV")):
    ax.plot(T, F, "-", color=col, lw=2, zorder=3)
    for t, f, k in zip(T, F, K):
        ax.plot(t, f, "o", ms=8, zorder=4,
                mfc=col if k else "white", mec=col, mew=1.8)
    ax.annotate(lab, xy=(T[-1], F[-1]), xytext=(4, 2),
                textcoords="offset points", color=col, fontsize=8.5,
                va="bottom", fontweight="bold")

ax.set_xlabel(r"irradiation temperature $T$ ($^\circ$C)")
ax.set_ylabel(r"$f_{\langle100\rangle}$  (content-weighted)")
ax.set_title("(a) Loop-character fraction vs temperature", loc="left")
ax.set_ylim(0, 1.05); ax.set_xlim(335, 585)
ax.plot([], [], "o", ms=7, mfc=MUTED, mec=MUTED, label="converged (3 dpa)")
ax.plot([], [], "o", ms=7, mfc="white", mec=MUTED, mew=1.6,
        label="dose-starved (transient)")
ax.legend(frameon=False, loc="upper left")

# ---------------------------------------------------------------- panel (b)
# peak-then-decay: f_100 against the growing mean 1/2<111> loop size.
tr = {
    r"400 $^\circ$C, $E_a^0$=1.6": (
        [5.3, 5.4, 5.7, 7.4, 8.6, 10.0, 13.9, 26.6, 60.4, 60.8],
        [0.0000, 0.0000, 0.0001, 0.0025, 0.0285, 0.1644, 0.3548, 0.3899,
         0.2776, 0.2767], C_VERM),
    r"450 $^\circ$C, $E_a^0$=1.6": (
        [5.4, 5.4, 6.3, 8.1, 9.2, 11.9, 18.4, 44.7, 91.9, 148.4],
        [0.0000, 0.0001, 0.0014, 0.0224, 0.1918, 0.5981, 0.7050, 0.5103,
         0.3456, 0.2562], C_ORANGE),
    r"500 $^\circ$C, Marian on": (
        [9.1, 11.9, 20.1, 50.1, 51.6, 51.6, 51.6],
        [0.1048, 0.4025, 0.5122, 0.3477, 0.3421, 0.3421, 0.3421], C_GREEN),
}
ax = axes[1]
for lab, (x, y, col) in tr.items():
    ax.plot(x, y, "-o", color=col, lw=2, ms=5, label=lab, zorder=3)
    ipk = int(np.argmax(y))
    ax.plot(x[ipk], y[ipk], "o", ms=10, mfc="none", mec=col, mew=2, zorder=4)
ax.annotate("peak", xy=(18.4, 0.705), xytext=(9, 3),
            textcoords="offset points", color=C_ORANGE, fontsize=8, ha="left")
ax.annotate("decay as loops outgrow\nthe unary window", xy=(95, 0.075),
            xytext=(0, 0), textcoords="offset points", color=MUTED,
            fontsize=7.5, ha="center", va="center")
ax.set_xscale("log")
ax.set_xlabel(r"mean $\frac{1}{2}\langle111\rangle$ loop size $\bar n_{111}$ (SIA)")
ax.set_ylabel(r"$f_{\langle100\rangle}$")
ax.set_title("(b) Transient: unary alone cannot hold large loops", loc="left")
ax.set_ylim(0, 0.8)
ax.legend(frameon=False, loc="upper left")

# ---------------------------------------------------------------- panel (c)
# Marian gate: magnitude vs temperature-sensitivity trade-off.
ax = axes[2]
TC = np.linspace(300, 600, 200)
kT = 8.617e-5 * (TC + 273.15)
for dH2, col in ((1.00, C_BLUE), (0.70, C_PURPLE), (0.55, C_GREEN),
                 (0.40, C_ORANGE)):
    P = 1.0 / (1.0 + np.exp((dH2 - 0.30) / kT))
    ax.plot(TC, P, "-", color=col, lw=2, zorder=3)
    ax.annotate(rf"$\Delta H_2$={dH2:.2f}", xy=(600, P[-1]), xytext=(3, 0),
                textcoords="offset points", color=col, fontsize=8, va="center")
ax.axhspan(1e-6, 1e-3, color=MUTED, alpha=0.07, lw=0)
ax.annotate("channel effectively OFF", xy=(430, 3e-5), color=MUTED,
            fontsize=7.5, ha="center")
ax.set_yscale("log")
ax.set_xlabel(r"irradiation temperature $T$ ($^\circ$C)")
ax.set_ylabel(r"$P_{\rm succ}(T)$  (Marian gate)")
ax.set_title("(c) Marian gate: magnitude vs. T-sensitivity", loc="left")
ax.set_xlim(300, 640); ax.set_ylim(1e-6, 1)

for a in axes:
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "loop_fraction_results.pdf", bbox_inches="tight")
fig.savefig(OUT / "loop_fraction_results.png", bbox_inches="tight", dpi=200)
print("wrote", OUT / "loop_fraction_results.pdf")
