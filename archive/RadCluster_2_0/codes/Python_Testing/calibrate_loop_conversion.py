"""
calibrate_loop_conversion.py — Phase-6 calibration of the ½⟨111⟩→⟨100⟩
conversion against the experimental ½⟨111⟩-loop fraction f₁₁₁(T).

Strategy
--------
The crossover temperature T* (where ⟨100⟩ overtakes ½⟨111⟩) is the dominant
knob; the other five (E_a0, γ_a, φ_max, σ_s, n_j_min) are anchored by Marian /
Dudarev and held at their physics defaults.  We:

  1. load the experimental f₁₁₁(T) for EUROFER-97 (the target material) from
     the microstructure database via loop_burgers_fraction;
  2. run the C++ solver (conversion ON) over a temperature sweep for each
     candidate T*, extracting the modeled f₁₁₁ = S_I^{111} / (S_I^{111}+S_I^{100})
     at a fixed dose;
  3. pick the T* minimising the squared error to the experimental points;
  4. write a model-vs-experiment plot and report the calibrated T*.

This is intentionally a 1-D calibration of the dominant knob: the data are
scattered across materials/ion-vs-neutron, so fitting all six knobs would
overfit.  Run cost is a few dozen short solver runs.

Usage:
  python codes/Python_Testing/calibrate_loop_conversion.py [--quick]
"""
import io
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "py_utils"))

from py_utils.simulation import RadClusterSimulation
from loop_burgers_fraction import extract_f111_table

DB = ROOT / "input" / "Ferritics_Irradiated_Microstructure_Data.xlsx"
G_DPA = 1.0e-6          # displacement rate [dpa/s] (default)


# ── experimental target ───────────────────────────────────────────────────────
def f111_experimental(material_filter="EUROFER"):
    """Return (T_C, f111) experimental points for the target material."""
    df = extract_f111_table(str(DB))
    sub = df[df["material"].str.contains(material_filter, case=False, na=False)]
    sub = sub.dropna(subset=["temperature_C", "f_111"])
    return sub["temperature_C"].to_numpy(), sub["f_111"].to_numpy()


# ── model evaluation ──────────────────────────────────────────────────────────
def f111_model(T_C, t_end_s, overrides=None, I=60):
    """Run the solver at temperature T_C to t_end_s; return modeled f₁₁₁."""
    _s = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=I, solver_mode="full_system",
            physics_option="full_CD_fission", C_floor=1e-25,
            he_kinetics="quasi_steady_state", i_mobile=5, v_mobile=2)
        inp = sim.input_data
        inp.reactions["T"] = T_C + 273.15
        for k, v in (overrides or {}).items():
            inp.reactions[k] = v
        inp._calculate_derived()
        sim.rebuild_rates()
        cfg = {"t_span": (1e-6, t_end_s), "n_points": 6, "log_time": True,
               "rtol": 1e-6, "atol": 1e-26,
               "solver_method": {"linsol": "dense"}, "loop_conversion": 1}
        r = sim.run_adaptive(solver_config=cfg, save_output=False)
    finally:
        sys.stdout, sys.stderr = _s
    y = r["y"]; y100 = r["y_sia100"]; sz = np.arange(1, I + 1)
    S111 = float(np.dot(sz, np.maximum(y[:I, -1], 0.0)))
    S100 = float(np.dot(sz, np.maximum(y100[:, -1], 0.0)))
    return S111 / (S111 + S100) if (S111 + S100) > 0 else 1.0


# ── calibration sweep ─────────────────────────────────────────────────────────
def main():
    quick = "--quick" in sys.argv
    Tgrid = np.array([250, 300, 350, 400, 450] if quick
                     else [250, 300, 330, 360, 400, 450, 500], dtype=float)
    Tstar_candidates = [400.0] if quick else [300.0, 350.0, 400.0, 450.0, 500.0]
    dose = 0.3                      # dpa target for the sweep (conversion saturates early)
    t_end = dose / G_DPA            # [s]

    Texp, fexp = f111_experimental("EUROFER")
    print(f"Experimental EUROFER points: {len(Texp)}  "
          f"(T {Texp.min():.0f}-{Texp.max():.0f} C)")

    # interpolate experiment onto the model grid for the loss
    order = np.argsort(Texp)
    fexp_on_grid = np.interp(Tgrid, Texp[order], fexp[order])

    results = {}   # Tstar -> modeled f111(Tgrid)
    best = (None, np.inf, None)
    for Tstar in Tstar_candidates:
        curve = np.array([
            f111_model(T, t_end, overrides={"T_star_conv_C": Tstar})
            for T in Tgrid])
        loss = float(np.mean((curve - fexp_on_grid) ** 2))
        results[Tstar] = curve
        print(f"  T*={Tstar:5.0f} C  ->  f111(T)="
              f"[{', '.join('%.2f' % c for c in curve)}]  loss={loss:.4f}")
        if loss < best[1]:
            best = (Tstar, loss, curve)

    print(f"\nBEST T* = {best[0]:.0f} C   (loss={best[1]:.4f})")

    # ── plot ──────────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(Texp, 100 * fexp, s=80, c="k", marker="o",
                   edgecolors="white", zorder=5, label="EUROFER-97 (exp.)")
        for Tstar, curve in results.items():
            lw = 2.5 if Tstar == best[0] else 1.2
            ax.plot(Tgrid, 100 * curve, "-o", lw=lw,
                    label=f"model T*={Tstar:.0f} C"
                          + (" (best)" if Tstar == best[0] else ""))
        ax.set_xlabel("Irradiation temperature (C)")
        ax.set_ylabel("1/2<111> loop fraction (%)")
        ax.set_title(f"Loop-conversion calibration ({dose:.0f} dpa)")
        ax.set_ylim(-5, 105); ax.grid(True, ls=":", alpha=0.4)
        ax.legend(fontsize=9, frameon=False)
        out = ROOT / "output" / "loop_conversion_calibration.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout(); fig.savefig(out, dpi=150)
        print(f"Saved plot: {out}")
    except Exception as e:
        print(f"(plot skipped: {e})")

    return best[0]


if __name__ == "__main__":
    main()
