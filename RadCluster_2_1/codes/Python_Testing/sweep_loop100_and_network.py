"""
sweep_loop100_and_network.py — exercise the two newest RadCluster_2_1 features
against experiment:

  (1) ½⟨111⟩ → ⟨100⟩ loop conversion — content-weighted f₁₁₁(T) crossover.
      Sweeps T ∈ {450, 500, 550} °C with the design-note-calibrated crossover
      T_star_conv_C = 340 °C, so the Dudarev unary driving force turns on and
      the large loops convert as T rises (⟨100⟩ becomes dominant at high T —
      the experimental high-T regime).

  (2) Loop → network-dislocation loss — dynamic ρ_net(dose).
      Enables LOOP_NETWORK_LOSS and drives ρ_net through run_adaptive's
      operator-split so it rises ~1 order of magnitude and saturates near the
      experimental EUROFER network density (~5×10¹⁴ m⁻²).

Run (from module root, needs a built solver.exe):
  python codes/Python_Testing/sweep_loop100_and_network.py [--quick] [--T 500]
"""
import argparse
import io
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]          # RadCluster_2_1/
sys.path.insert(0, str(ROOT))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from py_utils.simulation import RadClusterSimulation

# ── ρ_net channel calibration (target: rise from rho_d0 → ~5e14, saturate) ────
# The geometric capture switch P_ℓd turns on only when the elastic interaction
# range χ·d_loop reaches the network spacing L_ℓd = ρ_net^(-1/2).  For the few-nm
# loops here that needs χ ≳ 30 (at the realistic as-tempered ρ_net≈5e13–1e14,
# spacing ~100 nm); above ~30 it saturates.  With the switch on, w_c=b_111 (the
# physical capture width) and K_rec=0 give a monotonic rise that saturates at the
# rho_max ceiling — the classic irradiation network-saturation signature.
# Conversion crossover: T_star=250C makes the unary dF>0 cover the loop content
# (peaks ~n=200) by 500C → f_100≈1, while staying partial at 400C (cutoff ~n=98)
# → "significant by 400C, near-100% by 500C" (the requested experimental target).
T_STAR_CONV = 250.0
# E_a0 is THE efficiency lever: it is size-INDEPENDENT, so it speeds conversion at
# the small/mid sizes (n~10-50) where the loops actually live and where the
# growth-vs-conversion race is decided.  1.8 -> 1.6 lifts f_100(500C) 0.32 -> 0.95
# while leaving 350C low (0.16).  (gamma_a scales with P(n) and only moved large-n
# rates where no loops exist; T_star only moves the dF cutoff, which loops never
# reach — both were dead ends.  See docs/design_notes/loop_conversion_calibration.md.)
E_A0_CONV   = 1.6           # unary direct-rotation barrier [eV]
GAMMA_A_CONV = 0.02         # unary barrier size-slope [eV/segment]: lowered from the
                            # 0.03 default so n~200 loops convert by ~2 dpa at 500C
                            # (f_100->1) while staying partial at 400C (dF cutoff ~98).
DH2_CONV     = None         # Marian two-step gate dH2 [eV]; None = workbook (1.0 = OFF)
PHI_MAX_JUNC = None         # junction peak yield; None = workbook default (0.5)
RHO_D0      = 5.0e13        # initial network density [m^-2] (~as-tempered EUROFER)
RHO_MAX     = 5.0e14        # ceiling = experimental saturation plateau [m^-2]
# Capture width: the PHYSICAL value is w_c≈b_111≈2.5e-10 m, but the net climb
# velocity v_net is small (SIA/vacancy fluxes nearly cancel), so the physical
# channel moves ρ_net only ~0.1 %/dpa (saturates at ~160 dpa).  w_c is the
# designated amplification knob (loop_network_loss.tex; check_loop_network_loss
# amplifies it likewise).  50×b_111 makes ρ_net rise and saturate at the ceiling
# by ~3–6 dpa — a clear demonstration; precise w_c/K_rec are an offline calibration.
LOOP_NET_WC = 3.72e-8       # = 150×b_111 (amplified demo; reaches the 5e14 ceiling
                            # by ~6 dpa so ρ_net rises then saturates within 10 dpa)
LOOP_NET_KREC = 0.0         # no recovery → monotonic rise to the ceiling
LOOP_NET_CHI  = 50.0        # elastic interaction range (loop-diameters); ≥30 = on
LOOP_NET_NINC = 10          # incorporation onset size (sessile loops dominate)


def run_one(T_C, t_end, n_points, I=1000, V=1000, i_mobile=20, v_mobile=5,
            rtol=1e-5, network_loss=True, conversion=True):
    """One temperature: conversion ON + (optionally) network loss, via
    run_adaptive so ρ_net advances between segments. Returns the results dict.

    Uses the efficient full_system + bin_moment + Woodbury path (the original
    16-min config) so loops reach n≈200 — the size the network-loss channel
    needs to activate (the capture switch needs χ·d_loop ≳ network spacing)."""
    _s = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=V, solver_mode="full_system",
            equations="bin_moment", cascade="fission", C_floor=1e-25,
            he_kinetics="quasi_steady_state",
            i_mobile=i_mobile, v_mobile=v_mobile)
        re = sim.input_data.reactions
        re["T"]             = float(T_C) + 273.15
        re["rho_d"]         = RHO_D0
        re["T_star_conv_C"] = T_STAR_CONV
        re["E_a0_conv"]     = E_A0_CONV
        re["gamma_a_conv"]  = GAMMA_A_CONV
        re["n_loop_min"]    = 4
        # bin-moment discrete/bin split (matches the original I=1000 run)
        re["i_discrete"] = 50; re["v_discrete"] = 20
        re["I_bin"] = 20;      re["V_bin"] = 20
        re["shape_function"] = "linear"
        if network_loss:
            re["LOOP_NETWORK_LOSS"] = 1
            re["loop_net_chi"]      = LOOP_NET_CHI
            re["loop_net_w_c"]      = LOOP_NET_WC
            re["loop_net_K_rec"]    = LOOP_NET_KREC
            re["loop_net_n_inc"]    = LOOP_NET_NINC
            re["loop_net_rho_max"]  = RHO_MAX
        sim.input_data._calculate_derived()
        sim.rebuild_rates()
        cfg = {"t_span": (1e-6, t_end), "n_points": n_points, "log_time": True,
               "rtol": rtol, "atol": 1e-20, "timeout_s": 2400,
               "solver_method": {"linsol": "gmres", "preconditioner": "Woodbury",
                                 "window_width": 10, "window_pad": 20,
                                 "concentration_threshold": 1e-22},
               "loop_conversion": int(conversion)}
        r = sim.run_adaptive(solver_config=cfg, save_output=False,
                             points_per_segment=max(4, n_points // 10),
                             max_doublings=2, boundary_threshold=0.05)
    finally:
        sys.stdout, sys.stderr = _s
    return r


def run_conv(T_C, t_end=3e6, n_points=60, I=200, i_mobile=10, v_mobile=3,
             rtol=3e-3, timeout_s=4500, progress_path=None):
    """f_111(T) CROSSOVER point: loop conversion ON, network loss OFF, on a
    reduced active_window/discrete grid so the (stiff) conversion runs complete.
    f_100 is a growth-vs-conversion competition gated by the dF cutoff (T_star):
    large loops convert fast at high T (gamma_a tuned) -> f_100->1, but stay
    partial at low T where the dF cutoff caps the convertible size.

    If ``progress_path`` is given, per-SEGMENT f_100 / loop concentrations are
    appended live (via a wrapper on sim._merge_results) so an in-flight run's
    progress is visible without waiting for the point to finish."""
    _s = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=I, solver_mode="active_window",
            physics_option="full_CD_fission", C_floor=1e-25,
            he_kinetics="quasi_steady_state",
            i_mobile=i_mobile, v_mobile=v_mobile)
        re = sim.input_data.reactions
        re["T"] = float(T_C) + 273.15
        re["T_star_conv_C"] = T_STAR_CONV
        re["E_a0_conv"] = E_A0_CONV
        re["gamma_a_conv"] = GAMMA_A_CONV
        re["n_loop_min"] = 4
        # Marian kinetic channels (junction + absorption). Their T-gate is
        # P_success = 1/(1+exp((dH2-dH_rev)/kT)); at the shipped dH2=1.0 it is
        # ~1e-5, i.e. the channels are OFF and the unary does all the work
        # (which is why f_100 peaks then decays once loops outgrow the unary
        # window). Lower dH2 to engage them on the LARGE loops.
        if DH2_CONV is not None:
            re["dH2_conv"] = DH2_CONV
        if PHI_MAX_JUNC is not None:
            re["phi_max_junc"] = PHI_MAX_JUNC
        sim.input_data._calculate_derived(); sim.rebuild_rates()

        # Live per-segment progress: wrap _merge_results (called once per
        # completed segment with post-processed results) to append a log line.
        if progress_path is not None:
            _orig_merge = sim._merge_results
            def _logged_merge(acc, res, _o=_orig_merge, _T=T_C, _p=progress_path):
                try:
                    def last(k, d=0.0):
                        a = res.get(k)
                        return float(np.asarray(a)[-1]) if a is not None and len(a) else d
                    with open(_p, "a", encoding="utf-8") as fh:
                        fh.write(f"  T={_T:3.0f}C  dose={last('dose'):6.3f}dpa  "
                                 f"f_100={1-last('f_111_loop',1.0):.4f}  "
                                 f"N_111={last('N_loops_111'):.3e}  "
                                 f"N_100={last('N_loops_100'):.3e}  "
                                 f"mean_n111={last('mean_n_111'):.1f}  "
                                 f"mean_n100={last('mean_n_100'):.1f}\n")
                except Exception:
                    pass
                return _o(acc, res)
            sim._merge_results = _logged_merge

        cfg = {"t_span": (1e-6, t_end), "n_points": n_points, "log_time": True,
               "rtol": rtol, "atol": 1e-20, "timeout_s": timeout_s,
               "solver_method": {"linsol": "dense"},
               "loop_conversion": 1}
        r = sim.run_adaptive(solver_config=cfg, save_output=False,
                             points_per_segment=5,
                             max_doublings=2, boundary_threshold=0.05)
    finally:
        sys.stdout, sys.stderr = _s
    return r


def summarize(T_C, r):
    f111 = float(np.asarray(r["f_111_loop"])[-1])
    rho0 = float(np.asarray(r["rho_net"])[0])
    rho1 = float(np.asarray(r["rho_net"])[-1])
    dose = float(np.asarray(r["dose"])[-1])
    print(f"  T={T_C:3.0f}C dose={dose:5.2f}dpa | "
          f"f_111={f111:.4f} f_100={1-f111:.4f} | "
          f"rho_net {rho0:.2e} -> {rho1:.2e}  (x{rho1/rho0:.1f}) "
          f"mean_n_100={float(np.asarray(r['mean_n_100'])[-1]):.1f}")
    return dict(T_C=T_C, dose=dose, f111=f111, rho0=rho0, rho1=rho1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="fast probe: low dose, smaller domain")
    ap.add_argument("--T", type=float, default=None,
                    help="single temperature (°C) instead of the 3-point sweep")
    ap.add_argument("--no-conv", action="store_true",
                    help="disable loop conversion (isolate/verify the ρ_net channel)")
    ap.add_argument("--crossover", action="store_true",
                    help="f_111(T) crossover sweep: conversion ON, network OFF, "
                         "reduced grid (350/400/450/500/550 C)")
    ap.add_argument("--progress", type=str, default=None,
                    help="path to a live per-segment progress log file")
    ap.add_argument("--gamma-a", type=float, default=None,
                    help="override GAMMA_A_CONV (unary barrier size-slope)")
    ap.add_argument("--e-a0", type=float, default=None,
                    help="override E_A0_CONV (unary barrier offset)")
    ap.add_argument("--t-end", type=float, default=None,
                    help="override crossover t_end [s] (dose = t_end*1e-6 dpa)")
    ap.add_argument("--t-star", type=float, default=None,
                    help="override T_STAR_CONV (conversion crossover T, C)")
    ap.add_argument("--timeout", type=float, default=None,
                    help="override per-point solver timeout [s]")
    ap.add_argument("--dh2", type=float, default=None,
                    help="override dH2_conv [eV] — the Marian two-step gate. "
                         "Lower => P_success up => junction+absorption engage.")
    ap.add_argument("--phi-max", type=float, default=None,
                    help="override phi_max_junc (junction peak yield)")
    args = ap.parse_args()
    if args.gamma_a is not None:
        GAMMA_A_CONV = args.gamma_a
    if args.e_a0 is not None:
        E_A0_CONV = args.e_a0
    if args.t_star is not None:
        T_STAR_CONV = args.t_star
    if args.dh2 is not None:
        DH2_CONV = args.dh2
    if args.phi_max is not None:
        PHI_MAX_JUNC = args.phi_max

    if args.crossover:
        temps = ([args.T] if args.T is not None
                 else [350.0, 400.0, 450.0, 500.0, 550.0])
        print(f"f_111(T) crossover sweep: conversion ON, T_star={T_STAR_CONV}C, "
              f"E_a0={E_A0_CONV}eV, gamma_a={GAMMA_A_CONV}, network OFF, T={temps} C",
              flush=True)
        rows = []
        import time
        for T_C in temps:
            if args.progress:
                with open(args.progress, "a", encoding="utf-8") as fh:
                    fh.write(f"\n=== T={T_C:.0f}C starting ===\n")
            t0 = time.perf_counter()
            _kw = {}
            if args.t_end:
                _kw["t_end"] = args.t_end
            if args.timeout:
                _kw["timeout_s"] = args.timeout
            r = run_conv(T_C, progress_path=args.progress, **_kw)
            dt = time.perf_counter() - t0
            if r is None:
                print(f"  T={T_C:3.0f}C  FAILED", flush=True); continue
            f111 = float(np.asarray(r["f_111_loop"])[-1])
            dose = float(np.asarray(r["dose"])[-1])
            print(f"  T={T_C:3.0f}C dose={dose:5.2f}dpa | f_100={1-f111:.4f} | "
                  f"mean_n_100={float(np.asarray(r['mean_n_100'])[-1]):.1f} "
                  f"mean_n_111={float(np.asarray(r['mean_n_111'])[-1]):.1f} (wall {dt:.0f}s)",
                  flush=True)
            rows.append((T_C, 1 - f111))
        print("\nf_111(T) crossover   T[C]   f_100")
        for T_C, f100 in rows:
            print(f"  {T_C:5.0f}   {f100:.4f}")
        sys.exit(0)

    # Efficient full_system + bin_moment + Woodbury config (loops reach n≈200).
    if args.quick:
        t_end, n_points, I = 3e6, 60, 1000     # ~3 dpa verification
    else:
        t_end, n_points, I = 1e7, 80, 1000     # ~10 dpa

    temps = [args.T] if args.T is not None else [450.0, 500.0, 550.0]
    print(f"Sweep T={temps} °C  (t_end={t_end:.0e}s ~ {t_end*1e-6:.0f} dpa, "
          f"I={I} bin_moment/full_system, T_star=340C, LOOP_NETWORK_LOSS on)")
    rows = []
    import time
    for T_C in temps:
        t0 = time.perf_counter()
        r = run_one(T_C, t_end, n_points, I=I, V=I, conversion=not args.no_conv)
        dt = time.perf_counter() - t0
        if r is None:
            print(f"  T={T_C:3.0f}C  FAILED (no results)")
            continue
        row = summarize(T_C, r)
        row["wall_s"] = dt
        print(f"       (wall {dt:.0f}s)")
        rows.append(row)

    print("\nSummary  T[C]  f_100   rho_net_final")
    for x in rows:
        print(f"  {x['T_C']:5.0f}  {1-x['f111']:.4f}  {x['rho1']:.2e}")
