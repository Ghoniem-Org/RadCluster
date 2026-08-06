#!/usr/bin/env python
"""run_ensemble - the campaign worker.  One instance per machine.

    machine k of M  ->  runs every design row with  row_id % M == k

Deterministic, restartable, needs no shared filesystem and no coordination: a
dead machine leaves an identifiable GAP rather than corrupting the design.
Results are appended one JSON object per line, so a crash costs one row rather
than the file, and merging is plain concatenation.

    python run_ensemble.py --design design/T2_design_v1.csv \
                           --machine 0 --of 4 --workers 10

WHY THE ADMISSIBILITY BLOCK EXISTS
----------------------------------
delta_FP is BLIND to grid truncation.  On 2026-08-01 it held at 1e-8 while
99.96 % of the <100> content was stacked against the top of the grid, because
the `if (n < I)` growth guard in rate_kernels.cpp halts growth WITHOUT losing
atoms.  Conservation is a health check, not a grid-adequacy check.

Across ~1250 rows with theta varying over its full prior box, some cells WILL
saturate the grid.  Without pile/occupancy recorded per row those cells enter
the Sobol estimator as ordinary data and their sensitivity is an artefact of
the grid, not of the physics.  So every row carries:

    occ_111, occ_100     mean size / I          (project rule: > ~0.1 suspect)
    pile_111, pile_100   content fraction in the top 2 % of the grid
    pile_v, occ_v        the SAME test on the VACANCY axis (added 2026-08-02)
    d_over_ceiling_100   d_100 / d(n=I)
    dose_reached, starved, delta_FP, delta_He, solver_rc, wall_s

The vacancy axis went unchecked until 2026-08-02 and is not a minor addition:
at V=120 the check_machine probe reported mean_n_v = 18.7 against a converged
153.3 and N_voids 10.6x low.  Note the asymmetry -- delta_FP is blind to
SIA-axis truncation but IS sensitive to the vacancy axis, because the swelling
identity S = S_I + Delta J^d has the vacancy inventory on the left.  So on the
vacancy side delta_FP is a real guard, and the admissibility bar was tightened
from 1e-3 to 1e-6 to make it bite.

and `admissible`, which merge_and_sobol uses to decide what may be estimated
from.  A row that is inadmissible is NOT a failed row -- it ran fine, it just
cannot answer the question asked of it.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
import traceback
from pathlib import Path

# Force, do NOT setdefault: an inherited OMP_NUM_THREADS (the
# workstation exports 24) would give every worker its own thread
# pool -- 10 workers x 24 threads on 24 cores -- and would also
# let reduction order vary between machines.  Many serial workers
# beat few threaded ones for an ensemble anyway.
os.environ["OMP_NUM_THREADS"] = "1"   # many serial jobs > few threaded

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO))

STOP_FILE = HERE / "CAMPAIGN_STOP"   # graceful-halt sentinel (campaign_ops)
B111_M = 2.482e-10       # <111> Burgers magnitude [m], for loop_net_w_c units
PILE_TOP_FRAC = 0.02     # "top 2 % of the grid"
PILE_TOL = 0.05          # above this the size readout is a ceiling artefact
OCC_TOL = 0.10           # RECORDED ONLY since 2026-08-03 -- no longer a gate.
                         # Probe: row 24 at occ_100=0.171 / pile_100=2.3e-10 is
                         # bit-identical at I=3200 and I=12800.  Occupancy is not
                         # a truncation test (plan S11(c)-2).
TOPBIN_TOL = 0.02        # bin_moment: max content fraction in the TOP BIN.
                         # `pile` is measured on the reconstructed per-size
                         # distribution and inherits the closure's smoothing;
                         # this reads the raw moments instead.

# CLAUDE.md S8's own standard: delta_FP below ~1e-6 is good, above 1e-3 signals
# a coding error.  Admissibility used 1e-3, i.e. it only rejected the "coding
# error" band.  That is too loose for the VACANCY axis, where delta_FP is the
# primary guard: the 2026-08-02 vacancy study passed V=240 at delta_FP=8.3e-4
# (under the old bar) while mean_n_v was still 11 % from converged, and reached
# 1.4e-7 at V=480 where it had converged.  1e-6 separates those two cleanly and
# is comfortably met by good configs (production rows run 1.5e-8 to 2.9e-12).
#
# 2026-08-03: 1e-6 -> 1e-2 (plan S11(h)).  Three reasons the tightening above
# was wrong: (i) 1e-6 IS the solver's own rtol, so it rejected for solver noise
# -- only 12 % of grid-converged v2 rows cleared it; (ii) it was tightened to
# catch VACANCY truncation, which V=10000 now designs out, and at 1e-2 it still
# rejects the known-bad case by 53x (row 229 at V=600 sits at 0.530, dropping to
# 7.1e-6 at V=2400); (iii) bin_moment carries a closure error discrete does not,
# so a bar set from discrete experience rejects the method rather than the run.
#
# THIS VALUE IS ADVISORY.  It sets the worker's own `admissible` field, which
# feeds only the live monitor.  The gate that decides what enters a Sobol index
# is applied post-hoc in merge_and_sobol.usable(), from the per-row quantities
# recorded here -- so a threshold can be revised without re-running anything.
# That is not a nicety: re-scoring the stopped v2 campaign under the corrected
# rules moved it from 12 admissible rows to 52 at zero CPU cost.
DFP_TOL = 1e-2

# ---------------------------------------------------------------------------
# TRUNCATION IS NO LONGER A GATE -- author decision, 2026-08-05 (plan S11(q)).
#
# The campaign's question is the RELATIVE RANKING of parameters, not the
# absolute value of any observable, and a ranking survives a truncated tail.
# Every truncation quantity below is still MEASURED AND RECORDED on every row;
# what changed is that none of them may reject a row.
#
# delta_FP is included in this withdrawal, and that is the non-obvious part.
# Under `discrete` it is an honest conservation check.  Under `bin_moment` it is
# not independent of truncation: the top bin's closure cannot absorb overflow,
# so delta_FP tracks <100> truncation almost exactly -- 1.71e-04 on the
# converged I=50000 row against 0.907 on the truncated I=10000 one, the same row
# and the same theta.  Leaving it as a gate at 1e-2 would go on rejecting
# precisely the truncated rows this decision keeps, under a different name.
#
# What still rejects a row: the solver failed outright.  Everything else is a
# recorded covariate.
TRUNCATION_GATES = False

# ---------------------------------------------------------------------------
# A TIMED-OUT ROW IS KEPT -- author decision, 2026-08-05 (plan S12(s)).
#
# `--timeout-s` does not kill the solver: cpp_bridge sends a graceful interrupt
# and the solver flushes the trajectory it has, so a row that runs out of time
# still returns everything up to the dose it reached.  Discarding that was
# throwing away finished work -- 57 of the 516 pooled v2 rows were dropped for
# this alone.
#
# It is NOT, however, sound to read such a row at whatever dose it happened to
# reach: dose moves the observables far harder than any numerical choice we have
# measured, and how far a row gets is a function of theta.  So the observables
# are recorded on a fixed DOSE LADDER (below) and screened at a common dose.
# `starved` stays recorded, and stops rejecting.
STARVED_GATE = False

# Dose ladder for the rev-6 1 dpa campaign.  Must be ascending and should end at
# the campaign dose.  A row contributes to every rung it reached.
# Low rungs added 2026-08-06.  Dose accrues very late in wall-clock terms: a
# 300 s cut of a row that reaches 1.0 dpa in 1700 s sat at 0.003 dpa and filled
# ZERO of the old 0.1-1.0 rungs, so the partial trajectory the deadline exists
# to preserve contributed nothing. Rungs are free -- the trajectory is already
# computed -- so cover the decades a cut row can actually land in.
DOSE_CHECKPOINTS = (0.005, 0.01, 0.02, 0.05,
                    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


# ------------------------------------------------------------------- utilities
def sha256_file(p: Path) -> str:
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except Exception:
        return "unavailable"


def git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(REPO),
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return "unknown"


def read_design(path: Path) -> tuple[list[dict], dict]:
    meta_p = path.with_suffix(".meta.json")
    meta = json.loads(meta_p.read_text(encoding="utf-8")) if meta_p.exists() else {}
    got = sha256_file(path)
    if meta.get("design_sha256") and meta["design_sha256"] != got:
        raise SystemExit(
            f"DESIGN HASH MISMATCH\n  expected {meta['design_sha256']}\n"
            f"  got      {got}\nThe design file has been edited since it was "
            f"generated. Saltelli pairing cannot be trusted; regenerate or "
            f"restore it. Refusing to run.")
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    cols = lines[0].split(",")
    rows = []
    for ln in lines[1:]:
        vals = ln.split(",")
        d = {}
        for c, v in zip(cols, vals):
            if c in ("row_id", "cond_row_id", "base_idx", "param_j"):
                d[c] = int(v)
            elif c in ("matrix", "condition"):
                d[c] = v
            else:
                d[c] = float(v)
        rows.append(d)
    return rows, meta


# ------------------------------------------------------- theta -> InputData
def apply_theta(sim, spec: dict, row: dict, cond: dict):
    """Write one design row into the simulation's InputData.

    Returns the dict of values actually written, for the provenance record.
    Constructor-level parameters (i_mobile, v_mobile) are handled by the caller
    BEFORE construction -- they cannot be set here.
    """
    d = sim.input_data
    by_key = {p["key"]: p for p in spec["parameters"]}
    written = {}

    def put(sheet, key, val):
        tgt = getattr(d, sheet, None)
        if tgt is None:
            raise KeyError(f"InputData has no sheet {sheet!r}")
        tgt[key] = val
        written[f"{sheet}.{key}"] = val

    for key, val in row.items():
        if key not in by_key:
            continue
        p = by_key[key]
        sheets = p["sheet"]
        if sheets in ("ctor", "postproc"):
            continue                      # handled elsewhere / post-processing
        if isinstance(sheets, str):
            sheets = [sheets]
        v = val
        if key == "loop_net_w_c":
            v = val * B111_M              # spec is in units of b_111
        for sh in sheets:
            put(sh, key, v)

        # derived partners: the <100> binding law is a FIXED offset from <111>
        # so the character ratio cannot absorb the f_100 data (plan S2.1).
        if key == "B_111":
            put("dissociation", "B_100", val * 0.9246)

    # A_111 is DERIVED from E_b_i2 (plan S2.2c), never sampled directly.
    if "E_b_i2" in row and "B_111" in row:
        # E_b_loop(2) = A_111 * 2^-B_111  =>  A_111 = E_b_i2 * 2^B_111
        a111 = float(row["E_b_i2"]) * (2.0 ** float(row["B_111"]))
        put("dissociation", "A_111", a111)
        put("dissociation", "A_100", a111 * 0.9545)

    # fixed values, written explicitly so a stale workbook cannot leak in
    for f in spec.get("fixed", []):
        k, v = f["key"], f["value"]
        if k in ("psucc_abs_pref", "dH_rev_conv", "gamma_a_conv", "n_j_min_junc",
                 "n_j_min_frac", "T_star_conv_C", "nu0_conv", "loop_net_rho_max",
                 "n_ref_conv", "absorb_boost_100", "grow_boost_100"):
            put("reactions", k, v)

    # condition: the ACTUAL T and dose rate of the experiment.
    # T_eq = T - 60 C is a data-side plotting device (plan S3.2) and must never
    # reach the simulator -- asserted here rather than trusted.
    assert "T_eq" not in cond, "T_eq must not be passed to the simulator"
    put("reactions", "T", float(cond["T_K"]))
    # null in conditions.json means "keep the workbook value" -- test for None,
    # not for key presence, or the null is coerced and blows up.
    if cond.get("G") is not None:
        put("reactions", "G", float(cond["G"]))
    if cond.get("rho_d") is not None:
        put("reactions", "rho_d", float(cond["rho_d"]))
    return written


# --------------------------------------------------------- bin-moment config
def apply_bin_config(sim, cfg):
    """Write the bin-moment layout onto the workbook dict.

    MUST be called after the constructor and before _calculate_derived() +
    rebuild_rates().  Bin configuration is NOT a constructor argument -- passing
    i_discrete=/I_bin=/shape_function= to RadClusterSimulation() raises
    TypeError since 2026-08-02, and before that it was silently dropped, which
    produced nine bit-identical "refinement" runs (plan S10(i)).

    No-op in discrete mode, so it is safe to call unconditionally.
    """
    if cfg["equations"] != "bin_moment":
        return
    r = sim.input_data.reactions
    r["i_discrete"]    = int(cfg["i_discrete"])
    r["v_discrete"]    = int(cfg["v_discrete"])
    r["I_bin"]         = int(cfg["I_bin"])
    r["V_bin"]         = int(cfg["V_bin"])
    r["shape_function"] = str(cfg["shape_function"])


def bin_layout(sim, cfg):
    """Read the REALISED bin layout back off sim.rate_equations and verify it.

    Two reasons this is not bookkeeping:

    1.  The layout must be read from `rate_equations` (the BinMomentRateEquations
        that owns i_discrete/I_bin/r/shape_function), NOT `reaction_rates` --
        querying the wrong object returns nan rather than raising (plan S10(i)).
    2.  The realised bin count is NOT the requested one.  `r` is derived as
        (I/i_discrete)**(1/I_bin) and build_bins() then walks floor(edge*r) with
        a minimum increment of 1, so the walk overshoots or undershoots the
        target by a bin or two.  Recording the request and silently running a
        different layout is exactly the class of error S10(i) is about.

    Raises if the request did not take effect at all -- a silent fallback to the
    coarse default (i_discrete=10, I_bin=6) is the failure mode that must not
    reach a production row.
    """
    out = {}
    if cfg["equations"] != "bin_moment":
        return out
    re_obj = getattr(sim, "rate_equations", None)
    if re_obj is None:
        raise RuntimeError("bin_moment requested but sim.rate_equations is None")
    got_id = int(getattr(re_obj, "i_discrete", -1))
    got_ib = int(getattr(re_obj, "I_bin", -1))
    got_vb = int(getattr(re_obj, "V_bin", -1))
    got_sf = str(getattr(re_obj, "shape_function", "?"))
    out["bin_i_discrete"] = got_id
    out["bin_v_discrete"] = int(getattr(re_obj, "v_discrete", -1))
    out["bin_I_bin"] = got_ib
    out["bin_V_bin"] = got_vb
    out["bin_r"] = float(getattr(re_obj, "r", float("nan")))
    out["bin_shape"] = got_sf
    out["bin_n_mom"] = int(getattr(re_obj, "n_mom", -1))
    out["N_eq"] = int(getattr(re_obj, "N_eq", -1))
    # The request must have LANDED.  i_discrete and shape_function are set
    # verbatim, so they must match exactly; the bin counts are allowed to drift
    # by the floor() walk but not to collapse to the default.
    if got_id != int(cfg["i_discrete"]) or got_sf != str(cfg["shape_function"]):
        raise RuntimeError(
            f"bin config did not take effect: asked i_discrete="
            f"{cfg['i_discrete']} shape={cfg['shape_function']}, "
            f"rate_equations reports i_discrete={got_id} shape={got_sf}. "
            f"This is the silent-default failure of plan S10(i).")
    for name, want, got in (("I_bin", int(cfg["I_bin"]), got_ib),
                            ("V_bin", int(cfg["V_bin"]), got_vb)):
        if want > 0 and abs(got - want) > max(3, 0.2 * want):
            raise RuntimeError(
                f"bin config {name}: asked {want}, realised {got} -- too far "
                f"off to be the floor() walk. Check i_discrete < I.")
    return out


# ------------------------------------------------------------- one evaluation
def evaluate(args):
    row, spec, cond, cfg = args
    t0 = time.time()
    rec = {"row_id": int(row["row_id"]), "cond_row_id": int(row["cond_row_id"]),
           "matrix": row["matrix"], "base_idx": int(row["base_idx"]),
           "param_j": int(row["param_j"]), "condition": row["condition"]}
    # theta hash: catches design drift without storing 20 floats per row
    theta_keys = sorted(k for k in row
                        if k not in ("row_id", "cond_row_id", "matrix",
                                     "base_idx", "param_j", "condition"))
    rec["theta_hash"] = hashlib.sha256(
        json.dumps({k: row[k] for k in theta_keys}, sort_keys=True).encode()
    ).hexdigest()[:16]

    try:
        from RadCluster_2_1.py_utils.simulation import RadClusterSimulation
        i_mob = int(row.get("i_mobile", cfg["i_mobile_default"]))
        v_mob = int(row.get("v_mobile", cfg["v_mobile_default"]))
        _s = sys.stdout, sys.stderr
        buf = io.StringIO()
        try:
            sys.stdout = sys.stderr = buf
            sim = RadClusterSimulation(
                I=cfg["I"], V=cfg["V"], solver_mode=cfg["solver_mode"],
                equations=cfg["equations"], cascade=cond.get("cascade", "fission"),
                C_floor=cfg["C_floor"], he_kinetics="quasi_steady_state",
                i_mobile=i_mob, v_mobile=v_mob)
            written = apply_theta(sim, spec, row, cond)
            apply_bin_config(sim, cfg)
            sim.input_data._calculate_derived()
            sim.rebuild_rates()
            rec.update(bin_layout(sim, cfg))
            G = float(sim.input_data.reactions["G"])
            scfg = {"t_span": (1e-6, cfg["dose"] / G), "n_points": cfg["n_points"],
                    "log_time": True, "rtol": cfg["rtol"], "atol": cfg["atol"],
                    "timeout_s": cfg["timeout_s"],
                    "solver_method": {"linsol": "gmres",
                                      "preconditioner": cfg["preconditioner"],
                                      "concentration_threshold": 1e-22},
                    "loop_conversion": 1}
            res = sim.run(solver_config=scfg, save_output=False)
        finally:
            sys.stdout, sys.stderr = _s
        rec["solver_rc"] = 0
        rec.update(observe(res, sim, cfg, float(row.get("d_min_tem", 1.0))))
        rec["n_written"] = len(written)
    except Exception as exc:
        rec["solver_rc"] = 1
        rec["error"] = f"{type(exc).__name__}: {exc}"[:300]
        rec["traceback_tail"] = traceback.format_exc()[-500:]
        rec["admissible"] = False
    rec["wall_s"] = round(time.time() - t0, 1)
    return rec


def per_size_populations(res, sim, cfg):
    """Return (c111, c_v, binned, topbin) at the FINAL time, per size.

    One extraction that works in both modes, so every downstream gate is
    computed the same way and is comparable across the discrete/bin_moment
    bridge (plan S11(e)-2).

      discrete    : `y` is (I + V + extras, nt) -- plain slices.
      bin_moment  : `y` is (N_eq, nt) with the layout
                    [0 : i_d]                       discrete SIA
                    [i_d : i_d + P*I_bin]           SIA bin moments, P per bin
                    [i_VAC : i_VAC + v_d]           discrete vacancies
                    [i_VAC+v_d : ... + P*V_bin]     vacancy bin moments
                    mirroring RadClusterSimulation._size_distributions.

    `topbin` maps axis -> content fraction in the top bin, computed from the RAW
    mu1 moments (not the reconstruction).  Empty dict in discrete mode.

    Returns c111 / c_v as None when an axis genuinely cannot be extracted; the
    caller treats None as FAILING the gate, never as passing.
    """
    I, V = int(cfg["I"]), int(cfg["V"])
    y = res.get("y")
    topbin = {}
    if y is None:
        return None, None, False, topbin
    y = np.asarray(y)
    if y.ndim != 2:
        return None, None, False, topbin
    yj = np.maximum(y[:, -1], 0.0)

    re_obj = getattr(sim, "rate_equations", None)
    binned = cfg.get("equations") == "bin_moment" and hasattr(re_obj, "bins")

    if not binned:
        c111 = yj[:I] if yj.shape[0] >= I else None
        c_v = yj[I:I + V] if yj.shape[0] >= I + V else None
        return c111, c_v, False, topbin

    from RadCluster_2_1.py_utils.bin_moment_rates import reconstruct_distribution
    P = int(getattr(re_obj, "n_mom", 2))
    sf = str(getattr(re_obj, "shape_function", "linear"))
    i_d = int(getattr(re_obj, "i_discrete", 0))
    v_d = int(getattr(re_obj, "v_discrete", 0))
    I_bin = int(getattr(re_obj, "I_bin", 0))
    V_bin = int(getattr(re_obj, "V_bin", 0))
    i_vac = int(getattr(re_obj, "i_VAC", i_d + P * I_bin))

    def _axis(start, n_bins, n_disc, bins, n_max, disc_start):
        """Reconstruct one axis and its top-bin content fraction."""
        if n_bins <= 0:
            c = np.zeros(n_max)
            c[:n_disc] = yj[disc_start:disc_start + n_disc]
            return c, None
        mom = yj[start:start + P * n_bins]
        if mom.size < P * n_bins:
            return None, None
        mu0 = mom[0::P][:n_bins]
        mu1 = mom[1::P][:n_bins] if P >= 2 else None
        mu2 = mom[2::P][:n_bins] if P >= 3 else None
        c = reconstruct_distribution(sf, mu0, mu1, mu2, bins, n_max)
        c[:n_disc] = yj[disc_start:disc_start + n_disc]
        # Top-bin content fraction from the RAW moments.  mu1 is content
        # (sum n*c_n); with P==1 there is no mu1, so fall back to mu0 weighted
        # by the bin's geometric midpoint.
        if mu1 is not None:
            tot = float(np.sum(np.maximum(mu1, 0.0)))
            frac = float(max(mu1[-1], 0.0)) / tot if tot > 0 else None
        else:
            mid = np.array([(lo + hi - 1) / 2.0 for lo, hi in bins])
            w = mid * np.maximum(mu0, 0.0)
            tot = float(w.sum())
            frac = float(w[-1]) / tot if tot > 0 else None
        return c, frac

    c111, topbin["111"] = _axis(i_d, I_bin, i_d, getattr(re_obj, "bins", []),
                                I, 0)
    c_v, topbin["v"] = _axis(i_vac + v_d, V_bin, v_d,
                             getattr(re_obj, "vac_bins", []), V, i_vac)

    # <100> is emitted per-size in both modes (y_sia100), so its top-bin figure
    # is taken on the SIA bin partition it shares with <111>.
    y100 = res.get("y_sia100")
    if y100 is not None:
        y100 = np.asarray(y100)
        if y100.ndim == 2 and y100.shape[0] >= I and I_bin > 0:
            c100 = np.maximum(y100[:I, -1], 0.0)
            n_ax = np.arange(1, I + 1)
            k100 = n_ax * c100
            tot = float(k100.sum())
            if tot > 0:
                lo, hi = getattr(re_obj, "bins", [(I, I + 1)])[-1]
                topbin["100"] = float(k100[lo - 1:min(hi - 1, I)].sum()) / tot
    return c111, c_v, True, topbin


def observe(res, sim, cfg, d_min_nm):
    """Observation operator + the admissibility block."""
    out = {}
    d = sim.input_data
    I = cfg["I"]
    Om = float(d.derived["Omega"])
    b100 = float(d.energetics.get("b_100", 0.2867)) * 1e-9
    b111 = B111_M
    t = np.asarray(res["t"], float)

    def ser(k):
        if k not in res:
            return np.zeros_like(t)
        v = np.asarray(res[k], float).ravel()
        return v if v.size == t.size else np.zeros_like(t)

    G = float(d.reactions["G"])
    out["dose_reached"] = float(t[-1] * G)
    out["starved"] = bool(out["dose_reached"] < 0.95 * cfg["dose"])
    out["delta_FP"] = float(ser("delta_FP")[-1])
    out["delta_He"] = float(ser("delta_He")[-1])

    N111, N100 = float(ser("N_loops_111")[-1]), float(ser("N_loops_100")[-1])
    n111, n100 = float(ser("mean_n_111")[-1]), float(ser("mean_n_100")[-1])
    out["N_loops_111"], out["N_loops_100"] = N111, N100
    out["mean_n_111"], out["mean_n_100"] = n111, n100
    out["d_111_nm"] = float(2*np.sqrt(max(n111, 0)*Om/(np.pi*b111))*1e9)
    out["d_100_nm"] = float(2*np.sqrt(max(n100, 0)*Om/(np.pi*b100))*1e9)
    out["N_voids"] = float(ser("N_voids")[-1])
    # Cavity SIZE, not just density.  The database reports both, and with
    # swelling withdrawn as an observable (it is N_voids x mean_n_v exactly, so
    # it carries no independent information -- author, 2026-08-02) the void
    # channel would otherwise be density-only.  Spherical cavity:
    #   d = 2 (3 m Omega / 4 pi)^(1/3)          [CLAUDE.md S8 / plan S3.1-2]
    nv = float(ser("mean_n_v")[-1])
    out["mean_n_v"] = nv
    out["d_cavity_nm"] = float(
        2.0 * (3.0 * max(nv, 0.0) * Om / (4.0 * np.pi)) ** (1.0 / 3.0) * 1e9)
    # Still emitted (free, and a useful health check), but NOT screened or
    # calibrated on -- see OBSERVABLES in merge_and_sobol.py.  Note this is the
    # vacancy-INVENTORY swelling of the conservation identity, which includes
    # sub-visible clusters; the database reports S = (pi/6) N_c d_c^3 from the
    # VISIBLE cavity population.  They are different numbers (plan S3.1-3).
    out["S_inventory"] = float(ser("swelling")[-1])

    # ── DOSE LADDER ──────────────────────────────────────────────────────────
    # Accepting a timed-out row AT WHATEVER DOSE IT REACHED would be a mistake
    # of a different order than accepting a truncated tail.  Row 24 gives
    # d_100 = 5.36 nm at 0.1 dpa and 26.505 nm at 1.0 dpa -- a factor of 5 in
    # the primary observable from dose alone, against the ~2 % the whole
    # I=50000 -> 200000 grid refinement moved it.  And because a row is
    # expensive for PHYSICS reasons, `dose_reached` is correlated with theta:
    # it is a confound aligned with the quantity being estimated, not noise.
    #
    # The fix costs nothing, because the trajectory is already here (n_points).
    # Record the observables on a fixed dose ladder as well as at the end, so
    # merge_and_sobol can screen every row at a COMMON dose.  A timed-out row
    # is then not discarded -- it contributes to every checkpoint it actually
    # reached, and only drops out above that, pairwise, per checkpoint.
    #
    # f_100_* is deliberately absent: it needs the per-size reconstruction at
    # that time index, and every f_100 variant is already withheld under
    # bin_moment by BIN_MOMENT_BLOCKED, so there is nothing to screen with it.
    dose_t = t * G
    ladder = {}
    for dck in DOSE_CHECKPOINTS:
        if out["dose_reached"] < dck * (1.0 - 1e-9):
            continue                      # never got here; not an error
        j = int(np.searchsorted(dose_t, dck, side="right")) - 1
        if j < 0:
            continue
        n111c = float(ser("mean_n_111")[j])
        n100c = float(ser("mean_n_100")[j])
        nvc = float(ser("mean_n_v")[j])
        ladder[f"{dck:g}"] = {
            "dose": float(dose_t[j]),
            "N_loops_111": float(ser("N_loops_111")[j]),
            "N_loops_100": float(ser("N_loops_100")[j]),
            "d_111_nm": float(2*np.sqrt(max(n111c, 0)*Om/(np.pi*b111))*1e9),
            "d_100_nm": float(2*np.sqrt(max(n100c, 0)*Om/(np.pi*b100))*1e9),
            "N_voids": float(ser("N_voids")[j]),
            "mean_n_v": nvc,
            "d_cavity_nm": float(
                2.0*(3.0*max(nvc, 0.0)*Om/(4.0*np.pi))**(1.0/3.0)*1e9),
            "delta_FP": float(ser("delta_FP")[j]),
        }
    out["at_dose"] = ladder or None

    # --- admissibility -----------------------------------------------------
    out["occ_111"] = n111 / I
    out["occ_100"] = n100 / I
    d_ceil = float(2*np.sqrt(I*Om/(np.pi*b100))*1e9)
    out["d_ceiling_100_nm"] = d_ceil
    out["d_over_ceiling_100"] = out["d_100_nm"] / d_ceil if d_ceil > 0 else None

    # PER-SIZE EXTRACTION.  In discrete mode `y` is the per-size state and this
    # is a slice.  In bin_moment mode `y` is (N_eq, nt) with N_eq ~ 150 -- the
    # old code sliced y[:I] out of it, got N_eq rows back, and died on the
    # broadcast against arange(1, I+1).  Worse, the vacancy guard below tested
    # `y.shape[0] >= I + V`, which is false for a binned state, so pile_v fell
    # through to None and the gate read None as PASSING (plan S11(h)).
    c111, c_v, binned, topbin = per_size_populations(res, sim, cfg)
    y100 = res.get("y_sia100")
    pile111 = pile100 = None
    f_num = f_cont = f_tem = None
    if c111 is not None:
        n_ax = np.arange(1, I + 1)
        top = int((1 - PILE_TOP_FRAC) * I)
        c111 = np.maximum(c111 - cfg["C_floor"], 0.0)
        k111 = n_ax * c111
        if k111.sum() > 0:
            pile111 = float(k111[top:].sum() / k111.sum())
        if y100 is not None:
            y100 = np.asarray(y100)
            # y_sia100 is emitted PER-SIZE (I, nt) in BOTH modes -- verified
            # 2026-08-03 on a binned run (N_eq=55, y_sia100 shape (200, 20)) --
            # so the <100> path needs no reconstruction.
            if y100.ndim == 2 and y100.shape[0] >= I:
                    c100 = np.maximum(y100[:I, -1] - cfg["C_floor"], 0.0)
                    k100 = n_ax * c100
                    if k100.sum() > 0:
                        pile100 = float(k100[top:].sum() / k100.sum())
                    # loop-fraction conventions (plan S3.1-4): all three
                    s0, s1 = float(c100.sum()), float(c111.sum())
                    f_num = s0 / (s0 + s1) if (s0 + s1) > 0 else None
                    q0, q1 = float(k100.sum()), float(k111.sum())
                    f_cont = q0 / (q0 + q1) if (q0 + q1) > 0 else None
                    # The cutoff is applied to BOTH characters: they have
                    # different b, so the same n maps to a different d.
                    dd100 = 2*np.sqrt(n_ax*Om/(np.pi*b100))*1e9
                    dd111 = 2*np.sqrt(n_ax*Om/(np.pi*b111))*1e9
                    # Emit the whole cutoff curve rather than one sampled value:
                    # d_min is post-processing, so this is free, and the spread
                    # across cutoffs is itself the floor on sigma_model for
                    # f_100 (plan S3.1-4; convention alone moved one measurement
                    # 0.824 -> 0.408).
                    for dm in cfg.get("d_min_sweep", [0.8, 1.0, 1.25, 1.5]):
                        a0 = float(c100[dd100 >= dm].sum())
                        a1 = float(c111[dd111 >= dm].sum())
                        tag = f"{dm:g}".replace(".", "p")
                        out[f"f_100_tem_{tag}"] = ((a0 / (a0 + a1))
                                                   if (a0 + a1) > 0 else None)
                        out[f"N_100_vis_{tag}"] = a0 / Om
                        out[f"N_111_vis_{tag}"] = a1 / Om
                    v0 = float(c100[dd100 >= d_min_nm].sum())
                    v1 = float(c111[dd111 >= d_min_nm].sum())
                    f_tem = v0 / (v0 + v1) if (v0 + v1) > 0 else None
                    out["N_100_visible"] = v0 / Om
                    out["N_111_visible"] = v1 / Om
    out["pile_111"], out["pile_100"] = pile111, pile100
    out["f_100_number"], out["f_100_content"] = f_num, f_cont
    out["f_100_tem"], out["d_min_tem_nm"] = f_tem, d_min_nm

    # --- VACANCY axis -------------------------------------------------------
    # The SIA-only check above is one-sided.  delta_FP is built on the swelling
    # identity S = S_I + Delta J^d, and S is the VACANCY inventory, so a
    # truncated vacancy grid shows up there -- which is how this was found: the
    # check_machine probe (V=120) sits at delta_FP = 4.1e-3 while mean cavity
    # size is ~168 vacancies, i.e. the typical cavity does not fit on its grid.
    #
    # This matters more since swelling was withdrawn and d_cavity_nm added:
    # N_voids / mean_n_v / d_cavity_nm are the observables a truncated V
    # corrupts, and nothing was checking that axis.
    #
    # State layout (discrete): SIA 0..I-1, vacancies I..I+V-1.
    # (bin_moment): i_discrete + P*I_bin + v_discrete + P*V_bin + He extras,
    # handled by per_size_populations() above.
    V = int(cfg["V"])
    pile_v = None
    if c_v is not None:
        m_ax = np.arange(1, V + 1)
        c_v = np.maximum(c_v - cfg["C_floor"], 0.0)
        k_v = m_ax * c_v
        if k_v.sum() > 0:
            pile_v = float(k_v[int((1 - PILE_TOP_FRAC) * V):].sum() / k_v.sum())
    out["pile_v"] = pile_v
    out["occ_v"] = (out["mean_n_v"] / V) if V > 0 else None

    # occ_v is recorded but deliberately NOT a reject criterion.  The occupancy
    # heuristic assumes a converged mean sits far below the ceiling, which holds
    # for loops but not for cavities: the probe converges at mean_n_v = 153 with
    # occ_v = 0.319 (V=480).  Rejecting on occ_v > 0.10 would demand V >~ 1600 to
    # accept an answer that is already converged at 480.
    #
    # occ_111 / occ_100 are now recorded on the SAME footing -- WITHDRAWN as
    # reject criteria on 2026-08-03 (plan S11(c)-2).  Probe: row 24 at
    # occ_100 = 0.171 with pile_100 = 2.3e-10 returns d100 = 5.3573 and
    # N100 = 3.99268e21 at BOTH I=3200 and I=12800.  Occupancy is not a
    # truncation test; pile is.  Under the old rule 71 of 275 v2 rows with no
    # measurable truncation were discarded.
    bad_pile = (pile100 is not None and pile100 > PILE_TOL) or \
               (pile111 is not None and pile111 > PILE_TOL) or \
               (pile_v is not None and pile_v > PILE_TOL)

    # TOP-BIN GATE (bin_moment only).  `pile` is measured on the RECONSTRUCTED
    # per-size distribution, so it inherits the closure's smoothing and can look
    # clean while the top BIN -- the thing that actually holds the overflow --
    # is loaded.  Gate on the raw moments as well.
    bad_topbin = False
    for axis in ("111", "100", "v"):
        val = topbin.get(axis)
        out[f"topbin_{axis}"] = val
        if val is not None and val > TOPBIN_TOL:
            bad_topbin = True

    # FAIL CLOSED.  A missing truncation measure is not a pass.  Before
    # 2026-08-03 an uncomputable pile_v silently satisfied the gate, which is
    # exactly how a binned run would have reported an unguarded vacancy axis as
    # clean.  If we could not measure an axis, we do not certify it.
    unmeasured = [nm for nm, val in (("pile_111", pile111), ("pile_100", pile100),
                                     ("pile_v", pile_v)) if val is None]
    if binned:
        unmeasured += [f"topbin_{ax}" for ax in ("111", "100", "v")
                       if topbin.get(ax) is None]
    out["unmeasured_gates"] = unmeasured or None

    # Still COMPUTED and RECORDED -- it is the covariate the ranking-stability
    # check needs (merge_and_sobol --require-converged) and the only way to say
    # afterwards how deep the truncation ran.  It simply no longer rejects.
    out["grid_limited"] = bool(bad_pile or bad_topbin or unmeasured)
    out["grid_converged"] = not out["grid_limited"]      # explicit, for clarity

    if TRUNCATION_GATES:
        out["admissible"] = bool((not out["starved"]) and (not out["grid_limited"])
                                 and abs(out["delta_FP"]) < DFP_TOL)
    elif STARVED_GATE:
        out["admissible"] = not out["starved"]
    else:
        # A row is admissible if it EXISTS.  The solver returned a trajectory;
        # how far up the dose ladder that trajectory reaches is recorded in
        # `at_dose` and decided per-checkpoint downstream, not here.
        out["admissible"] = True
    return out


_PATTERN_CACHE = {}


def _weighted_pattern(of, weights):
    """Deterministic machine-assignment pattern of length sum(weights).

    Sainte-Lague / highest-averages: repeatedly hand the next slot to the
    machine with the largest w_k / (2 a_k + 1), where a_k is how many slots it
    already holds.  This spreads each machine's slots as evenly as the weights
    allow, rather than giving it one contiguous block.

    Interleaving is the point, not an aesthetic: the design rows are Saltelli
    A / B / AB_j matrices, so a contiguous block belongs disproportionately to
    a few base samples.  If one machine dies, scattered gaps cost pairwise
    deletion a few rows spread over many base samples; a contiguous gap would
    wipe out whole base samples and take their partners with them.

    Pure function of (of, weights) -- every machine computes the same pattern
    from the same CLI arguments, so no coordination is needed.
    """
    key = (of, tuple(weights))
    if key in _PATTERN_CACHE:
        return _PATTERN_CACHE[key]
    # integerise so the pattern is finite and exactly periodic
    scale = 1
    while any(abs(w * scale - round(w * scale)) > 1e-9 for w in weights) and scale < 10 ** 6:
        scale *= 10
    ints = [max(1, int(round(w * scale))) for w in weights]
    g = 0
    for v in ints:
        g = math.gcd(g, v)
    ints = [v // max(g, 1) for v in ints]

    total = sum(ints)
    assigned = [0] * of
    pattern = []
    for _ in range(total):
        k = max(range(of), key=lambda j: (ints[j] / (2 * assigned[j] + 1), -j))
        pattern.append(k)
        assigned[k] += 1
    _PATTERN_CACHE[key] = pattern
    return pattern


MACHINES_FILE = HERE / "machines.json"


def _host_facts() -> dict:
    """Everything the registry is allowed to match on."""
    import subprocess as _sp
    hw_model = ""
    if platform.system() == "Darwin":
        try:
            hw_model = _sp.run(["sysctl", "-n", "hw.model"], capture_output=True,
                               text=True, timeout=5).stdout.strip()
        except Exception:
            pass
    return {"system": platform.system(), "node": platform.node(),
            "hw_model": hw_model, "cpu_count": os.cpu_count() or 0,
            "env": set(os.environ)}


def _matches(rule: dict, f: dict) -> bool:
    if "system" in rule and rule["system"] != f["system"]:
        return False
    if "hw_model_prefix" in rule and not f["hw_model"].startswith(rule["hw_model_prefix"]):
        return False
    if "cpu_count" in rule and int(rule["cpu_count"]) != f["cpu_count"]:
        return False
    if "node_regex" in rule and not re.search(rule["node_regex"], f["node"], re.I):
        return False
    if "env_any" in rule and not any(v in f["env"] for v in rule["env_any"]):
        return False
    # A rule with no clause at all would match everything -- refuse it.
    return bool(rule)


def load_registry(path: Path = None) -> dict:
    p = path or MACHINES_FILE
    if not p.exists():
        raise SystemExit(f"--machine auto needs {p}, which does not exist.")
    return json.loads(p.read_text(encoding="utf-8"))


def detect_machine(reg: dict, facts: dict = None):
    """Which registry entry is this host?  Returns the entry, or raises.

    FAILS LOUDLY ON NO MATCH, and equally loudly on an AMBIGUOUS match.  Both
    are the same underlying hazard: an index chosen wrongly means two machines
    compute the same rows and some other rows are computed by nobody, and that
    surfaces at merge time as "rows MISSING" -- which reads like a machine that
    never reported rather than the misconfiguration it is.  A campaign is
    cheaper to start late than to re-run.
    """
    f = facts or _host_facts()
    hits = [m for m in reg["machines"]
            if any(_matches(rule, f)
                   for rule in [m.get("match", {})] + m.get("match_alternatives", []))]
    if len(hits) == 1:
        return hits[0]
    where = (f"node={f['node']!r} system={f['system']!r} "
             f"hw_model={f['hw_model']!r} cpu_count={f['cpu_count']}")
    if not hits:
        raise SystemExit(
            f"\n*** machine detection FAILED -- this host matches no entry in "
            f"{MACHINES_FILE.name}.\n    {where}\n"
            f"    Known: " + ", ".join(f"{m['index']}={m['name']}"
                                       for m in reg["machines"]) + "\n"
            f"    Fix it in ONE of two ways, never by guessing an index:\n"
            f"      python run_ensemble.py --register <index>   # pin this host\n"
            f"      python run_ensemble.py --machine <index> ... # one-off override\n")
    raise SystemExit(
        f"\n*** machine detection AMBIGUOUS -- this host matches "
        f"{[m['index'] for m in hits]}.\n    {where}\n"
        f"    Tighten the match rules in {MACHINES_FILE.name}; two machines "
        f"sharing an index compute the same rows.\n")


def register_host(index: int, path: Path = None) -> None:
    """Pin THIS host to a registry index, by writing a fingerprint it will match.

    Deliberately narrow: node + system + cpu_count.  A loose rule is how two
    hosts end up sharing an index.
    """
    p = path or MACHINES_FILE
    reg = load_registry(p)
    f = _host_facts()
    ent = next((m for m in reg["machines"] if m["index"] == index), None)
    if ent is None:
        raise SystemExit(f"--register {index}: no such index in {p.name}")
    rule = {"system": f["system"], "node_regex": "^" + re.escape(f["node"]) + "$",
            "cpu_count": f["cpu_count"]}
    if f["hw_model"]:
        rule["hw_model_prefix"] = f["hw_model"]
    ent.setdefault("match_alternatives", []).append(rule)
    p.write_text(json.dumps(reg, indent=2), encoding="utf-8")
    print(f"registered this host to machine {index} ({ent['name']}): {rule}")
    print(f"*** COMMIT AND PUSH {p.name} -- every participant reads it.")


def assign_machine(row_id, of, weights=None):
    """Which machine owns this row_id.

    weights=None reproduces the original even split (row_id % of) EXACTLY, so
    an existing campaign resumed without --weights keeps its assignment and
    every already-computed row stays where it was.
    """
    if not weights:
        return row_id % of
    pattern = _weighted_pattern(of, weights)
    return pattern[row_id % len(pattern)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", type=Path, required=True)
    # DETECTED, not typed.  `auto` reads machines.json, matches this host, and
    # takes index/of/weights/workers from it, so the one value whose mistyping
    # silently corrupts a campaign is no longer typed at all.  An integer still
    # works and still overrides.
    ap.add_argument("--machine", default="auto",
                    help="this machine's index k, or 'auto' (default) to detect "
                         "it from machines.json")
    ap.add_argument("--of", type=int, default=0,
                    help="total machines M; taken from machines.json under --machine auto")
    ap.add_argument("--registry", type=Path, default=MACHINES_FILE,
                    help="the machine registry (default machines.json)")
    ap.add_argument("--register", type=int, default=None, metavar="INDEX",
                    help="pin THIS host to a registry index and exit. Use once "
                         "per new participant, then commit machines.json.")
    # Hoffman2 is ONE machine index (3) but runs as a 16-task array job, so its
    # share is split a second time INSIDE the index.  Each subtask takes every
    # Kth row of machine 3's assignment and writes its own results file, so
    # tasks never collide and a task killed at h_rt costs only its own rows.
    ap.add_argument("--subtask", type=int, default=0,
                    help="this array task's 0-based id within the machine")
    ap.add_argument("--of-subtasks", type=int, default=1,
                    help="number of array tasks sharing this machine index")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--weights", default=None,
                    help="relative capacity of EACH machine, either a comma list "
                         "('6,14,22,22', normally each machine's --workers) or "
                         "'@path' to a campaign_layout.json written by "
                         "campaign_layout.py. Splits the design in proportion to "
                         "capacity instead of evenly, so heterogeneous machines "
                         "finish together. MUST be identical on every machine -- "
                         "prefer '@' form: a 66-entry list retyped on six hosts is "
                         "a silent row collision waiting to happen, whereas the "
                         "file is committed once and its hash is recorded per row. "
                         "Omit for the default even split (row_id %% M).")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--conditions", type=Path, default=HERE / "conditions.json")
    # FROZEN GRID, author 2026-08-06 (plan S12(u)).  The single source of truth
    # is machines.json["grid"]; these defaults mirror it so a bare invocation
    # cannot silently run something else.  Measured at 10-way concurrency
    # against the I=30000/V=5000/b25/i_mobile=50 base (4171 s on design row 21):
    #   I=10000            1838 s (2.3x), screened observables move <=6%
    #   V stays 5000       V=2000 collapses mean_n_v by 70%
    #   dose stays 1.0     halving it saves 17%, not 50%
    # i_mobile 40 and bins 20 are the author's call; i_mobile is PHYSICS, so
    # these rows are NOT comparable with anything run at i_mobile=50.
    ap.add_argument("--I", type=int, default=10000)
    ap.add_argument("--V", type=int, default=5000)
    ap.add_argument("--dose", type=float, default=1.0)
    ap.add_argument("--equations", default="discrete",
                    choices=["discrete", "bin_moment"])
    # --- bin_moment layout (ignored unless --equations bin_moment) -----------
    # These are NOT constructor kwargs; apply_bin_config() writes them to the
    # workbook dict and bin_layout() verifies they landed.  Defaults are the
    # plan S11(f) rev-6 campaign values, NOT the library defaults -- the library
    # default (i_discrete=10, I_bin=6) is the coarse binning plan S10(i) warns
    # about, and silently running it is the failure this whole block prevents.
    ap.add_argument("--i-discrete", type=int, default=40,
                    help="max individually-tracked SIA size; must be >= i_mobile")
    ap.add_argument("--v-discrete", type=int, default=5,
                    help="max individually-tracked vacancy size; must be >= v_mobile")
    ap.add_argument("--i-bin", type=int, default=20,
                    help="target SIA bin count; r=(I/i_discrete)^(1/I_bin) is DERIVED")
    ap.add_argument("--v-bin", type=int, default=20,
                    help="target vacancy bin count; r derived the same way")
    ap.add_argument("--shape-function", default="linear",
                    choices=["constant", "linear", "lognormal"],
                    help="intra-bin closure: P=1/2/3 moments per bin")
    # Used only when the DESIGN does not carry an i_mobile / v_mobile column.
    # Rev 6 withdraws both from theta (plan S11(f)) precisely because a sampled
    # i_mobile would change the bin layout row to row, so under rev 6 these are
    # the campaign values, not fallbacks.
    # 40, frozen 2026-08-06. This is PHYSICS, not a numerical knob: it is the
    # maximum mobile SIA cluster size, and it was a Tier-2 theta with a prior
    # before bin_moment required it fixed. Measured at 30 vs 50 it moved
    # N_100 +109%, N_111 -18%, d_111 -9%, so rows produced here CANNOT be
    # pooled with anything run at i_mobile=50 -- run_cfg_sha enforces that.
    ap.add_argument("--i-mobile-default", type=int, default=40)
    ap.add_argument("--v-mobile-default", type=int, default=5)
    ap.add_argument("--solver-mode", default="active_window",
                    choices=["active_window", "full_system"],
                    help="active_window (default) keeps the ACTIVE system at "
                         "50-200 unknowns regardless of I. Verified 2026-08-02 to "
                         "reproduce full_system to printed precision at I=1600 "
                         "(d100/occ/pile/N100 all identical) at 2.4x lower cost "
                         "and 3 orders better delta_FP. Preconditioner follows "
                         "this choice automatically.")
    ap.add_argument("--rtol", type=float, default=1e-6)
    # A BUDGET, NOT A CLIFF, since 2026-08-05 (plan S12(s)).  On expiry the
    # solver is asked to finalize gracefully and its partial trajectory is KEPT;
    # the row then contributes to every rung of the dose ladder it reached.  So
    # this no longer needs to cover the p99 to avoid losing work -- pick it from
    # the cost distribution so the COMMON RUNG lands high.  3600 s is a smoke-
    # test default and is far below the ~3130 s median / long tail of the rev-6
    # baseline grid; the campaign scripts set it explicitly.
    ap.add_argument("--timeout-s", type=float, default=3600)
    ap.add_argument("--limit", type=int, default=0, help="run only the first n rows (smoke test)")
    ap.add_argument("--stop-after-s", type=float, default=0.0,
                    help="stop SUBMITTING new rows after this many seconds and exit "
                         "cleanly (0 = no limit). For batch schedulers: set it to "
                         "the job's walltime minus one row's p99, so the worker "
                         "parks itself before the scheduler SIGKILLs it mid-row. "
                         "Rows already running are allowed to finish and are "
                         "written; a resubmission resumes by skipping them.")
    ap.add_argument("--allow-mixed", action="store_true",
                    help="resume even though git/solver/workbook/design changed "
                         "since the existing rows (see RESTART SAFETY)")
    a = ap.parse_args(argv)

    # ── who am I? ────────────────────────────────────────────────────────────
    if a.register is not None:
        register_host(a.register, a.registry)
        return 0
    detected = None
    if str(a.machine).lower() == "auto":
        reg = load_registry(a.registry)
        detected = detect_machine(reg)
        a.machine = int(detected["index"])
        a.of = a.of or int(reg["of"])
        if a.weights is None:
            a.weights = "@" + str(a.registry)
        if a.of_subtasks == 1 and detected.get("subtasks"):
            a.of_subtasks = int(detected["subtasks"])
        # Slots come from the registry too: a machine whose --workers does not
        # match the slot count its WEIGHT was computed from will finish out of
        # step with everyone else, which is the thing the weights exist to stop.
        a.workers = int(detected.get("slots", a.workers))
        if a.of_subtasks > 1:
            a.workers = max(1, a.workers // a.of_subtasks)
        print(f"  detected machine {a.machine} = {detected['name']} "
              f"({detected.get('speed_source','?')} speed "
              f"{detected.get('speed')}), {a.workers} worker(s)"
              + (f", subtask {a.subtask}/{a.of_subtasks}" if a.of_subtasks > 1 else ""))
    else:
        a.machine = int(a.machine)
        if not a.of:
            raise SystemExit("--of is required when --machine is given explicitly")
    if not 0 <= a.subtask < a.of_subtasks:
        raise SystemExit(f"--subtask {a.subtask} outside 0..{a.of_subtasks-1}")

    spec = json.loads((HERE / "parameters.json").read_text(encoding="utf-8"))
    rows, meta = read_design(a.design)
    conds = (json.loads(a.conditions.read_text(encoding="utf-8"))
             if a.conditions.exists() else
             {"N2": {"T_K": 573.15}, "N5": {"T_K": 623.15},
              "I1": {"T_K": 623.15, "cascade": "fission"}})

    weights = None
    weights_sha = ""
    if a.weights:
        if a.weights.startswith("@"):
            lay_p = Path(a.weights[1:])
            if not lay_p.exists():
                raise SystemExit(f"--weights {a.weights}: {lay_p} does not exist. "
                                 f"Generate it with campaign_layout.py and commit "
                                 f"it, so every machine reads the same file.")
            lay = json.loads(lay_p.read_text(encoding="utf-8"))
            a.weights = lay["weights"]
            weights_sha = sha256_file(lay_p)[:16]
            if int(lay.get("of", 0)) != a.of:
                raise SystemExit(
                    f"--of {a.of} but {lay_p.name} was built for "
                    f"of={lay['of']}. They must agree, or the machines disagree "
                    f"about how many participants exist and rows are both "
                    f"duplicated and dropped. Re-run campaign_layout.py.")
            n_part = len(lay.get("participants") or lay.get("machines") or [])
            print(f"  weights from {lay_p.name} (sha {weights_sha}): "
                  f"{n_part} participants")
        weights = [float(w) for w in a.weights.split(",")]
        if len(weights) != a.of:
            raise SystemExit(f"--weights has {len(weights)} entries but --of is "
                             f"{a.of}; they must match, and the SAME list must be "
                             f"passed on every machine.")
        if min(weights) <= 0:
            raise SystemExit("--weights entries must all be > 0.")
    mine = [r for r in rows if assign_machine(r["row_id"], a.of, weights) == a.machine]
    # SECOND-LEVEL SPLIT, for a machine that runs as a scheduler array job.
    # Deterministic and stateless: task t takes every of_subtasks-th row of this
    # machine's assignment, in design order.  Nothing is shared between tasks,
    # so they cannot collide, and a task killed at the walltime limit costs only
    # its own in-flight rows.  of_subtasks=1 is the identity.
    if a.of_subtasks > 1:
        n_before = len(mine)
        mine = mine[a.subtask::a.of_subtasks]
        print(f"  subtask {a.subtask}/{a.of_subtasks}: {len(mine)} of this "
              f"machine's {n_before} rows")
    if a.limit:
        mine = mine[:a.limit]
    _sfx = f"_machine{a.machine}" + (f"_t{a.subtask}" if a.of_subtasks > 1 else "")
    out = a.out or (HERE / "results" / f"{a.design.stem}{_sfx}.jsonl")
    out.parent.mkdir(parents=True, exist_ok=True)

    # resume: skip row_ids already present
    done = set()
    prior = []
    if out.exists():
        for ln in out.read_text(encoding="utf-8").splitlines():
            try:
                r = json.loads(ln)
            except Exception:
                continue          # truncated last line after a hard kill
            done.add(r["row_id"])
            prior.append(r)
    todo = [r for r in mine if r["row_id"] not in done]

    # Provenance only -- the run itself resolves the binary via py_utils.cpp_bridge.
    # Must match cpp_bridge's search order, otherwise solver_sha256 records
    # 'unavailable' on every non-Windows machine and the whole point of hashing
    # the binary (proving four machines ran the SAME solver) is lost.
    _build = REPO / "RadCluster_2_1" / "build"
    for _cand in (_build / "Release" / "solver.exe", _build / "solver.exe",
                  _build / "solver"):
        if _cand.exists():
            solver = _cand
            break
    else:
        solver = _build / "solver"
    # RUN CONFIGURATION IS PROVENANCE.  I/V/dose/equations/rtol are NOT part of
    # theta, so without recording them here nothing in the result row would show
    # that the numerics moved between two batches — and a grid change silently
    # invalidates comparability far more than most theta changes do.  Hashed
    # into one short field so it can be compared cheaply on resume, and kept in
    # full in the manifest.
    # Preconditioner follows the solver mode and is NOT independently selectable.
    # CLAUDE.md S9: Woodbury is a full_system optimisation -- its 58-RHS setup
    # cost is counterproductive against active_window's 50-200 unknown active
    # system, where Jacobi+GMRES converges fine.  Pairing them the other way
    # round is always a mistake, so the code does not offer the choice.
    precond = "Woodbury" if a.solver_mode == "full_system" else "Jacobi"
    run_cfg = {"I": a.I, "V": a.V, "dose": a.dose, "equations": a.equations,
               "rtol": a.rtol, "atol": 1e-20, "C_floor": 1e-25,
               "solver_mode": a.solver_mode, "linsol": "gmres",
               "preconditioner": precond, "loop_conversion": 1,
               "i_mobile_default": a.i_mobile_default,
               "v_mobile_default": a.v_mobile_default,
               "n_points": 40, "timeout_s": a.timeout_s}
    # run_cfg is BOTH the runtime config handed to evaluate() AND the source of
    # run_cfg_sha.  timeout_s must stay IN it -- cpp_bridge needs it -- but must
    # be left OUT of the hash: it is a wall-clock budget, not a model setting,
    # and the budgets have to differ between a Mac and a 4.67x-slower cluster
    # node whose h_rt caps a task at 24 h.  Hashing it made identical physics
    # produce different hashes, which is the false alarm that teaches people to
    # ignore the check.  Excluded at the hash, not removed from the dict:
    # deleting the key made every row die with KeyError: 'timeout_s' (2026-08-05).
    SHA_EXCLUDE = {"timeout_s"}
    # The bin layout changes the numerics as much as I/V do, so it belongs in
    # run_cfg_sha.  Added only in bin_moment mode so an existing discrete
    # campaign's hash -- and therefore its resume check -- is unchanged.
    if a.equations == "bin_moment":
        run_cfg.update({"i_discrete": a.i_discrete, "v_discrete": a.v_discrete,
                        "I_bin": a.i_bin, "V_bin": a.v_bin,
                        "shape_function": a.shape_function})
        if a.i_discrete >= a.I or a.v_discrete >= a.V:
            raise SystemExit(
                f"bin_moment needs i_discrete < I and v_discrete < V; got "
                f"i_discrete={a.i_discrete} I={a.I}, v_discrete={a.v_discrete} "
                f"V={a.V}. Otherwise I_bin collapses to 0 and the run is "
                f"silently fully discrete.")
        # Mobile sizes must be individually resolved -- a mobile cluster living
        # inside a bin has no per-size concentration for the coalescence and
        # annihilation sums to use.
        if a.i_discrete < a.i_mobile_default or a.v_discrete < a.v_mobile_default:
            raise SystemExit(
                f"bin_moment needs i_discrete >= i_mobile and v_discrete >= "
                f"v_mobile; got i_discrete={a.i_discrete} vs i_mobile="
                f"{a.i_mobile_default}, v_discrete={a.v_discrete} vs v_mobile="
                f"{a.v_mobile_default}.")
        # A SAMPLED i_mobile/v_mobile changes i_discrete, hence the bin layout
        # and N_eq, from row to row.  A variance decomposition across rows that
        # do not share a discretisation is not interpretable, which is why plan
        # S11(f) withdraws both from theta for the bin_moment campaign.  Catch
        # the stale-design case here rather than after 1000 rows.
        sampled = sorted({k for k in ("i_mobile", "v_mobile")
                          if any(k in r for r in rows)})
        if sampled and not a.allow_mixed:
            raise SystemExit(
                f"design carries sampled {'/'.join(sampled)} but --equations "
                f"bin_moment fixes the bin layout from --i-discrete/--v-discrete. "
                f"Rows would not share a discretisation and their indices would "
                f"not be comparable (plan S11(f)).\n"
                f"  Use a design with {'/'.join(sampled)} withdrawn from theta, "
                f"or pass --allow-mixed if you have a specific reason.")
    run_cfg_sha = hashlib.sha256(json.dumps(
        {k: v for k, v in run_cfg.items() if k not in SHA_EXCLUDE},
        sort_keys=True).encode()).hexdigest()[:16]
    prov = {"git_sha": git_sha(), "machine_id": platform.node(),
            "solver_sha256": sha256_file(solver)[:16],
            "workbook_sha256": sha256_file(
                REPO / "RadCluster_2_1" / "input" / "input_parameters.xlsx")[:16],
            "design_sha256": meta.get("design_sha256", "")[:16],
            "run_cfg_sha": run_cfg_sha,
            # Recorded, not hashed -- resource policy, not physics.
            "timeout_s": a.timeout_s, "stop_after_s": a.stop_after_s,
            "workers": a.workers,
            # The row->machine map is a function of (of, weights).  If two
            # machines disagree about it they silently compute overlapping rows
            # AND leave a hole, which looks like "some rows missing" at merge
            # time rather than like the configuration error it is.  Carried per
            # row so merge_and_sobol can say which it was.
            "weights_sha": weights_sha,
            # weights_sha above is sha256 of the machines.json FILE, so it moves
            # when anything in that file moves -- timeout_s, a match rule, a
            # comment -- none of which change which machine owns which row.
            # Setting timeout_s 12000 -> 3600 duly reported a PROVENANCE SPLIT
            # between two machines whose partition was verified identical over
            # all 1008 rows.  A checker that cries wolf is worse than no checker:
            # the real split it exists to catch arrives in the same list as the
            # noise.  So hash the MAP -- (of, weights) and nothing else -- and
            # leave the file hash beside it as a record of which file was read.
            "weights_map_sha": hashlib.sha256(json.dumps(
                {"of": a.of, "weights": weights}, sort_keys=True
            ).encode()).hexdigest()[:16],
            "of": a.of,
            # carried per row (not only in run_cfg_sha) so merge_and_sobol can
            # apply the per-mode observable-fidelity restriction without having
            # to resolve the hash back to a configuration
            "equations": a.equations,
            "python": platform.python_version()}
    print(f"machine {a.machine}/{a.of}  rows {len(mine)} "
          f"({len(done)} done, {len(todo)} to run)  workers {a.workers}")
    print(f"  provenance {json.dumps(prov)}")

    # RESTART SAFETY.  Resuming is only a benefit if the new rows are
    # comparable to the old ones.  If the code, solver or workbook moved since
    # the existing rows were written, appending to the same file silently mixes
    # two populations into one Sobol estimate.  Refuse unless told otherwise.
    if prior:
        drift = {}
        for f in ("git_sha", "solver_sha256", "workbook_sha256", "design_sha256",
                  "run_cfg_sha"):
            old = {r.get(f) for r in prior if r.get(f)}
            if old and prov[f] not in old:
                drift[f] = (sorted(old), prov[f])
        if drift and not a.allow_mixed:
            print("\n  *** REFUSING TO RESUME - the environment changed since the "
                  "existing rows were written:")
            for f, (old, new) in drift.items():
                print(f"        {f}: was {old} -> now {new}")
            print("\n  Those rows are not comparable to the ones this process "
                  "would produce.")
            print("  Choose one:")
            print("    - restore the previous state (git checkout / rebuild), or")
            print(f"    - archive the old results and start clean:")
            print(f"        mv {out.name} {out.stem}_pre-change.jsonl")
            print("    - or, only if you are certain the change cannot affect "
                  "results\n      (e.g. a comment or README edit), re-run with "
                  "--allow-mixed.")
            return 3
        if drift:
            print(f"  WARNING: resuming across an environment change "
                  f"({', '.join(drift)}) because --allow-mixed was given.")
    if meta.get("revision_pending"):
        print(f"  *** design carries REVISION_PENDING parameters: "
              f"{meta['revision_pending']}")

    # DERIVED from run_cfg, never rebuilt by hand.  These were two separate
    # literals -- run_cfg for the provenance hash, cfg for the workers -- and
    # they had to be edited in lockstep.  Adding solver_mode/preconditioner to
    # run_cfg alone made every row die with KeyError: 'solver_mode'.  Deriving
    # one from the other makes that class of drift impossible, and also
    # guarantees that what the workers actually ran is exactly what run_cfg_sha
    # hashed -- which is the whole point of the restart-safety check.
    cfg = dict(run_cfg)

    # GRACEFUL STOP.  submit/as_completed rather than ex.map: map() has no way
    # to stop feeding work, so a stop request could only be honoured by killing
    # the pool, which loses every in-flight row.  Here we simply stop SUBMITTING
    # once the sentinel appears, let the running rows finish and be written, and
    # exit 0.  Because rows are appended as they complete and the resume filter
    # above skips row_ids already present, a restart picks up exactly where this
    # left off -- no row lost, none recomputed.
    from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
    t0 = time.time()

    # LIVE START MARKER.  The manifest is written at the END of a run, so while
    # a campaign is in flight there is nothing on disk that says when it began
    # -- and without a true elapsed, the only throughput a status tool can form
    # is sum(row_wall)/workers over the rows that HAVE landed.  That is biased
    # low exactly when it matters: early, with the slow rows still running and
    # contributing nothing to the numerator OR the denominator.  Every wrong
    # ETA this campaign produced came from that estimator.  Six lines here make
    # completions-per-hour a division by real wall instead.
    started_p = out.with_suffix(".started.json")
    started_p.write_text(json.dumps({
        "unix": t0, "at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "machine": a.machine, "workers": a.workers,
        "rows_assigned": len(mine), "rows_already_done": len(mine) - len(todo),
        "rows_to_run": len(todo), "machine_id": platform.node(),
    }, indent=2), encoding="utf-8")

    n_ok = n_bad = n_inadm = 0
    stopped = None
    pending = list(todo)
    # WRITE BY PATH, NOT BY HELD HANDLE.  This used to be
    #     with ProcessPoolExecutor(...) as ex, out.open("a") as fh:
    # holding one file object open for the whole multi-day run.  A file handle
    # names an INODE, not a path, so anything that replaces the file underneath
    # -- and `git pull --rebase --autostash` does exactly that when the results
    # file is dirty, which it always is while a worker is running -- leaves the
    # worker appending to an unlinked inode that nothing can ever read.  It is
    # silent: the log keeps printing completed rows, the .jsonl stops growing,
    # and the rows are unrecoverable once the process exits (macOS volfs cannot
    # open an unlinked inode).  Measured 2026-08-06: 5 rows, 6.2 core-hours,
    # written into a dead inode over three hours before anyone noticed.
    #
    # Re-opening per row costs one open()/close() against ~3000 s of solve, and
    # makes a swapped file cost nothing: the next row simply appends to whatever
    # now lives at that path.
    def append_row(rec):
        with out.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            os.fsync(fh.fileno())

    with ProcessPoolExecutor(max_workers=a.workers) as ex:

        def submit(r):
            return ex.submit(evaluate,
                             (r, spec, conds.get(r["condition"], {"T_K": 623.15}),
                              cfg))

        inflight = set()
        while pending and len(inflight) < a.workers:
            inflight.add(submit(pending.pop(0)))

        while inflight:
            finished, inflight = wait(inflight, return_when=FIRST_COMPLETED)
            for fut in finished:
                rec = fut.result()
                rec.update(prov)
                append_row(rec)
                if rec.get("solver_rc"):
                    n_bad += 1
                elif not rec.get("admissible"):
                    n_inadm += 1
                else:
                    n_ok += 1
                print(f"  row {rec['row_id']:6d} {rec['condition']:>4s} "
                      f"{rec['matrix']:>2s} "
                      f"{'FAIL' if rec.get('solver_rc') else ('INADM' if not rec.get('admissible') else ' ok  ')} "
                      f"d100={rec.get('d_100_nm', float('nan')):6.2f} "
                      f"N100={rec.get('N_loops_100', float('nan')):9.3e} "
                      f"pile={rec.get('pile_100')} "
                      f"dFP={rec.get('delta_FP', float('nan')):8.1e} "
                      f"{rec['wall_s']:.0f}s", flush=True)
            if stopped is None and STOP_FILE.exists():
                try:
                    stopped = json.loads(STOP_FILE.read_text(encoding="utf-8"))
                except Exception:
                    stopped = {"reason": "unreadable STOP file"}
                print(f"\n  STOP requested ({stopped.get('reason')}) -- not "
                      f"submitting further rows; {len(inflight)} in flight will "
                      f"finish and be written.\n", flush=True)
            # SELF-IMPOSED DEADLINE, for batch schedulers.  A job killed at its
            # walltime loses every in-flight row -- with 4 workers that is 4 rows
            # of up to an hour each, silently, on every task of an array job.
            # Parking one row's-worth of time early costs a fraction of that and
            # makes the shortfall visible in the manifest instead.  Same path as
            # the STOP sentinel: stop SUBMITTING, let the running rows land.
            if stopped is None and a.stop_after_s > 0:
                el = time.time() - t0
                if el >= a.stop_after_s:
                    stopped = {"reason": f"--stop-after-s {a.stop_after_s:.0f} "
                                         f"reached at {el:.0f}s (scheduler "
                                         f"walltime guard)"}
                    print(f"\n  DEADLINE reached ({el:.0f}s of "
                          f"{a.stop_after_s:.0f}s) -- not submitting further "
                          f"rows; {len(inflight)} in flight will finish and be "
                          f"written. Resubmit to resume.\n", flush=True)
            if stopped is None:
                while pending and len(inflight) < a.workers:
                    inflight.add(submit(pending.pop(0)))

    # ── run-level manifest ───────────────────────────────────────────────────
    # The per-row block carries hashes so a mixed file is detectable; this
    # carries the FULL configuration behind those hashes, so a result set is
    # interpretable years later without re-deriving anything.  Mirrors the
    # repo's output/<stamp>/provenance.md convention (CLAUDE.md "Output
    # format"), in JSON because it is also read back by campaign_ops.
    all_recs = []
    if out.exists():
        for ln in out.read_text(encoding="utf-8").splitlines():
            try:
                all_recs.append(json.loads(ln))
            except Exception:
                pass
    ws = [r["wall_s"] for r in all_recs if r.get("wall_s")]
    man = {
        "schema_version": "1",
        "written_at_unix": time.time(),
        "written_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "campaign": {"design_file": str(a.design), "design_sha256":
                     meta.get("design_sha256"), "tier": meta.get("tier"),
                     "N": meta.get("N"), "p": meta.get("p"),
                     "conditions": meta.get("conditions"),
                     "parameters": meta.get("parameters"),
                     "parameters_version": meta.get("parameters_version"),
                     "revision_pending": meta.get("revision_pending"),
                     "rows_total": meta.get("rows_total")},
        "machine": {"index": a.machine, "of": a.of,
                    "weights": weights, "rows_assigned": len(mine),
                    "machine_id": platform.node(), "platform": platform.platform(),
                    "python": platform.python_version(),
                    "cpu_count": os.cpu_count(), "workers": a.workers,
                    "omp_num_threads": os.environ.get("OMP_NUM_THREADS")},
        "code": {"git_sha": prov["git_sha"], "solver_sha256": prov["solver_sha256"],
                 "solver_path": str(solver),
                 "workbook_sha256": prov["workbook_sha256"]},
        "run_config": run_cfg,
        "run_cfg_sha": run_cfg_sha,
        "conditions_file": str(a.conditions),
        "conditions": conds,
        "rows": {"assigned": len(mine), "completed": len(all_recs),
                 "admissible": sum(1 for r in all_recs
                                   if not r.get("solver_rc") and r.get("admissible")),
                 "inadmissible": sum(1 for r in all_recs
                                     if not r.get("solver_rc") and not r.get("admissible")),
                 "failed": sum(1 for r in all_recs if r.get("solver_rc")),
                 "not_started": len(pending)},
        "timing": {"session_wall_s": round(time.time() - t0, 1),
                   "row_wall_mean_s": (float(np.mean(ws)) if ws else None),
                   "row_wall_median_s": (float(np.median(ws)) if ws else None),
                   "row_wall_p90_s": (float(np.percentile(ws, 90)) if ws else None),
                   "core_hours": (sum(ws) / 3600.0 if ws else 0.0)},
        "stopped": stopped,
        "admissibility_rule": {
            "pile_tol": PILE_TOL, "occ_tol": OCC_TOL,
            "delta_FP_tol": DFP_TOL,
            "pile_top_frac": PILE_TOP_FRAC,
            "axes": ["sia_111", "sia_100", "vacancy"],
            "note": "delta_FP is blind to SIA-axis truncation (the <100> pile-up "
                    "held delta_FP at 1e-8 while 96% of content sat at n=I), so "
                    "grid adequacy is judged by pile/occupancy. It is NOT blind "
                    "to the VACANCY axis -- delta_FP rests on the swelling "
                    "identity S = S_I + Delta J^d and S is the vacancy "
                    "inventory. The vacancy axis was unchecked until 2026-08-02; "
                    "at the check_machine probe (V=120) mean_n_v moved 18.7->136 "
                    "and swelling 58x on going to V=240."},
    }
    man_p = out.with_suffix(".manifest.json")
    man_p.write_text(json.dumps(man, indent=2), encoding="utf-8")

    tag = "STOPPED" if stopped else "done"
    print(f"\n{tag} in {time.time()-t0:.0f}s: {n_ok} admissible, "
          f"{n_inadm} inadmissible, {n_bad} failed -> {out}")
    print(f"  manifest -> {man_p.name}")
    if stopped:
        print(f"  {len(pending)} row(s) were never started and remain assigned "
              f"to this machine.")
        if "walltime guard" in str(stopped.get("reason", "")):
            # No sentinel was written, so telling the operator to clear one sends
            # them looking for a file that does not exist -- and on a cluster
            # this message is read by whoever is debugging a resubmission loop.
            print("  No STOP flag was set: this was the --stop-after-s deadline. "
                  "Resubmit the SAME command to resume; completed rows are "
                  "skipped automatically.")
        else:
            print("  Clear the flag (campaign_ops.clear_stop()) and re-run the same "
                  "command to resume; completed rows are skipped automatically.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
