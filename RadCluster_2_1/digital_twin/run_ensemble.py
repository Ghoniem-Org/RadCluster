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
OCC_TOL = 0.10           # project occupancy rule (SIA axes only -- see observe())

# CLAUDE.md S8's own standard: delta_FP below ~1e-6 is good, above 1e-3 signals
# a coding error.  Admissibility used 1e-3, i.e. it only rejected the "coding
# error" band.  That is too loose for the VACANCY axis, where delta_FP is the
# primary guard: the 2026-08-02 vacancy study passed V=240 at delta_FP=8.3e-4
# (under the old bar) while mean_n_v was still 11 % from converged, and reached
# 1.4e-7 at V=480 where it had converged.  1e-6 separates those two cleanly and
# is comfortably met by good configs (production rows run 1.5e-8 to 2.9e-12).
DFP_TOL = 1e-6


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
            sim.input_data._calculate_derived()
            sim.rebuild_rates()
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

    # --- admissibility -----------------------------------------------------
    out["occ_111"] = n111 / I
    out["occ_100"] = n100 / I
    d_ceil = float(2*np.sqrt(I*Om/(np.pi*b100))*1e9)
    out["d_ceiling_100_nm"] = d_ceil
    out["d_over_ceiling_100"] = out["d_100_nm"] / d_ceil if d_ceil > 0 else None

    y, y100 = res.get("y"), res.get("y_sia100")
    pile111 = pile100 = None
    f_num = f_cont = f_tem = None
    if y is not None:
        y = np.asarray(y)
        if y.ndim == 2:
            n_ax = np.arange(1, I + 1)
            c111 = np.maximum(y[:I, -1] - cfg["C_floor"], 0.0)
            top = int((1 - PILE_TOP_FRAC) * I)
            k111 = n_ax * c111
            if k111.sum() > 0:
                pile111 = float(k111[top:].sum() / k111.sum())
            if y100 is not None:
                y100 = np.asarray(y100)
                if y100.ndim == 2:
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
    V = int(cfg["V"])
    pile_v = None
    if y is not None and getattr(y, "ndim", 0) == 2 and y.shape[0] >= I + V:
        m_ax = np.arange(1, V + 1)
        c_v = np.maximum(y[I:I + V, -1] - cfg["C_floor"], 0.0)
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
    # The vacancy axis is guarded instead by pile_v and by delta_FP, which --
    # unlike on the SIA axis -- is genuinely sensitive here.
    bad_pile = (pile100 is not None and pile100 > PILE_TOL) or \
               (pile111 is not None and pile111 > PILE_TOL) or \
               (pile_v is not None and pile_v > PILE_TOL)
    out["grid_limited"] = bool(bad_pile or out["occ_111"] > OCC_TOL
                               or out["occ_100"] > OCC_TOL)
    out["admissible"] = bool((not out["starved"]) and (not out["grid_limited"])
                             and abs(out["delta_FP"]) < DFP_TOL)
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
    ap.add_argument("--machine", type=int, required=True, help="this machine's index k")
    ap.add_argument("--of", type=int, required=True, help="total machines M")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--weights", default=None,
                    help="comma-separated relative capacity of EACH machine, "
                         "e.g. '6,14,22,22' (normally each machine's --workers). "
                         "Splits the design in proportion to capacity instead of "
                         "evenly, so heterogeneous machines finish together. "
                         "MUST be given identically on every machine. "
                         "Omit for the default even split (row_id %% M).")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--conditions", type=Path, default=HERE / "conditions.json")
    ap.add_argument("--I", type=int, default=800)
    ap.add_argument("--V", type=int, default=600)
    ap.add_argument("--dose", type=float, default=0.1)
    ap.add_argument("--equations", default="discrete")
    ap.add_argument("--solver-mode", default="active_window",
                    choices=["active_window", "full_system"],
                    help="active_window (default) keeps the ACTIVE system at "
                         "50-200 unknowns regardless of I. Verified 2026-08-02 to "
                         "reproduce full_system to printed precision at I=1600 "
                         "(d100/occ/pile/N100 all identical) at 2.4x lower cost "
                         "and 3 orders better delta_FP. Preconditioner follows "
                         "this choice automatically.")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--timeout-s", type=float, default=3600)
    ap.add_argument("--limit", type=int, default=0, help="run only the first n rows (smoke test)")
    ap.add_argument("--allow-mixed", action="store_true",
                    help="resume even though git/solver/workbook/design changed "
                         "since the existing rows (see RESTART SAFETY)")
    a = ap.parse_args(argv)

    spec = json.loads((HERE / "parameters.json").read_text(encoding="utf-8"))
    rows, meta = read_design(a.design)
    conds = (json.loads(a.conditions.read_text(encoding="utf-8"))
             if a.conditions.exists() else
             {"N2": {"T_K": 573.15}, "N5": {"T_K": 623.15},
              "I1": {"T_K": 623.15, "cascade": "fission"}})

    weights = None
    if a.weights:
        weights = [float(w) for w in a.weights.split(",")]
        if len(weights) != a.of:
            raise SystemExit(f"--weights has {len(weights)} entries but --of is "
                             f"{a.of}; they must match, and the SAME list must be "
                             f"passed on every machine.")
        if min(weights) <= 0:
            raise SystemExit("--weights entries must all be > 0.")
    mine = [r for r in rows if assign_machine(r["row_id"], a.of, weights) == a.machine]
    if a.limit:
        mine = mine[:a.limit]
    out = a.out or (HERE / "results" / f"{a.design.stem}_machine{a.machine}.jsonl")
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

    solver = REPO / "RadCluster_2_1" / "build" / "Release" / "solver.exe"
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
               "i_mobile_default": 10, "v_mobile_default": 5,
               "n_points": 40, "timeout_s": a.timeout_s}
    run_cfg_sha = hashlib.sha256(
        json.dumps(run_cfg, sort_keys=True).encode()).hexdigest()[:16]
    prov = {"git_sha": git_sha(), "machine_id": platform.node(),
            "solver_sha256": sha256_file(solver)[:16],
            "workbook_sha256": sha256_file(
                REPO / "RadCluster_2_1" / "input" / "input_parameters.xlsx")[:16],
            "design_sha256": meta.get("design_sha256", "")[:16],
            "run_cfg_sha": run_cfg_sha,
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
    n_ok = n_bad = n_inadm = 0
    stopped = None
    pending = list(todo)
    with ProcessPoolExecutor(max_workers=a.workers) as ex, \
            out.open("a", encoding="utf-8") as fh:

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
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
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
        print("  Clear the flag (campaign_ops.clear_stop()) and re-run the same "
              "command to resume; completed rows are skipped automatically.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
