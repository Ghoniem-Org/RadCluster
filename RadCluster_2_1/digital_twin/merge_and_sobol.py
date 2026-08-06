#!/usr/bin/env python
"""merge_and_sobol - concatenate the machines' results and estimate S_i, S_i^T.

    python merge_and_sobol.py --design design/T2_design_v1.csv \
                              --results results/ --out report/

Merging is plain concatenation keyed on row_id, so it is order-independent and
idempotent: re-running after a machine catches up simply adds the new rows.

THE SUBTLETY THAT MAKES THIS FILE WORTH READING
-----------------------------------------------
Saltelli estimators are PAIRWISE.  For parameter j and base index i they need
the triple ( f(A_i), f(B_i), f(AB_i^(j)) ).  The plan says "discard and flag"
failed runs, which is right for the physics but incomplete here:

  * dropping a row globally biases EVERY index, because the A_i and B_i rows are
    shared across all p parameters;
  * dropping nothing and imputing is worse.

The correct handling is PAIRWISE DELETION: if AB_i^(j) is unusable, base index i
is excluded from parameter j's estimator ONLY.  If A_i or B_i is unusable, base
index i is excluded from all p estimators (they genuinely share those rows).
Every index therefore has its own effective sample size, reported as n_eff --
an index computed from 3 of 16 base points is not comparable to one computed
from 16 and must not be read as a screening result.

"Unusable" means failed OR inadmissible.  An inadmissible row (grid-limited or
dose-starved) ran to completion and conserves, but its observables are an
artefact of the numerics; feeding it to the estimator measures the grid, not
the physics.  See run_ensemble.py for why delta_FP cannot detect this.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# f_100_tem is emitted at every cutoff in the parameters.json d_min sweep, so the
# convention spread shows up as the spread ACROSS these columns rather than being
# buried in a variance decomposition (d_min is post-processing, not a design axis).
# SWELLING IS WITHDRAWN AS AN OBSERVABLE (author, 2026-08-02).  Two reasons,
# the second the stronger:
#   (a) in ferritic/martensitic steels swelling is small and not a quantity the
#       campaign needs to reproduce accurately; and
#   (b) it carries NO INDEPENDENT INFORMATION.  S_inventory is the vacancy
#       content, i.e. N_voids x mean_n_v exactly -- verified numerically:
#       (1-0.115)(1-0.269) = 0.647, reproducing the measured 35.3 % bin_moment
#       error to three digits.  Feeding a deterministic function of two other
#       observables into a variance decomposition double-counts the void
#       channel and inflates the apparent number of constraints on theta.
# It is still EMITTED by run_ensemble (free, and a useful health check); it is
# simply not screened or calibrated on.  d_cavity_nm/mean_n_v replace it so the
# void channel keeps its size information rather than becoming density-only.
OBSERVABLES = ["N_loops_100", "d_100_nm", "f_100_number", "f_100_content",
               "f_100_tem_0p8", "f_100_tem_1", "f_100_tem_1p25", "f_100_tem_1p5",
               "N_loops_111", "d_111_nm", "N_voids", "d_cavity_nm"]

# ── OBSERVABLE FIDELITY BY EQUATIONS MODE ────────────────────────────────────
# An observable is only screenable in a mode where it has been VALIDATED against
# a converged discrete reference.  This is not bookkeeping: an index computed
# from an observable the closure gets wrong is a confident number about the
# numerics, not the physics — the same failure mode as the grid saturation that
# delta_FP could not see.
#
# bin_moment: the <100> block is a SESSILE, ONE-WAY population — born small by
# conversion, then advected upward in size space with no back-diffusion.  Moment
# closures suffer numerical diffusion on sharply-peaked advecting distributions,
# which broadens the reconstruction and biases the mean UP.  Measured 2026-08-02
# against a converged discrete reference (pile = 0.0000, delta_FP = 3.5e-07):
# mean_n_100 648.0 vs 265.2, +144 %, and the error GREW with refinement of the
# comparison (+58 % at the coarser setting).  The 1/2<111> block is unaffected
# (mean_n_111 within 6-7 % in two independent comparisons) because it is
# populated across its whole range rather than peaked and advecting.
#
# Author's decision 2026-08-02: if the <100> closure problem persists, screen
# with bin_moment on the 1/2<111>/void observables only.  Encoded here so the
# restriction is enforced by the tooling rather than remembered.
BIN_MOMENT_BLOCKED = {"N_loops_100", "d_100_nm", "f_100_number", "f_100_content",
                      "f_100_tem_0p8", "f_100_tem_1", "f_100_tem_1p25",
                      "f_100_tem_1p5"}


def observables_for_mode(equations: str, blocked_override=None) -> tuple[list, set]:
    """Observables screenable in this equations mode, and those withheld."""
    if equations != "bin_moment":
        return list(OBSERVABLES), set()
    blocked = set(blocked_override if blocked_override is not None
                  else BIN_MOMENT_BLOCKED)
    return [o for o in OBSERVABLES if o not in blocked], blocked


def equations_mode_of(recs: dict) -> str | None:
    """What mode produced these rows?  Read from the manifest-backed rows."""
    modes = {r.get("equations") for r in recs.values() if r.get("equations")}
    if len(modes) == 1:
        return modes.pop()
    return None if not modes else "MIXED"


def load_results(res_dir: Path) -> dict[int, dict]:
    recs: dict[int, dict] = {}
    files = sorted(res_dir.glob("*.jsonl"))
    if not files:
        raise SystemExit(f"no .jsonl files in {res_dir}")
    dup = 0
    for f in files:
        for ln in f.read_text(encoding="utf-8").splitlines():
            if not ln.strip():
                continue
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue                       # truncated final line after a crash
            rid = r["row_id"]
            if rid in recs:
                dup += 1
                # keep the later record: a resumed run supersedes a crashed one
            recs[rid] = r
    print(f"  loaded {len(recs)} unique rows from {len(files)} file(s)"
          + (f"  ({dup} duplicate row_id(s) superseded)" if dup else ""))
    return recs


def check_provenance(recs: dict[int, dict]) -> None:
    """Four machines must have run the SAME code against the SAME design."""
    # weights_sha / of describe the row->machine MAP rather than the physics.
    # A split on those is a different and nastier failure: the machines did not
    # partition the design, so some rows were computed twice and others by
    # nobody. That shows up at merge time as "rows MISSING", which reads like a
    # machine that never reported rather than the misconfiguration it is.
    for field in ("git_sha", "solver_sha256", "workbook_sha256", "design_sha256",
                  "weights_sha", "of"):
        vals = defaultdict(list)
        for r in recs.values():
            vals[r.get(field, "missing")].append(r.get("machine_id", "?"))
        if len(vals) > 1:
            print(f"  *** PROVENANCE SPLIT on {field}: results are NOT "
                  f"comparable across machines")
            for v, ms in vals.items():
                print(f"        {v}  <- {sorted(set(ms))}")


# Admissibility is applied HERE, not in the worker.  The worker records every
# gate quantity per row; this module decides what passes.  That separation is
# the whole reason the 2026-08-03 re-scoring of the stopped v2 campaign cost no
# CPU: moving delta_FP from 1e-6 to 1e-3 took its admissible count from 12 to 38
# without re-running a single row (plan S11(h)).  Never move a threshold into
# the worker -- it makes every past row unusable at the new threshold.
DFP_TOL_DEFAULT = 1e-2   # plan S11(h)


PILE_TOL = 0.05
TOPBIN_TOL = 0.02


def grid_ok(r: dict):
    """Re-derive the grid verdict from the recorded pile / top-bin quantities.

    Deliberately does NOT trust the row's stored `grid_limited`.  Rows written
    before 2026-08-03 had the withdrawn `occ > 0.10` rule folded into that flag,
    and on the stopped v2 campaign it alone condemned 71 of 275 rows that show
    no measurable truncation (plan S11(c)-2).  Recomputing here is what lets a
    superseded scoring rule be undone without re-running anything.

    Returns (ok, reason).  reason is None when ok.
    """
    if r.get("starved"):
        return False, "dose-starved"
    piles = {k: r.get(k) for k in ("pile_111", "pile_100", "pile_v")}
    # Fail closed: an axis we could not measure is not an axis we certified.
    missing = [k for k, v in piles.items() if v is None]
    if missing:
        return False, f"unmeasured:{','.join(sorted(missing))}"
    hot = [k for k, v in piles.items() if v > PILE_TOL]
    if hot:
        return False, f"pile>{PILE_TOL:g}:{','.join(sorted(hot))}"
    if r.get("equations") == "bin_moment":
        tb = {k: r.get(f"topbin_{k}") for k in ("111", "100", "v")}
        miss_tb = [k for k, v in tb.items() if v is None]
        if miss_tb or r.get("unmeasured_gates"):
            return False, "unmeasured:topbin"
        hot_tb = [k for k, v in tb.items() if v > TOPBIN_TOL]
        if hot_tb:
            return False, f"topbin>{TOPBIN_TOL:g}:{','.join(sorted(hot_tb))}"
    return True, None


def obs_value(r: dict, observable: str, at_dose=None):
    """The observable, read either at end-of-run or off the dose ladder.

    at_dose=None reads the end of the trajectory -- correct only when every row
    ran to the same dose.  Passing a rung reads `at_dose[rung]`, which is what
    makes a timed-out row comparable with a complete one instead of merely
    present alongside it.
    """
    if at_dose is None:
        return r.get(observable)
    lad = r.get("at_dose") or {}
    cell = lad.get(f"{at_dose:g}")
    return None if cell is None else cell.get(observable)


def dose_ladder_coverage(recs, design_row_ids):
    """rung -> how many of the design's rows reached it. Picks the common dose."""
    out = {}
    for rid in design_row_ids:
        r = recs.get(rid)
        if not r or r.get("solver_rc"):
            continue
        for rung in (r.get("at_dose") or {}):
            out[rung] = out.get(rung, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: float(kv[0])))


def usable(r: dict, dfp_tol: float = DFP_TOL_DEFAULT,
           require_converged: bool = False, at_dose=None) -> bool:
    """Does this row enter a Sobol index?

    DEFAULT POLICY since 2026-08-05 (author, plan S11(q)): a row is usable if it
    ran and reached the dose.  Truncation of the size tail does not disqualify
    it, because the campaign estimates a RANKING of parameters, not an absolute
    observable, and the ranking tolerates a truncated tail.

    require_converged=True restores the pre-2026-08-05 rule (grid_ok + delta_FP).
    Keep it working: re-estimating the indices both ways is the cheap empirical
    answer to "did truncation change the ranking", and it costs no CPU because
    every gate quantity is still recorded per row.  See --require-converged.
    """
    if not r or r.get("solver_rc"):
        return False
    if at_dose is not None:
        # A timed-out row is not rejected -- it is simply absent from the rungs
        # it never reached, pairwise, exactly like a missing AB row.
        if f"{at_dose:g}" not in (r.get("at_dose") or {}):
            return False
    elif r.get("starved"):
        # No common dose was requested, so this row's observables sit at a
        # DIFFERENT dose from everyone else's.  Refuse rather than silently
        # pool doses; --at-dose is the way to include it (plan S12(s)).
        return False
    if not require_converged:
        return True
    # Rows too old to carry the gate quantities cannot be re-scored; fall back
    # to the worker's verdict rather than guessing, and let the caller report
    # them separately instead of silently pooling two scoring rules.
    if r.get("pile_100") is None and r.get("grid_limited") is None:
        return bool(r.get("admissible"))
    ok, _ = grid_ok(r)
    if not ok:
        return False
    dfp = r.get("delta_FP")
    return dfp is not None and abs(dfp) < dfp_tol


def sobol_indices(recs, design_rows, p_keys, observable, n_boot=200, seed=0,
                  dfp_tol=None, require_converged=False, at_dose=None):
    """Saltelli/Jansen estimators with pairwise deletion.

    S_j    = mean_i f(B_i) [ f(AB_i^j) - f(A_i) ] / V        (Saltelli 2010)
    S_T_j  = mean_i [ f(A_i) - f(AB_i^j) ]^2 / (2 V)         (Jansen 1999)
    """
    # index the design: (base_idx, matrix, param_j) -> row_id
    idx = {}
    for d in design_rows:
        idx[(d["base_idx"], d["matrix"], d["param_j"])] = d["row_id"]
    bases = sorted({d["base_idx"] for d in design_rows})

    def val(b, m, j):
        rid = idx.get((b, m, j))
        if rid is None:
            return None
        r = recs.get(rid)
        if not usable(r, dfp_tol if dfp_tol is not None else DFP_TOL_DEFAULT,
                      require_converged, at_dose):
            return None
        v = obs_value(r, observable, at_dose)
        if v is None:
            return None
        v = float(v)
        return v if np.isfinite(v) else None

    fA, fB, ok_base = {}, {}, []
    for b in bases:
        a_, b_ = val(b, "A", -1), val(b, "B", -1)
        if a_ is not None and b_ is not None:
            fA[b], fB[b] = a_, b_
            ok_base.append(b)
    if len(ok_base) < 3:
        return None

    allf = np.array([fA[b] for b in ok_base] + [fB[b] for b in ok_base])
    V = float(np.var(allf, ddof=1))
    rng = np.random.default_rng(seed)
    res = {}
    for j, key in enumerate(p_keys):
        rows = [(fA[b], fB[b], val(b, "AB", j)) for b in ok_base]
        rows = [(a_, b_, ab) for (a_, b_, ab) in rows if ab is not None]
        n = len(rows)
        if n < 3 or V <= 0:
            res[key] = {"S": None, "ST": None, "n_eff": n,
                        "S_lo": None, "S_hi": None}
            continue
        A = np.array([r[0] for r in rows])
        B = np.array([r[1] for r in rows])
        AB = np.array([r[2] for r in rows])
        S = float(np.mean(B * (AB - A)) / V)
        ST = float(np.mean((A - AB) ** 2) / (2.0 * V))
        # bootstrap over base indices
        bs = []
        for _ in range(n_boot):
            k = rng.integers(0, n, n)
            Vb = np.var(np.concatenate([A[k], B[k]]), ddof=1)
            if Vb > 0:
                bs.append(np.mean(B[k] * (AB[k] - A[k])) / Vb)
        lo, hi = (float(np.percentile(bs, 2.5)),
                  float(np.percentile(bs, 97.5))) if bs else (None, None)
        res[key] = {"S": S, "ST": ST, "n_eff": n, "S_lo": lo, "S_hi": hi}
    return {"V": V, "n_base_ok": len(ok_base), "indices": res}


def main(argv=None):
    global DFP_TOL_DEFAULT
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=HERE / "results")
    ap.add_argument("--out", type=Path, default=HERE / "report")
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--dfp-tol", type=float, default=DFP_TOL_DEFAULT,
                    help="delta_FP admissibility bar, applied HERE not in the "
                         "worker, so it can be changed without re-running "
                         "anything (plan S11(h)). Only has an effect together "
                         "with --require-converged. Default 1e-2.")
    ap.add_argument("--require-converged", action="store_true",
                    help="restore the pre-2026-08-05 rule: reject rows whose "
                         "size tail is truncated (pile/top-bin) or whose "
                         "delta_FP exceeds --dfp-tol. OFF by default -- "
                         "truncation is recorded, not gated (plan S11(q)). Run "
                         "the report BOTH ways: if the parameter RANKING is the "
                         "same, truncation did not buy the ranking anything, "
                         "which is the claim the default policy rests on.")
    ap.add_argument("--at-dose", type=float, default=None,
                    help="screen every row at this COMMON dose, read off the "
                         "per-row `at_dose` ladder, instead of at the end of "
                         "each trajectory. This is what lets a timed-out row "
                         "contribute: it enters every rung it reached and is "
                         "absent above that, pairwise. Without it, rows that "
                         "stopped early are excluded, because their "
                         "observables sit at a different dose from everyone "
                         "else's and dose moves them far harder than any "
                         "numerical choice (plan S12(s)).")
    a = ap.parse_args(argv)
    DFP_TOL_DEFAULT = a.dfp_tol
    RC = a.require_converged
    AD = a.at_dose

    meta = json.loads(a.design.with_suffix(".meta.json").read_text(encoding="utf-8"))
    p_keys = meta["parameters"]
    lines = a.design.read_text(encoding="utf-8").strip().splitlines()
    cols = lines[0].split(",")
    design = []
    for ln in lines[1:]:
        v = dict(zip(cols, ln.split(",")))
        design.append({"row_id": int(v["row_id"]), "condition": v["condition"],
                       "base_idx": int(v["base_idx"]), "matrix": v["matrix"],
                       "param_j": int(v["param_j"])})

    recs = load_results(a.results)
    check_provenance(recs)

    total = len(design)
    have = sum(1 for d in design if d["row_id"] in recs)
    ok = sum(1 for d in design if usable(recs.get(d["row_id"]), a.dfp_tol, RC, AD))
    failed = sum(1 for d in design
                 if d["row_id"] in recs and recs[d["row_id"]].get("solver_rc"))
    inadm = have - ok - failed
    print(f"\n  policy: {'require-converged (pre-2026-08-05)' if RC else 'truncation RECORDED, NOT GATED (plan S11(q))'}"
          f"{f'; screened at COMMON dose {AD:g} dpa' if AD is not None else '; screened at END OF RUN'}")
    print(f"  coverage: {have}/{total} rows present, {ok} admissible, "
          f"{inadm} inadmissible, {failed} failed, {total-have} MISSING")

    # The dose ladder decides which common dose is worth screening at, so print
    # it before anything is estimated.  A rung nearly every row reaches costs
    # nothing to use; one reached by half of them halves n_eff everywhere.
    cov = dose_ladder_coverage(recs, [d["row_id"] for d in design])
    if cov:
        print(f"  dose ladder (rows reaching each rung, of {have} present):")
        print("    " + "  ".join(f"{k}:{v}" for k, v in cov.items()))
    n_starved = sum(1 for r in recs.values() if r.get("starved"))
    if AD is None and n_starved:
        if cov:
            best = max(float(k) for k, v in cov.items() if v == max(cov.values()))
            print(f"    *** {n_starved} timed-out row(s) EXCLUDED. Re-run with "
                  f"--at-dose {best:g} to include them at a common dose.")
        else:
            # Rows predating the dose ladder (plan S12(s)) carry no `at_dose`
            # block, so --at-dose cannot rescue them; only a re-run can.
            print(f"    *** {n_starved} timed-out row(s) EXCLUDED, and these rows "
                  f"carry NO dose ladder -- they predate plan S12(s), so "
                  f"--at-dose cannot recover them.")
    if RC:
        # The delta_FP bar is the one gate whose right value is a judgement
        # call, so show what the alternatives would buy BEFORE any index is
        # read.  A large spread here means the result depends on the threshold
        # and the threshold needs an argument, not a default.
        alts = [t for t in (1e-2, 1e-3, 1e-4, 1e-6) if t != a.dfp_tol]
        counts = "   ".join(
            f"{t:.0e}: {sum(1 for d in design if usable(recs.get(d['row_id']), t, True))}"
            for t in sorted([a.dfp_tol] + alts, reverse=True))
        print(f"  admissible vs delta_FP bar   {counts}      (in use: {a.dfp_tol:.0e})")
    else:
        # Truncation is not a gate, but it is not invisible either: say how much
        # of the pool carries it, so a ranking read off this report is read in
        # the knowledge of how much of it came from truncated rows.
        pool = [recs[d["row_id"]] for d in design
                if usable(recs.get(d["row_id"]), a.dfp_tol, False, AD)]
        if pool:
            trunc = sum(1 for r in pool if not grid_ok(r)[0])
            hi_dfp = sum(1 for r in pool if abs(r.get("delta_FP") or 0) >= a.dfp_tol)
            print(f"  of the {len(pool)} rows used: {trunc} ({100*trunc/len(pool):.0f}%) "
                  f"have a truncated size tail, {hi_dfp} ({100*hi_dfp/len(pool):.0f}%) "
                  f"have delta_FP >= {a.dfp_tol:g}")
            print(f"  -> re-run with --require-converged and compare the RANKING, "
                  f"not the values")
    if total - have:
        miss = sorted(d["row_id"] for d in design if d["row_id"] not in recs)
        by_mach = defaultdict(int)
        for m in miss:
            by_mach[m % 4] += 1
        print(f"    missing row_ids mod 4: {dict(by_mach)}  "
              f"(a whole residue class means a machine never reported)")

    # why runs were inadmissible -- this is diagnostic, not bookkeeping.
    # Reasons come from the SAME re-derivation usable() applies, so the counts
    # add up against the coverage line above; reading the row's stored flags
    # here would report the superseded rule.
    reasons = defaultdict(int)
    for r in recs.values():
        if r.get("solver_rc") or usable(r, a.dfp_tol, RC, AD):
            continue
        if not RC:
            # Under the default policy the only way to be inadmissible without
            # having failed is to have not reached the dose.
            reasons["dose-starved (use --at-dose)" if AD is None
                    else f"did not reach {AD:g} dpa"] += 1
            continue
        ok, why = grid_ok(r)
        if not ok:
            reasons[why] += 1
        elif abs(r.get("delta_FP") or 0) >= a.dfp_tol:
            reasons[f"delta_FP>={a.dfp_tol:g}"] += 1
    if reasons:
        print("    inadmissible because:")
        for k, v in sorted(reasons.items(), key=lambda t: -t[1]):
            print(f"       {v:5d}  {k}")

    a.out.mkdir(parents=True, exist_ok=True)

    # Restrict to observables validated in the mode that produced these rows.
    mode = equations_mode_of(recs)
    obs_list, blocked = observables_for_mode(mode or "discrete")
    if mode == "MIXED":
        print("\n  *** rows come from MORE THAN ONE equations mode - they are "
              "not poolable; separate them before estimating.")
    if blocked:
        print(f"\n  equations mode = {mode}: withholding {len(blocked)} "
              f"observable(s) whose closure is not validated in this mode:")
        print(f"    {', '.join(sorted(blocked))}")
        print("  Screening on them would produce indices about the numerics, "
              "not the physics.")
        print("  (<100> in bin_moment: mean_n_100 +144% vs a converged discrete "
              "reference, 2026-08-02.)")

    rows_out = []
    for cond in meta["conditions"]:
        dc = [d for d in design if d["condition"] == cond]
        for obs in obs_list:
            r = sobol_indices(recs, dc, p_keys, obs, n_boot=a.n_boot,
                              dfp_tol=a.dfp_tol, require_converged=RC,
                              at_dose=AD)
            if r is None:
                print(f"  {cond:>4s} {obs:16s} -- too few usable base points")
                continue
            print(f"\n  {cond} / {obs}   Var={r['V']:.4g}  "
                  f"base points usable {r['n_base_ok']}/{meta['N']}")
            ranked = sorted(r["indices"].items(),
                            key=lambda kv: -(kv[1]["ST"] or -1))
            for key, v in ranked[:6]:
                if v["S"] is None:
                    continue
                flag = "  <-- n_eff LOW" if v["n_eff"] < 0.5 * meta["N"] else ""
                print(f"      {key:18s} S={v['S']:7.3f} "
                      f"[{v['S_lo']:6.3f},{v['S_hi']:6.3f}]  "
                      f"ST={v['ST']:7.3f}  n_eff={v['n_eff']:3d}{flag}")
            for key, v in r["indices"].items():
                rows_out.append({"condition": cond, "observable": obs,
                                 "parameter": key, **v})

    csv = a.out / "T2_sobol_indices.csv"
    with csv.open("w", encoding="utf-8", newline="") as fh:
        fh.write("condition,observable,parameter,S,S_lo,S_hi,ST,n_eff\n")
        for r in rows_out:
            def g(k):
                return "" if r[k] is None else f"{r[k]:.6g}"
            fh.write(f"{r['condition']},{r['observable']},{r['parameter']},"
                     f"{g('S')},{g('S_lo')},{g('S_hi')},{g('ST')},{r['n_eff']}\n")
    print(f"\n  wrote {csv}  ({len(rows_out)} rows)")
    print("\n  Read n_eff before reading any index: an S_i^T computed from a")
    print("  handful of base points is not a screening result. Parameters whose")
    print("  n_eff collapsed are telling you their prior region is unintegrable")
    print("  at this tolerance -- itself a finding (plan Tier-2 failure handling).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
