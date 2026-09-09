"""
runs.py — the run manifest for the paper-revision verification study.

The ~25 runs of `docs/Formulation/paper_revision_verification_plan.md` as DATA,
not as two dozen hand-copied command lines.  One entry per table cell; the
runner (`campaign.py`) turns an entry into a `RadClusterSimulation` and a
timestamped output directory.

WHY A MANIFEST.  Three of the plan's rules are only enforceable if the runs are
enumerable in one place:

  * S3.2 -- every Table 1-3 column must share `t_span` and `n_points`, or
    "15.72 dpa" is a different point in different columns and the comparison
    drifts.  Here that is one constant (`GRID`) that every such entry inherits,
    rather than a number retyped 15 times.
  * S1.3 -- hold the bin ratio `r` at production and vary only `i_discrete`.
    The rungs below carry the plan's `I_bin`/`V_bin` for exactly that reason;
    `bin_layout()` then reads the REALISED layout back and refuses a silent
    fallback.
  * S3.4.1 -- no run that missed the comparison dose may be reported as a
    metric.  `dose_read` travels with the entry, so the reporting step can
    check it mechanically instead of by discipline.

Entries are ordered by the plan's execution order (S5): the Table 4 discrete
arm is the cost probe and comes first among the runs that must actually be
computed.
"""
from __future__ import annotations

# ── Shared reference configuration (S1.4, from the reference provenance) ────
# The anchor run these tables are read against.  Deviating from any of this in
# an entry is a deliberate act and shows up as an explicit key.
REFERENCE = "20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5"

# S3.2: ONE output grid for Tables 1-3.  15.72 dpa is point 35 of this grid.
# t_span end is dose/G with G = 1e-7 dpa/s, i.e. 40 dpa -> 4e8 s.
GRID = {"t_span_end_s": 4.0e8, "n_points": 37, "dose": 40.0}
DOSE_READ_MAIN = 15.72          # the dose Tables 1-3 are scored at
DOSE_READ_T4 = (0.2, 1.0, 2.0)  # S2.4: three doses, not one

BASE = {
    "I": 80000, "V": 20000,
    "equations": "bin_moment", "cascade": "fission",
    "shape_function": "linear",
    "solver_mode": "full_system",
    "i_mobile": 5, "v_mobile": 5,
    "rtol": 1e-5, "atol": 1e-20,
    "preconditioner": "woodbury",
    "he_kinetics": "quasi_steady_state",
    "lnl": 1,
    "omp_threads": 12,
    "dose": GRID["dose"],
    "n_points": GRID["n_points"],
    "dose_read": DOSE_READ_MAIN,
}

# TIMEOUT for the non-linear closures (P=1 constant, P=3 lognormal).  The plan
# says of these "may be unstable / may time out -- that is a result, not a
# failure" (S6), but the default 86400 s means collecting that result costs the
# claiming machine a full day, during which it takes nothing else off the board.
#
# Measured 2026-09-06: T4_C4_P3_4k (lognormal, 3 moments/bin) ran 66 minutes and
# reached 7.1e-4 of 2 dpa -- 0.04% -- while the SAME rung under the linear
# closure finished in 54 s.  Thirty minutes is already 33x the linear cost, so a
# run still under a percent at that point has demonstrated what the table needs
# to record.  The cap turns a 24 h occupancy into a 30 min measurement of the
# same fact; the record carries dose_reached and starved either way.
TIMEOUT_NONLINEAR_CLOSURE = 1800

# Production rung (= the reference run).  Its grid recurs in several tables.
B4_GRID = {"i_discrete": 100, "I_bin": 18, "v_discrete": 5, "V_bin": 20}


# P = 3 (lognormal, 3 moments per bin) DOES NOT ADVANCE on this model.  Measured
# 2026-09-06 on both domains: T4_C4_P3_4k reached 0.0022 of 2 dpa, T2_P3 reached
# 0.00054 of 40 dpa with a single output step running 2.9 h.  Plan S6 anticipates
# exactly this ("P = 3 may time out -- that is a result"), so the two rows stay
# in the manifest as "did not reach" entries with the dose they achieved; they
# are simply never CLAIMED again.  Re-enable by clearing this field once the
# lognormal closure's stiffness is addressed.
# The TD pair was the decisive test for a path difference that is now FIXED
# (cfce995).  check_bm_vs_discrete.py re-asks it on every run in 20 s, so the
# pair would spend ~2 days confirming what a 20 s test already shows.  Domain
# and dose coverage is bought instead by raising the test's own grid.
TD_SUPERSEDED = ("superseded by cfce995 + check_bm_vs_discrete.py; the pair "
                 "cost ~2 days to re-confirm a 20 s result")

P3_EXCLUDED = ("P=3 lognormal does not advance; measured 2026-09-06 on both "
               "domains. Recorded as 'did not reach' per plan S3.4.1.")


def _r(run_id, table, label, notes="", **over):
    """One manifest entry: BASE, overridden, plus identity and provenance."""
    e = dict(BASE)
    e.update(over)
    e.update({"run_id": run_id, "table": table, "label": label, "notes": notes})
    return e


def manifest() -> list[dict]:
    """The full study, in execution order (plan S5)."""
    runs: list[dict] = []

    # ── Table 4 — binning verification against the exact solution ───────────
    # I = V = 4000, MODEL-TO-MODEL.  Runs FIRST: the discrete arm is the only
    # genuinely uncertain cost (plan S5.2).
    #
    # WHY 4000 AND NOT THE PLAN'S I=8000/V=20000.  The plan's revision note
    # rejected a reduced domain because "V = 1000 sits below the mean cavity
    # (1071 vacancies at 15.72 dpa), so d_void reads the ceiling, not the
    # model".  That is a VALIDATION criterion and Table 4 is not validation:
    # both arms solve the same equations on the same domain, so a truncation
    # affects them identically and cancels out of the closure error being
    # measured.  It limits what may be CLAIMED (nothing in Table 4 may be
    # compared to experiment) rather than whether the measurement is sound.
    # Tables 1-3 and 5, which ARE scored against experiment, stay at the
    # production domain for exactly that reason.
    #
    # Measured, not assumed: the plan's own I=8000/V=20000 arm was launched
    # 2026-09-06 and reached 2.2e-3 of 2 dpa in 14 minutes, with per-step cost
    # growing ~3.3x per output step -- an extrapolated 12-30 days.  The cost is
    # in V, not I: N_eq = 8000 SIA + 20000 VAC + 6, plus an appended 8000
    # <100> block, and the vacancy block is 71% of the main system.  The plan's
    # S6 fallback ("cut I to 4000") removes ~22% of the equations where ~100x
    # was needed.  The existing G5b run (discrete, I=4000, V=4000) reached
    # 6.06 dpa in 6.8 h SINGLE-THREADED, which brackets this table's doses.
    T4 = {"I": 4000, "V": 4000, "dose": 2.0, "n_points": 37,
          "dose_read": DOSE_READ_T4}
    runs.append(_r("T4_D_4k", "T4", "D (exact)",
                   "Fully discrete: the exact solution of the same equations. "
                   "Cost probe -- run first (plan S5.2).",
                   equations="discrete", **T4))
    # The ladder holds r at the production value and varies only the resolution
    # (plan S1.3).  I_bin/V_bin below are ln(domain/discrete)/ln(r) rounded, so
    # r_i stays 1.44-1.47 and r_v 1.51-1.53 across every rung.
    for rid, lbl, ratio, g in (
        ("T4_C1_4k", "C1", 10,  {"i_discrete": 400, "I_bin": 6,  "v_discrete": 400, "V_bin": 6}),
        ("T4_C2_4k", "C2", 40,  {"i_discrete": 100, "I_bin": 10, "v_discrete": 100, "V_bin": 9}),
        ("T4_C3_4k", "C3", 160, {"i_discrete": 25,  "I_bin": 14, "v_discrete": 25,  "V_bin": 12}),
        ("T4_C4_4k", "C4", 800, {"i_discrete": 5,   "I_bin": 18, "v_discrete": 5,   "V_bin": 16}),
    ):
        note = f"resolution 1/{ratio} on both axes"
        if rid == "T4_C4_4k":
            note = ("SIA axis is EXACTLY production: i_discrete/I = 5/4000 = "
                    "1/800, and holding r_i returns I_bin = 18, the production "
                    "value.  The vacancy axis is 1/800, FINER than production's "
                    "5/20000 = 1/4000 -- matching it would need v_discrete = 1, "
                    "below v_mobile = 5.  So this column if anything UNDERSTATES "
                    "the production vacancy closure error; state that rather "
                    "than claiming an exact two-axis analogue.")
        runs.append(_r(rid, "T4", lbl, note, **T4, **g))
    # S2.4: closure ORDER scored against the exact answer, not against B4.
    C4 = {"i_discrete": 5, "I_bin": 18, "v_discrete": 5, "V_bin": 16}
    runs.append(_r("T4_C4_P1_4k", "T4", "C4, P=1 constant",
                   "Closure order scored against the exact column.",
                   shape_function="constant", timeout_s=TIMEOUT_NONLINEAR_CLOSURE,
                   **T4, **C4))
    runs.append(_r("T4_C4_P3_4k", "T4", "C4, P=3 lognormal",
                   "Closure order scored against the exact column.  RAN and "
                   "did not advance: 0.0022 of 2 dpa in 1866 s, all three "
                   "scoring rungs MISSING.  Kept as a 'did not reach' row.",
                   shape_function="lognormal", excluded=P3_EXCLUDED,
                   timeout_s=TIMEOUT_NONLINEAR_CLOSURE,
                   **T4, **C4))

    # ── DIAGNOSTIC PAIR — is the gap the closure, or the two RHS paths? ─────
    # Table 4 shows N_111 flat at ~3.6e21 across i_discrete 400 -> 5 (an 80x
    # change) while the exact arm gives 1.26e22.  C1 integrates sizes 1-400
    # individually and the exact distribution holds 89% of its loops in n =
    # 21-100, entirely inside that discrete region -- so a CLOSURE error cannot
    # explain it.  Either the binned large-size tail (77% of C4's content sits
    # above n = 500, against 0.1% in the exact solution) starves the small
    # loops through the coupled dynamics, or the two RHS implementations do not
    # solve the same equations.
    #
    # CLAUDE.md S1: "When I_bin = 0 and i_discrete = I, all equations are
    # discrete -> recovers full_CD".  So BM_FULL is the bin_moment CODE PATH
    # carrying NO closure at all.  Run against DISC_REF on an identical grid:
    #
    #   BM_FULL == DISC_REF  -> the paths agree; the fault is the closure/tail
    #   BM_FULL != DISC_REF  -> the paths differ; that is a solver defect and
    #                           every bin_moment result in the study inherits it
    #
    # 0.2 dpa, not 2: the gap is already wide by ~1e-3 dpa (N_111 3.36e21 vs
    # 2.30e21) and this costs ~1.5 h instead of ~17 h.  Both share t_span and
    # n_points so no interpolation is needed to compare them.
    DIAG = {"I": 4000, "V": 4000, "dose": 0.2, "n_points": 37,
            "dose_read": (0.05, 0.1, 0.2), "table": "TD"}
    runs.append(_r("TD_BM_FULL", "TD", "bin_moment, i_discrete = I",
                   "The bin_moment code path with NO closure (I_bin = V_bin = 0). "
                   "SUPERSEDED: it was queued to find the path difference, which "
                   "cfce995 then found and fixed, and "
                   "codes/Python_Testing/check_bm_vs_discrete.py now settles the "
                   "same question in 20 s to 7e-8 rather than in ~2 days. Stopped "
                   "at 0.063 of 0.2 dpa with steps at 148 min and lengthening.",
                   excluded=TD_SUPERSEDED,
                   equations="bin_moment", i_discrete=4000, I_bin=0,
                   v_discrete=4000, V_bin=0,
                   **{k: v for k, v in DIAG.items() if k != "table"}))
    runs.append(_r("TD_DISC_REF", "TD", "discrete, same grid",
                   "Matched-grid discrete reference for TD_BM_FULL. SUPERSEDED "
                   "with its pair: half a comparison proves nothing.",
                   excluded=TD_SUPERSEDED,
                   equations="discrete",
                   **{k: v for k, v in DIAG.items() if k != "table"}))

    # ── Table 4-NC — closure verification WITHOUT loop coarsening ───────────
    # Salvages the 17.6 h discrete arm of 2026-09-07
    # (output/20260907_093333_T4_D_4k_...), which ran before cfce995 and so had
    # loop_coal silently absent.  Re-running the closure rungs with LOOP_COAL = 0
    # puts both arms back on the same equations, and the comparison becomes a
    # legitimate verification again -- at ~2 min a rung instead of days.
    #
    # WHY THIS IS SOUND.  Verification asks whether the closure reproduces the
    # exact solution OF THE SAME EQUATIONS.  It does not require those equations
    # to be physically right; that is validation, and plan S1.2 separates the
    # two deliberately.  Bin partitioning, moment closure and intra-bin shape --
    # the approximations the reviewer actually asked about -- are all exercised
    # here.
    #
    # WHAT MUST BE STATED.  The no-coarsening variant is known-unphysical:
    # rate_kernels.cpp:1930 says without the channel "the 1/2<111> mean size
    # pins AT the cutoff" and "the density runs 20-60x over the EUROFER97 data",
    # which is exactly what the 2026-09-07 arm shows (d_111 frozen at 1.76 nm,
    # N_111 3.5x the closure's).  So S1.3's TRANSFER argument weakens: the
    # closure error measured here is not demonstrably the production run's
    # closure error, because production runs WITH coarsening and that changes
    # the distribution the closure has to represent.  Say so in the paper rather
    # than letting a reader assume otherwise.
    NC = {"I": 4000, "V": 4000, "dose": 2.0, "n_points": 37,
          "dose_read": (2.0,), "loop_coal": 0}
    for rid, lbl, g in (
        # C0 is finer than C1 -- N_eq 3214 against 830 -- so the size-effect
        # figures have a rung ABOVE the previous finest and the trend is not
        # read off an endpoint.
        ("T4_C0_nc", "C0 (no coarsening)", {"i_discrete": 1600, "I_bin": 2, "v_discrete": 1600, "V_bin": 2}),
        # C0b is the insurance rung: still finer than C1 (N_eq 1614 vs 830) but
        # ~4x cheaper than C0, because the bin_moment discrete-discrete
        # coalescence block is O(i_discrete^2).  Whichever of C0/C0b lands is
        # the "finer than C1" point in the size-effect figures.
        ("T4_C0b_nc", "C0b (no coarsening)", {"i_discrete": 800, "I_bin": 4, "v_discrete": 800, "V_bin": 4}),
        # C0c: the smallest step above C1 that still satisfies "more equations
        # than C1" (1226 vs 830).  Cost in this variant grows far faster than
        # O(i_discrete^2) -- C1 ran in 12 s, C0b was at 1.1% after 7 min --
        # because with LOOP_COAL = 0 the frozen spectrum is stiffer the more of
        # it is resolved discretely rather than smoothed into bins.
        ("T4_C0c_nc", "C0c (no coarsening)", {"i_discrete": 600, "I_bin": 5, "v_discrete": 600, "V_bin": 5}),
        ("T4_C1_nc", "C1 (no coarsening)", {"i_discrete": 400, "I_bin": 6,  "v_discrete": 400, "V_bin": 6}),
        ("T4_C2_nc", "C2 (no coarsening)", {"i_discrete": 100, "I_bin": 10, "v_discrete": 100, "V_bin": 9}),
        ("T4_C3_nc", "C3 (no coarsening)", {"i_discrete": 25,  "I_bin": 14, "v_discrete": 25,  "V_bin": 12}),
        ("T4_C4_nc", "C4 (no coarsening)", {"i_discrete": 5,   "I_bin": 18, "v_discrete": 5,   "V_bin": 16}),
    ):
        runs.append(_r(rid, "T4NC", lbl,
                       "Scored against the 2026-09-07 discrete arm at 2.0 dpa.",
                       **NC, **g))
    runs.append(_r("T4_C4_P1_nc", "T4NC", "C4, P=1 (no coarsening)",
                   "Closure ORDER against the same discrete arm.",
                   shape_function="constant", **NC,
                   **{"i_discrete": 5, "I_bin": 18, "v_discrete": 5, "V_bin": 16}))

    # ── Table 1 — closure convergence at the production domain ──────────────
    # Hold r at production, vary only i_discrete (S1.3).  B4 already exists as
    # the reference run and is not recomputed.
    for rid, lbl, g in (
        ("T1_B1", "B1", {"i_discrete": 6400, "I_bin": 7, "v_discrete": 1600, "V_bin": 6}),
        ("T1_B2", "B2", {"i_discrete": 1600, "I_bin": 11, "v_discrete": 400, "V_bin": 9}),
        ("T1_B3", "B3", {"i_discrete": 400, "I_bin": 14, "v_discrete": 100, "V_bin": 13}),
    ):
        runs.append(_r(rid, "T1", lbl,
                       "B1 is the most-resolved rung; deviations in Table 1 are "
                       "quoted relative to it." if rid == "T1_B1" else "", **g))

    # ── Table 2 — intra-bin closure and tolerance, ONE knob at a time ───────
    # Deviations relative to B4, since these sit inside the closure rather than
    # on the discrete/binned axis.
    runs.append(_r("T2_P1", "T2", "P=1 constant",
                   "May be unstable -- that is a result, not a failure (S6).",
                   shape_function="constant", timeout_s=TIMEOUT_NONLINEAR_CLOSURE,
                   **B4_GRID))
    runs.append(_r("T2_P3", "T2", "P=3 lognormal",
                   "May time out -- that is a result, not a failure (S6).  It "
                   "did: 0.00054 of 40 dpa in 174 min, 2.9 h on a single output "
                   "step, stopped by hand rather than left to burn a core to "
                   "the 24 h cap.",
                   shape_function="lognormal", excluded=P3_EXCLUDED,
                   timeout_s=TIMEOUT_NONLINEAR_CLOSURE,
                   **B4_GRID))
    runs.append(_r("T2_IBIN10", "T2", "I_bin 10",
                   **{**B4_GRID, "I_bin": 10}))
    runs.append(_r("T2_IBIN40", "T2", "I_bin 40",
                   **{**B4_GRID, "I_bin": 40}))
    runs.append(_r("T2_RTOL4", "T2", "rtol 1e-4", rtol=1e-4, **B4_GRID))
    runs.append(_r("T2_RTOL6", "T2", "rtol 1e-6", rtol=1e-6, **B4_GRID))

    # ── Table 3 — helium ───────────────────────────────────────────────────
    # Model-to-model; no experimental column.  The fusion arms are NOT optional
    # (S2.3.2): a He block evaluated only at fission is flat across all six
    # observables and reads as evasion.  The flatness is the result, but only
    # with the fusion columns beside it.
    # CONTROL COLUMN (added 2026-09-07).  The plan writes Table 3's first
    # column as "fission QSS (= B4)" and S1.3 says B4 "already exists as the
    # reference run and is not recomputed".  That makes the one QSS/dynamic
    # contrast in the study uncontrolled: the reference was produced by a
    # different code state on a different output grid, so subtracting it from
    # T3_FISS_DYN measures the code and the grid as well as the helium closure.
    #
    # This row is the controlled twin.  It differs from T3_FISS_DYN in exactly
    # one key -- he_kinetics -- and BASE already supplies cascade="fission" and
    # he_kinetics="quasi_steady_state", so the difference is the closure and
    # nothing else.  It is also, by construction, the production rung, so it
    # doubles as the B4 column Tables 1 and 2 quote.
    runs.append(_r("T3_FISS_QSS", "T3", "fission QSS (control)",
                   "Controlled twin of T3_FISS_DYN: identical in every respect "
                   "except he_kinetics.  Replaces the plan's uncomputed "
                   "'= B4' column, which was not a controlled comparator.",
                   **B4_GRID))
    runs.append(_r("T3_FISS_DYN", "T3", "fission dynamic",
                   he_kinetics="dynamic", **B4_GRID))
    runs.append(_r("T3_FUS_QSS", "T3", "fusion QSS",
                   cascade="fusion", **B4_GRID))
    runs.append(_r("T3_FUS_DYN", "T3", "fusion dynamic",
                   cascade="fusion", he_kinetics="dynamic", **B4_GRID))
    # The last two Table 3 columns need the S3.1 he_model unwelding, which this
    # pass does NOT implement.  They are carried here, blocked, so the manifest
    # is the whole study rather than the runnable part of it.
    runs.append(_r("T3_FUS_CASE1", "T3", "fusion Case 1",
                   "BLOCKED on plan S3.1 (unweld he_model from cascade).",
                   cascade="fusion", he_model="case1", blocked_on="S3.1", **B4_GRID))
    runs.append(_r("T3_FUS_CASE2", "T3", "fusion Case 2",
                   "BLOCKED on plan S3.1 (unweld he_model from cascade).",
                   cascade="fusion", he_model="case2", blocked_on="S3.1", **B4_GRID))

    # ── Table 5 — null knobs (the direct reply to the reviewer) ─────────────
    # Performance-only knobs get a NULL TEST, not a curve (S3.4.4).  If any of
    # these moves a converged answer, that is a bug report.
    runs.append(_r("T5_WOODBURY", "T5", "Woodbury",
                   "Null test: preconditioner enters no rate kernel.",
                   preconditioner="woodbury", **B4_GRID))
    runs.append(_r("T5_JACOBI", "T5", "Jacobi",
                   "Null test: must agree with Woodbury to many sig figs.",
                   preconditioner="jacobi", **B4_GRID))
    runs.append(_r("T5_THREADS1", "T5", "1 thread",
                   "Null test: thread count must not move the answer.",
                   omp_threads=1, **B4_GRID))
    runs.append(_r("T5_THREADS12", "T5", "12 threads",
                   "Null test: paired with T5_THREADS1.",
                   omp_threads=12, **B4_GRID))

    _check_unique(runs)
    return runs


def _check_unique(runs):
    seen = {}
    for e in runs:
        if e["run_id"] in seen:
            raise ValueError(f"duplicate run_id {e['run_id']!r} in manifest")
        seen[e["run_id"]] = e


# TABLE 5 IS NOT CLAIMED BY DEFAULT.  Its columns measure WALL-CLOCK as the
# dependent variable -- Woodbury vs Jacobi, 1 thread vs 12 -- so a T5 run that
# shares a machine with other solves measures the contention, not the knob.
# Plan S3.4.4: "performance-only knobs get a null test, not a curve"; a null
# test taken under load is worthless.  It happened: the unrestricted loop
# finished the exact arm, moved on to the next free run, and started
# T5_WOODBURY alongside three other solvers.  Claim these only with an explicit
# `--only T5` on an otherwise idle machine.
TIMING_SENSITIVE_TABLES = {"T5"}


def runnable(runs=None, tables=None) -> list[dict]:
    """Entries a worker may CLAIM: not blocked on unimplemented code, and not
    excluded by a measured failure.  Excluded entries stay in manifest() so the
    board still displays their recorded outcome."""
    want = set(tables or ())
    return [e for e in (runs or manifest())
            if not e.get("blocked_on") and not e.get("excluded")
            and (e["table"] not in TIMING_SENSITIVE_TABLES or e["table"] in want)]


if __name__ == "__main__":
    import json, sys
    m = manifest()
    if "--json" in sys.argv:
        print(json.dumps(m, indent=2))
    else:
        blocked = [e for e in m if e.get("blocked_on")]
        print(f"{len(m)} runs ({len(blocked)} blocked, {len(m) - len(blocked)} runnable)\n")
        hdr = f"{'run_id':16s} {'tbl':4s} {'label':22s} {'eqs':9s} " \
              f"{'I':>6s} {'V':>6s} {'i_d':>5s} {'I_b':>4s} {'v_d':>5s} {'V_b':>4s} {'dose':>5s}"
        print(hdr); print("-" * len(hdr))
        for e in m:
            print(f"{e['run_id']:16s} {e['table']:4s} {e['label'][:22]:22s} "
                  f"{e['equations']:9s} {e['I']:6d} {e['V']:6d} "
                  f"{e.get('i_discrete', 0):5d} {e.get('I_bin', 0):4d} "
                  f"{e.get('v_discrete', 0):5d} {e.get('V_bin', 0):4d} "
                  f"{e['dose']:5.1f}" + ("   BLOCKED" if e.get("blocked_on") else ""))
