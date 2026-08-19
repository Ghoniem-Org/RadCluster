#!/usr/bin/env python3
"""Propose the next calibration stage for whatever machine this is.

READS  calibration_ledger.json (written by learn.py)
WRITES design/<stage>_calib.csv, design/<stage>_labels.json, and the `next`
       block of the ledger + CALIBRATION_GUIDE.md

IT DOES NOT LAUNCH ANYTHING.  It prints the exact run_ensemble command, sized
to this machine slot count and row budget, and stops.  A design nobody read
can burn 20 h of a machine time, and the machine is usually one the author is
not sitting at.

HOW IT CHOOSES.  Greedy, and deliberately simple enough to argue with:

  1. Pin every parameter to the BEST VALID ROW in the ledger.  That row is the
     only point known to reach full dose grid-clean, so it is the only safe
     place to stand while varying something.
  2. Rank the residuals |log10(model/target)| over the six observables.
  3. Choose levers for the biggest residual, preferring in this order:
        never-varied  >  inconclusive  >  live
     and never proposing a lever the ledger calls `dead`.
  4. Size the factorial to the machine slot count so every row runs
     concurrently, and refuse if the estimated cost exceeds the row budget.

The lever->observable map below is the ONE piece of physics judgement in this
file.  It is explicit so it can be corrected as the campaign learns, rather
than buried in a heuristic.

Usage:
    python plan.py                       # propose for this machine
    python plan.py --machine 0           # propose as if machine 0
    python plan.py --levers eta,f_cl_i   # override the choice
    python plan.py --stage S15           # override the name
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
LEDGER = HERE / "calibration_ledger.json"
MACHINES = HERE / "machines.json"

# ---------------------------------------------------------------------------
# WHICH LEVERS PLAUSIBLY MOVE WHICH OBSERVABLE.
#
# Physics judgement, stated openly.  `candidates` are ordered by expected
# strength.  A lever is only ever proposed if the ledger has not called it
# dead, so this table degrades gracefully as verdicts accumulate.
# ---------------------------------------------------------------------------
LEVER_MAP = {
    # Too many loops / too much SIA locked in loops -> attack production and
    # the clustered fraction born in cascades before touching sink geometry.
    "N_100": ["eta", "f_cl_i", "E_b_i2", "i_mobile", "s_i"],
    "N_111": ["eta", "f_cl_i", "E_b_i2", "B_111", "i_mobile"],
    # Loops too small -> growth vs nucleation balance, then mobility.
    "d_100": ["f_cl_i", "eta", "E_b_i2", "L_hat", "phi_max_junc"],
    "d_111": ["B_111", "E_b_i2", "f_cl_i", "E_m_i", "L_hat"],
    # Cavities absent -> vacancy clustering, surface energy, He.
    "N_void": ["f_cl_v", "gamma_s", "E_b_v2", "E_b_hV_1", "s_v"],
    "d_void": ["gamma_s", "E_b_v2", "f_cl_v", "E_m_v", "E_b_hV_1"],
}

# SIGN of d(observable)/d(lever), where it is unambiguous.  Used to spend rows
# on the corrective side instead of bracketing symmetrically: N_100 is 12x too
# HIGH and eta raises it, so eta is swept DOWNWARD from the base rather than
# across the whole box.  Omit a pair whose sign is genuinely uncertain -- a
# missing entry falls back to a symmetric bracket, which is never wrong, only
# less efficient.
DIRECTION = {
    ("eta", "N_100"): +1, ("eta", "N_111"): +1,      # more production, more nuclei
    ("eta", "d_100"): -1, ("eta", "d_111"): -1,      # ... spread over more loops
    ("f_cl_i", "N_100"): +1, ("f_cl_i", "N_111"): +1,  # born clustered = born nucleated
    ("f_cl_i", "d_100"): -1, ("f_cl_i", "d_111"): -1,
    ("E_b_i2", "N_100"): +1, ("E_b_i2", "N_111"): +1,  # small clusters survive
    ("E_b_i2", "d_100"): -1, ("E_b_i2", "d_111"): -1,
    ("B_111", "d_111"): +1,
    ("f_cl_v", "N_void"): +1,
    # V1 measured both of these directly, one lever at a time:
    #   gamma_s 3.0 -> 1.5   N_void 0.00036x -> 47860x   (lower capillary
    #                        penalty, small cavities survive)
    #   E_m_v   0.50 -> 0.90 N_void 0.165x   -> 6793x    (slower vacancies,
    #                        less recombination before they cluster)
    ("gamma_s", "N_void"): -1, ("gamma_s", "d_void"): -1,
    ("E_m_v", "N_void"): +1, ("E_m_v", "d_void"): +1,
}

# Prior box.  A design must never leave this; values are physical bounds, not
# preferences.  (lo, hi, kind) where kind picks linear or log spacing.
PRIOR_BOX = {
    "eta":          (0.05, 0.60, "log"),
    "f_cl_i":       (0.02, 0.50, "log"),
    "f_cl_v":       (0.05, 0.80, "linear"),
    "E_b_i2":       (0.40, 1.20, "linear"),
    "E_b_v2":       (0.10, 0.45, "linear"),
    "B_111":        (0.20, 0.60, "linear"),
    "L_hat":        (10.0, 3000.0, "log"),
    "i_mobile":     (5, 60, "linear"),
    "gamma_s":      (1.5, 3.0, "linear"),
    "s_i":          (1.0, 3.0, "linear"),
    "s_v":          (1.0, 3.0, "linear"),
    "E_m_i":        (0.25, 0.55, "linear"),
    "E_m_v":        (0.50, 0.90, "linear"),
    "phi_max_junc": (0.05, 0.60, "linear"),
    "E_b_hV_1":     (1.0, 3.0, "linear"),
}

# NUMERICAL-VALIDITY CEILINGS, distinct from the physical prior box above.
# A value can be perfectly physical and still make the run unusable -- not by
# timing out (that is `affordability`, learned from results) but by driving the
# solver grid into pile-up, which invalidates the whole row including the
# observables the stage was actually about.  Entries are curated, and each one
# cites the evidence.  Intersected with PRIOR_BOX by box_for().
VALIDITY_BOX = {
    # S14: f_cl_v = 0.5 ran grid-clean to 15 dpa (occ_100 = 0.34); f_cl_v = 0.7
    # piled at the ceiling (occ_100 = 0.90, grid_limited).  The boundary inside
    # (0.5, 0.7) is unresolved, so 0.5 is the highest value known to work.
    # NOTE the caveat: S14 moved f_cl_v, E_b_v2 and s_v together, so the
    # pile-up is not attributable to f_cl_v alone.  The cap is conservative for
    # that reason and can be lifted once a stage varies f_cl_v on its own.
    "f_cl_v": (None, 0.5),
}

# A lever whose corrective direction retains less than this fraction of its
# prior box -- after affordability clipping -- is not worth a stage.
MIN_USABLE_SPAN = 0.10

# Columns the solver reads as counts, not reals.  A design that writes
# i_mobile = 32.5 is not a smaller step than 33 -- it is an invalid input.
INTEGER_COLS = {"i_mobile", "v_mobile", "i_discrete", "v_discrete"}

BOOKKEEPING = {"row_id", "condition", "cond_row_id", "matrix", "base_idx",
               "param_j", "theta_id"}

# Levers that are COUNTS, not continuous parameters.  run_ensemble takes
# i_mobile through int() at construction (line ~413), so an evenly spaced
# float grid would run 32 while the design recorded 32.5 -- the row's own
# provenance would then disagree with what was computed, and learn.py would
# pair rows on a value no solver ever saw.  Round here instead.
INTEGER_LEVERS = {"i_mobile"}


def load_ledger() -> dict:
    if not LEDGER.exists():
        raise SystemExit("no calibration_ledger.json -- run `python learn.py` first")
    return json.loads(LEDGER.read_text(encoding="utf-8"))


def detect(machine_index=None) -> tuple:
    reg = json.loads(MACHINES.read_text(encoding="utf-8"))
    if machine_index is not None:
        for m in reg["machines"]:
            if m["index"] == machine_index:
                return m, reg
        raise SystemExit("no machine with index %s in machines.json" % machine_index)
    sys.path.insert(0, str(HERE))
    import run_ensemble as RE
    return RE.detect_machine(reg), reg


def build_commands(stage: str, machine: dict, reg: dict, slots: int,
                   budget_h: float, share: list, grid_v: int = 2000,
                   grid_vbin: int = 20) -> list:
    """The run_ensemble invocation(s) for this stage.

    THE SHARD INDEX IS NOT THE REGISTRY INDEX.  run_ensemble assigns a row with
    `assign_machine(row_id, of) == row_id % of`, so with `--of 1` only
    `--machine 0` is ever given any rows.  Passing this host registry index
    (2 for the workstation) with `--of 1` silently assigns it ZERO rows and the
    run exits having done nothing.  Shard index is therefore always the
    position within the participating set, while the registry index is used
    only to name the output file so two machines cannot collide.
    """
    design_rel = "design/%s_calib.csv" % stage

    def one(shard: int, of: int, reg_index: int, workers: int) -> list:
        return [
            "python run_ensemble.py \\",
            "    --design %s \\" % design_rel,
            "    --conditions conditions_S8.json \\",
            "    --spec parameters_S4.json \\",
            "    --out results/%s_calib_machine%d.jsonl \\" % (stage, reg_index),
            "    --machine %d --of %d \\" % (shard, of),
            "    --equations bin_moment --i-discrete 100 --i-bin 36 \\",
            "    --v-discrete 5 --v-bin %d --allow-mixed \\" % grid_vbin,
            "    --I 80000 --V %d --dose 15.0 --lnl 1 --rtol 1e-5 \\" % grid_v,
            "    --solver-mode full_system \\",
            "    --timeout-s %d --workers %d --omp-threads 1"
            % (int(budget_h * 3600), workers),
        ]

    if not share:
        return one(0, 1, machine["index"], slots)

    lines = []
    for shard, reg_index in enumerate(share):
        m = next(x for x in reg["machines"] if x["index"] == reg_index)
        lines.append("# shard %d/%d -- run this ON %s"
                     % (shard, len(share), m["name"]))
        lines.extend(one(shard, len(share), reg_index, m.get("slots", 1)))
        lines.append("")
    return lines


def residual_rank(led: dict) -> tuple:
    """Observables ordered by how badly they miss, worst first.

    Returns (chaseable, deferred).  Deferred observables are still ranked and
    reported -- they are just not allowed to select levers.  See
    `policy.deprioritized_observables` in the ledger for why each is deferred.
    """
    r = led.get("residuals") or {}
    dep = (led.get("policy") or {}).get("deprioritized_observables") or {}
    out = [(k, abs(math.log10(v))) for k, v in r.items() if v]
    out.sort(key=lambda kv: -kv[1])
    return ([kv for kv in out if kv[0] not in dep],
            [kv for kv in out if kv[0] in dep])


def pick_levers(led: dict, worst: list, n_wanted: int = 3) -> tuple:
    """Choose levers for the worst residual, honouring the ledger verdicts."""
    verdict = {}
    for name, L in led["levers"].items():
        for c in L["columns"]:
            # A column inherits the strongest claim made about any group it was
            # in: live beats inconclusive beats dead is NOT the order -- dead is
            # only earned by a group that isolated it, so a column that is dead
            # in one group and live in another is treated as live.
            cur = verdict.get(c)
            rank = {"live": 3, "inconclusive": 2, "dead": 1}
            if cur is None or rank[L["verdict"]] > rank[cur]:
                verdict[c] = L["verdict"]
    never = set(led.get("untested_columns") or [])
    unwired = set(((led.get("policy") or {}).get("unwired_parameters") or {}))

    def priority(col):
        if col in unwired:
            return 99                     # physics never reads it -- see policy
        if col in never:
            return 0                      # never varied -- most information
        v = verdict.get(col)
        if v == "inconclusive":
            return 1
        if v == "live":
            return 2
        return 99                         # dead -- never propose

    params = (led.get("best") or {}).get("params") or {}
    ratios = led.get("residuals") or {}
    chosen, why, skipped = [], [], []
    for obs, miss in worst:
        for col in LEVER_MAP.get(obs, []):
            if col in chosen or priority(col) == 99:
                continue
            if col not in PRIOR_BOX:
                continue
            try:
                base = float(params.get(col))
            except (TypeError, ValueError):
                base = sum(PRIOR_BOX[col][:2]) / 2.0
            pref = prefer_for(col, obs, ratios.get(obs))
            span = usable_span(col, led, base, pref)
            if span < MIN_USABLE_SPAN:
                skipped.append((col, obs, span))
                continue
            chosen.append(col)
            why.append((col, obs, miss, priority(col)))
            if len(chosen) >= n_wanted:
                return chosen, why, skipped
    return chosen, why, skipped


def box_for(col: str, led: dict, base: float):
    """The prior box for `col`, clipped to the region measured as affordable.

    THE CLIFF LOCATION IS UNKNOWN.  S14 shows E_b_i2 = 0.60 never reaches dose
    and 0.75 always does; the boundary is somewhere in between and nothing
    measures where.  Placing a design point inside that unresolved interval is
    a coin flip that costs a full row budget when it loses, so the box is
    clipped to the nearest value KNOWN to complete -- not to just above the
    value known to fail.
    """
    lo, hi, kind = PRIOR_BOX[col]
    vlo, vhi = VALIDITY_BOX.get(col, (None, None))
    if vlo is not None:
        lo = max(lo, vlo)
    if vhi is not None:
        hi = min(hi, vhi)
    rec = (led.get("affordability") or {}).get(col)
    if not rec:
        return lo, hi, kind
    bad = rec.get("unaffordable_numeric") or []
    good = rec.get("affordable_numeric") or []
    below = [u for u in bad if u < base]
    above = [u for u in bad if u > base]
    if below:
        safe = [g for g in good if g > max(below)]
        lo = max(lo, min(safe) if safe else max(below))
    if above:
        safe = [g for g in good if g < min(above)]
        hi = min(hi, max(safe) if safe else min(above))
    return lo, hi, kind


def usable_span(col: str, led: dict, base: float, prefer: int) -> float:
    """Fraction of the prior box still available in the corrective direction."""
    lo0, hi0, kind = PRIOR_BOX[col]
    lo, hi, _ = box_for(col, led, base)
    full = (math.log10(hi0 / lo0) if kind == "log" and lo0 > 0 else hi0 - lo0)
    if prefer < 0:
        a, b = lo, min(base, hi)
    elif prefer > 0:
        a, b = max(base, lo), hi
    else:
        a, b = lo, hi
    if b <= a:
        return 0.0
    got = (math.log10(b / a) if kind == "log" and a > 0 else b - a)
    return max(0.0, got / full) if full else 0.0


def levels_for(col: str, base: float, n: int, prefer: int = 0,
               box=None) -> list:
    """n levels for `col`, spanning the prior box.

    `prefer` is the corrective direction: -1 sweeps from the box floor up to
    the base, +1 from the base up to the ceiling, 0 spans the whole box.  The
    base is always one of the endpoints when prefer != 0, so the stage stays
    anchored to the point already known to be grid-clean.
    """
    lo, hi, kind = box or PRIOR_BOX[col]
    base = min(max(base, lo), hi)
    if prefer < 0:
        lo, hi = lo, base
    elif prefer > 0:
        lo, hi = base, hi
    if hi <= lo:                     # base sits on the box edge -- fall back
        lo, hi, _ = box or PRIOR_BOX[col]
    if kind == "log" and lo > 0:
        a, b = math.log10(lo), math.log10(hi)
        vals = [10 ** (a + (b - a) * i / (n - 1)) for i in range(n)]
    else:
        vals = [lo + (hi - lo) * i / (n - 1) for i in range(n)]
    if col in INTEGER_LEVERS:
        # sorted(set(...)) can return FEWER levels than asked when the box is
        # narrow; that is correct -- duplicate rows would just be re-running
        # the same point -- and build_design sizes the factorial off the list.
        vals = sorted({int(round(v)) for v in vals})
    return vals


def prefer_for(col: str, obs: str, ratio) -> int:
    """Which way to sweep `col` to shrink the residual on `obs`."""
    sign = DIRECTION.get((col, obs))
    if not sign or not ratio:
        return 0
    # ratio > 1 means the model is too HIGH, so move the lever the way that
    # lowers the observable.
    return -sign if ratio > 1 else +sign


def parse_box_override(spec: str) -> dict:
    """--box i_mobile=5:32 -> {"i_mobile": (5.0, 32.0)}."""
    out = {}
    for part in (spec or "").split(","):
        part = part.strip()
        if not part:
            continue
        col, _, rng = part.partition("=")
        lo, _, hi = rng.partition(":")
        col = col.strip()
        if col not in PRIOR_BOX:
            raise SystemExit("--box %s: unknown lever (not in PRIOR_BOX)" % col)
        try:
            out[col] = (float(lo), float(hi))
        except ValueError:
            raise SystemExit("--box %s: expected lever=lo:hi" % part)
        if out[col][0] >= out[col][1]:
            raise SystemExit("--box %s: lo must be < hi" % part)
    return out


def build_design(led: dict, levers: list, n_levels: list, stage: str,
                 row_base: int, why: list = None, box_override: dict = None) -> tuple:
    best = led["best"]
    if not best:
        raise SystemExit("ledger has no valid best row -- nothing safe to pin to")
    params = dict(best["params"])
    target_obs = {c: o for c, o, _, _ in (why or [])}
    ratios = led.get("residuals") or {}

    grids = []
    for col, n in zip(levers, n_levels):
        try:
            base = float(params.get(col))
        except (TypeError, ValueError):
            base = sum(PRIOR_BOX[col][:2]) / 2.0
        obs = target_obs.get(col)
        lo, hi, kind = box_for(col, led, base)
        ov = (box_override or {}).get(col)
        if ov:
            # INTERSECT, never replace.  An override is a statement about where
            # this stage should spend its rows, not a licence to leave the
            # physical box or re-enter a span the ledger measured as
            # unaffordable.
            lo, hi = max(lo, ov[0]), min(hi, ov[1])
            if lo >= hi:
                raise SystemExit(
                    "--box %s=%g:%g leaves nothing after intersecting the "
                    "admissible span (%g, %g) for that lever."
                    % (col, ov[0], ov[1], *box_for(col, led, base)[:2]))
        grids.append(levels_for(col, base, n,
                                prefer_for(col, obs, ratios.get(obs)),
                                (lo, hi, kind)))

    combos = [[]]
    for g in grids:
        combos = [c + [v] for c in combos for v in g]

    header = ["row_id", "condition", "cond_row_id", "matrix", "base_idx",
              "param_j", "theta_id"] + [k for k in params if k not in BOOKKEEPING]
    rows, labels = [], {}
    for i, combo in enumerate(combos):
        rid = row_base + i
        rec = {"row_id": rid, "condition": best.get("condition") or "S330",
               "cond_row_id": rid, "matrix": "S", "base_idx": i,
               "param_j": -1, "theta_id": i}
        for k in header[7:]:
            rec[k] = params[k]
        tag = []
        for col, v in zip(levers, combo):
            if col in INTEGER_COLS:
                v = int(round(v))
            rec[col] = repr(v) if isinstance(v, float) else v
            tag.append("%s%.3g" % (col.replace("_", ""), v))
        rows.append(rec)
        labels[str(rid)] = "%s%02d_%s" % (stage, i, "_".join(tag))
    return header, rows, labels


def cost_key(cm: dict, machine: dict):
    """Which cost_model entry is this registry machine?

    THE TWO FILES KEY THE SAME MACHINE DIFFERENTLY.  learn.py keys cost_model by
    platform.node() -- the value rows actually carry -- while machines.json keys
    by human name.  'MATRIX-PC2' never equals 'Matrix-PC', so entry 1 matched
    nothing and fell through to the "scaled from another machine" branch, which
    then picked MATRIX-PC2's OWN numbers (it has the worst median) and inflated
    them by the 0.5/0.4 speed ratio: a measured 14.02 h worst row was reported
    as 17.5 h and tripped the budget refusal on a design that fits.  Match on
    the registry's own node_regex, which is the one field that already means
    "what this host calls itself".
    """
    for k in (machine["name"], machine["name"].replace(" ", "-")):
        if k in cm:
            return k
    pats = [r for r in [(machine.get("match") or {}).get("node_regex")] if r]
    pats += [a["node_regex"] for a in machine.get("match_alternatives", [])
             if a.get("node_regex")]
    for k in cm:
        if any(re.match(p, k, re.I) for p in pats):
            return k
    return None


def estimate_cost(led: dict, machine: dict, n_rows: int) -> dict:
    cm = led.get("cost_model", {})
    # Prefer this machine own measured rows; fall back to the campaign worst
    # case scaled by the registry speed ratio.
    key = cost_key(cm, machine)
    mine = cm.get(key) if key else None
    if mine:
        med, mx = mine["median_row_h"], mine["max_row_h"]
        src = "measured on this machine"
    elif any(v.get("n") for v in cm.values()):
        worst = max((v for v in cm.values() if v.get("n")),
                    key=lambda v: v["median_row_h"])
        ref_speed = 0.5
        scale = ref_speed / max(machine.get("speed", 1.0), 1e-6)
        med, mx = worst["median_row_h"] * scale, worst["max_row_h"] * scale
        src = "scaled from another machine by registry `speed`"
    else:
        return {"known": False}
    slots = machine.get("slots", 1)
    waves = math.ceil(n_rows / max(slots, 1))
    return {"known": True, "source": src, "median_row_h": round(med, 2),
            "max_row_h": round(mx, 2), "slots": slots, "waves": waves,
            "est_wall_h": round(mx * waves, 1)}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Propose the next calibration stage.")
    ap.add_argument("--machine", type=int, default=None)
    ap.add_argument("--stage", default=None)
    ap.add_argument("--levers", default=None,
                    help="comma-separated override, e.g. eta,f_cl_i,E_b_i2")
    ap.add_argument("--for", dest="for_obs", default=None,
                    help="observable an explicit --levers override is aimed at, "
                         "e.g. --for N_void; without it the sweep is symmetric "
                         "because the corrective direction is unknown")
    ap.add_argument("--levels", default=None,
                    help="comma-separated level count per lever, e.g. 3,2,2")
    ap.add_argument("--row-base", type=int, default=None)
    ap.add_argument("--box", default=None,
                    help="narrow a lever's box for this stage only, e.g. "
                         "--box gamma_s=1.5:1.8. Use to target a region a "
                         "sibling stage does not cover, instead of duplicating "
                         "it. Clamps into PRIOR_BOX; never widens past it.")
    ap.add_argument("--grid-v", type=int, default=2000,
                    help="max vacancy cluster size (--V). Raise for cavity "
                         "stages: at V=2000 every setting that grows realistic "
                         "cavities hits the ceiling and returns GRID+NOCONV.")
    ap.add_argument("--grid-vbin", type=int, default=20,
                    help="vacancy bin count. Scale with --grid-v to hold the "
                         "bin ratio: r = (V/v_discrete)^(1/v_bin), which is "
                         "1.349 at the V=2000/v_bin=20 default.")
    ap.add_argument("--budget-h", type=float, default=None,
                    help="row budget in hours; overrides machines.json")
    ap.add_argument("--share", default=None,
                    help="split this design across registry machine indices, "
                         "e.g. --share 0,2 for MacBook Pro + workstation")
    ap.add_argument("--write", action="store_true",
                    help="write the design files and update the ledger/guide")
    a = ap.parse_args(argv)

    led = load_ledger()
    if a.box:
        for item in a.box.split(","):
            col, rng = item.split("=")
            lo, hi = (float(x) for x in rng.split(":"))
            col = col.strip()
            if col not in PRIOR_BOX:
                raise SystemExit("--box names unknown lever %r" % col)
            plo, phi, kind = PRIOR_BOX[col]
            # Clamp INTO the prior box: a stage may narrow its search, never
            # widen it past a bound that exists for physical reasons.
            PRIOR_BOX[col] = (max(lo, plo), min(hi, phi), kind)
            print("box override: %s -> %.4g .. %.4g"
                  % (col, PRIOR_BOX[col][0], PRIOR_BOX[col][1]))
    machine, reg = detect(a.machine)
    worst, deferred = residual_rank(led)

    if a.levers:
        levers = [s.strip() for s in a.levers.split(",") if s.strip()]
        # An explicit override still may not sweep a parameter the physics does
        # not read: those rows come back bit-identical, so the stage is pure
        # waste no matter who asked for it.
        unwired = ((led.get("policy") or {}).get("unwired_parameters") or {})
        bad = [c for c in levers if c in unwired]
        if bad:
            for c in bad:
                print("REFUSING --levers %s: the physics never reads it." % c)
                print("  %s" % unwired[c])
            raise SystemExit(2)
        tgt = a.for_obs or "operator override"
        why = [(c, tgt, 0.0, -1) for c in levers]
        skipped = []
    else:
        levers, why, skipped = pick_levers(led, worst)
    if not levers:
        raise SystemExit("no admissible levers -- every candidate is dead in "
                         "the ledger. Widen LEVER_MAP or PRIOR_BOX.")

    slots = machine.get("slots", 1)
    if a.levels:
        n_levels = [int(x) for x in a.levels.split(",")]
    else:
        # Largest full factorial that fits the slot count in one wave.
        n_levels = [2] * len(levers)
        while math.prod(n_levels) * 3 // 2 <= slots and max(n_levels) < 4:
            n_levels[n_levels.index(min(n_levels))] += 1
    n_rows = math.prod(n_levels)

    # Stage name: next S-number after the highest already present.
    if a.stage:
        stage = a.stage
    else:
        nums = [int(s[1:].split("_")[0]) for s in led["stages"]
                if s.startswith("S") and s[1:].split("_")[0].isdigit()]
        stage = "S%d" % ((max(nums) if nums else 14) + 1)
    row_base = a.row_base if a.row_base is not None else 4000 + \
        (int(stage[1:]) - 10) * 100

    header, rows, labels = build_design(led, levers, n_levels, stage, row_base,
                                        why, parse_box_override(a.box))
    share = [int(x) for x in a.share.split(",")] if a.share else []
    cost = estimate_cost(led, machine,
                         math.ceil(n_rows / len(share)) if share else n_rows)

    if a.budget_h:
        budget_h, budget_src = a.budget_h, "--budget-h"
    elif machine.get("timeout_s"):
        budget_h, budget_src = machine["timeout_s"] / 3600.0, "machines.json (this machine)"
    elif reg.get("timeout_s"):
        budget_h, budget_src = reg["timeout_s"] / 3600.0, "machines.json (global default)"
    else:
        budget_h, budget_src = 20.0, "built-in fallback"
    over = (cost.get("known") and budget_h and cost["max_row_h"] > budget_h)

    b = led["best"]
    print("machine        : %d  %s  (%d slots)" % (machine["index"], machine["name"], slots))
    print("pinned to      : %s  (%s, %d/6 in range)"
          % (b["label"] or b["row_id"], b["stage"], b["n_in_range"]))
    print("worst residuals:")
    for k, m in worst[:4]:
        print("    %-7s %8.4gx off  (log10 = %.2f)" % (k, led["residuals"][k], m))
    for k, m in deferred:
        print("    %-7s %8.4gx off  (log10 = %.2f)  DEFERRED by policy"
              % (k, led["residuals"][k], m))
    print("levers chosen  :")
    tagname = {-1: "operator override", 0: "never varied",
               1: "inconclusive", 2: "live"}
    for col, obs, miss, pr in why:
        lv = next((n for c, n in zip(levers, n_levels) if c == col), "?")
        print("    %-12s %-2s levels   for %-7s  [%s]"
              % (col, lv, obs, tagname.get(pr, "?")))
    for col, obs, span in skipped:
        print("    %-12s SKIPPED for %-7s -- only %.0f%% of its box is "
              "affordable in the corrective direction" % (col, obs, span * 100))
    print("design         : %s, %d rows (%s)"
          % (stage, n_rows, " x ".join(str(n) for n in n_levels)))
    if cost.get("known"):
        print("est. cost      : %.1f h wall  (%d wave(s) x %.1f h worst row, %s)"
              % (cost["est_wall_h"], cost["waves"], cost["max_row_h"], cost["source"]))
    if share:
        print("sharded across : %s"
              % ", ".join("%d=%s" % (i, next(m["name"] for m in reg["machines"]
                                             if m["index"] == i)) for i in share))
    print("row budget     : %.1f h  (from %s)%s"
          % (budget_h, budget_src, "   *** EXCEEDED ***" if over else ""))
    print()

    if over:
        print("REFUSING to recommend as-is: the worst measured row (%.1f h) exceeds "
              "the row budget (%.1f h)." % (cost["max_row_h"], budget_h))
        print("Either pass --budget-h, raise timeout_s for machine %d in "
              "machines.json, or cut the design with --levels."
              % machine["index"])
        print()

    design_rel = "design/%s_calib.csv" % stage
    labels_rel = "design/%s_labels.json" % stage
    cmd = build_commands(stage, machine, reg, slots, budget_h, share,
                         a.grid_v, a.grid_vbin)

    if not a.write:
        print("DRY RUN -- nothing written. Re-run with --write to emit:")
        print("    %s   (%d rows)" % (design_rel, n_rows))
        print("    %s" % labels_rel)
        print()
        print("Then run:")
        print("\n".join(cmd))
        return 0

    dpath = HERE / design_rel
    with dpath.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    (HERE / labels_rel).write_text(json.dumps(labels, indent=1), encoding="utf-8")

    claim = {
        "stage": stage,
        "machine": machine["index"], "machine_name": machine["name"],
        "levers": levers, "n_levels": n_levels, "n_rows": n_rows,
        "design": design_rel, "labels": labels_rel,
        "pinned_to": {"stage": b["stage"], "row_id": b["row_id"],
                      "label": b["label"]},
        "rationale": ("worst residuals are %s; sweeping %s, which the ledger "
                      "has not retired" % (
                          ", ".join("%s (%.3gx)" % (k, led["residuals"][k])
                                    for k, _ in worst[:3]),
                          ", ".join(levers))),
        "cost_estimate": cost,
        "commands": cmd,
    }
    # PER-MACHINE CLAIMS.  A single `next` field is one slot for a fleet that
    # now runs three stages at once: writing S18 here silently erased
    # Matrix-PC's S17 claim.  `next` is kept as the claim made by THIS run so
    # the guide has something to render, but the durable record is per machine.
    led.setdefault("next_by_machine", {})[str(machine["index"])] = claim
    led["next"] = claim
    LEDGER.write_text(json.dumps(led, indent=2), encoding="utf-8")

    sys.path.insert(0, str(HERE))
    import learn
    (HERE / "CALIBRATION_GUIDE.md").write_text(learn.render_guide(led),
                                               encoding="utf-8")

    print("wrote %s (%d rows), %s" % (design_rel, n_rows, labels_rel))
    print("updated calibration_ledger.json[next] and CALIBRATION_GUIDE.md")
    print()
    print("Review the design, then run:")
    print("\n".join(cmd))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
