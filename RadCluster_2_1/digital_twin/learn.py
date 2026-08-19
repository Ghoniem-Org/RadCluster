#!/usr/bin/env python3
"""Ingest every calibration stage and re-derive what the campaign has learned.

WHY THIS EXISTS.  S1-S14 each answered one question, and the answer lived in a
commit message.  A machine that pulls the repo cannot read commit messages into
a plan, and a human returning after a crash re-derives the same verdicts by
hand (this happened on 2026-08-18).  This script turns the result files back
into two artefacts that DO travel:

    calibration_ledger.json   machine-readable state -- what was swept, what it
                              did, which levers are dead, where the best point
                              is, what each machine costs per row.
    CALIBRATION_GUIDE.md      the same thing for a human, regenerated so it can
                              never drift from the results it summarises.

`plan.py` reads the ledger and proposes the next stage.  Nothing here writes a
design or launches a run -- learning and planning are kept separate so that
re-ingesting results can never, by itself, start a job.

THE ONE RULE THIS FILE ENFORCES.  A row is a MEASUREMENT only if it reached the
target dose, is not grid-limited, and its grid converged.  Rows that piled at
the ceiling still carry information about DIRECTION, and are used for that, but
their size readouts are never scored or fitted.  S12 and the V2 block of S14
are the reason: both look like "big loops" and neither is a number.

Usage:
    python learn.py                 # re-ingest, rewrite ledger + guide
    python learn.py --quiet         # no stdout summary
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
LEDGER = HERE / "calibration_ledger.json"
GUIDE = HERE / "CALIBRATION_GUIDE.md"
TARGETS = HERE / "targets_330C_15dpa.json"

# Bookkeeping columns -- never a physics lever.
BOOKKEEPING = {"row_id", "condition", "cond_row_id", "matrix", "base_idx",
               "param_j", "theta_id"}

OBS = ["N_loops_100", "d_100_nm", "N_loops_111", "d_111_nm",
       "N_voids", "d_cavity_nm"]
SHORT = {"N_loops_100": "N_100", "d_100_nm": "d_100", "N_loops_111": "N_111",
         "d_111_nm": "d_111", "N_voids": "N_void", "d_cavity_nm": "d_void"}

# A lever whose FULL tested span moves every observable by less than this is
# declared dead.  5 % is well above the solver run-to-run reproducibility and
# well below anything that could close a 12x miss.
DEAD_THRESHOLD = 0.05

# Observables the planner should NOT chase, with the reason.  This is a
# deliberate, reviewable override of "worst residual first" -- it is seeded
# once and then preserved across regenerations like `notes`, so changing your
# mind means editing the ledger, not the code.
DEFAULT_POLICY = {
    "deprioritized_observables": {
        "N_void": ("Missed by ~300x AND nearly unresponsive to its own "
                   "governing parameters: the S14 vacancy triple moved "
                   "N_void only 2.57e18 -> 9.08e18 (3.5x) while moving loop "
                   "content 160x. A residual that large with a response that "
                   "small is evidence of a structural defect in the cavity "
                   "channel, not a parameter that needs tuning. Re-enable "
                   "once cavity nucleation is shown to respond at all."),
        "d_void": ("Pinned at 0.56-0.57 nm across every row of every stage, "
                   "including a 160x swing in loop content. Same reasoning as "
                   "N_void."),
    },
}


# --------------------------------------------------------------- loading ----
def load_targets() -> dict:
    return json.loads(TARGETS.read_text(encoding="utf-8"))


def discover_stages() -> dict:
    """Pair every results/*.jsonl with its design/*.csv and labels, if present."""
    stages = {}
    for res in sorted((HERE / "results").glob("*.jsonl")):
        stem = res.stem
        if stem.startswith("archive"):
            continue
        # S14_calib.jsonl -> S14 ; T9_screen_partial.jsonl -> T9
        tag = stem.split("_")[0]
        design = None
        for cand in (HERE / "design" / (stem + ".csv"),
                     HERE / "design" / (tag + "_calib.csv"),
                     HERE / "design" / (tag + "_probe.csv")):
            if cand.exists():
                design = cand
                break
        labels = {}
        for cand in (HERE / "design" / (tag + "_labels.json"),
                     HERE / "design" / (stem + "_labels.json")):
            if cand.exists():
                labels = json.loads(cand.read_text(encoding="utf-8"))
                break
        stages[stem] = {"tag": tag, "results": res, "design": design,
                        "labels": labels}
    return stages


def read_rows(path: Path) -> list:
    out = []
    for ln in path.read_text(encoding="utf-8").splitlines():
        ln = ln.strip()
        if not ln:
            continue
        try:
            out.append(json.loads(ln))
        except json.JSONDecodeError:
            pass
    return out


def read_design(path) -> dict:
    if not path or not Path(path).exists():
        return {}
    with Path(path).open(encoding="utf-8") as fh:
        return {r["row_id"]: r for r in csv.DictReader(fh)}


# -------------------------------------------------------------- validity ----
def classify(row: dict, dose_target: float = 15.0):
    """Is this row a measurement?  Returns (valid, flags)."""
    flags = []
    if row.get("grid_limited"):
        flags.append("GRID")
    if not row.get("grid_converged", True):
        flags.append("NOCONV")
    reached = row.get("dose_reached", 0.0) or 0.0
    if abs(reached - dose_target) > 0.02 * dose_target:
        flags.append("DOSE")
    if row.get("starved"):
        flags.append("STARV")
    if not row.get("admissible", True):
        flags.append("INADM")
    if not row.get("conserving", True):
        flags.append("NOCONS")
    # NOCONS is recorded but does NOT disqualify: every <100> row in the
    # campaign carries it (see commit ddcb203).  Disqualifying on it would
    # leave zero measurements and hide the rest of the signal.
    hard = {"GRID", "NOCONV", "DOSE", "STARV", "INADM"}
    return (not (hard & set(flags))), flags


def score(row: dict, targets: dict) -> dict:
    o = targets["observables"]
    per, n_ok = {}, 0
    for k in OBS:
        v = row.get(k)
        if v is None:
            per[k] = None
            continue
        ok = o[k]["lo"] <= v <= o[k]["hi"]
        n_ok += ok
        per[k] = {"value": v, "in_range": bool(ok),
                  "ratio": v / o[k]["target"] if o[k]["target"] else None}
    return {"n_in_range": n_ok, "per": per}


def log_distance(row: dict, targets: dict) -> float:
    """RMS of log10(model/target) over the six.  Scale-free, so a 300x miss on
    voids is not drowned out by a 1.2x miss on a diameter."""
    o = targets["observables"]
    acc, n = 0.0, 0
    for k in OBS:
        v, t = row.get(k), o[k]["target"]
        if v and t and v > 0:
            acc += math.log10(v / t) ** 2
            n += 1
    return math.sqrt(acc / n) if n else float("inf")


# ---------------------------------------------------- lever sensitivity -----
def covarying_groups(design: dict) -> list:
    """Partition swept columns into blocks that always move together.

    V0/V1/V2 changes f_cl_v, E_b_v2 and s_v as one knob; the S12 sink block
    moves five Z-factors at once.  Treating those as independent columns would
    find no single-column pairs at all and report every lever as untested.
    """
    rows = list(design.values())
    if not rows:
        return []
    swept = [c for c in rows[0]
             if c not in BOOKKEEPING and len({r[c] for r in rows}) > 1]
    sig = {}
    for c in swept:
        seen, pattern = {}, []
        for r in rows:                      # value -> class index, by first use
            pattern.append(seen.setdefault(r[c], len(seen)))
        sig.setdefault(tuple(pattern), []).append(c)
    return list(sig.values())


def group_key(design_row: dict, group: list) -> tuple:
    return tuple(design_row[c] for c in group)


def sensitivity(stage: dict, rows: list, targets: dict,
                dose_target: float = 15.0) -> dict:
    """Response of each observable to each co-varying lever group.

    Only pairs that differ in EXACTLY ONE group are used, so the response is
    attributable.  Magnitude is measured on valid rows only; a pair where one
    side went invalid is recorded as a direction with `drove_invalid`.
    """
    design = read_design(stage["design"])
    if not design:
        return {}
    groups = covarying_groups(design)
    by_id = {str(r["row_id"]): r for r in rows if r.get("row_id") is not None}
    out = {}

    for group in groups:
        others = [g for g in groups if g is not group]
        pairs, drove_invalid = [], False
        ids = [i for i in by_id if i in design]
        for a in ids:
            for b in ids:
                if a >= b:
                    continue
                da, db = design[a], design[b]
                if group_key(da, group) == group_key(db, group):
                    continue
                if any(group_key(da, g) != group_key(db, g) for g in others):
                    continue
                ra, rb = by_id[a], by_id[b]
                va, _ = classify(ra, dose_target)
                vb, _ = classify(rb, dose_target)
                if va and vb:
                    pairs.append((da, db, ra, rb))
                else:
                    drove_invalid = True

        resp = {}
        for k in OBS:
            best = 0.0
            for da, db, ra, rb in pairs:
                x, y = ra.get(k), rb.get(k)
                if x and y and x > 0 and y > 0:
                    best = max(best, abs(y / x - 1.0))
            resp[k] = best

        span = {}
        for c in group:
            vals = sorted({float(design[i][c]) for i in ids})
            span[c] = [vals[0], vals[-1]]

        peak = max(resp.values()) if resp else 0.0
        out["+".join(group)] = {
            "columns": group, "span": span, "n_pairs": len(pairs),
            "response": resp, "peak_response": peak,
            "drove_invalid": drove_invalid,
            "verdict": ("dead" if (pairs and peak < DEAD_THRESHOLD)
                        else "live" if pairs
                        else "inconclusive"),
        }
    return out


# ---------------------------------------------------------- affordability ---
def affordability(stage: dict, rows: list, dose_target: float) -> dict:
    """Which lever VALUES prevent a row from ever reaching dose.

    This is a COST verdict, not a physics one, and it is measured with a
    controlled contrast: within one stage, a level counts as unaffordable only
    if it produced zero full-dose rows while ANOTHER level of the same
    co-varying group -- everything else held equal -- produced at least one.

    S14 is the motivating case.  All six rows at E_b_i2 = 0.60 burned the full
    20 h budget and reached 0.004-6.06 of 15 dpa; all six at 0.75 finished in
    5.3-19.8 h.  Nothing about that is visible in the physics verdicts, because
    a row that never reaches dose is excluded from scoring -- so without this
    the planner cheerfully proposes another stage in the same dead direction.
    """
    design = read_design(stage["design"])
    if not design:
        return {}
    by_id = {str(r["row_id"]): r for r in rows if r.get("row_id") is not None}
    out = {}
    for group in covarying_groups(design):
        tally = {}
        for rid, drow in design.items():
            r = by_id.get(rid)
            if r is None:
                continue
            reached = abs((r.get("dose_reached") or 0.0) - dose_target) \
                <= 0.02 * dose_target
            key = group_key(drow, group)
            slot = tally.setdefault(key, [0, 0])
            slot[0] += 1
            slot[1] += reached
        if len(tally) < 2:
            continue
        any_ok = any(v[1] > 0 for v in tally.values())
        # ATTRIBUTE ONLY WHAT IS ATTRIBUTABLE.  A co-varying group that fails
        # says nothing about which member caused it: S14's V2 level moved
        # f_cl_v, E_b_v2 and s_v together, so blaming each individually would
        # cap E_b_v2 at 0.35 on evidence that never isolated it.  Multi-column
        # groups are recorded under the joined name and never clip a single
        # column's box.
        if len(group) > 1:
            continue
        for key, (n, ok) in tally.items():
            for col, val in zip(group, key):
                rec = out.setdefault(col, {"unaffordable": [], "affordable": [],
                                           "evidence": []})
                if ok == 0 and n >= 3 and any_ok:
                    if val not in rec["unaffordable"]:
                        rec["unaffordable"].append(val)
                        rec["evidence"].append(
                            "%s: 0 of %d rows reached dose at %s=%s"
                            % (stage.get("stem", "?"), n, col, val))
                elif ok > 0 and val not in rec["affordable"]:
                    rec["affordable"].append(val)
    # Drop columns with nothing adverse to say.
    return {c: r for c, r in out.items() if r["unaffordable"]}


# ------------------------------------------------------------ inventory -----
def inventory(row: dict, targets: dict):
    """Compare SIA content locked in loops against the measurement.

    Uses the ROW OWN n<->d relation (n / d^2 from its own mean_n and mean d),
    so no lattice constant or Burgers vector is assumed here and the comparison
    cannot drift from whatever the solver actually used.
    """
    o = targets["observables"]
    out = {}
    for char, nkey, dkey, Nkey, tN, td in (
        ("100", "mean_n_100", "d_100_nm", "N_loops_100", "N_loops_100", "d_100_nm"),
        ("111", "mean_n_111", "d_111_nm", "N_loops_111", "N_loops_111", "d_111_nm"),
    ):
        n, d, N = row.get(nkey), row.get(dkey), row.get(Nkey)
        if not (n and d and N and d > 0):
            return None
        n_per_d2 = n / (d ** 2)
        n_exp = n_per_d2 * (o[td]["target"] ** 2)
        out[char] = {
            "model_N": N, "model_mean_n": n, "model_content": N * n,
            "exp_N": o[tN]["target"], "exp_mean_n_at_target_d": n_exp,
            "exp_content": o[tN]["target"] * n_exp,
            "ratio": (N * n) / (o[tN]["target"] * n_exp),
        }
    tot_m = sum(v["model_content"] for v in out.values())
    tot_e = sum(v["exp_content"] for v in out.values())
    out["total"] = {"model_content": tot_m, "exp_content": tot_e,
                    "ratio": tot_m / tot_e if tot_e else None}
    return out


# ----------------------------------------------------------- cost model -----
def cost_model(all_rows: list) -> dict:
    """Wall time per row, per machine -- from COMPLETED rows only.

    A row cut at the budget records wall_s == the budget, which measures the
    timeout and not the cost.  Counting those pins max_row_h to the budget, and
    since plan.py refuses a design whose worst row exceeds the budget, one
    timed-out stage would make the planner refuse every subsequent design on
    that machine.  S14 did exactly this: six rows at 20.0 h that never reached
    dose.  Timeouts are counted separately, because a machine that times out
    often is still worth knowing about.
    """
    done, timed_out = {}, {}
    for r in all_rows:
        w, m = r.get("wall_s"), r.get("machine_id")
        if not (w and m):
            continue
        tgt = r.get("dose_target")
        reached = r.get("dose_reached")
        finished = (tgt is None or reached is None
                    or abs(reached - tgt) <= 0.02 * tgt)
        (done if finished else timed_out).setdefault(m, []).append(float(w))
    out = {}
    # SORTED, or the key order differs per machine and every learn.py run
    # produces a spurious ledger diff -- the exact churn 72e1c9f fixed
    # elsewhere.  Unordered iteration here reintroduced it through cost_model.
    for m in sorted(set(done) | set(timed_out)):
        ws = sorted(done.get(m, []))
        rec = {"n": len(ws), "n_timed_out": len(timed_out.get(m, []))}
        if ws:
            rec.update({"median_row_h": round(ws[len(ws) // 2] / 3600.0, 2),
                        "max_row_h": round(ws[-1] / 3600.0, 2),
                        "min_row_h": round(ws[0] / 3600.0, 2)})
        out[m] = rec
    return out


# --------------------------------------------------------------- ingest -----
def ingest() -> dict:
    targets = load_targets()
    stages_meta = discover_stages()
    prev = json.loads(LEDGER.read_text(encoding="utf-8")) if LEDGER.exists() else {}

    goal_dose = float(targets["condition"]["dose_dpa"])
    stages, all_rows, best = {}, [], None
    for stem, meta in sorted(stages_meta.items()):
        rows = read_rows(meta["results"])
        if not rows:
            continue
        all_rows.extend(rows)

        # A stage is IN SCOPE only if it ran to the dose the targets describe.
        # The T-series ran to 1 dpa; scoring those against a 15 dpa target set
        # would flag every row DOSE and then, because no valid pair survives,
        # silently report their levers as "dead".  That is how f_cl_i -- never
        # validly tested at all -- was reported dead on the first ingest.
        doses = [r.get("dose_target") for r in rows if r.get("dose_target")]
        stage_dose = max(set(doses), key=doses.count) if doses else goal_dose
        in_scope = abs(stage_dose - goal_dose) < 1e-9

        design = read_design(meta["design"])
        recs, n_valid = [], 0
        for r in rows:
            valid, flags = classify(r, stage_dose)
            n_valid += valid
            sc = score(r, targets)
            drow = design.get(str(r.get("row_id")), {})
            rec = {"row_id": r.get("row_id"),
                   "label": meta["labels"].get(str(r.get("row_id"))),
                   "valid": valid, "flags": flags,
                   "n_in_range": sc["n_in_range"],
                   "log_distance": round(log_distance(r, targets), 3),
                   "dose_reached": r.get("dose_reached"),
                   "observables": {SHORT[k]: (sc["per"][k]["value"]
                                              if sc["per"][k] else None)
                                   for k in OBS}}
            recs.append(rec)
            if valid and in_scope:
                cand = {"stage": stem, "row_id": r.get("row_id"),
                        "label": rec["label"], "n_in_range": sc["n_in_range"],
                        "log_distance": rec["log_distance"],
                        "flags": flags,
                        "observables": rec["observables"],
                        "ratios": {SHORT[k]: (round(sc["per"][k]["ratio"], 4)
                                              if sc["per"][k] and sc["per"][k]["ratio"]
                                              else None) for k in OBS},
                        "params": {c: v for c, v in drow.items()
                                   if c not in BOOKKEEPING},
                        "inventory": inventory(r, targets)}
                if best is None or (cand["n_in_range"], -cand["log_distance"]) > \
                                   (best["n_in_range"], -best["log_distance"]):
                    best = cand

        swept = ["+".join(g) for g in covarying_groups(design)] if design else []
        stages[stem] = {
            # as_posix(), never str(): str() of a Path emits os.sep, so the same
            # results rewrite every stored path when the ingesting machine
            # changes OS.  Forward slashes read back correctly on Windows too.
            "design": (meta["design"].relative_to(HERE).as_posix()
                       if meta["design"] else None),
            "results": meta["results"].relative_to(HERE).as_posix(),
            "dose_target": stage_dose, "in_scope": in_scope,
            "n_rows": len(rows), "n_valid": n_valid,
            "swept": swept,
            "sensitivity": sensitivity(meta, rows, targets, stage_dose),
            "affordability": (affordability(dict(meta, stem=stem), rows, stage_dose)
                              if in_scope else {}),
            "rows": recs,
        }

    # Roll per-stage sensitivity up into one verdict per lever group.  Only
    # in-scope stages vote: a 1 dpa response is not evidence about a 15 dpa
    # target, and mixing them would let a stage that never reached dose
    # manufacture a verdict.
    levers = {}
    for stem, st in stages.items():
        if not st["in_scope"]:
            continue
        for name, s in st["sensitivity"].items():
            L = levers.setdefault(name, {"columns": s["columns"], "seen_in": [],
                                         "peak_response": 0.0, "response": {},
                                         "spans": {}, "n_pairs": 0,
                                         "drove_invalid": False})
            L["seen_in"].append(stem)
            L["peak_response"] = max(L["peak_response"], s["peak_response"])
            L["n_pairs"] += s["n_pairs"]
            L["drove_invalid"] = L["drove_invalid"] or s["drove_invalid"]
            for k, v in s["response"].items():
                L["response"][k] = max(L["response"].get(k, 0.0), v)
            for c, sp in s["span"].items():
                old = L["spans"].get(c)
                L["spans"][c] = ([min(sp[0], old[0]), max(sp[1], old[1])]
                                 if old else sp)
    for name, L in levers.items():
        # ZERO valid pairs means UNTESTED, never dead.  A dead verdict retires
        # a lever from every future design, so it must be earned by a pair of
        # real measurements that bracket the span and did not move.
        if L["n_pairs"] == 0:
            L["verdict"] = "inconclusive"
        elif L["peak_response"] < DEAD_THRESHOLD:
            L["verdict"] = "dead"
        else:
            L["verdict"] = "live"
        moving = [k for k, v in L["response"].items() if v >= DEAD_THRESHOLD]
        moving.sort(key=lambda k: -L["response"][k])
        L["moves"] = [SHORT[k] for k in moving]

    # REACHABILITY.  How big a change each observable needs, against the
    # biggest change any single lever has actually produced across its full
    # tested span.  S20 is why this exists: it swept E_m_i, L_hat and B_111
    # over their whole boxes and moved d_111 from 1.05 to at most 1.33 nm,
    # when reaching the measured 6.2 nm needs +490 %.  Without this the
    # planner keeps proposing those levers for d_111 forever, because "live"
    # only means "moves more than 5 %".
    #
    # NOT A PROOF OF IMPOSSIBILITY.  Levers can compound, so a residual no
    # SINGLE lever can close may still be reachable by several together.  It
    # is recorded as a strong hint and ranked on, never as a hard block.
    reach = {}
    if best:
        for k in OBS:
            s_ = SHORT[k]
            ratio = (best["ratios"] or {}).get(s_)
            if not ratio:
                continue
            need = (1.0 / ratio) if ratio < 1 else ratio
            top, who = 0.0, None
            for name, L in levers.items():
                r = L["response"].get(k, 0.0)
                if r > top:
                    top, who = r, name
            reach[s_] = {
                "ratio": ratio,
                "needed_factor": round(need, 3),
                "needed_pct": round((need - 1) * 100, 1),
                "best_single_lever_pct": round(top * 100, 1),
                "best_single_lever": who,
                "single_lever_sufficient": bool(top >= (need - 1.0)),
            }

    # Which physics columns has the campaign never varied AT THE TARGET DOSE?
    # "Never varied" is literal -- the column took one value in every in-scope
    # design.  A column that WAS varied but yielded no valid pair is not listed
    # here; it carries an `inconclusive` verdict in the lever table instead,
    # which is a different and weaker kind of ignorance.
    tested_cols = {c for L in levers.values() for c in L["columns"]}
    all_cols = set()
    for stem, st in stages.items():
        if not st["in_scope"]:
            continue
        d = read_design(HERE / st["design"]) if st["design"] else {}
        if d:
            all_cols |= {c for c in next(iter(d.values())) if c not in BOOKKEEPING}
    untested = sorted(all_cols - tested_cols)

    # Roll affordability up across in-scope stages.
    afford = {}
    for stem, st in stages.items():
        for col, rec in (st.get("affordability") or {}).items():
            a = afford.setdefault(col, {"unaffordable": [], "affordable": [],
                                        "evidence": []})
            for kind in ("unaffordable", "affordable"):
                for v in rec.get(kind, []):
                    if v not in a[kind]:
                        a[kind].append(v)
            a["evidence"].extend(rec["evidence"])
    for col, a in afford.items():
        for kind in ("unaffordable", "affordable"):
            try:
                a[kind + "_numeric"] = sorted(float(v) for v in a[kind])
            except (TypeError, ValueError):
                pass

    ledger = {
        "schema": 1,
        "updated": None,                     # set below, from CONTENT not clock
        "goal": {
            "statement": ("Find one parameter vector that puts all six "
                          "observables inside their experimental ranges at "
                          "EUROFER97 330 C / 15 dpa, neutron."),
            "targets_file": TARGETS.name,
            "condition": targets["condition"],
        },
        "stages": stages,
        "levers": levers,
        "untested_columns": untested,
        "reachability": reach,
        "affordability": afford,
        "best": best,
        "residuals": best["ratios"] if best else {},
        "cost_model": cost_model(all_rows),
        "policy": prev.get("policy", DEFAULT_POLICY),
        "notes": prev.get("notes", []),      # hand-written insight survives
        "next": prev.get("next"),            # plan.py owns these two fields
        "next_by_machine": prev.get("next_by_machine", {}),
    }

    # RE-INGESTING THE SAME RESULTS MUST BE A BYTE-FOR-BYTE NO-OP.  The ledger
    # travels between machines through git, and step 1 of the protocol is
    # `git pull --rebase`.  Any field that moves on every run -- a wall-clock
    # stamp, or an os.sep-flavoured path -- makes two machines rewrite the same
    # lines of a 30k-line JSON file to different values, which is a rebase
    # conflict; and a machine that cannot rebase cannot ingest anyone else's
    # results, which is exactly the cross-machine learning this file exists to
    # carry.  So `updated` records when the DERIVED CONTENT last changed, not
    # when the script last ran.  Compare as serialized JSON: `prev` came back
    # through json.loads, so a tuple built here would never equal its list.
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")
    if prev:
        def body(d):
            return json.dumps({k: v for k, v in d.items() if k != "updated"},
                              sort_keys=True, indent=2)
        if body(ledger) == body(prev):
            stamp = prev.get("updated") or stamp
    ledger["updated"] = stamp
    return ledger


# --------------------------------------------------------------- render -----
def fmt(v, sig=3):
    if v is None:
        return "-"
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, (int, float)):
        return "%.*g" % (sig, v)
    return str(v)


def render_guide(led: dict) -> str:
    T = load_targets()["observables"]
    b = led["best"]
    L = []
    A = L.append

    A("# Calibration Guide - EUROFER97 digital twin")
    A("")
    A("*Derived by `learn.py`; content last changed %s (re-running with no new "
      "results leaves this file untouched). Do not hand-edit: edits are "
      "overwritten. Durable notes belong in `calibration_ledger.json` under "
      "`notes`, which is preserved across regenerations.*" % led["updated"])
    A("")
    A("## Goal")
    A("")
    A(led["goal"]["statement"])
    A("")
    c = led["goal"]["condition"]
    A("Condition: **%s, %s C, %s dpa, %s** - targets in `%s`."
      % (c["material"], c["T_C"], c["dose_dpa"], c["irradiation"],
         led["goal"]["targets_file"]))
    A("")

    A("## Where the campaign stands")
    A("")
    if not b:
        A("No valid (full-dose, grid-clean, converged) row yet.")
    else:
        A("Best valid row: **%s** from `%s` - **%d/6** observables in range, "
          "log-distance %s."
          % (b["label"] or b["row_id"], b["stage"], b["n_in_range"],
             b["log_distance"]))
        A("")
        A("| observable | model | target | range | ratio | in range |")
        A("|---|---|---|---|---|---|")
        for k in OBS:
            s = SHORT[k]
            v = b["observables"].get(s)
            t = T[k]
            inr = (t["lo"] <= v <= t["hi"]) if v is not None else False
            A("| %s | %s | %s | %s - %s | %s x | %s |"
              % (s, fmt(v), fmt(t["target"]), fmt(t["lo"]), fmt(t["hi"]),
                 fmt(b["ratios"].get(s)), "yes" if inr else "**no**"))
        A("")
        inv = b.get("inventory")
        if inv:
            A("### Defect inventory (SIA content locked in loops)")
            A("")
            A("Computed from the row own n<->d relation, so no lattice "
              "constant is assumed.")
            A("")
            A("| character | model N x nbar | experiment N x nbar | ratio |")
            A("|---|---|---|---|")
            for ch in ("100", "111"):
                d = inv[ch]
                A("| <%s> | %s | %s | %s x |"
                  % (ch, fmt(d["model_content"]), fmt(d["exp_content"]),
                     fmt(d["ratio"])))
            A("| **total** | %s | %s | **%s x** |"
              % (fmt(inv["total"]["model_content"]),
                 fmt(inv["total"]["exp_content"]),
                 fmt(inv["total"]["ratio"])))
            A("")

    rch = led.get("reachability") or {}
    if rch:
        A("## Can any lever still close this?")
        A("")
        A("How far each observable must move from the best row, against the "
          "largest change any SINGLE lever has actually produced across its "
          "full tested span. Levers can compound, so a `no` is a strong hint "
          "that the residual is structural rather than a proof that it is.")
        A("")
        A("| observable | ratio | needs | best single lever | that lever | single lever enough? |")
        A("|---|---|---|---|---|---|")
        for s_ in ("N_100", "d_100", "N_111", "d_111", "N_void", "d_void"):
            r = rch.get(s_)
            if not r:
                continue
            A("| %s | %s x | %+.0f%% | %+.0f%% | `%s` | %s |"
              % (s_, fmt(r["ratio"]), r["needed_pct"],
                 r["best_single_lever_pct"], r["best_single_lever"] or "-",
                 "yes" if r["single_lever_sufficient"] else "**no**"))
        A("")

    A("## What each lever does")
    A("")
    A("A lever is **dead** when its full tested span moves every observable by "
      "less than %d%%. Only pairs of rows differing in exactly one lever are "
      "used, and only rows that are real measurements (full dose, grid-clean, "
      "converged)." % int(DEAD_THRESHOLD * 100))
    A("")
    A("`inconclusive` means the span was attempted but no pair of valid rows "
      "bracketed it - usually because the rows piled at the grid ceiling or "
      "never reached dose. An inconclusive lever is an OPEN question, not a "
      "closed one.")
    A("")
    A("| lever | tested span | verdict | valid pairs | peak response | moves | stages |")
    A("|---|---|---|---|---|---|---|")
    for name, v in sorted(led["levers"].items(),
                          key=lambda kv: (kv[1]["verdict"] != "live",
                                          -kv[1]["peak_response"])):
        span = "; ".join("%s %s->%s" % (c, fmt(s[0]), fmt(s[1]))
                         for c, s in v["spans"].items())
        moves = ", ".join(v["moves"]) or "-"
        note = " (drives rows off-grid)" if v["drove_invalid"] else ""
        peak = ("%.1f%%" % (v["peak_response"] * 100)) if v["n_pairs"] else "-"
        A("| `%s` | %s | **%s**%s | %d | %s | %s | %s |"
          % (name, span, v["verdict"], note, v["n_pairs"], peak, moves,
             ", ".join(sorted(set(v["seen_in"])))))
    A("")
    if led["untested_columns"]:
        A("### Never varied")
        A("")
        A("These columns exist in the design but have never taken more than one "
          "value, so the campaign has no evidence about them:")
        A("")
        A("`" + "`, `".join(led["untested_columns"]) + "`")
        A("")

    A("## Stage history")
    A("")
    A("Only stages that ran to %s dpa are in scope; the rest are listed for "
      "provenance but cast no vote on any lever." % c["dose_dpa"])
    A("")
    A("| stage | dose | scope | rows | valid | swept | best in-range |")
    A("|---|---|---|---|---|---|---|")
    for stem, st in sorted(led["stages"].items()):
        bi = max([r["n_in_range"] for r in st["rows"] if r["valid"]] or [None])
        A("| `%s` | %s | %s | %d | %d | %s | %s |"
          % (stem, fmt(st["dose_target"]),
             "in" if st["in_scope"] else "out",
             st["n_rows"], st["n_valid"],
             ", ".join("`%s`" % s for s in st["swept"]) or "-",
             bi if bi is not None else "-"))
    A("")

    A("## Cost model (measured)")
    A("")
    A("Completed rows only. A row cut at the budget measures the timeout, not "
      "the cost, so timeouts are counted in their own column.")
    A("")
    A("| machine | completed | median row | min | max | timed out |")
    A("|---|---|---|---|---|---|")
    for m, cm in sorted(led["cost_model"].items()):
        if cm["n"]:
            A("| %s | %d | %s h | %s h | %s h | %d |"
              % (m, cm["n"], cm["median_row_h"], cm["min_row_h"],
                 cm["max_row_h"], cm["n_timed_out"]))
        else:
            A("| %s | 0 | - | - | - | %d |" % (m, cm["n_timed_out"]))
    A("")
    A("`plan.py` sizes a stage from this table and the machine slot count, and "
      "refuses to propose a design whose estimated cost exceeds the machine "
      "row budget.")
    A("")

    aff = led.get("affordability") or {}
    if aff:
        A("## Unaffordable lever values")
        A("")
        A("Measured, not assumed: a level counts here only if it produced ZERO "
          "full-dose rows while another level of the same lever - everything "
          "else held equal - produced at least one. `plan.py` will not place a "
          "design point on these values.")
        A("")
        A("| lever | unaffordable at | evidence |")
        A("|---|---|---|")
        for col, rec in sorted(aff.items()):
            A("| `%s` | %s | %s |"
              % (col, ", ".join(str(v) for v in rec["unaffordable"]),
                 "; ".join(rec["evidence"])))
        A("")

    dep = (led.get("policy") or {}).get("deprioritized_observables") or {}
    if dep:
        A("## Deferred observables")
        A("")
        A("The planner will not propose levers for these until the entry is "
          "removed from `calibration_ledger.json` under "
          "`policy.deprioritized_observables`. They still appear in every "
          "score above - they are deferred, not ignored.")
        A("")
        for k, reason in sorted(dep.items()):
            A("- **%s** - %s" % (k, reason))
        A("")

    if led.get("notes"):
        A("## Curated notes")
        A("")
        for n in led["notes"]:
            A("- %s" % n)
        A("")

    claims = led.get("next_by_machine") or {}
    if claims:
        A("## Claimed stages")
        A("")
        A("One row per machine. A machine claims a stage by running `plan.py "
          "--write` on it; the claim is not a lock, only a record of what that "
          "machine was last told to run.")
        A("")
        A("| machine | stage | levers | rows | design |")
        A("|---|---|---|---|---|")
        for idx in sorted(claims, key=lambda k: int(k)):
            c = claims[idx]
            A("| %s (%s) | **%s** | %s | %s | `%s` |"
              % (idx, c.get("machine_name", "?"), c.get("stage"),
                 ", ".join("`%s`" % l for l in c.get("levers", [])),
                 c.get("n_rows", "?"), c.get("design", "")))
        A("")

    A("## Next stage")
    A("")
    nx = led.get("next")
    if not nx:
        A("None proposed yet - run `python plan.py` on the machine that will "
          "execute it.")
    else:
        A("**%s** - %s" % (nx["stage"], nx["rationale"]))
        A("")
        A("Sweeping: %s" % ", ".join("`%s`" % s for s in nx["levers"]))
        A("")
        if nx.get("design"):
            A("Design: `%s` (%s rows)" % (nx["design"], nx.get("n_rows", "?")))
            A("")
        if nx.get("commands"):
            A("Run it with:")
            A("")
            A("```bash")
            for line in nx["commands"]:
                A(line)
            A("```")
            A("")

    A("## Multi-machine protocol")
    A("")
    A("Results travel through git; `campaign_ops.sync_results()` is safe to "
      "call while a run is in flight (rows are appended, so the newest lands "
      "in the next sync).")
    A("")
    A("```bash")
    A("git pull --rebase origin main            # 1. take everyone results")
    A("python learn.py                          # 2. re-derive ledger + this guide")
    A("python plan.py                           # 3. propose this machine next stage")
    A("#    ... review the printed command, then run it ...")
    A("python -c \"import campaign_ops as c; c.sync_results()\"   # 4. publish rows")
    A("```")
    A("")
    A("Two files must never be committed - both are in `.gitignore` and both "
      "have caused an outage before: `CAMPAIGN_STOP` (every machine that "
      "pulled it halted at startup, commit 8a600c7) and `results/*.pid`.")
    A("")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description="Re-ingest calibration results.")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)

    led = ingest()
    LEDGER.write_text(json.dumps(led, indent=2), encoding="utf-8")
    GUIDE.write_text(render_guide(led), encoding="utf-8")

    if not a.quiet:
        b = led["best"]
        print("ingested %d stage files, %d rows (%d valid)"
              % (len(led["stages"]),
                 sum(s["n_rows"] for s in led["stages"].values()),
                 sum(s["n_valid"] for s in led["stages"].values())))
        if b:
            print("best: %s (%s) %d/6 in range, log-distance %s"
                  % (b["label"] or b["row_id"], b["stage"], b["n_in_range"],
                     b["log_distance"]))
        by = {}
        for k, v in led["levers"].items():
            by.setdefault(v["verdict"], []).append(k)
        print("in-scope stages: %d of %d"
              % (sum(1 for s in led["stages"].values() if s["in_scope"]),
                 len(led["stages"])))
        print("levers: %d live, %d dead, %d inconclusive"
              % (len(by.get("live", [])), len(by.get("dead", [])),
                 len(by.get("inconclusive", []))))
        for verdict in ("live", "dead", "inconclusive"):
            if by.get(verdict):
                print("  %-13s %s" % (verdict + ":", ", ".join(sorted(by[verdict]))))
        if led["untested_columns"]:
            print("  never varied: " + ", ".join(led["untested_columns"]))
        print("wrote %s and %s" % (LEDGER.name, GUIDE.name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
