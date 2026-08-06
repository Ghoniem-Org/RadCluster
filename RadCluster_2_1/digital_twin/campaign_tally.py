#!/usr/bin/env python3
"""Tally a multi-machine campaign: who has reported what, and is it poolable.

Run this after `git pull` on any participating machine.  It answers the three
questions that matter between pulls:

  1. COVERAGE  - how many of the design's rows exist, per machine, against the
     count that machine OWNS under the frozen weights.  A machine that has not
     pushed is indistinguishable from a machine that has not started, so the
     table says "no file" rather than "0 done" for the former.

  2. POOLABILITY - whether every reported row carries the same git/solver/
     workbook/design/run_cfg hashes.  A split here means the rows cannot be
     pooled into one Sobol estimate, and it is silent unless something checks.
     This reuses merge_and_sobol.check_provenance so the two agree by
     construction.

  3. THROUGHPUT - completions per hour, computed from the rows themselves.
     NOT mean wall over completed rows: with W workers running concurrently,
     mean wall over the rows that HAVE finished is biased low early on, because
     the slow rows are precisely the ones still running.  Every ETA this
     campaign got wrong was got wrong that way.  Here, throughput is

         (rows done) / (elapsed wall on that machine)

     where elapsed comes from the machine's own first-to-last row span plus the
     wall of the longest row still in flight.  It converges from below and
     cannot flatter itself.

Usage:
    python campaign_status.py                      # frozen design in machines.json
    python campaign_status.py --design design/T3_rev6.csv
    python campaign_status.py --at-dose 0.1        # ladder view at one rung
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent

import run_ensemble as RE
from merge_and_sobol import check_provenance


def load_by_machine(res_dir: Path, stem: str):
    """{machine_index: {row_id: rec}}, plus the files each came from.

    Keyed on the machine index recorded IN the row, not on the filename: a
    file renamed or copied between hosts would otherwise be miscounted, and
    the subtask suffix (_t3) already makes one machine write many files.
    """
    per = defaultdict(dict)
    files = defaultdict(list)
    for f in sorted(res_dir.glob(f"{stem}_machine*.jsonl")):
        n = 0
        for ln in f.read_text(encoding="utf-8").splitlines():
            if not ln.strip():
                continue
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue                    # truncated last line after a kill
            k = r.get("machine")
            if k is None:                   # older rows predate the field
                k = _machine_from_name(f.name, stem)
            per[k][r["row_id"]] = r
            n += 1
        files[_machine_from_name(f.name, stem)].append((f.name, n))
    return per, files


def load_starts(res_dir: Path, stem: str):
    """{machine_index: marker} from the .started.json a live worker drops."""
    out = {}
    for f in sorted(res_dir.glob(f"{stem}_machine*.started.json")):
        try:
            m = json.loads(f.read_text())
        except Exception:
            continue
        k = m.get("machine")
        if k is None:
            continue
        # Many subtasks of one machine each drop a marker; the earliest is the
        # one that bounds the elapsed wall for all of them.
        if k not in out or float(m["unix"]) < float(out[k]["unix"]):
            out[k] = m
    return out


def _machine_from_name(name: str, stem: str) -> int:
    tail = name[len(stem):].lstrip("_")
    digits = ""
    for ch in tail[len("machine"):] if tail.startswith("machine") else "":
        if ch.isdigit():
            digits += ch
        else:
            break
    return int(digits) if digits else -1


def owned_counts(design: Path, of: int, weights):
    """rows each machine index is responsible for, and the full row_id set."""
    ids = []
    with design.open() as fh:
        for row in csv.DictReader(fh):
            ids.append(int(row["row_id"]))
    own = defaultdict(list)
    for rid in ids:
        own[RE.assign_machine(rid, of, weights)].append(rid)
    return own, set(ids)


def fmt_eta(hours):
    if hours is None:
        return "     -"
    if hours > 999:
        return "  >999h"
    return f"{hours:6.1f}h"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", type=Path, default=None,
                    help="default: the design named in machines.json")
    ap.add_argument("--registry", type=Path, default=HERE / "machines.json")
    ap.add_argument("--results", type=Path, default=HERE / "results")
    ap.add_argument("--at-dose", type=float, default=None,
                    help="also report how many rows reached this dose rung")
    ap.add_argument("--sync", action="store_true",
                    help="first send THIS machine's results and pull everyone "
                         "else's, then tally (git is the transport)")
    ap.add_argument("--no-push", action="store_true",
                    help="with --sync: pull and commit, but do not push")
    a = ap.parse_args(argv)

    if a.sync:
        import campaign_ops as ops
        print("sync")
        try:
            ops.sync_results(push=not a.no_push)
        except Exception as e:
            # A sync failure must not hide the tally: the local rows are still
            # worth seeing, and "could not reach the remote" is a different
            # problem from "the campaign is not progressing".
            print(f"  *** sync failed ({e}); tallying LOCAL results only")

    reg = json.loads(a.registry.read_text())
    design = a.design or (HERE / reg["design"])
    of = int(reg["of"])
    weights = [float(x) for x in str(reg["weights"]).split(",")]
    names = {int(m["index"]): m["name"] for m in reg["machines"]}

    own, all_ids = owned_counts(design, of, weights)
    per, files = load_by_machine(a.results, design.stem)
    starts = load_starts(a.results, design.stem)

    # A start marker is only a clock for the machine we are ON.  Another
    # machine's marker arrives over git with ITS unix epoch; subtracting our
    # local now from it measures the pull latency, not that machine's runtime.
    # So exact throughput is offered for this host only, and everyone else
    # falls back to the ceiling estimator -- which the table labels as such.
    now = time.time()
    try:
        here_k = RE.detect_machine(reg)["index"]
    except SystemExit:
        here_k = None
    starts = {k: v for k, v in starts.items() if k == here_k}

    print(f"\ndesign   {design.name}   {len(all_ids)} rows")
    print(f"registry {a.registry.name}   of={of}  weights={reg['weights']}  "
          f"timeout_s={reg.get('timeout_s')}")
    print()

    hdr = (f"{'#':>2} {'machine':<18} {'owns':>5} {'done':>5} {'%':>6} "
           f"{'cut':>4} {'rc!=0':>6} {'med_s':>7} {'rows/h':>7} {'ETA':>7}  status")
    print(hdr)
    print("-" * len(hdr))

    tot_done = tot_owned = 0
    eta_hours = []
    for k in range(of):
        n_own = len(own.get(k, []))
        tot_owned += n_own
        recs = per.get(k, {})
        n_done = len(recs)
        tot_done += n_done

        if not files.get(k) and not recs:
            print(f"{k:>2} {names.get(k,'?'):<18} {n_own:>5} {'-':>5} {'-':>6} "
                  f"{'-':>4} {'-':>6} {'-':>7} {'-':>7} {'-':>7}  NO FILE - not pushed / not started")
            continue

        walls = sorted(r.get("wall_s", 0.0) for r in recs.values())
        med = walls[len(walls) // 2] if walls else 0.0
        cut = sum(1 for r in recs.values() if r.get("starved"))
        bad = sum(1 for r in recs.values() if r.get("solver_rc"))
        nw = max((int(r.get("workers") or 1)) for r in recs.values()) if recs else 1

        # THROUGHPUT.  Completions divided by REAL elapsed wall, taken from the
        # .started.json marker the worker drops at launch.  The tempting
        # alternative -- sum(row_wall)/workers -- counts only rows that have
        # landed, so early in a run, when the slow rows are still occupying
        # their workers and contributing to neither numerator nor denominator,
        # it overstates the rate by however bimodal the cost distribution is.
        # On this campaign that was better than 2x.  Where no marker exists
        # (a run started before this field, or another machine's file pulled
        # without it) the rate is printed with a '~' and IS that biased
        # estimator; treat it as a ceiling, not an estimate.
        st = starts.get(k)
        exact = st is not None and now is not None
        if exact:
            elapsed_h = max(now - float(st["unix"]), 1.0) / 3600.0
            base = int(st.get("rows_already_done") or 0)
            n_this = max(n_done - base, 0)      # rows THIS session produced
            rate = (n_this / elapsed_h) if elapsed_h > 0 else None
        else:
            burned_h = (sum(walls) / nw) / 3600.0 if walls else 0.0
            rate = (n_done / burned_h) if burned_h > 0 else None
        eta = ((n_own - n_done) / rate) if rate else None
        if eta is not None and exact:
            eta_hours.append(eta)

        pct = 100.0 * n_done / n_own if n_own else 0.0
        note = f"{nw}w"
        # STARTUP TRANSIENT.  For the first couple of row-times the workers are
        # all still on their first row, so completions/elapsed reads far below
        # the steady state it will settle to -- the mirror image of the
        # sum(wall)/workers bias, and just as misleading if quoted as a rate.
        # Flagged until elapsed covers ~2 median rows, after which the pipeline
        # is full and the division means what it says.
        if exact and med > 0 and elapsed_h * 3600.0 < 2.0 * med:
            note += (f"  TRANSIENT: {elapsed_h*3600/med:.1f} row-times in, "
                     f"rate not yet meaningful")
        if not exact and rate:
            note += "  ~rate: no start marker, CEILING only"
        if bad:
            note += f"  {bad} FAIL"
        if n_done >= n_own:
            note += "  COMPLETE"
        mark = "" if exact else "~"
        print(f"{k:>2} {names.get(k,'?'):<18} {n_own:>5} {n_done:>5} {pct:>5.1f}% "
              f"{cut:>4} {bad:>6} {med:>7.0f} "
              f"{(f'{mark}{rate:.1f}' if rate else '-'):>7} {fmt_eta(eta)}  {note}")

    print("-" * len(hdr))
    print(f"{'':>2} {'TOTAL':<18} {tot_owned:>5} {tot_done:>5} "
          f"{100.0*tot_done/tot_owned if tot_owned else 0:>5.1f}%")

    # Partition sanity: every design row owned exactly once, and nothing
    # reported that the design does not contain.
    seen = {}
    dup = []
    for k, recs in per.items():
        for rid in recs:
            if rid in seen:
                dup.append((rid, seen[rid], k))
            seen[rid] = k
    missing = len(all_ids - set(seen))
    print(f"\ncoverage  {len(seen)}/{len(all_ids)} distinct rows   "
          f"{missing} not yet reported"
          + (f"   *** {len(dup)} row(s) computed by TWO machines" if dup else ""))
    if dup[:5]:
        for rid, k1, k2 in dup[:5]:
            print(f"    row {rid}: machine {k1} and machine {k2}")

    pool = {}
    for recs in per.values():
        pool.update(recs)
    if pool:
        print("\nprovenance")
        check_provenance(pool)
        print("  (no output above = all reported rows are poolable)")

        rungs = defaultdict(int)
        for r in pool.values():
            for d in (r.get("at_dose") or {}):
                rungs[float(d)] += 1
        if rungs:
            print("\ndose ladder (rows reaching each rung, of "
                  f"{len(pool)} reported)")
            line = "  " + "  ".join(f"{d:g}:{rungs[d]}"
                                    for d in sorted(rungs))
            print(line)
            full = max(rungs) if rungs else 0
            common = [d for d in sorted(rungs) if rungs[d] == len(pool)]
            print(f"  highest rung ALL reported rows reached: "
                  f"{max(common):g}" if common else
                  "  no rung reached by every reported row")

        if a.at_dose is not None:
            key = f"{a.at_dose:g}"
            n = sum(1 for r in pool.values() if key in (r.get("at_dose") or {}))
            print(f"\n  at dose {key}: {n}/{len(pool)} reported rows have a value")

    if eta_hours:
        print(f"\nETA to finish the design: {max(eta_hours):.1f} h "
              f"(slowest machine; optimistic - excludes rows in flight)")
    print()


if __name__ == "__main__":
    main()
