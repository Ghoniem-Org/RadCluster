#!/usr/bin/env python
"""campaign_layout - size a heterogeneous campaign and emit the exact commands.

WHY THIS EXISTS
---------------
run_ensemble already takes --weights (Sainte-Lague, interleaved).  What it does
NOT do is tell you what the weights should be, and getting them wrong is not
cosmetic: the v2 campaign split 1104 rows evenly over four machines whose real
throughput differed by ~3x, so the fast machines idled while the slow one held
a quarter of the design.  Worse, --weights MUST be byte-identical on every
participant -- a typo on one machine silently gives two machines the same rows
and leaves a hole somewhere else.

So the layout is computed ONCE, from measured throughput, written to
campaign_layout.json, and every participant is handed a command line generated
from that one file.

THROUGHPUT, NOT WALL TIME PER ROW
---------------------------------
The obvious measure -- mean wall_s per row -- is wrong under oversubscription.
A machine running 16 workers on 4 cores shows a 4x inflated wall_s per row while
its throughput is unchanged, because wall_s is elapsed time, not CPU time.  The
v2 pool shows exactly this ambiguity: MATRIX-PC2 looks 12-15x slower than the
Mac on the same three design rows, but only 2.7x slower per row averaged over
its whole set.  Those cannot both be capacity.

The measure that is invariant to oversubscription is

    throughput = rows_completed / session_wall_s        [rows/second]

which is what this module uses, read from each machine's manifest.  A machine
with no manifest cannot be measured and must declare its capacity instead --
loudly, so an assumption is never mistaken for a measurement.

    python campaign_layout.py --results results/ \
        --declare hoffman2=64 --declare MATRIX-PC2=4
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_manifests(res_dir: Path) -> dict:
    """Measured throughput per participant, from the manifests they wrote."""
    out = {}
    for f in sorted(glob.glob(str(res_dir / "*.manifest.json"))):
        try:
            m = json.loads(Path(f).read_text(encoding="utf-8"))
        except Exception:
            continue
        mach, tim = m.get("machine", {}), m.get("timing", {})
        rows = m.get("rows", {})
        done = rows.get("completed") or 0
        wall = tim.get("session_wall_s") or 0
        if not done or not wall:
            continue
        out[int(mach.get("index", -1))] = {
            "machine_id": mach.get("machine_id", "?"),
            "workers": mach.get("workers"),
            "cpu_count": mach.get("cpu_count"),
            "completed": done,
            "session_wall_s": wall,
            "rows_per_hour": 3600.0 * done / wall,
            "core_s_per_row": tim.get("row_wall_mean_s"),
            "source": "MEASURED",
        }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, default=HERE / "results")
    ap.add_argument("--design", default="design/T3_rev6.csv")
    ap.add_argument("--participant", action="append", default=[],
                    help="name:slots  -- one per participant, IN ORDER. 'slots' is "
                         "the number of concurrent single-core rows it can run. "
                         "Repeat for each. Example: --participant mac:14")
    ap.add_argument("--speed", action="append", default=[],
                    help="name:factor -- per-CORE speed relative to the reference "
                         "machine (default 1.0). Use when a participant's cores "
                         "are not the same speed, NOT to describe core count.")
    ap.add_argument("--split", action="append", default=[],
                    help="name:K -- expand this participant into K machine indices "
                         "of equal weight, for a scheduler ARRAY JOB. Each index "
                         "is an independent worker with its own results file, so "
                         "a task that is killed at the walltime limit costs only "
                         "its own in-flight rows. Required for a cluster: one "
                         "64-core index means one 64-core job, which queues badly "
                         "and dies whole; 16 x 4-core tasks queue fast and die "
                         "one at a time.")
    ap.add_argument("--out", type=Path, default=HERE / "campaign_layout.json")
    ap.add_argument("--n-rows", type=int, default=1008)
    ap.add_argument("--row-cost-s", type=float, default=3600.0,
                    help="measured seconds per row on the REFERENCE machine")
    a = ap.parse_args(argv)

    meas = read_manifests(a.results)
    if meas:
        print("MEASURED throughput from manifests (rows completed / session wall):")
        for k, v in sorted(meas.items()):
            print(f"   machine {k}  {v['machine_id']:16s} "
                  f"{v['workers']}w/{v['cpu_count']}c  "
                  f"{v['completed']:4d} rows in {v['session_wall_s']/3600:6.2f} h "
                  f"= {v['rows_per_hour']:6.2f} rows/h "
                  f"({v['rows_per_hour']/max(v['workers'] or 1,1):.2f} per worker)")
    else:
        print("no manifests found -- every capacity below is DECLARED, not measured")
    print()

    if not a.participant:
        ap.error("give at least one --participant name:slots")

    speeds = {}
    for s in a.speed:
        n, _, f = s.partition(":")
        speeds[n] = float(f)

    splits = {}
    for s in a.split:
        n, _, k = s.partition(":")
        splits[n] = int(k)

    parts = []
    for p in a.participant:
        name, _, slots = p.partition(":")
        slots = int(slots)
        sp = speeds.get(name, 1.0)
        k = splits.get(name, 1)
        if k < 1:
            ap.error(f"--split {name}:{k} must be >= 1")
        if slots % k:
            ap.error(f"--split {name}:{k} does not divide {slots} slots evenly; "
                     f"pick a K that divides {slots}")
        per = slots // k
        for i in range(k):
            parts.append({"name": f"{name}[{i+1}/{k}]" if k > 1 else name,
                          "group": name, "task": i + 1, "n_tasks": k,
                          "slots": per, "speed": sp,
                          "weight": round(per * sp, 4)})

    total_w = sum(p["weight"] for p in parts)
    # Wall-clock estimate: the campaign finishes when the LAST participant does,
    # and weights are chosen precisely so that is all of them at once.
    total_core_s = a.n_rows * a.row_cost_s
    eta_h = total_core_s / (total_w * 3600.0) if total_w else float("inf")

    print(f"{'participant':>16s} {'slots':>6s} {'speed':>6s} {'weight':>8s} "
          f"{'share':>7s} {'rows':>6s} {'ETA h':>7s}")
    for p in parts:
        share = p["weight"] / total_w
        rows = share * a.n_rows
        print(f"{p['name']:>16s} {p['slots']:6d} {p['speed']:6.2f} {p['weight']:8.2f} "
              f"{100*share:6.1f}% {rows:6.0f} {rows*a.row_cost_s/(p['slots']*p['speed']*3600):7.1f}")
    wstr = ",".join(str(p["weight"]) for p in parts)
    print(f"\n  total capacity {total_w:.2f} slot-equivalents; "
          f"{a.n_rows} rows x {a.row_cost_s/3600:.2f} h  ->  ETA {eta_h:.1f} h "
          f"({eta_h/24:.1f} days) if all participants run continuously")

    layout = {"design": a.design, "of": len(parts), "weights": wstr,
              "participants": parts, "n_rows": a.n_rows,
              "row_cost_s": a.row_cost_s, "eta_hours": round(eta_h, 2),
              "measured": {str(k): v for k, v in meas.items()}}
    a.out.write_text(json.dumps(layout, indent=2), encoding="utf-8")
    print(f"\n  wrote {a.out}")

    print("\n" + "=" * 78)
    print("  COMMANDS -- --weights and --of MUST be byte-identical on every")
    print("  participant. Copy from here; do not retype. A single typo gives two")
    print("  machines the same rows and leaves a hole somewhere else.")
    print("=" * 78)
    phys = " ".join([
        "--equations bin_moment --I <SET_BY_GATE_i3> --V 10000",
        "--i-discrete 50 --v-discrete 5 --i-bin 25 --v-bin 25",
        "--i-mobile-default 50 --v-mobile-default 5",
        "--dose 1.0 --rtol 1e-6 --solver-mode full_system"])
    seen = set()
    for k, p in enumerate(parts):
        grp = p.get("group", p["name"])
        if p.get("n_tasks", 1) > 1:
            if grp in seen:
                continue
            seen.add(grp)
            base = k          # index of task 1 of this group
            print(f"\n# {grp}: ARRAY JOB, {p['n_tasks']} tasks x {p['slots']} core(s), "
                  f"machine indices {base}..{base + p['n_tasks'] - 1}")
            print(f"#   in the job script, with SGE_TASK_ID running 1..{p['n_tasks']}:")
            print(f"MACHINE=$(( {base} + SGE_TASK_ID - 1 ))")
            print(f"python run_ensemble.py --design {a.design} \\")
            print(f"    --machine $MACHINE --of {len(parts)} --weights {wstr} \\")
            print(f"    --workers {p['slots']} \\")
            print(f"    {phys} \\")
            print(f"    --timeout-s <p99_HERE> --stop-after-s <h_rt MINUS one row>")
        else:
            print(f"\n# {p['name']}  (machine {k} of {len(parts)}, weight {p['weight']})")
            print(f"python run_ensemble.py --design {a.design} \\")
            print(f"    --machine {k} --of {len(parts)} --weights {wstr} \\")
            print(f"    --workers {p['slots']} \\")
            print(f"    {phys} \\")
            print(f"    --timeout-s <p99_ON_THIS_MACHINE>")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
