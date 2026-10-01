#!/usr/bin/env python
"""Close the self-learning loop: turn UNVERIFIED candidates into measurements.

`learn.py` will not score an observable that has not survived a change of grid
EXTENT (see learn.EXTENT_TOL).  That is the right rule -- refinement at fixed
extent proved nothing, and occupancy proved less -- but on its own it only
*labels* the problem.  This script is the other half: it finds the rows the
ledger most wants verified and emits the runs that would verify them.

    python verify.py                  # show what needs verifying, and the commands
    python verify.py --run --workers 6

THE ONE INVARIANT THAT MATTERS.  `run_ensemble` builds `theta_hash` from every
design column EXCEPT row_id / cond_row_id / matrix / base_idx / param_j /
condition -- and `theta_id` IS included.  A verification row must therefore be a
byte-for-byte copy of the original design row with ONLY row_id and cond_row_id
changed.  Renumber `theta_id` and the hash changes, the pair never forms, and
the run is silently wasted.  This is not hypothetical: the V12/V13/V15 designs
renumber theta_id and so cannot pair with their own parents.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEDGER = HERE / "calibration_ledger.json"
RESULTS = HERE / "results"

# Columns run_ensemble excludes from theta_hash; only these may be rewritten.
REWRITABLE = ("row_id", "cond_row_id")


def recover_V(row: dict):
    occ, mn = row.get("occ_v"), row.get("mean_n_v")
    if not occ or mn is None or occ <= 0:
        return None
    return int(round(mn / occ))


def load_result_rows() -> dict:
    """row_id -> raw result row (last one wins; they are replicates)."""
    out = {}
    for f in sorted(RESULTS.glob("*.jsonl")):
        for line in f.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            out.setdefault(str(r.get("row_id")), r)
    return out


def candidates(led: dict, top: int) -> list:
    """Unverified rows, best first.

    Ranked by how much a verification would TELL US: rows already sitting in
    several bands are the ones whose status the campaign is actually resting on,
    so they are the ones worth the CPU.
    """
    out = []
    for stem, st in led.get("stages", {}).items():
        if not st.get("in_scope"):
            continue
        for r in st.get("rows", []):
            if not r.get("valid") or r.get("n_verified", 0) > 0:
                continue
            out.append({"stage": stem, "design": st.get("design"),
                        "row_id": str(r.get("row_id")), "label": r.get("label"),
                        "n_in_range": r.get("n_in_range", 0),
                        "log_distance": r.get("log_distance", 1e9)})
    out.sort(key=lambda c: (-c["n_in_range"], c["log_distance"]))
    return out[:top]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=8,
                    help="how many unverified rows to verify (default 8)")
    ap.add_argument("--factor", type=int, default=4,
                    help="grid-extent multiplier (default 4)")
    ap.add_argument("--run", action="store_true", help="launch, don't just print")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tag", default="AUTOVERIFY")
    a = ap.parse_args(argv)

    if not LEDGER.exists():
        print("no ledger - run learn.py first", file=sys.stderr)
        return 2
    led = json.loads(LEDGER.read_text(encoding="utf-8"))
    res = load_result_rows()
    cands = candidates(led, a.top)
    if not cands:
        print("nothing to verify: every valid in-scope row already has an "
              "extent pair.")
        return 0

    # Group by (design file, original grid config).  Rows run on different grids
    # cannot share one verification command.
    groups: dict = {}
    skipped = []
    for c in cands:
        raw = res.get(c["row_id"])
        if raw is None or not c["design"]:
            skipped.append((c["row_id"], "no raw result or no design file"))
            continue
        V = recover_V(raw)
        if not V:
            skipped.append((c["row_id"], "cannot recover V (no occ_v)"))
            continue
        key = (c["design"], V, raw.get("bin_V_bin"), raw.get("bin_I_bin"),
               raw.get("bin_i_discrete"), raw.get("bin_v_discrete"),
               raw.get("equations"))
        groups.setdefault(key, []).append((c, raw))

    print(f"{len(cands)} unverified candidates -> {len(groups)} verification run(s)\n")
    cmds = []
    for gi, ((design, V, vbin, ibin, idisc, vdisc, eqs), items) in \
            enumerate(sorted(groups.items(), key=lambda kv: str(kv[0]))):
        src = HERE / design
        if not src.exists():
            skipped.append((design, "design file missing"))
            continue
        by_id = {r["row_id"]: r for r in csv.DictReader(src.open(encoding="utf-8"))}
        hdr = list(next(iter(by_id.values())).keys())
        newV = V * a.factor
        # Hold the bin RATIO fixed while the extent grows, so the test isolates
        # EXTENT from RESOLUTION.  r = (V/v_discrete)**(1/v_bin).
        import math
        vd = max(int(vdisc or 5), 1)
        if vbin and V > vd:
            r_bin = (V / vd) ** (1.0 / int(vbin))
            newvbin = max(2, int(round(math.log(newV / vd) / math.log(r_bin))))
        else:
            newvbin = vbin
        out_rows, picked = [], []
        for off, (c, _raw) in enumerate(items):
            drow = by_id.get(c["row_id"])
            if drow is None:
                skipped.append((c["row_id"], f"row not in {design}"))
                continue
            nr = dict(drow)                     # VERBATIM copy ...
            nid = str(900000 + gi * 100 + off)  # ... except the two safe fields
            nr["row_id"] = nid
            if "cond_row_id" in nr:
                nr["cond_row_id"] = nid
            out_rows.append(nr)
            picked.append(f"{c['row_id']}({c['label'] or '-'}, {c['n_in_range']}/6)")
        if not out_rows:
            continue
        name = f"{a.tag}_g{gi}"
        dpath = HERE / "design" / f"{name}.csv"
        with dpath.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=hdr)
            w.writeheader()
            w.writerows(out_rows)
        cmd = ["python", "run_ensemble.py",
               "--design", f"design/{name}.csv",
               "--conditions", "conditions_S8.json", "--spec", "parameters_S4.json",
               "--out", f"results/{name}.jsonl",
               "--machine", "0", "--of", "1",
               "--equations", str(eqs or "bin_moment"),
               "--i-discrete", str(idisc or 100), "--i-bin", str(ibin or 36),
               "--v-discrete", str(vd), "--v-bin", str(newvbin),
               "--allow-mixed", "--I", "80000", "--V", str(newV),
               "--dose", "15.0", "--lnl", "1", "--rtol", "1e-5",
               "--solver-mode", "full_system", "--timeout-s", "43200",
               "--workers", str(a.workers), "--omp-threads", "1"]
        cmds.append(cmd)
        print(f"[{name}] V {V} -> {newV} (v_bin {vbin} -> {newvbin}, ratio held)")
        for p in picked:
            print(f"    {p}")
        print("    " + " ".join(cmd) + "\n")

    if skipped:
        print("skipped:")
        for rid, why in skipped:
            print(f"    {rid}: {why}")

    if a.run:
        for cmd in cmds:
            print("launching:", " ".join(cmd[-6:]))
            subprocess.Popen(cmd, cwd=HERE,
                             stdout=open(HERE / "results" / "verify.log", "a"),
                             stderr=subprocess.STDOUT)
        print(f"\n{len(cmds)} run(s) launched.  Re-run learn.py when they finish; "
              "pairs form automatically.")
    else:
        print("dry run - pass --run to launch.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
