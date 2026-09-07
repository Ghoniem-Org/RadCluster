#!/usr/bin/env python3
"""
rescore.py — re-apply the observation operator to runs already computed.

    python rescore.py --dry-run     # show what would change
    python rescore.py               # rewrite and publish the corrected records

WHY.  The dose ladder scored the wrong dose and said it had not (b5e9a51): the
end point rounded past the final grid point, and a requested dose that was not
itself a grid point silently took the last point BELOW it, however far below.
On the exact arm `T4_D_4k` the "0.2 dpa" column carried 0.0665 dpa numbers -- 3x
low -- while `reached_all_scoring_doses` reported True.  Every run computed
before that fix carries columns labelled with doses they were not taken at.

The trajectories are sound; only the extraction was wrong.  So this re-reads
each run's saved solution and re-runs `run_ensemble.observe` under the current
code, rather than recomputing 17 hours of solver time.

SCOPE.  A run can only be re-scored on the machine that produced it: the
trajectory lives in `output/<stamp>_<run_id>_.../plots/plot_data.pkl`, and
`output/` is gitignored (a 38 MB pickle per run, and one of live objects that
will not survive a class refactor -- see .gitignore).  What crossed the wire is
the compact record, which is the thing being corrected.  Run this on EVERY
machine that computed rows; each fixes its own and pushes.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
DT = HERE.parent
MOD = DT.parent
REPO = MOD.parent
for p in (str(REPO), str(DT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import runs as manifest_mod                     # noqa: E402
from gitsync import CampaignSync                # noqa: E402


def find_runs():
    """Map run_id -> newest output directory produced on this machine."""
    ids = {e["run_id"] for e in manifest_mod.manifest()}
    found = {}
    for d in sorted((MOD / "output").glob("*/"), reverse=True):
        if not (d / "plots" / "plot_data.pkl").exists():
            continue
        name = d.name
        # <stamp>_<run_id>_<solver_mode>_... -- match the longest id that fits,
        # so T4_C4_P1_4k is not mistaken for T4_C4_4k.
        for rid in sorted(ids, key=len, reverse=True):
            if f"_{rid}_" in name and rid not in found:
                found[rid] = d
                break
    return found


def rescore(rid, run_dir, entry):
    """Recompute the observables for one run from its saved solution."""
    from RadCluster_2_1.py_utils import visualization as viz
    import run_ensemble as re_mod

    d = viz.load_plot_data(str(run_dir / "plots" / "plot_data.pkl"))
    res, idata = d["results"], d["input_data"]
    sim = types.SimpleNamespace(input_data=idata)

    cfg = {"I": entry["I"], "V": entry["V"], "dose": float(entry["dose"]),
           "C_floor": float(idata.reactions.get("C_floor", 1e-15))}
    out = re_mod.observe(res, sim, cfg, 1.0)

    reads = entry["dose_read"]
    reads = reads if isinstance(reads, (list, tuple)) else [reads]
    lad = out.get("at_dose") or {}
    scored = {}
    for target in reads:
        key = f"{float(target):g}"
        got = lad.get(key)
        scored[key] = got if got else {"missing": True,
                                       "dose_reached": out.get("dose_reached")}
    out["scored"] = scored
    out["reached_all_scoring_doses"] = all(
        not v.get("missing") for v in scored.values())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would change; publish nothing")
    ap.add_argument("--only", default=None, help="restrict to one run_id")
    ap.add_argument("--branch", default="campaign-verification")
    a = ap.parse_args(argv)

    entries = {e["run_id"]: e for e in manifest_mod.manifest()}
    found = find_runs()
    if a.only:
        found = {k: v for k, v in found.items() if k == a.only}
    if not found:
        print("no re-scoreable runs on this machine "
              "(need output/<...>/plots/plot_data.pkl)")
        return 0

    sync = CampaignSync(REPO, branch=a.branch)
    sync.ensure()
    claims = sync.claims()

    print(f"{len(found)} run(s) with a saved solution on this machine\n")
    changed = []
    for rid in sorted(found):
        entry, run_dir = entries[rid], found[rid]
        try:
            new = rescore(rid, run_dir, entry)
        except Exception as exc:
            print(f"  {rid:14s} RESCORE FAILED: {type(exc).__name__}: {exc}")
            continue
        old = claims.get(rid, {})
        old_sc, new_sc = old.get("scored") or {}, new["scored"]
        moved = []
        for k in sorted(new_sc, key=float):
            o, n = old_sc.get(k) or {}, new_sc[k]
            od, nd = o.get("dose"), n.get("dose")
            if od is None or nd is None:
                moved.append(f"{k}: {'MISSING' if o.get('missing') else od} -> "
                             f"{'MISSING' if n.get('missing') else nd}")
            elif abs(od - nd) > 1e-9 * max(abs(od), 1.0):
                moved.append(f"{k} dpa: scored at {od:.4g} -> {nd:.4g}"
                             f"  (d100 {o.get('d_100_nm', float('nan')):.3f}"
                             f" -> {n.get('d_100_nm', float('nan')):.3f} nm)")
        flag = "" if new["reached_all_scoring_doses"] else "   [did not reach all]"
        print(f"  {rid:14s} {'CHANGED' if moved else 'unchanged'}{flag}")
        for m in moved:
            print(f"                   {m}")
        if moved:
            changed.append((rid, new))

    if a.dry_run:
        print(f"\ndry run: {len(changed)} record(s) would be republished")
        return 0

    for rid, new in changed:
        rec = dict(claims.get(rid, {}))
        rec.update(new)
        rec["rescored_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        rec["rescore_reason"] = ("dose ladder fix b5e9a51: end-point rounding "
                                 "and off-grid checkpoints scored the wrong dose")
        ok = sync.publish(rid, rec, {"observables.json":
                                     json.dumps(rec, indent=2, sort_keys=True,
                                                default=str) + "\n"})
        print(f"  republished {rid}: {'ok' if ok else 'PUSH FAILED'}")
    print(f"\n{len(changed)} record(s) corrected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
