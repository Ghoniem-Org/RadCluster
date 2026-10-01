#!/usr/bin/env python
"""Validate the Saltelli estimator against Ishigami, whose indices are analytic.

    f(x) = sin(x1) + a sin^2(x2) + b x3^4 sin(x1),   x ~ U[-pi, pi]^3

with a = 7, b = 0.1 the exact indices are

    S1  = 0.3139   S2  = 0.4424   S3  = 0.0000
    ST1 = 0.5576   ST2 = 0.4424   ST3 = 0.2437

x3 is the interesting one: S3 = 0 exactly while ST3 = 0.244, i.e. it acts ONLY
through its interaction with x1.  An estimator that got S right but ST wrong
would pass a naive check and then silently screen out exactly the parameters
that matter through interactions -- which for this campaign would mean dropping
a loop-conversion parameter that is inert alone but live in combination.

Also tests PAIRWISE DELETION, the part that is easy to get wrong: some AB rows
are deliberately marked unusable and the indices must stay near truth, with
n_eff reporting the reduced count for the affected parameter only.

    python test_sobol.py
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from merge_and_sobol import sobol_indices   # noqa: E402

EXACT_S = {"x1": 0.3139, "x2": 0.4424, "x3": 0.0000}
EXACT_ST = {"x1": 0.5576, "x2": 0.4424, "x3": 0.2437}


def ishigami(x1, x2, x3, a=7.0, b=0.1):
    return np.sin(x1) + a * np.sin(x2) ** 2 + b * (x3 ** 4) * np.sin(x1)


def build(N=4096, seed=7):
    spec = {"version": "ishigami", "parameters": [
        {"id": i + 1, "key": f"x{i+1}", "sheet": "reactions", "tier": 2,
         "lo": -np.pi, "hi": np.pi, "nominal": 0.0, "prior": "U", "units": "-"}
        for i in range(3)], "fixed": []}
    td = Path(tempfile.mkdtemp())
    (td / "parameters.json").write_text(json.dumps(spec), encoding="utf-8")
    out = td / "d.csv"
    subprocess.run([sys.executable, str(HERE / "design.py"), "--tier", "2",
                    "--N", str(N), "--seed", str(seed), "--conditions", "C0",
                    "--spec", str(td / "parameters.json"), "--out", str(out)],
                   check=True, capture_output=True)
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    cols = lines[0].split(",")
    design, recs = [], {}
    for ln in lines[1:]:
        v = dict(zip(cols, ln.split(",")))
        rid = int(v["row_id"])
        design.append({"row_id": rid, "condition": v["condition"],
                       "base_idx": int(v["base_idx"]), "matrix": v["matrix"],
                       "param_j": int(v["param_j"])})
        recs[rid] = {"row_id": rid, "solver_rc": 0, "admissible": True,
                     "y": float(ishigami(float(v["x1"]), float(v["x2"]),
                                         float(v["x3"])))}
    return design, recs


def report(tag, r, drop_note=""):
    print(f"\n  {tag}{drop_note}")
    print(f"    {'param':>5s} {'S':>8s} {'exact':>8s} {'err':>8s} "
          f"{'ST':>8s} {'exact':>8s} {'err':>8s} {'n_eff':>6s}")
    worst = 0.0
    for k in ("x1", "x2", "x3"):
        v = r["indices"][k]
        eS, eT = abs(v["S"] - EXACT_S[k]), abs(v["ST"] - EXACT_ST[k])
        worst = max(worst, eS, eT)
        print(f"    {k:>5s} {v['S']:8.4f} {EXACT_S[k]:8.4f} {eS:8.4f} "
              f"{v['ST']:8.4f} {EXACT_ST[k]:8.4f} {eT:8.4f} {v['n_eff']:6d}")
    return worst


def main():
    N = 4096
    print(f"Ishigami validation of the Saltelli/Jansen estimators (N={N})")
    design, recs = build(N=N)
    r = sobol_indices(recs, design, ["x1", "x2", "x3"], "y", n_boot=0)
    w1 = report("all rows usable", r)

    # --- pairwise deletion: knock out 30 % of the AB rows for x2 only --------
    rng = np.random.default_rng(3)
    killed = 0
    for d in design:
        if d["matrix"] == "AB" and d["param_j"] == 1 and rng.random() < 0.30:
            recs[d["row_id"]]["admissible"] = False   # ran, but inadmissible
            killed += 1
    r2 = sobol_indices(recs, design, ["x1", "x2", "x3"], "y", n_boot=0)
    w2 = report("after dropping 30% of x2's AB rows", r2,
                f"  ({killed} rows marked inadmissible)")

    n = {k: r2["indices"][k]["n_eff"] for k in ("x1", "x2", "x3")}
    print(f"\n    n_eff after deletion: {n}")

    ok = True
    if w1 > 0.02:
        print(f"\n  FAIL: worst error {w1:.4f} > 0.02 with all rows usable"); ok = False
    else:
        print(f"\n  clean-case worst error {w1:.4f} < 0.02   OK")
    if w2 > 0.03:
        print(f"  FAIL: worst error {w2:.4f} > 0.03 after deletion"); ok = False
    else:
        print(f"  deletion-case worst error {w2:.4f} < 0.03   OK")
    if not (n["x2"] < n["x1"] and n["x1"] == n["x3"] == N):
        print(f"  FAIL: pairwise deletion leaked - only x2 should lose rows"); ok = False
    else:
        print(f"  pairwise deletion confined to x2 ({n['x2']}/{N}); "
              f"x1 and x3 untouched   OK")
    if abs(r["indices"]["x3"]["S"]) > 0.02 <= r["indices"]["x3"]["ST"]:
        pass
    if r["indices"]["x3"]["ST"] < 0.15:
        print(f"  FAIL: ST3 = {r['indices']['x3']['ST']:.3f}, interaction missed"); ok = False
    else:
        print(f"  x3 caught as interaction-only (S3~0, ST3="
              f"{r['indices']['x3']['ST']:.3f})   OK")
    print("\n  " + ("ESTIMATOR VALIDATED" if ok else "*** ESTIMATOR IS WRONG ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
