#!/usr/bin/env python3
"""Generate a SPACE-FILLING calibration design for the T8 experiment match.

WHY NOT SALTELLI.  T2/T3 used a Saltelli A/B/AB_j design because the question
was the RELATIVE RANKING of parameters (Sobol indices), which needs that exact
pairing and costs N*(p+2) rows to buy it.  The question now is different: find
parameter vectors whose loop character and mean size match the measurement.
That is estimation, not attribution -- the pairing buys nothing, and a Latin
hypercube covers the box far better per row.  At p = 7 a Saltelli design would
need 9 rows per base sample; this needs 1.

WHY THESE SEVEN.  Screened on the T3 campaign (1008 rows) against the two
observables the campaign never screened -- f_100_tem_1 and the loop diameters --
because <100> was withheld under bin_moment:

    parameter      rho(d_100)  rho(d_111)  rho(pile_100)  rho(f_100)
    L_hat            -0.540      +0.011       -0.699        +0.004
    E_a0_conv        +0.650      +0.268       +0.119        -0.583
    phi_max_junc     +0.322      +0.267       -0.019        -0.556
    Z_i              +0.420      +0.270       +0.015        -0.590
    f_cl_v           -0.337      -0.251       +0.136        +0.605
    f_cl_i           +0.343      -0.031       +0.284        +0.013
    dH2_abs_conv     +0.323      +0.113       +0.430        -0.183

L_hat is the one lever that is strong on the <100> axis and ORTHOGONAL to
d_111 (rho = +0.011), so it de-saturates the grid without disturbing the
1/2<111> size that sets the character.  Every other parameter trades the two
against each other.  It is sampled HIGH by default: all four grid-clean
full-dose rows in T3 had L_hat >= 1016, and the L_hat < 100 half of the campaign
has median pile_100 = 1.000, i.e. totally saturated.

THE TARGET is the 85 grid-clean rows (pile_100 < 0.05, 8.4 % of T3), where the
model already reproduces the experiment qualitatively: median f_100_tem_1 =
0.331 against a "< 0.5" datum, median d_111 = 2.31 nm against 3.4 nm, and
d_100 down to 1.61 nm.  This design searches that region rather than the full
prior box.

Every non-sampled parameter is pinned to a BASE ROW rather than the workbook, so
the design sits at a point already known to be grid-clean, conserving and cheap
(T3 row 38: pile_100 = 1e-17, delta_FP = 4.6e-3, 1852 s to full dose).

Usage:
    python make_calib_design.py --n 64 --condition C_330_G6 --out design/T8_calib.csv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent

# key -> (lo, hi, scale).  "log" samples uniformly in log10.
SAMPLED = {
    "L_hat":        (300.0,  3500.0, "log"),
    "E_a0_conv":    (1.6,    2.4,    "lin"),
    "phi_max_junc": (0.1,    1.0,    "lin"),
    "Z_i":          (1.02,   1.15,   "lin"),
    "f_cl_v":       (0.2,    0.7,    "lin"),
    "f_cl_i":       (0.05,   0.25,   "lin"),
    "dH2_abs_conv": (0.26,   0.45,   "lin"),
}


def latin_hypercube(n: int, d: int, rng: random.Random) -> list[list[float]]:
    """Stratified sample in [0,1)^d: one point per stratum per dimension.

    Plain uniform sampling leaves gaps and clumps that matter at n = 64 over 7
    dimensions; the LHS guarantees each parameter's marginal is covered evenly,
    which is what makes a small design usable for estimation.
    """
    cols = []
    for _ in range(d):
        strata = [(i + rng.random()) / n for i in range(n)]
        rng.shuffle(strata)
        cols.append(strata)
    return [[cols[j][i] for j in range(d)] for i in range(n)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=64)
    ap.add_argument("--condition", required=True,
                    help="condition key; must exist in the --conditions file "
                         "passed to run_ensemble")
    ap.add_argument("--seed", type=int, default=20260813)
    ap.add_argument("--base-design", type=Path, default=HERE / "design/T3_rev6.csv")
    ap.add_argument("--base-row", default="38")
    ap.add_argument("--base-condition", default="N2")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    base_rows = [r for r in csv.DictReader(a.base_design.open())
                 if r["row_id"] == a.base_row and r["condition"] == a.base_condition]
    if not base_rows:
        raise SystemExit(f"base row {a.base_row}/{a.base_condition} not found in "
                         f"{a.base_design}")
    base = base_rows[0]
    params = json.loads((HERE / "design/T3_rev6.meta.json").read_text())["parameters"]
    missing = [k for k in SAMPLED if k not in params]
    if missing:
        raise SystemExit(f"sampled keys absent from the parameter list: {missing}")

    rng = random.Random(a.seed)
    keys = list(SAMPLED)
    pts = latin_hypercube(a.n, len(keys), rng)

    cols = ["row_id", "condition", "cond_row_id", "matrix", "base_idx",
            "param_j"] + params
    out = []
    for i, pt in enumerate(pts):
        row = {"row_id": 800 + i, "condition": a.condition, "cond_row_id": 800 + i,
               "matrix": "L", "base_idx": i, "param_j": -1}
        for p in params:
            row[p] = base[p]                      # pinned unless sampled below
        for u, k in zip(pt, keys):
            lo, hi, scale = SAMPLED[k]
            if scale == "log":
                row[k] = repr(10.0 ** (math.log10(lo) + u * (math.log10(hi)
                                                             - math.log10(lo))))
            else:
                row[k] = repr(lo + u * (hi - lo))
        out.append(row)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out)

    sha = hashlib.sha256(a.out.read_bytes()).hexdigest()
    meta = {
        "tier": "T8_calibration",
        "kind": "latin_hypercube",
        "n": a.n,
        "seed": a.seed,
        "condition": a.condition,
        "sampled": {k: {"lo": v[0], "hi": v[1], "scale": v[2]}
                    for k, v in SAMPLED.items()},
        "base_row": f"{a.base_design.name}:{a.base_row}/{a.base_condition}",
        "parameters": params,
        "design_sha256": sha,
        "note": ("Space-filling, NOT Saltelli: this design is for ESTIMATION "
                 "(match the measurement), not attribution (Sobol indices), so "
                 "it carries no A/B/AB_j pairing and merge_and_sobol must not be "
                 "pointed at it."),
    }
    a.out.with_suffix(".meta.json").write_text(json.dumps(meta, indent=1))
    print(f"wrote {a.out}  ({a.n} rows, {len(keys)} sampled, condition={a.condition})")
    print(f"  sha256 {sha[:16]}")
    for k, (lo, hi, sc) in SAMPLED.items():
        vals = [float(r[k]) for r in out]
        print(f"  {k:14s} [{lo:g}, {hi:g}] {sc:3s} -> "
              f"min={min(vals):.4g} max={max(vals):.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
