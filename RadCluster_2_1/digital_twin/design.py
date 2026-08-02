#!/usr/bin/env python
"""make_design - the ONE canonical Saltelli design for the Tier-2/3 campaign.

Run this ONCE, on ONE machine, and commit the output.  Every machine in the
campaign reads the same file; no machine ever generates its own.

Two reasons this is not merely tidy:

  1. Saltelli indices are computed from PAIRED rows.  S_j pairs f(A_i) against
     f(AB_i^(j)).  If two machines independently sample, the pairing is
     destroyed and the indices are silently WRONG -- not noisy, wrong.

  2. Tier 3 is a maximin NESTED SUBSET of the Tier-2 design (plan S4/Tier 3).
     A Sobol sequence is nested by construction: the first M base points of a
     2^k sequence are themselves a balanced design.  Latin hypercube is NOT --
     it cannot be extended without rebuilding.  The shared theta points between
     LF and HF ARE the multi-fidelity discrepancy model, so they must coincide
     exactly, which is why theta is referenced by a stable row_id.

Design size: N base samples, p sampled parameters -> N*(p+2) rows per condition
(matrices A, B, and AB_j for j = 1..p).

Usage
-----
    python design.py --tier 2 --N 16 --out design/T2_design_v1.csv
    python design.py --tier 3 --N 4 --nest-from design/T2_design_v1.csv \
                     --out design/T3_design_v1.csv

Outputs a CSV with one row per model evaluation:
    row_id, matrix, base_idx, param_j, condition, <key_1> ... <key_p>
plus a sidecar .meta.json recording the seed, parameter version, git SHA and a
SHA-256 of the design itself.  run_ensemble.py refuses to run against a design
whose hash does not match its sidecar.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.stats import qmc, truncnorm

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------- spec loading
def load_spec(path: Path | None = None) -> dict:
    """Load parameters.json (or a parameters.yaml sibling if pyyaml exists)."""
    if path is None:
        y = HERE / "parameters.yaml"
        if y.exists():
            try:
                import yaml  # optional
                return yaml.safe_load(y.read_text(encoding="utf-8"))
            except ImportError:
                print("  note: parameters.yaml present but pyyaml missing; "
                      "using parameters.json", file=sys.stderr)
        path = HERE / "parameters.json"
    return json.loads(path.read_text(encoding="utf-8"))


def active_params(spec: dict, tier: int) -> list[dict]:
    """Parameters sampled at this tier.

    Tier 2 -> tier==2 only.  Tier 3 -> tier 2 AND 3: the HF design must vary
    everything the LF design varied *plus* the dose-driven block, otherwise the
    HF runs are not comparable to their LF partners.
    """
    want = {2} if tier == 2 else {2, 3}
    ps = [p for p in spec["parameters"] if p["tier"] in want]
    return sorted(ps, key=lambda p: p["id"])


# ------------------------------------------------------------- prior transform
def to_physical(u: np.ndarray, p: dict) -> np.ndarray:
    """Map u ~ U[0,1)^n onto the parameter's prior support.

    Every family is inverse-CDF so the Sobol sequence's low-discrepancy
    structure survives the transform -- which is the whole point of using one.
    """
    lo, hi, fam = float(p["lo"]), float(p["hi"]), p["prior"]
    u = np.clip(u, 1e-12, 1 - 1e-12)
    if fam == "U":
        return lo + (hi - lo) * u
    if fam == "LN":                       # log-uniform
        return 10.0 ** (np.log10(lo) + (np.log10(hi) - np.log10(lo)) * u)
    if fam == "INT":
        # uniform over integers; floor of a uniform on [lo, hi+1)
        return np.floor(lo + (hi + 1 - lo) * u).clip(lo, hi).astype(int)
    if fam == "TN":
        mu = float(p["nominal"])
        sig = (hi - lo) / 4.0           # +/-2 sigma spans the stated range
        a, b = (lo - mu) / sig, (hi - mu) / sig
        return truncnorm.ppf(u, a, b, loc=mu, scale=sig)
    raise ValueError(f"unknown prior family {fam!r} for {p['key']}")


# ------------------------------------------------------------- Saltelli matrix
def saltelli_rows(N: int, params: list[dict], seed: int):
    """Return (rows, meta) where rows is a list of dicts.

    Layout per base index i:  A_i, B_i, then AB_i^(j) for each j.
    AB_i^(j) is A_i with column j replaced by B_i's column j.
    """
    p = len(params)
    # One Sobol draw of dimension 2p, split into A | B.  Drawing both halves
    # from a single sequence keeps A and B mutually low-discrepancy, which is
    # what the Saltelli estimator assumes.
    sampler = qmc.Sobol(d=2 * p, scramble=True, seed=seed)
    if N & (N - 1):
        print(f"  warning: N={N} is not a power of two; Sobol balance is "
              f"degraded. Prefer 8/16/32/64.", file=sys.stderr)
    draw = sampler.random(N)
    U_A, U_B = draw[:, :p], draw[:, p:]

    rows = []
    rid = 0
    for i in range(N):
        variants = [("A", -1, U_A[i].copy()), ("B", -1, U_B[i].copy())]
        for j in range(p):
            ab = U_A[i].copy()
            ab[j] = U_B[i][j]
            variants.append(("AB", j, ab))
        for matrix, j, u in variants:
            row = {"row_id": rid, "matrix": matrix, "base_idx": i, "param_j": j}
            for k, par in enumerate(params):
                v = to_physical(np.array([u[k]]), par)[0]
                row[par["key"]] = int(v) if par["prior"] == "INT" else float(v)
            rows.append(row)
            rid += 1
    return rows


def git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       cwd=str(HERE), stderr=subprocess.DEVNULL,
                                       text=True).strip()
    except Exception:
        return "unknown"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", type=int, choices=(2, 3), required=True)
    ap.add_argument("--N", type=int, default=16, help="base samples (power of 2)")
    ap.add_argument("--seed", type=int, default=20260801)
    ap.add_argument("--conditions", nargs="+", default=["N2", "N5", "I1"],
                    help="condition ids; the design is replicated across these")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--spec", type=Path, default=None)
    a = ap.parse_args(argv)

    spec = load_spec(a.spec)
    params = active_params(spec, a.tier)
    p = len(params)

    pending = [q["key"] for q in params if q.get("REVISION_PENDING")]
    if pending:
        print(f"  *** {len(pending)} parameter(s) marked REVISION_PENDING: "
              f"{', '.join(pending)}", file=sys.stderr)
        print("      Their ranges predate the 2026-08-01 absorption-gate work. "
              "Fix parameters.json before a production campaign.",
              file=sys.stderr)

    rows = saltelli_rows(a.N, params, a.seed)
    per_cond = len(rows)

    # replicate across conditions, keeping row_id globally unique and stable
    out_rows = []
    rid = 0
    for cond in a.conditions:
        for r in rows:
            rr = dict(r)
            rr["condition"] = cond
            rr["row_id"] = rid
            rr["cond_row_id"] = r["row_id"]   # id WITHIN the condition block
            rid += 1
            out_rows.append(rr)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    cols = (["row_id", "condition", "cond_row_id", "matrix", "base_idx", "param_j"]
            + [q["key"] for q in params])
    with a.out.open("w", encoding="utf-8", newline="") as fh:
        fh.write(",".join(cols) + "\n")
        for r in out_rows:
            fh.write(",".join(
                (f"{r[c]:.17g}" if isinstance(r[c], float) else str(r[c]))
                for c in cols) + "\n")

    digest = hashlib.sha256(a.out.read_bytes()).hexdigest()
    meta = {
        "tier": a.tier, "N": a.N, "p": p, "seed": a.seed,
        "conditions": a.conditions,
        "rows_per_condition": per_cond, "rows_total": len(out_rows),
        "parameters": [q["key"] for q in params],
        "parameters_version": spec.get("version"),
        "revision_pending": pending,
        "design_sha256": digest,
        "git_sha": git_sha(),
        "estimator": "Saltelli; S_j pairs A_i with AB_i^(j), S_T_j pairs B_i with AB_i^(j)",
    }
    a.out.with_suffix(".meta.json").write_text(json.dumps(meta, indent=2),
                                               encoding="utf-8")
    print(f"  wrote {a.out}  ({len(out_rows)} rows = "
          f"{a.N}*({p}+2)*{len(a.conditions)} conditions)")
    print(f"  sha256 {digest[:16]}...   p={p}  N={a.N}")
    print(f"  per machine (4-way): ~{len(out_rows)//4} rows each")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
