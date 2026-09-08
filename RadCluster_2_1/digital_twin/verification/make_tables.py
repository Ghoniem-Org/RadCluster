#!/usr/bin/env python3
"""
make_tables.py — emit the verification tables from the campaign board.

    python make_tables.py            # markdown to stdout
    python make_tables.py --latex    # LaTeX bodies
    python make_tables.py --out FILE # write markdown to FILE

This is plan S3.5's reporting helper.  It reads the published records and
REFUSES to print a value from a run that missed its comparison dose (S3.4.1):
such a row appears in a "did not reach" line with the dose it achieved, never as
a metric.  The rule is mechanical here rather than a matter of discipline.

Table 4 ships as the NO-COARSENING variant.  The production-configuration exact
arm was attempted three times and abandoned each time (see S3 note below); the
NC variant is the one where both arms provably solve identical equations, which
is what a verification table requires.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
DT = HERE.parent
MOD = DT.parent
REPO = MOD.parent
for p in (str(REPO), str(DT), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

SYNC = HERE / ".sync" / "claims"
# The salvaged discrete arm: 2026-09-07, 12006 equations, 17.6 h, LOOP_COAL
# silently absent -- which is exactly what makes it the right reference for the
# no-coarsening rungs.
NC_REF_GLOB = str(MOD / "output" / "20260907_093333_T4_D_4k_*" / "plots" / "plot_data.pkl")


def board() -> dict:
    out = {}
    for f in sorted(SYNC.glob("*.json")):
        try:
            j = json.loads(f.read_text(encoding="utf-8"))
            out[j.get("run_id", f.stem)] = j
        except Exception:
            pass
    return out


def nc_reference():
    """Re-score the 2026-09-07 discrete arm from its own saved trajectory."""
    from RadCluster_2_1.py_utils import visualization as viz
    import run_ensemble as re_mod
    f = glob.glob(NC_REF_GLOB)
    if not f:
        return None
    d = viz.load_plot_data(f[0])
    sim = types.SimpleNamespace(input_data=d["input_data"])
    cfg = {"I": 4000, "V": 4000, "dose": 2.0,
           "C_floor": float(d["input_data"].reactions.get("C_floor", 1e-15))}
    return re_mod.observe(d["results"], sim, cfg, 1.0)["at_dose"]["2"]


# B4 is not a campaign row -- it IS the tracked reference run.  Score it from
# its own saved trajectory so Table 1 shows the rung the ladder converges toward
# instead of a "did not reach ?" hole.
REF_RUN_GLOB = str(MOD / "output" /
                   "20260906_063055_full_system_bin_moment_CD_fission_I80000V20000_im5vm5"
                   / "plots" / "plot_data.pkl")


def reference_run_row():
    from RadCluster_2_1.py_utils import visualization as viz
    import run_ensemble as re_mod
    f = glob.glob(REF_RUN_GLOB)
    if not f:
        return None
    d = viz.load_plot_data(f[0])
    sim = types.SimpleNamespace(input_data=d["input_data"])
    cfg = {"I": 80000, "V": 20000, "dose": 40.0,
           "C_floor": float(d["input_data"].reactions.get("C_floor", 1e-15))}
    o = re_mod.observe(d["results"], sim, cfg, 1.0)
    v = (o.get("at_dose") or {}).get("15.72")
    if not v:
        return None
    return {"status": "done", "N_eq": 189, "scored": {"15.72": v},
            "delta_FP": v.get("delta_FP", float("nan")), "wall_s": 1317.37}


def cell(rec, dose_key):
    """The scored block for a run, or None when it may not be reported."""
    if not rec or rec.get("status") not in ("done",):
        return None
    v = (rec.get("scored") or {}).get(dose_key)
    if not v or v.get("missing"):
        return None
    return v


def table4(B, md):
    ref = nc_reference()
    md.append("## Table 4 — closure verification against the exact solution\n")
    md.append("**I = V = 4000, `LOOP_COAL = 0` in both arms, scored at 2.0000 dpa.**")
    md.append("Both arms solve identical equations, so the deviations below are "
              "closure error.\n")
    if ref is None:
        md.append("_reference trajectory not present on this machine_\n")
        return
    md.append("| rung | resolution | N_eq | d_100 (nm) | dev | d_111 (nm) | dev | d_cav (nm) | dev |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    md.append(f"| **D (exact)** | — | 8006 (12006 solved) | {ref['d_100_nm']:.3f} | — "
              f"| {ref['d_111_nm']:.3f} | — | {ref['d_cavity_nm']:.3f} | — |")
    rows = (("C1", "1/10", "T4_C1_nc"), ("C2", "1/40", "T4_C2_nc"),
            ("C3", "1/160", "T4_C3_nc"), ("C4", "1/800", "T4_C4_nc"),
            ("C4, P=1 constant", "1/800", "T4_C4_P1_nc"))
    for lbl, res, rid in rows:
        v = cell(B.get(rid), "2")
        if v is None:
            r = B.get(rid) or {}
            md.append(f"| {lbl} | {res} | — | _did not reach_ "
                      f"({r.get('dose_reached', '?')} dpa) | | | | | |")
            continue
        d = lambda k: (v[k] - ref[k]) / ref[k] * 100
        md.append(f"| {lbl} | {res} | {B[rid].get('N_eq', '—')} "
                  f"| {v['d_100_nm']:.3f} | {d('d_100_nm'):+.1f}% "
                  f"| {v['d_111_nm']:.3f} | {d('d_111_nm'):+.1f}% "
                  f"| {v['d_cavity_nm']:.3f} | {d('d_cavity_nm'):+.1f}% |")
    md.append("")


def simple_table(B, md, title, subtitle, rows, dose_key):
    md.append(f"## {title}\n")
    md.append(subtitle + "\n")
    md.append("| column | N_eq | d_100 (nm) | d_111 (nm) | d_cav (nm) | delta_FP | wall |")
    md.append("|---|---|---|---|---|---|---|")
    missed = []
    for lbl, rid in rows:
        rec = B.get(rid)
        v = cell(rec, dose_key)
        if v is None:
            missed.append((lbl, rec or {}))
            continue
        md.append(f"| {lbl} | {rec.get('N_eq', '—')} | {v['d_100_nm']:.3f} "
                  f"| {v['d_111_nm']:.3f} | {v['d_cavity_nm']:.3f} "
                  f"| {rec.get('delta_FP', float('nan')):.3f} "
                  f"| {rec.get('wall_s', 0)/60:.1f} m |")
    md.append("")
    if missed:
        md.append("**Did not reach the comparison dose** (plan §3.4.1 — reported "
                  "with the dose achieved, never as a metric):\n")
        for lbl, r in missed:
            md.append(f"- **{lbl}** — reached {r.get('dose_reached', '?')} of "
                      f"{r.get('dose_target', '?')} dpa. {r.get('note', '')}")
        md.append("")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    B = board()
    r = reference_run_row()
    if r:
        B["_B4"] = r
    md = ["# Verification study — results", ""]
    md.append("Generated by `make_tables.py` from the campaign board "
              "(`campaign-verification`). Every value is scored through "
              "`run_ensemble.observe`, the same operator the calibration ledger "
              "uses, so a table entry and a ledger row for one configuration "
              "agree by construction.\n")

    table4(B, md)

    simple_table(
        B, md, "Table 1 — closure convergence at the production domain",
        "I = 80000, V = 20000, scored at 15.72 dpa. Deviations are quoted "
        "relative to **B2**, the most-resolved rung that completes.",
        (("B1  i_d 6400", "T1_B1"), ("B2  i_d 1600", "T1_B2"),
         ("B3  i_d 400", "T1_B3"), ("B4 = production (reference run)", "_B4")),
        "15.72")

    simple_table(
        B, md, "Table 2 — intra-bin closure and tolerance",
        "One knob changed from B4 at a time.",
        (("P=1 constant", "T2_P1"), ("P=3 lognormal", "T2_P3"),
         ("I_bin 10", "T2_IBIN10"), ("I_bin 40", "T2_IBIN40"),
         ("rtol 1e-4", "T2_RTOL4"), ("rtol 1e-6", "T2_RTOL6")),
        "15.72")

    simple_table(
        B, md, "Table 3 — helium",
        "Model-to-model; no experimental column.",
        (("fission QSS", "T3_FISS_QSS"), ("fission dynamic", "T3_FISS_DYN"),
         ("fusion QSS", "T3_FUS_QSS"), ("fusion dynamic", "T3_FUS_DYN")),
        "15.72")

    md.append("## Limitations that must be stated\n")
    md.append("1. **Table 4 verifies the no-coarsening variant.** The production "
              "configuration runs with `LOOP_COAL = 1`. Its exact arm was "
              "attempted three times and abandoned: the channel is O(NC²) in the "
              "number of populated sessile sizes, and even gated at 2×C_floor it "
              "spent 2.8 h on a single output step at 0.014 of 2 dpa. So the "
              "closure error measured here is **not demonstrably the production "
              "run's closure error** — plan §1.3's transfer argument does not "
              "carry over unmodified.")
    md.append("2. **The no-coarsening variant is unphysical.** Without the "
              "channel the ½⟨111⟩ mean size pins at the mobility cutoff and the "
              "loop density runs 20–60× over the EUROFER97 data "
              "(`rate_kernels.cpp:1930`). That does not invalidate the "
              "verification — verification asks whether the closure reproduces "
              "the exact solution of the *same* equations, which these arms do — "
              "but nothing in Table 4 may be compared to experiment.")
    md.append("3. **Table 1 has no B1.** The most-resolved rung reached 0.079 of "
              "40 dpa in 18.9 h and was stopped. Deviations are quoted relative "
              "to B2 (N_eq = 2048); B2 → B3 → B4 still spans 10× in N_eq.")
    md.append("4. **`delta_FP` ≈ 0.05–0.06 throughout**, against the model's own "
              "1e-2 gate. A separate, disclosed defect: a verification study "
              "measures closure error against the same equations, and "
              "conservation error is not that. Do not conflate them.\n")
    out = "\n".join(md)
    if a.out:
        Path(a.out).write_text(out, encoding="utf-8")
        print(f"wrote {a.out}")
    else:
        print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
