#!/usr/bin/env python3
"""Score T4 rows against the experimental bands in targets_T4.json.

WHY THIS IS NOT A ONE-LINER.  The database reports what TEM measures: ONE loop
density and ONE mean diameter, over BOTH Burgers characters, above a detection
limit.  The model reports the two populations separately and unfiltered.  So the
comparison has to be built, not read off:

    N_loop_total  <-  N_100_vis_1 + N_111_vis_1      (TEM-visible, d_min = 1 nm)
    d_loop_mean   <-  number-weighted mean of d_100_nm and d_111_nm using those
                      SAME visible densities as weights
    dominant BV   <-  f_100_tem_1 > 0.5  vs  < 0.5

Using the raw N_loops_100 / N_loops_111 instead of the _vis_1 fields would
compare an unfiltered model against a filtered measurement, and the model's
1/2<111> population sits at ~0.8 nm -- entirely below the detection limit -- so
that error is worth a factor of several, not a few percent.

The diameter weighting is an APPROXIMATION: it weights each character's mean by
its visible count, which is right only if the within-character distributions do
not overlap the cutoff too strongly.  It is the best available from the row
fields; the exact quantity would need the size distribution.

Usage:
    python score_T4.py --rows <file.jsonl> [--targets targets_T4.json] [--top N]
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def visible(rec: dict) -> tuple[float, float]:
    """TEM-visible <100> and 1/2<111> number densities at the 1.0 nm cutoff."""
    n100 = rec.get("N_100_vis_1")
    n111 = rec.get("N_111_vis_1")
    if n100 is None:
        n100 = rec.get("N_100_visible", rec.get("N_loops_100", 0.0)) or 0.0
    if n111 is None:
        n111 = rec.get("N_111_visible", rec.get("N_loops_111", 0.0)) or 0.0
    return float(n100), float(n111)


def derived(rec: dict) -> dict:
    n100, n111 = visible(rec)
    tot = n100 + n111
    d100 = rec.get("d_100_nm") or 0.0
    d111 = rec.get("d_111_nm") or 0.0
    dmean = ((n100 * d100 + n111 * d111) / tot) if tot > 0 else float("nan")
    return {
        "N_loop_total_m3": tot,
        "d_loop_mean_nm": dmean,
        "f_100_tem_1": rec.get("f_100_tem_1"),
        "N_cavity_m3": rec.get("N_voids"),
        "d_cavity_nm": rec.get("d_cavity_nm"),
    }


def misfit(model: dict, target: dict) -> tuple[float, list[str]]:
    """Sum of |log10(model/target)| over the fields the target actually pins.

    Log ratio, not relative error: densities span decades, and a factor-3 miss
    should cost the same whether it is high or low.  Diameters use the same
    measure so one channel cannot dominate purely by units.
    """
    terms, notes = [], []

    # DENSITY IS DIAGNOSTIC ONLY, NOT SCORED (2026-08-12).  The model carries
    # two independent <100> number densities that disagree: N_100_vis_1 (from
    # the per-size y_sia100 array) exceeds N_loops_100 (from post_process, same
    # m^-3 units after post_process.py:471) in 95.3 % of 1173 T3 rows, median
    # 626x, max 65540x.  A detection-limit filter can only REMOVE loops, so one
    # of the two is wrong by ~3 decades and it is not yet established which.
    #
    # Character and mean size are UNAFFECTED, and that is not a hope -- it is
    # exact.  Both are built from the same two per-size sums a0 (<100>) and a1
    # (1/2<111>), so the 1/Omega normalisation cancels identically:
    #     f_100  = a0 / (a0 + a1)
    #     d_mean = (a0 d100 + a1 d111) / (a0 + a1) = f_100 d100 + (1-f_100) d111
    # Whatever scales a0 and a1 together cannot move either number.  So the
    # calibration scores on these two and reports density beside them.
    for key, tkey in (("d_loop_mean_nm", "d_loop_mean_nm"),):
        t = target.get(tkey)
        m = model.get(key)
        if t is None or m is None or not (m > 0) or not (t > 0):
            continue
        r = math.log10(m / t)
        terms.append(abs(r))
        notes.append(f"{tkey}: model {m:.4g} vs {t:.4g}  ({10**r:.2f}x)")

    tN, mN = target.get("N_loop_total_m3"), model.get("N_loop_total_m3")
    if tN and mN and mN > 0:
        notes.append(f"[diagnostic, NOT scored] N_loop_total_m3: model {mN:.4g} "
                     f"vs {tN:.4g}  ({mN/tN:.2f}x)")
    # character is a direction, not a number: score it as a hit/miss penalty
    exp_f = target.get("f_100_tem_1_expected")
    f = model.get("f_100_tem_1")
    if exp_f and f is not None:
        want_hi = exp_f.strip().startswith(">")
        ok = (f > 0.5) if want_hi else (f < 0.5)
        if not ok:
            terms.append(0.5)
        notes.append(f"character: f_100_tem_1 {f:.3f}, data says {exp_f} "
                     f"-> {'OK' if ok else 'INVERTED'}")
    return (sum(terms) if terms else float("nan")), notes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=Path, required=True)
    ap.add_argument("--targets", type=Path, default=HERE / "targets_T4.json")
    ap.add_argument("--top", type=int, default=10)
    a = ap.parse_args(argv)

    tg = json.loads(a.targets.read_text(encoding="utf-8"))
    by_id = {}
    for fam in ("neutron", "ion"):
        for p in tg[fam]["points"]:
            by_id[p["id"]] = p
    conds = json.loads((HERE / "conditions_T4.json").read_text(encoding="utf-8"))

    rows = [json.loads(l) for l in a.rows.read_text(encoding="utf-8").splitlines()
            if l.strip()]
    if not rows:
        print("no rows yet")
        return 0

    scored = []
    for r in rows:
        c = conds.get(r.get("condition"), {})
        t = by_id.get(c.get("_target_id"))
        if t is None:
            continue
        m = derived(r)
        s, notes = misfit(m, t)
        scored.append((s, r, m, t, notes))
    scored.sort(key=lambda x: (math.isnan(x[0]), x[0]))

    print(f"scored {len(scored)} rows against targets_T4.json\n")
    for s, r, m, t, notes in scored[:a.top]:
        chan = (f"chi={r.get('loop_net_chi')} w_c={r.get('loop_net_w_c')} "
                f"K_rec={r.get('loop_net_K_rec')}")
        print(f"row {r['row_id']:>4}  {r.get('condition')}  misfit={s:.3f}   {chan}")
        print(f"     dose_reached={r.get('dose_reached')}  starved={r.get('starved')}  "
              f"topbin_100={r.get('topbin_100')}  pile_100={r.get('pile_100')}")
        for n in notes:
            print(f"     {n}")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
