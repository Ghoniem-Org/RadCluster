#!/usr/bin/env python3
"""Score T9 on the EUROFER97 character CROSSOVER, per parameter vector.

WHY NOT score_T4.  score_T4 compares ONE row to ONE datum.  The character target
is not a value at a temperature -- it is a trend across two:

    300 C / 15 dpa    ->  f_100_tem_1 < 0.5   (1/2<111> dominant)
    350 C / 16.3 dpa  ->  f_100_tem_1 > 0.5   (<100> dominant)

Neither end alone is sufficient.  A theta that converts uniformly too little
passes the cold end and fails the hot one; a theta that converts uniformly too
much does the reverse.  Only the PAIR pins the crossover temperature, which is
the quantity E_a0_conv actually sets.  So rows are grouped by theta_id (both
conditions carry the same parameter vector) and scored jointly.

The 330 C datum is deliberately NOT used for character -- its database cell is a
crossover statement, not an assignment (see targets_T4.json
_character_correction).  It still pins size and density and is scored separately
when available.

REPORTED, NOT SCORED: whether d_111 clears the 1.0 nm TEM cutoff.  That is the
diagnostic the whole widened box exists to move -- if d_111 stays below 1 nm
everywhere, the 1/2<111> population is invisible by construction and no amount
of parameter search can make the character come out right, which is a structural
finding rather than a bad fit.

Usage:
    python score_crossover.py --rows t9.jsonl [--design design/T9_crossover.csv]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

COLD, HOT = "C_300_G7", "C_350_G7"
D_MIN_TEM = 1.0          # nm; the cutoff the model's 1/2<111> has to clear
# hot end: 350 C / 16.3 dpa, d = 10-15 nm, N = 0.8e22
HOT_D, HOT_D_RANGE, HOT_N = 12.5, (10.0, 15.0), 8.0e21
# cold end: 300 C / 15 dpa, d < 5 nm, N > 1.4e22 (reported as a lower bound)
COLD_D_MAX, COLD_N_MIN = 5.0, 1.4e22


def d_mean(r: dict):
    f, a, b = r.get("f_100_tem_1"), r.get("d_100_nm"), r.get("d_111_nm")
    if f is None or not a or not b:
        return None
    return f * a + (1.0 - f) * b


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=Path, required=True)
    ap.add_argument("--design", type=Path, default=HERE / "design/T9_crossover.csv")
    ap.add_argument("--top", type=int, default=15)
    a = ap.parse_args(argv)

    theta_of, param_of = {}, {}
    for d in csv.DictReader(a.design.open()):
        theta_of[int(d["row_id"])] = int(d["theta_id"])
        param_of[int(d["theta_id"])] = d

    rows = [json.loads(l) for l in a.rows.read_text().splitlines() if l.strip()]
    by_theta: dict[int, dict] = {}
    for r in rows:
        t = theta_of.get(r.get("row_id"))
        if t is None:
            continue
        by_theta.setdefault(t, {})[r.get("condition")] = r

    print(f"{len(rows)} rows -> {len(by_theta)} theta vectors "
          f"({sum(1 for v in by_theta.values() if len(v) == 2)} with BOTH ends)\n")

    # --- the structural question, asked of every row we have ----------------
    d111 = [r["d_111_nm"] for r in rows if isinstance(r.get("d_111_nm"), (int, float))]
    if d111:
        over = [x for x in d111 if x >= D_MIN_TEM]
        print(f"d_111 vs the {D_MIN_TEM} nm TEM cutoff over {len(d111)} rows: "
              f"min={min(d111):.2f}  median={sorted(d111)[len(d111)//2]:.2f}  "
              f"max={max(d111):.2f} nm")
        print(f"  rows clearing the cutoff: {len(over)}/{len(d111)} "
              f"({100*len(over)/len(d111):.0f} %)"
              + ("   <-- the widened box DOES move it" if over else
                 "   <-- STRUCTURAL: 1/2<111> invisible everywhere, character"
                 " cannot be fitted"))
        print()

    scored = []
    for t, ends in by_theta.items():
        c, h = ends.get(COLD), ends.get(HOT)
        if not c or not h:
            continue
        fc, fh = c.get("f_100_tem_1"), h.get("f_100_tem_1")
        if fc is None or fh is None:
            continue
        cold_ok, hot_ok = fc < 0.5, fh > 0.5
        # crossover margin: how decisively both ends land on the right side.
        # Negative if either end is wrong; the pair is what is being scored.
        margin = min(0.5 - fc, fh - 0.5)
        terms = []
        dh = d_mean(h)
        if dh and dh > 0:
            terms.append(abs(math.log10(dh / HOT_D)))
        nh = (h.get("N_100_vis_1") or 0) + (h.get("N_111_vis_1") or 0)
        if nh > 0:
            terms.append(abs(math.log10(nh / HOT_N)))
        dc = d_mean(c)
        if dc and dc > COLD_D_MAX:                 # only a bound: penalise excess
            terms.append(abs(math.log10(dc / COLD_D_MAX)))
        nc = (c.get("N_100_vis_1") or 0) + (c.get("N_111_vis_1") or 0)
        if 0 < nc < COLD_N_MIN:                    # lower bound, penalise shortfall
            terms.append(abs(math.log10(nc / COLD_N_MIN)))
        size_misfit = sum(terms) if terms else float("nan")
        penalty = 0.0 if (cold_ok and hot_ok) else 1.0
        scored.append((penalty + size_misfit, cold_ok, hot_ok, margin,
                       size_misfit, t, c, h))

    scored.sort(key=lambda x: (math.isnan(x[0]), x[0]))
    both = [s for s in scored if s[1] and s[2]]
    print(f"theta vectors reproducing the CROSSOVER (cold<0.5 AND hot>0.5): "
          f"{len(both)}/{len(scored)}\n")

    print(f"{'theta':>6} {'score':>7} {'cold f100':>10} {'hot f100':>9} "
          f"{'d_mean C':>9} {'d_mean H':>9} {'N_hot':>10} {'cross?':>7}")
    for s in scored[:a.top]:
        tot, ck, hk, mg, sm, t, c, h = s
        nh = (h.get("N_100_vis_1") or 0) + (h.get("N_111_vis_1") or 0)
        dc, dh = d_mean(c), d_mean(h)
        print(f"{t:6d} {tot:7.3f} {c['f_100_tem_1']:10.4f} {h['f_100_tem_1']:9.4f} "
              f"{(dc or float('nan')):9.2f} {(dh or float('nan')):9.2f} {nh:10.2e} "
              f"{'YES' if (ck and hk) else ('cold' if ck else ('hot' if hk else 'no')):>7}")

    if both:
        print("\nCROSSOVER-CORRECT vectors, parameter values "
              "(* = outside the original elicited prior):")
        meta = json.loads(a.design.with_suffix(".meta.json").read_text())
        samp = meta["sampled"]
        for s in both[:5]:
            t = s[5]
            d = param_of[t]
            bits = []
            for k, v in samp.items():
                val = float(d[k])
                op = v.get("original_prior")
                star = "*" if (op and not (op[0] <= val <= op[1])) else ""
                bits.append(f"{k}={val:.4g}{star}")
            print(f"  theta {t}: " + "  ".join(bits))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
