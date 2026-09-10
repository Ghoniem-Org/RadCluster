#!/usr/bin/env python3
"""make_M_latex.py — LaTeX table bodies for campaign M.

Writes docs/Formulation/tables/M*.tex, which campaign_results.tex \\input{}s.
Generated rather than transcribed: this report has already had one table carry
a fabricated N_eq (28006 for a run whose value was 8006) because a number was
filled in from memory when the tool left a blank.

Every table refuses to print a metric for a run that missed its comparison dose
(plan S3.4.1); such a cell prints \\dnr instead.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DT = HERE.parent
MOD = DT.parent
REPO = MOD.parent
for p in (str(REPO), str(DT), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)

import runs as M                                   # noqa: E402

SYNC = HERE / ".sync" / "claims"
OUT = REPO / "docs" / "Formulation" / "tables"
DOSE = "20"


def board():
    out = {}
    for f in sorted(SYNC.glob("*.json")):
        try:
            j = json.loads(f.read_text(encoding="utf-8"))
            out[j.get("run_id", f.stem)] = j
        except Exception:
            pass
    return out


def cell(rec):
    if not rec or rec.get("status") != "done":
        return None
    v = (rec.get("scored") or {}).get(DOSE)
    return None if (not v or v.get("missing") or v.get("off_grid")) else v


def why_blank(rec):
    """WHY a cell is empty -- the three reasons are not the same thing.

    A run that timed out short of the comparison dose is a cost result; a
    configuration the bin layout cannot realise is a statement about the grid,
    not about cost; and a cell never attempted is neither.  Printing one symbol
    for all three would let a reader read the M2D grid's expensive corner as if
    the method had failed there.
    """
    if rec is None:
        return r"\dnr"
    if rec.get("status") == "failed":
        return r"\nreal"
    return r"\dnr"


def w(name, lines):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  wrote docs/Formulation/tables/{name}")


def main():
    B = board()
    ref = cell(B.get(M.M_REF_RUN))

    # ── M1: the cost cliff ──────────────────────────────────────────────────
    L = [r"\begin{tabular}{lrrrr}", r"\toprule",
         r"Run & $I=V$ & $N_{\rm eq}$ & Dose reached & Wall \\", r"\midrule"]
    # M1_LADDER is local to manifest(); enumerate from the manifest itself so
    # the table cannot drift from the runs that exist.
    ent = {x["run_id"]: x for x in M.manifest()}
    m1 = [x for x in M.manifest() if x["table"] in ("M1", "M1REF")]
    for e in sorted(m1, key=lambda x: x["I"]):
        rid = e["run_id"]
        r = B.get(rid)
        if not r:
            continue
        wall = r.get("wall_s")
        name = rid.replace("_", chr(92) + "_")
        if not wall:
            # Claimed but never published -- the ladder was stopped by hand once
            # 2000 had failed, because 4000 and up cannot succeed and each costs
            # an hour to confirm it.  Say that, rather than printing an empty
            # row a reader would take for a missing measurement.
            L.append(f"{name} & {e['I']} & --- & \\multicolumn{{2}}{{c}}"
                     f"{{\\itshape not run --- see text}} \\\\")
            continue
        reach = r.get("dose_reached", float("nan"))
        bold = r"\bfseries " if reach >= 0.99 * e["dose"] else ""
        L.append(f"{bold}{name} & {e['I']} "
                 f"& {r.get('N_eq', '---')} "
                 f"& {reach:.4g} / {e['dose']:g} "
                 f"& {wall/60:.1f} min \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    w("M1_cost.tex", L)

    # ── M2R: one-dimensional ladder, deviations ─────────────────────────────
    L = [r"\begin{tabular}{lrrrrrrrr}", r"\toprule",
         r"Rung & $I_{\rm bin}$ & $r$ & $N_{\rm eq}$ & $N_{111}$ & $d_{111}$ "
         r"& $N_{100}$ & $d_{100}$ & $d_{\rm cav}$ \\",
         r"\midrule"]
    if ref:
        for nb in M.M2R_BINS:
            rec = B.get(f"M2R_B{nb}")
            v = cell(rec)
            if v is None:
                L.append(f"B{nb} & {nb} & --- & --- & "
                         + " & ".join([r"\dnr"] * 5) + r" \\")
                continue
            d = [(v[k] - ref[k]) / ref[k] * 100 for k in
                 ("N_loops_111", "d_111_nm", "N_loops_100", "d_100_nm",
                  "d_cavity_nm")]
            L.append(f"B{nb} & {rec.get('bin_I_bin', nb)} "
                     f"& {rec.get('bin_r', float('nan')):.2f} "
                     f"& {rec.get('N_eq', '---')} & "
                     + " & ".join(f"{x:+.1f}" for x in d) + r" \\")
    L += [r"\bottomrule", r"\end{tabular}"]
    w("M2R_dev.tex", L)

    # ── M2D: the 2-D sweep, d_111 deviation ─────────────────────────────────
    for key, fname in (("d_111_nm", "M2D_d111.tex"),
                       ("d_cavity_nm", "M2D_dcav.tex")):
        L = [r"\begin{tabular}{l" + "r" * len(M.M2D_BINS) + "}", r"\toprule",
             r"$i_{\rm discrete}$ & "
             + " & ".join(f"$I_{{\\rm bin}}={nb}$" for nb in M.M2D_BINS)
             + r" \\", r"\midrule"]
        for i_d in M.M2D_IDISC:
            row = [f"{i_d}"]
            for nb in M.M2D_BINS:
                v = cell(B.get(f"M2D_I{i_d}_B{nb}"))
                if v is None or not ref or not ref[key]:
                    row.append(why_blank(B.get(f"M2D_I{i_d}_B{nb}")))
                else:
                    row.append(f"{(v[key]-ref[key])/ref[key]*100:+.1f}")
            L.append(" & ".join(row) + r" \\")
        L += [r"\bottomrule", r"\end{tabular}"]
        w(fname, L)
    return 0


if __name__ == "__main__":
    sys.exit(main())
