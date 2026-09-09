#!/usr/bin/env python3
"""make_M_report.py — campaign M: binning convergence against the exact arm.

    python make_M_report.py            # markdown to stdout
    python make_M_report.py --out FILE

Campaign M asks one question that Table 4 could not answer cleanly: at a FIXED
domain and a FIXED discrete core (i_discrete = v_discrete = 50), how do the six
observables converge as the bin count rises?  Table 4 varied the discrete core
and the bin count together, so its rungs moved along two axes at once.  Here
only the binning moves, and the reference is a fully discrete solve of the same
domain -- so a deviation is closure error and nothing else.

The reference is M1_D10000: I = V = 10000 discrete, one core, no wall-clock
cap.  Until it finishes, the rungs are still reported -- their mutual spread IS
self-convergence evidence -- but the deviation columns are withheld rather than
computed against the finest rung.  Quoting convergence against the finest
AVAILABLE approximation is how a ladder appears to converge to the wrong
answer; plan S3.4.1 exists for the same reason.
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

import runs as manifest_mod          # noqa: E402

SYNC = HERE / ".sync" / "claims"
SCORE_DOSE = "20"                    # every M row shares this grid (plan S3.2)

# The six observables, in the order the paper reports them.
OBS = [("N_loops_111", "N_111 (m^-3)", "{:.3e}"),
       ("d_111_nm",    "d_111 (nm)",   "{:.3f}"),
       ("N_loops_100", "N_100 (m^-3)", "{:.3e}"),
       ("d_100_nm",    "d_100 (nm)",   "{:.3f}"),
       ("N_voids",     "N_cav (m^-3)", "{:.3e}"),
       ("d_cavity_nm", "d_cav (nm)",   "{:.3f}")]


def board() -> dict:
    out = {}
    for f in sorted(SYNC.glob("*.json")):
        try:
            j = json.loads(f.read_text(encoding="utf-8"))
            out[j.get("run_id", f.stem)] = j
        except Exception:
            pass
    return out


def cell(rec):
    """The scored block at the comparison dose, or None if it may not be used."""
    if not rec or rec.get("status") != "done":
        return None
    v = (rec.get("scored") or {}).get(SCORE_DOSE)
    if not v or v.get("missing") or v.get("off_grid"):
        return None
    return v


def ladder(B, md, tbl, prefix, bins, dom, refid):
    """One binning ladder, with deviations only if its exact arm completed."""
    ref_rec = B.get(refid) if refid else None
    ref = cell(ref_rec)
    md.append(f"## {tbl} — binning convergence at I = V = {dom}\n")
    md.append("`i_discrete = v_discrete = 50` throughout; only the bin count "
              "moves, so a deviation below is closure error and nothing else.\n")
    if ref is None:
        why = (f"`{refid}` did not complete" if refid else
               "the discrete solve at this domain diverges")
        md.append(f"> **No exact arm: {why}.**  The rungs are listed so their "
                  "mutual spread can be read, but the deviation columns are "
                  "WITHHELD.  Scoring a ladder against its own finest rung is "
                  "how one converges confidently to the wrong answer.\n")

    md.append("| rung | I_bin=V_bin | r | N_eq | "
              + " | ".join(o[1] for o in OBS) + " |")
    md.append("|" + "---|" * (4 + len(OBS)))
    if ref is not None:
        md.append(f"| **{refid} (exact)** | — | — | {ref_rec.get('N_eq', '—')} | "
                  + " | ".join(f[2].format(ref[f[0]]) for f in OBS) + " |")
    rows = []
    for nb in bins:
        rid = f"{prefix}{nb}"
        rec = B.get(rid)
        v = cell(rec)
        if v is None:
            r = rec or {}
            got = r.get("dose_reached")
            why = (f"reached {got:.4g} dpa" if isinstance(got, (int, float))
                   else r.get("status", "not run"))
            md.append(f"| {rid} | {nb} | — | — | _{why}_ | "
                      + " | ".join(["—"] * (len(OBS) - 1)) + " |")
            continue
        md.append(f"| {rid} | {rec.get('bin_I_bin', nb)} "
                  f"| {rec.get('bin_r', float('nan')):.3f} "
                  f"| {rec.get('N_eq', '—')} | "
                  + " | ".join(f[2].format(v[f[0]]) for f in OBS) + " |")
        rows.append((nb, rec, v))
    md.append("")

    if ref is not None and rows:
        md.append(f"### {tbl}: deviation from the exact arm (%)\n")
        md.append("| rung | N_eq | r | "
                  + " | ".join(o[1].split(" ")[0] for o in OBS) + " | max |")
        md.append("|" + "---|" * (4 + len(OBS)))
        for nb, rec, v in rows:
            devs = [((v[k] - ref[k]) / ref[k] * 100.0) if ref[k] else float("nan")
                    for k, _, _ in OBS]
            worst = max(abs(d) for d in devs)
            md.append(f"| {prefix}{nb} | {rec.get('N_eq', '—')} "
                      f"| {rec.get('bin_r', float('nan')):.2f} | "
                      + " | ".join(f"{d:+.1f}" for d in devs)
                      + f" | **{worst:.1f}** |")
        md.append("")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    B = board()
    ent = {e["run_id"]: e for e in manifest_mod.manifest()}

    ref_id = f"M1_D{manifest_mod.M_REF_DOMAIN}"
    ref_rec = B.get(ref_id)
    ref = cell(ref_rec)

    md = ["# Campaign M — system size and binning", ""]
    md.append(f"Production mobility (`i_mobile = v_mobile = 5`), fission "
              f"cascade, `LOOP_COAL` and `prec_bw` at their defaults in BOTH "
              f"arms.  All rows share one 45-point log grid to "
              f"{SCORE_DOSE} dpa, so the comparison dose is the same point in "
              f"every column (plan S3.2).\n")

    # ── M1: the cost curve ──────────────────────────────────────────────────
    md.append("## M1 — cost against domain (fully discrete)\n")
    md.append("| rung | I = V | N_eq | dose reached | wall |")
    md.append("|---|---|---|---|---|")
    for rid, e in ent.items():
        if e["table"] not in ("M1", "M1REF"):
            continue
        r = B.get(rid)
        if not r:
            md.append(f"| {rid} | {e['I']} | — | _not run_ | — |")
            continue
        st = r.get("status", "?")
        w = r.get("wall_s")
        md.append(f"| {rid} | {e['I']} | {r.get('N_eq', '—')} "
                  f"| {r.get('dose_reached', float('nan')):.4g} / {e['dose']:g} dpa "
                  f"| {(w/60 if w else float('nan')):.1f} min ({st}) |")
    md.append("")

    # ── The binning ladders ─────────────────────────────────────────────────
    # Two of them: M2 at I = 10000 (no exact arm -- the discrete solve there
    # diverges) and M2R at I = 1000 (M1_D1000 completes, so it can be scored).
    for tbl, prefix, bins, dom, refid in (
            ("M2R", "M2R_B", manifest_mod.M2R_BINS, manifest_mod.M2R_DOMAIN,
             manifest_mod.M_REF_RUN),
            ("M2", "M2_B", manifest_mod.M2_BINS, manifest_mod.M_BIN_DOMAIN,
             None)):
        ladder(B, md, tbl, prefix, bins, dom, refid)

    md.append("## What must be stated\n")
    md.append(f"1. **The domain is truncated.** V = {manifest_mod.M_BIN_DOMAIN} "
              "sits below the mean cavity size at production dose, so d_cav "
              "reads the ceiling at the coarse rungs rather than the model. "
              "Both arms see the identical truncation, so it cancels out of the "
              "closure error — but nothing here may be compared to experiment.")
    md.append("2. **`delta_FP` runs ~0.05 throughout**, against the model's own "
              "1e-2 gate.  A separate, disclosed defect: closure error and "
              "conservation error are different quantities.  Do not conflate "
              "them.")
    md.append("3. **The reference is a single trajectory, not a converged "
              "domain.**  M1 shows how cost grows with I; it does not show that "
              "I = V = 10000 is itself converged in domain.  A rung agreeing "
              "with it agrees with the exact solution ON THIS DOMAIN.\n")

    out = "\n".join(md)
    if a.out:
        Path(a.out).write_text(out, encoding="utf-8")
        print(f"wrote {a.out}")
    else:
        print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
