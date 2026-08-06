#!/usr/bin/env python
"""make_tables - the analysis tables for a campaign, joined into one place.

    python make_tables.py --design design/T2_design_v1.csv --results results/

Writes a single named artifact set into report/, default stem
`provisional_twin_map` (override with --name):

  provisional_twin_map.md             the report: what was run, model against
                                      experiment, the sigma_model floor, and an
                                      explicit statement of what is NOT included.

  provisional_twin_map_runs.csv       ONE ROW PER MODEL EVALUATION: identity,
                                      every theta component, every observable,
                                      the admissibility block, and provenance
                                      hashes.  The master table -- everything
                                      else is a view of it.  The .jsonl stores
                                      only theta_hash, so this join against the
                                      design file is what puts theta beside the
                                      result it produced.

  ..._runs_admissible.csv             the same, usable rows only.

  ..._observables.csv                 per observable per condition: model
                                      p05/median/p95 over admissible rows
                                      against the experimental band, with the
                                      in-band fraction.

  ..._targets.csv                     the experimental table, standalone.

"PROVISIONAL" is meant literally: see the three exclusions below.  This is a
map of where the model sits relative to the data under the current prior box,
not a calibrated twin.

WHAT THESE TABLES DO NOT CONTAIN, AND WHY
-----------------------------------------
1. RESULTS EXTRAPOLATED TO EXPERIMENTAL CONDITIONS.  Deliberately absent.  The
   Tier-2 ensemble runs at d_LF = 0.1 dpa; the database is at 15-30 dpa.  There
   is no defensible closed-form extrapolation across that gap -- <100> loops in
   this model have no size cap (loop->network is off and its P_ld is ~0 at these
   loop sizes anyway), so a run that matches at 0.1 dpa keeps growing and any
   naive scaling would be confidently wrong.  Bridging LF to experimental
   conditions is precisely the job of the multi-fidelity emulator (plan Tier 4,
   surrogate.py + calibrate.py), which is not built.  Emitting a column labelled
   "extrapolated" before that exists would be inventing data.

2. POSTERIOR UNCERTAINTY.  Also absent.  Tier 2 produces SENSITIVITY (Sobol
   S_i, S_i^T -- which parameters matter), not UNCERTAINTY (posterior spread on
   theta given the data).  What these tables can honestly show is the ENSEMBLE
   SPREAD -- p05/median/p95 of each observable over the prior box -- which is a
   prior-predictive spread, not a posterior.  They are different quantities and
   conflating them would overstate what the campaign knows.

3. sigma_model.  The plan expects model error to dominate sigma_exp for f_100.
   The convention spread IS computable here (f_100 is emitted at four TEM
   cutoffs) and is reported as a floor on sigma_model, but the full model-error
   term belongs to calibration.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

OBSERVABLES = ["N_loops_100", "d_100_nm", "d_100_content_nm", "mean_n_100",
               "N_loops_111", "d_111_nm", "mean_n_111",
               "f_100_number", "f_100_content",
               "f_100_tem_0p8", "f_100_tem_1", "f_100_tem_1p25", "f_100_tem_1p5",
               "N_100_visible", "N_111_visible", "N_voids", "mean_n_v",
               "d_cavity_nm", "S_inventory"]
ADMIN = ["admissible", "grid_limited", "starved", "dose_reached",
         "occ_100", "occ_111", "pile_100", "pile_111",
         "d_over_ceiling_100", "d_ceiling_100_nm", "delta_FP", "delta_He",
         "solver_rc", "wall_s", "error"]
PROV = ["machine_id", "git_sha", "solver_sha256", "workbook_sha256",
        "design_sha256", "run_cfg_sha", "theta_hash"]
TEM_COLS = ["f_100_tem_0p8", "f_100_tem_1", "f_100_tem_1p25", "f_100_tem_1p5"]


def read_design(design: Path):
    meta = json.loads(design.with_suffix(".meta.json").read_text(encoding="utf-8"))
    lines = design.read_text(encoding="utf-8").strip().splitlines()
    cols = lines[0].split(",")
    rows = {}
    for ln in lines[1:]:
        v = dict(zip(cols, ln.split(",")))
        rows[int(v["row_id"])] = v
    return rows, meta, cols


def load_results(results_dir: Path):
    recs = {}
    for f in sorted(Path(results_dir).glob("*.jsonl")):
        for ln in f.read_text(encoding="utf-8").splitlines():
            if not ln.strip():
                continue
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue
            recs[r["row_id"]] = r
    return recs


def _s(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        return "" if not np.isfinite(v) else f"{v:.10g}"
    return str(v).replace(",", ";").replace("\n", " ")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=HERE / "results")
    ap.add_argument("--targets", type=Path, default=HERE / "targets.json")
    ap.add_argument("--out", type=Path, default=HERE / "report")
    ap.add_argument("--name", default="provisional_twin_map",
                    help="stem for the artifact set (default provisional_twin_map)")
    a = ap.parse_args(argv)

    design, meta, dcols = read_design(a.design)
    recs = load_results(a.results, meta.get('design_sha256'))
    targets = json.loads(a.targets.read_text(encoding="utf-8"))
    theta_keys = meta["parameters"]
    a.out.mkdir(parents=True, exist_ok=True)

    # ── 1. master table: theta beside the result it produced ────────────────
    head = (["row_id", "condition", "matrix", "base_idx", "param_j"]
            + [f"theta.{k}" for k in theta_keys]
            + [f"obs.{k}" for k in OBSERVABLES]
            + [f"adm.{k}" for k in ADMIN]
            + [f"prov.{k}" for k in PROV])
    full, adm_rows = [], []
    for rid in sorted(design):
        d = design[rid]
        r = recs.get(rid)
        line = [rid, d["condition"], d["matrix"], d["base_idx"], d["param_j"]]
        line += [d.get(k, "") for k in theta_keys]
        line += [(r or {}).get(k) for k in OBSERVABLES]
        line += [(r or {}).get(k) for k in ADMIN]
        line += [(r or {}).get(k) for k in PROV]
        full.append(line)
        if r and not r.get("solver_rc") and r.get("admissible"):
            adm_rows.append(line)

    def write(path, rows):
        with path.open("w", encoding="utf-8", newline="") as fh:
            fh.write(",".join(head) + "\n")
            for ln in rows:
                fh.write(",".join(_s(v) for v in ln) + "\n")

    stem = a.name
    write(a.out / f"{stem}_runs.csv", full)
    write(a.out / f"{stem}_runs_admissible.csv", adm_rows)

    # ── 2. observables vs the experimental table, per condition ─────────────
    obs_rows = []
    conds = meta["conditions"]
    for cond in conds + ["ALL"]:
        sel = [r for rid, r in recs.items()
               if not r.get("solver_rc") and r.get("admissible")
               and (cond == "ALL" or design[rid]["condition"] == cond)]
        for key, t in targets.get("observables", {}).items():
            vals = np.array([r[key] for r in sel
                             if r.get(key) is not None
                             and np.isfinite(r.get(key, np.nan))], dtype=float)
            lo, hi, tgt = t.get("lo"), t.get("hi"), t.get("target")
            row = {"condition": cond, "observable": key, "units": t.get("units", ""),
                   "n_admissible": int(vals.size),
                   "model_p05": None, "model_median": None, "model_p95": None,
                   "exp_target": tgt, "exp_lo": lo, "exp_hi": hi,
                   "frac_in_band": None, "median_over_target": None}
            if vals.size:
                row["model_p05"] = float(np.percentile(vals, 5))
                row["model_median"] = float(np.median(vals))
                row["model_p95"] = float(np.percentile(vals, 95))
                if lo is not None and hi is not None:
                    row["frac_in_band"] = float(np.mean((vals >= lo) & (vals <= hi)))
                if tgt:
                    row["median_over_target"] = float(np.median(vals) / tgt)
            obs_rows.append(row)
    ohead = list(obs_rows[0].keys()) if obs_rows else []
    with (a.out / f"{stem}_observables.csv").open("w", encoding="utf-8",
                                                  newline="") as fh:
        fh.write(",".join(ohead) + "\n")
        for r in obs_rows:
            fh.write(",".join(_s(r[k]) for k in ohead) + "\n")

    # ── 3. the experimental table, standalone ───────────────────────────────
    with (a.out / f"{stem}_targets.csv").open("w", encoding="utf-8",
                                              newline="") as fh:
        fh.write("observable,target,lo,hi,units,note\n")
        for k, t in targets.get("observables", {}).items():
            fh.write(",".join([k, _s(t.get("target")), _s(t.get("lo")),
                               _s(t.get("hi")), _s(t.get("units")),
                               _s(t.get("note", ""))]) + "\n")

    # ── 4. sigma_model floor from the TEM-convention spread ─────────────────
    adm = [r for r in recs.values()
           if not r.get("solver_rc") and r.get("admissible")]
    spread = None
    if adm:
        per_run = []
        for r in adm:
            v = [r[c] for c in TEM_COLS
                 if r.get(c) is not None and np.isfinite(r.get(c, np.nan))]
            if len(v) >= 2:
                per_run.append(max(v) - min(v))
        if per_run:
            spread = {"n": len(per_run), "median": float(np.median(per_run)),
                      "p95": float(np.percentile(per_run, 95))}

    # ── 5. human-readable summary ───────────────────────────────────────────
    n_tot = len(design)
    n_done = len(recs)
    md = [f"# Provisional twin map — {a.design.name}", "",
          "*Provisional*: this maps where the model sits relative to the data "
          "under the current prior box. It is not a calibrated twin — see "
          "\"What is deliberately NOT here\" at the end.", "",
          f"- design `{meta.get('design_sha256','')[:16]}`, tier {meta.get('tier')}, "
          f"N={meta.get('N')}, p={meta.get('p')}, conditions {', '.join(conds)}",
          f"- rows {n_done}/{n_tot} evaluated, {len(adm)} admissible", ""]
    if meta.get("revision_pending"):
        md += [f"> **Design carries REVISION_PENDING parameters:** "
               f"{', '.join(meta['revision_pending'])}. Their ranges predate the "
               f"2026-08-01 absorption-gate work.", ""]
    md += ["## Files", "",
           "| file | contents |",
           "|---|---|",
           f"| `{stem}_runs.csv` | one row per evaluation: θ, observables, admissibility, provenance |",
           f"| `{stem}_runs_admissible.csv` | the same, usable rows only |",
           f"| `{stem}_observables.csv` | model p05/median/p95 vs the experimental band, per condition |",
           f"| `{stem}_targets.csv` | the experimental table |",
           "| `T2_sobol_indices.csv` | S_i, S_i^T (from `merge_and_sobol.py`) |", ""]
    md += ["## Model vs experiment (all conditions, admissible rows)", "",
           "| observable | n | p05 | median | p95 | target | in band |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for r in [x for x in obs_rows if x["condition"] == "ALL"]:
        if not r["n_admissible"]:
            continue
        fb = "n/a" if r["frac_in_band"] is None else f"{100*r['frac_in_band']:.0f}%"
        tg = "n/a" if r["exp_target"] is None else f"{r['exp_target']:.3g}"
        md.append(f"| `{r['observable']}` | {r['n_admissible']} | "
                  f"{r['model_p05']:.3g} | {r['model_median']:.3g} | "
                  f"{r['model_p95']:.3g} | {tg} | {fb} |")
    md += ["", "p05/median/p95 are the **prior-predictive** spread of the "
           "ensemble over the prior box. They are *not* a posterior — see below.", ""]
    if spread:
        md += ["## σ_model floor from the TEM-cutoff convention", "",
               f"Per-run spread of `f_100_tem` across d_min ∈ {{0.8, 1.0, 1.25, "
               f"1.5}} nm, over {spread['n']} admissible runs:", "",
               f"- median **{spread['median']:.3f}**, p95 **{spread['p95']:.3f}**",
               "", "This is a *floor* on σ_model for the loop-fraction "
               "observable: it is the variation produced by the reporting "
               "convention alone, with θ held fixed. If a calibration returns a "
               "posterior on the conversion parameters tighter than this, σ_model "
               "was set too small.", ""]
    md += ["## What is deliberately NOT here", "",
           "**Results extrapolated to experimental conditions.** The Tier-2 "
           "ensemble runs at 0.1 dpa; the database is at 15–30 dpa. There is no "
           "defensible closed-form extrapolation across that gap — ⟨100⟩ loops "
           "have no size cap in the present model, so a run matching at 0.1 dpa "
           "keeps growing and any naive scaling would be confidently wrong. "
           "Bridging LF to experimental conditions is the multi-fidelity "
           "emulator's job (plan Tier 4, `surrogate.py` + `calibrate.py`), which "
           "is not built. A column labelled *extrapolated* before that exists "
           "would be invented data.", "",
           "**Posterior uncertainty.** Tier 2 produces *sensitivity* (which "
           "parameters matter), not *uncertainty* (posterior spread on θ given "
           "the data). The spread reported above is prior-predictive. Posterior "
           "UQ requires Tier 4.", ""]
    (a.out / f"{stem}.md").write_text("\n".join(md), encoding="utf-8")

    print(f"  wrote {a.out}/")
    for f in (f"{stem}.md", f"{stem}_runs.csv", f"{stem}_runs_admissible.csv",
              f"{stem}_observables.csv", f"{stem}_targets.csv"):
        p = a.out / f
        print(f"    {f:32s} {p.stat().st_size:>9,d} bytes")
    print(f"\n  {n_done}/{n_tot} rows evaluated, {len(adm)} admissible")
    print("  NOTE: no extrapolation to experimental conditions and no posterior")
    print(f"  UQ - both require Tier 4. See {stem}.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
