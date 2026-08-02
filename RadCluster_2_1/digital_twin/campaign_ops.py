#!/usr/bin/env python
"""campaign_ops - build, monitor and stop control for the digital-twin campaign.

Everything the control notebook drives lives here so the notebook cells stay
thin and this logic stays testable from a shell.

Three groups:

  BUILD        ensure_solver()    - verify or compile the C++ solver for THIS
                                    machine.  The binary is not in git, so each
                                    machine must build its own.

  STOP CONTROL request_stop()     - graceful halt.  Writes a sentinel that
                                    run_ensemble checks between rows: it stops
                                    submitting new work, lets in-flight rows
                                    finish and be written, then exits cleanly.
                                    Nothing is lost and the campaign resumes
                                    from where it stopped.

  STATUS       campaign_status()  - coverage, per-machine progress, timing and
                                    remaining-compute estimate, observables
                                    against the experimental bands, and the
                                    reasons rows were rejected.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
MODULE = REPO / "RadCluster_2_1"
STOP_FILE = HERE / "CAMPAIGN_STOP"


# ══════════════════════════════════════════════════════════ build ════════════
def find_solver() -> Path | None:
    for p in (MODULE / "build" / "Release" / "solver.exe",
              MODULE / "build" / "solver.exe",
              MODULE / "build" / "solver"):
        if p.exists():
            return p
    return None


def sha256_file(p) -> str:
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except Exception:
        return "unavailable"


def _newest_source_mtime() -> float:
    newest = 0.0
    for pat in ("**/*.cpp", "**/*.h", "**/CMakeLists.txt"):
        for f in (MODULE / "cpp_utils").glob(pat):
            newest = max(newest, f.stat().st_mtime)
    cml = MODULE / "CMakeLists.txt"
    if cml.exists():
        newest = max(newest, cml.stat().st_mtime)
    return newest


def solver_info() -> dict:
    p = find_solver()
    src = _newest_source_mtime()
    info = {"path": str(p) if p else None, "exists": p is not None,
            "sha256": sha256_file(p)[:16] if p else None,
            "mtime": p.stat().st_mtime if p else None,
            "newest_source_mtime": src}
    info["stale"] = bool(p and src > p.stat().st_mtime)
    return info


def ensure_solver(force: bool = False, verbose: bool = True) -> dict:
    """Verify the solver, and build it if missing, stale, or force=True.

    'Stale' means a .cpp/.h/CMakeLists is newer than the binary -- the exact
    situation that silently produces results from code you think you replaced.
    """
    info = solver_info()
    need = force or (not info["exists"]) or info["stale"]
    if verbose:
        if not info["exists"]:
            print("  solver: NOT FOUND -> building")
        elif info["stale"]:
            print("  solver: STALE (C++ source is newer than the binary) -> rebuilding")
        elif force:
            print("  solver: rebuild forced")
        else:
            print(f"  solver: OK  {info['path']}  sha {info['sha256']}")
    if not need:
        info["built"] = False
        return info

    cfg = subprocess.run(["cmake", "-S", str(MODULE), "-B", str(MODULE / "build")],
                         capture_output=True, text=True)
    if cfg.returncode != 0 and not (MODULE / "build").exists():
        info["built"] = False
        info["error"] = (cfg.stdout + cfg.stderr)[-1500:]
        if verbose:
            print("  cmake configure FAILED:\n" + info["error"])
        return info
    bld = subprocess.run(["cmake", "--build", str(MODULE / "build"),
                          "--config", "Release"], capture_output=True, text=True)
    out = bld.stdout + bld.stderr
    if bld.returncode != 0:
        info["built"] = False
        info["error"] = out[-2500:]
        if verbose:
            print("  build FAILED:\n" + info["error"])
            if "LNK1104" in out or "Permission denied" in out or "Text file busy" in out:
                print("\n  The linker could not write solver.exe. A campaign run is")
                print("  probably still holding it -- stop the workers first")
                print("  (campaign_ops.request_stop('rebuild')), then retry.")
        return info
    info = solver_info()
    info["built"] = True
    if verbose:
        print(f"  solver: BUILT  {info['path']}  sha {info['sha256']}")
    return info


def environment() -> dict:
    try:
        gs = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(REPO),
                                     stderr=subprocess.DEVNULL, text=True).strip()
        br = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"],
                                     cwd=str(REPO), stderr=subprocess.DEVNULL,
                                     text=True).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"],
                                             cwd=str(REPO), text=True).strip())
    except Exception:
        gs, br, dirty = "unknown", "unknown", False
    return {"machine_id": platform.node(), "platform": platform.platform(),
            "python": platform.python_version(), "cpu_count": os.cpu_count(),
            "git_sha": gs, "branch": br, "worktree_dirty": dirty,
            "workbook_sha256": sha256_file(MODULE / "input" /
                                           "input_parameters.xlsx")[:16]}


# ═══════════════════════════════════════════════════ stop control ════════════
def request_stop(reason: str = "user requested") -> Path:
    """Ask every running worker on THIS machine to stop gracefully.

    Workers finish the rows already in flight, write them, and exit 0.  No row
    is lost and no partial row is written, so a later restart resumes cleanly.
    """
    STOP_FILE.write_text(json.dumps(
        {"reason": reason, "requested_at": time.time(),
         "by": platform.node()}, indent=2), encoding="utf-8")
    print(f"  STOP requested: {reason}")
    print(f"  wrote {STOP_FILE}")
    print("  Workers finish their in-flight rows, then exit. Clear it with")
    print("  campaign_ops.clear_stop() before restarting.")
    return STOP_FILE


def clear_stop() -> None:
    if STOP_FILE.exists():
        STOP_FILE.unlink()
        print("  STOP cleared - the campaign may be restarted.")
    else:
        print("  no STOP flag set.")


def stop_requested() -> dict | None:
    if not STOP_FILE.exists():
        return None
    try:
        return json.loads(STOP_FILE.read_text(encoding="utf-8"))
    except Exception:
        return {"reason": "unreadable STOP file"}


# ══════════════════════════════════════════════════════════ status ═══════════
def load_targets(path: Path | None = None) -> dict:
    p = path or (HERE / "targets.json")
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def read_design(design: Path) -> tuple[list[dict], dict]:
    meta_p = Path(design).with_suffix(".meta.json")
    meta = json.loads(meta_p.read_text(encoding="utf-8")) if meta_p.exists() else {}
    lines = Path(design).read_text(encoding="utf-8").strip().splitlines()
    cols = lines[0].split(",")
    rows = []
    for ln in lines[1:]:
        v = dict(zip(cols, ln.split(",")))
        rows.append({"row_id": int(v["row_id"]), "condition": v["condition"],
                     "matrix": v["matrix"], "base_idx": int(v["base_idx"]),
                     "param_j": int(v["param_j"])})
    return rows, meta


def load_results(results_dir: Path) -> dict[int, dict]:
    recs: dict[int, dict] = {}
    for f in sorted(Path(results_dir).glob("*.jsonl")):
        for ln in f.read_text(encoding="utf-8").splitlines():
            if not ln.strip():
                continue
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue          # truncated last line after a hard kill
            recs[r["row_id"]] = r
    return recs


def campaign_status(design: Path, results_dir: Path, n_machines: int = 4,
                    workers_per_machine: int | None = None) -> dict:
    design_rows, meta = read_design(design)
    recs = load_results(results_dir)
    total = len(design_rows)
    done = [r for r in recs.values()]
    failed = [r for r in done if r.get("solver_rc")]
    admis = [r for r in done if not r.get("solver_rc") and r.get("admissible")]
    inadm = [r for r in done if not r.get("solver_rc") and not r.get("admissible")]

    walls = [r["wall_s"] for r in done if r.get("wall_s")]
    mean_w = float(np.mean(walls)) if walls else None
    med_w = float(np.median(walls)) if walls else None
    p90_w = float(np.percentile(walls, 90)) if walls else None

    per_machine = {}
    for k in range(n_machines):
        mine = [d for d in design_rows if d["row_id"] % n_machines == k]
        got = [d for d in mine if d["row_id"] in recs]
        ids = {recs[d["row_id"]].get("machine_id") for d in got}
        per_machine[k] = {"assigned": len(mine), "done": len(got),
                          "remaining": len(mine) - len(got),
                          "machine_id": sorted(x for x in ids if x) or None}

    W = workers_per_machine or max(1, (os.cpu_count() or 4) - 2)
    eta = {}
    if mean_w:
        for k, v in per_machine.items():
            eta[k] = v["remaining"] * mean_w / W
    core_h_used = sum(walls) / 3600.0 if walls else 0.0
    core_h_left = (sum(v["remaining"] for v in per_machine.values()) * mean_w
                   / 3600.0) if mean_w else None

    reasons = {}
    for r in inadm:
        if r.get("starved"):
            reasons["dose-starved"] = reasons.get("dose-starved", 0) + 1
        if r.get("grid_limited"):
            reasons["grid-limited"] = reasons.get("grid-limited", 0) + 1
        if abs(r.get("delta_FP") or 0) >= 1e-3:
            reasons["delta_FP>=1e-3"] = reasons.get("delta_FP>=1e-3", 0) + 1

    prov = {}
    for f in ("git_sha", "solver_sha256", "workbook_sha256", "design_sha256"):
        vals = {}
        for r in done:
            vals.setdefault(r.get(f, "missing"), []).append(r.get("machine_id", "?"))
        prov[f] = {k: sorted(set(v)) for k, v in vals.items()}

    return {"design": str(design), "meta": meta, "total": total,
            "done": len(done), "admissible": len(admis),
            "inadmissible": len(inadm), "failed": len(failed),
            "missing": total - len(done),
            "pct": 100.0 * len(done) / total if total else 0.0,
            "wall_mean_s": mean_w, "wall_median_s": med_w, "wall_p90_s": p90_w,
            "core_hours_used": core_h_used, "core_hours_remaining": core_h_left,
            "per_machine": per_machine, "eta_s": eta,
            "inadmissible_reasons": reasons, "provenance": prov,
            "workers_assumed": W,
            "stop": stop_requested(),
            "_admis": admis}


def _fmt_hms(s):
    if s is None:
        return "n/a"
    s = int(s)
    return f"{s//3600:d}h{(s%3600)//60:02d}m" if s >= 3600 else f"{s//60:d}m{s%60:02d}s"


def observable_summary(admis: list[dict], targets: dict) -> list[dict]:
    """Where the ensemble sits relative to each experimental band."""
    out = []
    for key, t in targets.get("observables", {}).items():
        vals = np.array([r[key] for r in admis
                         if r.get(key) is not None and np.isfinite(r.get(key, np.nan))],
                        dtype=float)
        if vals.size == 0:
            out.append({"observable": key, "n": 0, **t})
            continue
        lo, hi = t.get("lo"), t.get("hi")
        frac_in = (float(np.mean((vals >= lo) & (vals <= hi)))
                   if lo is not None and hi is not None else None)
        out.append({"observable": key, "n": int(vals.size),
                    "p05": float(np.percentile(vals, 5)),
                    "median": float(np.median(vals)),
                    "p95": float(np.percentile(vals, 95)),
                    "target": t.get("target"), "lo": lo, "hi": hi,
                    "units": t.get("units", ""), "frac_in_band": frac_in,
                    "brackets": (bool(vals.min() <= (t.get("target") or np.inf)
                                      <= vals.max())
                                 if t.get("target") is not None else None)})
    return out


def render_status(st: dict, targets: dict | None = None) -> None:
    """Human-readable campaign report."""
    m = st["meta"]
    print("=" * 78)
    print(f"CAMPAIGN  {Path(st['design']).name}   "
          f"p={m.get('p','?')} N={m.get('N','?')} "
          f"conditions={','.join(m.get('conditions', []))}")
    print("=" * 78)
    if st["stop"]:
        print(f"  *** STOP FLAG SET: {st['stop'].get('reason')} "
              f"(by {st['stop'].get('by')}) - workers are winding down ***\n")

    bar_n = 46
    filled = int(bar_n * st["pct"] / 100)
    print(f"  progress  [{'#'*filled}{'.'*(bar_n-filled)}] {st['pct']:5.1f}%  "
          f"{st['done']}/{st['total']} rows")
    print(f"            admissible {st['admissible']}   "
          f"inadmissible {st['inadmissible']}   failed {st['failed']}   "
          f"missing {st['missing']}")
    if st["inadmissible_reasons"]:
        print(f"            rejected because: {st['inadmissible_reasons']}")

    print(f"\n  timing    per row  mean {_fmt_hms(st['wall_mean_s'])}   "
          f"median {_fmt_hms(st['wall_median_s'])}   "
          f"p90 {_fmt_hms(st['wall_p90_s'])}")
    print(f"            core-hours used {st['core_hours_used']:.1f}"
          + (f"   remaining ~{st['core_hours_remaining']:.1f}"
             if st["core_hours_remaining"] is not None else ""))

    print(f"\n  {'machine':>8s} {'id':>18s} {'assigned':>9s} {'done':>6s} "
          f"{'left':>6s} {'ETA @'+str(st['workers_assumed'])+'w':>12s}")
    for k, v in st["per_machine"].items():
        mid = (v["machine_id"] or ["-"])[0]
        print(f"  {k:>8d} {str(mid)[:18]:>18s} {v['assigned']:>9d} "
              f"{v['done']:>6d} {v['remaining']:>6d} "
              f"{_fmt_hms(st['eta_s'].get(k)):>12s}")

    for f, vals in st["provenance"].items():
        if len(vals) > 1:
            print(f"\n  *** PROVENANCE SPLIT on {f} - results are NOT comparable:")
            for v, ms in vals.items():
                print(f"        {str(v)[:16]}  <- {ms}")

    if targets and st["_admis"]:
        rows = observable_summary(st["_admis"], targets)
        print(f"\n  {'observable':>16s} {'n':>5s} {'p05':>11s} {'median':>11s} "
              f"{'p95':>11s} {'target':>11s} {'in band':>8s}")
        for r in rows:
            if r["n"] == 0:
                print(f"  {r['observable']:>16s} {0:>5d}  (no admissible rows yet)")
                continue
            fb = ("n/a" if r["frac_in_band"] is None
                  else f"{100*r['frac_in_band']:.0f}%")
            tg = (f"{'n/a':>11s}" if r["target"] is None
                  else f"{r['target']:11.3g}")
            print(f"  {r['observable']:>16s} {r['n']:>5d} {r['p05']:11.3g} "
                  f"{r['median']:11.3g} {r['p95']:11.3g} {tg} {fb:>8s}")
        print("\n  'in band' = fraction of admissible runs inside the experimental")
        print("  range. Early in a campaign a low value is normal - the prior box")
        print("  is wide on purpose. A value stuck at 0% once several hundred rows")
        print("  are in means the box does not contain the data, which is a result.")


def watch(design: Path, results_dir: Path, n_machines: int = 4,
          workers_per_machine: int | None = None, interval: float = 60.0,
          iterations: int = 10_000, targets: dict | None = None) -> dict:
    """Live monitor. Interrupt the cell (kernel stop) to leave it."""
    from IPython.display import clear_output
    targets = targets if targets is not None else load_targets()
    st = None
    try:
        for _ in range(iterations):
            st = campaign_status(design, results_dir, n_machines,
                                 workers_per_machine)
            clear_output(wait=True)
            print(f"  refreshed {time.strftime('%H:%M:%S')}  "
                  f"(every {interval:.0f}s; interrupt the kernel to stop watching)\n")
            render_status(st, targets)
            if st["missing"] == 0:
                print("\n  CAMPAIGN COMPLETE.")
                break
            time.sleep(interval)
    except KeyboardInterrupt:
        print("\n  (stopped watching - the campaign itself is unaffected)")
    return st
