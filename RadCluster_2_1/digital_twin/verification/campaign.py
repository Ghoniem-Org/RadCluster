#!/usr/bin/env python3
"""
campaign.py — run the verification study, one claimed run at a time.

    # take runs until none are left, syncing after each
    python campaign.py --loop

    # a single named run
    python campaign.py --run T4_C4

    # what is the board doing
    python campaign.py --status

Each run produces the SAME artifact set as the tracked reference run --
`output/<stamp>_<tag>_<mode>_<physics>_I..V.._im..vm../` with provenance.md,
summary.csv, diagnostics.txt, results_*.npy and plots/ -- because it goes
through `RadClusterSimulation._save_output`, the same code path the reference
came from.  The `<tag>` is the run_id (T1_B3, T4_D, ...); without it the four
Table-1 rungs differ only in `i_discrete` and would produce four directory
names distinguishable only by timestamp.

LIVE PROGRESS.  The C++ solver emits a `[diag]` line per output step; the
progress callback turns that into two local files that update DURING a run:

    verification/logs/<machine>.jsonl   append-only event log
    verification/logs/<machine>.status  a rendered snapshot, rewritten per step

so `cat`/`tail -f` answers "where is this run" without waiting for it to end.
The six observables are NOT available mid-run -- they need the full state, and
the callback carries only the solver's own diagnostic scalars -- so the live
view reports position (dose reached of target, elapsed, step) and the final
record carries the physics.  Saying otherwise would be inventing numbers.

MULTI-MACHINE.  Runs are claimed through `gitsync.CampaignSync`: the git remote
is the lock, one file per run, push races resolved by retry.  After every run
the result is committed and pushed, so any machine (or the user, anywhere) sees
it immediately.  Local output directories stay local -- `output/` is gitignored
apart from the pinned reference -- so what travels is the compact record:
provenance.md, summary.csv, diagnostics.txt and the observables JSON.
"""
from __future__ import annotations

import argparse
import io
import json
import os
import sys
import threading
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
DT = HERE.parent                     # digital_twin/
MOD = DT.parent                      # RadCluster_2_1/
REPO = MOD.parent
for p in (str(REPO), str(DT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import runs as manifest_mod                     # noqa: E402
from gitsync import CampaignSync, machine_id    # noqa: E402

LOG_DIR = HERE / "logs"


# ── live progress log ────────────────────────────────────────────────────────
class ProgressLog:
    """Append-only event log plus a rendered status snapshot.

    Both are local and gitignored: they update per solver output step, which is
    far too often to commit.  What reaches the remote is the per-run record,
    pushed once the run ends.
    """

    def __init__(self, machine: str, run_id: str, entry: dict):
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        self.jsonl = LOG_DIR / f"{machine}.jsonl"
        # PER-RUN status file, not per-machine.  A single <machine>.status is
        # overwritten by whatever ran most recently -- a smoke test clobbered a
        # live T4_D status the moment it was introduced, and two concurrent runs
        # on one machine would fight over it continuously.  `tail -f *.status`
        # and the --status view both glob, so nothing is lost by splitting.
        self.status = LOG_DIR / f"{machine}.{run_id}.status"
        self.machine, self.run_id, self.entry = machine, run_id, entry
        self.t0 = time.time()
        self.G = None
        self.last = {}
        self.n_steps = 0
        self._t_step = time.time()
        self._stop = threading.Event()
        self._beat = None

    # HEARTBEAT.  render() otherwise fires only on a solver output step, and the
    # steps are log-spaced: T4_D reached step 18 in 4 s and then spent a long
    # time on step 19.  A status file frozen at "elapsed 0.1 min" is
    # indistinguishable from a dead process at exactly the moment the user most
    # wants to know the difference -- which defeats the point of a live log.
    # A daemon thread re-renders on a timer, so `elapsed` and `updated` keep
    # moving while `steps` stands still: working-but-slow now looks different
    # from stopped.
    def start_heartbeat(self, period_s=30.0):
        def _loop():
            while not self._stop.wait(period_s):
                try:
                    self.render(state="running")
                except Exception:
                    pass          # a status file must never kill a run
        self._beat = threading.Thread(target=_loop, daemon=True)
        self._beat.start()

    def stop_heartbeat(self):
        self._stop.set()
        if self._beat is not None:
            self._beat.join(timeout=2.0)

    def event(self, kind, **kw):
        rec = {"ts": time.strftime("%Y-%m-%dT%H:%M:%S"), "machine": self.machine,
               "run_id": self.run_id, "event": kind,
               "elapsed_s": round(time.time() - self.t0, 1), **kw}
        # Append by PATH and fsync, never through a held handle: the campaign
        # runs for days and anything that replaces the file underneath would
        # otherwise leave us writing into an unlinked inode (the failure
        # run_ensemble.py documents at length).
        with self.jsonl.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(rec, default=str) + "\n")
            fh.flush()
            os.fsync(fh.fileno())
        return rec

    def tick(self, row: dict):
        """One C++ solver output step."""
        self.n_steps += 1
        t = float(row.get("t", 0.0) or 0.0)
        dose = t * self.G if self.G else None
        self._t_step = time.time()
        self.last = {"t_s": t, "dose": dose, "step": self.n_steps}
        self.render(state="running")
        # The event log gets every step; the rendered status is what a human
        # tails.  37 steps per run makes this cheap.
        self.event("tick", t_s=t, dose=dose, step=self.n_steps)

    def render(self, state, extra=""):
        e = self.entry
        el = time.time() - self.t0
        tgt = float(e.get("dose", 0) or 0)
        dose = self.last.get("dose")
        frac = (dose / tgt * 100.0) if (dose and tgt) else 0.0
        bar_w = 34
        filled = max(0, min(bar_w, int(bar_w * frac / 100.0)))
        bar = "#" * filled + "." * (bar_w - filled)
        lines = [
            f"machine   : {self.machine}",
            f"run       : {self.run_id}   ({e.get('table')} / {e.get('label')})",
            f"state     : {state}",
            f"grid      : {e.get('equations')}  I={e.get('I')} V={e.get('V')}  "
            f"i_d={e.get('i_discrete', '-')} I_bin={e.get('I_bin', '-')} "
            f"v_d={e.get('v_discrete', '-')} V_bin={e.get('V_bin', '-')}  "
            f"shape={e.get('shape_function')}",
            f"dose      : [{bar}] "
            + (f"{dose:.4g} / {tgt:g} dpa ({frac:.1f}%)" if dose else
               f"-- / {tgt:g} dpa"),
            # Not "of n_points": run_adaptive drives the solver in segments and
            # the callback fires per output step of EACH, so the running total
            # legitimately exceeds the grid size (it read "9 of 8" the first
            # time).  Dose is the honest progress measure; this is a liveness
            # indicator -- if it stops advancing, the solver is stuck.
            f"steps     : {self.n_steps} solver output steps"
            f"  (grid: {e.get('n_points')} points)",
            f"step age  : {time.time() - self._t_step:.0f} s since the last step"
            f"{'  <-- long step, still working' if time.time() - self._t_step > 300 else ''}",
            f"elapsed   : {el/60:.1f} min",
            f"updated   : {time.strftime('%Y-%m-%d %H:%M:%S')}",
        ]
        if extra:
            lines.append("")
            lines.append(extra)
        tmp = self.status.with_suffix(".tmp")
        tmp.write_text("\n".join(lines) + "\n", encoding="utf-8")
        os.replace(tmp, self.status)     # atomic: a tailer never sees a partial


# ── building and running one manifest entry ──────────────────────────────────
def build_sim(e: dict):
    """A RadClusterSimulation configured from one manifest entry."""
    from RadCluster_2_1.py_utils.simulation import RadClusterSimulation
    sim = RadClusterSimulation(
        I=e["I"], V=e["V"], solver_mode=e["solver_mode"],
        equations=e["equations"], cascade=e["cascade"],
        he_kinetics=e["he_kinetics"],
        i_mobile=e["i_mobile"], v_mobile=e["v_mobile"])
    # Bin layout through run_ensemble's own helper, so the workbook key names
    # live in exactly one place.  It MUST run after the constructor and before
    # _calculate_derived()/rebuild_rates(): passing i_discrete= to the
    # constructor raises, and before 2026-08-02 was silently dropped, which is
    # how nine bit-identical "refinement" runs got produced.
    import run_ensemble as _re
    _re.apply_bin_config(sim, {
        "equations": e["equations"],
        "i_discrete": e.get("i_discrete", 0), "v_discrete": e.get("v_discrete", 0),
        "I_bin": e.get("I_bin", 0), "V_bin": e.get("V_bin", 0),
        "shape_function": e["shape_function"]})
    sim.input_data.reactions["LOOP_NETWORK_LOSS"] = int(e.get("lnl", 1))
    sim.input_data._calculate_derived()
    sim.rebuild_rates()
    return sim


def solver_config(e: dict, sim) -> dict:
    G = float(sim.input_data.reactions["G"])
    return {
        "t_span": (1e-6, float(e["dose"]) / G),
        "n_points": int(e["n_points"]),
        "log_time": True,
        "rtol": float(e["rtol"]), "atol": float(e["atol"]),
        "timeout_s": float(e.get("timeout_s", 86400)),
        "solver_method": {"linsol": "gmres",
                          "preconditioner": e.get("preconditioner", "woodbury"),
                          "concentration_threshold": 1e-22},
        "loop_conversion": 1,
    }


def execute(e: dict, log: ProgressLog, save_plots=True) -> dict:
    """Run one entry to completion.  Returns the record to publish."""
    import run_ensemble as re_mod

    rec = {"run_id": e["run_id"], "table": e["table"], "label": e["label"],
           "machine": log.machine}
    t0 = time.time()
    if e.get("omp_threads"):
        os.environ["OMP_NUM_THREADS"] = str(e["omp_threads"])

    sim = build_sim(e)
    sim.run_tag = e["run_id"]            # -> output/<stamp>_<run_id>_<...>/
    G = float(sim.input_data.reactions["G"])
    log.G = G
    scfg = solver_config(e, sim)

    cfg = {"I": e["I"], "V": e["V"], "dose": float(e["dose"]),
           "C_floor": float(sim.input_data.reactions.get("C_floor", 1e-15)),
           "equations": e["equations"],
           "i_discrete": e.get("i_discrete", 0), "I_bin": e.get("I_bin", 0),
           "v_discrete": e.get("v_discrete", 0), "V_bin": e.get("V_bin", 0),
           "shape_function": e["shape_function"]}

    # Realised bin layout, read back off rate_equations -- plan S3.3 requires
    # the tables to report this, not the requested values.  It also RAISES on a
    # silent fallback to the coarse default, which is the failure this study
    # cannot afford to discover after the fact.
    if e["equations"] == "bin_moment":
        rec.update(re_mod.bin_layout(sim, cfg))

    log.render(state="starting")
    log.start_heartbeat()
    log.event("start", **{k: e.get(k) for k in
                          ("table", "label", "equations", "I", "V", "dose",
                           "i_discrete", "I_bin", "v_discrete", "V_bin",
                           "shape_function", "rtol", "preconditioner",
                           "he_kinetics", "cascade", "omp_threads")})

    # The solver's own chatter goes to a buffer; the progress callback is what
    # reaches the log.  Without this the terminal is unusable for a long loop.
    buf = io.StringIO()
    saved = sys.stdout, sys.stderr
    try:
        sys.stdout = sys.stderr = buf
        if int(e.get("lnl", 1)):
            # LOOP_NETWORK_LOSS is DEAD under plain run(): v_net is built from
            # the segment-frozen monomers, which only run_adaptive writes.
            # max_doublings=0 takes the operator splitting and refuses the
            # domain growth, which would otherwise leave rungs of one table on
            # different grids.
            res = sim.run_adaptive(solver_config=scfg, save_output=save_plots,
                                   progress_callback=log.tick,
                                   timeout_s=scfg["timeout_s"], max_doublings=0)
        else:
            res = sim.run(solver_config=scfg, save_output=save_plots,
                          progress_callback=log.tick)
    finally:
        sys.stdout, sys.stderr = saved
        log.stop_heartbeat()

    if res is None:
        rec.update({"status": "failed", "error": "solver returned None"})
        rec["wall_s"] = round(time.time() - t0, 1)
        return rec

    # Six observables through the campaign's own operator, so a verification
    # table and a ledger row for the same configuration agree by construction.
    rec.update(re_mod.observe(res, sim, cfg, 1.0))
    rec["status"] = "done"
    rec["wall_s"] = round(time.time() - t0, 1)
    rec["dose_target"] = float(e["dose"])

    # S3.4.1: a run that missed its comparison dose is reported as "did not
    # reach", never as a metric.  Recorded on the row so the reporting step can
    # enforce it mechanically rather than by discipline.
    reads = e["dose_read"]
    reads = reads if isinstance(reads, (list, tuple)) else [reads]
    lad = rec.get("at_dose") or {}
    rec["scored"] = {}
    for target in reads:
        key = f"{float(target):g}"
        got = lad.get(key)
        rec["scored"][key] = got if got else {"missing": True,
                                              "dose_reached": rec.get("dose_reached")}
    rec["reached_all_scoring_doses"] = all(
        not v.get("missing") for v in rec["scored"].values())

    out_dir = getattr(sim, "_last_output_dir", None)
    if out_dir:
        rec["output_dir"] = str(out_dir)
    return rec


def artifacts_for(rec: dict) -> dict:
    """The compact record that travels to the remote (not the PNGs, not the
    38 MB plot_data.pkl -- .gitignore excludes the latter for good reasons)."""
    art = {"observables.json": json.dumps(rec, indent=2, sort_keys=True,
                                          default=str) + "\n"}
    od = rec.get("output_dir")
    if od:
        d = Path(od)
        for name in ("provenance.md", "summary.csv", "diagnostics.txt"):
            f = d / name
            if f.exists():
                try:
                    art[name] = f.read_text(encoding="utf-8", errors="replace")
                except Exception:
                    pass
    return art


# ── driver ───────────────────────────────────────────────────────────────────
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", help="run this run_id and stop")
    ap.add_argument("--loop", action="store_true",
                    help="keep claiming runs until the board is empty")
    ap.add_argument("--status", action="store_true", help="print the board")
    ap.add_argument("--list", action="store_true", help="print the manifest")
    ap.add_argument("--release", help="drop a claim so another machine can take it")
    ap.add_argument("--branch", default="campaign-verification")
    ap.add_argument("--only", default=None,
                    help="restrict to a table, e.g. --only T4")
    ap.add_argument("--offline", action="store_true",
                    help="no fetch/push; claims stay local (smoke tests)")
    ap.add_argument("--no-plots", action="store_true",
                    help="skip the output directory entirely (diagnostics only)")
    ap.add_argument("--max-runs", type=int, default=0,
                    help="stop after N runs (0 = no limit)")
    ap.add_argument("--smoke", action="store_true",
                    help="tiny end-to-end run (seconds) exercising the real "
                         "code path: output directory, live log, claim, sync")
    a = ap.parse_args(argv)

    entries = {e["run_id"]: e for e in manifest_mod.manifest()}
    order = [e["run_id"] for e in manifest_mod.runnable()]

    if a.smoke:
        # Deliberately NOT in the manifest: it is a pipeline test, and a study
        # run_id on the shared board would have to be explained forever after.
        smoke = dict(manifest_mod.BASE)
        smoke.update({
            "run_id": "SMOKE", "table": "--", "label": "pipeline smoke test",
            "I": 200, "V": 100, "i_discrete": 10, "I_bin": 6,
            "v_discrete": 5, "V_bin": 6,
            # 0.01 dpa clears DOSE_CHECKPOINTS[0] = 0.005, so the ladder has a
            # rung and the six observables are actually exercised.  At 1e-3 the
            # scoring block returned {"missing": True} and every observable
            # printed nan -- which is correct behaviour, but tests nothing.
            "dose": 1e-2, "n_points": 10, "dose_read": 0.005,
            "timeout_s": 600, "omp_threads": 2,
        })
        entries["SMOKE"] = smoke
        order = ["SMOKE"]
        a.run = "SMOKE"
    if a.only:
        order = [r for r in order if entries[r]["table"] == a.only]

    if a.list:
        manifest_mod.__dict__["__name__"] = "__main__"
        os.system(f"{sys.executable} {HERE/'runs.py'}")
        return 0

    sync = CampaignSync(REPO, branch=a.branch, offline=a.offline)
    sync.ensure()

    if a.release:
        ok = sync.release(a.release)
        print(f"release {a.release}: {'ok' if ok else 'nothing to release'}")
        return 0

    if a.status:
        print_status(sync, entries, order)
        return 0

    machine = machine_id()
    done_count = 0
    while True:
        rid = a.run if a.run else sync.next_unclaimed(order)
        if rid is None:
            print("board empty — every run is claimed or finished.")
            break
        if rid not in entries:
            print(f"unknown run_id {rid!r}", file=sys.stderr)
            return 2
        e = entries[rid]
        # The smoke test never touches the shared board.  It is a pipeline
        # check, and a "SMOKE done" row sitting among the study runs would have
        # to be explained to every later reader of the board.
        if a.smoke:
            print(f"\n=== SMOKE (pipeline test, not published) on {machine} ===",
                  flush=True)
            log = ProgressLog(machine, "SMOKE", e)
            rec = execute(e, log, save_plots=not a.no_plots)
            log.render(state=rec.get("status", "?").upper(), extra=summary_line(rec))
            print(f"  {summary_line(rec)}", flush=True)
            print("  board not touched (smoke)", flush=True)
            return 0
        if not sync.claim(rid, meta={"table": e["table"], "label": e["label"]}):
            print(f"  {rid}: claimed by another machine, skipping")
            if a.run:
                return 0
            continue

        print(f"\n=== {rid}  ({e['table']} / {e['label']}) on {machine} ===",
              flush=True)
        log = ProgressLog(machine, rid, e)
        try:
            rec = execute(e, log, save_plots=not a.no_plots)
        except BaseException as exc:      # a failed run must still be published
            rec = {"run_id": rid, "table": e["table"], "label": e["label"],
                   "status": "failed", "machine": machine,
                   "error": f"{type(exc).__name__}: {exc}"[:300],
                   "traceback_tail": traceback.format_exc()[-800:]}
            log.event("error", error=rec["error"])
            log.render(state="FAILED", extra=rec["error"])
            if isinstance(exc, KeyboardInterrupt):
                sync.publish(rid, rec, artifacts_for(rec))
                raise

        log.event("end", status=rec.get("status"),
                  wall_s=rec.get("wall_s"),
                  dose_reached=rec.get("dose_reached"))
        log.render(state=rec.get("status", "?").upper(), extra=summary_line(rec))
        print(f"  {summary_line(rec)}", flush=True)

        # AUTOMATIC SYNC AFTER EVERY RUN, success or failure.
        pushed = sync.publish(rid, rec, artifacts_for(rec))
        print(f"  synced to {a.branch}: {'yes' if pushed else 'NO (will retry next run)'}",
              flush=True)

        done_count += 1
        if a.run or (a.max_runs and done_count >= a.max_runs) or not a.loop:
            break
    return 0


def summary_line(rec: dict) -> str:
    if rec.get("status") != "done":
        return f"{rec.get('status', '?')}: {rec.get('error', '')}"[:200]
    sc = (rec.get("scored") or {})
    first = next((v for v in sc.values() if not v.get("missing")), {})
    return (f"dose {rec.get('dose_reached', 0):.3g}/{rec.get('dose_target', 0):g} dpa  "
            f"N_eq={rec.get('N_eq', '-')}  "
            f"d100={first.get('d_100_nm', float('nan')):.3g} "
            f"d111={first.get('d_111_nm', float('nan')):.3g} "
            f"dcav={first.get('d_cavity_nm', float('nan')):.3g} nm  "
            f"dFP={rec.get('delta_FP', float('nan')):.2g}  "
            f"{rec.get('wall_s', 0)/60:.1f} min")


def print_status(sync, entries, order):
    held = sync.claims()
    print(f"\ncampaign board — branch {sync.branch}   "
          f"(this machine: {sync.machine})\n")
    hdr = f"{'run_id':16s} {'tbl':4s} {'state':9s} {'machine':12s} {'result':52s}"
    print(hdr); print("-" * len(hdr))
    n_done = n_run = n_free = 0
    for rid in order:
        c = held.get(rid)
        if not c:
            n_free += 1
            print(f"{rid:16s} {entries[rid]['table']:4s} {'free':9s} {'-':12s}")
            continue
        st = c.get("status", "?")
        n_done += st == "done"
        n_run += st == "running"
        print(f"{rid:16s} {entries[rid]['table']:4s} {st:9s} "
              f"{str(c.get('machine', '-'))[:12]:12s} {summary_line(c)[:52]}")
    print(f"\n{n_done} done, {n_run} running, {n_free} free, "
          f"of {len(order)} runnable\n")
    for f in sorted(LOG_DIR.glob("*.status")):
        print(f"--- live: {f.name} " + "-" * 30)
        print(f.read_text(encoding="utf-8"))


if __name__ == "__main__":
    sys.exit(main())
