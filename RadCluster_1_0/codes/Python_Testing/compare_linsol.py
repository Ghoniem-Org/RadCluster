#!/usr/bin/env python3
"""
compare_linsol.py — Linear-solver / preconditioner comparison.

Six configurations are timed back-to-back, all with identical physics
overrides (taken verbatim from the current cell of
RadCluster_1_0.ipynb). Each run is capped at TIMEOUT_S wall-clock
seconds; if the cap is hit, the C++ solver finalises the current
output step gracefully and the partial trajectory is kept.

Group 1  EQUATIONS='discrete'    (full per-size, I_bin=V_bin=0)
   D-WB   gmres + Woodbury
   D-JC   gmres + Jacobi
   D-KLU  klu (sparse direct, full_CD only)

Group 2  EQUATIONS='bin_moment'  (i_discrete=v_discrete=50, I_bin=V_bin=20)
   B-WB   gmres + Woodbury
   B-JC   gmres + Jacobi
   B-KLU  klu                       <-- C++ rejects: full_CD only

After the 6 runs the script writes
    compare_linsol_results.json
    compare_linsol_report.md
in the same directory.
"""
from __future__ import annotations

import sys, os, io, time, json, subprocess, importlib, traceback
from pathlib import Path
from datetime import datetime

# Force UTF-8 stdout so bridge/solver prints with non-ASCII (any "->" etc.)
# never crash the parse path on Windows cp1252 consoles.
os.environ.setdefault('PYTHONIOENCODING', 'utf-8')
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    sys.stderr.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

# ── Repo paths ────────────────────────────────────────────────────────────────
MODULE_ROOT = Path(__file__).resolve().parent.parent.parent  # RadCluster_1_0/
REPO_ROOT   = MODULE_ROOT.parent
HERE        = Path(__file__).resolve().parent
for p in (str(REPO_ROOT), str(MODULE_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

# ── Build C++ solver (skip if up-to-date) ─────────────────────────────────────
build_dir = MODULE_ROOT / 'build'
exe_name  = 'solver.exe' if sys.platform == 'win32' else 'solver'
exe_paths = [build_dir / 'Release' / exe_name,
             build_dir / 'Debug'   / exe_name,
             build_dir              / exe_name]
if not any(p.exists() for p in exe_paths):
    print("Building C++ solver ...", flush=True)
    build_dir.mkdir(exist_ok=True)
    cmake_src = MODULE_ROOT / 'cpp_utils'
    for cmd in (['cmake', '-S', str(cmake_src), '-B', str(build_dir),
                 '-DCMAKE_BUILD_TYPE=Release'],
                ['cmake', '--build', str(build_dir), '--config', 'Release']):
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print(res.stderr[-1500:]); sys.exit(1)
    print("Build OK\n", flush=True)

# ── Reload py_utils so any in-progress edits land ─────────────────────────────
import RadCluster_1_0.py_utils.defect_production as _dp
import RadCluster_1_0.py_utils.binding_energies  as _be
import RadCluster_1_0.py_utils.bin_moment_rates  as _bmr
import RadCluster_1_0.py_utils.input_data        as _inp
import RadCluster_1_0.py_utils.reaction_rates    as _rr
import RadCluster_1_0.py_utils.rate_equations    as _re
import RadCluster_1_0.py_utils.cpp_bridge        as _cb
import RadCluster_1_0.py_utils.post_process      as _pp
import RadCluster_1_0.py_utils.simulation        as _sim
for _m in (_dp, _be, _bmr, _inp, _rr, _re, _cb, _pp, _sim):
    importlib.reload(_m)
from RadCluster_1_0.py_utils.simulation import RadClusterSimulation

# ══════════════════════════════════════════════════════════════════════════════
# 1.  Parameters — mirror the notebook cell exactly
# ══════════════════════════════════════════════════════════════════════════════
TIMEOUT_S = 3000.0          # per-run wall-clock cap

# Domain (kept identical across both groups so physics is comparable)
I = V              = int(1e4)
i_mobile, v_mobile = 10, 2

# Time horizon and tolerances — copied from the notebook cell
T_SPAN   = (1e-6, 1e4)
N_POINTS = 200
RTOL     = 1e-6
ATOL     = 1e-20
C_FLOOR  = 1e-25

# Physics overrides — copied verbatim from the notebook cell
PHYSICS_OVERRIDES = {
    'T':         673,
    'eta':       0.3,
    'f_cl_i':    0.5,
    'f_cl_v':    0.15,
    'E_m_1D':    0.4,
    'L_hat':     5,
    'c_C':       1e-3,
    'E_b_C_SIA': 0.65,
    'rho_d':     1e14,
    'Z_i':       1.08,
    'Z_ii':      1.01,
}

# Bin-moment hybrid grid (user spec: I_discrete=v_discrete=50, I_bin=V_bin=20)
BIN_MOMENT_GRID = dict(i_discrete=50, v_discrete=50,
                       I_bin=20,     V_bin=20,
                       shape_function='linear')

# Six cases.  prec_type only consulted when linsol='gmres'.
CASES = [
    # tag,    equations,     linsol,   preconditioner
    ('D-WB',  'discrete',    'gmres',  'Woodbury'),
    ('D-JC',  'discrete',    'gmres',  'Jacobi'  ),
    ('D-KLU', 'discrete',    'klu',    None      ),
    ('B-WB',  'bin_moment',  'gmres',  'Woodbury'),
    ('B-JC',  'bin_moment',  'gmres',  'Jacobi'  ),
    ('B-KLU', 'bin_moment',  'klu',    None      ),  # expected: rejected by C++
]


# ══════════════════════════════════════════════════════════════════════════════
# 2.  Helpers
# ══════════════════════════════════════════════════════════════════════════════
def build_sim(equations: str) -> RadClusterSimulation:
    """Construct a fresh RadClusterSimulation with the right equation form."""
    saved_out, saved_err = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = io.StringIO()
    try:
        sim = RadClusterSimulation(
            I=I, V=V,
            solver_mode='full_system',          # Woodbury requires full_system
            equations=equations,
            cascade='fission',
            C_floor=C_FLOOR,
            he_kinetics='quasi_steady_state',
            i_mobile=i_mobile, v_mobile=v_mobile,
        )
        # Inject overrides
        ov = dict(PHYSICS_OVERRIDES)
        if equations == 'discrete':
            ov['i_discrete'] = I        # full per-size
            ov['v_discrete'] = V
            ov['I_bin']      = 0
            ov['V_bin']      = 0
            ov['shape_function'] = 'linear'
        else:
            ov.update(BIN_MOMENT_GRID)

        inp = sim.input_data
        for key, val in ov.items():
            placed = False
            for d in (inp.production_fission, inp.production_fusion,
                      inp.diffusion, inp.reactions,
                      inp.energetics, inp.dissociation):
                if key in d:
                    d[key] = val
                    placed = True
            if not placed:
                inp.reactions[key] = val
        for mob_key in ('i_mobile', 'v_mobile'):
            v_int = int(ov.get(mob_key, inp.diffusion.get(mob_key, 0)))
            inp.diffusion[mob_key] = v_int
            inp.reactions[mob_key] = v_int
        inp._calculate_derived()
        sim.rebuild_rates()
    finally:
        sys.stdout, sys.stderr = saved_out, saved_err
    return sim


def make_solver_config(linsol: str, preconditioner) -> dict:
    method = {
        'linsol':                  linsol,
        'window_width':            10,
        'concentration_threshold': 1e-22,
        'window_pad':              20,
    }
    if linsol == 'gmres' and preconditioner is not None:
        method['preconditioner'] = preconditioner
    return {
        't_span':   T_SPAN,
        'n_points': N_POINTS,
        'log_time': True,
        'rtol':     RTOL,
        'atol':     ATOL,
        'solver_method': method,
    }


def _last(arr):
    try:
        return float(arr[-1])
    except Exception:
        return None


def run_one(tag, equations, linsol, preconditioner):
    print(f"\n{'='*72}\n  CASE {tag}: equations={equations!r}  linsol={linsol!r}"
          f"  prec={preconditioner!r}\n{'='*72}", flush=True)

    record = dict(tag=tag, equations=equations, linsol=linsol,
                  preconditioner=preconditioner,
                  N_eq=None, wall_s=None, status='unknown',
                  t_final=None, dose_final=None, n_pts=None,
                  swelling_pct=None, delta_FP=None, delta_He=None,
                  mean_n_i=None, mean_n_v=None, error=None)

    try:
        sim = build_sim(equations)
        record['N_eq'] = int(sim.rate_equations.N_eq)
        cfg            = make_solver_config(linsol, preconditioner)

        t0 = time.perf_counter()
        results = sim.run(solver_config=cfg, save_output=False,
                          timeout_s=TIMEOUT_S)
        wall = time.perf_counter() - t0
        record['wall_s'] = wall

        if results is None:
            record['status'] = 'failed'
            print(f"  -> no results returned ({wall:.1f} s)", flush=True)
            return record

        # Authoritative timeout flag set by cpp_bridge when wall cap or Ctrl+C
        # tripped the graceful shutdown path; fall back to a wall-clock heuristic.
        record['status'] = ('timeout' if (results.get('partial')
                                          or wall >= TIMEOUT_S - 1)
                            else 'completed')
        t_arr = results.get('t', [])
        record['n_pts']        = int(len(t_arr))
        record['t_final']      = _last(t_arr)
        record['dose_final']   = _last(results.get('dose'))
        sw = _last(results.get('swelling'))
        record['swelling_pct'] = sw * 100.0 if sw is not None else None
        record['delta_FP']     = _last(results.get('delta_FP'))
        record['delta_He']     = _last(results.get('delta_He'))
        record['mean_n_i']     = _last(results.get('mean_n_i'))
        record['mean_n_v']     = _last(results.get('mean_n_v'))

        tf = record['t_final']
        tf_s = '—' if tf is None else f"{tf:.3e}"
        print(f"  -> status={record['status']}  wall={wall:.2f} s  "
              f"n_pts={record['n_pts']}  t_final={tf_s}", flush=True)
        return record

    except Exception as exc:
        record['status'] = 'error'
        record['error']  = f"{type(exc).__name__}: {exc}"
        traceback.print_exc()
        return record


# ══════════════════════════════════════════════════════════════════════════════
# 3.  Run the matrix
# ══════════════════════════════════════════════════════════════════════════════
records = [run_one(*c) for c in CASES]

# ── persist raw results ───────────────────────────────────────────────────────
out_json = HERE / 'compare_linsol_results.json'
with out_json.open('w') as f:
    json.dump({
        'timestamp':        datetime.now().isoformat(timespec='seconds'),
        'I': I, 'V': V, 'i_mobile': i_mobile, 'v_mobile': v_mobile,
        't_span': T_SPAN, 'n_points': N_POINTS, 'rtol': RTOL, 'atol': ATOL,
        'physics_overrides': PHYSICS_OVERRIDES,
        'bin_moment_grid':   BIN_MOMENT_GRID,
        'timeout_s':         TIMEOUT_S,
        'records':           records,
    }, f, indent=2, default=str)
print(f"\nRaw results -> {out_json}")


# ══════════════════════════════════════════════════════════════════════════════
# 4.  Auto-generate the markdown report
# ══════════════════════════════════════════════════════════════════════════════
def fmt_t(x):     return '—' if x is None else f"{x:.3e}"
def fmt_w(x):     return '—' if x is None else f"{x:8.1f}"
def fmt_pct(x):   return '—' if x is None else f"{x:.4f}"
def fmt_n(x):     return '—' if x is None else f"{int(x):d}"

def find(tag):
    for r in records:
        if r['tag'] == tag: return r
    return None

def ratio(a, b):
    if a is None or b is None or b == 0: return '—'
    return f"{a/b:.2f}×"

D_WB,  D_JC,  D_KLU = find('D-WB'),  find('D-JC'),  find('D-KLU')
B_WB,  B_JC,  B_KLU = find('B-WB'),  find('B-JC'),  find('B-KLU')

# Pick the fastest *successful* run in each group as the in-group reference
def fastest(group):
    ok = [r for r in group if r['status'] == 'completed' and r['wall_s']]
    return min(ok, key=lambda r: r['wall_s']) if ok else None

ref_d = fastest([D_WB, D_JC, D_KLU])
ref_b = fastest([B_WB, B_JC, B_KLU])

lines = []
lines.append("# Linear-solver / preconditioner comparison\n")
lines.append(f"_Generated {datetime.now().isoformat(timespec='seconds')}_\n")
lines.append("## 1. Setup\n")
lines.append(f"- **Domain:** `I = V = {I}`,  `i_mobile = {i_mobile}`,  `v_mobile = {v_mobile}`")
lines.append(f"- **Time horizon:** `t_span = {T_SPAN}`,  `n_points = {N_POINTS}`,  "
             f"`rtol = {RTOL:.0e}`,  `atol = {ATOL:.0e}`")
lines.append(f"- **Solver mode:** `full_system`,  cascade = `fission`,  he_kinetics = `quasi_steady_state`")
lines.append(f"- **Wall-clock cap per run:** `TIMEOUT_S = {TIMEOUT_S:.0f} s`  "
             f"(C++ child finalises gracefully on cap, partial bin file kept)")
lines.append("- **Group 1 (`discrete`):** full per-size — `i_discrete = I`, `v_discrete = V`, "
             "`I_bin = V_bin = 0`")
lines.append(f"  - `N_eq` = {D_WB['N_eq'] if D_WB else '—'}")
lines.append("- **Group 2 (`bin_moment`):** hybrid — "
             "`i_discrete = v_discrete = 50`, `I_bin = V_bin = 20`, `shape_function = 'linear'` (P=2)")
lines.append(f"  - `N_eq` = {B_WB['N_eq'] if B_WB else '—'}\n")
lines.append("**Physics overrides (identical for all 6 runs):**\n")
lines.append("```python")
lines.append("PHYSICS_OVERRIDES = " + json.dumps(PHYSICS_OVERRIDES, indent=4))
lines.append("```\n")

# ── Results table ────────────────────────────────────────────────────────────
lines.append("## 2. Raw results\n")
lines.append("| Tag   | Equations    | linsol | precond  | N_eq  | status     | wall (s) | n_pts | t_final (s) | dose (dpa) | swelling (%) | δ_FP    | δ_He    |")
lines.append("|-------|--------------|--------|----------|------:|------------|---------:|------:|-------------|------------|--------------|---------|---------|")
for r in records:
    lines.append(
        f"| {r['tag']:<5} | {r['equations']:<12} | {r['linsol']:<6} | "
        f"{(r['preconditioner'] or '—'):<8} | "
        f"{r['N_eq'] if r['N_eq'] else '—':>5} | {r['status']:<10} | "
        f"{fmt_w(r['wall_s']):>8} | {fmt_n(r['n_pts']):>5} | "
        f"{fmt_t(r['t_final'])} | {fmt_t(r['dose_final'])} | "
        f"{fmt_pct(r['swelling_pct']):>12} | "
        f"{fmt_t(r['delta_FP']):>7} | {fmt_t(r['delta_He']):>7} |"
    )
lines.append("")

# ── Intra-group analysis ─────────────────────────────────────────────────────
lines.append("## 3. Intra-group comparison\n")
lines.append("### 3.1 `discrete` (Group 1)\n")
lines.append(f"Reference (fastest completed): **{ref_d['tag'] if ref_d else 'none'}** "
             f"({fmt_w(ref_d['wall_s']) if ref_d else '—'} s).\n")
lines.append("| Tag   | precond  | wall (s) | speed-up vs ref | status     |")
lines.append("|-------|----------|---------:|-----------------|------------|")
for r in (D_WB, D_JC, D_KLU):
    if r is None: continue
    rel = ratio(r['wall_s'], ref_d['wall_s']) if (ref_d and r['wall_s']) else '—'
    lines.append(f"| {r['tag']:<5} | {(r['preconditioner'] or r['linsol']):<8} | "
                 f"{fmt_w(r['wall_s'])} | {rel} | {r['status']} |")
lines.append("")
lines.append(
    "_Notes._  Within Group 1 the C++ KLU path (linsol=`klu`) builds the exact "
    "sparsity pattern of the per-size Jacobian and uses a colored finite-difference "
    "Jacobian, giving sparse-direct factorisation. GMRES + Woodbury exploits the "
    "bordered-banded structure $J = T + UV^T$ with band half-width "
    f"`{2*max(i_mobile, v_mobile)+1}` and rank `{i_mobile+v_mobile}`. GMRES + Jacobi "
    "is the legacy diagonal preconditioner — included as a baseline.")
lines.append("")

lines.append("### 3.2 `bin_moment` (Group 2)\n")
lines.append(f"Reference (fastest completed): **{ref_b['tag'] if ref_b else 'none'}** "
             f"({fmt_w(ref_b['wall_s']) if ref_b else '—'} s).\n")
lines.append("| Tag   | precond  | wall (s) | speed-up vs ref | status     |")
lines.append("|-------|----------|---------:|-----------------|------------|")
for r in (B_WB, B_JC, B_KLU):
    if r is None: continue
    rel = ratio(r['wall_s'], ref_b['wall_s']) if (ref_b and r['wall_s']) else '—'
    lines.append(f"| {r['tag']:<5} | {(r['preconditioner'] or r['linsol']):<8} | "
                 f"{fmt_w(r['wall_s'])} | {rel} | {r['status']} |")
lines.append("")
lines.append(
    "_Notes._  KLU is **not wired for `bin_moment`** in the C++ solver — the "
    "sparsity-pattern builder is `full_CD`-only (`solver.cpp` lines 386-426), so "
    "B-KLU is expected to fail at startup with the message "
    "`[KLU] linsol=klu only supports full_CD modes`. The two GMRES variants are "
    "the practical options for the bin-moment formulation.")
lines.append("")

# ── Inter-group analysis ─────────────────────────────────────────────────────
lines.append("## 4. Inter-group comparison\n")
lines.append("Compares the cheapest successful run in each group at matched physics.\n")
if ref_d and ref_b:
    spd = ratio(ref_d['wall_s'], ref_b['wall_s'])
    lines.append(f"| Group        | best tag | wall (s) | N_eq | swelling (%) | δ_FP    | δ_He    |")
    lines.append(f"|--------------|----------|---------:|-----:|--------------|---------|---------|")
    for label, r in (('discrete',   ref_d), ('bin_moment', ref_b)):
        lines.append(
            f"| {label:<12} | {r['tag']:<8} | {fmt_w(r['wall_s'])} | "
            f"{r['N_eq']:>4} | {fmt_pct(r['swelling_pct']):>12} | "
            f"{fmt_t(r['delta_FP']):>7} | {fmt_t(r['delta_He']):>7} |"
        )
    lines.append("")
    lines.append(f"**Speed ratio (discrete / bin_moment):** {spd}.\n")
    lines.append("**Solution agreement:** swelling and conservation diagnostics should "
                 "match between the two formulations to within the bin-moment closure "
                 "error (`linear` shape, P=2 → truncation error $O((r-1)^3)$). "
                 "A swelling discrepancy above ~1 % (relative) or a δ_FP / δ_He gap "
                 "above $10^{-3}$ is a red flag for the binning grid.\n")
else:
    lines.append("_Inter-group comparison skipped — at least one group had no completed run._\n")

# ── Recommendations ──────────────────────────────────────────────────────────
lines.append("## 5. Recommendations (data-driven, fill in after inspection)\n")
lines.append("- **Best linsol for `discrete`:** _e.g._ KLU when the Jacobian fits in cache; "
             "Woodbury when N_eq exceeds a few thousand.")
lines.append("- **Best linsol for `bin_moment`:** Woodbury or Jacobi GMRES — pick from the "
             "Group-2 wall-clock column above.")
lines.append("- **Recommended production setup:** the formulation/linsol pair with the "
             "lowest wall-clock that also keeps δ_FP, δ_He < $10^{-6}$.")
lines.append("")

report_path = HERE / 'compare_linsol_report.md'
report_path.write_text('\n'.join(lines), encoding='utf-8')
print(f"Report      -> {report_path}")
