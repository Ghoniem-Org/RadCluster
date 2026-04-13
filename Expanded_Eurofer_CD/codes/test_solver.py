"""Quick solver test — runs 2 short segments and prints CVODE diagnostics."""
import sys, os, io
from pathlib import Path

MODULE_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT   = MODULE_ROOT.parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(MODULE_ROOT))

import importlib
import Expanded_Eurofer_CD.py_utils.defect_production as _dp
import Expanded_Eurofer_CD.py_utils.binding_energies as _be
import Expanded_Eurofer_CD.py_utils.bin_moment_rates as _bmr
import Expanded_Eurofer_CD.py_utils.input_data as _inp
import Expanded_Eurofer_CD.py_utils.reaction_rates as _rr
import Expanded_Eurofer_CD.py_utils.rate_equations as _re
import Expanded_Eurofer_CD.py_utils.cpp_bridge as _cb
import Expanded_Eurofer_CD.py_utils.post_process as _pp
import Expanded_Eurofer_CD.py_utils.simulation as _sim
for _m in [_dp, _be, _bmr, _inp, _rr, _re, _cb, _pp, _sim]:
    importlib.reload(_m)
from Expanded_Eurofer_CD.py_utils.simulation import ExpandedEuroferCDSimulation

# ── Configuration (matches notebook) ──
I, V = 100_000, 100_000
i_mobile, v_mobile = 10, 5

_saved = sys.stdout, sys.stderr
sys.stdout = sys.stderr = io.StringIO()
try:
    sim = ExpandedEuroferCDSimulation(
        I=I, V=V, solver_mode='cpp_full',
        physics_option='bin_moment_CD_fission',
        C_floor=1e-25, he_options='quasi_steady_state',
        i_mobile=i_mobile, v_mobile=v_mobile)
finally:
    sys.stdout, sys.stderr = _saved

# ── Overrides ──
OVERRIDES = {
    'eta': 0.3, 'f_cl_i': 0.58, 'f_cl_v': 0.45,
    'E_m_1D': 0.34, 'i_mobile': i_mobile, 'L_hat': 71.8,
    'c_C': 1.94e-4, 'E_b_C_SIA': 0.45, 'rho_d': 1e14,
    'Z_i': 1.1, 'Z_ii': 1.1, 'shape_function': 'lognormal',
    'i_discrete': i_mobile, 'v_discrete': v_mobile,
    'I_bin': 30, 'V_bin': 30,
}
sys.stdout = sys.stderr = io.StringIO()
try:
    inp = sim.input_data
    for k, v in OVERRIDES.items():
        placed = False
        for d in [inp.production_fission, inp.production_fusion,
                  inp.diffusion, inp.reactions, inp.energetics, inp.dissociation]:
            if k in d:
                d[k] = v; placed = True
        if not placed:
            inp.reactions[k] = v
    if 'i_mobile' in OVERRIDES:
        inp.diffusion['i_mobile'] = int(OVERRIDES['i_mobile'])
        inp.reactions['i_mobile'] = int(OVERRIDES['i_mobile'])
    inp._calculate_derived()
    sim.rebuild_rates()
finally:
    sys.stdout, sys.stderr = _saved

# ── Run just 2 short segments ──
SOLVER_CONFIG = {
    't_span': (1e-6, 3.16),      # segments 1-5 time range (up to original failure)
    'n_points': 50,
    'log_time': True,
    'rtol': 1e-6, 'atol': 1e-25,
    'solver_method': {
        'backend': 'cvode', 'lmm': 'bdf', 'linsol': 'gmres',
        'window_prec': 1, 'window_gmres_maxl': 20,
    }
}

print(f'N_eq = {sim.rate_equations.N_eq}')
print(f'Running t=[{SOLVER_CONFIG["t_span"][0]:.1e}, {SOLVER_CONFIG["t_span"][1]:.1e}]...')
print('(CVODE diagnostics will appear on stderr)\n', flush=True)

results = sim.run_adaptive(
    solver_config=SOLVER_CONFIG,
    save_output=False,
    boundary_threshold=0.1,
    max_doublings=0,
    points_per_segment=10,
)

if results is not None:
    print(f'\nDone: {len(results["t"])} time points')
    print(f't = [{results["t"][0]:.2e}, {results["t"][-1]:.2e}]')
    print(f'delta_FP = {results["delta_FP"][-1]:.3e}')
else:
    print('\nSimulation FAILED')
