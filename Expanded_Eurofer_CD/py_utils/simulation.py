"""
simulation.py — Expanded_Eurofer_CD simulation orchestrator.

Responsibilities
----------------
1. Load inputs and initialise physics modules.
2. Dispatch to the appropriate solver mode and physics option.
3. Write timestamped output directory with provenance, CSV, and plots.

Solver modes
------------
cpp_full        → C++ SUNDIALS CVODE BDF, full system, via cpp_bridge.run_cpp_solver
cpp_sliding_win → C++ SUNDIALS CVODE BDF with sliding SIA window, via cpp_bridge
sliding_OpenMP  → C++ sliding window + OpenMP, via cpp_bridge

Physics options
---------------
full_CD_fission      → RateEquations he_mode='case2' (Eq. 175)
full_CD_fusion       → RateEquations he_mode='case1' (Eq. 174)
bin_moment_CD_fission → BinMomentRateEquations he_mode='case2'
bin_moment_CD_fusion  → BinMomentRateEquations he_mode='case1'
"""

import time as _time
import subprocess
import sys
from pathlib import Path
from datetime import datetime

import numpy as np

from .input_data      import InputData
from .reaction_rates  import ReactionRates
from .rate_equations  import RateEquations
from .bin_moment_rates import BinMomentRateEquations
from .                import post_process


BASE_DIR = Path(__file__).parent.parent


class ExpandedEuroferCDSimulation:
    """
    Orchestrates the Expanded_Eurofer_CD simulation workflow.

    Parameters
    ----------
    N              : int, optional  — override max SIA cluster size
    M              : int, optional  — override max vacancy cluster size
    solver_mode    : str, optional  — 'cpp_full' | 'cpp_sliding_win' | 'sliding_OpenMP'
    physics_option : str, optional  — 'full_CD_fission' | 'full_CD_fusion' |
                                      'bin_moment_CD_fission' | 'bin_moment_CD_fusion'
    excel_file     : path-like, optional
    C_floor        : float, optional — concentration floor (default 1e-15).
                                       Any state variable below this is clamped
                                       before evaluating rate terms; derivatives
                                       of clamped variables are also clamped ≥ 0.
    he_options     : str, optional   — 'dynamic' (default) integrates free He as a
                                       full ODE (Eq. 157); 'quasi_steady_state'
                                       eliminates it via dc_h/dt = 0 (valid because
                                       E_m_h = 0.06 eV gives rapid equilibration).
    """

    def __init__(self, N=None, M=None, solver_mode=None,
                 physics_option=None, excel_file=None,
                 C_floor=None, he_options=None):
        print("Initializing Expanded_Eurofer_CD simulation…")
        kwargs = {}
        if excel_file is not None:
            kwargs['excel_file'] = excel_file

        self.input_data = InputData(
            N=N, M=M,
            solver_mode=solver_mode,
            physics_option=physics_option,
            **kwargs
        )

        # Inject optional overrides before building rate equations
        if C_floor is not None:
            self.input_data.reactions['C_floor'] = float(C_floor)
        if he_options is not None:
            validated = str(he_options).lower()
            if validated not in ('dynamic', 'quasi_steady_state'):
                import warnings
                warnings.warn(f"Unknown he_options='{he_options}'. Using 'dynamic'.")
                validated = 'dynamic'
            self.input_data.reactions['he_options'] = validated

        self.reaction_rates = ReactionRates(self.input_data)

        po = self.input_data.physics_option
        if 'bin_moment' in po:
            self.rate_equations = BinMomentRateEquations(
                self.input_data, self.reaction_rates
            )
        else:
            self.rate_equations = RateEquations(
                self.input_data, self.reaction_rates
            )

        print(f"Simulation initialized: solver_mode='{self.input_data.solver_mode}'"
              f"  physics_option='{po}'")

    # ── Main entry point ──────────────────────────────────────────────────────

    def run(self, solver_config=None, save_output=True):
        """
        Run the simulation using the configured solver mode and physics option.

        Parameters
        ----------
        solver_config : dict, optional
            Keys: t_span, n_points, rtol, atol, log_time,
                  solver_method (dict with backend, linsol, window_*, etc.)
        save_output : bool
            Write timestamped output/ directory.

        Returns
        -------
        results : dict  (see post_process.calculate_derived_quantities)
        """
        if solver_config is None:
            solver_config = self._default_solver_config()

        sm = self.input_data.solver_mode
        print(f"\nLaunching solver_mode='{sm}' …")

        if sm in ('cpp_full', 'cpp_sliding_win', 'sliding_OpenMP'):
            results = self._run_cpp(solver_config)
        else:
            raise ValueError(f"Unknown solver_mode='{sm}'. "
                             "Use cpp_full, cpp_sliding_win, or sliding_OpenMP.")

        if results is not None and save_output:
            self._save_output(results, solver_config)

        return results

    # ── C++ solver dispatch ───────────────────────────────────────────────────

    def _run_cpp(self, solver_config):
        """Invoke the C++ solver via cpp_bridge."""
        from . import cpp_bridge
        results = cpp_bridge.run_cpp_solver(self, solver_config, base_dir=BASE_DIR)
        return results

    # ── Default solver config ─────────────────────────────────────────────────

    def _default_solver_config(self):
        re = self.input_data.reactions
        sm = self.input_data.solver_mode
        po = self.input_data.physics_option

        # Linear solver and window parameters from reactions sheet
        linsol = str(re.get('linsol', 'dense')).lower()
        w0_i   = int(float(re.get('window_w0_i',  100)))
        w_w    = int(float(re.get('window_width', 500)))
        C_exp  = float(re.get('window_C_exp',  1e-18))
        n_thr  = int(float(re.get('window_omp', 0)))

        # Map solver mode to window_mode integer
        window_mode_map = {
            'cpp_full':        0,
            'cpp_sliding_win': 3,
            'sliding_OpenMP':  4,
        }
        win_mode = window_mode_map.get(sm, 0)

        # For bin_moment, use gmres (larger state space often benefits)
        if 'bin_moment' in po:
            linsol = 'gmres'

        return {
            't_span':    (float(re.get('t_begin', 1e-8)),
                          float(re.get('t_end',   1e7))),
            'n_points':  int(float(re.get('n_points', 200))),
            'log_time':  bool(int(float(re.get('log_time', 1)))),
            'rtol':      float(re.get('rtol', 1e-8)),
            'atol':      float(re.get('atol', 1e-20)),
            'solver_method': {
                'backend':              'cvode',
                'lmm':                  'bdf',
                'linsol':               linsol,
                'window_mode':          win_mode,
                'window_w0_i':          w0_i,
                'window_width':         w_w,
                'window_C_expand':      C_exp,
                'window_expand_pad':    10,
                'window_omp_threads':   n_thr,
                'window_gmres_maxl':    20,
                'window_prec':          1,
            },
        }

    # ── Output writing ────────────────────────────────────────────────────────

    def _save_output(self, results, solver_config):
        """Write timestamped output directory."""
        try:
            git_hash = subprocess.check_output(
                ['git', 'rev-parse', '--short', 'HEAD'],
                cwd=str(BASE_DIR), stderr=subprocess.DEVNULL
            ).decode().strip()
        except Exception:
            git_hash = 'unknown'

        sm = self.input_data.solver_mode
        po = self.input_data.physics_option
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        label = f"{ts}_{sm}_{po}_{git_hash}"

        out_dir  = BASE_DIR / 'output' / label
        plot_dir = out_dir / 'plots'
        plot_dir.mkdir(parents=True, exist_ok=True)

        # Provenance
        d = self.input_data.derived
        with open(out_dir / 'provenance.md', 'w') as f:
            f.write(f"# Expanded_Eurofer_CD run\n\n")
            f.write(f"- timestamp:      {ts}\n")
            f.write(f"- git_hash:       {git_hash}\n")
            f.write(f"- solver_mode:    {sm}\n")
            f.write(f"- physics_option: {po}\n")
            f.write(f"- T:              {d['T']} K\n")
            f.write(f"- G:              {d['G']} dpa/s\n")
            f.write(f"- N:              {self.input_data.N}\n")
            f.write(f"- M:              {self.input_data.M}\n")
            f.write(f"- spectrum:       {d['spectrum']}\n")
            f.write(f"- t_span:         {solver_config['t_span']}\n")
            f.write(f"- rtol/atol:      {solver_config['rtol']} / {solver_config['atol']}\n")

        # Summary CSV
        import csv
        row = post_process.summary_csv_row(results, self.input_data,
                                           solver_label=f"{sm}/{po}")
        with open(out_dir / 'summary.csv', 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(row.keys()))
            writer.writeheader()
            writer.writerow(row)

        # Binary results (numpy)
        np.save(str(out_dir / 'results_t.npy'),  results['t'])
        np.save(str(out_dir / 'results_y.npy'),  results['y'])

        # Plots
        from . import visualization
        visualization.save_all_plots(results, self.input_data,
                                     str(plot_dir), label=f"{sm}/{po}")

        print(f"Output written to: {out_dir}")
        return out_dir
