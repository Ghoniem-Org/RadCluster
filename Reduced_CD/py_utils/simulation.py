# -*- coding: utf-8 -*-
"""
simulation.py — ClusterDynamics simulation orchestrator

Manages solver setup, time integration, result packaging, plotting,
and output saving (pkl + CSV).
"""

import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.integrate import solve_ivp

from py_utils.input_data import InputData
from py_utils.reaction_rates import ClusterRates
from py_utils.rate_equations import ClusterDynamicsODE


class ClusterDynamicsSimulation:
    """
    Orchestrates a single cluster-dynamics run:
      1. Load inputs
      2. Build rate coefficients and ODE system
      3. Integrate
      4. Post-process and save
    """

    def __init__(self, input_file=None):
        if input_file is None:
            input_file = Path(__file__).parent.parent / 'input' / 'input_parameters.xlsx'

        print("Loading input data...")
        self.input_data = InputData(input_file)

        print("Building rate coefficients...")
        self.rates = ClusterRates(self.input_data)

        print("Setting up ODE system...")
        self.ode = ClusterDynamicsODE(self.input_data, self.rates)

    # ------------------------------------------------------------------
    def _ode_wrapper(self, t, y):
        try:
            if np.any(np.isnan(y)) or np.any(np.isinf(y)):
                raise ValueError(f"Invalid state at t={t:.2e}")
            dydt = self.ode.ode_system(t, y)
            if np.any(np.isnan(dydt)) or np.any(np.isinf(dydt)):
                raise ValueError(f"Invalid derivatives at t={t:.2e}")
            return dydt
        except Exception as e:
            raise ValueError(f"ODE error at t={t:.2e}: {e}")

    # ------------------------------------------------------------------
    def run(self, t_begin=1e-1, t_end=1e6, n_points=500,
            method='LSODA', rtol=1e-6, atol=1e-20, log_time=True):
        """Integrate the ODE system and return a results dict."""

        if log_time:
            t_eval = np.logspace(np.log10(t_begin), np.log10(t_end), n_points)
        else:
            t_eval = np.linspace(t_begin, t_end, n_points)

        y0 = self.ode.get_initial_conditions()

        print(f"\nIntegrating: t=[{t_begin:.1e}, {t_end:.1e}] s  ({method})")
        sol = solve_ivp(
            self._ode_wrapper, [0, t_end], y0,
            method=method, t_eval=t_eval,
            rtol=rtol, atol=atol,
            dense_output=False,
            max_step=t_end / 100,
        )

        if not sol.success:
            raise RuntimeError(f"Integration failed: {sol.message}")

        print(f"Integration complete. nfev={sol.nfev}")
        return self._package_results(sol)

    # ------------------------------------------------------------------
    def _package_results(self, sol):
        t   = sol.t
        y   = np.maximum(sol.y, 0.0)
        G   = self.input_data.material_params['G']
        N_v = self.ode.N_v
        N_i = self.ode.N_i

        # Extract distributions
        fv = {n: y[self.ode.idx_fv + n - 1, :] for n in range(1, N_v + 1)}
        fi = {n: y[self.ode.idx_fi + n - 1, :] for n in range(1, N_i + 1)}

        # Loop radii at each time point
        r_111 = np.zeros(len(t))
        r_100 = np.zeros(len(t))
        for k in range(len(t)):
            r_111[k], r_100[k] = self.ode.loop_radii(y[:, k])

        results = {
            'time': t,
            'dpa':  t * G,
            'fv':   fv,
            'fi':   fi,
            'CiL_111':   y[self.ode.idx_iL_111,  :],
            'CiL_100':   y[self.ode.idx_iL_100,  :],
            'CiL_i_111': y[self.ode.idx_iLi_111, :],
            'CiL_i_100': y[self.ode.idx_iLi_100, :],
            'C_void':    y[self.ode.idx_void,     :],
            'r_void':    y[self.ode.idx_rvoid,    :],
            'r_iL_111':  r_111,
            'r_iL_100':  r_100,
            'metadata': {
                'T':    self.input_data.material_params['T'],
                'G':    G,
                'N_v':  N_v,
                'N_i':  N_i,
                'N_loop': self.ode.N_loop,
                'solver': {'success': sol.success, 'nfev': sol.nfev, 'message': sol.message},
            },
        }
        return results

    # ------------------------------------------------------------------
    def save(self, results, output_dir=None):
        """Save results as .pkl and summary .csv."""
        if output_dir is None:
            output_dir = Path.cwd() / 'output'
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        T  = int(self.input_data.material_params['T'] - 273)
        G  = self.input_data.material_params['G']
        tag = 'Ion' if G > 1e-5 else 'Neutron'
        stem = f'ClusterDynamics_{T}C_{tag}'

        pkl_file = output_dir / f'{stem}.pkl'
        with open(pkl_file, 'wb') as f:
            pickle.dump(results, f)
        print(f"Results saved → {pkl_file}")

        # Summary CSV
        csv_data = {
            'time':       results['time'],
            'dpa':        results['dpa'],
            'fv1':        results['fv'][1],
            'fi1':        results['fi'][1],
            'CiL_111':    results['CiL_111'],
            'CiL_100':    results['CiL_100'],
            'r_iL_111_nm': results['r_iL_111'] * 1e9,
            'r_iL_100_nm': results['r_iL_100'] * 1e9,
            'C_void':     results['C_void'],
            'r_void_nm':  results['r_void'] * 1e9,
        }
        csv_file = output_dir / f'{stem}_summary.csv'
        pd.DataFrame(csv_data).to_csv(csv_file, index=False)
        print(f"Summary saved  → {csv_file}")


# -----------------------------------------------------------------------
def run_cluster_dynamics_simulation(config, input_file=None):
    """
    Top-level entry point for notebooks.

    Parameters
    ----------
    config : dict
        Keys: t_begin, t_end, n_points, method, rtol, atol, log_time
    input_file : str or Path, optional

    Returns
    -------
    results : dict
    sim     : ClusterDynamicsSimulation
    """
    try:
        sim = ClusterDynamicsSimulation(input_file=input_file)
        results = sim.run(
            t_begin=config['t_begin'],
            t_end=config['t_end'],
            n_points=config['n_points'],
            method=config['method'],
            rtol=config['rtol'],
            atol=config['atol'],
            log_time=config['log_time'],
        )
        sim.save(results)
        return results, sim

    except Exception as e:
        import traceback
        print(f"Fatal error: {e}")
        traceback.print_exc()
        return None, None
