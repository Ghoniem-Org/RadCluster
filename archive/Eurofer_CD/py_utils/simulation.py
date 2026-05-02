"""
simulation.py — Eurofer_CD simulation orchestrator.

Responsibilities:
  1. Load inputs and initialise physics modules (input_data, reaction_rates,
     rate_equations).
  2. Run the ODE integration with x_max-gated segmented stepping.
  3. Delegate post-processing to post_process and plotting to visualization.

Adapted from Full_CD/py_utils/simulation.py for the Eurofer_CD state vector
(Ni SIA + Nv vacancy + 1 He species).
"""

import time as _time

import numpy as np
from scipy.integrate import solve_ivp

from .input_data     import InputData
from .reaction_rates import ReactionRates
from .rate_equations import RateEquations
from .               import post_process


class EuroferCDSimulation:
    """
    Orchestrates the Eurofer_CD simulation workflow:
    load inputs → initialise physics → run ODE → post-process → (optionally) plot.

    Parameters
    ----------
    Nv      : int, optional  — override max vacancy cluster size
    Ni      : int, optional  — override max interstitial cluster size
    he_mode : str, optional  — override He-vacancy reduction mode
    excel_file : path-like, optional  — override Excel input path
    """

    def __init__(self, Nv=None, Ni=None, he_mode=None, excel_file=None):
        print("Initializing Eurofer_CD simulation…")
        kwargs = {}
        if excel_file is not None:
            kwargs['excel_file'] = excel_file
        try:
            self.input_data     = InputData(Nv=Nv, Ni=Ni, he_mode=he_mode, **kwargs)
            self.reaction_rates = ReactionRates(self.input_data)
            self.rate_equations = RateEquations(self.input_data, self.reaction_rates)
            print("Simulation initialized successfully.")
        except Exception as e:
            print(f"Error during initialization: {e}")
            raise

    # ── Python ODE solver ────────────────────────────────────────────────────

    def run_simulation(self, t_span=None, n_segments=None,
                       rtol=None, atol=None, verbose=True):
        """
        Integrate using LSODA with x_max-gated segmented stepping.

        Parameters read from Model_Parameters sheet by default; keyword
        arguments override them.

        Parameters
        ----------
        t_span     : (t0, tf), optional
        n_segments : int, optional
        rtol       : float, optional
        atol       : float, optional
        verbose    : bool

        Returns
        -------
        results : dict  (see post_process.calculate_derived_quantities)
        """
        mp = self.input_data.model_params

        t0  = float(mp.get('t_begin', 1e-8))
        tf  = float(mp.get('t_end',   1e6))
        if t_span is not None:
            t0, tf = t_span
        nseg = int(mp.get('n_segments', 60)) if n_segments is None else int(n_segments)
        _rtol = float(mp.get('rtol', 1e-8))  if rtol is None else float(rtol)
        # atol=1e-50 was too tight: late-time quasi-steady-state derivatives are
        # small differences of large terms; LSODA cannot achieve 1e-50 accuracy
        # and hits convergence failures.  1e-16 is the practical floor for
        # double-precision ODE integration.
        _atol = float(mp.get('atol', 1e-16)) if atol is None else float(atol)

        Ni = self.rate_equations.Ni
        re = self.rate_equations

        t_edges = np.logspace(np.log10(t0), np.log10(tf), nseg + 1)

        y0    = re.get_initial_conditions()
        y_cur = y0.copy()

        t_out        = [t_edges[0]]
        y_out        = [y0.copy()]
        xmax_history = []
        x_max_cur    = 3

        wall_t0 = _time.perf_counter()

        print(f"\nRunning LSODA solver  ({nseg} segments,"
              f"  t=[{t0:.1e}, {tf:.1e}] s,"
              f"  he_mode='{re.he_mode}')")

        for seg in range(nseg):
            t0_seg, t1_seg = t_edges[seg], t_edges[seg + 1]

            # Update x_max gating estimate
            x_max_next = self._predict_xmax(
                t0_seg, float(y_cur[re.i_SIA]),
                float(y_cur[re.i_VAC]), x_max_cur, Ni, margin=3,
            )
            xmax_history.append(x_max_next)

            gated = self._make_gated_rhs(x_max_next)
            seg_span = t1_seg - t0_seg
            sol   = solve_ivp(
                gated, (t0_seg, t1_seg), y_cur,
                method='LSODA', rtol=_rtol, atol=_atol,
                max_step=seg_span,   # prevent LSODA jumping past segment boundary
            )

            if not sol.success and verbose:
                print(f"  [seg {seg}] {sol.message}")

            y_cur     = np.maximum(sol.y[:, -1], 0.0)
            x_max_cur = x_max_next

            t_out.append(t1_seg)
            y_out.append(y_cur.copy())

            if verbose:
                Ci1 = float(y_cur[re.i_SIA])
                Cv1 = float(y_cur[re.i_VAC])
                pct = 100 * x_max_next / Ni
                print(f"  seg {seg:3d}  t={t1_seg:.2e}s  x_max={x_max_next:3d}"
                      f" ({pct:.0f}%)  Ci1={Ci1:.2e}  Cv1={Cv1:.2e}")

        elapsed = _time.perf_counter() - wall_t0
        print(f"\nSolver finished  —  wall time: {elapsed:.2f} s")

        t_arr = np.array(t_out)
        y_arr = np.column_stack(y_out)

        results = post_process.calculate_derived_quantities(
            t_arr, y_arr, self.input_data, self.rate_equations,
            xmax_history=xmax_history,
        )
        results['metadata'] = {
            'solver_stats': {
                'success':       True,
                'message':       f'LSODA segmented ({nseg} segments)',
                'wall_time':     elapsed,
                'n_segments':    nseg,
                'n_time_points': len(t_arr),
                'he_mode':       re.he_mode,
            },
        }
        return results

    # ── x_max gating helpers ─────────────────────────────────────────────────

    @staticmethod
    def _predict_xmax(t_now, Ci1, Cv1, x_max_now, Ni, margin=3):
        """
        Heuristic x_max estimate: linear growth on a log-time axis.
        x_max never shrinks; capped at Ni.
        """
        t_frac    = np.clip(np.log10(max(t_now, 1e-8) / 1e-8) / 14.0, 0.0, 1.0)
        x_max_new = int(3 + (Ni - 3) * t_frac) + margin
        return min(max(x_max_new, x_max_now), Ni)

    def _make_gated_rhs(self, x_max):
        """Return an ODE callable that freezes SIA cluster sizes > x_max."""
        re     = self.rate_equations
        freeze = slice(re.i_SIA + x_max, re.i_VAC)  # freeze SIA clusters > x_max

        def gated(t, y, fh=freeze):
            dydt = re.ode_system(t, y)
            dydt[fh] = 0.0
            return dydt

        return gated
