"""
simulation.py – ClusterDynamics simulation orchestrator.

Responsibilities:
  1. Load inputs and initialise physics modules (input_data, reaction_rates,
     rate_equations).
  2. Run the ODE integration with x_max-gated segmented stepping (Python solver).
  3. Delegate post-processing to post_process and plotting to visualization.

The segmented approach mimics CD.ipynb: the time span [t0, tf] is divided
into n_segments logarithmically spaced intervals.  Before each segment an
x_max estimate is computed (heuristic or NN-based) and only interstitial
cluster sizes ≤ x_max are integrated — sizes above are frozen at zero.

The ODE interface is kept as a thin, isolated layer to facilitate future
migration of the heavy-lifting to a C++ back-end (see cpp_bridge.py).
"""

import time as _time
import traceback

import numpy as np
from scipy.integrate import solve_ivp

from py_utils.input_data import InputData
from py_utils.reaction_rates import ReactionRates
from py_utils.rate_equations import RateEquations
from py_utils import pre_process, post_process


class ClusterDynamicsSimulation:
    """
    Orchestrates the ClusterDynamics simulation workflow:
    load inputs → initialise physics → run ODE → post-process → plot.
    """

    def __init__(self, Nv=None, Ni=None):
        print("Initializing ClusterDynamics simulation…")
        try:
            self.input_data      = InputData(Nv=Nv, Ni=Ni)
            self.reaction_rates  = ReactionRates(self.input_data)
            self.rate_equations  = RateEquations(self.input_data, self.reaction_rates)
            pre_process.validate_setup(self.input_data)
            print("✓ Simulation initialized successfully!")
        except Exception as e:
            print(f"❌ Error during initialization: {e}")
            raise

    # ── Python ODE solver (segmented, x_max-gated) ───────────────────────────

    def run_simulation(self, t_span=(1e-6, 1e5), n_segments=60,
                       rtol=1e-8, atol=1e-50, verbose=True):
        """
        Integrate using LSODA with x_max-gated segmented stepping.

        Parameters
        ----------
        t_span      : (t0, tf) integration window [s]
        n_segments  : number of log-spaced time segments
        rtol, atol  : ODE solver tolerances
        verbose     : print per-segment progress

        Returns
        -------
        results : dict or None
        """
        Ni   = self.rate_equations.Ni
        Nv   = self.rate_equations.Nv
        t_edges = np.logspace(np.log10(t_span[0]), np.log10(t_span[1]), n_segments + 1)

        y0      = pre_process.get_initial_conditions(self.rate_equations)
        y_cur   = y0.copy()

        t_out       = [t_edges[0]]
        y_out       = [y0.copy()]
        xmax_history = []

        x_max_cur = 3      # start with only 3 interstitial sizes active
        wall_t0   = _time.perf_counter()

        print(f"\nRunning Python LSODA solver  ({n_segments} segments,"
              f"  t ∈ [{t_span[0]:.1e}, {t_span[1]:.1e}] s)")

        for seg in range(n_segments):
            t0_seg, t1_seg = t_edges[seg], t_edges[seg + 1]

            # ── Update x_max estimate ────────────────────────────────────────
            x_max_next = self._predict_xmax(
                t0_seg, float(y_cur[0]), float(y_cur[Nv]),
                x_max_cur, Ni, margin=3,
            )
            xmax_history.append(x_max_next)

            # ── Integrate this segment ────────────────────────────────────────
            gated = self._make_gated_rhs(x_max_next)
            sol   = solve_ivp(
                gated, (t0_seg, t1_seg), y_cur,
                method='LSODA', rtol=rtol, atol=atol,
            )

            if not sol.success and verbose:
                print(f"  [seg {seg}] ⚠️  {sol.message}")

            y_cur     = np.maximum(sol.y[:, -1], 0.0)
            x_max_cur = x_max_next

            t_out.append(t1_seg)
            y_out.append(y_cur.copy())

            if verbose and (seg % 10 == 0 or seg == n_segments - 1):
                pct = 100 * x_max_next / Ni
                print(f"  seg {seg:3d}  t={t1_seg:.2e}s  x_max={x_max_next:3d}"
                      f" ({pct:.0f}% of Ni)  Ci={float(y_cur[Nv]):.2e}")

        elapsed = _time.perf_counter() - wall_t0
        print(f"\n✓ Python solver finished  —  wall time: {elapsed:.2f} s")

        t_arr = np.array(t_out)              # (n_segments+1,)
        y_arr = np.column_stack(y_out)       # (N, n_segments+1)

        results = post_process.calculate_derived_quantities(
            t_arr, y_arr, self.input_data, self.rate_equations,
            xmax_history=xmax_history,
        )
        results['metadata'] = {
            'solver_stats': {
                'success':  True,
                'message':  f'Python LSODA segmented ({n_segments} segments)',
                'wall_time': elapsed,
                'n_segments': n_segments,
                'n_time_points': len(t_arr),
            },
        }
        return results

    # ── x_max gating helpers ─────────────────────────────────────────────────

    @staticmethod
    def _predict_xmax(t_now, Cv, Ci, x_max_now, Ni, margin=3):
        """
        Heuristic x_max estimate: linear growth on a log-time axis.

        Mirrors the heuristic branch in CD.ipynb (no NN required).
        x_max never shrinks; capped at Ni.
        """
        t_frac    = np.clip(np.log10(max(t_now, 1e-6) / 1e-6) / 11.0, 0.0, 1.0)
        x_max_new = int(3 + (Ni - 3) * t_frac) + margin
        return min(max(x_max_new, x_max_now), Ni)

    def _make_gated_rhs(self, x_max):
        """Return an ODE callable that freezes interstitial sizes > x_max."""
        re = self.rate_equations
        freeze_hi = slice(re.Nv + x_max, re.N)

        def gated(t, y, fh=freeze_hi):
            dydt = re._rhs_full(y)
            dydt[fh] = 0.0
            return dydt

        return gated
