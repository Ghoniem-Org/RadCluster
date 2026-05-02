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
    I              : int, optional  — override max SIA cluster size
    V              : int, optional  — override max vacancy cluster size
    solver_mode    : str, optional  — 'cpp_full' | 'cpp_sliding_win' | 'sliding_OpenMP'
    physics_option : str, optional  — 'full_CD_fission' | 'full_CD_fusion' |
                                      'bin_moment_CD_fission' | 'bin_moment_CD_fusion'
    excel_file     : path-like, optional
    C_floor        : float, optional — concentration floor (default 1e-15).
                                       Any state variable below this is clamped
                                       before evaluating rate terms.
    he_options     : str, optional   — 'dynamic' (default) integrates free He as a
                                       full ODE (Eq. 157); 'quasi_steady_state'
                                       eliminates it via dc_h/dt = 0 (valid because
                                       E_m_h = 0.06 eV gives rapid equilibration).
    i_mobile       : int, optional   — max mobile SIA cluster size (default 1).
                                       1 = mono-defect form; set > 1 to enable 1D
                                       glide and SIA–SIA coalescence.
    v_mobile       : int, optional   — max mobile vacancy cluster size (default 1).
                                       1 = mono-defect form; set > 1 to enable
                                       vacancy–vacancy coalescence.
    """

    def __init__(self, I=None, V=None, solver_mode=None,
                 physics_option=None, excel_file=None,
                 C_floor=None, he_options=None,
                 i_mobile=None, v_mobile=None,
                 **legacy_kw):
        # Backward compatibility: accept old N/M and n_max_i/m_max_v kwargs
        if I is None and 'N' in legacy_kw:
            I = legacy_kw.pop('N')
        if V is None and 'M' in legacy_kw:
            V = legacy_kw.pop('M')
        if i_mobile is None and 'n_max_i' in legacy_kw:
            i_mobile = legacy_kw.pop('n_max_i')
        if v_mobile is None and 'm_max_v' in legacy_kw:
            v_mobile = legacy_kw.pop('m_max_v')

        print("Initializing Expanded_Eurofer_CD simulation…")
        kwargs = {}
        if excel_file is not None:
            kwargs['excel_file'] = excel_file

        self.input_data = InputData(
            I=I, V=V,
            solver_mode=solver_mode,
            physics_option=physics_option,
            **kwargs
        )

        # Inject optional overrides before building rate equations
        if i_mobile is not None:
            self.input_data.diffusion['i_mobile'] = int(i_mobile)
            self.input_data.derived['i_mobile']   = int(i_mobile)
        if v_mobile is not None:
            self.input_data.diffusion['v_mobile'] = int(v_mobile)
            self.input_data.derived['v_mobile']   = int(v_mobile)
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

        # Always defined so notebook fallbacks like
        # `sim._accumulated_results` don't AttributeError when run_adaptive
        # is interrupted before it can initialize this attribute itself.
        self._accumulated_results = None

        print(f"Simulation initialized: solver_mode='{self.input_data.solver_mode}'"
              f"  physics_option='{po}'")

    def rebuild_rates(self):
        """Rebuild reaction rates and rate equations from current input_data.

        Call this after modifying input_data parameter dicts and calling
        input_data._calculate_derived() to ensure the solver uses the
        updated rate constants.
        """
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

    def _resize_domain(self, I_new, V_new):
        """Resize the cluster domain and rebuild all rate objects."""
        self.input_data.reactions['I'] = int(I_new)
        self.input_data.reactions['V'] = int(V_new)
        # InputData caches I, V as properties from reactions dict
        self.rebuild_rates()

    def _boundary_fraction_at(self, results, idx):
        """
        Compute tail fractions at a specific time index.

        Returns (frac_I, frac_V) — fraction of the total SIA/vacancy content
        in the last 10% of sizes.  Values > ~0.01 indicate boundary pileup.
        """
        y = results.get('y')
        if y is None:
            return 0.0, 0.0

        re = self.rate_equations
        is_bin = hasattr(re, 'bins')
        I, V = self.input_data.I, self.input_data.V
        ns = np.arange(1, I + 1, dtype=float)
        ms = np.arange(1, V + 1, dtype=float)

        yj = np.maximum(y[:, idx], 0.0)

        if is_bin:
            from .bin_moment_rates import reconstruct_distribution
            I_bin = getattr(re, 'I_bin', getattr(re, 'K', 0))
            V_bin = getattr(re, 'V_bin', getattr(re, 'K_v', 0))
            P     = getattr(re, 'n_mom', 2)
            sf    = getattr(re, 'shape_function', 'linear')
            i_d   = getattr(re, 'i_discrete', 0)
            v_d   = getattr(re, 'v_discrete', 0)
            iv    = getattr(re, 'i_VAC', i_d + P * I_bin)

            # SIA: discrete + binned reconstruction
            if I_bin > 0:
                mom = yj[i_d:i_d + P * I_bin]
                mu0 = mom[0::P][:I_bin]
                mu1 = mom[1::P][:I_bin] if P >= 2 else None
                mu2 = mom[2::P][:I_bin] if P >= 3 else None
                c_n = reconstruct_distribution(sf, mu0, mu1, mu2, re.bins, I)
                c_n[:i_d] = yj[:i_d]
            else:
                c_n = np.zeros(I)
                c_n[:i_d] = yj[:i_d]

            # VAC: discrete + binned reconstruction
            if V_bin > 0:
                vac_start = iv + v_d
                vmom = yj[vac_start:vac_start + P * V_bin]
                vmu0 = vmom[0::P][:V_bin]
                vmu1 = vmom[1::P][:V_bin] if P >= 2 else None
                vmu2 = vmom[2::P][:V_bin] if P >= 3 else None
                c_v = reconstruct_distribution(sf, vmu0, vmu1, vmu2,
                                               re.vac_bins, V)
                c_v[:v_d] = yj[iv:iv + v_d]
            else:
                # All discrete vacancies: only v_discrete entries stored
                c_v = np.zeros(V)
                c_v[:v_d] = yj[iv:iv + v_d]
        else:
            c_n = yj[:I]
            c_v = yj[I:I + V]

        # Fraction of content in last 10% of sizes
        tail_n = int(max(I * 0.1, 1))
        tail_m = int(max(V * 0.1, 1))
        total_sia = np.dot(ns, np.maximum(c_n, 0.0))
        total_vac = np.dot(ms, np.maximum(c_v, 0.0))
        tail_sia  = np.dot(ns[-tail_n:], np.maximum(c_n[-tail_n:], 0.0))
        tail_vac  = np.dot(ms[-tail_m:], np.maximum(c_v[-tail_m:], 0.0))

        # Guard: if total content is at floor level (all concentrations ≈ C_floor),
        # the tail fraction is a floor artifact, not real boundary pileup.
        # Require total content to be above I * C_floor * I (roughly I^2 * C_floor)
        # to distinguish real content from uniform-floor noise.
        C_floor = getattr(self.input_data, 'C_floor', 1e-20)
        sia_threshold = I * I * C_floor * 10.0
        vac_threshold = V * V * C_floor * 10.0
        frac_I = tail_sia / max(total_sia, 1e-300) if total_sia > sia_threshold else 0.0
        frac_V = tail_vac / max(total_vac, 1e-300) if total_vac > vac_threshold else 0.0
        return frac_I, frac_V

    def _boundary_fraction(self, results):
        """Check tail fractions at the last time step (backward compat)."""
        y = results.get('y')
        if y is None:
            return 0.0, 0.0
        return self._boundary_fraction_at(results, -1)

    def _find_first_exceedance(self, results, threshold):
        """
        Scan all output time points and return the index of the first one
        where either the SIA or vacancy tail fraction exceeds *threshold*.

        Returns None if no exceedance is found.
        """
        y = results.get('y')
        if y is None:
            return None
        n_t = y.shape[1]
        for j in range(n_t):
            frac_I, frac_V = self._boundary_fraction_at(results, j)
            if frac_I > threshold or frac_V > threshold:
                return j
        return None

    def _expand_state(self, y_old, I_new, V_new):
        """
        Map ODE state from the current domain to a larger (I_new, V_new) domain.

        1. Reconstruct per-size distributions from old state.
        2. Pad with C_floor for new sizes beyond old I, V.
        3. Resize domain (rebuilds bins and rate equations).
        4. Project padded distributions onto new bins.

        Returns y0_new for the enlarged domain.
        """
        re_old = self.rate_equations
        is_bin = hasattr(re_old, 'bins')
        I_old  = self.input_data.I
        V_old  = self.input_data.V
        C_floor = float(self.input_data.reactions.get('C_floor', 1e-25))

        yj = np.maximum(y_old, 0.0)

        if is_bin:
            from .bin_moment_rates import (reconstruct_distribution,
                                           moments_from_distribution)
            K_i_old = re_old.K_i
            K_v_old = re_old.K_v
            iv_old  = re_old.i_VAC
            P_old   = getattr(re_old, 'n_mom', 2)
            sf_old  = getattr(re_old, 'shape_function', 'linear')
            i_d_old = getattr(re_old, 'i_discrete', 0)
            v_d_old = getattr(re_old, 'v_discrete', 0)

            # Reconstruct SIA per-size from moments
            sia_mom = yj[i_d_old:i_d_old + P_old * K_i_old]
            sia_mu0 = sia_mom[0::P_old][:K_i_old]
            sia_mu1 = sia_mom[1::P_old][:K_i_old] if P_old >= 2 else None
            sia_mu2 = sia_mom[2::P_old][:K_i_old] if P_old >= 3 else None
            c_n_old = reconstruct_distribution(sf_old, sia_mu0, sia_mu1,
                                               sia_mu2, re_old.bins, I_old)
            c_n_old[:i_d_old] = yj[:i_d_old]

            # Reconstruct vacancy per-size from moments
            vac_mom = yj[iv_old + v_d_old:iv_old + v_d_old + P_old * K_v_old]
            vmu0 = vac_mom[0::P_old][:K_v_old]
            vmu1 = vac_mom[1::P_old][:K_v_old] if P_old >= 2 else None
            vmu2 = vac_mom[2::P_old][:K_v_old] if P_old >= 3 else None
            c_v_old = reconstruct_distribution(sf_old, vmu0, vmu1, vmu2,
                                               re_old.vac_bins, V_old)
            c_v_old[:v_d_old] = yj[iv_old:iv_old + v_d_old]

            # Extract He state
            if re_old.he_mode == 'case2':
                Q_tot = yj[re_old.i_Qtot]
                Q_k_old = None
            else:
                Q_k_old = yj[re_old.i_Q:re_old.i_Q + K_v_old]
                Q_tot = np.sum(Q_k_old)

            c_h = None if re_old.qss_He else yj[re_old.i_He]

            # Extend distributions to new domain (pad with C_floor)
            c_n_new = np.full(I_new, C_floor)
            c_n_new[:I_old] = c_n_old
            c_v_new = np.full(V_new, C_floor)
            c_v_new[:V_old] = c_v_old

            # Resize domain — rebuilds bins and rate equations
            self._resize_domain(I_new, V_new)
            re_new = self.rate_equations
            P_new = re_new.n_mom

            # Project onto new bins
            sia_mu0_new, sia_mu1_new, sia_mu2_new = moments_from_distribution(
                c_n_new, re_new.bins, n_mom=P_new)
            vac_mu0_new, vac_mu1_new, vac_mu2_new = moments_from_distribution(
                c_v_new, re_new.vac_bins, n_mom=P_new)

            # Build new state vector
            y0 = np.full(re_new.N_eq, C_floor)
            # Discrete SIA
            i_d_new = re_new.i_discrete
            y0[:i_d_new] = c_n_new[:i_d_new]
            # Binned SIA moments
            for k in range(re_new.K_i):
                y0[i_d_new + P_new * k] = sia_mu0_new[k]
                if P_new >= 2 and sia_mu1_new is not None:
                    y0[i_d_new + P_new * k + 1] = sia_mu1_new[k]
                if P_new >= 3 and sia_mu2_new is not None:
                    y0[i_d_new + P_new * k + 2] = sia_mu2_new[k]
            # Discrete VAC
            iv_new = re_new.i_VAC
            v_d_new = re_new.v_discrete
            y0[iv_new:iv_new + v_d_new] = c_v_new[:v_d_new]
            # Binned VAC moments
            vac_start_new = iv_new + v_d_new
            for k in range(re_new.K_v):
                y0[vac_start_new + P_new * k] = vac_mu0_new[k]
                if P_new >= 2 and vac_mu1_new is not None:
                    y0[vac_start_new + P_new * k + 1] = vac_mu1_new[k]
                if P_new >= 3 and vac_mu2_new is not None:
                    y0[vac_start_new + P_new * k + 2] = vac_mu2_new[k]

            # He state
            if re_new.he_mode == 'case2':
                y0[re_new.i_Qtot] = Q_tot
            else:
                # Distribute Q_tot proportionally among new vacancy bins
                vac_total = np.sum(vac_mu0_new)
                if vac_total > 0:
                    for k in range(re_new.K_v):
                        y0[re_new.i_Q + k] = Q_tot * vac_mu0_new[k] / vac_total
            if not re_new.qss_He and c_h is not None:
                y0[re_new.i_He] = c_h

            return y0

        else:
            # Full per-size mode — pad state vector directly
            c_n_old = yj[:I_old]
            c_v_old = yj[I_old:I_old + V_old]

            if re_old.he_mode == 'case2':
                Q_tot = yj[re_old.i_Qtot]
                Q_m_old = None
            else:
                Q_m_old = yj[re_old.i_Q:re_old.i_Q + V_old]
                Q_tot = np.sum(Q_m_old)
            c_h = None if re_old.qss_He else yj[re_old.i_He]

            # Resize domain
            self._resize_domain(I_new, V_new)
            re_new = self.rate_equations

            y0 = np.full(re_new.N_eq, C_floor)
            y0[:I_old] = c_n_old
            y0[I_new:I_new + V_old] = c_v_old

            if re_new.he_mode == 'case2':
                y0[re_new.i_Qtot] = Q_tot
            else:
                y0[re_new.i_Q:re_new.i_Q + V_old] = Q_m_old
            if not re_new.qss_He and c_h is not None:
                y0[re_new.i_He] = c_h

            return y0

    # ── Time-series keys that should be concatenated when merging segments ────

    _TS_KEYS = {'t', 'dose', 'C_SIA_tot', 'C_VAC_tot', 'C_He_tot',
                'C_He_free', 'mean_n_i', 'mean_n_v', 'N_loops', 'N_voids',
                'swelling', 'C_i1', 'C_v1', 'delta_FP', 'delta_He',
                'J_SIA_fixed', 'J_SIA_mutual', 'J_VAC_fixed', 'J_VAC_mutual',
                'J_He_sink'}

    def _slice_results(self, results, start, end):
        """Slice a results dict to time indices [start:end)."""
        sliced = {}
        n_t = len(results['t'])
        for key, val in results.items():
            if key in self._TS_KEYS and isinstance(val, np.ndarray):
                sliced[key] = val[start:end]
            elif key == 'y' and isinstance(val, np.ndarray) and val.ndim == 2:
                sliced[key] = val[:, start:end]
            else:
                sliced[key] = val
        return sliced

    def _merge_results(self, accumulated, new_segment):
        """
        Merge two results dicts by concatenating time-series arrays.

        The first time point of *new_segment* is skipped if it duplicates
        the last time point of *accumulated* (continuation overlap).
        The raw state ``y`` is decomposed into SIA / VAC / He blocks,
        each zero-padded independently to the larger layout, then
        reassembled and concatenated so that plots cover the full dose
        range even after domain-doubling.
        """
        if accumulated is None:
            return new_segment
        if new_segment is None:
            return accumulated

        # Detect overlap: first point of new matches last of accumulated
        skip = 0
        if (len(accumulated['t']) > 0 and len(new_segment['t']) > 0 and
                abs(new_segment['t'][0] - accumulated['t'][-1]) /
                max(abs(accumulated['t'][-1]), 1e-30) < 1e-8):
            skip = 1

        merged = {}
        for key in new_segment:
            old_val = accumulated.get(key)
            new_val = new_segment[key]
            if key in self._TS_KEYS and isinstance(new_val, np.ndarray):
                if old_val is not None and isinstance(old_val, np.ndarray):
                    merged[key] = np.concatenate([old_val, new_val[skip:]])
                else:
                    merged[key] = new_val
            elif key == 'y':
                if old_val is not None and isinstance(old_val, np.ndarray):
                    merged[key] = self._merge_y_blocks(
                        accumulated, new_segment, skip)
                else:
                    merged[key] = new_val
            else:
                merged[key] = new_val
        return merged

    @staticmethod
    def _merge_y_blocks(accumulated, new_segment, skip):
        """
        Remap old y to the new layout and concatenate.

        Uses ``_y_i_VAC`` and ``_y_i_He`` stored in each results dict to
        decompose y into SIA / VAC / He blocks, zero-pad each to the
        larger size, then reassemble into a single array with the new
        layout.
        """
        old_y  = accumulated['y']
        new_y  = new_segment['y']

        iv_old = accumulated.get('_y_i_VAC', old_y.shape[0])
        ih_old = accumulated.get('_y_i_He',  old_y.shape[0])
        iv_new = new_segment.get('_y_i_VAC', new_y.shape[0])
        ih_new = new_segment.get('_y_i_He',  new_y.shape[0])

        # If layout indices aren't stored, fall back to keeping latest only
        if '_y_i_VAC' not in accumulated or '_y_i_VAC' not in new_segment:
            return new_y

        # Decompose into blocks
        sia_old = old_y[:iv_old, :]
        vac_old = old_y[iv_old:ih_old, :]
        he_old  = old_y[ih_old:, :]

        sia_new = new_y[:iv_new, :]
        vac_new = new_y[iv_new:ih_new, :]
        he_new  = new_y[ih_new:, :]

        # Zero-pad each block of old to match new sizes
        def _pad_rows(arr, target_rows):
            if arr.shape[0] >= target_rows:
                return arr[:target_rows, :]
            pad = np.zeros((target_rows - arr.shape[0], arr.shape[1]))
            return np.vstack([arr, pad])

        sia_old = _pad_rows(sia_old, iv_new)
        vac_old = _pad_rows(vac_old, ih_new - iv_new)
        he_n    = he_new.shape[0]
        he_old  = _pad_rows(he_old, he_n)

        # Reassemble old y in new layout
        old_remapped = np.vstack([sia_old, vac_old, he_old])
        return np.concatenate([old_remapped, new_y[:, skip:]], axis=1)

    def run_adaptive(self, solver_config=None, save_output=True,
                     progress_callback=None, timeout_s=None,
                     boundary_threshold=0.01, max_doublings=6,
                     points_per_segment=10):
        """
        Truly adaptive domain doubling via short time-segments.

        The total time span is divided into short segments of
        *points_per_segment* output points each.  After every segment
        the tail fraction is checked.  If it exceeds *boundary_threshold*:

        1. The first exceedance point within the segment is located.
        2. Results up to that point are kept.
        3. I and/or V are doubled; the state vector is mapped to the
           larger domain via ``_expand_state``.
        4. Integration resumes from the exceedance time.

        This allows starting with a small domain (e.g. I = V = 1000) and
        growing it on demand as the distribution spreads, without ever
        restarting from t = 0.

        Parameters
        ----------
        boundary_threshold : float
            Adapt if tail fraction > this (default 0.01 = 1 %).
        max_doublings : int
            Maximum number of doublings (default 6; 1000 → 64 000).
        points_per_segment : int
            Output points per solver call (default 10).  Smaller values
            detect exceedance earlier but add subprocess overhead.

        Returns
        -------
        results : dict
        """
        if solver_config is None:
            solver_config = self._default_solver_config()

        t_begin, t_end = solver_config['t_span']
        n_points_total = solver_config.get('n_points', 100)
        log_time       = solver_config.get('log_time', True)

        # Precompute checkpoint times (segment boundaries)
        n_segments = max(1, n_points_total // points_per_segment)
        if log_time and t_begin > 0:
            checkpoints = np.geomspace(t_begin, t_end, n_segments + 1)
        else:
            checkpoints = np.linspace(t_begin, t_end, n_segments + 1)

        accumulated  = None
        self._accumulated_results = None   # expose for graceful interrupt
        current_t    = t_begin
        y0_override  = None
        n_doublings_I = 0     # independent doubling counters
        n_doublings_V = 0
        seg_count    = 0
        interrupted  = False
        # Allow extra iterations for expansion restarts within segments
        max_iters    = n_segments + 2 * max_doublings + 1

        try:
          while current_t < t_end * (1 - 1e-10) and seg_count < max_iters:
            seg_count += 1
            I_cur = self.input_data.I
            V_cur = self.input_data.V

            # Next checkpoint after current_t
            future = checkpoints[checkpoints > current_t * (1 + 1e-10)]
            seg_t_end = float(future[0]) if len(future) > 0 else t_end

            print(f"\n[Segment {seg_count}] I={I_cur}  V={V_cur}"
                  f"  t=[{current_t:.2e}, {seg_t_end:.2e}]")

            seg_config = dict(solver_config)
            seg_config['t_span']   = (current_t, seg_t_end)
            seg_config['n_points'] = points_per_segment

            results = self.run(
                solver_config=seg_config,
                save_output=False,
                progress_callback=progress_callback,
                timeout_s=timeout_s,
                y0_override=y0_override,
            )

            if results is None:
                print("  Solver failed — returning accumulated results.")
                if accumulated is not None and save_output:
                    self._diag_text = self.reaction_rates.format_diagnostic(
                        mean_n_i=accumulated['mean_n_i'][-1]
                        if 'mean_n_i' in accumulated else None)
                    self._save_output(accumulated, solver_config)
                return accumulated

            # Check tail at the last point of this segment
            frac_I, frac_V = self._boundary_fraction_at(results, -1)
            can_double_I = frac_I > boundary_threshold and n_doublings_I < max_doublings
            can_double_V = frac_V > boundary_threshold and n_doublings_V < max_doublings

            if can_double_I or can_double_V:
                # Find the first exceedance point within this segment
                exceed_idx = self._find_first_exceedance(
                    results, boundary_threshold)
                if exceed_idx is None:
                    exceed_idx = results['y'].shape[1] - 1

                frac_I_ex, frac_V_ex = self._boundary_fraction_at(
                    results, exceed_idx)
                t_exceed = results['t'][exceed_idx]

                # Keep results up to and including exceedance
                partial = self._slice_results(results, 0, exceed_idx + 1)
                accumulated = self._merge_results(accumulated, partial)
                self._accumulated_results = accumulated

                # Double each dimension independently if it exceeds AND
                # still has budget remaining
                do_I = frac_I_ex > boundary_threshold and n_doublings_I < max_doublings
                do_V = frac_V_ex > boundary_threshold and n_doublings_V < max_doublings
                I_new = I_cur * 2 if do_I else I_cur
                V_new = V_cur * 2 if do_V else V_cur
                which = []
                if do_I:
                    which.append(f'I: {I_cur}→{I_new} [{n_doublings_I+1}/{max_doublings}]')
                    n_doublings_I += 1
                if do_V:
                    which.append(f'V: {V_cur}→{V_new} [{n_doublings_V+1}/{max_doublings}]')
                    n_doublings_V += 1
                print(f"  Exceedance at t={t_exceed:.3e}: "
                      f"SIA={frac_I_ex:.3f}  VAC={frac_V_ex:.3f}")
                print(f"  Doubling: {', '.join(which)}")

                # Map state to enlarged domain
                y_at = results['y'][:, exceed_idx]
                y0_override = self._expand_state(y_at, I_new, V_new)
                current_t   = t_exceed
            else:
                # Segment OK (or both dimensions maxed out) — accumulate
                exceeded_but_maxed = []
                if frac_I > boundary_threshold and n_doublings_I >= max_doublings:
                    exceeded_but_maxed.append(f'SIA={frac_I:.3f}')
                if frac_V > boundary_threshold and n_doublings_V >= max_doublings:
                    exceeded_but_maxed.append(f'VAC={frac_V:.3f}')
                if exceeded_but_maxed:
                    print(f"  Tail: {', '.join(exceeded_but_maxed)} "
                          f"(max doublings reached)")
                else:
                    print(f"  Tail OK: SIA={frac_I:.3f}  VAC={frac_V:.3f}")

                accumulated = self._merge_results(accumulated, results)
                self._accumulated_results = accumulated
                current_t   = results['t'][-1]
                y0_override = results['y'][:, -1]

        except KeyboardInterrupt:
            interrupted = True
            print("\n\n*** KeyboardInterrupt — stopping gracefully. ***")
            # Merge the last completed segment if it wasn't merged yet
            if results is not None and results is not accumulated:
                try:
                    accumulated = self._merge_results(accumulated, results)
                    self._accumulated_results = accumulated
                except Exception:
                    pass   # keep whatever was accumulated before

        # ── Finished (or interrupted) — save and return ───────────────────
        if accumulated is not None:
            I_final = self.input_data.I
            V_final = self.input_data.V
            tot = n_doublings_I + n_doublings_V
            status = "INTERRUPTED" if interrupted else "complete"
            n_pts  = len(accumulated['t'])
            print(f"\nAdaptive run {status}: {seg_count} segments, "
                  f"{tot} doublings (I×{n_doublings_I}, V×{n_doublings_V}), "
                  f"final domain I={I_final} V={V_final}, "
                  f"{n_pts} time points saved")
            if save_output:
                self._diag_text = self.reaction_rates.format_diagnostic(
                    mean_n_i=accumulated['mean_n_i'][-1]
                    if 'mean_n_i' in accumulated else None)
                self._save_output(accumulated, solver_config)
        elif interrupted:
            print("\nNo completed segments to save.")

        return accumulated

    # ── Main entry point ──────────────────────────────────────────────────────

    def run(self, solver_config=None, save_output=True, progress_callback=None,
            timeout_s=None, y0_override=None):
        """
        Run the simulation using the configured solver mode and physics option.

        Parameters
        ----------
        solver_config     : dict, optional
            Keys: t_span, n_points, rtol, atol, log_time,
                  solver_method (dict with backend, linsol, window_*, etc.)
        save_output       : bool
            Write timestamped output/ directory.
        progress_callback : callable or None
            Called once per output time step during the C++ solver run with a
            dict of concentrations and rate-breakdown values (atom fraction).
            Enables verbose mode in the C++ solver automatically.
        timeout_s         : float or None
            Maximum wall-clock seconds for the C++ solver.
        y0_override       : ndarray or None
            Custom initial conditions (adaptive continuation).

        Returns
        -------
        results : dict  (see post_process.calculate_derived_quantities)
        """
        if solver_config is None:
            solver_config = self._default_solver_config()

        sm = self.input_data.solver_mode
        print(f"\nLaunching solver_mode='{sm}' …")

        if sm in ('cpp_full', 'cpp_sliding_win', 'sliding_OpenMP'):
            results = self._run_cpp(solver_config, progress_callback,
                                    timeout_s=timeout_s,
                                    y0_override=y0_override)
        else:
            raise ValueError(f"Unknown solver_mode='{sm}'. "
                             "Use cpp_full, cpp_sliding_win, or sliding_OpenMP.")

        if results is not None:
            # Store diagnostic text for file output (no inline printing)
            self._diag_text = self.reaction_rates.format_diagnostic(
                mean_n_i=results['mean_n_i'][-1] if 'mean_n_i' in results else None
            )
            if save_output:
                self._save_output(results, solver_config)

        return results

    # ── C++ solver dispatch ───────────────────────────────────────────────────

    def _run_cpp(self, solver_config, progress_callback=None, timeout_s=None,
                 y0_override=None):
        """Invoke the C++ solver via cpp_bridge."""
        from . import cpp_bridge
        results = cpp_bridge.run_cpp_solver(
            self, solver_config, base_dir=BASE_DIR,
            progress_callback=progress_callback,
            timeout_s=timeout_s,
            y0_override=y0_override,
        )
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
            f.write(f"- I:              {self.input_data.I}\n")
            f.write(f"- V:              {self.input_data.V}\n")
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

        # Diagnostics file
        diag_path = out_dir / 'diagnostics.txt'
        with open(diag_path, 'w') as f:
            f.write(f"# Diagnostics for {label}\n\n")
            f.write(f"## Rate constants\n")
            f.write(getattr(self, '_diag_text', '') + '\n\n')
            f.write(f"## Conservation diagnostics (final time step)\n")
            f.write(f"delta_FP = {results['delta_FP'][-1]:.6e}  (Frenkel pair)\n")
            f.write(f"delta_He = {results['delta_He'][-1]:.6e}  (He balance)\n\n")
            f.write(f"## Key results\n")
            f.write(f"Final dose:       {results['dose'][-1]:.4e} dpa\n")
            f.write(f"Swelling (final): {results['swelling'][-1]*100:.6f} %\n")
            f.write(f"C_He_tot (final): {results['C_He_tot'][-1]:.3e} m^-3\n")
            f.write(f"mean_n_i (final): {results['mean_n_i'][-1]:.2f}\n")
            f.write(f"mean_n_v (final): {results['mean_n_v'][-1]:.2f}\n\n")
            # Write progress table if collected
            prog = getattr(self, '_progress_rows', None)
            if prog:
                f.write(f"## Time-step diagnostics ({len(prog)} rows)\n")
                keys = list(prog[0].keys())
                f.write('\t'.join(keys) + '\n')
                for row in prog:
                    f.write('\t'.join(f"{row.get(k, 0.0):.6e}" for k in keys) + '\n')
        print(f"Diagnostics written to: {diag_path}")

        # Binary results (numpy)
        np.save(str(out_dir / 'results_t.npy'),  results['t'])
        np.save(str(out_dir / 'results_y.npy'),  results['y'])

        # Plots
        from . import visualization
        visualization.save_all_plots(results, self.input_data,
                                     str(plot_dir), label=f"{sm}/{po}",
                                     rate_eq_obj=self.rate_equations)

        print(f"Output written to: {out_dir}")
        return out_dir
