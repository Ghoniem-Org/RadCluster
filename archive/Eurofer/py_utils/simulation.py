# -*- coding: utf-8 -*-
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from py_utils.input_data import InputData
from py_utils.reaction_rates import ReactionRates
from py_utils.rate_equations import RateEquations

class EuroferMicroSimulation:
    def __init__(self):
        INPUT_FILE = Path(__file__).parent.parent / 'input' / 'input_parameters_Neutron.xlsx'
        print(f"Loading input data from {INPUT_FILE}...")
        self.input_data = InputData(INPUT_FILE)

        print("Initializing reaction rates...")
        self.reaction_rates = ReactionRates(self.input_data)

        print("Setting up rate equations...")
        self.rate_equations = RateEquations(self.input_data, self.reaction_rates)

    def _ode_wrapper(self,t, y):
        try:
            # Check for unphysical values
            if np.any(y < -1e-15):  # Allow small numerical errors
                neg_indices = np.where(y < -1e-15)[0]
                if len(neg_indices) < 5:  # Only show first few to avoid spam
                    neg_names = [self.rate_equations.concentration_names[i] for i in neg_indices[:3]]
                    print(f"Warning: Negative concentrations at t={t:.2e}: {neg_names}")

            if np.any(np.isnan(y)) or np.any(np.isinf(y)):
                raise ValueError(f"Invalid concentrations at t={t:.2e}")

            # Call the actual ODE system
            dydt = self.rate_equations.ode_system(t, y)

            # Check derivatives for sanity
            if np.any(np.isnan(dydt)) or np.any(np.isinf(dydt)):
                raise ValueError(f"Invalid derivatives at t={t:.2e}")

            return dydt

        except Exception as e:
            raise ValueError(f"Error in ODE evaluation at t={t:.2e}: {e}")

    def run_simulation(self, t_begin=1e-6, t_end=1e8, n_points=1000, method='LSODA', rtol=1e-6, atol=1e-20, log_time=True, max_step=None):
        print(f"Running simulation from 0 to {t_end:.1e} seconds...")
        print(f"Using {method} solver with rtol={rtol:.1e}, atol={atol:.1e}")

        # Set up time points
        if log_time:
            t_eval = np.logspace(np.log10(t_begin), np.log10(t_end), n_points)
        else:
            t_eval = np.linspace(0, t_end, n_points)

        # Get initial conditions
        y0 = self.rate_equations.get_initial_conditions()

        # Set default max_step if not provided
        if max_step is None:
            max_step = t_end / 100  # Prevent overly large steps

        try:
            print(f"\nIntegrating with {method} solver...")
            sol = solve_ivp(self._ode_wrapper, [0, t_end], y0, method=method, t_eval=t_eval, rtol=rtol, atol=atol, dense_output=True, max_step=max_step)

            if sol.success:
                print("Integration successful!")
                self._check_solution_quality(sol)
                return self._process_solution(sol)
            else:
                raise RuntimeError(f"Integration failed: {sol.message}")

        except Exception as e:
            print(f"Error during integration: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _check_solution_quality(self, sol):
        """
        Check the quality of the numerical solution

        Parameters:
        -----------
        sol : scipy.integrate.OdeResult
            Solution object from solve_ivp
        """
        print("Checking solution quality...")
        warnings_found = []

        # Check for negative concentrations
        min_values = np.min(sol.y, axis=1)
        if np.any(min_values < -1e-12):
            neg_indices = np.where(min_values < -1e-12)[0]
            neg_names = [self.rate_equations.concentration_names[i] for i in neg_indices]
            warnings_found.append(f"Significant negative concentrations: {neg_names}")

        # Check for conservation issues (simplified check)
        # Total point defects should be roughly conserved over short times
        if len(sol.t) > 10:
            early_total = np.sum(sol.y[:4, :10], axis=0)  # First 4 concentrations, first 10 points
            late_total = np.sum(sol.y[:4, -10:], axis=0)  # First 4 concentrations, last 10 points

            relative_change = np.abs(np.mean(late_total) - np.mean(early_total)) / np.mean(early_total)
            if relative_change > 1.0:  # More than 100% change might indicate issues
                warnings_found.append(f"Large total defect change: {relative_change * 100:.1f}%")

        # Check for numerical oscillations in final solution
        if len(sol.t) > 20:
            for i in range(sol.y.shape[0]):
                recent_values = sol.y[i, -20:]
                if len(recent_values) > 5 and np.mean(recent_values) > 1e-20:
                    relative_std = np.std(recent_values) / np.mean(recent_values)
                    if relative_std > 0.5:  # 50% relative standard deviation
                        name = self.rate_equations.concentration_names[i]
                        warnings_found.append(f"Oscillations in {name} (rel_std={relative_std:.2f})")

        # Check final time reached
        if hasattr(sol, 't_events') and sol.t_events:
            if sol.t[-1] < 0.9 * sol.t_events[0]:
                warnings_found.append("Simulation may have terminated early")

        # Report findings
        if warnings_found:
            print("Solution quality warnings:")
            for warning in warnings_found:
                print(f"    - {warning}")
        else:
            print("Solution quality: Good")

        return warnings_found

    def _process_solution(self, sol):
        """
        Process the raw solution into useful results

        Parameters:
        -----------
        sol : scipy.integrate.OdeResult
            Solution object from solve_ivp

        Returns:
        --------
        results : dict
            Processed results dictionary
        """
        print("Processing results...")

        time = sol.t
        concentrations = sol.y

        # Calculate derived quantities
        results = self._calculate_derived_quantities(time, concentrations)

        # Add raw output data
        results['raw'] = sol.y

        # Add metadata
        results['metadata'] = {
            'solver_stats': {
                'success': sol.success,
                'message': sol.message,
                'nfev': sol.nfev,
                'njev': getattr(sol, 'njev', None),
                'nlu': getattr(sol, 'nlu', None),
                'status': sol.status,
                't_events': getattr(sol, 't_events', None)
            },
            'parameters': {
                'temperature': self.input_data.material_params['T'],
                'dose_rate': self.input_data.material_params['G'],
                'dislocation_density': self.input_data.material_params['rho'],
                'end_time': time[-1],
                'n_time_points': len(time)
            }
        }
        print("Results processing complete!")
        return results

    def _calculate_derived_quantities(self, time, concentrations):
        """
        Calculate all derived physical quantities

        Parameters:
        -----------
        time : numpy.ndarray
            Time points
        concentrations : numpy.ndarray
            Concentration matrix [n_species x n_time_points]

        Returns:
        --------
        results : dict
            Results dictionary with derived quantities
        """
        n_points = len(time)

        # Initialize arrays for derived quantities
        loop_sizes = {
            'ril_111': np.zeros(n_points),
            'ril_100': np.zeros(n_points)
        }

        # Calculate quantities at each time point
        print("Calculating derived quantities...")
        progress_points = max(1, n_points // 10)

        Cveq = np.zeros(n_points)
        for i, t in enumerate(time):
            # Progress indicator
            if i % progress_points == 0:
                print(f"  Processing {100 * i / n_points:.0f}%")

            y_current = concentrations[:, i]
            y_current = np.maximum(y_current, 1e-20)

            try:
                # Update reaction rates for current state
                self.reaction_rates.update_state(y_current, t)

                # Calculate loop sizes
                loop_sizes['ril_111'][i] = self.rate_equations.calculate_ril_111(y_current)
                loop_sizes['ril_100'][i] = self.rate_equations.calculate_ril_100(y_current)
                Cveq[i] = self.reaction_rates.Cal_Cveq(y_current)

            except Exception as e:
                print(f"Warning: Error calculating derived quantities at t={t:.2e}: {e}")
                # Use previous values or defaults
                if i > 0:
                    loop_sizes['ril_111'][i] = loop_sizes['ril_111'][i - 1]
                    loop_sizes['ril_100'][i] = loop_sizes['ril_100'][i - 1]
                else:
                    # Set defaults for first point
                    loop_sizes['ril_111'][i] = 1e-9
                    loop_sizes['ril_100'][i] = 1e-9

        # Package results
        results = {
            'time': time,
            'concentrations': {
                name: concentrations[i, :]
                for i, name in enumerate(self.rate_equations.concentration_names)
            },
            'loop_sizes': loop_sizes,
            'Cveq': Cveq
        }

        return results

    def plot_results(self, results, save_plots=True, output_dir=None,
                     dpa_range=None, use_dpa=True):
        """
        Create comprehensive plots of simulation results

        Parameters:
        -----------
        results : dict
            Simulation results dictionary
        save_plots : bool
            Whether to save plots to file
        output_dir : str or Path
            Directory to save plots
        dpa_range : tuple, optional
            (dpa_min, dpa_max) range for plotting. If None, uses full range
        use_dpa : bool
            If True, plot vs DPA. If False, plot vs time
        """

        print(f"\n PLOT DEBUG INFO:")
        print(f"   use_dpa = {use_dpa}")
        print(f"   dpa_range = {dpa_range}")

        # Set default output directory relative to EuroferExperiments/
        if output_dir is None:
            output_dir = Path(__file__).parent.parent / 'output' / 'plots'
        else:
            output_dir = Path(output_dir)

        if save_plots:
            output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize ALL variables at the start to avoid scope issues
        x_min_limit, x_max_limit = None, None
        x_scale = 'log'

        time = results['time']
        conc = results['concentrations']
        loops = results['loop_sizes']

        # Calculate DPA from time and dose rate
        if use_dpa:
            dose_rate = self.input_data.material_params['G']
            x_data = time * dose_rate  # DPA = G * t
            x_label = 'DPA (displacements per atom)'

            print(f"   Dose rate G = {dose_rate:.2e} dpa/s")
            print(f"   Calculated DPA range: {x_data[0]:.2e} to {x_data[-1]:.2e}")

            # Apply DPA range filter if specified
            if dpa_range is not None:
                dpa_min, dpa_max = dpa_range
                x_min_limit, x_max_limit = dpa_min, dpa_max  # Store for axis limits
                print(f"   Setting x-axis limits to: {x_min_limit:.2e} to {x_max_limit:.2e}")

                mask = (x_data >= dpa_min) & (x_data <= dpa_max)

                if np.sum(mask) == 0:
                    print(f"Warning: No data points in DPA range {dpa_min:.2e} to {dpa_max:.2e}")
                    print(f"Available DPA range: {x_data[0]:.2e} to {x_data[-1]:.2e}")
                    # Don't filter if no data in range
                    mask = np.ones(len(x_data), dtype=bool)
                    x_min_limit, x_max_limit = None, None
                    print(f"   Reverting to full range")
                else:
                    print(f"   Filtering data: {len(x_data)} → {np.sum(mask)} points")
                    x_data = x_data[mask]

                    # Apply mask to all data arrays
                    for key in conc:
                        conc[key] = conc[key][mask]
                    for key in loops:
                        loops[key] = loops[key][mask]

                    print(f"Plotting DPA range: {dpa_min:.2e} to {dpa_max:.2e}")
                    print(f"   Final data range: {x_data[0]:.2e} to {x_data[-1]:.2e}")
            else:
                print(f"Plotting full DPA range: {x_data[0]:.2e} to {x_data[-1]:.2e}")
        else:
            x_data = time
            x_label = 'Time (s)'

            # Apply time range filter if specified and dpa_range is given
            if dpa_range is not None:
                dose_rate = self.input_data.material_params['G']
                t_min, t_max = dpa_range[0] / dose_rate, dpa_range[1] / dose_rate
                x_min_limit, x_max_limit = t_min, t_max
                mask = (x_data >= t_min) & (x_data <= t_max)

                if np.sum(mask) == 0:
                    print(f"Warning: No data points in time range {t_min:.2e} to {t_max:.2e} s")
                    mask = np.ones(len(x_data), dtype=bool)
                    x_min_limit, x_max_limit = None, None
                else:
                    x_data = x_data[mask]

                    # Apply mask to all data arrays
                    for key in conc:
                        conc[key] = conc[key][mask]
                    for key in loops:
                        loops[key] = loops[key][mask]

                    print(f"Plotting time range: {t_min:.2e} to {t_max:.2e} s")

        # Create comprehensive figure
        fig, axes = plt.subplots(3, 3, figsize=(18, 12))
        plot_title = f"EuroferMicro Simulation Results (G = {self.input_data.material_params['G']:.2e} dpa/s, temp = {self.input_data.material_params['T']} K)"
        fig.suptitle(plot_title, fontsize=16, fontweight='bold')
        Omega = self.input_data.physical_props['Omega'] * 1e6  # Convert to cm³

        # Plot 1: Point defect concentrations
        ax = axes[0, 0]
        # ax.loglog(x_data, conc['Cv'] / Omega, 'b-', label='Vacancies', linewidth=2)
        # ax.loglog(x_data, conc['Ci'] / Omega, 'r-', label='Interstitials', linewidth=2)
        # ax.loglog(x_data, conc['C2i'] / Omega, 'g--', label='Di-interstitials', linewidth=2)
        # ax.loglog(x_data, conc['C3i'] / Omega, 'm--', label='Tri-interstitials', linewidth=2)
        # ax.loglog(x_data, conc['CiL_i_111'] / Omega, 'y--', label='1/2<111> i-Loop Interstitials', linewidth=2)
        # ax.loglog(x_data, conc['CiL_i_100'] / Omega, 'k--', label='<100> i-Loop Interstitials', linewidth=2)
        ax.loglog(x_data, conc['Cv'] , 'b-', label='Vacancies', linewidth=2)
        ax.loglog(x_data, conc['Ci'] , 'r-', label='Interstitials', linewidth=2)
        ax.loglog(x_data, conc['C2i'] , 'g--', label='Di-interstitials', linewidth=2)
        ax.loglog(x_data, conc['C3i'] , 'm--', label='Tri-interstitials', linewidth=2)
        ax.loglog(x_data, conc['CiL_i_111'], 'y--', label='1/2<111> i-Loop Interstitials', linewidth=2)
        ax.loglog(x_data, conc['CiL_i_100'], 'k--', label='<100> i-Loop Interstitials', linewidth=2)
        ax.set_xlabel(x_label)
        # ax.set_ylabel('Concentration')
        ax.set_title('Point Defect Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 2: Loop concentrations (convert to number density)
        ax = axes[0, 1]
        Omega = self.input_data.physical_props['Omega']
        ax.loglog(x_data, conc['CiL_111'] / Omega, 'b-', label='1/2<111> Interstitial Loops', linewidth=2)
        ax.loglog(x_data, conc['CiL_100'] / Omega, 'r-', label='<100> Interstitial Loops', linewidth=2)
        ax.set_xlabel(x_label)
        ax.set_ylabel('Number Density ($m^{-3}$)')
        ax.set_title('Loop Density Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 3: Loop sizes
        ax = axes[0, 2]
        ax.plot(x_data, loops['ril_111'] * 1e9, 'b--', label='1/2<111> Interstitial Loops', linewidth=2)
        ax.plot(x_data, loops['ril_100'] * 1e9, 'r--', label='<100> Interstitial Loops', linewidth=2)
        ax.set_xlabel(x_label)
        ax.set_ylabel('Loop Radius (nm)')
        ax.set_title('Loop Size Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 4: Void density
        ax = axes[1, 0]
        Omega = self.input_data.physical_props['Omega']
        ax.loglog(x_data,conc['C_void'] / Omega, 'b-', label='Void', linewidth=2)
        ax.set_xlabel(x_label)
        ax.set_ylabel("Number Density ($m^{-3}$)")
        ax.set_title('Void Density Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 5: Void radius
        ax = axes[1, 1]
        ax.plot(x_data, conc['r_void'] * 1e9, 'b-', label='Void', linewidth=2)
        ax.set_xlabel(x_label)
        ax.set_ylabel('Void Radius (nm)')
        ax.set_title('Void Size Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 6:
        ax = axes[1, 2]
        ax.loglog(x_data, self.reaction_rates.omega_v*(conc['Cv']-results['Cveq']), 'b-', label='Vacancies flux', linewidth=2)
        ax.loglog(x_data, self.reaction_rates.omega_i*conc['Ci'], 'r-', label='Interstitials flux', linewidth=2)
        ax.set_xlabel(x_label)
        ax.set_ylabel('Relative Flux')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 7:
        ax = axes[2, 0]
        ax.loglog(x_data, conc['Ctrap_i'], 'r-', label='Ctrap_i', linewidth=2)
        ax.loglog(x_data, conc['Ctrap_v'], 'b-', label='Ctrap_v', linewidth=2)
        ax.set_xlabel(x_label)
        ax.set_ylabel('Trap Occupancy (at/at)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_plots:
            if use_dpa:
                if self.input_data.material_params['G'] > 1e-5:
                    plot_file = output_dir / f'Eurofermicro_results_dpa_{int(self.input_data.material_params["T"]-273)}C_Ion.png'
                else:
                    plot_file = output_dir / f'Eurofermicro_results_dpa_{int(self.input_data.material_params["T"] - 273)}C_Neutron.png'
            else:
                plot_file = output_dir / f'Eurofermicro_results_time_{int(self.input_data.material_params["T"]-273)}C.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            print(f"Plots saved to {plot_file}")

        plt.show()

    def save_results(self, results, filename=None):
        """Save results to file in the output directory"""

        # Set default filename in output directory
        if filename is None:
            output_dir = Path(__file__).parent.parent / 'output'
            if self.input_data.material_params['G'] > 1e-5:
                filename = output_dir / f'Eurofermicro_results_{int(self.input_data.material_params["T"] - 273)}C_Ion.pkl'
            else:
                filename = output_dir / f'Eurofermicro_results_{int(self.input_data.material_params["T"] - 273)}C_Neutron.pkl'

        else:
            filename = Path(filename)

        # Ensure output directory exists
        filename.parent.mkdir(parents=True, exist_ok=True)

        try:
            with open(filename, 'wb') as f:
                pickle.dump(results, f)
            print(f"Results saved to {filename}")

            # Also save a summary CSV with proper path handling
            # Create CSV filename by replacing .pkl with _summary.csv
            csv_filename = str(filename).replace('.pkl', '_summary.csv')
            csv_file = Path(csv_filename)
            self._save_summary_csv(results, csv_file)

        except Exception as e:
            print(f"Error saving results: {e}")

    def _save_summary_csv(self, results, filename):
        """Save a summary CSV file"""
        try:
            # Convert Path object to ensure proper handling
            filename = Path(filename)
            df_data = {
                'time': results['time'],
                'Cv': results['concentrations']['Cv'],
                'Ci': results['concentrations']['Ci'],
                # 'N_void': results['concentrations']['N_void'],
                'C2i': results['concentrations']['C2i'],
                'C3i': results['concentrations']['C2i'],
                'C_void': results['concentrations']['C_void'],
                'CiL_111': results['concentrations']['CiL_111'],
                'CiL_100': results['concentrations']['CiL_100'],
                'ril_111_nm': results['loop_sizes']['ril_111'] * 1e9,
                'ril_100_nm': results['loop_sizes']['ril_100'] * 1e9,
                'r_void_nm': results['concentrations']['r_void'] * 1e9,
            }
            df = pd.DataFrame(df_data)
            df.to_csv(filename, index=False)
            print(f"Summary saved to {filename}")

        except Exception as e:
            print(f"Error saving summary CSV: {e}")


def run_eurofermicro_simulation(SIMULATION_CONFIG):
    try:
        # Initialize simulation
        print("\n" + "=" * 50)
        print("INITIALIZING SIMULATION")
        print("=" * 50)

        sim = EuroferMicroSimulation()

        # Run simulation
        print("\n" + "=" * 50)
        print("RUNNING SIMULATION")
        print("=" * 50)

        results = sim.run_simulation(
            t_begin=SIMULATION_CONFIG['t_begin'],
            t_end=SIMULATION_CONFIG['t_end'],
            n_points=SIMULATION_CONFIG['n_points'],
            method=SIMULATION_CONFIG['method'],
            rtol=SIMULATION_CONFIG['rtol'],
            atol=SIMULATION_CONFIG['atol'],
            log_time=SIMULATION_CONFIG['log_time']
        )

        if results is not None:
            # Create plots
            print("\n" + "=" * 50)
            print("GENERATING PLOTS")
            print("=" * 50)

            sim.plot_results(
                results,
                save_plots=True,
                output_dir=str(Path(__file__).parent.parent / 'output' / 'plots')
            )

            # Save results
            print("\n" + "=" * 50)
            print("SAVING RESULTS")
            print("=" * 50)

            sim.save_results(results)

            # Print final summary
            print_simulation_summary(results, sim, SIMULATION_CONFIG)

            return results, sim
        else:
            print("Simulation failed!")
            return None, None

    except Exception as e:
        print(f"Fatal error during simulation: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def print_simulation_summary(results, sim, SIMULATION_CONFIG):
    """Print a comprehensive simulation summary"""
    print("\n" + "=" * 60)
    print("SIMULATION SUMMARY")
    print("=" * 60)

    # Simulation parameters
    print("Simulation Parameters:")
    print(f"  Temperature: {sim.input_data.material_params['T']} K")
    print(f"  Dose rate: {sim.input_data.material_params['G']:.2e} dpa/s")
    print(f"  Total time: {results['time'][-1]:.2e} seconds")
    print(f"  Solver: {SIMULATION_CONFIG['method']}")
    print(f"  Tolerances: rtol={SIMULATION_CONFIG['rtol']:.1e}, atol={SIMULATION_CONFIG['atol']:.1e}")

    # Final concentrations
    print("\nFinal Concentrations:")
    for name, values in results['concentrations'].items():
        print(f"  {name}: {values[-1]:.2e}")

    # Final loop sizes
    print("\nFinal Loop Sizes:")
    for name, values in results['loop_sizes'].items():
        print(f"  {name}: {values[-1] * 1e9:.2f} nm")

    # Solver statistics
    if 'metadata' in results:
        stats = results['metadata']['solver_stats']
        print(f"\nSolver Statistics:")
        print(f"  Success: {stats['success']}")
        print(f"  Function evaluations: {stats['nfev']}")
        if stats['njev']:
            print(f"  Jacobian evaluations: {stats['njev']}")

    print("=" * 60)
