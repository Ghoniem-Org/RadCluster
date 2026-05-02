"""
pre_process.py – Initial conditions and parameter validation.

Called before the ODE solver is invoked.
"""

import warnings
import numpy as np


def validate_setup(input_data):
    """Warn about parameter values outside typical physical ranges."""
    p = input_data.material_params
    T     = p['T']
    P     = p['P']
    rho_d = p['rho_d']

    if not (300 <= T <= 1200):
        warnings.warn(f"Temperature {T} K outside typical range [300–1200 K]")
    if not (1e-8 <= P <= 1e-2):
        warnings.warn(f"Production rate {P} dpa/s outside typical range [1e-8 – 1e-2]")
    if not (1e7 <= rho_d <= 1e14):
        warnings.warn(f"Dislocation density {rho_d} cm/cm³ outside typical range [1e7 – 1e14]")

    Cv_eq = input_data.derived['Cv_eq']
    if Cv_eq > 1e-4:
        warnings.warn(f"Cv_eq={Cv_eq:.2e} seems high — check temperature/E_f_v")

    print("✓ Setup validation complete.")


def get_initial_conditions(rate_equations, custom_values=None):
    """
    Build and validate the initial condition vector.

    Parameters
    ----------
    rate_equations : RateEquations
    custom_values  : dict, optional – overrides for specific species names

    Returns
    -------
    y0 : numpy.ndarray  shape (N,)
    """
    y0 = rate_equations.get_initial_conditions()

    if custom_values:
        names = rate_equations.concentration_names
        for name, val in custom_values.items():
            if name in names:
                y0[names.index(name)] = val

    if np.any(y0 < 0):
        raise ValueError("Negative initial concentrations detected")
    if np.any(np.isnan(y0)) or np.any(np.isinf(y0)):
        raise ValueError("Invalid initial concentrations (NaN or Inf)")

    return y0
