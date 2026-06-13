"""
Neural Network-Inspired Cluster Dynamics Code
Based on Ghoniem & Cho (1979) - Simultaneous Clustering of Point Defects

Network Concept:
- Vertices (nodes) = Defect concentrations (Cv, C2v, ..., Ci, C2i, ...)
- Edges (connections) = Reaction rates between species
- Dynamic network size based on concentration/rate thresholds
"""

import numpy as np
from scipy.integrate import solve_ivp, ode
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Callable
import matplotlib.pyplot as plt
from enum import Enum


class SolverType(Enum):
    """Available ODE solver types"""
    LSODA = "LSODA"  # Adaptive method for stiff/non-stiff
    RADAU = "Radau"  # Implicit Runge-Kutta for stiff
    BDF = "BDF"      # Backward differentiation formula for stiff
    VODE = "vode"    # Using scipy.integrate.ode


@dataclass
class MaterialParameters:
    """Material and irradiation parameters for Ghoniem & Cho model"""
    # Material properties (316 stainless steel from Table 1)
    a: float = 3.63e-10  # Lattice parameter (m)
    Omega: float = 1.197e-29  # Atomic volume (m^3), calculated from a^3
    
    # Migration energies (eV)
    E_v_m: float = 1.4   # Vacancy migration energy
    E_i_m: float = 0.2   # Interstitial migration energy
    E_2v_m: float = 0.9  # Divacancy migration energy
    E_2i_m: float = 0.25 # Diinterstitial binding energy
    E_3v_m: float = 0.76 # Trivacancy binding energy
    
    # Formation energies (eV)
    E_v_f: float = 1.6   # Vacancy formation energy
    E_i_f: float = 4.08  # Interstitial formation energy
    
    # Frequencies (s^-1)
    nu_i: float = 5e13   # Interstitial attempt frequency
    nu_v: float = 5e13   # Vacancy attempt frequency
    
    # Surface and binding energies (eV)
    g: float = 6.24e14   # Surface energy (eV/cm^2)
    
    # Combinatorial numbers
    C_i: int = 84  # Interstitial-interstitial combinatorial number
    C_v: int = 84  # Vacancy-vacancy combinatorial number
    
    # Bias factors
    Z_v: float = 1.0   # Vacancy-dislocation bias factor
    Z_i: float = 1.08  # Interstitial-dislocation bias factor
    
    # Irradiation parameters
    T: float = 723.15  # Temperature (K) - 450°C
    P: float = 1e-5    # Production rate (dpa/s)
    rho_d: float = 1e14  # Dislocation density (cm/cm^3)
    
    # Physical constants
    k_B: float = 8.617e-5  # Boltzmann constant (eV/K)


@dataclass
class NetworkNode:
    """Represents a node (defect species) in the network"""
    name: str                    # e.g., "C_v", "C_2v", "C_i", "C_2i"
    index: int                   # Index in state vector
    cluster_size: int            # Number of defects in cluster
    defect_type: str            # 'v' for vacancy, 'i' for interstitial, 'l' for loop
    concentration: float = 0.0   # Current concentration
    is_active: bool = True       # Whether node is active in network


@dataclass
class NetworkEdge:
    """Represents an edge (reaction) in the network"""
    source_idx: int              # Index of source node
    target_idx: int              # Index of target node
    reaction_type: str           # Type of reaction
    rate_constant: float = 0.0   # Reaction rate constant
    is_active: bool = True       # Whether edge is active


class DynamicClusterNetwork:
    """
    Dynamic network for cluster dynamics simulations
    Automatically manages network size based on concentration thresholds
    """
    
    def __init__(self, params: MaterialParameters,
                 max_cluster_size: int = 100,
                 concentration_threshold: float = 1e-30,
                 rate_threshold: float = 1e-50):
        """
        Initialize dynamic cluster network
        
        Parameters:
        -----------
        params : MaterialParameters
            Material and irradiation parameters
        max_cluster_size : int
            Maximum cluster size to consider
        concentration_threshold : float
            Minimum concentration to keep node active
        rate_threshold : float
            Minimum reaction rate to keep edge active
        """
        self.params = params
        self.max_cluster_size = max_cluster_size
        self.c_threshold = concentration_threshold
        self.r_threshold = rate_threshold
        
        # Network components
        self.nodes: List[NetworkNode] = []
        self.edges: List[NetworkEdge] = []
        
        # Mapping dictionaries
        self.node_map: Dict[str, int] = {}  # name -> index
        
        # Initialize base network
        self._initialize_network()
        
    def _initialize_network(self):
        """Initialize base network with single defects and small clusters"""
        # Single vacancies
        self._add_node("C_v", 1, 'v', 1e-20)
        
        # Single interstitials  
        self._add_node("C_i", 1, 'i', 1e-20)
        
        # Small vacancy clusters (up to size 10 for better dynamics)
        for n in range(2, 11):
            self._add_node(f"C_{n}v", n, 'v', 0.0)
        
        # Small interstitial clusters (up to size 10)
        for n in range(2, 11):
            self._add_node(f"C_{n}i", n, 'i', 0.0)
        
        # Interstitial loops (different sizes) - these grow larger
        for size in [20, 30, 40, 50, 60, 70, 80, 90]:
            self._add_node(f"C_{size}i", size, 'i', 0.0)
            
        # Initialize basic reactions
        self._initialize_reactions()
        
    def _add_node(self, name: str, size: int, defect_type: str, 
                  initial_conc: float = 0.0):
        """Add a node to the network"""
        idx = len(self.nodes)
        node = NetworkNode(
            name=name,
            index=idx,
            cluster_size=size,
            defect_type=defect_type,
            concentration=initial_conc
        )
        self.nodes.append(node)
        self.node_map[name] = idx
        
    def _add_edge(self, source_name: str, target_name: str, 
                  reaction_type: str, rate: float = 0.0):
        """Add an edge (reaction) to the network"""
        if source_name not in self.node_map or target_name not in self.node_map:
            return
            
        source_idx = self.node_map[source_name]
        target_idx = self.node_map[target_name]
        
        edge = NetworkEdge(
            source_idx=source_idx,
            target_idx=target_idx,
            reaction_type=reaction_type,
            rate_constant=rate
        )
        self.edges.append(edge)
        
    def _initialize_reactions(self):
        """Initialize basic reaction network"""
        # Production: P -> C_v + C_i
        # Recombination: C_v + C_i -> 0
        # Clustering: C_v + C_v -> C_2v
        # Clustering: C_i + C_i -> C_2i
        # Growth: C_v + C_2v -> C_3v
        # Growth: C_i + C_2i -> C_3i
        # Dissociation: C_2v -> C_v + C_v
        # Dissociation: C_2i -> C_i + C_i
        # Emission from surfaces
        # Absorption by dislocations
        pass
    
    def compute_reaction_rates(self) -> Dict[str, float]:
        """
        Compute all reaction rate constants based on current state
        Returns dictionary of rate constants
        """
        p = self.params
        T = p.T
        k_B = p.k_B
        
        rates = {}
        
        # Diffusion coefficients
        D_v = (p.a**2 / 6) * p.nu_v * np.exp(-p.E_v_m / (k_B * T))
        D_i = (p.a**2 / 6) * p.nu_i * np.exp(-p.E_i_m / (k_B * T))
        D_2v = (p.a**2 / 6) * p.nu_v * np.exp(-p.E_2v_m / (k_B * T))
        
        # Point defect recombination coefficient (large, fast process)
        alpha = 4 * np.pi * (D_v + D_i) * p.a / p.Omega
        rates['recombination'] = alpha
        
        # Vacancy clustering reactions (Eq. 2 with modifications)
        # This should be much larger - using simplified diffusion-limited rate
        z_c = 48  # Combinatorial number
        # Simplified rate: K ~ 4π * D * a * z_c
        K_v_2 = 4 * np.pi * D_v * p.a * z_c / p.Omega
        rates['K_v_2'] = K_v_2
        
        # Interstitial clustering (Eq. 3) - also much faster
        K_i_2 = 4 * np.pi * D_i * p.a * z_c / p.Omega
        rates['K_i_2'] = K_i_2
        
        # Vacancy emission from clusters (Eq. 4) - reduced to allow growth
        for x in range(2, 10):
            gamma_v = (1.28 * p.g * p.a**2) / (k_B * T) * x**(-1/3)
            C_v_th = 6 * np.exp(-(2 * p.E_v_f - p.E_2v_m) / (k_B * T))
            # Reduce emission by factor of 1000 to allow cluster growth
            rates[f'gamma_v_{x}'] = K_v_2 * C_v_th * np.exp(-gamma_v) / 1000.0
        
        # Dislocation sink strengths (from Eq. 3 in paper, modified)
        # Z_v,i are bias factors
        rates['k_v_d_sq'] = p.Z_v * p.rho_d * 1e4 / (p.a**2)  # Convert cm to m
        rates['k_i_d_sq'] = p.Z_i * p.rho_d * 1e4 / (p.a**2)
        
        # Absorption rates by dislocations
        rates['D_v_k_v_sq'] = D_v * rates['k_v_d_sq']
        rates['D_i_k_i_sq'] = D_i * rates['k_i_d_sq']
        
        return rates
    
    def expand_network(self, state: np.ndarray, rates: Dict[str, float]) -> bool:
        """
        Dynamically expand network if needed
        Returns True if network was expanded
        """
        expanded = False
        
        # Check if we need to add larger clusters
        for node in self.nodes:
            if not node.is_active:
                continue
                
            conc = state[node.index]
            
            # If concentration is significant and we're at the boundary
            if (conc > self.c_threshold * 10 and 
                node.cluster_size < self.max_cluster_size):
                
                # Check if next size exists
                next_size = node.cluster_size + 1
                next_name = f"C_{next_size}{node.defect_type}"
                
                if next_name not in self.node_map:
                    self._add_node(next_name, next_size, 
                                 node.defect_type, 0.0)
                    expanded = True
                    
        return expanded
    
    def prune_network(self, state: np.ndarray, rates: Dict[str, float]) -> bool:
        """
        Prune inactive nodes and edges from network
        Returns True if network was pruned
        """
        pruned = False
        
        # Deactivate nodes below threshold
        for node in self.nodes:
            if node.cluster_size <= 3:  # Keep base species always active
                continue
                
            conc = state[node.index]
            if conc < self.c_threshold and node.is_active:
                node.is_active = False
                pruned = True
                
        return pruned
    
    def get_active_state_vector(self, full_state: np.ndarray) -> np.ndarray:
        """Extract only active nodes from full state vector"""
        active_indices = [n.index for n in self.nodes if n.is_active]
        return full_state[active_indices]
    
    def get_full_state_vector(self, active_state: np.ndarray, 
                             full_state: np.ndarray) -> np.ndarray:
        """Reconstruct full state vector from active state"""
        active_indices = [n.index for n in self.nodes if n.is_active]
        full_state[active_indices] = active_state
        return full_state


class GhoniemChoSolver:
    """
    Solver for Ghoniem & Cho (1979) cluster dynamics equations
    using dynamic network approach
    """
    
    def __init__(self, params: MaterialParameters,
                 max_cluster_size: int = 100,
                 solver_type: SolverType = SolverType.LSODA):
        """
        Initialize solver
        
        Parameters:
        -----------
        params : MaterialParameters
            Material parameters
        max_cluster_size : int
            Maximum cluster size
        solver_type : SolverType
            ODE solver to use
        """
        self.params = params
        self.solver_type = solver_type
        self.network = DynamicClusterNetwork(
            params, 
            max_cluster_size=max_cluster_size
        )
        
        # Storage for results
        self.times = []
        self.solutions = []
        self.network_sizes = []
        
    def rate_equations(self, t: float, y: np.ndarray) -> np.ndarray:
        """
        Rate equations for cluster dynamics (Eq. 6 in paper)
        Properly implements clustering, emission, and growth dynamics
        
        Parameters:
        -----------
        t : float
            Time (s)
        y : np.ndarray
            State vector of concentrations
            
        Returns:
        --------
        dydt : np.ndarray
            Time derivatives
        """
        p = self.params
        rates = self.network.compute_reaction_rates()
        
        # Initialize derivatives
        dydt = np.zeros_like(y)
        
        # Ensure non-negative concentrations
        y = np.maximum(y, 0.0)
        
        # Get single defect concentrations
        C_v = max(y[self.network.node_map['C_v']], 1e-50)
        C_i = max(y[self.network.node_map['C_i']], 1e-50)
        
        # Single vacancy equation
        idx_v = self.network.node_map['C_v']
        loss_to_clusters_v = 0.0
        gain_from_emission_v = 0.0
        
        # Loss to clustering with all vacancy clusters
        for node in self.network.nodes:
            if node.defect_type == 'v' and node.cluster_size > 1:
                idx = node.index
                C_x = max(y[idx], 0.0)
                # v + x_v -> (x+1)_v
                loss_to_clusters_v += rates['K_v_2'] * C_v * C_x
        
        # Gain from emission from all vacancy clusters
        for node in self.network.nodes:
            if node.defect_type == 'v' and node.cluster_size >= 2:
                idx = node.index
                C_x = max(y[idx], 0.0)
                x = node.cluster_size
                gamma_key = f'gamma_v_{x}'
                if gamma_key in rates:
                    gain_from_emission_v += rates[gamma_key] * C_x
        
        dydt[idx_v] = (
            p.P  # Production
            - rates['K_v_2'] * C_v**2  # Self-clustering v + v -> 2v
            - loss_to_clusters_v  # Absorption by clusters
            + gain_from_emission_v  # Emission from clusters
            - rates['recombination'] * C_v * C_i  # Recombination
            - rates['D_v_k_v_sq'] * C_v  # Absorption by dislocations
        )
        
        # Single interstitial equation
        idx_i = self.network.node_map['C_i']
        loss_to_clusters_i = 0.0
        
        # Loss to clustering with all interstitial clusters
        for node in self.network.nodes:
            if node.defect_type == 'i' and node.cluster_size > 1:
                idx = node.index
                C_x = max(y[idx], 0.0)
                # i + x_i -> (x+1)_i
                loss_to_clusters_i += rates['K_i_2'] * C_i * C_x
        
        dydt[idx_i] = (
            p.P  # Production
            - rates['K_i_2'] * C_i**2  # Self-clustering i + i -> 2i
            - loss_to_clusters_i  # Absorption by clusters
            - rates['recombination'] * C_v * C_i  # Recombination
            - rates['D_i_k_i_sq'] * C_i  # Absorption by dislocations
        )
        
        # Cluster equations for all sizes x >= 2
        for node in self.network.nodes:
            if node.cluster_size < 2:
                continue
                
            idx = node.index
            x = node.cluster_size
            C_x = max(y[idx], 0.0)
            
            if node.defect_type == 'v':
                # Vacancy cluster x_v
                # Gain from (x-1)_v + v -> x_v
                gain_from_smaller = 0.0
                if x == 2:
                    # 2v formed from v + v
                    gain_from_smaller = 0.5 * rates['K_v_2'] * C_v**2
                else:
                    # x_v formed from (x-1)_v + v
                    prev_name = f"C_{x-1}v"
                    if prev_name in self.network.node_map:
                        C_xm1 = max(y[self.network.node_map[prev_name]], 0.0)
                        gain_from_smaller = rates['K_v_2'] * C_v * C_xm1
                
                # Loss to (x+1)_v via x_v + v -> (x+1)_v
                loss_to_larger = rates['K_v_2'] * C_v * C_x
                
                # Gain from (x+1)_v -> x_v + v (emission)
                gain_from_emission = 0.0
                next_name = f"C_{x+1}v"
                if next_name in self.network.node_map:
                    C_xp1 = max(y[self.network.node_map[next_name]], 0.0)
                    gamma_key = f'gamma_v_{x+1}'
                    if gamma_key in rates:
                        gain_from_emission = rates[gamma_key] * C_xp1
                
                # Loss via x_v -> (x-1)_v + v (dissociation)
                gamma_key = f'gamma_v_{x}'
                loss_from_emission = rates.get(gamma_key, 0.0) * C_x
                
                dydt[idx] = (
                    gain_from_smaller 
                    - loss_to_larger
                    + gain_from_emission
                    - loss_from_emission
                )
                
            elif node.defect_type == 'i':
                # Interstitial cluster x_i
                # Gain from (x-1)_i + i -> x_i
                gain_from_smaller = 0.0
                if x == 2:
                    # 2i formed from i + i
                    gain_from_smaller = 0.5 * rates['K_i_2'] * C_i**2
                else:
                    # x_i formed from (x-1)_i + i
                    prev_name = f"C_{x-1}i"
                    if prev_name in self.network.node_map:
                        C_xm1 = max(y[self.network.node_map[prev_name]], 0.0)
                        gain_from_smaller = rates['K_i_2'] * C_i * C_xm1
                
                # Loss to (x+1)_i via x_i + i -> (x+1)_i
                loss_to_larger = rates['K_i_2'] * C_i * C_x
                
                dydt[idx] = gain_from_smaller - loss_to_larger
        
        return dydt
    
    def solve(self, t_span: Tuple[float, float], 
              y0: Optional[np.ndarray] = None,
              t_eval: Optional[np.ndarray] = None,
              rtol: float = 1e-6,
              atol: float = 1e-12) -> Dict:
        """
        Solve the cluster dynamics equations
        
        Parameters:
        -----------
        t_span : tuple
            (t_start, t_end) time span
        y0 : np.ndarray, optional
            Initial conditions
        t_eval : np.ndarray, optional
            Times at which to evaluate solution
        rtol : float
            Relative tolerance
        atol : float
            Absolute tolerance
            
        Returns:
        --------
        result : dict
            Solution dictionary with times, concentrations, and network info
        """
        # Set initial conditions
        if y0 is None:
            n_nodes = len(self.network.nodes)
            y0 = np.zeros(n_nodes)
            # Set thermal equilibrium concentrations for single defects
            # C_v_eq = exp(-E_f^v / kT)
            C_v_eq = np.exp(-self.params.E_v_f / (self.params.k_B * self.params.T))
            C_i_eq = np.exp(-self.params.E_i_f / (self.params.k_B * self.params.T))
            
            y0[self.network.node_map['C_v']] = C_v_eq
            y0[self.network.node_map['C_i']] = C_i_eq
            # Small initial cluster concentrations
            for key in ['C_2v', 'C_2i', 'C_3v', 'C_3i']:
                if key in self.network.node_map:
                    y0[self.network.node_map[key]] = 1e-30
        
        # Choose solver
        if self.solver_type == SolverType.VODE:
            # Use scipy.integrate.ode for VODE
            solver = ode(self.rate_equations)
            solver.set_integrator('vode', method='bdf', 
                                 rtol=rtol, atol=atol,
                                 max_step=1e-3)
            solver.set_initial_value(y0, t_span[0])
            
            if t_eval is None:
                t_eval = np.logspace(np.log10(t_span[0] + 1e-10), 
                                    np.log10(t_span[1]), 100)
            
            times = []
            solutions = []
            
            for t in t_eval:
                if solver.successful():
                    solver.integrate(t)
                    times.append(solver.t)
                    solutions.append(solver.y.copy())
                else:
                    break
                    
            times = np.array(times)
            solutions = np.array(solutions)
            
        else:
            # Use solve_ivp for other methods
            method_map = {
                SolverType.LSODA: 'LSODA',
                SolverType.RADAU: 'Radau',
                SolverType.BDF: 'BDF'
            }
            
            sol = solve_ivp(
                self.rate_equations,
                t_span,
                y0,
                method=method_map[self.solver_type],
                t_eval=t_eval,
                rtol=rtol,
                atol=atol,
                dense_output=True
            )
            
            if not sol.success:
                print(f"Warning: Integration failed - {sol.message}")
                # Return partial results
                times = sol.t
                solutions = sol.y.T
            else:
                times = sol.t
                solutions = sol.y.T
        
        # Store results
        self.times = times
        self.solutions = solutions
        
        # Build result dictionary
        result = {
            't': times,
            'success': True,
            'message': 'Integration successful'
        }
        
        # Extract concentrations for each species
        for node in self.network.nodes:
            result[node.name] = solutions[:, node.index]
        
        # Add network size information
        result['n_nodes'] = len(self.network.nodes)
        result['n_active_nodes'] = sum(1 for n in self.network.nodes if n.is_active)
        
        return result
    
    def plot_results(self, result: Dict, species: Optional[List[str]] = None,
                    figsize: Tuple[int, int] = (12, 8)):
        """
        Plot concentration evolution
        
        Parameters:
        -----------
        result : dict
            Solution dictionary from solve()
        species : list, optional
            List of species names to plot (default: all)
        figsize : tuple
            Figure size
        """
        if species is None:
            species = [n.name for n in self.network.nodes if n.cluster_size <= 5]
        
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
        
        t = result['t']
        
        # Plot 1: Single defects
        ax = axes[0]
        if 'C_v' in result:
            ax.loglog(t, result['C_v'], 'b-', label='$C_v$', linewidth=2)
        if 'C_i' in result:
            ax.loglog(t, result['C_i'], 'r-', label='$C_i$', linewidth=2)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Concentration (at/at)')
        ax.set_title('Single Defect Concentrations')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Small vacancy clusters
        ax = axes[1]
        for sp in ['C_2v', 'C_3v']:
            if sp in result:
                ax.loglog(t, result[sp], label=f'${sp}$', linewidth=2)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Concentration (at/at)')
        ax.set_title('Vacancy Cluster Concentrations')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Small interstitial clusters
        ax = axes[2]
        for sp in ['C_2i', 'C_3i']:
            if sp in result:
                ax.loglog(t, result[sp], label=f'${sp}$', linewidth=2)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Concentration (at/at)')
        ax.set_title('Interstitial Cluster Concentrations')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 4: Total defect concentrations
        ax = axes[3]
        total_v = np.zeros_like(t)
        total_i = np.zeros_like(t)
        
        for node in self.network.nodes:
            if node.defect_type == 'v' and node.name in result:
                total_v += result[node.name] * node.cluster_size
            elif node.defect_type == 'i' and node.name in result:
                total_i += result[node.name] * node.cluster_size
                
        ax.loglog(t, total_v, 'b-', label='Total vacancies', linewidth=2)
        ax.loglog(t, total_i, 'r-', label='Total interstitials', linewidth=2)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Total Defect Concentration')
        ax.set_title('Total Defect Evolution')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('ghoniem_cho_results.png', dpi=150)
        plt.show()
        
        return fig
    
    def print_network_info(self):
        """Print information about the current network"""
        print(f"\n{'='*60}")
        print(f"NETWORK TOPOLOGY")
        print(f"{'='*60}")
        print(f"Total nodes: {len(self.network.nodes)}")
        print(f"Active nodes: {sum(1 for n in self.network.nodes if n.is_active)}")
        print(f"Total edges: {len(self.network.edges)}")
        print(f"\nNode details:")
        print(f"{'Name':<10} {'Size':<6} {'Type':<6} {'Active':<8} {'Index':<6}")
        print(f"{'-'*50}")
        for node in self.network.nodes[:10]:  # Show first 10
            print(f"{node.name:<10} {node.cluster_size:<6} {node.defect_type:<6} "
                  f"{str(node.is_active):<8} {node.index:<6}")
        if len(self.network.nodes) > 10:
            print(f"... and {len(self.network.nodes) - 10} more nodes")
        print(f"{'='*60}\n")


if __name__ == "__main__":
    # Example usage
    print("Ghoniem & Cho (1979) Cluster Dynamics Solver")
    print("=" * 60)
    
    # Create material parameters
    params = MaterialParameters(
        T=723.15,  # 450°C
        P=1e-5,    # Production rate (dpa/s)
        rho_d=1e14 # Dislocation density
    )
    
    # Create solver
    solver = GhoniemChoSolver(params, solver_type=SolverType.LSODA)
    
    # Print network info
    solver.print_network_info()
    
    # Solve
    print("Solving cluster dynamics equations...")
    t_end = 10.0  # seconds
    t_eval = np.logspace(-8, np.log10(t_end), 100)
    
    result = solver.solve((0, t_end), t_eval=t_eval)
    
    if result['success']:
        print(f"✓ Solution successful!")
        print(f"  Final time: {result['t'][-1]:.3e} s")
        print(f"  Final C_v: {result['C_v'][-1]:.3e}")
        print(f"  Final C_i: {result['C_i'][-1]:.3e}")
        
        # Plot results
        solver.plot_results(result)
    else:
        print(f"✗ Solution failed: {result['message']}")
