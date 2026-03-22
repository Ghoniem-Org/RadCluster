# ClusterDynamics — Claude Code Instructions

## Objectives

Model the time evolution of radiation-induced defect clusters in **316 stainless steel at 450 °C** under constant displacement damage. The goal is to predict:

- Transient and steady-state concentrations of vacancy clusters (up to size Nv) and interstitial clusters (up to size Ni)
- Total vacancy and interstitial content as a function of time
- Mean cluster sizes and size distributions
- Cluster front propagation for very large systems (Nv ~ 10⁴, Ni ~ 10⁶)

The theoretical framework is **Ghoniem & Cho (1979)** cluster dynamics for metallic alloys under irradiation.

---

## Physical Model & Equations

### State Vector

```
y[N] = [Cv1, Cv2, ..., Cv_Nv, Ci1, Ci2, ..., Ci_Ni]
```

- `Cvx` — concentration of vacancy cluster containing `x` vacancies (atom fraction)
- `Cix` — concentration of interstitial cluster containing `x` interstitials (atom fraction)
- Total equations: `N = Nv + Ni`  (default: 50 + 100 = 150; production: up to 10,000 + 1,000,000)

### Rate Equations

**Monomers (x = 1):**

```
dCv1/dt = P - α·Cv1·Ci1 - Kc_v(1)·Cv1·Cv1 + Γcv(2)·Cv2 - Σ[Kc_v(x)·Cv1·Cvx] - ρd·Dv·Cv1
dCi1/dt = P - α·Cv1·Ci1 - Kc_i(1)·Ci1·Ci1 + Γci(2)·Ci2 - Σ[Kc_i(x)·Ci1·Cix] - K_nuc_i·Ci1²
```

**Small clusters (x = 2):**

```
dCv2/dt = ½·Kc_v(1)·Cv1² - [Kc_v(2)·Cv1 + Kc_i(2)·Ci1 + Γcv(2)]·Cv2 + Γcv(3)·Cv3 - ρd·Dv·Cv2
```

**Large clusters (x ≥ 3):**

```
dCvx/dt = Kc_v(x-1)·Cv1·Cv(x-1) - [Kc_v(x)·Cv1 + Kc_i(x)·Ci1 + Γcv(x)]·Cvx + Γcv(x+1)·Cv(x+1)
dCix/dt = Kl_i(x-1)·Ci1·Ci(x-1) - [Kl_v(x)·Cv1 + Kl_i(x)·Ci1 + Γci(x)]·Cix + Γci(x+1)·Ci(x+1)
```

### Rate Constants

**Spherical vacancy cluster capture (size x):**
```
Kc_v(x) = 2.216 · Zv² · x^(2/3) / (1 + 0.1128·Zv·x^(1/3)) · νv · exp(-E_m_v / kT)
Kc_i(x) = 2.216 · Zi² · x^(2/3) / (1 + 0.1128·Zi·x^(1/3)) · νi · exp(-E_m_i / kT)
```

**Circular interstitial loop capture (size x):**
```
Kl_v(x) = 1.555 · Zv · x^0.5 · νv · exp(-E_m_v / kT)
Kl_i(x) = 1.555 · Zi · x^0.5 · νi · exp(-E_m_i / kT)
```

**Thermal emission (vacancy from size-x cluster):**
```
Γcv(x) = Kc_v(x) · Cv_eq · exp((1.28 · g · a² / kT) · x^(-1/3))
```

### Material Parameters (316 SS, Ghoniem & Cho 1979)

| Symbol | Value | Description |
|---|---|---|
| T | 723.15 K | Temperature (450 °C) |
| P | 1×10⁻⁶ dpa/s | Displacement damage rate |
| E_m_v | 1.6 eV | Vacancy migration energy |
| E_m_i | 0.2 eV | Interstitial migration energy |
| a | 3.63 Å | Lattice parameter |
| Cv_eq | ~7×10⁻¹² | Thermal equilibrium vacancy concentration |
| ρd | network dislocation density | Sink strength for point defects |

---

## Solution Algorithm

### Python Path (LSODA, segmented)

1. `InputData` loads parameters and computes derived quantities (diffusion coefficients, equilibrium concentrations, attempt frequencies)
2. `ReactionRates` pre-computes all six rate-constant arrays `KCV`, `KCI`, `KLV`, `KLI`, `GCV`, `GLV` at initialization
3. `RateEquations` provides `ode_system(t, y)` using pre-computed arrays (no recomputation per call)
4. `ClusterDynamicsSimulation.run_simulation()` integrates using `scipy.integrate.solve_ivp(method='LSODA')`:
   - Time span divided into `n_segments` (default 60) logarithmically-spaced segments
   - An `x_max` heuristic gates the upper interstitial equations: `dCix/dt = 0` for `x > x_max`
   - `x_max` is predicted from linear growth on the log-time axis and updated between segments
5. `post_process.calculate_derived_quantities()` computes totals, mean sizes, and validates solution quality

### C++ Path (CVODE, continuous)

1. `cpp_bridge.write_param_file()` serializes all parameters to a text file (`key=value` format, arrays as `KCV_0`, `KCV_1`, ...)
2. `solver.exe --param_file=<path>` is launched as a subprocess
3. CVODE BDF integrates all N equations in one continuous run (no gating needed)
4. Output: `n_points` rows × `(1 + N)` columns, space-separated scientific notation
5. `cpp_bridge.run_cpp_solver()` parses stdout into NumPy arrays and invokes post-processing

### C++ Window Solver Modes

| Mode | Name | Description |
|---|---|---|
| 0 | Full | All N equations simultaneously (baseline) |
| 1 | Phase I | Upper truncation: active window `[1..x_hi]`; expands when `C[x_hi] > C_expand` |
| 2 | Phase II | Sliding window: upper expansion + lower contraction (QSS criterion); Jacobi preconditioner |
| 3 | Phase III | Constant-width window `W`; slides upward; nucleation guard at `t < t_start` |
| 4 | Phase IV | Same as Phase III + OpenMP intra-RHS parallelism + pre-allocated scratch buffers |

---

## Development Methodology: Phase I → II → III

### Phase I — Dynamic Upper-Truncation Window

**Goal:** Enable large Nv/Ni (up to 1000/10,000) by starting with a small active window and expanding only when needed.

**Method:**
- Start with `[1..w0_v]` × `[1..w0_i]` active (default `w0_v = w0_i = 100`)
- At each output point, check if `C[x_hi] > window_C_expand` (e.g., 1×10⁻²⁰)
- If triggered, expand `x_hi += window_expand_pad` (e.g., +20 equations) and reinitialize CVODE
- Equations above `x_hi` are zeroed (frozen); full Nv+Ni output row is always written
- GMRES linear solver used (matrix-free, adapts to changing system size)

**Achieved:** ~40× speedup vs full solver for Ni=10,000 (6 s vs 239 s); 0.013% error on Cv1 steady state.

**Key parameters:** `window_mode=1`, `window_w0_v`, `window_w0_i`, `window_C_expand`, `window_expand_pad`

---

### Phase II — Sliding Window with Lower Contraction (QSS)

**Goal:** Further reduce active equations by also freezing small clusters that have reached quasi-steady state, enabling Ni up to 100,000.

**Method:**
- Retains the Phase I upper expansion
- **Geometric expansion:** `x_hi *= window_expand_factor` (e.g., ×2) instead of additive `+pad`, reducing CVODE reinitializations from ~50,000 to ~20
- **Lower contraction:** A cluster `x` is frozen (dropped from CVODE) when its relative rate of change satisfies `|dC/dt|/C < window_C_contract` (e.g., 1×10⁻³); minimum `window_min_active_i` equations always kept active
- **Jacobi diagonal preconditioner** (`window_prec=1`): preconditiones GMRES iterations using the diagonal of the Jacobian, improving convergence for the large banded system
- Frozen concentrations contribute to the RHS of active equations as constant source/sink corrections (`recompute_frozen_sums()`)
- Active window at time t: `[Ci1] ∪ [Ci_{x_lo}..Ci_{x_hi}]`

**Achieved:** ~0.01% accuracy vs full reference; enables Ni = 100,000.

**Key parameters:** `window_mode=2`, `window_C_contract`, `window_min_active_i`, `window_expand_factor`, `window_prec`

---

### Phase III — Constant-Width Sliding Window

**Goal:** Scale to Ni = 10⁶ (1,010,000 total equations) while maintaining accuracy with minimal active equations.

**Method:**
- Window of fixed width `W` (e.g., W = 1000 interstitial equations) slides rigidly upward with the cluster front
- Lower bound: `x_lo = max(2, x_hi - W + 1)` — always exactly W equations active
- **Nucleation guard:** Lower sliding suppressed until `t > window_t_start` (e.g., 10 s); prevents destroying the small-cluster nucleation front during early transient
- **Auto-activation:** Phase III only engages when `N_EQ > window_N_thresh` (default 2000); smaller systems fall back to Phase II
- GMRES + Jacobi preconditioner required at this scale
- Output always writes all N rows; inactive equations hold their last frozen value

**Achieved (Nv=10⁴, Ni=10⁶):** 7 CVODE reinitializations total; final active window Ci[15001..16000] (1.5% of Ni); wall time ~534 s (old build, -O2, no SIMD); Cv1(t_end) = 3.264×10⁻⁶.

**Key parameters:** `window_mode=3`, `window_width`, `window_t_start`, `window_N_thresh`

---

### Phase IV — Multithread-OpenMP

**Goal:** Exploit the M3 Max's 12 P-cores to reduce wall time beyond what Phase III achieves with the sliding-window algorithm alone.

**Method:**
Three stacked speedup layers, all activated by `window_mode=4`:

1. **Compiler (`-O3 -march=native`):** ARM NEON SIMD auto-vectorisation of all loops; benefits Phase III too. CMake detects Apple Homebrew libomp automatically.
   - **Note:** `-ffast-math` is intentionally omitted — it alters IEEE-754 behaviour and causes CVODE's floating-point error estimator to diverge on stiff problems.
2. **Pre-allocated scratch buffers (`WindowDataOMP::Cv_buf`, `Ci_buf`):** Eliminates `malloc`/`free` on every RHS call (millions of calls over a full run). Buffers are resized on each CVODE reinitialisation (`resize_buffers()`).
3. **OpenMP intra-RHS parallelism (`rhs_window_omp()`):** Uses a **single** `#pragma omp parallel` region per RHS call (one fork-join, not many), with `OMP_MIN_WORK=20000` threshold:
   - If `n_ci_win < 20000` (current W=1000 case): falls back to the serial buffer-optimised path automatically
   - If `n_ci_win ≥ 20000` (larger windows): parallel loops activate for buffer fill, KLV/KLI reductions, dCvx, dCix, and floor enforcement

**Code structure:** `WindowDataOMP` (rate_equations.h), `rhs_window_omp` / `prec_setup_window_omp` / `prec_solve_window_omp` / `recompute_frozen_sums_omp` (rate_equations.cpp). All guarded by `#ifdef CD_HAVE_OPENMP` so the code compiles without OpenMP too (graceful fallback to Phase III).

**Key parameters:** `window_mode=4`, `window_omp_threads` (0 = `OMP_NUM_THREADS`; 10 recommended for M3 Max)

**Benchmark results (Nv=10⁴, Ni=10⁶, n_points=20):**

| Configuration | Wall time | Speedup vs original |
|---|---|---|
| Phase III (old build, -O2, no SIMD) | ~534 s | 1× (reference) |
| Phase III (new build, -O3 -march=native) | ~78 s | **6.9×** |
| Phase IV (4/8/10 threads, W=1000) | ~80 s | **6.7×** |

At W=1000 the `OMP_MIN_WORK=20000` threshold forces the serial fallback, so Phase IV matches Phase III speed. The 6.9× gain comes entirely from `-O3 -march=native` (ARM NEON SIMD), shared by both phases. Phase IV OpenMP will activate and provide additional speedup when the window width reaches ≥ 20,000 (larger Ni problems).

Cv1(t_end) = 3.2636×10⁻⁶ — identical across Phase III and all Phase IV thread counts (0.000% relative error).

**Notes:**
- At Nv=10⁴, Ni=10⁶ with n_points=1000, the stdout I/O pipeline (1000 rows × 1,010,001 numbers ≈ 20 GB text through subprocess pipe) can dominate total wall time, masking pure ODE integration speedup. Use fewer output points for pure timing benchmarks.
- `brew install libomp` required before cmake on macOS; cmake auto-detects via `brew --prefix libomp`.

---

## File Map

```
ClusterDynamics/
├── codes/
│   ├── cluster_dynamics.ipynb    Main notebook: Cells 0–5 (config, Python solver,
│   │                             C++ full, Phase I, Phase II, Phase III)
│   ├── CD.ipynb                  Legacy/alternative notebook
│   ├── solver.cpp                C++ entry point: CLI parsing, CVODE/ARKODE setup,
│   │                             window mode dispatch, stdout output
│   └── test_window_solver.py     Python test utilities for window solver validation
│
├── cpp_utils/
│   ├── CMakeLists.txt            CMake build: C++17, SUNDIALS 7.1.1, links cvode/arkode/gmres/band
│   ├── parameters.h              Parameters struct: Nv/Ni/N_EQ, rate arrays, window config,
│   │                             solver settings, initial conditions y0[]
│   ├── rate_equations.h          RHS function declarations: rhs_cd() (full),
│   │                             rhs_window() (Phases I–III), WindowData struct,
│   │                             recompute_frozen_sums(), prec_setup_window()
│   └── rate_equations.cpp        C++ ODE RHS implementation (mirrors Python _rhs_full)
│
├── py_utils/
│   ├── input_data.py             InputData class: 316 SS parameters (Ghoniem & Cho 1979),
│   │                             Nv/Ni/C_floor, derived Dv/Di/Cv_eq/alpha
│   ├── reaction_rates.py         ReactionRates class: arrays KCV/KCI/KLV/KLI/GCV/GLV
│   ├── rate_equations.py         RateEquations class: ode_system(), gated_ode_system(),
│   │                             _rhs_full() (150-equation NumPy implementation)
│   ├── simulation.py             ClusterDynamicsSimulation: segmented LSODA, x_max gating,
│   │                             _predict_xmax(), _make_gated_rhs()
│   ├── cpp_bridge.py             write_param_file() (key=value serialization),
│   │                             run_cpp_solver() (subprocess + output parsing)
│   ├── post_process.py           calculate_derived_quantities() (totals, mean sizes),
│   │                             check_solution_quality()
│   ├── pre_process.py            Input validation, parameter range warnings
│   └── visualization.py          plot_results(), CDVisualizer, create_run_directory()
│
├── utilities/
│   └── ghoniem_cho_network.py    Legacy network model (unused)
│
├── build/                        CMake artifacts (Debug/ and Release/ subdirs)
└── output/                       Timestamped run directories (<YYYYMMDD_HHMMSS_git-hash>/)
                                  Each contains: provenance.md + 13 PNG figures
```

### File Relationships

```
cluster_dynamics.ipynb
  ├── Cell 1 (Python LSODA)
  │     └── simulation.py → rate_equations.py → reaction_rates.py → input_data.py
  │                       → post_process.py → visualization.py
  │
  ├── Cells 2–5 (C++ solver paths)
  │     ├── simulation.py (init only) → input_data.py + reaction_rates.py
  │     └── cpp_bridge.py
  │           ├── write_param_file() → parameters.h (data contract)
  │           └── solver.exe (subprocess)
  │                 ├── parameters.h (struct)
  │                 ├── rate_equations.h/.cpp (C++ RHS)
  │                 └── [stdout] → cpp_bridge.py → post_process.py → visualization.py
  │
  └── Cell 0 (PLOT_CONFIG) → visualization.py (axis limits)
```

---

## Build and Execution

```bash
# Build C++ solver
cd ClusterDynamics/cpp_utils
cmake -S . -B ../build -DCMAKE_BUILD_TYPE=Release
cmake --build ../build --config Release
# Output: ClusterDynamics/build/Release/solver.exe

# Run notebook
.Zr_venv/Scripts/python.exe -m nbconvert --to notebook --execute \
    --ExecutePreprocessor.kernel_name=fluor_zr codes/cluster_dynamics.ipynb
```

---

## Known Notes

- `np.maximum(y, 1e-100)` applied in Python RHS to prevent negative concentrations
- C++ uses `C_floor = 1e-100` parameter to clamp concentrations
- Parameter file I/O used (not CLI args) to avoid Windows command-line length limits for large Nv/Ni
- Phase III at Nv=10⁴, Ni=10⁶ requires ~8 GB RAM for the full output array
