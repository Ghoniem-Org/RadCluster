# CUDA C++ Plan

**Project:** RadCluster_2_0 — Graph-Based Cluster Dynamics for Irradiated Materials
**Document:** Development plan for a GPU (CUDA C++) integrator channel
**Status:** Draft / design
**Scope:** Adds a third solver channel alongside the existing CPU C++ CVODE channel; defines the optimal-selection conditions for each channel and the switch-over rule.

---

## 1. Background — the current C++ channel (the baseline)

The existing full_system channel (cpp_utils/core/) is a single-CPU, serial-vector SUNDIALS CVODE BDF integrator:

- **State vector:** N_VNew_Serial (host memory, one CPU).
- **Linear solvers (solver.cpp):** SUNLinSol_Dense, SUNLinSol_Band, SUNLinSol_KLU (sparse direct on a CSC SUNSparseMatrix), SUNLinSol_SPGMR (matrix-free Krylov).
- **Sparse Jacobian (sparse_jacobian.cpp):** sparsity pattern derived from the mobile-species coupling (i_mobile / v_mobile), greedy column coloring (color_columns_greedy -> build_color_groups), then colored finite-difference perturbation (sparse_fd_jac) assembled into a CSC matrix for KLU.
- **Woodbury preconditioner (rhs_dispatch.cpp):** bordered-banded Sherman-Morrison-Woodbury. The banded core is factorized with LAPACK band LU (dgbtrf/dgbtrs); the dense capacitance/border block is factorized with dense LU (dgetrf/dgetrs). The band half-widths (kl, ku) and the border rank are set by the number of mobile species.
- **Parallelism today:** OpenMP only, and only in the RHS assembly (auto-picked from N_eq). Time integration and linear algebra are serial.

The two costs that dominate a step are (a) RHS evaluations (nfe) and (b) linear-solve work (nli for SPGMR, factor/solve for direct). These are the GPU targets.

---

## 2. Objective

Add a CUDA C++ channel (e.g. SOLVER_MODE = 'full_system_cuda') that moves the **RHS evaluation, the finite-difference Jacobian, and the linear solve** onto the GPU, while preserving the existing Python -> cpp_bridge.py -> subprocess -> .bin-results contract. The CPU channel remains the default and the correctness reference.

---

## 3. Architecture of the new channel

### 3.1 Integration boundary (Python side)
- Add 'full_system_cuda' as a SOLVER_MODE value; extend the linsol map in cpp_bridge.py with GPU variants (e.g. cuda_sparse, cuda_gmres).
- Add a CUDA-availability probe; if no GPU/toolkit is present, fall back to the CPU channel transparently.
- Reuse the existing input-file marshalling and .bin result parsing unchanged.

### 3.2 Device data + vector ops
- Move the state to SUNDIALS CUDA N_Vector (N_VNew_Cuda) so all CVODE vector operations (norms, AXPY, BDF history) run on device and per-step host<->device copies are eliminated.

### 3.3 GPU RHS kernel
- Port the per-size growth / shrinkage / emission accumulation (currently OpenMP-parallel) to a CUDA kernel parallelized over the size axis.
- The coalescence convolution (size a + size b -> a+b) is the highest arithmetic-intensity term and the strongest GPU candidate.

### 3.4 GPU Jacobian
- The colored finite-difference Jacobian maps naturally: each color group is one batched RHS evaluation, so coloring + batched perturbation becomes a batched GPU kernel feeding a device CSC matrix.

### 3.5 GPU linear solvers (two paths matched to Jacobian structure)
- **Sparse-direct path:** cuSOLVER/cuSPARSE sparse or batched factorization as the KLU analog (heavy-coalescence / broadly sparse Jacobian).
- **Iterative path:** GPU SPGMR + GPU Woodbury preconditioner. The banded core solve maps to cuSPARSE batched-band / gtsv2 routines; the low-rank border correction maps to small cuBLAS GEMM/GETRF (arrow / bordered-banded Jacobian).

---

## 4. Phased development

1. **Build infrastructure.** Add CUDA toolkit + SUNDIALS-CUDA detection to CMakeLists.txt, guarded so the CPU build is unaffected when no GPU is present. Add the 'full_system_cuda' plumbing in cpp_bridge.py.
2. **Vector + RHS port.** Switch to N_VNew_Cuda and a CUDA RHS kernel, keeping CPU linear solvers initially. Isolates and measures RHS speedup; validates correctness.
3. **GPU linear solvers.** Add the cuSOLVER sparse-direct path (KLU analog), then GPU SPGMR + GPU Woodbury preconditioner (SPGMR/Woodbury analog).
4. **Calibration + auto-dispatch.** Run the benchmark harness of Section 4A across N_eq and the two Jacobian regimes on the target GPU; fit the launch-amortization threshold W* and the regime boundaries; encode them as dispatcher constants implementing the decision rule (Section 7).
5. **Validation.** Compare delta_FP / delta_He conservation diagnostics and size distributions against the CPU channel across the full test matrix.

---

## 4A. Phase 4 in detail — calibration harness and auto-dispatch

The purpose of Phase 4 is to replace the order-of-magnitude estimates in Sections 5-7 with empirically fitted constants for a specific GPU + CPU pair, and to wire those constants into the runtime dispatcher. It has four parts: a benchmark matrix, an instrumented runner, a fitting step, and the dispatcher integration.

### 4A.1 Benchmark matrix (the N_eq x {arrow, coalescent} sweep)

Two structural regimes are swept independently, each across a geometric ladder of system sizes.

- **Size ladder.** N_eq in {2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6}, realized by scaling the cluster domain (I, V) and, where needed, the bin-moment expansion. The lower rungs overlap the current validated im5vm2 / I=V=1000 case (N_eq ~ 2e3) so the harness is anchored to known-good output.
- **ARROW regime.** i_mobile = v_mobile = 1 (single mobile monomers). Produces the bordered-banded Jacobian: narrow band + thin (rank ~ number of mobile species) border. Coalescence edges present but with only monomers mobile, so off-band fill is minimal.
- **COALESCENT regime.** i_mobile, v_mobile raised so a large fraction of sizes are mobile (e.g. mobile cutoffs at 10%, 25%, 50% of I, V), activating SIA_SIA / VAC_VAC coalescence across many sizes. Produces high off-band fill and a wide/high-rank border. Sweep the mobile fraction as a secondary axis to locate where the structure transitions from arrow-like to coalescent.
- **Solver cross-product.** For each (regime, N_eq) cell, run every applicable solver:
  - CPU: band, gmres+Woodbury, klu (full_CD discrete only), and dense at the smallest sizes as a reference.
  - CUDA: cuda_gmres+GPU-Woodbury and cuda_sparse (cuSOLVER).
- **Fixed controls.** Identical t_span, n_points, rtol, atol, physics_option, temperature, and PARAM_OVERRIDES across all cells so only the solver and structure vary. Pin OMP_NUM_THREADS for the CPU runs and record it (the CPU baseline must be its best honest configuration, not a single thread).

### 4A.2 Instrumented runner and metrics

A driver script (e.g. codes/Notebooks or a bench/ harness) loops over the matrix, launches each solver through cpp_bridge.py, and records per cell:

- **Primary:** wall-clock solve time (median of >= 3 repeats; discard a warm-up run to exclude JIT/allocation and first-touch GPU context cost).
- **CVODE work counters** (already emitted): nsteps, nfe (RHS evals), nni, nli, nli_per_nni, npe (preconditioner setups), nps (preconditioner solves), ncfn, netf.
- **Structure descriptors:** N_eq, measured Jacobian nnz and average nonzeros per row (from sparse_jacobian.cpp), band half-widths kl/ku, and Woodbury border rank.
- **Derived per-RHS flop estimate:** W = N_eq x (avg nonzeros per row).
- **GPU-only:** host<->device transfer time, kernel time vs total, peak device memory; flag any cell that exceeds device memory (these define the CUDA feasibility ceiling).
- **Correctness gate:** delta_FP, delta_He and final macroscopic outputs must match the CPU reference within tolerance, else the cell is marked invalid and excluded from fitting.

Results are written to a tidy CSV (one row per regime x N_eq x solver x repeat) plus a provenance record (git SHA, GPU model, driver/toolkit versions, CPU model, thread count).

### 4A.3 Fitting W* and the regime boundaries

- **Crossover W\*.** For each regime, plot best-CPU vs best-CUDA wall time against W (log-log). The crossover is where the two curves intersect; W* is taken as that intersection (optionally shifted by a safety margin so CUDA is only chosen when it wins by more than the measurement noise). Expect a lower W* in the COALESCENT regime than in ARROW.
- **Arrow-vs-coalescent boundary.** Along the mobile-fraction axis, locate where off-band fill (avg nonzeros per row, or border rank / N_eq) crosses the level at which cuda_sparse begins to beat cuda_gmres+Woodbury and at which CPU klu begins to beat CPU band/Woodbury. This fixes the structural classifier threshold used in Section 7.
- **Per-solver winner map.** Record, for each (regime, W) region, which single solver was fastest, so the dispatcher can pick not just CPU-vs-CUDA but the specific linsol/preconditioner.
- **Robustness.** Refit if the GPU, toolkit, or CPU thread count changes; store the calibration context alongside the constants so a stale calibration is detectable.

### 4A.4 Dispatcher integration

- Encode W*, the structural-classifier threshold, and the per-solver winner map as named constants (e.g. a small JSON or header consumed by cpp_bridge.py / solver.cpp), not magic numbers scattered in code.
- At run start, compute N_eq and the structure descriptors, estimate W, classify ARROW vs COALESCENT, and apply the Section 7 rule to select channel + linsol + preconditioner. Always fall back to CPU when no GPU is present or device memory is insufficient.
- Honor an explicit user override (SOLVER_MODE / solver_method) so the auto-dispatch can be bypassed for testing.
- Emit the chosen channel, the inputs to the decision (N_eq, W, regime, border rank), and the reason, into provenance.md for every run.

### 4A.5 Deliverables

A reproducible benchmark CSV + plots (wall time vs W per regime, CPU vs CUDA), the fitted constants file, and a short calibration report recording W*, the structural boundary, the per-solver winner map, and the hardware/software context they were measured on.

## 5. Optimal-selection conditions for the CPU C++ solver

Prefer the existing CPU channel when ANY of the following hold:

- **Small to moderate system:** N_eq below roughly 1e4. Per-step work is microseconds; serial CVODE vector ops are effectively free and GPU launch/transfer latency dominates.
- **Arrow / bordered-banded Jacobian with a thin border:** single (or few) mobile species (small i_mobile / v_mobile). The banded core + thin low-rank border is exactly the regime the LAPACK band + dense Woodbury preconditioner handles near-optimally. Banded triangular solves are inherently sequential and do not favor the GPU at modest size.
- **Latency-sensitive / many short runs:** parameter sweeps, adaptive-doubling segments, smoke tests, where kernel-launch and host<->device transfer overhead would dominate wall time.
- **No GPU available**, or GPU memory insufficient for the device state + Jacobian + Krylov basis.
- **Linear-solver choice within the CPU channel:**
  - dense -> only as a correctness check or for very small N_eq.
  - band -> arrow/banded structure, narrow bandwidth.
  - gmres + Woodbury -> arrow/bordered-banded; near-optimal while the border stays low-rank.
  - klu (full_CD discrete only) -> moderate heavy-coalescence sparsity with limited fill-in.

---

## 6. Optimal-selection conditions for the CUDA C++ solver

Prefer the CUDA channel when the per-step arithmetic is large enough to amortize launch/transfer overhead AND the structure exposes parallelism:

- **Large system:** N_eq in the ~1e4-1e5 break-even band and dominant above ~1e5-1e6 (reached via the 111/100 SIA loop split, multiple populations, large I/V domains, or full bin-moment expansions with many shape moments).
- **Heavy-coalescence / broadly-sparse or wide-bordered Jacobian:** many mobile sizes (large i_mobile / v_mobile). Coalescence produces dense off-band blocks and a high-rank border. This is where (a) the RHS convolution has high arithmetic intensity, (b) CPU sparse-direct suffers heavy fill-in but cuSOLVER batched/sparse factorization maps well, and (c) the Woodbury border correction becomes a sizeable GEMM (a GPU strength). The crossover N_eq drops here, possibly to a few x 1e3-1e4.
- **High RHS cost per step:** when nfe and (for SPGMR) nli are large and the RHS/Jacobian assembly dominates, the batched GPU RHS + colored-FD Jacobian win even before the linear solve does.
- **GPU linear-solver choice within the CUDA channel:**
  - cuda_sparse (cuSOLVER) -> heavy-coalescence, high-fill Jacobian.
  - cuda_gmres + GPU Woodbury -> large banded/bordered systems where the band solve and border GEMM are large enough to fill the device.

A practical dispatch metric: prefer CUDA when N_eq x (avg nonzeros per row) — i.e. the per-RHS flop count — exceeds a calibrated threshold, rather than N_eq alone.

---

## 7. Switch-over decision rule

Let N = N_eq, let W = per-RHS flop estimate (N x avg nonzeros/row), and let the calibrated launch-amortization threshold be W*. Classify the Jacobian as ARROW (thin border, few mobile species) or COALESCENT (high fill / wide border, many mobile species).

- **ARROW and N < ~1e4:** CPU, band or gmres+Woodbury.
- **ARROW and N >> 1e4 (W > W*):** CUDA, gmres + GPU Woodbury (gain comes mainly from batched RHS + Jacobian throughput, since the band solve is partly sequential).
- **COALESCENT and N moderate-to-large (W > W*):** CUDA — cuda_sparse (cuSOLVER) if fill is high, else cuda_gmres + GPU Woodbury.
- **COALESCENT and N small (W <= W*):** CPU, klu.
- **Any case with no GPU or insufficient device memory:** CPU (auto-fallback).

Summary: switch to CUDA when **W > W\*** AND (**system is large-banded** OR **coalescence has produced a high-fill / high-rank-border Jacobian**); otherwise stay on the CPU channel. Thresholds W* and the regime boundaries are fixed empirically in Phase 4 on the target GPU and stored as dispatcher constants.

---

## 8. Validation criteria

- Conservation diagnostics delta_FP (Eq. 96) and delta_He (Eq. 97) must match the CPU channel within tolerance across the test matrix.
- Matched size distributions and macroscopic outputs (swelling, mean sizes, loop/void densities) versus the CPU reference.
- Benchmark report: wall time vs N_eq for ARROW and COALESCENT regimes on CPU vs CUDA, used to fit W* and the crossover boundaries.

---

## 9. Open implementation details (confirmed from source)

- Band structure uses LAPACK kl/ku half-widths and ldab leading dimension; the Woodbury border is a dense block factorized with dgetrf — its rank is set by the mobile-species count, which sizes the GPU GEMM in the preconditioner.
- The sparse Jacobian sparsity pattern is derived from i_mobile / v_mobile coupling, then greedy-colored and built as CSC for KLU — this determines whether cuSOLVER batched-tridiagonal or general-sparse routines are the right GPU target.
