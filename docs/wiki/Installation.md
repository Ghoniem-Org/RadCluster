# Installation

There are two ways to get a working RadCluster: the published container image, which needs no build, and a build from source, which is faster for long runs.

## Container image

The image carries the compiled solver, the Python package, the input workbook and the scripts.

```bash
docker pull ghcr.io/ghoniem/radcluster:v2.2.0
docker run --rm -it ghcr.io/ghoniem/radcluster:v2.2.0
```

Inside the container:

- `solver` is on `PATH`;
- `/work` holds `py_utils/`, `input/` and `codes/`, and `PYTHONPATH` is `/work`;
- output is written inside the container and is lost when it exits.

Mount a host directory to keep the output:

```bash
docker run --rm -v "$PWD/output:/work/output" ghcr.io/ghoniem/radcluster:v2.2.0 \
    python3 codes/Python_Testing/check_eurofer_rag.py
```

The solver in the image is compiled for `x86-64-v2`, so it runs on any x86-64 processor from about 2009 onwards. A local build with `-march=native` is somewhat faster; use it for production runs.

The image is built and published by `.github/workflows/release.yml` when a `v*` tag is pushed. The workflow then starts the published image and runs `check_eurofer_rag.py` in it.

## Python environment

```bash
git clone https://github.com/Ghoniem/RadCluster.git
cd RadCluster
python -m pip install -r requirements.txt
python -m ipykernel install --user --name radcluster --display-name "RadCluster"   # optional
```

- Python 3.9 or later, with `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl`, `jupyter` and `ipykernel`.
- `networkx` is optional. Only the graph export in `py_utils/core/to_networkx.py` and the two graph figure scripts use it.

The Python package alone is enough to read the workbook, build a material graph, inspect the model and run the reference graph walker. A simulation needs the C++ solver.

## C++ solver

Required:

- a C++17 compiler (GCC, Clang or MSVC) and CMake 3.15 or later;
- SUNDIALS 7.1.1 (CVODE). The solver calls `SUNContext_Create(SUN_COMM_NULL, ...)`, which does not exist before version 7, so the version 6 packages of some Linux distributions do not compile.

Build:

```bash
cd radcluster_code
cmake -S cpp_utils -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

The executable is `build/solver`, or `build/Release/solver.exe` with a multi-configuration generator. `py_utils/cpp_bridge.py` looks for it in `build/Release/`, `build/Debug/` and `build/`, in that order.

**Where SUNDIALS is found.** CMake looks first in a `Libraries/sundials-7.1.1/` folder beside the repository (`mac_install/` on macOS, `win_install/` on Windows), then wherever `find_package(SUNDIALS)` resolves: a system install, a vcpkg install, or a prefix given with `-DCMAKE_PREFIX_PATH=<sundials install>`.

**Optional dependencies.** Each is detected at configure time, and the solver builds without it.

| Dependency | Enables | If missing |
|---|---|---|
| LAPACK and BLAS (Accelerate on macOS; OpenBLAS, reference LAPACK or MKL elsewhere) | the Woodbury preconditioner | it is compiled out, and the solver uses the Jacobi preconditioner |
| OpenMP (`brew install libomp` on macOS; included with MSVC) | the parallel sweeps of the right-hand side, in both solver modes | the solver runs on one thread |
| SuiteSparse/KLU, with SUNDIALS built `-DENABLE_KLU=ON` | `linsol='klu'`, the sparse direct solver | that option stops with a message |

**Compiler flags.** `-O3 -march=native`, or `/arch:AVX2 /GL /fp:precise` with MSVC. Fast-math is not enabled: it changes IEEE-754 arithmetic and can make CVODE's error estimate diverge on these stiff systems. For a build that must run on other machines, pass `-DRADCLUSTER_ARCH_FLAGS=-march=x86-64-v2`.

**From the notebook.** The first code cell of `RadCluster_2_1.ipynb` builds the solver when no executable for the current platform is found. After editing C++ sources, set `FORCE_REBUILD = True` in that cell.

## Windows notes

- The build copies the SUNDIALS, LAPACK or MKL, and KLU runtime DLLs beside the executable, so the binary runs from a fresh checkout.
- A vcpkg SUNDIALS built with KLU (triplet `x64-windows-rel`, under `%VCPKG_ROOT%` or `C:/vcpkg`) is searched first, before the `Libraries/` folder and before a system install.

## UCLA Hoffman2 cluster

- [docs/Hoffman2_Guide.md](https://github.com/Ghoniem/RadCluster/blob/main/docs/Hoffman2_Guide.md): logging in, interactive and batch jobs, and JupyterLab through an SSH tunnel.
- `scripts/h2jupynb` starts a Jupyter session on the cluster.
- `radcluster_code/digital_twin/hoffman2_setup.sh`, `hoffman2_probe.sh` and `hoffman2_array.sh` set up and submit the campaign as an array job. See [Calibration and verification campaign](Calibration-and-Verification-Campaign).

## Check the installation

```bash
cd radcluster_code
PYTHONPATH=$PWD python codes/Python_Testing/check_eurofer_rag.py
```

It reads the workbook, declares the EUROFER-97 graph for both equation forms and ends with `ALL CHECKS PASSED`. It does not need the C++ solver. Checks that do are listed in [Tests and conservation checks](Tests-and-Conservation-Checks).
