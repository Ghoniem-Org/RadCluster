# Code reference (Doxygen)

Doxygen builds a browsable reference of the code from its docstrings and header comments. It is not published on the web yet; build it locally.

## Build it

Install [Doxygen](https://www.doxygen.nl/download.html) (1.9 or later) and Python 3, then from the repository root:

```bash
doxygen docs/doxygen/Doxyfile
```

Open `docs/doxygen/build/html/index.html`. The build takes a few seconds.

## What it covers

| Covered | Source |
|---|---|
| the Python package: abstract core, EUROFER-97 declaration, inputs, rates, orchestration, post-processing | `radcluster_code/py_utils/` |
| the C++ solver: CVODE driver, preconditioners, sparse Jacobian, EUROFER-97 kernels | `radcluster_code/cpp_utils/` |
| the physics and numerics checks, and the figure generators | `radcluster_code/codes/` |
| the calibration and verification campaign | `radcluster_code/digital_twin/` |
| the physics and solver reference, as a page with its equations typeset | `radcluster_code/CLAUDE.md` |

Not covered: the notebooks, the input workbook, the experimental database, the derivations and papers in `docs/Formulation/`, and the campaign's design and result files.

## Finding things

| To find | Look in |
|---|---|
| a Python class | *Classes*, under its package: `py_utils.core.rag.ReactionAdmissibilityGraph`, `py_utils.simulation.RadClusterSimulation` |
| a Python module and its functions | *Namespaces*: `py_utils.reaction_rates`, `py_utils.binding_energies` |
| a script | *Namespaces*, under the script's name: `check_he_conservation`, `run_ensemble` |
| the C++ parameter structure | *Classes*: `Parameters`, `UserData` |
| a C++ file | *Files*: `solver.cpp`, `rate_kernels.cpp`, `rhs_dispatch.cpp` |
| the equations | *Related Pages*: Physics and Solver Reference |

Every page links to the annotated source.

## How the code is documented

The code is not written with Doxygen commands. Python docstrings are shown as written, and every class, function and file is listed whether or not it has a comment. The C++ headers use `/** ... */` comments, which Doxygen reads directly.

When adding code, keep to that:

- give every Python module a docstring that says what it does and which section or equations of the manuscript it implements;
- give every public class and function a docstring, with its parameters and units;
- give every C++ function declared in a header a `/** ... */` comment.

## From GitHub

The workflow `.github/workflows/documentation.yml` builds the reference on every push to `main` and keeps the HTML as a downloadable artifact of the run (Actions tab). It publishes to GitHub Pages only when it is run by hand, and only after Pages has been enabled in the repository settings with "GitHub Actions" as the source.

## These pages

The wiki pages are kept in the repository, in `docs/wiki/`, one Markdown file per page. Edit them there and copy them to the wiki.

Details: [docs/doxygen/README.md](https://github.com/Ghoniem/RadCluster/blob/main/docs/doxygen/README.md).
