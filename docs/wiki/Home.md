# RadCluster: graph-based cluster dynamics for irradiated materials

RadCluster simulates the evolution of radiation-induced defect populations in crystalline solids under neutron or ion irradiation: self-interstitial atom (SIA) dislocation loops, vacancy and helium–vacancy cavities, and free interstitial helium. The reference host material is the reduced-activation ferritic–martensitic steel **EUROFER-97**.

**Current version:** 2.2.0. See the [README](https://github.com/Ghoniem/RadCluster/blob/main/README.md) for the release notes.

## Pages

| To | Read |
|---|---|
| install and run | [Installation](Installation), [Quick start](Quick-Start), [Running a simulation](Running-a-Simulation) |
| set a case and read its results | [Inputs and parameters](Inputs-and-Parameters), [Outputs and provenance](Outputs-and-Provenance) |
| understand the design | [Architecture](Architecture), [Physics and equations](Physics-and-Equations), [Solver](Solver), [Glossary](Glossary) |
| understand the EUROFER-97 model | [EUROFER-97 reaction graph](EUROFER-97-Reaction-Graph), [Loop character and network](Loop-Character-and-Network), [Calibration and verification campaign](Calibration-and-Verification-Campaign) |
| work on the code | [Tests and conservation checks](Tests-and-Conservation-Checks), [Adding a host material](Adding-a-Host-Material), [Code reference](Code-Reference) |

> **Name.** The active module is `radcluster_code/`. It was called `RadCluster_2_1` before version 2.1.1, and that name remains in run labels, in file headers and in the driver notebook `RadCluster_2_1.ipynb`.

## The governing idea

> A material is declared. Adding a host material adds a declaration and its rate kernels; it does not rewrite the solver.

A material is declared as a **Reaction Admissibility Graph** (RAG):

- its **vertices** are the admissible clusters, each identified by polarity (vacancy or interstitial type), size, population and trapped-solute composition;
- its **edges** are the admissible elementary reactions, each of one of ten abstract classes;
- its **kernels** are the size-resolved rate constants.

A reaction that is physically forbidden is never an edge, and the sign pattern of each edge class is fixed, so the signed defect content is conserved by construction.

The graph is used in two ways:

- **In Python**, the reference right-hand side is assembled by walking the graph, edge by edge. It is the executable specification, and the substrate of the bin-moment reduction.
- **In C++**, the production right-hand side is the kernel file of the material, written for the same reactions and integrated by CVODE. Every simulation runs through it, and carries the graph as the description of what was run.

See [Architecture](Architecture) for what each side does today.

## Two layers

| Layer | Location | Contents |
|---|---|---|
| 1: core | `py_utils/core/`, `cpp_utils/core/` | in Python: the cluster identifier, the ten edge classes, the graph, the graph walker, the state layout, the state-space reductions. In C++: the CVODE driver and its preconditioners |
| 2: material | `py_utils/materials/eurofer97/`, `cpp_utils/materials/eurofer97/` | the EUROFER-97 populations, edges and rate kernels (processes P1–P8) |

The Python core carries no assumption about bcc iron or EUROFER-97. The C++ core is the solver layer that a second material would share.

## What a run chooses

| Choice | Values | Meaning |
|---|---|---|
| equations | `discrete`, `bin_moment` | one equation per cluster size, or logarithmic bins carrying moments |
| cascade | `fission`, `fusion` | the helium reduction: a decoupled inventory for fission, a mean-field loading for fusion. The cascade source itself follows the workbook's `spectrum` |
| solver mode | `full_system`, `active_window` | CVODE on the full state vector, or on two sliding size windows |
| linear solver | `dense`, `band`, `gmres`, `klu` | with a Jacobi or a Woodbury preconditioner for GMRES |
| helium kinetics | `dynamic`, `quasi_steady_state` | free helium as an unknown, or eliminated algebraically |

## Getting started

Without building anything, with the published container image:

```bash
docker pull ghcr.io/ghoniem/radcluster:v2.2.0
docker run --rm -it ghcr.io/ghoniem/radcluster:v2.2.0
```

From source:

```bash
git clone https://github.com/Ghoniem/RadCluster.git
cd RadCluster
python -m pip install -r requirements.txt
cd radcluster_code
cmake -S cpp_utils -B build -DCMAKE_BUILD_TYPE=Release     # needs SUNDIALS 7.1.1
cmake --build build --config Release
```

Then open `radcluster_code/codes/Notebooks/RadCluster_2_1.ipynb`. Details: [Installation](Installation).

## Where to read more

- [README](https://github.com/Ghoniem/RadCluster/blob/main/README.md): capabilities, methodology, repository structure, branches and citation.
- [docs/Formulation](https://github.com/Ghoniem/RadCluster/tree/main/docs/Formulation): the framework manuscript, derivations and study reports, by topic.
- [radcluster_code/CLAUDE.md](https://github.com/Ghoniem/RadCluster/blob/main/radcluster_code/CLAUDE.md): the working physics and solver reference, with every equation.
- [radcluster_code/digital_twin](https://github.com/Ghoniem/RadCluster/tree/main/radcluster_code/digital_twin): the calibration and verification campaign.

## License

MIT License, © 2026 Nasr M. Ghoniem.

Primary reference: N. M. Ghoniem, *A Generalized Graph-Based Cluster Dynamics Framework for Irradiated Materials*, Journal of Nuclear Materials (submitted, 2026). Section and equation numbers in the code and in these pages refer to that document: [Generalized_Cluster_Dynamics.pdf](https://github.com/Ghoniem/RadCluster/blob/main/docs/Formulation/cluster_dynamics_framework/Generalized_Cluster_Dynamics.pdf).
