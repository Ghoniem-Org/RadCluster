# docs/Formulation — index

Derivations, formulation notes, paper sections and study reports for RadCluster,
grouped by topic. LaTeX sources sit next to their compiled PDFs and `.bib` files;
compile each document from inside its own folder.

| Folder | Topic | Contents |
|---|---|---|
| `cluster_dynamics_framework/` | The generalized graph-based cluster-dynamics framework and its paper drafts | `Generalized_Cluster_Dynamics.pdf` (main paper), `Defect_Cluster_Dynamics.pdf` + `Abstract_and_Introduction.tex/.bib`, `generalized_defect_cluster_dynamics_section.tex`, `simulation_methodology.tex`, novelty analyses (`radcluster_novelty*.pdf/.tex`) |
| `reaction_admissibility_graph/` | The RAG: implementation, architecture, solver coupling | `supplement_S1/` (Supplementary S1: CRN lineage, NetworkX use, three-layer architecture, solver–RAG relation; drop-in `.tex`, TikZ figures, `.bib`, preview PDF); `figures/` (EUROFER-97 RAG figures, GraphML exports, structural report — written by `radcluster_code/codes/make_rag_figure.py` and `make_rag_hyperedge_figure.py`) |
| `reaction_kernels/` | Rate kernels and thermodynamic consistency | 1D-migration reaction frequencies (`one_dimensional_reaction_kernels.tex`, `Reaction Kernels.pdf`); detailed balance for effective kernels (`detailed_balance_effective_kernels.tex`, `detailed balance.pdf`) |
| `state_space_reduction_and_solvers/` | Bin-moment reduction and integrator development | `bin_moments.tex`, `Cuda_C++_Plan.md` (GPU integrator plan) |
| `dislocation_loops/` | ½⟨111⟩→⟨100⟩ loop conversion and loop→network loss | `loop_111_to_100_conversion.tex/.bib`, `radcluster_two_channel_loop_conversion.tex/.pdf`, `derivation_loop_conversion_kernels.tex/.pdf`, `loop_character_analysis.tex` + `_1.pdf`, `loop_network_loss.tex/.bib` |
| `segregation_and_precipitation/` | Radiation-induced segregation and solute precipitation (RadCluster_2_1 objective b) | `radcluster_2_1_RIS_plan.tex/.bib` |
| `steel_applications_and_calibration/` | F/M-steel application, digital twin, calibration | `Ferritic_Martensitic_Steel_Microstructure.pdf`, `radcluster_digital_twin_appendix.tex`, `RadCluster_AI_Calibration_Adaptive_Solver_Roadmap.md` |
| `verification_campaign/` | Verification and approximation study for the paper revision | `paper_revision_verification_plan.md`, `campaign_results.tex/.pdf`, `campaign_M_results.md`, `verification_results.md`, `d111_investigation.md`, `fill_c0_row.py`; `figs/` and `tables/` (written by `radcluster_code/codes/make_size_effect_figures.py` and `digital_twin/verification/make_M_*.py`) |

Generator scripts in `radcluster_code/` point at these sub-folders; if you move a
folder, update the `OUT` / `--out` paths in those scripts as well.
