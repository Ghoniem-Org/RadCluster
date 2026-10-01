# RadCluster_2_1 — Python utilities for cluster dynamics in bcc Fe / EUROFER97.
# Physics: Ghoniem (2026), Rate_Equations.pdf.
# Solver modes: full_system | active_window
# Physics axes: equations ∈ {discrete, bin_moment}  ×  cascade ∈ {fission, fusion}
#   Combined legacy strings: full_CD_fission | full_CD_fusion |
#                            bin_moment_CD_fission | bin_moment_CD_fusion

from .input_data import make_physics_option, split_physics_option  # noqa: F401
