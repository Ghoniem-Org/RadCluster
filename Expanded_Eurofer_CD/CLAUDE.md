# Expanded_Eurofer_CD — Physics and Solver Reference

## Source Document

All equations cite:
> Ghoniem, N.M. (2026), *"A Cluster Dynamics Model for Radiation Damage Evolution
> in Ferritic-Martensitic Steels"* (Rate_Equations.pdf)

---

## 1. State Vector

### full_CD modes  (Eq. 152, 155, 157)

| Segment | Indices | Variable |
|---|---|---|
| SIA clusters | 0 … N−1 | c_n, n = 1..N |
| Vacancy/He grid | N … N + M·(L_max+1) − 1 | c_{m,ℓ}, m=1..M, ℓ=0..ℓ_max(m) |
| Free He | last | c_h |

Total equations: N + M·(ℓ_max+1) + 1.

**fission** physics option uses He-reduction **Case 2** (Eq. 175, Section 8.2):
decoupled scalar ℓ_tot → reduced system size N + M + 2.

**fusion** physics option uses He-reduction **Case 1** (Eq. 174, Section 8.1):
mean-field ℓ̄(m) per void class → reduced system size N + 2M + 1.

### bin_moment_CD modes  (Chapter 9, Eq. 188–208)

Logarithmic bin partition (Eq. 188–191): bin k spans [n_k, n_{k+1}).
Per bin, two moments tracked (Eq. 192):
  μ_k^(0) = Σ_{n∈bin_k} c_n      (zeroth moment — bin density)
  μ_k^(1) = Σ_{n∈bin_k} n·c_n    (first moment — bin content)

System size ≈ 2·K_bins (Table 33: ~113 bins for N=5000, r=2).

---

## 2. Cascade Production Spectrum  (Eqs. 1–13, Tables 2, 5)

Power-law cluster number distribution:
  ε_m = C · m^{−s}   (Eq. 7)

Normalisation constants (Eq. 9–10):
  C_i = f_i^cl / Σ_{m=2}^{m1} m^{1−s_i}
  C_v = f_v^cl / Σ_{n=2}^{n1} n^{1−s_v}

Production rates (Eq. 152 source term):
  P_n^i = η · G · ε_n^i    [atom frac / s]
  P_m^v = η · G · ε_m^v

| Parameter | Fission | Fusion |
|---|---|---|
| η (survival) | 0.30 | 0.28 |
| f_i^cl | 0.58 | 0.65 |
| s_i | 1.6 | 1.5 |
| m1 | 20 | 50 |
| f_v^cl | 0.15 | 0.20 |
| s_v | 2.5 | 2.3 |
| n1 | 10 | 20 |

---

## 3. Defect Diffusion and Solute Trapping  (Eqs. 17–52)

Point-defect diffusivity (Eq. 17):
  D_α = a² · ν_α · exp(−E_m^α / k_B T)

EUROFER solute-trapping effective frequency (Eq. 42):
  ω_α^eff = ω_α / [1 + Σ_s z_s · c_s · exp(E_b^{s,α} / k_B T)]

This applies to SIA (Eq. 42), vacancy (Eq. 48), and SIA clusters via loop
trapping (Eq. 52).

SIA cluster 1D glide (Eq. 33):
  D_n^{1D} = (3a²ν_0^{1D}) / (2n^{s_1D}) · exp(−E_m^{1D} / k_B T)

Mobility cutoffs: n_max^i = 100 SIAs, m_max^v = 5 vacancies.

---

## 4. Reaction Rate Constants  (Eqs. 109–143, Tables 28, 30)

Geometric prefactors (Eq. 128):
  A_sph  = (48π²)^{1/3} ≈ 7.818   — spherical clusters
  A_loop = 8√(π/√3)    ≈ 10.78   — dislocation loops
  A_1D   = 9/(8π^{2/3}) ≈ 2.632   — pure 1D glide
  B_rot  = (4/π)(8π/3)^{1/3} ≈ 2.627  — 1D/3D crossover

Jump frequency: ω_α = D_α / a²  (Eq. 22–24)

3D spherical capture rate (Eq. 109):
  K_sph(α, m) = A_sph · m^{1/3} · ω_α^eff

Loop capture rate (Eq. 113):
  K_loop(i, n) = A_loop · n^{1/2} · Z_i^loop · ω_i^eff

V–SIA recombination (Eq. 130):
  K_iv = 4√3 · π · ω_v^eff ≈ 21.77 · ω_v^eff

Mixed 1D/3D effective rate (Eq. 121, 141):
  K_n,m^eff = A_sph · m^{1/3} · ω_n^{1D} / (1 + B_rot · L̂² · m^{−1/3})
  where L̂ = L/a (normalised mean free path)

Thermal emission (Eq. 122, 138–140):
  α_α(n) = A_sph · (n−1)^{1/3} · ω_α^eff · exp(−E_b / k_B T)

Trap mutation (Eq. 142):
  Γ_TM(m, ℓ) = ν_0 · exp(−E_TM(m, ℓ) / k_B T)

Radiation re-solution (Eq. 143):
  Γ_res(m, ℓ) = b_0 · ℓ · φ̇

Fixed sink strengths (Eq. 134–137) — dislocation network + grain boundaries + precipitates:
  D_α^d  = Z_α · ρ_d · D_α^eff
  D_α^gb = π² D_α^eff / d_g²
  D_α^p  = Z_p · ρ_p · r_p · D_α^eff

---

## 5. Master Equations  (Eqs. 152, 155, 157)

**SIA clusters** (Eq. 152), n = 1..N:
  dc_n/dt = P_n^i
           + K_sph(i,n−1)·c_1·c_{n−1}  [n−1 grows]
           − K_sph(i,n)·c_1·c_n        [n grows]
           + K_loop(i,n+1)·c_1·c_{n+1} [n+1 emits]
           − K_loop(i,n)·c_1·c_n       [n emits → captured at loops]
           + α_i(n+1)·c_{n+1}          [thermal SIA emission from n+1]
           − α_i(n)·c_n                [thermal SIA emission]
           − K_sph(v,n)·c_v·c_n        [vacancy annihilation]
           + K_sph(v,n+1)·c_v·c_{n+1}  [n+1 shrinks]
           − Σ_m K_n,m^eff·c_n·c_{-m}  [1D glide recombination]
           − D_i^{sink}·c_n            [fixed sinks]

**Vacancy/He clusters** (Eq. 155), m = 1..M, ℓ = 0..ℓ_max(m):
  dc_{m,ℓ}/dt = P_m^v · δ_{ℓ,0}
               + K_sph(v,m−1)·c_v·c_{m−1,ℓ}
               − K_sph(v,m)·c_v·c_{m,ℓ}
               + α_v(m+1,ℓ)·c_{m+1,ℓ}
               − α_v(m,ℓ)·c_{m,ℓ}
               − K_sph(i,m)·c_i·c_{m,ℓ}
               + K_sph(i,m+1)·c_i·c_{m+1,ℓ}
               + K_sph(h,m)·c_h·c_{m,ℓ−1}   [He capture]
               − K_sph(h,m)·c_h·c_{m,ℓ}
               + α_h(m,ℓ+1)·c_{m,ℓ+1}       [He emission]
               − α_h(m,ℓ)·c_{m,ℓ}
               + Γ_TM(m,ℓ)·c_{m,ℓ}→c_{m+1,ℓ−1} [trap mutation]
               + Γ_res·ℓ·c_{m,ℓ}→c_{m,ℓ−1}  [re-solution]
               − D_v^{sink}·c_{m,ℓ}

**Free He** (Eq. 157):
  dc_h/dt = G_He
           − Σ_m Σ_ℓ K_sph(h,m)·c_{m,ℓ}·c_h
           + Σ_m α_h(m,1)·c_{m,1}
           + Σ_m Γ_res(m,1)·c_{m,1}
           − D_h^{sink}·c_h

---

## 6. He-Vacancy State-Space Reduction  (Section 8)

### Case 1 — Fusion (fast equilibration, Eq. 174)
Track marginal c_m^tot = Σ_ℓ c_{m,ℓ} and mean loading ℓ̄_m = Q_m / c_m^tot.
Auxiliary ODE for He content per class Q_m.
System size: N + 2M + 1.

### Case 2 — Fission (decoupled, Eq. 175)
Track only c_m (pure void class) and scalar total He Q_tot = Σ_{m,ℓ} ℓ·c_{m,ℓ}.
He modifies void emission via effective binding energy.
System size: N + M + 2.

---

## 7. Size-Bin Moment Reduction  (Chapter 9, Eqs. 188–211)

Logarithmic bin edges (Eq. 188):
  n_k = floor(n_1 · r^k),  r > 1  (typically r = 2)

Bin moments (Eq. 192):
  μ_k^(0) = Σ_{n∈B_k} c_n
  μ_k^(1) = Σ_{n∈B_k} n·c_n

Bin-integrated RHS (Eq. 193–197):
  dμ_k^(q)/dt = [sources + sinks integrated over B_k]

Piecewise-constant closure (Eq. 198–200):
  c_n ≈ μ_k^(0) / |B_k|  for n ∈ B_k

Hat-function (Galerkin) closure (Eq. 201–206):
  c_n = φ_{k,0}(n)·μ_k^(0) + φ_{k,1}(n)·μ_k^(1)

Inter-bin upwind flux (Eq. 207–208).

Conservation diagnostic (Eq. 211): δ_FP^bin.

---

## 8. Post-Processing  (Eqs. 161–165)

Swelling identity (Eq. 161):
  S(t) = S_I(t) + ΔJ^d(t)

Frenkel pair conservation (Eq. 164):
  δ_FP = |Σ_n n·c_n − Σ_m m·c_m − ΔJ^d| / (G·t)

He conservation (Eq. 165):
  δ_He = |c_h + Σ_{m,ℓ} ℓ·c_{m,ℓ} − G_He·t| / (G_He·t)

---

## 9. Solver Modes

| Mode | Description |
|---|---|
| `cpp_full` | Full system, SUNDIALS CVODE BDF, dense/band/GMRES linear solver |
| `cpp_sliding_win` | Sliding-window SIA truncation, CVODE BDF GMRES |
| `sliding_OpenMP` | Sliding window + OpenMP intra-RHS parallelism |

Physics options (4 per solver mode):
- `full_CD_fission`: Eqs. 152/155/157 with He Case 2 (decoupled)
- `full_CD_fusion`:  Eqs. 152/155/157 with He Case 1 (mean-field)
- `bin_moment_CD_fission`: Chapter 9 bin-moment, fission cascade
- `bin_moment_CD_fusion`:  Chapter 9 bin-moment, fusion cascade

---

## 10. Binding Energies  (Eqs. 62–108, Tables 18–19)

Void capillary (Eq. 66–67):
  E_b^v(m) = E_f_v − A_void · [m^{2/3} − (m−1)^{2/3}]
  A_void = 4π · γ_s · r_0²  (J to eV conversion applied)

Bubble vacancy binding (Eq. 70–73):
  E_b^{bub}(m, ℓ) = E_b^v(m) + He-pressure correction via virial EOS

He virial EOS (Eq. 64–65):
  P_He V = N_He k_B T [1 + B2·(N_He/V) + B3·(N_He/V)²]
  B2 = 1.67×10⁻²⁹ m³/atom,  B3 = 1.84×10⁻⁵⁸ m⁶/atom²

He binding to bubble (Eq. 76–77, Table 19):
  E_b^He(m, ℓ) = E_s^He + P_He(m,ℓ)·Ω − f_blend · A^He(α)·exp(−μ(ℓ−α))

Interstitial loop (Eq. 106–108):
  E_b^{loop}(n) = A_111·n^{−B_111}  (power law, n ≤ n_tr)
  blended to continuum limit via tanh(n − n_tr) / σ_tr

---

## 11. Parameter Sources

| Table | Content |
|---|---|
| Table 2 | Cascade production parameters (fission/fusion) |
| Table 5 | bcc Fe lattice, energetics (updated DFT values) |
| Table 8 | He EOS virial coefficients |
| Table 16 | EUROFER dissolved solute concentrations and trapping energies |
| Table 18 | Atomistic fitting amplitudes A(m) for void correction |
| Table 19 | He binding amplitudes A^He(α) |
| Table 25 | Geometric rate constant prefactors |
| Table 26 | Dislocation bias factors and sink parameters |
| Table 27 | Trap mutation barriers E_TM(m, ℓ) |
| Table 28 | Full rate constant table (3D reactions) |
| Table 29 | Model parameters (solver settings) |
| Table 30 | Rate constant table (1D/mixed reactions) |
