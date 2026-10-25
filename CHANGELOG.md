# UQMaT — Changelog

---

## v1.0.0 — First stable release

### Summary
First public release of UQMaT. All features tested and documented.
Package runs on any MATLAB R2020b+ installation with no external dependencies.
Optional MEX C++ solvers provide 2–50× speedup for large systems.

### Core features
- **`UQMaT.Params`** — unified parameter class for all system types
- **`UQMaT.Params2D`** — 2D lattice shorthand with `toParams()` and `solve()`
- **`UQMaT.System`** — routes Params to the correct Hamiltonian builder
- **`UQMaT.EV`** — single-solve cache; all observables share E and V; explicit solver selection
- **`UQMaT.Config`** — global MEX on/off/auto control

### Models
- 1D: chain, SSH, Kitaev, HNKitaev (non-Hermitian), ChainBC
- 2D: square, Kagome, TMD (MoS₂/WS₂/MoSe₂/WSe₂)
- 3D: cubic (anisotropic), diamond (FCC)
- Heterostructures: NS/SS' junction, 1D wire on 2D substrate
- Continuum: FreeElectron, BdG
- Quasiclassical: Eilenberger, Usadel, Keldysh

### Solvers
- `'direct'` — dense `eig()` (MEX-accelerated if compiled)
- `'lanczos'` — sparse `eigs()`
- `'resolvent'` — G(ω) = (ω+iη−H)⁻¹
- `'recursive'` — layer-by-layer RGF, includes transmission T(ω)
- `'kkr'` — KKR t-matrix + structure constants, CPA self-energy
- `'chebyshev'` — KPM/Chebyshev, O(N·M·Nrand)

### Green's function Params classes
- `ParamsResolventGreen`, `ParamsRecursiveGreen`, `ParamsKKRGreen`, `ParamsChebyshev`

### Observables
- DOS, LDOS (spin-resolved, particle-hole)
- Spin texture, pair amplitude, probability density
- Self-consistent gap equation
- Bulk topological: Chern, Berry curvature, Wilson loop, entanglement spectrum, Z2, quantum metric, winding, Zak, Pfaffian
- Real-space: edge states, IPR, pair amplitude

### Majorana polarization
- Charge C(i) = |u|² − |v|²
- Local MP(i) = 2 Re[Σ_σ u_σ(i) v_σ*(i)]
- Left/right: MP_LR = MP_L · MP_R (real, signed)
- Nonlocal: MP_nl(i,j) = 2 Re[Σ_σ u_σ(i) v_σ*(j)]
- Ref: Sticlet et al. PRL 108, 096802 (2012); Awoga & Cayao PRB (2023)

### Order parameter symmetry
- SC (singlet + triplet), AM (altermagnet), CDW
- Site-resolved Δ(i) vector for inhomogeneous/proximity SC and vortex cores
- Symmetry labels: onsite, px/py, dx2y2, dxy, chiral_p, extended_s, gwave

### MEX native solvers (optional C++)
- `uqmat_eig.cpp` — LAPACK eigensolver (2–5× speedup)
- `uqmat_green.cpp` — batched GF inversion (5–20× speedup)
- `uqmat_kpm.cpp` — BLAS Chebyshev KPM (10–50× speedup)
- `build_mex()` — one-command compiler
- Auto-fallback to MATLAB if MEX not compiled

### Examples (examples/)
- `ex01` – `ex07`: chain, SSH, Kitaev, SC, 2D lattice, vortex, wire-on-substrate
- `ex08` – `ex10`: Green's function methods (resolvent, recursive, KKR)

### examplesMEX/
- `ex01_mex` – `ex07_mex`: MEX vs MATLAB timing comparison
- `benchmark_mex`: scaling benchmark across system sizes

### Bug fixes (pre-release)
- `EV` made handle class (lazy observable cache now persists)
- `EV` constructor ternary `? :` operator replaced
- `Params.sc` property restored; `am_components cell` validator fixed
- `Spectral.ldos`, `ldos_ph`, `pair_amplitude`: Hdim-aware strides for mixed heterostructures
- `ChainBC`: `zeros(1,3)(k)` invalid indexing fixed; `dvec_from_sym` added
- `RealSpace.Invariants.fromEV`: crash when `obj.sys=[]` fixed
- `Chain1D_Square2D.buildH`: `bfS`/`bfC` derived from actual H dimensions
- `Chain1D_Square2D`: `pS.system='square'` silent failure fixed
- `Plotter`: `RdBu`/`viridis`/`bwr` → `resolve_colormap()` (works without toolboxes)
- `Plotter`: `clim()` → `caxis()` for R2020b compatibility
- `uqmat_path`/`mex_path`: `clear classes` removed from functions (wiped return value)
- Path resolution: walk-up search for `+UQMaT/` handles double-nested zip extraction
- `Recursive.n_layers`: `{mustBeInteger}` only (allows 0 = auto)

---

## Pre-release development (v2.x internal)

### v2.7.0
- MEX C++ native solvers, Recursive GF, KKR GF, GF Params classes, ex08–10

### v2.6.0
- Solver selection interface, vector `Delta_s`, ex07 (1D wire on 2D substrate)

### v2.5.0
- `EV` single-solve design, HNKitaev, Chern/Berry/Wilson/Z2/quantum metric, Majorana polarization

### v2.4.0
- Heterostructures, cubic/diamond/TMD, quasiclassical (Eilenberger/Usadel/Keldysh)

### v2.3.0
- Kagome lattice, OPS framework (`SC`/`AM`/`CDW` struct interface), altermagnetism

### v2.2.0
- `B=[Bx,By,Bz]` vector, `EnergyGrid`, site-resolved `Delta_s`

### v2.1.0
- Unified observables, Plotter, topological invariants (winding, Zak, Pfaffian)

### v2.0.0
- Solver hierarchy, BdG assembly, 2D square lattice

### v1.x (initial)
- 1D chains (chain, SSH, Kitaev), DOS, LDOS, spectrum
