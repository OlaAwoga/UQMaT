# UQMaT — Unified Quantum Materials Toolbox

**Version:** 1.0.0 | **Language:** MATLAB R2020b+ | **License:** MIT

A modular MATLAB toolkit for building, solving, and analysing tight-binding
and Bogoliubov-de Gennes (BdG) Hamiltonians on 1D, 2D, and 3D lattices.
Supports Hermitian and non-Hermitian systems, heterostructures, quasiclassical
Green's function methods, and optional MEX C++ native solvers.

---

## Setup

```matlab
% Run once per session (or add to startup.m):
run('/path/to/UQMaT/uqmat_setup.m')

% Or from within examples/:
uqmat_path()
```

---

## Quick start

```matlab
% Solve once — share everywhere
p   = UQMaT.Params('system','chain','N',60,'t',1,'mu',0);
sys = UQMaT.System(p);
ev  = UQMaT.EV(sys, 'direct');        % choose solver explicitly
[w, D] = ev.dos('range',[-3,3]);

% s-wave SC
p  = UQMaT.Params('system','chain','N',60, ...
     'SC',struct('Delta_s',0.3,'sym','onsite'));
ev = UQMaT.EV(UQMaT.System(p), 'direct');

% Topological SC — Majorana polarization
p   = UQMaT.Params('system','kitaev','N',80,'mu',-1,'Delta_p',0.5,'bc','obc');
ev  = UQMaT.EV(UQMaT.System(p), 'direct');
idx = ev.edge_idx();
mp  = ev.mp();
[MP_L, MP_R, MP_LR] = mp.left_right(idx);   % MP_LR < 0 → topological

% 2D shorthand
p  = UQMaT.Params2D('Nx',20,'Ny',20,'SC',struct('Delta_s',0.3,'sym','dx2y2'));
ev = p.solve();

% Green's function solvers
p  = UQMaT.ParamsResolventGreen('system','chain','N',200, ...
     'EnergyGrid',-3:0.01:3,'eta',0.01);
ev = p.solve();

% 1D wire on 2D substrate
p  = UQMaT.Params('system','chain1d_square2d', ...
     'hetero',struct('substrate',struct('Nx',20,'Ny',20,'t',1,'mu',0, ...
         'SC',struct('Delta_s',0.35,'sym','onsite')), ...
         'chain',struct('N',20,'t',0.8,'mu',-0.5,'B',[0,0,0.5],'lambda_R',0.5), ...
         'Gamma',0.3,'direction','x','chain_row',10));
ev = UQMaT.EV(UQMaT.System(p), 'direct');
```

---

## Package structure

```
UQMaT/
├── +UQMaT/
│   ├── Params.m                        ← all parameters
│   ├── Params2D.m                      ← 2D shorthand
│   ├── ParamsResolventGreen.m          ← resolvent GF params
│   ├── ParamsRecursiveGreen.m          ← recursive GF params
│   ├── ParamsKKRGreen.m                ← KKR GF params
│   ├── ParamsChebyshev.m               ← KPM params
│   ├── System.m                        ← routes Params → Model → H
│   ├── EV.m                            ← single solve, all observables share
│   ├── Config.m                        ← MEX on/off/auto, verbosity
│   ├── Observables.m                   ← legacy compatibility
│   │
│   ├── +SC/BdG.m
│   ├── +Models/
│   │   ├── +Chains1D/ (Chain, SSH, Kitaev, HNKitaev, ChainBC)
│   │   ├── +Lattice/  (Regular, SSH, Kitaev, Nambu, Square2D)
│   │   ├── +Lattice2D/ (+Square, +Kagome, +TMD)
│   │   ├── +Lattice3D/ (+Cubic, +Diamond)
│   │   ├── +Heterostructure/ (Junction1D, Chain1D_Square2D)
│   │   ├── +Continuum/ (FreeElectron, BdG)
│   │   └── +Quasiclassical/ (Eilenberger, Usadel, Keldysh)
│   │
│   ├── +Solver/
│   │   ├── +Diagonalization/ (Direct, Lanczos)
│   │   ├── +Greens/ (Resolvent, Recursive, KKR)
│   │   ├── +Spectral/ (Chebyshev)
│   │   └── +Native/ (uqmat_eig.cpp, uqmat_green.cpp, uqmat_kpm.cpp, build_mex.m)
│   │
│   ├── +Observables/
│   │   ├── +Conventional/ (Spectral, Nambu, Disorder, SpinResolved)
│   │   ├── +SC/Spectral.m
│   │   ├── +Spectral/DOS.m
│   │   └── +Topological/
│   │       ├── +Bulk/Invariants.m
│   │       ├── +RealSpace/ (Invariants, MajoranaPolarization)
│   │       └── PhaseDiagram.m
│   │
│   ├── +OrderParameterSymmetry/ (OPS, Chain1D, Square2D, Kagome, Honeycomb)
│   ├── +Disorder/Disorder.m
│   └── +Plotter/ (Spectrum, DOS, LDOS, Bands, Phase, SpinTexture, Lattice, utils)
│
├── examples/
│   ├── uqmat_path.m
│   ├── ex01_chain_spectrum.m
│   ├── ex02_topological_sc.m
│   ├── ex03_ssh_chain.m
│   ├── ex04_kitaev_chain.m
│   ├── ex05_square_lattice.m
│   ├── ex06_inhomogeneous_delta.m
│   ├── ex07_chain1d_square2d.m
│   ├── ex08_greens_resolvent.m
│   ├── ex09_greens_recursive.m
│   └── ex10_greens_kkr.m
│
├── examplesMEX/
│   ├── mex_path.m
│   ├── ex01_mex.m … ex07_mex.m
│   └── benchmark_mex.m
│
├── uqmat_setup.m
├── build_mex.m
├── INSTALL.md
├── CHANGELOG.md
├── CITATION.cff
├── LICENSE
└── README.md
```

---

## `UQMaT.EV` — solve once, use everywhere

```matlab
% Choose solver
ev = UQMaT.EV(sys, 'direct')                        % dense eig(), all states
ev = UQMaT.EV(sys, 'lanczos', 'k',10, 'sigma',0)   % sparse, k states near E=0
ev = UQMaT.EV(sys, 'resolvent', 'eta',0.02)         % G(ω) via inversion
ev = UQMaT.EV(sys, 'chebyshev', 'M',512)            % KPM, no eig, very large N
ev = UQMaT.EV(V, E, N, p)                           % supply own eigenvectors

% Energy grid — all forms accepted
ev.dos(-2:0.005:2)
ev.dos('range',[-2,2])
ev.dos('range',[-2,2],'Npts',300)
ev.dos('eta',0.02)
ev.dos()

% Observables
[E, V]          = ev.spectrum()
idx             = ev.edge_idx()
[w, D]          = ev.dos()
[w, D_up, D_dn] = ev.dos_spin()
[w, D_e, D_h]   = ev.dos_ph()
[w, L]          = ev.ldos(sites)
[w, Sx, Sy, Sz] = ev.spin_texture(sites)
F               = ev.pair_amplitude()
obs             = ev.spectral()
inv             = ev.realspace_invariants()
mp              = ev.mp()
```

### Solver reference

| Name | Class | Best for |
|---|---|---|
| `'direct'` | `Diagonalization.Direct` | N ≲ 2000, all states |
| `'lanczos'` | `Diagonalization.Lanczos` | Large N, few states near E=σ |
| `'resolvent'` | `Greens.Resolvent` | Exact G(ω), O(N³) per ω |
| `'recursive'` | `Greens.Recursive` | Large N, transport T(ω) |
| `'kkr'` | `Greens.KKR` | t-matrix, CPA disorder |
| `'chebyshev'` | `Spectral.Chebyshev` | Very large N, KPM DOS |

---

## Green's function Params classes

```matlab
% Resolvent: G(ω) = (ω+iη−H)⁻¹
p  = UQMaT.ParamsResolventGreen('system','chain','N',200, ...
     'EnergyGrid',-3:0.01:3,'eta',0.01)
ev = p.solve()

% Recursive: layer-by-layer, O(N·M²)
p  = UQMaT.ParamsRecursiveGreen('system','chain','N',1000, ...
     'EnergyGrid',-3:0.01:3,'eta',0.01,'n_layers',1000)
[ev, slv] = p.solve()
[w, T]    = slv.transmission()   % Landauer T(ω)

% KKR: t-matrix + structure constants
p  = UQMaT.ParamsKKRGreen('system','square','Nx',10,'Ny',10, ...
     'EnergyGrid',-3:0.01:3,'eta',0.02,'use_cpa',true)
[ev, slv] = p.solve()
[w, Sig]  = slv.self_energy()    % CPA self-energy

% KPM Chebyshev: O(N·M·Nrand), no eig
p  = UQMaT.ParamsChebyshev('system','square','Nx',50,'Ny',50, ...
     'M',1024,'Nrand',20,'EnergyGrid',-4:0.005:4)
ev = p.solve()
```

---

## `UQMaT.Params2D` — 2D shorthand

```matlab
p  = UQMaT.Params2D('Nx',20,'Ny',20,'t',1,'mu',0, ...
     'SC',struct('Delta_s',0.3,'sym','dx2y2'))
ev = p.solve()
```

---

## Systems and H sizes

| `system` | Normal | Spinful | BdG |
|---|---|---|---|
| `'chain'` (N) | N×N | 2N×2N | 4N×4N |
| `'ssh'` (2N) | 2N×2N | 4N×4N | 8N×8N |
| `'kitaev'` (N) | — | — | 2N×2N |
| `'hnkitaev'` (N) | — | — | 2N×2N (NH) |
| `'square'` (Nx·Ny) | Ns×Ns | 2Ns×2Ns | 4Ns×4Ns |
| `'kagome'` (3·Nx·Ny) | Ns×Ns | 2Ns×2Ns | 4Ns×4Ns |
| `'tmd'` (3·Nx·Ny) | Ns×Ns | 2Ns×2Ns | 4Ns×4Ns |
| `'cubic'` (Nx·Ny·Nz) | Ns×Ns | 2Ns×2Ns | 4Ns×4Ns |
| `'diamond'` (2·Nx·Ny·Nz) | Ns×Ns | 2Ns×2Ns | 4Ns×4Ns |
| `'junction1d'` | via `hetero` | | |
| `'chain1d_square2d'` | via `hetero` | | |

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `system` | `'chain'` | system type |
| `N` | 50 | 1D sites |
| `Nx`, `Ny`, `Nz` | 10 | 2D/3D sites |
| `bc`, `bcx`, `bcy`, `bcz` | `'obc'` | boundary conditions |
| `t` | 1.0 | hopping (scalar or vector) |
| `mu` | 0.0 | chemical potential |
| `B` | `[0,0,0]` | Zeeman field |
| `lambda_R` | 0.0 | Rashba SOC |
| `SC` | `struct()` | `Delta_s` (scalar or Ns-vector), `sym`, `Delta_t`, `sym_t` |
| `AM` | `struct()` | altermagnetism |
| `CDW` | `struct()` | charge density wave |
| `Delta_p` | 0.0 | Kitaev p-wave |
| `t_nh` | 0.0 | non-Hermitian asymmetry |
| `W` | 0.0 | Anderson disorder |
| `eta` | 0.01 | Lorentzian broadening |
| `EnergyGrid` | `[]` | energy grid |

---

## Order parameter symmetry

```matlab
% Site-resolved SC (proximity / vortex)
Delta_vec = zeros(N,1); Delta_vec(N/2+1:end) = 0.4;
'SC', struct('Delta_s',Delta_vec,'sym','onsite')

% 2D vortex — flatten Ny×Nx to Ns×1
'SC', struct('Delta_s',Ds_vortex(:),'sym','onsite')
```

| Lattice | Symmetry labels |
|---|---|
| 1D | `onsite`, `px`, `extended_s`, `dx2` |
| Square 2D | `onsite`, `extended_s`, `px`, `py`, `chiral_p`, `dx2y2`, `dxy` |
| Kagome | `onsite`, `extended_s`, `px`, `py`, `chiral_p`, `dx2y2`, `d_chiral` |

---

## Majorana polarization

```matlab
idx = ev.edge_idx();
mp  = ev.mp();

C                     = mp.charge(idx)           % |u|² - |v|²
MP                    = mp.local(idx)            % 2Re[Σ_σ u_σ v_σ*]
[MP_L, MP_R, MP_LR]  = mp.left_right(idx)       % MP_LR < 0 → topological
[P, M]               = mp.total(idx)             % M→0: Majorana, M→1: fermionic
MP_nl                = mp.nonlocal(idx)          % N×N matrix
[NMP, MP_LR_nl]      = mp.nonlocal_total(idx)   % cross-block product
[prof, x]            = mp.profile()
disp(mp)
```

---

## Topological invariants

```matlab
% Bulk (PBC)
inv = UQMaT.Observables.Topological.Bulk.Invariants(sys)
C   = inv.chern()           % Chern number
phi = inv.zak()             % Zak phase
nu  = inv.winding()         % winding number
Z2  = inv.Z2_2d()           % Z2 (class AII)
disp(inv)

% Real-space (OBC)
inv = ev.realspace_invariants()
[E_edge, idx] = inv.edge_states()
mp            = inv.majorana_polarization()
```

---

## MEX native solvers (optional)

```matlab
% Compile once (requires C++ compiler — see INSTALL.md):
build_mex()

% Status and control
UQMaT.Config.mex_status()
UQMaT.Config.use_mex('auto')    % default: use if compiled
UQMaT.Config.use_mex(true)      % require MEX
UQMaT.Config.use_mex(false)     % MATLAB only

% Benchmark
cd examplesMEX
benchmark_mex
```

| MEX solver | Speedup | Replaces |
|---|---|---|
| `uqmat_eig_mex` | 2–5× | `eig()` in Direct solver |
| `uqmat_green_mex` | 5–20× | MATLAB loop in Resolvent |
| `uqmat_kpm_mex` | 10–50× | MATLAB loop in Chebyshev |

---

## Examples

| File | Demonstrates |
|---|---|
| `ex01_chain_spectrum.m` | Normal vs SC chain; DOS, LDOS |
| `ex02_topological_sc.m` | Rashba+Zeeman SC; MZM; Majorana profile |
| `ex03_ssh_chain.m` | SSH; Zak phase; phase diagram |
| `ex04_kitaev_chain.m` | Kitaev; full MP workflow (C, MP, MP_LR, nonlocal) |
| `ex05_square_lattice.m` | 2D square; d-wave SC; cylinder bands |
| `ex06_inhomogeneous_delta.m` | Site-resolved Δ(i); proximity effect; vortex core |
| `ex07_chain1d_square2d.m` | 1D wire on 2D substrate; proximity; MZM; Γ scan |
| `ex08_greens_resolvent.m` | Resolvent G(ω): DOS, LDOS, spectral function |
| `ex09_greens_recursive.m` | Recursive GF: large N, transmission T(ω), disorder |
| `ex10_greens_kkr.m` | KKR: t-matrix, LDOS, CPA self-energy |
| `examplesMEX/ex01_mex.m`–`ex07_mex.m` | Same examples with MEX timing comparison |
| `examplesMEX/benchmark_mex.m` | MEX vs MATLAB scaling benchmark |

---

## Citation

```bibtex
@software{UQMaT,
  title   = {{UQMaT}: Unified Quantum Materials Toolbox},
  version = {1.0.0},
  url     = {https://github.com/yourname/UQMaT}
}
```

Key references implemented:

| Feature | Reference |
|---|---|
| BdG formalism | Zhu J.-X., Springer (2016) |
| Eilenberger | Eilenberger G., Z. Phys. **214**, 195 (1968) |
| Usadel | Usadel K., PRL **25**, 507 (1970) |
| Keldysh | Keldysh L.V., JETP **20**, 1018 (1965) |
| TMD model | Liu G.-B. et al., PRB **88**, 085433 (2013) |
| Chern number | Fukui T. et al., JPSJ **74**, 1674 (2005) |
| Berry curvature | Berry M.V., Proc. R. Soc. A **392**, 45 (1984) |
| Wilson loop | Soluyanov A.A., Vanderbilt D., PRB **83**, 235401 (2011) |
| Entanglement spectrum | Li H., Haldane F.D.M., PRL **101**, 010504 (2008) |
| Z2 invariant | Fu L., Kane C.L., PRB **76**, 045302 (2007) |
| Hatano-Nelson | Hatano N., Nelson D.R., PRL **77**, 570 (1996) |
| Non-Hermitian topology | Bergholtz E.J. et al., Rev. Mod. Phys. **93**, 015005 (2021) |
| Majorana polarization | Sticlet D. et al., PRL **108**, 096802 (2012) |
| Nonlocal MP | Sedlmayr N. et al., PRB **93**, 245109 (2016) |
| MP_LR | Awoga O.A., Cayao J., PRB (2023) |
| MP_LR (2026) | Kumar H., Cayao J., Awoga O.A., PRB (2026) |
| RGF transport | López Sancho M.P. et al., J. Phys. F **14**, 1205 (1984) |
| KKR | Korringa J., Physica **13**, 392 (1947) |
| KPM | Weisse A. et al., Rev. Mod. Phys. **78**, 275 (2006) |
| LAPACK | Anderson E. et al., LAPACK Users' Guide, SIAM (1999) |

---

## License

MIT — see `LICENSE`.
