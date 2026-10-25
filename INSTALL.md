# UQMaT — Installation Guide

## Requirements

### Mandatory (pure MATLAB, no compilation needed)

| Requirement | Version | Notes |
|---|---|---|
| MATLAB | R2020b or later | Core language features |
| No toolboxes | — | All core functionality uses base MATLAB only |

UQMaT runs on any system with MATLAB R2020b+. No C++ compiler, no toolboxes, no external libraries are needed for the pure MATLAB mode.

---

## Quick start (pure MATLAB)

```matlab
% Option 1: run setup script (recommended)
run('/path/to/UQMaT/uqmat_setup.m')

% Option 2: add to path manually
addpath('/path/to/UQMaT')
clear classes
```

Verify installation:

```matlab
p   = UQMaT.Params('system','chain','N',20,'t',1,'mu',0);
ev  = UQMaT.EV(UQMaT.System(p), 'direct');
[w, D] = ev.dos();
disp('UQMaT working.')
```

---

## Optional: MEX native solvers (C++, for large systems)

MEX solvers provide significant speedups for large systems (N > 500):

| Solver | Speedup | Use case |
|---|---|---|
| `uqmat_eig_mex` | 2–5× | Dense eigensolver |
| `uqmat_green_mex` | 5–20× | Batched Green's function |
| `uqmat_kpm_mex` | 10–50× | KPM/Chebyshev DOS |

### Step 1: Configure a C++ compiler in MATLAB

**Windows (Microsoft Visual C++):**
```matlab
mex -setup C++
% Select: Microsoft Visual C++ 2019 or 2022
% Free: install "Build Tools for Visual Studio" from microsoft.com
```

**macOS:**
```bash
% Install Xcode Command Line Tools (free):
xcode-select --install
```
```matlab
mex -setup C++
```

**Linux:**
```bash
# Ubuntu/Debian
sudo apt-get install g++ build-essential

# CentOS/RHEL
sudo yum install gcc-c++ make
```
```matlab
mex -setup C++
```

### Step 2: Compile MEX solvers

```matlab
% From MATLAB, in the UQMaT directory:
build_mex()
% or
UQMaT.Solver.Native.build_mex()
```

Expected output:
```
Compiler: GNU G++ (or MSVC, Xcode Clang)
  Compiling uqmat_eig.cpp ... OK (uqmat_eig_mex.mexa64)
  Compiling uqmat_green.cpp ... OK (uqmat_green_mex.mexa64)
  Compiling uqmat_kpm.cpp ... OK (uqmat_kpm_mex.mexa64)

Result: 3 compiled, 0 failed.
MEX solvers ready. UQMaT will use them automatically.
```

### Step 3: Verify

```matlab
UQMaT.Config.mex_status()
% Expected output:
%   uqmat_eig_mex:   compiled ✓
%   uqmat_green_mex: compiled ✓
%   uqmat_kpm_mex:   compiled ✓
%   MEX setting: auto

% Run benchmark
benchmark_mex
```

---

## Optional toolboxes

These toolboxes unlock additional features but are **never required**:

| Toolbox | Feature unlocked |
|---|---|
| Parallel Computing Toolbox | `parfor` in disorder ensemble averaging |
| Statistics and Machine Learning Toolbox | `ksdensity` for smoother DOS |
| Signal Processing Toolbox | Advanced spectral filtering |

UQMaT detects these at runtime and uses them if present, falling back to base MATLAB otherwise.

---

## Platform notes

### Windows
- Use backslashes or `fullfile()` for paths
- MSVC compiler recommended; MinGW-w64 also works
- Pre-built `.mexw64` files may be included in future releases

### macOS (Intel and Apple Silicon)
- Xcode Command Line Tools provides `clang++`
- Apple Silicon (M1/M2/M3): use `.mexmaca64` extension
- Intel Mac: use `.mexmaci64` extension

### Linux (HPC clusters)
- Load GCC module: `module load gcc`
- For OpenMP support: `module load gcc/11` or later
- Compile once on the login node, use on compute nodes

---

## Persistent path setup

To avoid running `uqmat_setup.m` every MATLAB session, add to `startup.m`:

```matlab
% Find startup.m location:
edit(fullfile(userpath, 'startup.m'))

% Add this line:
run('/path/to/UQMaT/uqmat_setup.m')
```

---

## Troubleshooting

**"No C++ compiler configured"**
```matlab
mex -setup C++   % select a compiler
```

**"MEX compilation failed"**
```matlab
build_mex('verbose')   % see full compiler output
```

**"Unrecognised function UQMaT.Params"**
```matlab
run('/path/to/UQMaT/uqmat_setup.m')  % add to path
```

**"Invalid use of operator" or other syntax errors**
- Ensure MATLAB R2020b or later
- Run `clear classes` after updating `.m` files

**MEX files not found after compilation**
```matlab
addpath('/path/to/UQMaT/+UQMaT/+Solver/+Native')
UQMaT.Config.mex_status()
```

---

## Uninstall

```matlab
% Remove from path
rmpath('/path/to/UQMaT')

% Remove MEX files (optional)
build_mex('clean')
```
