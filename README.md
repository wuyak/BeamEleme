**English** | [简体中文](./README_zh-CN.md)

---

# BeamElem

> Timoshenko Beam Modal Analysis Finite Element Solver

---

## 📋 Dependencies

- **MATLAB** — Core computing environment
- **chebfun** — Characteristic equation root-finding tool (included in project)
  - Website: [https://www.chebfun.org/](https://www.chebfun.org/)
  - Features: Open-source package with ~15 significant digits precision, providing root-finding, integration, and differential equation solving

---

## 📍 Project Overview

### Timoshenko Beam Modal Analysis

Timoshenko beam theory extends Euler-Bernoulli theory by adding shear deformation and rotary inertia terms, making it more suitable for short thick beams and higher-order modal analysis. Standard finite element methods suffer from shear locking when solving Timoshenko beams: elements become overly stiff in the slender beam limit, leading to non-convergent results. While analytical solutions can provide exact frequencies, the characteristic equations vary significantly across boundary conditions, and numerical root-finding is prone to missing or spurious roots.

### Purpose of This Project

This project implements 1D Timoshenko beam finite element numerical solutions, analytical solution computation, and error verification, providing modal frequencies, mode shapes, and visualization results. It supports 10 boundary conditions (PP/CF/CC/FF/CP/PR/CR/RR/RF/PF) and 3 vibration types (axial/torsional/bending-shear), with built-in material library and section property calculation, automatically outputting frequency errors, MAC values, and mode shape plots.

**Workflow**:

1. **Parametric Modeling** — Select properties from material library, calculate section geometric parameters
2. **FEM Solution** — Construct stiffness and mass matrices, apply boundary constraints, solve eigenvalue problem
3. **Analytical Solution** — Compute theoretical frequencies and mode shapes for three vibration types
4. **Automatic Verification** — Compare frequency errors and MAC values, automatically match modes
5. **Result Caching** — Generate unique identifier based on parameters to avoid redundant computation
6. **Visualization** — Load cached data and plot mode shapes

**Extension Experiments**: Adding new materials or section types only requires adding configurations in the `+parameters` module. Extending to coupled models requires modifying `+matrix`, `+solvers`, and other modules. See [`docs/project_structure.md`](docs/project_structure.md) for details.

## 🚀 Quick Start

**For detailed usage, parameter configuration, and notes, refer to**: [`example_basic.m`](example_basic.m)

```matlab
% Add project paths
addpath(genpath('src'));
addpath(genpath('chebfun'));

% Define material
mat = parameters.MaterialLibrary.Steel();

% Define section (rectangular: 20cm width × 10cm height)
sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);

% Create beam model (length 1m)
beam = parameters.BeamModel(mat, sec, 'L', 1.0);

% Solve (pinned-pinned beam, 100 elements, 10 modes)
wf_result = workflow.solve('PP', beam, 100, 10);

% Visualize 1st bending mode
fem_result = wf_result.bending_shear.fem_result;
visualization.plotMode(fem_result, 1);
```

**Output Example**:

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【Axial】Axial Vibration
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Mode   FEM (Hz)   Analytical (Hz)   Error     MAC
  1     3184.71      3184.71       0.00%    1.0000
  2     6369.43      6369.43       0.00%    1.0000
  3     9554.14      9554.14       0.00%    1.0000
  ...

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【Bending-Shear】Bending-Shear Vibration
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Mode   FEM (Hz)   Analytical (Hz)   Error     MAC
  1      837.42       838.15       +0.09%   0.9998
  2     2149.33      2153.26       +0.18%   0.9997
  3     3765.98      3773.01       +0.19%   0.9997
  ...

======== Error Statistics ========
【Axial Vibration】
  Max Frequency Error: 0.0000%
  Average MAC: 1.000000

【Torsional Vibration】
  Max Frequency Error: 0.0000%
  Average MAC: 1.000000

【Bending-Shear Vibration】
  Max Frequency Error: 0.1832%
  Average MAC: 0.999712

💾 Results cached to: .cache/4BB0F885_PP_N100_M10/
```

### Mode Shape Visualization

![Mode Shapes](PP_N100_mode_1_2_bending_shear.png)

*Bending-shear mode shapes (Mode 1-2, PP boundary): FEM results (blue) vs. Analytical solution (red dashed)*

### Cache Structure

**Cache Path Format**: `.cache/<hash>_<BC>_N<NElem>_M<n_modes>/`

**Parameter Description**:
- `<hash>`: Hash of beam model parameters (material+section+length), example: `4BB0F885`
- `<BC>`: Boundary condition, example: `PP` (pinned-pinned)
- `N<NElem>`: Number of elements, example: `N100` (100 elements)
- `M<n_modes>`: Number of modes, example: `M10` (10 modes)

**Cache Contents**:

```text
.cache/4BB0F885_PP_N100_M10/
├── metadata.mat                        # Beam model parameters and solver config
├── axial/                              # Axial vibration results
│   ├── fem_result.mat                  # FEM results (frequencies, modes)
│   ├── ana_result.mat                  # Analytical solution results
│   └── comp_result.mat                 # Comparison results (errors, MAC values)
├── torsion/                            # Torsional vibration results
│   └── (same three files)
└── bending_shear/                      # Bending-shear vibration results
    └── (same three files)
```

### Convergence Analysis

The FEM solver demonstrates excellent convergence across all boundary conditions:

![Convergence](src/+workflow/tests/convergence/output/convergence_rectangular_0.2x0.1.png)

*Convergence test for rectangular beam (0.2m × 0.1m, 50 modes): Frequency error decreases with increasing element count. All 10 boundary conditions show consistent second-order convergence.*

---

## 📁 Project Structure

```text
BeamElem/
├── example_basic.m                     # Basic usage example
├── chebfun/                            # chebfun open-source library
├── .cache/                             # Computation cache
├── src/                                # Core source code
│   ├── +workflow/                      # Solution workflow control (with convergence tests)
│   ├── +parameters/                    # Parameter management
│   ├── +solvers/                       # Finite element solver
│   ├── +matrix/                        # Element matrix generation (stiffness/mass)
│   ├── +analytical/                    # Analytical solution module (axial/torsional/bending-shear)
│   ├── +comparison/                    # FEM vs analytical comparison
│   ├── +boundary_conditions/           # Boundary condition management
│   └── +visualization/                 # Visualization tools (with test examples)
├── reference/                          # Theoretical references
│   └── REFERENCES.md
└── docs/                               # Project documentation
    ├── code_corrections.md
    ├── analytical_solution_limitations.md
    └── project_structure.md            # Detailed directory structure
```

> 📖 For detailed directory structure, see: [`docs/project_structure.md`](docs/project_structure.md)

---

## 📚 Theoretical Foundations

The project is based on the following core references:
- **Friedman & Kosmatka (1993)** — Improved elements avoiding shear locking
- **Khasawneh & Segalman (2019)** — Analytical solutions for bending-shear vibration
- **Blevins (1979)** — Frequency formulas for axial and torsional vibration
- **Roark's Formulas** — Section geometric properties
- **Hutchinson (2001)** — Shear correction coefficients

> 📖 For all references and code citations, see: [`reference/REFERENCES.md`](reference/REFERENCES.md)
