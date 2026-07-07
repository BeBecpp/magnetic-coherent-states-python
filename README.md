# Magnetic Coherent States Python

A reproducible computational research project for exploring coherent-state dynamics under a two-dimensional magnetic Schrodinger flow.

The project is inspired by two mathematical references included in `references/`:

1. Gregory Boil and San Vu Ngoc, **Long-time dynamics of coherent states in strong magnetic fields**.
2. Gregory Boil, Nicolas Raymond, and San Vu Ngoc, **On the stability of the Schwartz class under the magnetic Schrodinger flow**, arXiv:1806.05389.

This repository is not an official implementation of the papers and it does not claim to prove their theorems numerically. Its goal is narrower and concrete: build a clean, runnable, well-documented Python experiment that turns the main objects from the theory into an inspectable numerical model.

## Research question

What does a localized coherent-state-like wave packet do under the magnetic Schrodinger equation

```math
i h \partial_t \psi = (-i h \nabla - A(q))^2 \psi,
```

when the magnetic field is spatially variable?

The reference paper predicts a subtle long-time phenomenon: a coherent state in a strong magnetic field can decompose into multiple coherent components whose guiding-center motions occur at quantized speeds. Directly reproducing that theorem is beyond this small numerical project, but the code here implements the computational primitives needed for a first-stage exploration:

- coherent-state initialization,
- a variable magnetic field,
- a compatible vector potential,
- gauge-covariant finite differences using Peierls phases,
- sparse matrix time evolution,
- density visualization,
- center-of-mass and spreading diagnostics,
- a toy classical guiding-center comparison.

## Magnetic field used in the demo

The default demo uses the magnetic field shown in the motivating figure of the long-time coherent-states paper:

```math
B(x,y)=2-\cos(x)+y^2.
```

A compatible vector potential is

```math
A_1(x,y)=(2-\cos(x))y+\frac{y^3}{3}, \qquad A_2(x,y)=0,
```

so that

```math
\partial_y A_1 - \partial_x A_2 = B.
```

## Numerical method

The discrete Hamiltonian approximates

```math
L_{h,A}=(-ih\nabla-A)^2
```

on a rectangular grid. Instead of expanding derivatives term by term, nearest-neighbor couplings are multiplied by magnetic link variables:

```math
U_x(i+1/2,j)=\exp\left(-i A_x(x_{i+1/2},y_j)\Delta x/h\right),
```

```math
U_y(i,j+1/2)=\exp\left(-i A_y(x_i,y_{j+1/2})\Delta y/h\right).
```

This is the standard gauge-covariant finite-difference idea: magnetic phase information is attached to lattice links. Time evolution is computed with `scipy.sparse.linalg.expm_multiply`, avoiding dense matrix exponentials.

## Repository structure

```text
magnetic-coherent-states-python-research/
├── README.md
├── pyproject.toml
├── requirements.txt
├── PROJECT_PLAN.md
├── paper_summary.md
├── run_demo.py
├── src/
│   └── magnetic_coherent_states/
│       ├── __init__.py
│       ├── core.py
│       └── plotting.py
├── examples/
│   ├── run_demo.py
│   └── classical_guiding_center.py
├── experiments/
│   └── hbar_sweep.py
├── tests/
│   └── test_core.py
├── docs/
│   ├── numerical_method.md
│   ├── paper_summary.md
│   └── research_report.md
├── references/
│   ├── boil-vn.pdf
│   └── 1806.05389v1.pdf
└── outputs/
    └── generated figures and CSV files
```

## Installation

```bash
git clone https://github.com/BeBecpp/magnetic-coherent-states-python.git
cd magnetic-coherent-states-python
python -m venv .venv
```

Activate the virtual environment.

Windows:

```powershell
.venv\Scripts\activate
```

macOS/Linux:

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
pip install -e .
```

## Run the main demo

Fast smoke run:

```bash
python run_demo.py --quick
```

Research-style default run:

```bash
python run_demo.py
```

Custom run:

```bash
python run_demo.py --nx 96 --ny 96 --hbar 0.16 --frames 9 --tmax-scale 1.4 --out-dir outputs
```

The output directory will contain:

- `magnetic_field.png`
- `density_000.png`, `density_001.png`, ...
- `center_trajectory.png`
- `stats.csv`
- `simulation_data.npz`
- `research_summary.md`

## Run the classical guiding-center toy model

```bash
python examples/classical_guiding_center.py
```

This creates `outputs/classical_guiding_center.png`.

## Run an hbar sweep

```bash
python experiments/hbar_sweep.py --quick
```

This compares several values of the semiclassical parameter `hbar` using the same diagnostic pipeline.

## Run tests

```bash
python -m unittest discover -s tests
```

The tests check:

- discrete normalization of coherent states,
- compatibility of the vector potential with the magnetic field,
- Hermitian structure of the discrete Hamiltonian,
- basic evolution mass preservation.

## What this project proves and what it does not prove

This project demonstrates a working numerical pipeline for magnetic coherent-state evolution. It is suitable as a computational research prototype, a portfolio project, and a foundation for deeper experiments.

It does **not** prove the long-time asymptotic splitting theorem. The theorem is a semiclassical result in Schwartz seminorms and involves adiabatic normal forms, Landau-level-like projectors, and remainders of order `O(h^infinity)`. The present code is a finite-grid simulation and should be interpreted as exploratory evidence, not as a mathematical proof.

## Next research extensions

The most valuable next steps are:

1. implement projection diagnostics for approximate adiabatic Landau levels;
2. compare the quantum center motion against guiding-center trajectories more systematically;
3. run hbar-refinement studies and grid-convergence studies;
4. add absorbing boundaries or larger domains for longer-time runs;
5. implement wave-packet decomposition diagnostics to search for quantized-speed separation.
