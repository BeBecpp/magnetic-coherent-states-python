# Research Report

## Title

Numerical exploration of coherent-state dynamics under a variable magnetic Schrodinger flow

## Abstract

This project implements a reproducible Python simulation for localized wave-packet evolution under a two-dimensional magnetic Schrodinger operator. Motivated by Boil and Vu Ngoc's work on long-time coherent-state dynamics in strong magnetic fields, the implementation focuses on the numerical construction of a magnetic Hamiltonian, coherent-state-like initial data, sparse time evolution, and diagnostic visualization. The project does not reproduce the paper's full theorem, but it provides a clean computational foundation for investigating the first visible features of magnetic coherent-state propagation.

## Objective

The objective is to convert the mathematical ingredients of the reference papers into a working numerical experiment:

- define a nonconstant magnetic field;
- build a compatible vector potential;
- discretize the magnetic Laplacian in a gauge-covariant way;
- initialize a localized Gaussian wave packet;
- evolve it under the Schrodinger equation;
- measure center motion and spreading.

## Model

The simulation uses

```math
B(x,y)=2-\cos(x)+y^2
```

and

```math
A_1(x,y)=(2-\cos(x))y+y^3/3, \qquad A_2(x,y)=0.
```

This gives `partial_y A_1 - partial_x A_2 = B`.

## Implementation

The project is implemented as a Python package under `src/magnetic_coherent_states`. The main numerical function is `build_magnetic_hamiltonian`, which returns a sparse Hermitian matrix. The simulation is launched through `run_demo.py`.

## Validation checks

The repository includes unit tests for:

- discrete `L2` normalization;
- numerical curl of the vector potential;
- Hermiticity of the Hamiltonian;
- mass preservation under short evolution.

These are engineering validation checks, not theorem proofs.

## Results produced by the demo

A successful run produces:

- a magnetic-field heatmap;
- probability-density frames;
- a center-of-mass trajectory plot;
- a CSV table of time-dependent diagnostics;
- a compressed NPZ file containing the numerical state frames;
- a Markdown summary of the run.

## Interpretation

The resulting figures show that a localized packet evolves nontrivially under a variable magnetic field. The center-of-mass trajectory and participation ratio give a first numerical window into transport and spreading. These results are consistent with the motivation of studying guiding-center effects, but they should not be interpreted as a numerical proof of quantized-speed coherent-state splitting.

## Future work

The strongest next research step is to implement approximate adiabatic Landau-level projections. That would allow the simulation to separate different mode components and test whether their centers drift at different quantized speeds, which is the central phenomenon in the long-time paper.
