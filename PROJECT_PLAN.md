# Project Plan and Research Status

## Completed in this version

- Clean Python package layout under `src/magnetic_coherent_states/`.
- Runnable main demo through `python run_demo.py`.
- Gauge-covariant finite-difference Hamiltonian using magnetic link variables.
- Coherent-state-like Gaussian initialization.
- Sparse matrix exponential time evolution.
- Diagnostic outputs: density frames, center-of-mass trajectory, participation ratio, CSV, NPZ archive.
- Classical guiding-center toy model.
- hbar sweep experiment.
- Unit tests for core numerical properties.
- Research documentation: paper summary, numerical method, and research report.

## Scientific scope

This repository is a computational prototype. It does not attempt to fully reproduce the rigorous long-time theorem from Boil and Vu Ngoc. The theorem involves semiclassical normal forms, adiabatic projectors, Schwartz-seminorm estimates, and `O(h^infinity)` remainders. Those objects are not directly implemented here.

## Roadmap

### Stage 1 - Working numerical core

Status: completed.

- Build a sparse magnetic Hamiltonian.
- Evolve coherent wave packets.
- Generate basic diagnostics.

### Stage 2 - Semiclassical validation

Status: partially prepared.

- Run hbar sweeps.
- Compare center-of-mass motion across hbar values.
- Track spreading indicators.
- Add grid-convergence checks.

### Stage 3 - Guiding-center comparison

Status: toy model completed, systematic comparison pending.

- Integrate guiding-center dynamics from the same initial point.
- Overlay quantum center-of-mass paths with classical drift curves.
- Measure trajectory mismatch.

### Stage 4 - Toward the paper's splitting phenomenon

Status: future work.

- Build approximate Landau-level projectors.
- Decompose the wave packet into mode components.
- Measure whether components drift at different quantized speeds.
- Test longer time scales of order `1/h` while controlling boundary effects.
