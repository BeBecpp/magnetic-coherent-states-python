# Simulation Summary

## Configuration

- Grid: `72 x 72`
- Domain: `[-5, 5] x [-5, 5]`
- hbar: `0.16`
- Frames: `9`
- Time interval: `0` to `8.75`

## Main diagnostics

| quantity | initial | final |
|---|---:|---:|
| center x | -1.800000 | -0.372677 |
| center y | 0.350000 | 0.456080 |
| variance x | 0.096800 | 1.969929 |
| variance y | 0.096800 | 0.326967 |
| participation ratio | 1.216425 | 3.190210 |

## Interpretation

The packet begins as a localized Gaussian-like state. Under the magnetic Hamiltonian, the probability density moves and spreads. The center-of-mass trajectory is a finite-grid diagnostic, not a proof of the asymptotic guiding-center theorem. The participation ratio gives a compact scalar proxy for spreading: larger values indicate that the probability density occupies a larger effective area.
