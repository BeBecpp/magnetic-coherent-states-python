# Numerical Method

## Discrete domain

The wavefunction is represented on a uniform rectangular grid

```math
(x_i,y_j) \in [x_{min},x_{max}]\times[y_{min},y_{max}].
```

The grid uses a Dirichlet exterior approximation: values outside the computational box are treated as zero. This is simple and reproducible, but long-time simulations can be affected by boundary reflections. Larger boxes or absorbing layers are recommended for more serious long-time studies.

## Magnetic Hamiltonian

The continuous model is

```math
L_{h,A}=(-ih\nabla-A)^2.
```

The finite-difference Hamiltonian is assembled as a sparse matrix. The diagonal is the usual second-difference diagonal, while off-diagonal nearest-neighbor entries carry Peierls phases.

For the x-direction:

```math
H_{(i,j),(i+1,j)} = -\frac{h^2}{\Delta x^2}\exp\left(-i A_x(x_{i+1/2},y_j)\Delta x/h\right).
```

The backward link is the complex conjugate, which makes the matrix Hermitian up to floating-point tolerance. The y-direction is analogous.

## Time evolution

For each requested time `t`, the code computes

```math
\psi(t)=\exp(-iHt/h)\psi(0)
```

using `scipy.sparse.linalg.expm_multiply`. This avoids forming dense matrix exponentials and is suitable for moderate grid sizes.

## Diagnostics

The code records:

- probability density `|psi|^2`,
- center of mass,
- coordinate variances,
- participation ratio,
- a compressed NPZ archive of frames.

The participation ratio is

```math
PR(\psi)=\frac{1}{\int |\psi|^4 dxdy}.
```

It is a simple scalar indicator of spreading. It is not a substitute for the Schwartz seminorms used in the theorem.

## Numerical caveats

1. Smaller `hbar` values are more semiclassical but require finer grids.
2. Long times of order `1/hbar` are expensive and boundary-sensitive.
3. The finite-grid model cannot prove `O(hbar^infinity)` remainder estimates.
4. The default coherent-state initialization is a practical Gaussian packet, not a full symbolic construction from the paper.
