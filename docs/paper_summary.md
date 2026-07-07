# Paper Summary and Computational Translation

## 1. Long-time dynamics of coherent states in strong magnetic fields

**Authors:** Gregory Boil and San Vu Ngoc  
**Local file:** `references/boil-vn.pdf`

### Mathematical setting

The paper studies the magnetic Schrodinger equation in the plane:

```math
i\hbar \partial_t \varphi_\hbar^t = L_{\hbar,A}\varphi_\hbar^t,
```

with magnetic Laplacian

```math
L_{\hbar,A}=(-i\hbar\partial_1-A_1(q))^2+(-i\hbar\partial_2-A_2(q))^2.
```

The magnetic field is a smooth nonvanishing function `B(q)` and the vector potential satisfies

```math
\partial_2 A_1 - \partial_1 A_2 = B.
```

The initial states are coherent states: wave packets localized in phase space at scale approximately `sqrt(hbar)`.

### Main idea

For finite times, semiclassical propagation says a coherent state follows the classical Hamiltonian trajectory. The paper asks what happens on the longer time scale `t = O(1/hbar)` in a strong magnetic field.

The physical picture is:

1. a charged particle rapidly rotates because of the magnetic field;
2. its guiding center drifts slowly along level lines of the magnetic field;
3. a quantum coherent state cannot stay as a single purely classical particle because of uncertainty and energy quantization;
4. the initial packet decomposes into several coherent components;
5. each component follows the same guiding-center dynamics but with a different quantized speed.

The paper calls this a quantum ubiquity or splitting phenomenon.

### Why the theorem is hard to reproduce directly

The theorem is not just a numerical statement about one wave packet. It uses:

- semiclassical wavefront sets,
- adiabatic normal forms,
- a quantum adiabatic invariant,
- Landau-level-like spectral projectors,
- long-time control on intervals of length `O(1/hbar)`,
- remainder estimates of order `O(hbar^infinity)`.

A finite-grid simulation cannot prove these estimates. It can only explore related behavior.

## 2. On the stability of the Schwartz class under the magnetic Schrodinger flow

**Authors:** Gregory Boil, Nicolas Raymond, and San Vu Ngoc  
**Local file:** `references/1806.05389v1.pdf`  
**arXiv:** 1806.05389

### Main result

This paper proves that, under suitable assumptions on the magnetic field, the Schwartz class is stable under the magnetic Schrodinger flow. In simpler language:

> If the initial wavefunction is smooth and rapidly decaying, then magnetic Schrodinger evolution keeps it smooth and rapidly decaying.

The result also gives explicit control of Schwartz seminorms in terms of the semiclassical parameter.

### Why it matters here

Coherent states are built from Schwartz-class profiles. Long-time coherent-state propagation requires control not only of `L2` mass but also of localization, decay, and derivatives. This stability result supports the more refined coherent-state theorem.

## 3. Translation into code

| Paper concept | Code representation |
|---|---|
| Magnetic field `B(q)` | `magnetic_field_reference(x, y)` |
| Vector potential `A(q)` | `vector_potential_reference(x, y)` |
| Magnetic Laplacian `(-ih nabla - A)^2` | `build_magnetic_hamiltonian(...)` |
| Coherent state | `coherent_state(...)` |
| Schrodinger evolution | `evolve(...)` using sparse matrix exponentials |
| Probability density | `density(psi)` |
| Localization diagnostics | center of mass, variance, participation ratio |
| Guiding-center intuition | `examples/classical_guiding_center.py` |

## 4. Honest limitation statement

The implemented project is a finished computational research prototype, not a finished mathematical reproduction of the theorem. The current code demonstrates the core numerical mechanism and creates reproducible outputs, but the following paper-level components remain future research:

- construction of adiabatic Landau-level projectors;
- numerical decomposition into several coherent components;
- proof-grade comparison with the guiding-center Hamiltonian flow;
- rigorous long-time error control;
- Schwartz-seminorm diagnostics on a finite grid.
