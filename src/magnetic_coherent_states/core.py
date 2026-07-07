"""
Core numerical routines for magnetic coherent-state experiments.

The model operator is

    L_{h,A} = (-i h nabla - A)^2

on a two-dimensional rectangular grid. The discretization uses magnetic link
variables, also called Peierls phases, to keep gauge information attached to
nearest-neighbor lattice links.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import coo_matrix, csr_matrix, spmatrix
from scipy.sparse.linalg import expm_multiply

Array = NDArray[np.complex128]
RealArray = NDArray[np.float64]
VectorPotential = Callable[[RealArray, RealArray], tuple[RealArray, RealArray]]


@dataclass(frozen=True)
class Grid2D:
    """Uniform rectangular grid with a Dirichlet exterior approximation."""

    xmin: float = -6.0
    xmax: float = 6.0
    ymin: float = -6.0
    ymax: float = 6.0
    nx: int = 96
    ny: int = 96

    def __post_init__(self) -> None:
        if self.nx < 3 or self.ny < 3:
            raise ValueError("Grid dimensions must satisfy nx >= 3 and ny >= 3.")
        if not self.xmax > self.xmin:
            raise ValueError("Expected xmax > xmin.")
        if not self.ymax > self.ymin:
            raise ValueError("Expected ymax > ymin.")

    @property
    def x(self) -> RealArray:
        return np.linspace(self.xmin, self.xmax, self.nx, dtype=np.float64)

    @property
    def y(self) -> RealArray:
        return np.linspace(self.ymin, self.ymax, self.ny, dtype=np.float64)

    @property
    def dx(self) -> float:
        return float((self.xmax - self.xmin) / (self.nx - 1))

    @property
    def dy(self) -> float:
        return float((self.ymax - self.ymin) / (self.ny - 1))

    @property
    def cell_area(self) -> float:
        return self.dx * self.dy

    @property
    def mesh(self) -> tuple[RealArray, RealArray]:
        return np.meshgrid(self.x, self.y, indexing="ij")

    @property
    def shape(self) -> tuple[int, int]:
        return self.nx, self.ny

    @property
    def size(self) -> int:
        return self.nx * self.ny

    def flatten_index(self, i: int, j: int) -> int:
        if not (0 <= i < self.nx and 0 <= j < self.ny):
            raise IndexError("Grid index out of bounds.")
        return i * self.ny + j


def _require_positive_hbar(hbar: float) -> None:
    if not np.isfinite(hbar) or hbar <= 0:
        raise ValueError("hbar must be a positive finite number.")


def magnetic_field_reference(x: RealArray, y: RealArray) -> RealArray:
    """Reference field: B(x,y) = 2 - cos(x) + y^2."""

    return 2.0 - np.cos(x) + y**2


def vector_potential_reference(x: RealArray, y: RealArray) -> tuple[RealArray, RealArray]:
    """
    Vector potential compatible with ``magnetic_field_reference``.

    A_1(x,y) = (2 - cos(x)) y + y^3 / 3,  A_2(x,y) = 0.
    Therefore partial_y A_1 - partial_x A_2 = 2 - cos(x) + y^2.
    """

    a1 = (2.0 - np.cos(x)) * y + (y**3) / 3.0
    a2 = np.zeros_like(a1)
    return a1.astype(np.float64), a2.astype(np.float64)


def discrete_inner_product(phi: Array, psi: Array, grid: Grid2D) -> complex:
    """Discrete approximation of <phi, psi> = integral conj(phi) psi dx dy."""

    if phi.shape != grid.shape or psi.shape != grid.shape:
        raise ValueError("Wavefunction shape does not match grid shape.")
    return complex(np.sum(np.conjugate(phi) * psi) * grid.cell_area)


def l2_mass(psi: Array, grid: Grid2D) -> float:
    """Discrete L2 mass integral |psi|^2 dx dy."""

    return float(np.real(discrete_inner_product(psi, psi, grid)))


def normalize(psi: Array, grid: Grid2D) -> Array:
    """Normalize a wavefunction in the discrete L2 norm."""

    mass = l2_mass(psi, grid)
    if mass <= 0 or not np.isfinite(mass):
        raise ValueError("Cannot normalize a zero or invalid wavefunction.")
    return (psi / np.sqrt(mass)).astype(np.complex128)


def coherent_state(
    grid: Grid2D,
    hbar: float,
    q0: tuple[float, float] = (0.0, 0.0),
    p0: tuple[float, float] = (0.0, 0.0),
    width_factor: float = 1.0,
) -> Array:
    """
    Construct a normalized coherent-state-like Gaussian wave packet.

    The packet has envelope width proportional to sqrt(hbar). It is not a
    theorem-level coherent-state object from the paper; it is the numerical
    analogue used for finite-grid experiments.
    """

    _require_positive_hbar(hbar)
    if width_factor <= 0:
        raise ValueError("width_factor must be positive.")

    x, y = grid.mesh
    x0, y0 = q0
    px, py = p0
    sigma2 = hbar * width_factor**2
    envelope = np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2.0 * sigma2))
    phase = np.exp(1j * (px * (x - x0) + py * (y - y0)) / hbar)
    return normalize((envelope * phase).astype(np.complex128), grid)


def build_magnetic_hamiltonian(
    grid: Grid2D,
    hbar: float,
    vector_potential: VectorPotential = vector_potential_reference,
) -> csr_matrix:
    """
    Build a sparse gauge-covariant finite-difference Hamiltonian.

    The boundary is handled as a Dirichlet exterior approximation: missing
    neighbors outside the grid are treated as zero values. The diagonal term is
    therefore the same on boundary and interior nodes.
    """

    _require_positive_hbar(hbar)
    x_values = grid.x
    y_values = grid.y
    dx, dy = grid.dx, grid.dy

    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []

    diag_value = 2.0 * hbar**2 / dx**2 + 2.0 * hbar**2 / dy**2

    def add(row: int, col: int, value: complex) -> None:
        rows.append(row)
        cols.append(col)
        data.append(value)

    for i, x in enumerate(x_values):
        for j, y in enumerate(y_values):
            k = grid.flatten_index(i, j)
            add(k, k, diag_value)

            if i + 1 < grid.nx:
                x_mid = np.array([[x + dx / 2.0]], dtype=np.float64)
                y_mid = np.array([[y]], dtype=np.float64)
                ax_mid, _ = vector_potential(x_mid, y_mid)
                ux = np.exp(-1j * float(ax_mid[0, 0]) * dx / hbar)
                add(k, grid.flatten_index(i + 1, j), -(hbar**2 / dx**2) * ux)

            if i - 1 >= 0:
                x_mid = np.array([[x - dx / 2.0]], dtype=np.float64)
                y_mid = np.array([[y]], dtype=np.float64)
                ax_mid, _ = vector_potential(x_mid, y_mid)
                ux = np.exp(+1j * float(ax_mid[0, 0]) * dx / hbar)
                add(k, grid.flatten_index(i - 1, j), -(hbar**2 / dx**2) * ux)

            if j + 1 < grid.ny:
                x_mid = np.array([[x]], dtype=np.float64)
                y_mid = np.array([[y + dy / 2.0]], dtype=np.float64)
                _, ay_mid = vector_potential(x_mid, y_mid)
                uy = np.exp(-1j * float(ay_mid[0, 0]) * dy / hbar)
                add(k, grid.flatten_index(i, j + 1), -(hbar**2 / dy**2) * uy)

            if j - 1 >= 0:
                x_mid = np.array([[x]], dtype=np.float64)
                y_mid = np.array([[y - dy / 2.0]], dtype=np.float64)
                _, ay_mid = vector_potential(x_mid, y_mid)
                uy = np.exp(+1j * float(ay_mid[0, 0]) * dy / hbar)
                add(k, grid.flatten_index(i, j - 1), -(hbar**2 / dy**2) * uy)

    return coo_matrix((data, (rows, cols)), shape=(grid.size, grid.size)).tocsr()


# Backward-compatible name from the first prototype.
build_magnetic_laplacian = build_magnetic_hamiltonian


def is_hermitian(matrix: spmatrix, tolerance: float = 1e-10) -> bool:
    """Return True when a sparse matrix is Hermitian within a max-entry tolerance."""

    diff = (matrix - matrix.getH()).tocoo()
    if diff.nnz == 0:
        return True
    return bool(np.max(np.abs(diff.data)) <= tolerance)


def evolve(
    hamiltonian: csr_matrix,
    psi0: Array,
    grid: Grid2D,
    hbar: float,
    times: Iterable[float],
    *,
    renormalize: bool = True,
) -> list[Array]:
    """Evolve ``psi0`` at requested times using sparse matrix exponentials."""

    _require_positive_hbar(hbar)
    if psi0.shape != grid.shape:
        raise ValueError("psi0 shape does not match grid shape.")
    if hamiltonian.shape != (grid.size, grid.size):
        raise ValueError("Hamiltonian shape does not match grid size.")

    time_values = np.asarray(list(times), dtype=np.float64)
    if time_values.ndim != 1:
        raise ValueError("times must be a one-dimensional iterable.")
    if np.any(~np.isfinite(time_values)):
        raise ValueError("times must be finite.")

    psi0_vec = psi0.reshape(-1).astype(np.complex128)
    generator = (-1j / hbar) * hamiltonian

    frames: list[Array] = []
    for t in time_values:
        psi_t = expm_multiply(generator * float(t), psi0_vec).reshape(grid.shape)
        if renormalize:
            psi_t = normalize(psi_t, grid)
        frames.append(psi_t.astype(np.complex128))
    return frames


def density(psi: Array) -> RealArray:
    """Probability density |psi|^2."""

    return (np.abs(psi) ** 2).astype(np.float64)


def center_of_mass(psi: Array, grid: Grid2D) -> tuple[float, float]:
    """Discrete probability center of mass."""

    rho = density(psi)
    mass = float(np.sum(rho) * grid.cell_area)
    if mass <= 0:
        raise ValueError("Cannot compute center of mass for zero mass.")
    x, y = grid.mesh
    cx = float(np.sum(x * rho) * grid.cell_area / mass)
    cy = float(np.sum(y * rho) * grid.cell_area / mass)
    return cx, cy


def variance(psi: Array, grid: Grid2D) -> tuple[float, float]:
    """Coordinate variances of the probability density."""

    cx, cy = center_of_mass(psi, grid)
    x, y = grid.mesh
    rho = density(psi)
    mass = float(np.sum(rho) * grid.cell_area)
    vx = float(np.sum((x - cx) ** 2 * rho) * grid.cell_area / mass)
    vy = float(np.sum((y - cy) ** 2 * rho) * grid.cell_area / mass)
    return vx, vy


def participation_ratio(psi: Array, grid: Grid2D) -> float:
    """Spreading diagnostic: PR = 1 / integral |psi|^4 dx dy."""

    rho = density(psi)
    denom = float(np.sum(rho**2) * grid.cell_area)
    if denom <= 0:
        raise ValueError("Cannot compute participation ratio for zero density.")
    return 1.0 / denom


def save_npz(
    path: str | Path,
    grid: Grid2D,
    times: Sequence[float] | RealArray,
    frames: Sequence[Array],
    metadata: Mapping[str, object] | None = None,
) -> None:
    """Save simulation data in a compressed NumPy archive."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        x=grid.x,
        y=grid.y,
        times=np.asarray(times, dtype=np.float64),
        psi=np.stack(list(frames), axis=0),
        metadata=np.array([dict(metadata or {})], dtype=object),
    )
