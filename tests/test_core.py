from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from magnetic_coherent_states import (
    Grid2D,
    build_magnetic_hamiltonian,
    coherent_state,
    evolve,
    is_hermitian,
    l2_mass,
    magnetic_field_reference,
    vector_potential_reference,
)


class CoreNumericsTest(unittest.TestCase):
    def test_coherent_state_is_normalized(self) -> None:
        grid = Grid2D(nx=32, ny=32)
        psi = coherent_state(grid, hbar=0.2, q0=(0.1, -0.2), p0=(0.3, 0.0))
        self.assertAlmostEqual(l2_mass(psi, grid), 1.0, places=10)

    def test_reference_vector_potential_has_expected_curl(self) -> None:
        grid = Grid2D(xmin=-2, xmax=2, ymin=-2, ymax=2, nx=81, ny=81)
        x, y = grid.mesh
        a1, a2 = vector_potential_reference(x, y)
        partial_y_a1 = np.gradient(a1, grid.dy, axis=1)
        partial_x_a2 = np.gradient(a2, grid.dx, axis=0)
        curl = partial_y_a1 - partial_x_a2
        expected = magnetic_field_reference(x, y)
        err = np.max(np.abs(curl[2:-2, 2:-2] - expected[2:-2, 2:-2]))
        self.assertLess(err, 3e-3)

    def test_hamiltonian_is_hermitian(self) -> None:
        grid = Grid2D(nx=16, ny=15)
        hamiltonian = build_magnetic_hamiltonian(grid, hbar=0.25)
        self.assertTrue(is_hermitian(hamiltonian, tolerance=1e-12))

    def test_short_evolution_preserves_mass(self) -> None:
        grid = Grid2D(nx=18, ny=18)
        hbar = 0.25
        hamiltonian = build_magnetic_hamiltonian(grid, hbar=hbar)
        psi0 = coherent_state(grid, hbar=hbar, q0=(-0.4, 0.3))
        frames = evolve(hamiltonian, psi0, grid, hbar, [0.0, 0.1, 0.2])
        for frame in frames:
            self.assertAlmostEqual(l2_mass(frame, grid), 1.0, places=9)


if __name__ == "__main__":
    unittest.main()
