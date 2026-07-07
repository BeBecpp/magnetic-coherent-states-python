"""Toy classical guiding-center model for the reference magnetic field."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from magnetic_coherent_states import Grid2D, magnetic_field_reference


def grad_b(x: float, y: float) -> tuple[float, float]:
    return float(np.sin(x)), float(2.0 * y)


def guiding_center_rhs(_t: float, z: np.ndarray, speed: float = 1.0) -> list[float]:
    """
    Toy drift along magnetic-field level lines.

    This is an intuition plot, not the full Hamiltonian flow from the theorem.
    """

    x, y = float(z[0]), float(z[1])
    bx, by = grad_b(x, y)
    b_value = float(2.0 - np.cos(x) + y**2)
    denom = b_value**2 + 1e-12
    return [speed * (-by) / denom, speed * bx / denom]


def main() -> None:
    out_dir = ROOT / "outputs"
    out_dir.mkdir(exist_ok=True)
    grid = Grid2D(xmin=-5, xmax=5, ymin=-5, ymax=5, nx=220, ny=220)
    x, y = grid.mesh
    b = magnetic_field_reference(x, y)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.contour(x, y, b, levels=25)

    q0 = np.array([-1.8, 0.35], dtype=np.float64)
    for speed in [1.0, 3.0, 5.0]:
        sol = solve_ivp(
            lambda t, z: guiding_center_rhs(t, z, speed=speed),
            t_span=(0.0, 12.0),
            y0=q0,
            max_step=0.05,
        )
        ax.plot(sol.y[0], sol.y[1], label=f"speed={speed}")

    ax.scatter([q0[0]], [q0[1]], marker="x", label="initial center")
    ax.set_title("Toy guiding-center drift along B-level lines")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.axis("equal")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "classical_guiding_center.png", dpi=180)
    plt.close(fig)
    print(f"Saved {out_dir / 'classical_guiding_center.png'}")


if __name__ == "__main__":
    main()
