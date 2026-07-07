"""Plotting and report helpers for magnetic coherent-state experiments."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .core import Array, Grid2D, center_of_mass, density, magnetic_field_reference, participation_ratio, variance


def plot_density(grid: Grid2D, psi: Array, title: str, out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rho = density(psi)
    plt.figure(figsize=(7, 6))
    plt.imshow(
        rho.T,
        origin="lower",
        extent=[grid.xmin, grid.xmax, grid.ymin, grid.ymax],
        aspect="equal",
    )
    plt.colorbar(label=r"$|\psi(x,y,t)|^2$")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close()


def plot_magnetic_field(grid: Grid2D, out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    x, y = grid.mesh
    b = magnetic_field_reference(x, y)
    plt.figure(figsize=(7, 6))
    plt.imshow(
        b.T,
        origin="lower",
        extent=[grid.xmin, grid.xmax, grid.ymin, grid.ymax],
        aspect="equal",
    )
    plt.colorbar(label="B(x,y)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(r"Magnetic field $B(x,y)=2-\cos(x)+y^2$")
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close()


def plot_center_trajectory(centers: Sequence[tuple[float, float]], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.asarray(centers, dtype=np.float64)
    plt.figure(figsize=(6, 6))
    plt.plot(arr[:, 0], arr[:, 1], marker="o")
    plt.xlabel("center x")
    plt.ylabel("center y")
    plt.title("Wave-packet center-of-mass trajectory")
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(out_path, dpi=180)
    plt.close()


def write_research_summary(
    out_path: str | Path,
    *,
    grid: Grid2D,
    hbar: float,
    times: Sequence[float],
    frames: Sequence[Array],
) -> None:
    """Write a compact Markdown summary of one simulation run."""

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for t, psi in zip(times, frames):
        cx, cy = center_of_mass(psi, grid)
        vx, vy = variance(psi, grid)
        pr = participation_ratio(psi, grid)
        rows.append((float(t), cx, cy, vx, vy, pr))

    initial = rows[0]
    final = rows[-1]
    text = f"""# Simulation Summary

## Configuration

- Grid: `{grid.nx} x {grid.ny}`
- Domain: `[{grid.xmin}, {grid.xmax}] x [{grid.ymin}, {grid.ymax}]`
- hbar: `{hbar}`
- Frames: `{len(rows)}`
- Time interval: `{times[0]:.6g}` to `{times[-1]:.6g}`

## Main diagnostics

| quantity | initial | final |
|---|---:|---:|
| center x | {initial[1]:.6f} | {final[1]:.6f} |
| center y | {initial[2]:.6f} | {final[2]:.6f} |
| variance x | {initial[3]:.6f} | {final[3]:.6f} |
| variance y | {initial[4]:.6f} | {final[4]:.6f} |
| participation ratio | {initial[5]:.6f} | {final[5]:.6f} |

## Interpretation

The packet begins as a localized Gaussian-like state. Under the magnetic Hamiltonian, the probability density moves and spreads. The center-of-mass trajectory is a finite-grid diagnostic, not a proof of the asymptotic guiding-center theorem. The participation ratio gives a compact scalar proxy for spreading: larger values indicate that the probability density occupies a larger effective area.
"""
    out_path.write_text(text, encoding="utf-8")
