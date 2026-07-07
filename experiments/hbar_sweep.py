"""Run a small hbar sweep and compare final diagnostics."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import numpy as np

from magnetic_coherent_states import (
    Grid2D,
    build_magnetic_hamiltonian,
    center_of_mass,
    coherent_state,
    evolve,
    participation_ratio,
    vector_potential_reference,
    variance,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="Use fewer hbar values and a smaller grid.")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    hbars = [0.24, 0.18, 0.14] if args.quick else [0.28, 0.22, 0.18, 0.14, 0.11]
    n = 42 if args.quick else 58
    frames = 4 if args.quick else 6

    rows: list[dict[str, float]] = []
    for hbar in hbars:
        grid = Grid2D(xmin=-5, xmax=5, ymin=-5, ymax=5, nx=n, ny=n)
        hamiltonian = build_magnetic_hamiltonian(grid, hbar, vector_potential_reference)
        psi0 = coherent_state(grid, hbar, q0=(-1.8, 0.35), width_factor=1.1)
        times = np.linspace(0.0, 0.7 / hbar, frames)
        final = evolve(hamiltonian, psi0, grid, hbar, times)[-1]
        cx, cy = center_of_mass(final, grid)
        vx, vy = variance(final, grid)
        pr = participation_ratio(final, grid)
        rows.append(
            {
                "hbar": float(hbar),
                "final_time": float(times[-1]),
                "center_x": cx,
                "center_y": cy,
                "variance_x": vx,
                "variance_y": vy,
                "participation_ratio": pr,
            }
        )

    with (args.out_dir / "hbar_sweep.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {args.out_dir / 'hbar_sweep.csv'}")


if __name__ == "__main__":
    main()
