"""Command-line interface for the magnetic coherent-state demo."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from tqdm import tqdm

from .core import (
    Grid2D,
    build_magnetic_hamiltonian,
    center_of_mass,
    coherent_state,
    evolve,
    is_hermitian,
    participation_ratio,
    save_npz,
    vector_potential_reference,
    variance,
)
from .plotting import (
    plot_center_trajectory,
    plot_density,
    plot_magnetic_field,
    write_research_summary,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a reproducible magnetic coherent-state simulation.")
    parser.add_argument("--quick", action="store_true", help="Use a small grid for a fast smoke run.")
    parser.add_argument("--nx", type=int, default=72, help="Number of grid points in x.")
    parser.add_argument("--ny", type=int, default=72, help="Number of grid points in y.")
    parser.add_argument("--hbar", type=float, default=0.16, help="Semiclassical parameter hbar.")
    parser.add_argument("--frames", type=int, default=9, help="Number of saved time frames.")
    parser.add_argument("--tmax-scale", type=float, default=1.4, help="Final time is tmax_scale / hbar.")
    parser.add_argument("--out-dir", type=Path, default=Path("outputs"), help="Output directory.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.quick:
        args.nx = min(args.nx, 48)
        args.ny = min(args.ny, 48)
        args.frames = min(args.frames, 5)
        args.tmax_scale = min(args.tmax_scale, 0.7)

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    grid = Grid2D(xmin=-5, xmax=5, ymin=-5, ymax=5, nx=args.nx, ny=args.ny)
    hbar = args.hbar
    times = np.linspace(0.0, args.tmax_scale / hbar, args.frames)

    print(f"Grid: {grid.nx} x {grid.ny}")
    print(f"hbar: {hbar}")
    print("Building Hamiltonian...")
    hamiltonian = build_magnetic_hamiltonian(grid, hbar, vector_potential_reference)
    print(f"Hamiltonian shape: {hamiltonian.shape}, nnz={hamiltonian.nnz}")
    print(f"Hermitian check: {is_hermitian(hamiltonian)}")

    print("Creating initial coherent state...")
    psi0 = coherent_state(
        grid=grid,
        hbar=hbar,
        q0=(-1.8, 0.35),
        p0=(0.0, 0.0),
        width_factor=1.1,
    )

    print("Evolving...")
    frames = []
    for psi_t in tqdm(evolve(hamiltonian, psi0, grid, hbar, times), total=len(times)):
        frames.append(psi_t)

    print("Saving figures and data...")
    plot_magnetic_field(grid, out_dir / "magnetic_field.png")

    centers: list[tuple[float, float]] = []
    rows: list[dict[str, float]] = []
    for idx, (t, psi_t) in enumerate(zip(times, frames)):
        cx, cy = center_of_mass(psi_t, grid)
        vx, vy = variance(psi_t, grid)
        pr = participation_ratio(psi_t, grid)
        centers.append((cx, cy))
        rows.append(
            {
                "time": float(t),
                "center_x": cx,
                "center_y": cy,
                "variance_x": vx,
                "variance_y": vy,
                "participation_ratio": pr,
            }
        )
        plot_density(
            grid,
            psi_t,
            title=fr"Density at $t={t:.2f}$, $h={hbar}$",
            out_path=out_dir / f"density_{idx:03d}.png",
        )

    plot_center_trajectory(centers, out_dir / "center_trajectory.png")

    with (out_dir / "stats.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    save_npz(
        out_dir / "simulation_data.npz",
        grid=grid,
        times=times,
        frames=frames,
        metadata={
            "hbar": hbar,
            "field": "B(x,y)=2-cos(x)+y^2",
            "operator": "gauge-covariant finite-difference magnetic Hamiltonian",
            "grid": f"{grid.nx}x{grid.ny}",
        },
    )
    write_research_summary(out_dir / "research_summary.md", grid=grid, hbar=hbar, times=times, frames=frames)

    print(f"Done. See {out_dir}.")


if __name__ == "__main__":
    main()
