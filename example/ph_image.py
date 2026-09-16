"""
Compute persistent homology of a 3D grayscale volume via a cubical complex
with V-filtration (filtration values given on vertices / pixels, extended by
lower-star to edges, faces, and cubes).

Requires: numpy, gudhi
  pip install numpy gudhi
"""

from __future__ import annotations

from pathlib import Path

import gudhi as gd
import numpy as np


def load_grayscale_volume(path: Path, shape: tuple[int, int, int] = (10, 10, 10)) -> np.ndarray:
    """Load a volume whose pixels are stored as pixel(x,y,z)=data[x*ny*nz+y*nz+z]."""
    tokens: list[float] = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.split("#", 1)[0]
            tokens.extend(float(v) for v in line.split() if v)

    nx, ny, nz = shape
    expected = nx * ny * nz
    if len(tokens) != expected:
        raise ValueError(f"Expected {expected} values, got {len(tokens)} in {path}")

    # Flat order matches C++ ImageOctreeSF indexing; reshape to (x, y, z).
    volume = np.asarray(tokens, dtype=np.float64).reshape(shape)
    return volume


def compute_v_filtration_persistence(
    volume: np.ndarray,
    *,
    homology_coeff_field: int = 2,
    min_persistence: float = 0.0,
) -> list[tuple[int, tuple[float, float]]]:
    """
    Build a Gudhi CubicalComplex from vertex filtration values (V-construction)
    and compute persistence.

    filt(cube) = max{ filt(v) : v vertex of cube }  (lower-star / sublevel sets)
    """
    cc = gd.CubicalComplex(vertices=volume)
    return cc.persistence(
        homology_coeff_field=homology_coeff_field,
        min_persistence=min_persistence,
    )


def format_persistence(diagram: list[tuple[int, tuple[float, float]]]) -> str:
    by_dim: dict[int, list[tuple[float, float]]] = {}
    for dim, (birth, death) in diagram:
        by_dim.setdefault(dim, []).append((birth, death))

    lines = [f"Total intervals: {len(diagram)}"]
    for dim in sorted(by_dim):
        pairs = by_dim[dim]
        lines.append(f"\nH_{dim}: {len(pairs)} intervals")
        for birth, death in pairs:
            death_str = "inf" if death == float("inf") else f"{death:g}"
            lines.append(f"  [{birth:g}, {death_str})")
    return "\n".join(lines) + "\n"


def print_persistence(diagram: list[tuple[int, tuple[float, float]]]) -> None:
    print(format_persistence(diagram), end="")


def save_persistence(diagram: list[tuple[int, tuple[float, float]]], path: Path) -> None:
    path.write_text(format_persistence(diagram), encoding="utf-8")


def main() -> None:
    image_path = Path(__file__).resolve().parent / "image_10x10x10_spherical_shell.txt"
    volume = load_grayscale_volume(image_path, shape=(10, 10, 10))

    print(f"Loaded volume {volume.shape} from {image_path.name}")
    print(f"Value range: [{volume.min():g}, {volume.max():g}]")
    print("Cubical complex: V-filtration (vertices = pixel values)\n")

    diagram = compute_v_filtration_persistence(volume, homology_coeff_field=2)
    print_persistence(diagram)

    out_path = image_path.with_name(f"{image_path.stem}_ph.txt")
    save_persistence(diagram, out_path)
    print(f"\nSaved persistence to {out_path}")

    try:
        import matplotlib.pyplot as plt

        gd.plot_persistence_diagram(diagram)
        plt.title("Persistence diagram (V-filtration)")
        plt.tight_layout()
        plt.show()
    except ImportError:
        pass


if __name__ == "__main__":
    main()
