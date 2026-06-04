"""Generate Toeplitz alternating-minimization figures for the manuscript.

The script runs the Toeplitz alternating least-squares method, optionally
refines the result by nonlinear least squares, and writes only PNG figures.
It is intentionally self-contained so that the manuscript figures can be
recreated without the older exploratory scripts used during drafting.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from scipy.optimize import least_squares


Matrix = np.ndarray
Vector = np.ndarray


RESIDUAL_FIGURE_NAMES = {
    "toeplitz_3x3_ye_lim": "toeplitz_3x3_residual.png",
    "toeplitz_5x5_ye_lim_example2_depth15": "toeplitz_5x5_r15_residual.png",
}

FACTOR_GRID_FIGURE_NAMES = {
    "toeplitz_5x5_ye_lim_example2_depth15": "toeplitz_5x5_r15_factors.png",
}


@dataclass(frozen=True)
class ToeplitzExperiment:
    key: str
    target: Matrix
    num_factors: int
    als_sweeps: int
    tol: float
    max_nfev: int
    initial_factors: tuple[Matrix, ...] | None = None
    make_factor_grid: bool = False


@dataclass(frozen=True)
class ToeplitzResult:
    factors: tuple[Matrix, ...]
    residual_history: tuple[float, ...]
    final_residual: float
    nfev: int
    message: str


def toeplitz(first_column: list[float], first_row: list[float]) -> Matrix:
    """Construct a Toeplitz matrix from its first column and first row."""
    n = len(first_column)
    if len(first_row) != n:
        raise ValueError("first column and first row must have the same length")
    if first_column[0] != first_row[0]:
        raise ValueError("first column and first row must agree at the first entry")

    matrix = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            matrix[i, j] = first_row[j - i] if j >= i else first_column[i - j]
    return matrix


def toeplitz_linear_map(n: int) -> Matrix:
    """Return the matrix S satisfying vec(T)=Sb for Toeplitz parameters b."""
    columns: list[Vector] = []
    for diagonal in range(-(n - 1), n):
        basis = np.zeros((n, n), dtype=np.float64)
        for row in range(n):
            col = row + diagonal
            if 0 <= col < n:
                basis[row, col] = 1.0
        columns.append(basis.reshape(-1, order="F"))
    return np.column_stack(columns)


def toeplitz_from_params(params: Vector, linear_map: Matrix, n: int) -> Matrix:
    return (linear_map @ params).reshape((n, n), order="F")


def project_to_toeplitz(matrix: Matrix, linear_map: Matrix) -> Vector:
    rhs = matrix.reshape(-1, order="F")
    params, _, _, _ = np.linalg.lstsq(linear_map, rhs, rcond=None)
    return np.asarray(params, dtype=np.float64)


def product_of_factors(factors: list[Matrix] | tuple[Matrix, ...]) -> Matrix:
    product = np.eye(factors[0].shape[0], dtype=np.float64)
    for factor in factors:
        product = product @ factor
    return product


def residual_norm(target: Matrix, factors: list[Matrix] | tuple[Matrix, ...]) -> float:
    return float(np.linalg.norm(product_of_factors(factors) - target, ord="fro"))


def make_initial_parameters(
    target: Matrix,
    num_factors: int,
    linear_map: Matrix,
    initial_factors: tuple[Matrix, ...] | None,
) -> list[Vector]:
    if initial_factors is not None:
        if len(initial_factors) != num_factors:
            raise ValueError("initial_factors must match num_factors")
        return [project_to_toeplitz(factor, linear_map) for factor in initial_factors]

    identity_params = project_to_toeplitz(np.eye(target.shape[0], dtype=np.float64), linear_map)
    params = [identity_params.copy() for _ in range(num_factors)]
    params[0] = project_to_toeplitz(target, linear_map)
    return params


def block_best_response(target: Matrix, params: list[Vector], factor_index: int, linear_map: Matrix) -> Vector:
    n = target.shape[0]
    factors = [toeplitz_from_params(block, linear_map, n) for block in params]

    left = np.eye(n, dtype=np.float64)
    for factor in factors[:factor_index]:
        left = left @ factor

    right = np.eye(n, dtype=np.float64)
    for factor in factors[factor_index + 1 :]:
        right = right @ factor

    system = np.kron(right.T, left) @ linear_map
    rhs = target.reshape(-1, order="F")
    solution, _, _, _ = np.linalg.lstsq(system, rhs, rcond=None)
    return np.asarray(solution, dtype=np.float64)


def run_als_sweeps(target: Matrix, params: list[Vector], linear_map: Matrix, sweeps: int) -> tuple[list[Vector], list[float]]:
    n = target.shape[0]

    def current_factors() -> tuple[Matrix, ...]:
        return tuple(toeplitz_from_params(block, linear_map, n) for block in params)

    history = [residual_norm(target, current_factors())]
    for _ in range(sweeps):
        for factor_index in range(len(params)):
            params[factor_index] = block_best_response(target, params, factor_index, linear_map)
        history.append(residual_norm(target, current_factors()))
    return params, history


def refine_by_nonlinear_least_squares(
    target: Matrix,
    params: list[Vector],
    linear_map: Matrix,
    *,
    tol: float,
    max_nfev: int,
) -> tuple[list[Vector], list[float], int, str]:
    n = target.shape[0]
    width = linear_map.shape[1]
    x0 = np.concatenate(params)
    callback_history: list[float] = []

    def unpack(vector: Vector) -> list[Vector]:
        return [np.asarray(vector[i * width : (i + 1) * width], dtype=np.float64) for i in range(len(params))]

    def factors_from_vector(vector: Vector) -> tuple[Matrix, ...]:
        return tuple(toeplitz_from_params(block, linear_map, n) for block in unpack(vector))

    def objective(vector: Vector) -> Vector:
        factors = factors_from_vector(vector)
        return (product_of_factors(factors) - target).reshape(-1, order="F")

    def callback(vector: Vector, *args: object, **kwargs: object) -> None:
        callback_history.append(float(np.linalg.norm(objective(vector))))

    kwargs = {
        "method": "trf",
        "ftol": tol,
        "xtol": tol,
        "gtol": tol,
        "max_nfev": max_nfev,
    }

    try:
        result = least_squares(objective, x0, callback=callback, **kwargs)
    except TypeError:
        result = least_squares(objective, x0, **kwargs)

    refined_params = unpack(np.asarray(result.x, dtype=np.float64))
    final_residual = float(np.linalg.norm(objective(result.x)))
    if not callback_history or abs(callback_history[-1] - final_residual) > 1e-15:
        callback_history.append(final_residual)
    return refined_params, callback_history, int(result.nfev), str(result.message)


def solve_experiment(experiment: ToeplitzExperiment) -> ToeplitzResult:
    n = experiment.target.shape[0]
    linear_map = toeplitz_linear_map(n)
    params = make_initial_parameters(
        experiment.target,
        experiment.num_factors,
        linear_map,
        experiment.initial_factors,
    )

    params, als_history = run_als_sweeps(experiment.target, params, linear_map, experiment.als_sweeps)
    params, refinement_history, nfev, message = refine_by_nonlinear_least_squares(
        experiment.target,
        params,
        linear_map,
        tol=experiment.tol,
        max_nfev=experiment.max_nfev,
    )

    factors = tuple(toeplitz_from_params(block, linear_map, n) for block in params)
    final_residual = residual_norm(experiment.target, factors)
    history = list(als_history) + list(refinement_history)
    if not history or abs(history[-1] - final_residual) > 1e-15:
        history.append(final_residual)

    return ToeplitzResult(
        factors=factors,
        residual_history=tuple(history),
        final_residual=final_residual,
        nfev=nfev,
        message=message,
    )


def write_residual_plot(result: ToeplitzResult, output_path: Path) -> None:
    iterations = np.arange(len(result.residual_history))
    residuals = np.maximum(np.asarray(result.residual_history, dtype=np.float64), 1e-16)

    plt.figure(figsize=(5.4, 3.4))
    plt.plot(iterations, residuals, marker="o", markersize=3.5, linewidth=1.6, color="#1f77b4")
    plt.xlabel("iteration")
    plt.ylabel("residual")
    plt.yscale("log")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close()


def write_factor_grid(factors: tuple[Matrix, ...], output_path: Path, *, columns: int = 5) -> None:
    rows = int(np.ceil(len(factors) / columns))
    max_abs = max(float(np.max(np.abs(factor))) for factor in factors)
    max_abs = max(max_abs, 1.0)
    norm = colors.TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)

    fig, axes = plt.subplots(rows, columns, figsize=(2.45 * columns, 2.45 * rows + 0.6))
    axes_array = np.asarray(axes, dtype=object).reshape(rows, columns)
    image = None

    for idx, ax in enumerate(axes_array.flat):
        if idx >= len(factors):
            ax.axis("off")
            continue

        factor = factors[idx]
        image = ax.imshow(factor, cmap="coolwarm", norm=norm)
        ax.set_title(f"$M_{{{idx + 1}}}$", fontsize=10, pad=4)
        ax.set_xticks([])
        ax.set_yticks([])
        n = factor.shape[0]
        ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.8)
        ax.tick_params(which="minor", bottom=False, left=False)
        for spine in ax.spines.values():
            spine.set_visible(False)

    if image is None:
        raise RuntimeError("no factors were available for plotting")

    fig.subplots_adjust(left=0.04, right=0.92, top=0.92, bottom=0.05, wspace=0.25, hspace=0.32)
    colorbar = fig.colorbar(image, ax=axes_array.ravel().tolist(), fraction=0.025, pad=0.02)
    colorbar.set_label("entry value", rotation=90)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_experiments() -> list[ToeplitzExperiment]:
    target_3x3 = np.asarray(
        [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ],
        dtype=np.float64,
    )
    initial_3x3 = (
        toeplitz([5, 6, 7], [5, 4, 3]),
        np.eye(3, dtype=np.float64),
    )

    target_5x5 = np.asarray(
        [
            [2, 5, 2, 5, 3],
            [4, 5, 5, 2, 2],
            [2, 3, 2, 1, 5],
            [3, 1, 5, 2, 3],
            [4, 1, 2, 4, 3],
        ],
        dtype=np.float64,
    )

    return [
        ToeplitzExperiment(
            key="toeplitz_3x3_ye_lim",
            target=target_3x3,
            num_factors=2,
            als_sweeps=16,
            tol=1e-14,
            max_nfev=6000,
            initial_factors=initial_3x3,
        ),
        ToeplitzExperiment(
            key="toeplitz_5x5_ye_lim_example2_depth15",
            target=target_5x5,
            num_factors=15,
            als_sweeps=12,
            tol=1e-14,
            max_nfev=12000,
            make_factor_grid=True,
        ),
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "output",
        help="Directory where PNG figures are written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for experiment in build_experiments():
        result = solve_experiment(experiment)
        residual_path = args.output_dir / RESIDUAL_FIGURE_NAMES[experiment.key]
        write_residual_plot(result, residual_path)
        print(f"{experiment.key}: residual={result.final_residual:.3e}, nfev={result.nfev}")
        print(f"Wrote {residual_path}")

        if experiment.make_factor_grid:
            grid_path = args.output_dir / FACTOR_GRID_FIGURE_NAMES[experiment.key]
            write_factor_grid(result.factors, grid_path)
            print(f"Wrote {grid_path}")


if __name__ == "__main__":
    main()
