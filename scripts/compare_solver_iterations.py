#!/usr/bin/env python3

import argparse
import csv
import subprocess
from pathlib import Path


SOLVERS = ("simple", "projection")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run SIMPLE/projection cavity cases and plot iteration counts."
    )
    parser.add_argument(
        "--solver-path",
        default="build/cfd_solver",
        help="Path to the compiled cfd_solver executable.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/iteration_comparison",
        help="Directory for run artifacts, CSV, and SVG plot.",
    )
    parser.add_argument(
        "--grids",
        type=int,
        nargs="+",
        default=[8, 16, 32],
        help="Square grid sizes to compare.",
    )
    parser.add_argument("--re", type=float, default=100.0, help="Reynolds number.")
    parser.add_argument("--min-iters", type=int, default=50, help="Minimum nonlinear iterations.")
    parser.add_argument("--simple-max-iters", type=int, default=2500, help="SIMPLE iteration limit.")
    parser.add_argument(
        "--projection-max-iters",
        type=int,
        default=1000,
        help="Projection iteration limit.",
    )
    parser.add_argument(
        "--projection-dt",
        type=float,
        default=0.1,
        help="Projection pseudo-time step.",
    )
    return parser.parse_args()


def read_key_values(path):
    values = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key] = value
    return values


def read_last_residual(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No residual rows in {path}")
    return rows[-1]


def run_case(solver_path, output_root, solver, grid, re_value, min_iters, simple_max, projection_max, projection_dt):
    run_dir = output_root / "runs" / f"{solver}_nx{grid}"
    run_dir.mkdir(parents=True, exist_ok=True)

    max_iters = simple_max if solver == "simple" else projection_max
    cmd = [
        str(solver_path),
        "--case",
        "cavity",
        "--solver",
        solver,
        "--nx",
        str(grid),
        "--ny",
        str(grid),
        "--re",
        f"{re_value:g}",
        "--max-iters",
        str(max_iters),
        "--min-iters",
        str(min_iters),
        "--output-dir",
        str(run_dir),
    ]

    if solver == "simple":
        cmd.extend(["--alpha-u", "0.5", "--alpha-v", "0.5", "--alpha-p", "0.3"])
    else:
        cmd.extend(["--projection-dt", f"{projection_dt:g}"])

    completed = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{solver} nx={grid} failed with code {completed.returncode}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )

    summary = read_key_values(run_dir / "summary.txt")
    residual = read_last_residual(run_dir / "residuals.csv")
    return {
        "solver": solver,
        "nx": grid,
        "ny": grid,
        "re": re_value,
        "iterations": int(summary["iterations"]),
        "converged": summary["converged"],
        "continuity_residual": float(residual["continuity"]),
        "u_momentum_residual": float(residual["u_momentum"]),
        "v_momentum_residual": float(residual["v_momentum"]),
        "max_velocity_correction": float(residual["max_velocity_correction"]),
        "output_dir": str(run_dir),
    }


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "solver",
        "nx",
        "ny",
        "re",
        "iterations",
        "converged",
        "continuity_residual",
        "u_momentum_residual",
        "v_momentum_residual",
        "max_velocity_correction",
        "output_dir",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def svg_escape(text):
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def draw_svg(rows, path):
    by_grid = {}
    for row in rows:
        by_grid.setdefault(row["nx"], {})[row["solver"]] = row

    grids = sorted(by_grid)
    max_iterations = max(row["iterations"] for row in rows)
    y_max = int((max_iterations * 1.12 + 99) // 100 * 100)
    if y_max <= 0:
        y_max = 100

    width = max(780, 170 * len(grids) + 180)
    height = 520
    left = 82
    top = 78
    plot_width = width - 140
    plot_height = 315
    baseline = top + plot_height
    group_width = plot_width / len(grids)
    bar_width = min(42, group_width * 0.24)
    colors = {"simple": "#2563eb", "projection": "#0f766e"}

    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#f8fafc" />',
        f'<text x="{width / 2:.1f}" y="34" text-anchor="middle" font-size="24" font-weight="700" fill="#0f172a">'
        'SIMPLE vs Projection Iteration Counts</text>',
        f'<text x="{width / 2:.1f}" y="58" text-anchor="middle" font-size="13" fill="#475569">'
        'Lid-driven cavity, Re=100, same convergence gates</text>',
    ]

    for tick_index in range(6):
        value = y_max * tick_index / 5
        y = baseline - (value / y_max) * plot_height
        elements.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_width}" y2="{y:.1f}" '
            'stroke="#d8dee9" stroke-width="1" />'
        )
        elements.append(
            f'<text x="{left - 12}" y="{y + 4:.1f}" text-anchor="end" font-size="12" fill="#475569">'
            f"{int(value)}</text>"
        )

    elements.append(
        f'<line x1="{left}" y1="{baseline}" x2="{left + plot_width}" y2="{baseline}" stroke="#334155" stroke-width="1.5" />'
    )
    elements.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{baseline}" stroke="#334155" stroke-width="1.5" />'
    )

    for grid_index, grid in enumerate(grids):
        center = left + group_width * (grid_index + 0.5)
        simple = by_grid[grid]["simple"]
        projection = by_grid[grid]["projection"]
        speedup = simple["iterations"] / projection["iterations"]

        for offset, solver in [(-0.55, "simple"), (0.55, "projection")]:
            row = by_grid[grid][solver]
            bar_height = row["iterations"] / y_max * plot_height
            x = center + offset * bar_width - bar_width / 2
            y = baseline - bar_height
            elements.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" '
                f'rx="4" fill="{colors[solver]}" />'
            )
            elements.append(
                f'<text x="{x + bar_width / 2:.1f}" y="{y - 7:.1f}" text-anchor="middle" '
                'font-size="12" font-weight="600" fill="#0f172a">'
                f'{row["iterations"]}</text>'
            )

        elements.append(
            f'<text x="{center:.1f}" y="{baseline + 26}" text-anchor="middle" font-size="13" fill="#334155">'
            f'{grid}x{grid}</text>'
        )
        elements.append(
            f'<text x="{center:.1f}" y="{baseline + 45}" text-anchor="middle" font-size="12" fill="#64748b">'
            f'S/P = {speedup:.2f}x</text>'
        )

    legend_x = left + plot_width - 220
    legend_y = 82
    for idx, solver in enumerate(SOLVERS):
        y = legend_y + idx * 24
        elements.append(
            f'<rect x="{legend_x}" y="{y - 12}" width="18" height="18" rx="3" fill="{colors[solver]}" />'
        )
        elements.append(
            f'<text x="{legend_x + 28}" y="{y + 2}" font-size="13" fill="#0f172a">'
            f'{svg_escape(solver)}</text>'
        )

    elements.append(
        f'<text x="{left + plot_width / 2:.1f}" y="{height - 32}" text-anchor="middle" font-size="13" fill="#334155">'
        'Grid size</text>'
    )
    elements.append(
        f'<text x="24" y="{top + plot_height / 2:.1f}" text-anchor="middle" font-size="13" fill="#334155" '
        f'transform="rotate(-90 24 {top + plot_height / 2:.1f})">Iterations to convergence</text>'
    )
    elements.append("</svg>")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(elements))


def main():
    args = parse_args()
    solver_path = Path(args.solver_path)
    output_dir = Path(args.output_dir)
    if not solver_path.exists():
        raise SystemExit(f"Solver executable not found: {solver_path}")

    rows = []
    for grid in args.grids:
        for solver in SOLVERS:
            rows.append(
                run_case(
                    solver_path,
                    output_dir,
                    solver,
                    grid,
                    args.re,
                    args.min_iters,
                    args.simple_max_iters,
                    args.projection_max_iters,
                    args.projection_dt,
                )
            )

    csv_path = output_dir / "iteration_counts.csv"
    svg_path = output_dir / "iteration_comparison.svg"
    write_csv(rows, csv_path)
    draw_svg(rows, svg_path)

    print(f"Wrote {csv_path}")
    print(f"Wrote {svg_path}")
    for grid in sorted({row["nx"] for row in rows}):
        simple = next(row for row in rows if row["nx"] == grid and row["solver"] == "simple")
        projection = next(row for row in rows if row["nx"] == grid and row["solver"] == "projection")
        ratio = simple["iterations"] / projection["iterations"]
        print(
            f"{grid}x{grid}: simple={simple['iterations']} iterations, "
            f"projection={projection['iterations']} iterations, ratio={ratio:.2f}x"
        )


if __name__ == "__main__":
    main()
