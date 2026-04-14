#!/usr/bin/env python3

import argparse
import csv
import json
import math
from pathlib import Path


def read_pairs(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        if len(headers) != 2:
            raise ValueError(f"Expected two columns in {path}, got {headers}")
        x_key, y_key = headers
        for row in reader:
            rows.append((float(row[x_key]), float(row[y_key])))
    rows.sort(key=lambda item: item[0])
    return rows


def interpolate(series, x):
  if x <= series[0][0]:
      return series[0][1]
  if x >= series[-1][0]:
      return series[-1][1]

  for left, right in zip(series, series[1:]):
      x0, y0 = left
      x1, y1 = right
      if x0 <= x <= x1:
          t = (x - x0) / (x1 - x0)
          return y0 + t * (y1 - y0)

  raise RuntimeError(f"Interpolation failed for x={x}")


def compute_errors(reference, numerical):
    diffs = []
    for x_ref, y_ref in reference:
        y_num = interpolate(numerical, x_ref)
        diffs.append(y_num - y_ref)

    linf = max(abs(value) for value in diffs)
    l2 = math.sqrt(sum(value * value for value in diffs) / len(diffs))
    return {
        "linf": linf,
        "l2": l2,
        "samples": len(diffs),
    }


def svg_escape(text):
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def scale_point(x, y, x_min, x_max, y_min, y_max, left, top, width, height):
    x_norm = 0.0 if x_max == x_min else (x - x_min) / (x_max - x_min)
    y_norm = 0.0 if y_max == y_min else (y - y_min) / (y_max - y_min)
    px = left + x_norm * width
    py = top + (1.0 - y_norm) * height
    return px, py


def polyline_points(series, x_min, x_max, y_min, y_max, left, top, width, height):
    points = []
    for x, y in series:
        px, py = scale_point(x, y, x_min, x_max, y_min, y_max, left, top, width, height)
        points.append(f"{px:.2f},{py:.2f}")
    return " ".join(points)


def draw_ticks(elements, x_min, x_max, y_min, y_max, left, top, width, height):
    x_ticks = [0.0, 0.25, 0.5, 0.75, 1.0]
    y_span = y_max - y_min
    if y_span <= 0:
        y_ticks = [y_min]
    else:
        step = y_span / 4.0
        y_ticks = [y_min + step * i for i in range(5)]

    for tick in x_ticks:
        px, py0 = scale_point(tick, y_min, x_min, x_max, y_min, y_max, left, top, width, height)
        _, py1 = scale_point(tick, y_max, x_min, x_max, y_min, y_max, left, top, width, height)
        elements.append(
            f'<line x1="{px:.2f}" y1="{py0:.2f}" x2="{px:.2f}" y2="{py1:.2f}" '
            'stroke="#d8dde6" stroke-width="1" />'
        )
        elements.append(
            f'<text x="{px:.2f}" y="{top + height + 22:.2f}" text-anchor="middle" '
            'font-size="12" fill="#374151">'
            f'{tick:.2f}</text>'
        )

    for tick in y_ticks:
        px0, py = scale_point(x_min, tick, x_min, x_max, y_min, y_max, left, top, width, height)
        px1, _ = scale_point(x_max, tick, x_min, x_max, y_min, y_max, left, top, width, height)
        elements.append(
            f'<line x1="{px0:.2f}" y1="{py:.2f}" x2="{px1:.2f}" y2="{py:.2f}" '
            'stroke="#d8dde6" stroke-width="1" />'
        )
        elements.append(
            f'<text x="{left - 10:.2f}" y="{py + 4:.2f}" text-anchor="end" '
            'font-size="12" fill="#374151">'
            f'{tick:.2f}</text>'
        )


def draw_series_points(elements, series, color, x_min, x_max, y_min, y_max, left, top, width, height):
    for x, y in series:
        px, py = scale_point(x, y, x_min, x_max, y_min, y_max, left, top, width, height)
        elements.append(
            f'<circle cx="{px:.2f}" cy="{py:.2f}" r="3.5" fill="{color}" stroke="white" stroke-width="1" />'
        )


def draw_panel(elements, title, x_label, y_label, reference, numerical, errors, bounds):
    left, top, width, height = bounds
    x_min = 0.0
    x_max = 1.0
    y_values = [value for _, value in reference] + [value for _, value in numerical]
    y_pad = 0.08 * max(1e-12, max(y_values) - min(y_values))
    y_min = min(y_values) - y_pad
    y_max = max(y_values) + y_pad

    elements.append(
        f'<rect x="{left}" y="{top}" width="{width}" height="{height}" rx="10" '
        'fill="#ffffff" stroke="#cbd5e1" stroke-width="1.5" />'
    )
    elements.append(
        f'<text x="{left + 16}" y="{top + 26}" font-size="18" font-weight="700" fill="#0f172a">'
        f'{svg_escape(title)}</text>'
    )
    elements.append(
        f'<text x="{left + 16}" y="{top + 48}" font-size="12" fill="#475569">'
        f'Linf={errors["linf"]:.4f}, L2={errors["l2"]:.4f}, samples={errors["samples"]}</text>'
    )

    plot_left = left + 60
    plot_top = top + 70
    plot_width = width - 90
    plot_height = height - 115

    draw_ticks(elements, x_min, x_max, y_min, y_max, plot_left, plot_top, plot_width, plot_height)
    elements.append(
        f'<rect x="{plot_left}" y="{plot_top}" width="{plot_width}" height="{plot_height}" '
        'fill="none" stroke="#334155" stroke-width="1.5" />'
    )

    numerical_line = polyline_points(
        numerical, x_min, x_max, y_min, y_max, plot_left, plot_top, plot_width, plot_height
    )
    elements.append(
        f'<polyline points="{numerical_line}" fill="none" stroke="#0f766e" stroke-width="2.5" />'
    )
    draw_series_points(
        elements, reference, "#dc2626", x_min, x_max, y_min, y_max, plot_left, plot_top, plot_width, plot_height
    )

    legend_x = left + width - 175
    legend_y = top + 22
    elements.append(
        f'<rect x="{legend_x}" y="{legend_y}" width="155" height="44" rx="8" '
        'fill="#f8fafc" stroke="#cbd5e1" stroke-width="1" />'
    )
    elements.append(
        f'<line x1="{legend_x + 12}" y1="{legend_y + 14}" x2="{legend_x + 42}" y2="{legend_y + 14}" '
        'stroke="#0f766e" stroke-width="2.5" />'
    )
    elements.append(
        f'<text x="{legend_x + 50}" y="{legend_y + 18}" font-size="12" fill="#0f172a">Simulation</text>'
    )
    elements.append(
        f'<circle cx="{legend_x + 27}" cy="{legend_y + 31}" r="3.5" fill="#dc2626" stroke="white" stroke-width="1" />'
    )
    elements.append(
        f'<text x="{legend_x + 50}" y="{legend_y + 35}" font-size="12" fill="#0f172a">Benchmark</text>'
    )

    elements.append(
        f'<text x="{plot_left + plot_width / 2:.2f}" y="{top + height - 18:.2f}" text-anchor="middle" '
        f'font-size="13" fill="#334155">{svg_escape(x_label)}</text>'
    )
    elements.append(
        f'<text x="{left + 18:.2f}" y="{plot_top + plot_height / 2:.2f}" text-anchor="middle" '
        f'font-size="13" fill="#334155" transform="rotate(-90 {left + 18:.2f} {plot_top + plot_height / 2:.2f})">'
        f'{svg_escape(y_label)}</text>'
    )


def write_validation_plot(plot_path, results_dir, u_reference, v_reference, u_numerical, v_numerical, summary):
    width = 1280
    height = 900
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
    ]
    elements.append('<rect width="100%" height="100%" fill="#eef2ff" />')
    elements.append(
        '<text x="48" y="56" font-size="30" font-weight="700" fill="#0f172a">'
        'CFD Cavity Validation: Simulation vs Benchmark</text>'
    )
    elements.append(
        f'<text x="48" y="86" font-size="14" fill="#475569">Results: {svg_escape(results_dir)}</text>'
    )
    elements.append(
        f'<text x="48" y="110" font-size="14" fill="#475569">'
        f'Overall verdict: {"PASS" if summary["passed"] else "FAIL"}'
        '</text>'
    )

    draw_panel(
        elements,
        "u centerline at x = 0.5",
        "y",
        "u",
        u_reference,
        u_numerical,
        summary["u"],
        (40, 140, 1200, 330),
    )
    draw_panel(
        elements,
        "v centerline at y = 0.5",
        "x",
        "v",
        v_reference,
        v_numerical,
        summary["v"],
        (40, 510, 1200, 330),
    )
    elements.append("</svg>")
    plot_path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--u-benchmark", default="data/benchmarks/ghia_u_centerline.csv")
    parser.add_argument("--v-benchmark", default="data/benchmarks/ghia_v_centerline.csv")
    parser.add_argument("--u-threshold", type=float, default=0.15)
    parser.add_argument("--v-threshold", type=float, default=0.15)
    parser.add_argument("--plot-out", default=None)
    args = parser.parse_args()

    results_dir = Path(args.results)
    u_reference = read_pairs(args.u_benchmark)
    v_reference = read_pairs(args.v_benchmark)
    u_numerical = read_pairs(results_dir / "centerline_u.csv")
    v_numerical = read_pairs(results_dir / "centerline_v.csv")

    u_errors = compute_errors(u_reference, u_numerical)
    v_errors = compute_errors(v_reference, v_numerical)
    passed = u_errors["linf"] <= args.u_threshold and v_errors["linf"] <= args.v_threshold

    summary = {
        "passed": passed,
        "u": u_errors,
        "v": v_errors,
        "thresholds": {
            "u_linf": args.u_threshold,
            "v_linf": args.v_threshold,
        },
    }
    plot_path = Path(args.plot_out) if args.plot_out else results_dir / "validation_plot.svg"
    write_validation_plot(
        plot_path, results_dir, u_reference, v_reference, u_numerical, v_numerical, summary
    )
    summary["artifacts"] = {
        "plot": str(plot_path),
    }

    summary_path = results_dir / "validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
