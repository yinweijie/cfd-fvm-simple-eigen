---
name: scripts
description: "Skill for the Scripts area of cfd-fvm-simple-eigen. 21 symbols across 2 files."
---

# Scripts

21 symbols | 2 files | Cohesion: 96%

## When to Use

- Working with code in `scripts/`
- Understanding how svg_escape, scale_point, polyline_points work
- Modifying scripts-related functionality

## Key Files

| File | Symbols |
|------|---------|
| `scripts/validate_cavity.py` | svg_escape, scale_point, polyline_points, draw_ticks, draw_series_points (+8) |
| `scripts/compare_solver_iterations.py` | parse_args, read_key_values, read_last_residual, run_case, write_csv (+3) |

## Entry Points

Start here when exploring this area:

- **`svg_escape`** (Function) — `scripts/validate_cavity.py:54`
- **`scale_point`** (Function) — `scripts/validate_cavity.py:64`
- **`polyline_points`** (Function) — `scripts/validate_cavity.py:72`
- **`draw_ticks`** (Function) — `scripts/validate_cavity.py:80`
- **`draw_series_points`** (Function) — `scripts/validate_cavity.py:116`

## Key Symbols

| Symbol | Type | File | Line |
|--------|------|------|------|
| `svg_escape` | Function | `scripts/validate_cavity.py` | 54 |
| `scale_point` | Function | `scripts/validate_cavity.py` | 64 |
| `polyline_points` | Function | `scripts/validate_cavity.py` | 72 |
| `draw_ticks` | Function | `scripts/validate_cavity.py` | 80 |
| `draw_series_points` | Function | `scripts/validate_cavity.py` | 116 |
| `draw_series_line` | Function | `scripts/validate_cavity.py` | 124 |
| `draw_legend_entry` | Function | `scripts/validate_cavity.py` | 132 |
| `draw_panel` | Function | `scripts/validate_cavity.py` | 148 |
| `write_validation_plot` | Function | `scripts/validate_cavity.py` | 260 |
| `parse_args` | Function | `scripts/compare_solver_iterations.py` | 11 |
| `read_key_values` | Function | `scripts/compare_solver_iterations.py` | 50 |
| `read_last_residual` | Function | `scripts/compare_solver_iterations.py` | 62 |
| `run_case` | Function | `scripts/compare_solver_iterations.py` | 70 |
| `write_csv` | Function | `scripts/compare_solver_iterations.py` | 124 |
| `svg_escape` | Function | `scripts/compare_solver_iterations.py` | 145 |
| `draw_svg` | Function | `scripts/compare_solver_iterations.py` | 155 |
| `main` | Function | `scripts/compare_solver_iterations.py` | 261 |
| `read_pairs` | Function | `scripts/validate_cavity.py` | 9 |
| `interpolate` | Function | `scripts/validate_cavity.py` | 23 |
| `compute_errors` | Function | `scripts/validate_cavity.py` | 39 |

## Execution Flows

| Flow | Type | Steps |
|------|------|-------|
| `Main → Scale_point` | cross_community | 6 |
| `Main → Svg_escape` | cross_community | 4 |
| `Main → Read_key_values` | intra_community | 3 |
| `Main → Read_last_residual` | intra_community | 3 |
| `Main → Svg_escape` | intra_community | 3 |
| `Main → Interpolate` | intra_community | 3 |

## How to Explore

1. `gitnexus_context({name: "svg_escape"})` — see callers and callees
2. `gitnexus_query({query: "scripts"})` — find related execution flows
3. Read key files listed above for implementation details
