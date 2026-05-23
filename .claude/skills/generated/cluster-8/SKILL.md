---
name: cluster-8
description: "Skill for the Cluster_8 area of cfd-fvm-simple-eigen. 12 symbols across 2 files."
---

# Cluster_8

12 symbols | 2 files | Cohesion: 83%

## When to Use

- Working with code in `src/`
- Understanding how apply_cavity_boundary_conditions, SimpleSolver, apply_boundary_conditions work
- Modifying cluster_8-related functionality

## Key Files

| File | Symbols |
|------|---------|
| `src/simple_solver.cpp` | SimpleSolver, apply_boundary_conditions, prepare_iteration, refresh_face_velocities, load_pressure_correction (+6) |
| `src/discretization.cpp` | apply_cavity_boundary_conditions |

## Entry Points

Start here when exploring this area:

- **`apply_cavity_boundary_conditions`** (Function) — `src/discretization.cpp:215`
- **`SimpleSolver`** (Method) — `src/simple_solver.cpp:9`
- **`apply_boundary_conditions`** (Method) — `src/simple_solver.cpp:16`
- **`prepare_iteration`** (Method) — `src/simple_solver.cpp:22`
- **`refresh_face_velocities`** (Method) — `src/simple_solver.cpp:29`

## Key Symbols

| Symbol | Type | File | Line |
|--------|------|------|------|
| `apply_cavity_boundary_conditions` | Function | `src/discretization.cpp` | 215 |
| `SimpleSolver` | Method | `src/simple_solver.cpp` | 9 |
| `apply_boundary_conditions` | Method | `src/simple_solver.cpp` | 16 |
| `prepare_iteration` | Method | `src/simple_solver.cpp` | 22 |
| `refresh_face_velocities` | Method | `src/simple_solver.cpp` | 29 |
| `load_pressure_correction` | Method | `src/simple_solver.cpp` | 68 |
| `solve_pressure_correction_step` | Method | `src/simple_solver.cpp` | 109 |
| `correct_pressure_and_velocity` | Method | `src/simple_solver.cpp` | 147 |
| `compute_max_velocity_correction` | Method | `src/simple_solver.cpp` | 176 |
| `correct_fields_and_collect_metrics` | Method | `src/simple_solver.cpp` | 195 |
| `has_converged` | Method | `src/simple_solver.cpp` | 217 |
| `run` | Method | `src/simple_solver.cpp` | 232 |

## Execution Flows

| Flow | Type | Steps |
|------|------|-------|
| `Correct_fields_and_collect_metrics → U_cells` | cross_community | 5 |
| `Correct_fields_and_collect_metrics → Pressure_gradient_x` | cross_community | 5 |
| `Correct_fields_and_collect_metrics → V_cells` | cross_community | 5 |
| `Correct_fields_and_collect_metrics → Pressure_gradient_y` | cross_community | 5 |
| `SimpleSolver → U_cells` | cross_community | 5 |
| `SimpleSolver → Pressure_gradient_x` | cross_community | 5 |
| `SimpleSolver → V_cells` | cross_community | 5 |
| `SimpleSolver → Pressure_gradient_y` | cross_community | 5 |
| `Solve_pressure_correction_step → U_cells` | cross_community | 5 |
| `Solve_pressure_correction_step → Pressure_gradient_x` | cross_community | 5 |

## Connected Areas

| Area | Connections |
|------|-------------|
| Cluster_6 | 3 calls |
| Cfd | 3 calls |

## How to Explore

1. `gitnexus_context({name: "apply_cavity_boundary_conditions"})` — see callers and callees
2. `gitnexus_query({query: "cluster_8"})` — find related execution flows
3. Read key files listed above for implementation details
