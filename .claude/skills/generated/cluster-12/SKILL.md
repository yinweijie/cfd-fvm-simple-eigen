---
name: cluster-12
description: "Skill for the Cluster_12 area of cfd-fvm-simple-eigen. 12 symbols across 1 files."
---

# Cluster_12

12 symbols | 1 files | Cohesion: 84%

## When to Use

- Working with code in `src/`
- Understanding how ProjectionSolver, apply_boundary_conditions, set_projection_correction_factors work
- Modifying cluster_12-related functionality

## Key Files

| File | Symbols |
|------|---------|
| `src/projection_solver.cpp` | ProjectionSolver, apply_boundary_conditions, set_projection_correction_factors, prepare_iteration, refresh_face_velocities (+7) |

## Entry Points

Start here when exploring this area:

- **`ProjectionSolver`** (Method) — `src/projection_solver.cpp:10`
- **`apply_boundary_conditions`** (Method) — `src/projection_solver.cpp:21`
- **`set_projection_correction_factors`** (Method) — `src/projection_solver.cpp:26`
- **`prepare_iteration`** (Method) — `src/projection_solver.cpp:33`
- **`refresh_face_velocities`** (Method) — `src/projection_solver.cpp:40`

## Key Symbols

| Symbol | Type | File | Line |
|--------|------|------|------|
| `ProjectionSolver` | Method | `src/projection_solver.cpp` | 10 |
| `apply_boundary_conditions` | Method | `src/projection_solver.cpp` | 21 |
| `set_projection_correction_factors` | Method | `src/projection_solver.cpp` | 26 |
| `prepare_iteration` | Method | `src/projection_solver.cpp` | 33 |
| `refresh_face_velocities` | Method | `src/projection_solver.cpp` | 40 |
| `load_pressure_correction` | Method | `src/projection_solver.cpp` | 74 |
| `solve_pressure_increment_step` | Method | `src/projection_solver.cpp` | 105 |
| `project_pressure_and_velocity` | Method | `src/projection_solver.cpp` | 116 |
| `compute_max_velocity_correction` | Method | `src/projection_solver.cpp` | 141 |
| `project_fields_and_collect_metrics` | Method | `src/projection_solver.cpp` | 157 |
| `has_converged` | Method | `src/projection_solver.cpp` | 179 |
| `run` | Method | `src/projection_solver.cpp` | 191 |

## Execution Flows

| Flow | Type | Steps |
|------|------|-------|
| `Run → U_cells` | cross_community | 5 |
| `Run → Pressure_gradient_x` | cross_community | 5 |
| `Run → V_cells` | cross_community | 5 |
| `ProjectionSolver → U_cells` | cross_community | 5 |
| `ProjectionSolver → Pressure_gradient_x` | cross_community | 5 |
| `ProjectionSolver → V_cells` | cross_community | 5 |
| `ProjectionSolver → Pressure_gradient_y` | cross_community | 5 |
| `Project_fields_and_collect_metrics → U_cells` | cross_community | 5 |
| `Project_fields_and_collect_metrics → Pressure_gradient_x` | cross_community | 5 |
| `Project_fields_and_collect_metrics → V_cells` | cross_community | 5 |

## Connected Areas

| Area | Connections |
|------|-------------|
| Cluster_6 | 3 calls |
| Cfd | 3 calls |
| Cluster_8 | 1 calls |

## How to Explore

1. `gitnexus_context({name: "ProjectionSolver"})` — see callers and callees
2. `gitnexus_query({query: "cluster_12"})` — find related execution flows
3. Read key files listed above for implementation details
