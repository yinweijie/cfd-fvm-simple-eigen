---
name: cluster-6
description: "Skill for the Cluster_6 area of cfd-fvm-simple-eigen. 6 symbols across 1 files."
---

# Cluster_6

6 symbols | 1 files | Cohesion: 60%

## When to Use

- Working with code in `src/`
- Understanding how update_face_velocities, assemble_pressure_correction, compute_continuity_residual work
- Modifying cluster_6-related functionality

## Key Files

| File | Symbols |
|------|---------|
| `src/discretization.cpp` | a_p, compute_pressure_correction_coefficients, compute_mass_imbalance, update_face_velocities, assemble_pressure_correction (+1) |

## Entry Points

Start here when exploring this area:

- **`update_face_velocities`** (Function) — `src/discretization.cpp:267`
- **`assemble_pressure_correction`** (Function) — `src/discretization.cpp:625`
- **`compute_continuity_residual`** (Function) — `src/discretization.cpp:703`

## Key Symbols

| Symbol | Type | File | Line |
|--------|------|------|------|
| `update_face_velocities` | Function | `src/discretization.cpp` | 267 |
| `assemble_pressure_correction` | Function | `src/discretization.cpp` | 625 |
| `compute_continuity_residual` | Function | `src/discretization.cpp` | 703 |
| `compute_pressure_correction_coefficients` | Function | `src/discretization.cpp` | 130 |
| `compute_mass_imbalance` | Function | `src/discretization.cpp` | 148 |
| `a_p` | Method | `src/discretization.cpp` | 34 |

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
| Cfd | 2 calls |

## How to Explore

1. `gitnexus_context({name: "update_face_velocities"})` — see callers and callees
2. `gitnexus_query({query: "cluster_6"})` — find related execution flows
3. Read key files listed above for implementation details
