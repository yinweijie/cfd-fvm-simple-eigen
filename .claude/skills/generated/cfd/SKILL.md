---
name: cfd
description: "Skill for the Cfd area of cfd-fvm-simple-eigen. 51 symbols across 12 files."
---

# Cfd

51 symbols | 12 files | Cohesion: 83%

## When to Use

- Working with code in `src/`
- Understanding how flow_solver_kind_name, extract_u_centerline, extract_v_centerline work
- Modifying cfd-related functionality

## Key Files

| File | Symbols |
|------|---------|
| `src/discretization.cpp` | pressure_gradient_x, pressure_gradient_y, rhie_chow_u_face_velocity, rhie_chow_v_face_velocity, positive_part (+9) |
| `src/output.cpp` | write_pairs, write_matrix_csv, pressure_cells, u_cells, v_cells (+3) |
| `include/cfd/discretization.hpp` | apply_cavity_boundary_conditions, assemble_u_projection_momentum, assemble_v_projection_momentum, assemble_pressure_correction, update_face_velocities (+3) |
| `src/simple_solver.cpp` | load_u_solution, load_v_solution, solve_u_predictor, solve_v_predictor, update_u_correction_factors (+1) |
| `src/projection_solver.cpp` | load_u_solution, load_v_solution, solve_u_predictor, solve_v_predictor |
| `tests/test_projection_coefficients.cpp` | expect_near, sparse_coeff, fields |
| `tests/test_coefficients.cpp` | expect_near, sparse_coeff, fields |
| `include/cfd/solver.hpp` | flow_solver_kind_name |
| `include/cfd/case.hpp` | viscosity |
| `include/cfd/linear_system.hpp` | solve_linear_system |

## Entry Points

Start here when exploring this area:

- **`flow_solver_kind_name`** (Function) — `include/cfd/solver.hpp:13`
- **`extract_u_centerline`** (Function) — `src/output.cpp:77`
- **`extract_v_centerline`** (Function) — `src/output.cpp:102`
- **`out_dir`** (Function) — `src/output.cpp:133`
- **`assemble_u_momentum`** (Function) — `src/discretization.cpp:308`

## Key Symbols

| Symbol | Type | File | Line |
|--------|------|------|------|
| `flow_solver_kind_name` | Function | `include/cfd/solver.hpp` | 13 |
| `extract_u_centerline` | Function | `src/output.cpp` | 77 |
| `extract_v_centerline` | Function | `src/output.cpp` | 102 |
| `out_dir` | Function | `src/output.cpp` | 133 |
| `assemble_u_momentum` | Function | `src/discretization.cpp` | 308 |
| `assemble_v_momentum` | Function | `src/discretization.cpp` | 392 |
| `assemble_u_projection_momentum` | Function | `src/discretization.cpp` | 474 |
| `assemble_v_projection_momentum` | Function | `src/discretization.cpp` | 550 |
| `solve_linear_system` | Function | `include/cfd/linear_system.hpp` | 30 |
| `apply_cavity_boundary_conditions` | Function | `include/cfd/discretization.hpp` | 21 |
| `assemble_u_projection_momentum` | Function | `include/cfd/discretization.hpp` | 51 |
| `assemble_v_projection_momentum` | Function | `include/cfd/discretization.hpp` | 57 |
| `assemble_pressure_correction` | Function | `include/cfd/discretization.hpp` | 63 |
| `fields` | Function | `tests/test_projection_coefficients.cpp` | 27 |
| `update_face_velocities` | Function | `include/cfd/discretization.hpp` | 27 |
| `assemble_u_momentum` | Function | `include/cfd/discretization.hpp` | 39 |
| `assemble_v_momentum` | Function | `include/cfd/discretization.hpp` | 45 |
| `compute_continuity_residual` | Function | `include/cfd/discretization.hpp` | 69 |
| `fields` | Function | `tests/test_coefficients.cpp` | 29 |
| `write_results` | Function | `include/cfd/output.hpp` | 13 |

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

## How to Explore

1. `gitnexus_context({name: "flow_solver_kind_name"})` — see callers and callees
2. `gitnexus_query({query: "cfd"})` — find related execution flows
3. Read key files listed above for implementation details
