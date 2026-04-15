# SIMPLE Solver Theory Notes

This note documents the current lid-driven cavity solver implementation in code-first notation.
The goal is to let a reader move directly between the equations and the code without renaming symbols in their head.

## Problem Setup

- Domain: 2D square cavity with `lx = ly = 1`
- Mesh: uniform structured finite-volume grid with `nx * ny` interior control volumes
- Flow model: steady, incompressible, laminar
- Boundary condition:
  - top wall moves with `lid_velocity`
  - other walls are stationary
  - `pressure` and `pressure_correction` use zero normal gradient at walls

The code stores one ghost-cell layer around every cell-centered field.
That is why `pressure`, `pressure_correction`, `u`, `v`, `u_star`, and `v_star` all have size `(nx + 2, ny + 2)`.

## Primary Fields

The solver uses these code identifiers as the canonical notation:

- `pressure`: corrected cell-centered pressure
- `pressure_correction`: pressure increment field, often written as `p'` in textbooks
- `u`, `v`: corrected cell-centered velocity components
- `u_star`, `v_star`: momentum-predicted velocity components, often written as `u*`, `v*`
- `d_u`, `d_v`: SIMPLE velocity-correction factors derived from the relaxed momentum diagonal
- `u_face`, `v_face`: Rhie-Chow reconstructed face velocities on physical faces

## Boundary Conditions And Ghost Cells

`apply_cavity_boundary_conditions()` updates the ghost layer so the interior stencil sees the cavity walls through algebraic mirror rules:

- side walls:
  - `u(0, j) = -u(1, j)`
  - `u(nx + 1, j) = -u(nx, j)`
  - `v(0, j) = -v(1, j)`
  - `v(nx + 1, j) = -v(nx, j)`
- bottom wall:
  - `u(i, 0) = -u(i, 1)`
  - `v(i, 0) = -v(i, 1)`
- moving lid:
  - `u(i, ny + 1) = 2 * lid_velocity - u(i, ny)`
  - `v(i, ny + 1) = -v(i, ny)`
- pressure-like fields:
  - `pressure` and `pressure_correction` copy the adjacent interior value at every wall

The same rules are applied to both corrected and predicted velocity fields so that every SIMPLE stage sees the same wall treatment.

## How Boundary Terms Enter The Algebra

This solver does use ghost cells for all cell-centered fields.
The momentum equations are assembled only for the `nx * ny` interior control volumes, so ghost cells are never unknowns in the linear system.
Instead, the ghost-cell values are substituted into the stencil before or during assembly.

For a stationary velocity wall, the ghost rule is a mirror rule:

- `u_G = -u_P`
- `v_G = -v_P`

If a boundary-side stencil contribution would normally look like `a_B * phi_G`, then substituting the ghost relation gives:

`a_B * phi_G = a_B * (-phi_P) = -a_B * phi_P`

Moving that term to the left-hand side adds `+a_B` to the diagonal.
That is why stationary wall boundaries in the velocity equations do not create an extra constant source term in this code.
They only fold the missing neighbor coefficient back into `a_p`.

For the moving lid, the `u` ghost rule is:

- `u_G = 2 * U_lid - u_P`

Substituting that into the north-side contribution gives:

`a_N * u_G = a_N * (2 * U_lid - u_P) = -a_N * u_P + 2 * a_N * U_lid`

Again, `-a_N * u_P` is moved to the left-hand side, so `a_N` is added to the diagonal.
The remaining constant piece becomes the explicit source contribution `2 * a_N * U_lid`.

So the code is not sending a boundary coefficient into the source term arbitrarily.
It is eliminating the ghost-cell unknown algebraically.
The final assembled form is exactly the ghost-cell-substituted form of the same stencil.

## Momentum Equations

`assemble_u_momentum()` and `assemble_v_momentum()` build the cell-centered momentum systems using an upwind convection plus central diffusion stencil.

For one control volume, the stencil coefficients are:

- `a_e`, `a_w`, `a_n`, `a_s`: east, west, north, south neighbor coefficients
- `a_p_base = a_e + a_w + a_n + a_s + (fe - fw + fn - fs)`
- `a_p = a_p_base / alpha`

where:

- `fe`, `fw`, `fn`, `fs` are the face mass fluxes reconstructed from `u_face` and `v_face`
- `alpha` is `alpha_u` or `alpha_v`, depending on the momentum equation

The source term contains:

- the pressure force, stored by `pressure_force_u()` or `pressure_force_v()`
- the under-relaxation correction `((1 - alpha) / alpha) * a_p_base * velocity_old`

For the top wall in the `u` equation, the moving lid contributes the extra source term `2 * a_n * lid_velocity`.

When a cell touches a cavity wall, the missing neighbor link is not assembled as an off-diagonal term.
Instead, the corresponding coefficient is folded back into the diagonal `a_p`.
For the north wall in the `u` equation, that fold-back is paired with the extra moving-lid source term.
This is the algebraic result of substituting the ghost-cell relations described above.

## Discretization Coefficient Expressions (Typst)

The standalone Typst document for these coefficient expressions lives at [`docs/momentum_boundary_tables.typ`](docs/momentum_boundary_tables.typ).
It now includes full momentum right-hand-side expressions, momentum boundary tables, equivalent `S_P / S_U` linearization tables, Rhie-Chow face-velocity formulas, and the pressure-correction coefficient expressions in a single file that can be copied directly into Typst-based notes or papers.

## Rhie-Chow Face Velocities

The solver uses a collocated mesh, so face velocities are not taken from a separate staggered storage.
Instead, `update_face_velocities()` reconstructs them from cell-centered values with a Rhie-Chow-style correction:

- `u_face = interp(u_cells) - d_face * (dp_face - grad_avg)`
- `v_face = interp(v_cells) - d_face * (dp_face - grad_avg)`

More explicitly:

- `interp(...)` is the arithmetic average of the two neighboring cell-centered velocities
- `dp_face` is the one-sided pressure difference across the face
- `grad_avg` is the average of the adjacent cell-centered pressure gradients
- `d_face` is the average of `d_u` or `d_v` from the two neighboring cells

This keeps the pressure and velocity coupling stable on the collocated grid while preserving the expected behavior for linear pressure fields.

## Pressure Correction Equation

After the predicted velocities are available, `assemble_pressure_correction()` builds the SIMPLE pressure-correction system.

The right-hand side is the cell mass defect:

`mass_imbalance = rho * dy * (u_face_star(i - 1, j) - u_face_star(i, j)) + rho * dx * (v_face_star(i, j - 1) - v_face_star(i, j))`

The matrix coefficients are assembled from `d_u` and `d_v`:

- `a_e = rho * dy * avg(d_u) / dx`
- `a_w = rho * dy * avg(d_u) / dx`
- `a_n = rho * dx * avg(d_v) / dy`
- `a_s = rho * dx * avg(d_v) / dy`
- `a_p = a_e + a_w + a_n + a_s`

Because only pressure gradients matter, one cell is pinned as a reference:

- the `(1, 1)` row is replaced by `1 * pressure_correction(1, 1) = 0`

The pressure-like fields use a different wall treatment from velocity.
Their ghost rule is zero normal gradient:

- `pressure_G = pressure_P`
- `pressure_correction_G = pressure_correction_P`

So for pressure and pressure correction, the wall does not inject a Dirichlet value into the stencil the way the lid velocity does for `u`.
That is why the pressure-correction assembly in this code simply omits nonexistent outer neighbors and pins one reference cell, instead of adding a moving-wall-type source term.

## SIMPLE Iteration Flow

`SimpleSolver::run()` follows this sequence each iteration:

1. Apply boundary conditions and rebuild `u_face`, `v_face` from the current corrected fields.
2. Solve the `u` momentum equation to get `u_star`.
3. Update `d_u`, then rebuild face velocities using `u_star` and the current `v`.
4. Solve the `v` momentum equation to get `v_star`.
5. Update `d_v`, reapply boundary conditions, and solve the pressure-correction equation.
6. Correct the fields:
   - `u = u_star - d_u * grad(pressure_correction)`
   - `v = v_star - d_v * grad(pressure_correction)`
   - `pressure = pressure + alpha_p * pressure_correction`
7. Subtract `pressure(1, 1)` from every interior cell so the corrected pressure keeps a fixed reference level.
8. Reapply boundary conditions, rebuild face velocities, and evaluate convergence.

## Convergence Metrics

The solver records one `IterationMetrics` entry per iteration:

- `continuity_residual`: average absolute cell mass defect computed from the corrected `u_face` and `v_face`
- `u_momentum_residual`: linear-solver relative residual for the `u` system
- `v_momentum_residual`: linear-solver relative residual for the `v` system
- `pressure_correction_residual`: linear-solver relative residual for the pressure-correction system
- `max_velocity_correction`: max absolute change between the pre-correction and post-correction `u`, `v`

The run is considered converged only after `min_iterations` and only if momentum, continuity, and velocity-correction thresholds all pass.

## Output Quantities

`write_results()` exports:

- `u.csv`, `v.csv`, `p.csv`: interior cell-centered fields only
- `centerline_u.csv`: the `u` profile interpolated to `x = 0.5`
- `centerline_v.csv`: the `v` profile interpolated to `y = 0.5`
- `residuals.csv`: iteration history
- `summary.txt`: final convergence summary

The centerline files are interpolated because the cell centers on an even grid do not lie exactly on the geometric midlines.

## Symbol Mapping

| Code Identifier | Theory Meaning | Producer | Consumer | File / Function |
| --- | --- | --- | --- | --- |
| `pressure` | Corrected cell-centered pressure | `correct_pressure_and_velocity()` | momentum assembly, face reconstruction, output | `src/simple_solver.cpp`, `src/discretization.cpp`, `src/output.cpp` |
| `pressure_correction` | SIMPLE pressure increment (`p'`) | `load_pressure_correction()` | pressure/velocity correction, wall ghost update | `src/simple_solver.cpp`, `src/discretization.cpp` |
| `u` | Corrected cell-centered x-velocity | `correct_pressure_and_velocity()` | next SIMPLE iteration, output | `src/simple_solver.cpp`, `src/output.cpp` |
| `v` | Corrected cell-centered y-velocity | `correct_pressure_and_velocity()` | next SIMPLE iteration, output | `src/simple_solver.cpp`, `src/output.cpp` |
| `u_star` | Predicted x-velocity (`u*`) | `load_u_solution()` | `v` predictor, final correction | `src/simple_solver.cpp` |
| `v_star` | Predicted y-velocity (`v*`) | `load_v_solution()` | pressure correction, final correction | `src/simple_solver.cpp` |
| `d_u` | x-velocity correction factor from relaxed momentum diagonal | `update_u_correction_factors()` | Rhie-Chow reconstruction, final correction, pressure correction | `src/simple_solver.cpp`, `src/discretization.cpp` |
| `d_v` | y-velocity correction factor from relaxed momentum diagonal | `update_v_correction_factors()` | Rhie-Chow reconstruction, final correction, pressure correction | `src/simple_solver.cpp`, `src/discretization.cpp` |
| `u_face` | Rhie-Chow x-face velocity | `update_face_velocities()` | momentum fluxes, continuity residual, predictor-face reconstruction seed | `src/discretization.cpp` |
| `v_face` | Rhie-Chow y-face velocity | `update_face_velocities()` | momentum fluxes, continuity residual, predictor-face reconstruction seed | `src/discretization.cpp` |
| `a_e`, `a_w`, `a_n`, `a_s` | Neighbor stencil coefficients | `compute_momentum_coefficients()`, `compute_pressure_correction_coefficients()` | matrix assembly | `src/discretization.cpp` |
| `a_p_base` | Unrelaxed momentum diagonal | `compute_momentum_coefficients()` | under-relaxation source and relaxed diagonal | `src/discretization.cpp` |
| `a_p` | Relaxed diagonal / central coefficient | momentum and pressure-correction assembly | linear solves, SIMPLE correction factors | `src/discretization.cpp`, `src/simple_solver.cpp` |
| `mass_imbalance` | Cell continuity defect used as the pressure-correction RHS | `assemble_pressure_correction()` | diagnostics and pressure-correction solve | `src/discretization.cpp` |
