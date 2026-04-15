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

## Documentation Ownership Boundary

This Markdown note is the canonical home for implementation-stage meaning: SIMPLE chronology, function-to-stage mapping, correction formulas, pressure-reference handling, and the distinction between the gradients used in Rhie-Chow reconstruction and the gradients used in the final correction step.

The reusable discrete algebra lives in [`docs/momentum_boundary_tables.typ`](docs/momentum_boundary_tables.typ), which is the canonical home for coefficient formulas, boundary linearization tables, face-averaged `d_face` definitions, and the pressure-correction coefficient and `mass_imbalance` equations.

Use the split below to avoid drift:

| Topic | Canonical Doc | Why |
| --- | --- | --- |
| SIMPLE call order and stage semantics | `docs/simple_solver_theory.md` | tied directly to `SimpleSolver::run()` and helper sequencing |
| `pressure_gradient_x()` / `pressure_gradient_y()` in Rhie-Chow reconstruction | `docs/simple_solver_theory.md` | needs code-stage interpretation, not just algebra |
| `grad_pc_x` / `grad_pc_y` in final field correction | `docs/simple_solver_theory.md` | belongs to the correction stage, not the reusable stencil tables |
| Momentum and pressure-correction coefficient formulas | `docs/momentum_boundary_tables.typ` | reusable equation/table material |
| Boundary fold-back and `S_P / S_U` style tables | `docs/momentum_boundary_tables.typ` | reusable discrete-algebra reference |
| Face-averaged `d_face`, pressure-correction RHS, and coefficient definitions | `docs/momentum_boundary_tables.typ` | reusable notation for papers/notes |

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

After each momentum solve, `update_u_correction_factors()` or `update_v_correction_factors()` converts the relaxed diagonal into the stored SIMPLE correction factors:

- `d_u(i, j) = dx * dy / a_p^u(i, j)`
- `d_v(i, j) = dx * dy / a_p^v(i, j)`

These are not extra transport coefficients.
They are the local pressure-to-velocity response scales reused in both Rhie-Chow reconstruction and the later pressure-correction update of `u` and `v`.

## Discretization Coefficient Expressions (Typst)

The standalone Typst document for these coefficient expressions lives at [`docs/momentum_boundary_tables.typ`](docs/momentum_boundary_tables.typ).
It now includes full momentum right-hand-side expressions, momentum boundary tables, equivalent `S_P / S_U` linearization tables, Rhie-Chow face-velocity formulas, and the pressure-correction coefficient expressions in a single file that can be copied directly into Typst-based notes or papers.

## Rhie-Chow Face Velocities

The solver uses a collocated mesh, so face velocities are not taken from a separate staggered storage.
Instead, `update_face_velocities()` reconstructs them from cell-centered values with a Rhie-Chow-style correction.
Its interior-face helpers are `rhie_chow_u_face_velocity()` and `rhie_chow_v_face_velocity()`, which both call `pressure_gradient_x()` or `pressure_gradient_y()` on the current corrected `pressure` field:

- `u_face = interp(u_cells) - d_face * (dp_face - grad_avg)`
- `v_face = interp(v_cells) - d_face * (dp_face - grad_avg)`

More explicitly:

- `interp(...)` is the arithmetic average of the two neighboring cell-centered velocities
- `dp_face` is the one-sided pressure difference across the face
- `grad_avg` is the average of the adjacent cell-centered pressure gradients returned by `pressure_gradient_x()` or `pressure_gradient_y()`
- `d_face` is the average of `d_u` or `d_v` from the two neighboring cells

For the x-face reconstruction, the code uses:

- `pressure_gradient_x(i, j) = (pressure(i + 1, j) - pressure(i - 1, j)) / (2 * dx)`
- `u_face(i, j) = 0.5 * (u_cells(i, j) + u_cells(i + 1, j)) - 0.5 * (d_u(i, j) + d_u(i + 1, j)) * (dp_face - grad_avg)`

with `dp_face = (pressure(i + 1, j) - pressure(i, j)) / dx` and `grad_avg = 0.5 * (pressure_gradient_x(i, j) + pressure_gradient_x(i + 1, j))`.

The y-face reconstruction is the analogous `pressure_gradient_y()` / `v_face` version.
This gradient pair belongs specifically to face reconstruction.
It is distinct from the final correction gradients `grad_pc_x` and `grad_pc_y`, which are central differences of `pressure_correction` and are used only inside `correct_pressure_and_velocity()`.

This keeps the pressure and velocity coupling stable on the collocated grid while preserving the expected behavior for linear pressure fields.

## Pressure Correction Equation

After the predicted velocities are available, `assemble_pressure_correction()` builds the SIMPLE pressure-correction system.
It does not reuse the earlier `(u_star, v)` face refresh from `run()`.
Instead, it creates local `u_face_star` and `v_face_star` arrays and rebuilds predictor faces inside the assembly call with:

- `u_cells = u_star`
- `v_cells = v_star`
- the current corrected `pressure`
- the latest `d_u` and `d_v`

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

That row pinning is only the linear-system reference fix for `pressure_correction`.
The corrected physical pressure field is anchored again later by subtracting `pressure(1, 1)` from every interior cell after `correct_pressure_and_velocity()`.
Both mechanisms are required in the implementation and they act at different stages.

## Pressure And Velocity Correction Stage

Once `load_pressure_correction()` has scattered the solved increment field, `correct_pressure_and_velocity()` performs the cell-centered SIMPLE update:

- `grad_pc_x = (pressure_correction(i + 1, j) - pressure_correction(i - 1, j)) / (2 * dx)`
- `grad_pc_y = (pressure_correction(i, j + 1) - pressure_correction(i, j - 1)) / (2 * dy)`
- `u(i, j) = u_star(i, j) - d_u(i, j) * grad_pc_x`
- `v(i, j) = v_star(i, j) - d_v(i, j) * grad_pc_y`
- `pressure(i, j) = pressure(i, j) + alpha_p * pressure_correction(i, j)`

The names `grad_pc_x` and `grad_pc_y` are intentionally separate from `pressure_gradient_x()` and `pressure_gradient_y()`:

- `pressure_gradient_x()` / `pressure_gradient_y()` operate on the corrected pressure field and belong to Rhie-Chow face reconstruction.
- `grad_pc_x` / `grad_pc_y` operate on `pressure_correction` and belong to the final cell-centered correction stage.

After this update sweep, the solver removes the arbitrary constant level from the corrected pressure field by subtracting `pressure(1, 1)` from every interior cell.
So the implementation uses two pressure-reference mechanisms:

1. a pinned `(1, 1)` pressure-correction row during assembly
2. a post-correction pressure shift on the corrected `pressure` field

## SIMPLE Iteration Flow

`SimpleSolver::run()` follows this sequence each iteration:

1. `prepare_iteration()` applies cavity ghost-cell boundary conditions and rebuilds `u_face`, `v_face` from the current corrected fields `u`, `v`.
2. `solve_u_predictor()` assembles and solves the `u` momentum system, scatters the interior result into `u_star`, and updates `d_u` from the relaxed momentum diagonal.
3. `run()` immediately calls `refresh_face_velocities(u_star, v)` so the `v` predictor sees the newest x-face fluxes together with the current corrected `v`.
4. `solve_v_predictor()` assembles and solves the `v` momentum system, scatters the interior result into `v_star`, and updates `d_v`.
5. `run()` reapplies `apply_boundary_conditions()` before pressure correction so both predictor fields have boundary-consistent ghost cells.
6. `solve_pressure_correction_step()` calls `assemble_pressure_correction()`, which rebuilds predictor faces internally from `(u_star, v_star)` using the current `pressure`, `d_u`, and `d_v`, then solves for `pressure_correction`.
7. `correct_fields_and_collect_metrics()` calls `correct_pressure_and_velocity()`, reapplies boundary conditions, rebuilds corrected face velocities from `(u, v)`, and only then forms the continuity residual and other iteration metrics.

Two face-velocity refreshes therefore matter for different reasons:

- `refresh_face_velocities(u_star, v)` is the staging bridge between the `u` and `v` predictor solves.
- the predictor-face reconstruction inside `assemble_pressure_correction()` is the one that actually feeds the pressure-correction right-hand side, and it uses `(u_star, v_star)` rather than `(u_star, v)`.

## Convergence Metrics

The solver records one `IterationMetrics` entry per iteration:

- `continuity_residual`: average absolute cell mass defect computed after rebuilding corrected face velocities from `u`, `v`, `pressure`, `d_u`, and `d_v`
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
