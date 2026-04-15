#include "cfd/simple_solver.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace cfd {

// Initialize the collocated fields and their boundary-consistent face velocities.
SimpleSolver::SimpleSolver(CavityCase config)
    : grid_(config.mesh_spec), config_(std::move(config)), fields_(grid_) {
  apply_boundary_conditions();
  refresh_face_velocities(fields_.u, fields_.v);
}

// Apply the cavity ghost-cell boundary conditions to all working fields.
void SimpleSolver::apply_boundary_conditions() {
  apply_cavity_boundary_conditions(grid_, config_, &fields_);
}

// Refresh boundary conditions and corrected-field face velocities for a new iteration.
// This makes the predictor solve start from a self-consistent corrected state.
void SimpleSolver::prepare_iteration() {
  apply_boundary_conditions();
  refresh_face_velocities(fields_.u, fields_.v);
}

// Rebuild face velocities from the supplied cell-centered velocity fields.
// The caller chooses whether the reconstruction should use corrected or predicted velocities.
void SimpleSolver::refresh_face_velocities(
    const Eigen::MatrixXd& u_cells,
    const Eigen::MatrixXd& v_cells) {
  update_face_velocities(
      grid_,
      config_,
      u_cells,
      v_cells,
      fields_.pressure,
      fields_.d_u,
      fields_.d_v,
      &fields_.u_face,
      &fields_.v_face);
}

// Scatter the solved u predictor into the collocated u_star storage.
// The loop walks every interior control volume and uses the grid index mapping to
// copy the flat linear-solver result back into the 2D field storage.
void SimpleSolver::load_u_solution(const Eigen::VectorXd& values) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.u_star(i, j) = values(static_cast<int>(grid_.u_index(i, j)));
    }
  }
}

// Scatter the solved v predictor into the collocated v_star storage.
// The loop walks every interior control volume and uses the grid index mapping to
// copy the flat linear-solver result back into the 2D field storage.
void SimpleSolver::load_v_solution(const Eigen::VectorXd& values) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.v_star(i, j) = values(static_cast<int>(grid_.v_index(i, j)));
    }
  }
}

// Scatter the solved pressure-correction field into the collocated storage.
// The field is reset first so the pinned reference row and ghost cells do not keep stale values.
void SimpleSolver::load_pressure_correction(const Eigen::VectorXd& values) {
  fields_.pressure_correction.setZero();
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.pressure_correction(i, j) =
          values(static_cast<int>(grid_.pressure_index(i, j)));
    }
  }
}

// Assemble and solve the relaxed u-momentum predictor system.
// After the solve, the interior predictor is scattered into u_star and the relaxed
// diagonal is converted into d_u for the later Rhie-Chow and correction steps.
LinearSolveResult SimpleSolver::solve_u_predictor() {
  const MomentumAssembly assembly = assemble_u_momentum(grid_, config_, fields_);
  const LinearSolveResult result =
      solve_linear_system(assembly.system, SolverKind::kBiCGSTAB, 1e-8, 10000);
  // Keep the previous ghost layer until the next boundary-condition refresh.
  fields_.u_star = fields_.u;
  load_u_solution(result.solution);
  update_u_correction_factors(assembly.diagonal);
  return result;
}

// Assemble and solve the relaxed v-momentum predictor system.
// After the solve, the interior predictor is scattered into v_star and the relaxed
// diagonal is converted into d_v for the later Rhie-Chow and correction steps.
LinearSolveResult SimpleSolver::solve_v_predictor() {
  const MomentumAssembly assembly = assemble_v_momentum(grid_, config_, fields_);
  const LinearSolveResult result =
      solve_linear_system(assembly.system, SolverKind::kBiCGSTAB, 1e-8, 10000);
  // Keep the previous ghost layer until the next boundary-condition refresh.
  fields_.v_star = fields_.v;
  load_v_solution(result.solution);
  update_v_correction_factors(assembly.diagonal);
  return result;
}

// Assemble and solve the SIMPLE pressure-correction system.
// The resulting flat solution is scattered back into pressure_correction for the
// subsequent pressure and velocity correction stage.
LinearSolveResult SimpleSolver::solve_pressure_correction_step() {
  const PressureCorrectionAssembly assembly =
      assemble_pressure_correction(grid_, config_, fields_);
  const LinearSolveResult result =
      solve_linear_system(assembly.system, SolverKind::kConjugateGradient, 1e-10, 10000);
  load_pressure_correction(result.solution);
  return result;
}

// Convert the relaxed u diagonal into SIMPLE correction factors.
// The loop visits every interior cell and stores the local pressure-to-velocity
// correction scale used by both Rhie-Chow reconstruction and the final correction.
void SimpleSolver::update_u_correction_factors(const Eigen::VectorXd& diagonal) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      // SIMPLE uses the relaxed momentum diagonal to scale the pressure correction.
      fields_.d_u(i, j) = (grid_.dx() * grid_.dy()) / diagonal(static_cast<int>(grid_.u_index(i, j)));
    }
  }
}

// Convert the relaxed v diagonal into SIMPLE correction factors.
// The loop visits every interior cell and stores the local pressure-to-velocity
// correction scale used by both Rhie-Chow reconstruction and the final correction.
void SimpleSolver::update_v_correction_factors(const Eigen::VectorXd& diagonal) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      // SIMPLE uses the relaxed momentum diagonal to scale the pressure correction.
      fields_.d_v(i, j) = (grid_.dx() * grid_.dy()) / diagonal(static_cast<int>(grid_.v_index(i, j)));
    }
  }
}

// Apply the solved pressure increment to the cell-centered velocity and pressure fields.
// The first grid sweep computes central pressure-correction gradients and updates u, v,
// and p cell-by-cell; the second sweep removes the arbitrary pressure constant by
// subtracting the first-cell reference from every interior pressure value.
void SimpleSolver::correct_pressure_and_velocity() {
  // Convert the solved pressure increment into corrected cell-centered velocity and pressure fields.
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      const double grad_pc_x =
          (fields_.pressure_correction(i + 1, j) - fields_.pressure_correction(i - 1, j)) /
          (2.0 * grid_.dx());
      const double grad_pc_y =
          (fields_.pressure_correction(i, j + 1) - fields_.pressure_correction(i, j - 1)) /
          (2.0 * grid_.dy());

      fields_.u(i, j) = fields_.u_star(i, j) - fields_.d_u(i, j) * grad_pc_x;
      fields_.v(i, j) = fields_.v_star(i, j) - fields_.d_v(i, j) * grad_pc_y;
      fields_.pressure(i, j) += config_.controls.alpha_p * fields_.pressure_correction(i, j);
    }
  }

  // Pressure is defined up to a constant; pin the corrected field to the first cell.
  const double reference_pressure = fields_.pressure(1, 1);
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.pressure(i, j) -= reference_pressure;
    }
  }
}

// Measure the maximum cell-wise velocity change produced by the correction step.
// The loop compares each corrected cell against its pre-correction value and keeps
// the largest absolute update across both velocity components.
double SimpleSolver::compute_max_velocity_correction(
    const Eigen::MatrixXd& u_before,
    const Eigen::MatrixXd& v_before) const {
  double max_velocity_correction = 0.0;
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      max_velocity_correction = std::max(
          max_velocity_correction,
          std::max(std::abs(fields_.u(i, j) - u_before(i, j)),
                   std::abs(fields_.v(i, j) - v_before(i, j))));
    }
  }
  return max_velocity_correction;
}

// Finalize one SIMPLE iteration and collect residual-style diagnostics.
// This stage applies the correction, restores boundary-consistent ghost cells, rebuilds
// corrected face velocities, and then packages both linear-solver residuals and field
// change metrics into one IterationMetrics record.
IterationMetrics SimpleSolver::correct_fields_and_collect_metrics(
    int iteration,
    const LinearSolveResult& u_result,
    const LinearSolveResult& v_result,
    const LinearSolveResult& p_result,
    const Eigen::MatrixXd& u_before,
    const Eigen::MatrixXd& v_before) {
  correct_pressure_and_velocity();
  apply_boundary_conditions();
  refresh_face_velocities(fields_.u, fields_.v);

  IterationMetrics metrics;
  metrics.iteration = iteration;
  metrics.continuity_residual = compute_continuity_residual(grid_, config_, fields_);
  metrics.u_momentum_residual = u_result.relative_residual;
  metrics.v_momentum_residual = v_result.relative_residual;
  metrics.pressure_correction_residual = p_result.relative_residual;
  metrics.max_velocity_correction = compute_max_velocity_correction(u_before, v_before);
  return metrics;
}

// Check the momentum, continuity, and velocity-correction convergence gates.
bool SimpleSolver::has_converged(int iteration, const IterationMetrics& metrics) const {
  const bool momentum_ok =
      metrics.u_momentum_residual < config_.controls.momentum_tolerance &&
      metrics.v_momentum_residual < config_.controls.momentum_tolerance;
  const bool continuity_ok =
      metrics.continuity_residual < config_.controls.continuity_tolerance;
  const bool velocity_ok = metrics.max_velocity_correction < 1e-8;
  return iteration >= config_.controls.min_iterations && momentum_ok && continuity_ok && velocity_ok;
}

// Execute the full SIMPLE loop until convergence or the iteration limit is reached.
// Each iteration refreshes the corrected state, solves u and v predictors in sequence,
// builds the pressure correction, applies the correction, records diagnostics, and
// finally checks the convergence gates.
SolveSummary SimpleSolver::run() {
  SolveSummary summary;

  for (int iteration = 1; iteration <= config_.controls.max_iterations; ++iteration) {
    // SIMPLE predictor step: start from the corrected cell-centered fields.
    prepare_iteration();

    // Solve the u predictor first so its updated d_u and face fluxes can feed the v step.
    const LinearSolveResult u_result = solve_u_predictor();

    // Rebuild face velocities with the latest u predictor before solving v.
    refresh_face_velocities(fields_.u_star, fields_.v);
    const LinearSolveResult v_result = solve_v_predictor();

    // Pressure correction uses both predictors together with boundary-consistent ghost cells.
    apply_boundary_conditions();
    const LinearSolveResult p_result = solve_pressure_correction_step();

    // Compare corrected fields against the pre-correction state when forming diagnostics.
    const Eigen::MatrixXd u_before = fields_.u;
    const Eigen::MatrixXd v_before = fields_.v;
    const IterationMetrics metrics = correct_fields_and_collect_metrics(
        iteration, u_result, v_result, p_result, u_before, v_before);
    summary.residual_history.push_back(metrics);

    if (has_converged(iteration, metrics)) {
      summary.converged = true;
      summary.iterations = iteration;
      return summary;
    }
  }

  summary.iterations = static_cast<int>(summary.residual_history.size());
  return summary;
}

}  // namespace cfd
