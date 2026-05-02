#include "cfd/projection_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfd {

// Initialize the collocated fields for pseudo-transient projection iterations.
ProjectionSolver::ProjectionSolver(CavityCase config)
    : grid_(config.mesh_spec), config_(std::move(config)), fields_(grid_) {
  if (config_.controls.projection_dt <= 0.0) {
    throw std::invalid_argument("Projection time step must be positive");
  }
  set_projection_correction_factors();
  apply_boundary_conditions();
  refresh_face_velocities(fields_.u, fields_.v);
}

// Apply cavity ghost-cell boundary conditions to all working fields.
void ProjectionSolver::apply_boundary_conditions() {
  apply_cavity_boundary_conditions(grid_, config_, &fields_);
}

// Store the projection velocity response d = dt / rho used by the pressure increment.
void ProjectionSolver::set_projection_correction_factors() {
  const double velocity_response = config_.controls.projection_dt / config_.density;
  fields_.d_u.setConstant(velocity_response);
  fields_.d_v.setConstant(velocity_response);
}

// Refresh the corrected state at the start of one projection iteration.
void ProjectionSolver::prepare_iteration() {
  set_projection_correction_factors();
  apply_boundary_conditions();
  refresh_face_velocities(fields_.u, fields_.v);
}

// Rebuild face velocities from the supplied cell-centered velocity fields.
void ProjectionSolver::refresh_face_velocities(
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

// Scatter the flat u predictor solution into cell-centered storage.
void ProjectionSolver::load_u_solution(const Eigen::VectorXd& values) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.u_star(i, j) = values(static_cast<int>(grid_.u_index(i, j)));
    }
  }
}

// Scatter the flat v predictor solution into cell-centered storage.
void ProjectionSolver::load_v_solution(const Eigen::VectorXd& values) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.v_star(i, j) = values(static_cast<int>(grid_.v_index(i, j)));
    }
  }
}

// Scatter the pressure increment and clear stale ghost/reference values first.
void ProjectionSolver::load_pressure_correction(const Eigen::VectorXd& values) {
  fields_.pressure_correction.setZero();
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.pressure_correction(i, j) =
          values(static_cast<int>(grid_.pressure_index(i, j)));
    }
  }
}

// Solve the projection u predictor with pseudo-time inertia.
LinearSolveResult ProjectionSolver::solve_u_predictor() {
  const MomentumAssembly assembly = assemble_u_projection_momentum(grid_, config_, fields_);
  const LinearSolveResult result =
      solve_linear_system(assembly.system, SolverKind::kBiCGSTAB, 1e-8, 10000);
  fields_.u_star = fields_.u;
  load_u_solution(result.solution);
  return result;
}

// Solve the projection v predictor with pseudo-time inertia.
LinearSolveResult ProjectionSolver::solve_v_predictor() {
  const MomentumAssembly assembly = assemble_v_projection_momentum(grid_, config_, fields_);
  const LinearSolveResult result =
      solve_linear_system(assembly.system, SolverKind::kBiCGSTAB, 1e-8, 10000);
  fields_.v_star = fields_.v;
  load_v_solution(result.solution);
  return result;
}

// Solve the projection pressure increment from predicted mass imbalance.
LinearSolveResult ProjectionSolver::solve_pressure_increment_step() {
  const PressureCorrectionAssembly assembly =
      assemble_pressure_correction(grid_, config_, fields_);
  const LinearSolveResult result =
      solve_linear_system(assembly.system, SolverKind::kConjugateGradient, 1e-10, 10000);
  load_pressure_correction(result.solution);
  apply_boundary_conditions();
  return result;
}

// Correct the predicted velocity with the pressure increment and update pressure.
void ProjectionSolver::project_pressure_and_velocity() {
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
      fields_.pressure(i, j) += fields_.pressure_correction(i, j);
    }
  }

  const double reference_pressure = fields_.pressure(1, 1);
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.pressure(i, j) -= reference_pressure;
    }
  }
}

// Measure the maximum projection update across both velocity components.
double ProjectionSolver::compute_max_velocity_correction(
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

// Apply projection, refresh boundary-consistent faces, and package residual metrics.
IterationMetrics ProjectionSolver::project_fields_and_collect_metrics(
    int iteration,
    const LinearSolveResult& u_result,
    const LinearSolveResult& v_result,
    const LinearSolveResult& p_result,
    const Eigen::MatrixXd& u_before,
    const Eigen::MatrixXd& v_before) {
  project_pressure_and_velocity();
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

// Check the pseudo-time steady-state gates.
bool ProjectionSolver::has_converged(int iteration, const IterationMetrics& metrics) const {
  const bool momentum_ok =
      metrics.u_momentum_residual < config_.controls.momentum_tolerance &&
      metrics.v_momentum_residual < config_.controls.momentum_tolerance;
  const bool continuity_ok =
      metrics.continuity_residual < config_.controls.continuity_tolerance;
  const bool velocity_ok =
      metrics.max_velocity_correction < config_.controls.velocity_tolerance;
  return iteration >= config_.controls.min_iterations && momentum_ok && continuity_ok && velocity_ok;
}

// Execute projection iterations until the steady-state gates pass or the limit is reached.
SolveSummary ProjectionSolver::run() {
  SolveSummary summary;
  summary.solver = FlowSolverKind::kProjection;

  for (int iteration = 1; iteration <= config_.controls.max_iterations; ++iteration) {
    prepare_iteration();

    const Eigen::MatrixXd u_before = fields_.u;
    const Eigen::MatrixXd v_before = fields_.v;

    const LinearSolveResult u_result = solve_u_predictor();
    apply_boundary_conditions();
    refresh_face_velocities(fields_.u_star, fields_.v);
    const LinearSolveResult v_result = solve_v_predictor();

    apply_boundary_conditions();
    set_projection_correction_factors();
    const LinearSolveResult p_result = solve_pressure_increment_step();

    const IterationMetrics metrics = project_fields_and_collect_metrics(
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
