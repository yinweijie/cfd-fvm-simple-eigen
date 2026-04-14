#include "cfd/simple_solver.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace cfd {

SimpleSolver::SimpleSolver(CavityCase config)
    : grid_(config.mesh_spec), config_(std::move(config)), fields_(grid_) {
  apply_boundary_conditions();
  update_face_velocities(
      grid_,
      config_,
      fields_.u,
      fields_.v,
      fields_.pressure,
      fields_.d_u,
      fields_.d_v,
      &fields_.u_face,
      &fields_.v_face);
}

void SimpleSolver::apply_boundary_conditions() {
  apply_cavity_boundary_conditions(grid_, config_, &fields_);
}

void SimpleSolver::load_u_solution(const Eigen::VectorXd& values) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.u_star(i, j) = values(static_cast<int>(grid_.u_index(i, j)));
    }
  }
}

void SimpleSolver::load_v_solution(const Eigen::VectorXd& values) {
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.v_star(i, j) = values(static_cast<int>(grid_.v_index(i, j)));
    }
  }
}

void SimpleSolver::load_pressure_correction(const Eigen::VectorXd& values) {
  fields_.pressure_correction.setZero();
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.pressure_correction(i, j) =
          values(static_cast<int>(grid_.pressure_index(i, j)));
    }
  }
}

void SimpleSolver::correct_pressure_and_velocity() {
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

  const double reference_pressure = fields_.pressure(1, 1);
  for (int j = 1; j <= grid_.ny(); ++j) {
    for (int i = 1; i <= grid_.nx(); ++i) {
      fields_.pressure(i, j) -= reference_pressure;
    }
  }
}

SolveSummary SimpleSolver::run() {
  SolveSummary summary;

  for (int iteration = 1; iteration <= config_.controls.max_iterations; ++iteration) {
    apply_boundary_conditions();
    update_face_velocities(
        grid_,
        config_,
        fields_.u,
        fields_.v,
        fields_.pressure,
        fields_.d_u,
        fields_.d_v,
        &fields_.u_face,
        &fields_.v_face);

    const MomentumAssembly u_assembly = assemble_u_momentum(grid_, config_, fields_);
    const LinearSolveResult u_result = solve_linear_system(
        u_assembly.system, SolverKind::kBiCGSTAB, 1e-8, 10000);
    fields_.u_star = fields_.u;
    load_u_solution(u_result.solution);
    for (int j = 1; j <= grid_.ny(); ++j) {
      for (int i = 1; i <= grid_.nx(); ++i) {
        fields_.d_u(i, j) =
            (grid_.dx() * grid_.dy()) / u_assembly.diagonal(static_cast<int>(grid_.u_index(i, j)));
      }
    }

    update_face_velocities(
        grid_,
        config_,
        fields_.u_star,
        fields_.v,
        fields_.pressure,
        fields_.d_u,
        fields_.d_v,
        &fields_.u_face,
        &fields_.v_face);

    const MomentumAssembly v_assembly = assemble_v_momentum(grid_, config_, fields_);
    const LinearSolveResult v_result = solve_linear_system(
        v_assembly.system, SolverKind::kBiCGSTAB, 1e-8, 10000);
    fields_.v_star = fields_.v;
    load_v_solution(v_result.solution);
    for (int j = 1; j <= grid_.ny(); ++j) {
      for (int i = 1; i <= grid_.nx(); ++i) {
        fields_.d_v(i, j) =
            (grid_.dx() * grid_.dy()) / v_assembly.diagonal(static_cast<int>(grid_.v_index(i, j)));
      }
    }

    apply_boundary_conditions();

    const PressureCorrectionAssembly p_assembly =
        assemble_pressure_correction(grid_, config_, fields_);
    const LinearSolveResult p_result = solve_linear_system(
        p_assembly.system, SolverKind::kConjugateGradient, 1e-10, 10000);
    load_pressure_correction(p_result.solution);

    const Eigen::MatrixXd u_before = fields_.u;
    const Eigen::MatrixXd v_before = fields_.v;

    correct_pressure_and_velocity();
    apply_boundary_conditions();
    update_face_velocities(
        grid_,
        config_,
        fields_.u,
        fields_.v,
        fields_.pressure,
        fields_.d_u,
        fields_.d_v,
        &fields_.u_face,
        &fields_.v_face);

    double max_velocity_correction = 0.0;
    for (int j = 1; j <= grid_.ny(); ++j) {
      for (int i = 1; i <= grid_.nx(); ++i) {
        max_velocity_correction = std::max(
            max_velocity_correction,
            std::max(std::abs(fields_.u(i, j) - u_before(i, j)),
                     std::abs(fields_.v(i, j) - v_before(i, j))));
      }
    }

    IterationMetrics metrics;
    metrics.iteration = iteration;
    metrics.continuity_residual = compute_continuity_residual(grid_, config_, fields_);
    metrics.u_momentum_residual = u_result.relative_residual;
    metrics.v_momentum_residual = v_result.relative_residual;
    metrics.pressure_correction_residual = p_result.relative_residual;
    metrics.max_velocity_correction = max_velocity_correction;
    summary.residual_history.push_back(metrics);

    const bool momentum_ok =
        metrics.u_momentum_residual < config_.controls.momentum_tolerance &&
        metrics.v_momentum_residual < config_.controls.momentum_tolerance;
    const bool continuity_ok =
        metrics.continuity_residual < config_.controls.continuity_tolerance;
    const bool velocity_ok = metrics.max_velocity_correction < 1e-8;

    if (iteration >= config_.controls.min_iterations && momentum_ok && continuity_ok && velocity_ok) {
      summary.converged = true;
      summary.iterations = iteration;
      return summary;
    }
  }

  summary.iterations = static_cast<int>(summary.residual_history.size());
  return summary;
}

}  // namespace cfd
