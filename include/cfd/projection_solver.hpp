#pragma once

#include "cfd/case.hpp"
#include "cfd/discretization.hpp"
#include "cfd/solver.hpp"

namespace cfd {

// Steady cavity solver using pseudo-transient incremental projection iterations.
class ProjectionSolver {
 public:
  // Initialize grid storage, projection scales, and boundary-consistent face velocities.
  explicit ProjectionSolver(CavityCase config);

  // Run projection iterations until the steady convergence gates pass or the iteration limit is hit.
  SolveSummary run();

  const StructuredGrid& grid() const noexcept { return grid_; }
  const CavityCase& config() const noexcept { return config_; }
  const FlowFields& fields() const noexcept { return fields_; }

 private:
  StructuredGrid grid_;
  CavityCase config_;
  FlowFields fields_;

  // Apply cavity boundary conditions to corrected and predicted fields.
  void apply_boundary_conditions();
  // Update the projection pressure-to-velocity response factors from the pseudo-time step.
  void set_projection_correction_factors();
  // Refresh boundaries, correction factors, and corrected-field face velocities.
  void prepare_iteration();
  // Rebuild Rhie-Chow face velocities from the supplied cell-centered velocity fields.
  void refresh_face_velocities(const Eigen::MatrixXd& u_cells, const Eigen::MatrixXd& v_cells);
  // Scatter the solved u predictor into the collocated u_star field.
  void load_u_solution(const Eigen::VectorXd& values);
  // Scatter the solved v predictor into the collocated v_star field.
  void load_v_solution(const Eigen::VectorXd& values);
  // Scatter the solved pressure increment into pressure_correction.
  void load_pressure_correction(const Eigen::VectorXd& values);
  // Assemble and solve the pseudo-transient u predictor.
  LinearSolveResult solve_u_predictor();
  // Assemble and solve the pseudo-transient v predictor.
  LinearSolveResult solve_v_predictor();
  // Assemble and solve the projection pressure-increment equation.
  LinearSolveResult solve_pressure_increment_step();
  // Project the predicted velocity field and update pressure.
  void project_pressure_and_velocity();
  // Measure the largest corrected-cell velocity change in the current projection iteration.
  double compute_max_velocity_correction(
      const Eigen::MatrixXd& u_before,
      const Eigen::MatrixXd& v_before) const;
  // Finish one projection step and collect convergence diagnostics.
  IterationMetrics project_fields_and_collect_metrics(
      int iteration,
      const LinearSolveResult& u_result,
      const LinearSolveResult& v_result,
      const LinearSolveResult& p_result,
      const Eigen::MatrixXd& u_before,
      const Eigen::MatrixXd& v_before);
  // Check whether the current pseudo-time state is steady enough.
  bool has_converged(int iteration, const IterationMetrics& metrics) const;
};

}  // namespace cfd
