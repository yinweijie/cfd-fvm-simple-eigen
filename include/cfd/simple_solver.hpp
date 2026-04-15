#pragma once

#include <vector>

#include "cfd/case.hpp"
#include "cfd/discretization.hpp"

namespace cfd {

// Convergence history for one SIMPLE iteration.
struct IterationMetrics {
  int iteration = 0;
  double continuity_residual = 0.0;
  double u_momentum_residual = 0.0;
  double v_momentum_residual = 0.0;
  double pressure_correction_residual = 0.0;
  double max_velocity_correction = 0.0;
};

// Solver-level result returned after the final SIMPLE iteration.
struct SolveSummary {
  bool converged = false;
  int iterations = 0;
  std::vector<IterationMetrics> residual_history;
};

// Steady lid-driven cavity solver using a collocated SIMPLE loop.
class SimpleSolver {
 public:
  // Initialize grid storage and boundary-consistent face velocities from the case setup.
  explicit SimpleSolver(CavityCase config);

  // Run SIMPLE iterations until convergence or the iteration limit.
  SolveSummary run();

  const StructuredGrid& grid() const noexcept { return grid_; }
  const CavityCase& config() const noexcept { return config_; }
  const FlowFields& fields() const noexcept { return fields_; }

 private:
  StructuredGrid grid_;
  CavityCase config_;
  FlowFields fields_;

  // Apply cavity boundary conditions to the corrected and predicted fields.
  void apply_boundary_conditions();
  // Refresh boundaries and face velocities at the start of one SIMPLE iteration.
  void prepare_iteration();
  // Rebuild Rhie-Chow face velocities from the supplied cell-centered velocity fields.
  void refresh_face_velocities(const Eigen::MatrixXd& u_cells, const Eigen::MatrixXd& v_cells);
  // Scatter the solved u predictor into the collocated u_star field.
  void load_u_solution(const Eigen::VectorXd& values);
  // Scatter the solved v predictor into the collocated v_star field.
  void load_v_solution(const Eigen::VectorXd& values);
  // Scatter the solved pressure-correction field into pressure_correction.
  void load_pressure_correction(const Eigen::VectorXd& values);
  // Assemble and solve the u-momentum predictor system.
  LinearSolveResult solve_u_predictor();
  // Assemble and solve the v-momentum predictor system.
  LinearSolveResult solve_v_predictor();
  // Assemble and solve the SIMPLE pressure-correction system.
  LinearSolveResult solve_pressure_correction_step();
  // Update d_u from the relaxed u-momentum diagonal.
  void update_u_correction_factors(const Eigen::VectorXd& diagonal);
  // Update d_v from the relaxed v-momentum diagonal.
  void update_v_correction_factors(const Eigen::VectorXd& diagonal);
  // Apply the pressure-correction update to the velocity and pressure fields.
  void correct_pressure_and_velocity();
  // Measure the largest corrected-cell velocity change in the current iteration.
  double compute_max_velocity_correction(
      const Eigen::MatrixXd& u_before,
      const Eigen::MatrixXd& v_before) const;
  // Finish one corrector stage and collect convergence diagnostics.
  IterationMetrics correct_fields_and_collect_metrics(
      int iteration,
      const LinearSolveResult& u_result,
      const LinearSolveResult& v_result,
      const LinearSolveResult& p_result,
      const Eigen::MatrixXd& u_before,
      const Eigen::MatrixXd& v_before);
  // Check whether the current iteration satisfies all convergence gates.
  bool has_converged(int iteration, const IterationMetrics& metrics) const;
};

}  // namespace cfd
