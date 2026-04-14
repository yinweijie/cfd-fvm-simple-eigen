#pragma once

#include <vector>

#include "cfd/case.hpp"
#include "cfd/discretization.hpp"

namespace cfd {

struct IterationMetrics {
  int iteration = 0;
  double continuity_residual = 0.0;
  double u_momentum_residual = 0.0;
  double v_momentum_residual = 0.0;
  double pressure_correction_residual = 0.0;
  double max_velocity_correction = 0.0;
};

struct SolveSummary {
  bool converged = false;
  int iterations = 0;
  std::vector<IterationMetrics> residual_history;
};

class SimpleSolver {
 public:
  explicit SimpleSolver(CavityCase config);

  SolveSummary run();

  const StructuredGrid& grid() const noexcept { return grid_; }
  const CavityCase& config() const noexcept { return config_; }
  const FlowFields& fields() const noexcept { return fields_; }

 private:
  StructuredGrid grid_;
  CavityCase config_;
  FlowFields fields_;

  void apply_boundary_conditions();
  void load_u_solution(const Eigen::VectorXd& values);
  void load_v_solution(const Eigen::VectorXd& values);
  void load_pressure_correction(const Eigen::VectorXd& values);
  void correct_pressure_and_velocity();
};

}  // namespace cfd
