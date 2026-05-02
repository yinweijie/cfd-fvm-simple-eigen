#pragma once

#include <vector>

namespace cfd {

// Solver algorithm selected for one cavity run.
enum class FlowSolverKind {
  kSimple,
  kProjection,
};

// Return the stable CLI/file spelling for a solver kind.
inline const char* flow_solver_kind_name(FlowSolverKind kind) {
  switch (kind) {
    case FlowSolverKind::kSimple:
      return "simple";
    case FlowSolverKind::kProjection:
      return "projection";
  }
  return "unknown";
}

// Convergence history for one nonlinear or pseudo-time solver iteration.
struct IterationMetrics {
  int iteration = 0;
  double continuity_residual = 0.0;
  double u_momentum_residual = 0.0;
  double v_momentum_residual = 0.0;
  double pressure_correction_residual = 0.0;
  double max_velocity_correction = 0.0;
};

// Solver-level result returned after the final iteration.
struct SolveSummary {
  FlowSolverKind solver = FlowSolverKind::kSimple;
  bool converged = false;
  int iterations = 0;
  std::vector<IterationMetrics> residual_history;
};

}  // namespace cfd
