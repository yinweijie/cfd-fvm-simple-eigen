#include "cfd/simple_solver.hpp"

#include <cassert>
#include <cmath>

namespace {

void expect_near(double actual, double expected, double tolerance = 1e-10) {
  assert(std::abs(actual - expected) <= tolerance);
}

}  // namespace

int main() {
  cfd::CavityCase config;
  config.mesh_spec = {8, 8, 1.0, 1.0};
  config.controls.max_iterations = 200;
  config.controls.min_iterations = 10;
  config.controls.alpha_u = 0.5;
  config.controls.alpha_v = 0.5;
  config.controls.alpha_p = 0.3;
  config.output_dir = "build/test_simple_smoke_output";

  cfd::SimpleSolver solver(config);
  const cfd::SolveSummary summary = solver.run();

  assert(summary.converged);
  assert(summary.iterations > 1);
  assert(!summary.residual_history.empty());
  const double first_continuity = summary.residual_history.front().continuity_residual;
  const double last_continuity = summary.residual_history.back().continuity_residual;
  assert(std::isfinite(first_continuity));
  assert(std::isfinite(last_continuity));
  assert(last_continuity < first_continuity);

  const cfd::IterationMetrics& final_metrics = summary.residual_history.back();
  assert(std::isfinite(final_metrics.u_momentum_residual));
  assert(std::isfinite(final_metrics.v_momentum_residual));
  assert(std::isfinite(final_metrics.pressure_correction_residual));
  assert(std::isfinite(final_metrics.max_velocity_correction));

  const cfd::FlowFields& fields = solver.fields();
  expect_near(fields.pressure(1, 1), 0.0);
  for (int j = 1; j <= solver.grid().ny(); ++j) {
    expect_near(fields.u(0, j), -fields.u(1, j));
    expect_near(fields.u(solver.grid().nx() + 1, j), -fields.u(solver.grid().nx(), j));
    expect_near(fields.v(0, j), -fields.v(1, j));
    expect_near(fields.v(solver.grid().nx() + 1, j), -fields.v(solver.grid().nx(), j));
    expect_near(fields.u_face(0, j), 0.0);
    expect_near(fields.u_face(solver.grid().nx(), j), 0.0);

    for (int i = 1; i <= solver.grid().nx(); ++i) {
      assert(std::isfinite(fields.pressure(i, j)));
      assert(std::isfinite(fields.u(i, j)));
      assert(std::isfinite(fields.v(i, j)));
    }
  }
  for (int i = 1; i <= solver.grid().nx(); ++i) {
    expect_near(fields.u(i, 0), -fields.u(i, 1));
    expect_near(fields.u(i, solver.grid().ny() + 1), 2.0 * config.lid_velocity - fields.u(i, solver.grid().ny()));
    expect_near(fields.v(i, 0), -fields.v(i, 1));
    expect_near(fields.v(i, solver.grid().ny() + 1), -fields.v(i, solver.grid().ny()));
    expect_near(fields.v_face(i, 0), 0.0);
    expect_near(fields.v_face(i, solver.grid().ny()), 0.0);
  }

  return 0;
}
