#include "cfd/simple_solver.hpp"

#include <cassert>
#include <cmath>

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

  assert(summary.iterations > 1);
  assert(!summary.residual_history.empty());
  const double first_continuity = summary.residual_history.front().continuity_residual;
  const double last_continuity = summary.residual_history.back().continuity_residual;
  assert(std::isfinite(first_continuity));
  assert(std::isfinite(last_continuity));
  assert(last_continuity < first_continuity);

  const cfd::FlowFields& fields = solver.fields();
  for (int j = 1; j <= solver.grid().ny(); ++j) {
    for (int i = 1; i <= solver.grid().nx(); ++i) {
      assert(std::isfinite(fields.pressure(i, j)));
      assert(std::isfinite(fields.u(i, j)));
      assert(std::isfinite(fields.v(i, j)));
    }
  }

  return 0;
}
