#pragma once

#include <string>

#include "cfd/mesh.hpp"

namespace cfd {

// SIMPLE stopping criteria and under-relaxation settings for one solve.
struct SolverControls {
  int max_iterations = 2000;
  int min_iterations = 50;
  double continuity_tolerance = 1e-6;
  double momentum_tolerance = 1e-6;
  double alpha_u = 0.5;
  double alpha_v = 0.5;
  double alpha_p = 0.3;
};

// Full configuration for the steady lid-driven cavity benchmark case.
struct CavityCase {
  MeshSpec mesh_spec{32, 32, 1.0, 1.0};
  double density = 1.0;
  double reynolds = 100.0;
  double lid_velocity = 1.0;
  SolverControls controls{};
  std::string output_dir = "results/cavity_re100";

  // Derive the dynamic viscosity from density, lid speed, cavity length, and Reynolds number.
  double viscosity() const {
    return density * lid_velocity * mesh_spec.lx / reynolds;
  }
};

}  // namespace cfd
