#include "cfd/discretization.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace cfd {
namespace {

// Cell-face mass fluxes surrounding one control volume.
struct FaceMassFluxes {
  double fe;
  double fw;
  double fn;
  double fs;
};

// Upwind-diffusion stencil coefficients for one momentum row.
struct MomentumCoefficients {
  double a_e;
  double a_w;
  double a_n;
  double a_s;
  double a_p_base;
};

// Neighbor coefficients for one pressure-correction row.
struct PressureCorrectionCoefficients {
  double a_e;
  double a_w;
  double a_n;
  double a_s;

  // Return the central pressure-correction diagonal assembled from the neighbors.
  double a_p() const {
    return a_e + a_w + a_n + a_s;
  }
};

// Return max(value, 0) for the upwind coefficient split.
double positive_part(double value) {
  return value > 0.0 ? value : 0.0;
}

// Estimate the cell-centered x pressure gradient with a central difference.
double pressure_gradient_x(const StructuredGrid& grid, const Eigen::MatrixXd& pressure, int i, int j) {
  return (pressure(i + 1, j) - pressure(i - 1, j)) / (2.0 * grid.dx());
}

// Estimate the cell-centered y pressure gradient with a central difference.
double pressure_gradient_y(const StructuredGrid& grid, const Eigen::MatrixXd& pressure, int i, int j) {
  return (pressure(i, j + 1) - pressure(i, j - 1)) / (2.0 * grid.dy());
}

// Convert the neighboring pressure difference into the u-equation source term.
double pressure_force_u(const FlowFields& fields, int i, int j, double area) {
  return (fields.pressure(i - 1, j) - fields.pressure(i + 1, j)) * 0.5 * area;
}

// Convert the neighboring pressure difference into the v-equation source term.
double pressure_force_v(const FlowFields& fields, int i, int j, double area) {
  return (fields.pressure(i, j - 1) - fields.pressure(i, j + 1)) * 0.5 * area;
}

// Gather the four face mass fluxes needed to assemble one momentum control volume.
FaceMassFluxes compute_face_mass_fluxes(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields,
    int i,
    int j) {
  return FaceMassFluxes{
      config.density * fields.u_face(i, j) * grid.dy(),
      config.density * fields.u_face(i - 1, j) * grid.dy(),
      config.density * fields.v_face(i, j) * grid.dx(),
      config.density * fields.v_face(i, j - 1) * grid.dx(),
  };
}

// Build the relaxed momentum stencil from the local face fluxes and diffusion terms.
MomentumCoefficients compute_momentum_coefficients(
    const FaceMassFluxes& fluxes,
    double diff_x,
    double diff_y) {
  const double a_e = diff_x + positive_part(-fluxes.fe);
  const double a_w = diff_x + positive_part(fluxes.fw);
  const double a_n = diff_y + positive_part(-fluxes.fn);
  const double a_s = diff_y + positive_part(fluxes.fs);
  return MomentumCoefficients{
      a_e,
      a_w,
      a_n,
      a_s,
      a_e + a_w + a_n + a_s + (fluxes.fe - fluxes.fw + fluxes.fn - fluxes.fs),
  };
}

// Reconstruct one vertical-face velocity with the Rhie-Chow correction.
double rhie_chow_u_face_velocity(
    const StructuredGrid& grid,
    const Eigen::MatrixXd& u_cells,
    const Eigen::MatrixXd& pressure,
    const Eigen::MatrixXd& d_u,
    int i,
    int j) {
  const double interp = 0.5 * (u_cells(i, j) + u_cells(i + 1, j));
  const double dp_face = (pressure(i + 1, j) - pressure(i, j)) / grid.dx();
  const double grad_avg = 0.5 * (pressure_gradient_x(grid, pressure, i, j) +
                                 pressure_gradient_x(grid, pressure, i + 1, j));
  const double d_face = 0.5 * (d_u(i, j) + d_u(i + 1, j));
  return interp - d_face * (dp_face - grad_avg);
}

// Reconstruct one horizontal-face velocity with the Rhie-Chow correction.
double rhie_chow_v_face_velocity(
    const StructuredGrid& grid,
    const Eigen::MatrixXd& v_cells,
    const Eigen::MatrixXd& pressure,
    const Eigen::MatrixXd& d_v,
    int i,
    int j) {
  const double interp = 0.5 * (v_cells(i, j) + v_cells(i, j + 1));
  const double dp_face = (pressure(i, j + 1) - pressure(i, j)) / grid.dy();
  const double grad_avg = 0.5 * (pressure_gradient_y(grid, pressure, i, j) +
                                 pressure_gradient_y(grid, pressure, i, j + 1));
  const double d_face = 0.5 * (d_v(i, j) + d_v(i, j + 1));
  return interp - d_face * (dp_face - grad_avg);
}

// Build the pressure-correction stencil coefficients for one interior cell.
PressureCorrectionCoefficients compute_pressure_correction_coefficients(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields,
    int i,
    int j) {
  const double rho = config.density;
  return PressureCorrectionCoefficients{
      (i < grid.nx()) ? rho * grid.dy() * 0.5 * (fields.d_u(i, j) + fields.d_u(i + 1, j)) / grid.dx()
                      : 0.0,
      (i > 1) ? rho * grid.dy() * 0.5 * (fields.d_u(i, j) + fields.d_u(i - 1, j)) / grid.dx() : 0.0,
      (j < grid.ny()) ? rho * grid.dx() * 0.5 * (fields.d_v(i, j) + fields.d_v(i, j + 1)) / grid.dy()
                      : 0.0,
      (j > 1) ? rho * grid.dx() * 0.5 * (fields.d_v(i, j) + fields.d_v(i, j - 1)) / grid.dy() : 0.0,
  };
}

// Compute the mass defect of one pressure control volume from its face fluxes.
double compute_mass_imbalance(
    const StructuredGrid& grid,
    double density,
    const Eigen::MatrixXd& u_face,
    const Eigen::MatrixXd& v_face,
    int i,
    int j) {
  return density * grid.dy() * (u_face(i - 1, j) - u_face(i, j)) +
         density * grid.dx() * (v_face(i, j - 1) - v_face(i, j));
}

// Either connect an interior neighbor into the sparse row or fold the missing
// boundary-side coefficient back into the current diagonal. Some boundaries also
// contribute an explicit source term, such as the moving lid in the u equation.
template <typename NeighborColumnFn>
void couple_neighbor_or_fold_boundary(
    std::vector<Eigen::Triplet<double>>* triplets,
    int row,
    bool has_neighbor,
    NeighborColumnFn&& neighbor_column,
    double neighbor_coefficient,
    double* diagonal,
    double* rhs = nullptr,
    double boundary_source = 0.0) {
  if (has_neighbor) {
    triplets->emplace_back(row, neighbor_column(), -neighbor_coefficient);
    return;
  }

  *diagonal += neighbor_coefficient;
  if (rhs != nullptr) {
    *rhs += boundary_source;
  }
}

}  // namespace

// Allocate all collocated cell and face fields for the requested grid.
FlowFields::FlowFields(const StructuredGrid& grid)
    : pressure(grid.nx() + 2, grid.ny() + 2),
      pressure_correction(grid.nx() + 2, grid.ny() + 2),
      u(grid.nx() + 2, grid.ny() + 2),
      v(grid.nx() + 2, grid.ny() + 2),
      u_star(grid.nx() + 2, grid.ny() + 2),
      v_star(grid.nx() + 2, grid.ny() + 2),
      d_u(grid.nx() + 2, grid.ny() + 2),
      d_v(grid.nx() + 2, grid.ny() + 2),
      u_face(grid.nx() + 1, grid.ny() + 2),
      v_face(grid.nx() + 2, grid.ny() + 1) {
  reset();
}

// Reset every stored field and correction factor to zero.
void FlowFields::reset() {
  pressure.setZero();
  pressure_correction.setZero();
  u.setZero();
  v.setZero();
  u_star.setZero();
  v_star.setZero();
  d_u.setZero();
  d_v.setZero();
  u_face.setZero();
  v_face.setZero();
}

// Update the ghost layer so every SIMPLE stage sees the cavity wall conditions.
void apply_cavity_boundary_conditions(
    const StructuredGrid& grid,
    const CavityCase& config,
    FlowFields* fields) {
  const int nx = grid.nx();
  const int ny = grid.ny();

  // Mirror the side-wall cell values into one ghost layer to enforce no-slip walls.
  for (int j = 1; j <= ny; ++j) {
    fields->u(0, j) = -fields->u(1, j);
    fields->u(nx + 1, j) = -fields->u(nx, j);
    fields->u_star(0, j) = -fields->u_star(1, j);
    fields->u_star(nx + 1, j) = -fields->u_star(nx, j);

    fields->v(0, j) = -fields->v(1, j);
    fields->v(nx + 1, j) = -fields->v(nx, j);
    fields->v_star(0, j) = -fields->v_star(1, j);
    fields->v_star(nx + 1, j) = -fields->v_star(nx, j);
  }

  // Bottom wall is stationary; the top wall uses the moving-lid ghost relation.
  for (int i = 1; i <= nx; ++i) {
    fields->u(i, 0) = -fields->u(i, 1);
    fields->u(i, ny + 1) = 2.0 * config.lid_velocity - fields->u(i, ny);
    fields->u_star(i, 0) = -fields->u_star(i, 1);
    fields->u_star(i, ny + 1) = 2.0 * config.lid_velocity - fields->u_star(i, ny);

    fields->v(i, 0) = -fields->v(i, 1);
    fields->v(i, ny + 1) = -fields->v(i, ny);
    fields->v_star(i, 0) = -fields->v_star(i, 1);
    fields->v_star(i, ny + 1) = -fields->v_star(i, ny);
  }

  // Pressure and pressure_correction keep zero normal gradient at all cavity walls.
  for (int j = 1; j <= ny; ++j) {
    fields->pressure(0, j) = fields->pressure(1, j);
    fields->pressure(nx + 1, j) = fields->pressure(nx, j);
    fields->pressure_correction(0, j) = fields->pressure_correction(1, j);
    fields->pressure_correction(nx + 1, j) = fields->pressure_correction(nx, j);
  }
  for (int i = 0; i <= nx + 1; ++i) {
    fields->pressure(i, 0) = fields->pressure(i, 1);
    fields->pressure(i, ny + 1) = fields->pressure(i, ny);
    fields->pressure_correction(i, 0) = fields->pressure_correction(i, 1);
    fields->pressure_correction(i, ny + 1) = fields->pressure_correction(i, ny);
  }
}

// Recompute physical-face velocities from cell values and pressure correction factors.
// The function first clears the face arrays, then sweeps the vertical and horizontal
// face strips on the structured grid and reconstructs each interior face from the
// two adjacent cells plus the Rhie-Chow pressure correction.
void update_face_velocities(
    const StructuredGrid& grid,
    const CavityCase&,
    const Eigen::MatrixXd& u_cells,
    const Eigen::MatrixXd& v_cells,
    const Eigen::MatrixXd& pressure,
    const Eigen::MatrixXd& d_u,
    const Eigen::MatrixXd& d_v,
    Eigen::MatrixXd* u_face,
    Eigen::MatrixXd* v_face) {
  const int nx = grid.nx();
  const int ny = grid.ny();

  u_face->setZero();
  v_face->setZero();

  // Loop over each interior row of vertical faces. Physical boundary faces stay at zero
  // for the cavity walls, and every interior face is reconstructed from its neighboring cells.
  for (int j = 1; j <= ny; ++j) {
    (*u_face)(0, j) = 0.0;
    (*u_face)(nx, j) = 0.0;
    for (int i = 1; i < nx; ++i) {
      (*u_face)(i, j) = rhie_chow_u_face_velocity(grid, u_cells, pressure, d_u, i, j);
    }
  }

  // Loop over each interior column of horizontal faces. As above, the wall faces are
  // clamped to zero normal velocity and only the interior faces use Rhie-Chow reconstruction.
  for (int i = 1; i <= nx; ++i) {
    (*v_face)(i, 0) = 0.0;
    (*v_face)(i, ny) = 0.0;
    for (int j = 1; j < ny; ++j) {
      (*v_face)(i, j) = rhie_chow_v_face_velocity(grid, v_cells, pressure, d_v, i, j);
    }
  }
}

// Assemble the relaxed u-momentum system and expose its diagonal for SIMPLE corrections.
// The assembly loops over every interior control volume, computes local convection and
// diffusion coefficients from the current face fluxes, then writes the row into sparse
// triplets together with the under-relaxed diagonal and pressure/lid source terms.
MomentumAssembly assemble_u_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields) {
  const int nx = grid.nx();
  const int ny = grid.ny();
  const double mu = config.viscosity();
  const double dx = grid.dx();
  const double dy = grid.dy();
  const double diff_x = mu * dy / dx;
  const double diff_y = mu * dx / dy;
  const double alpha = config.controls.alpha_u;

  MomentumAssembly assembly;
  assembly.system.matrix.resize(static_cast<int>(grid.u_unknown_count()),
                                static_cast<int>(grid.u_unknown_count()));
  assembly.system.rhs = Eigen::VectorXd::Zero(static_cast<int>(grid.u_unknown_count()));
  assembly.diagonal = Eigen::VectorXd::Zero(static_cast<int>(grid.u_unknown_count()));

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(grid.u_unknown_count() * 5);

  for (int j = 1; j <= ny; ++j) {
    for (int i = 1; i <= nx; ++i) {
      // Map the current control volume to one row in the flattened linear system.
      const int row = static_cast<int>(grid.u_index(i, j));

      // Build the local stencil from the surrounding face fluxes and diffusion terms.
      const FaceMassFluxes fluxes = compute_face_mass_fluxes(grid, config, fields, i, j);
      const MomentumCoefficients coeffs =
          compute_momentum_coefficients(fluxes, diff_x, diff_y);
      double a_p = coeffs.a_p_base / alpha;

      // Start the RHS from the pressure source and the SIMPLE under-relaxation correction.
      assembly.system.rhs(row) =
          pressure_force_u(fields, i, j, dy) +
          ((1.0 - alpha) / alpha) * coeffs.a_p_base * fields.u(i, j);

      // For each stencil side, either connect the interior neighbor or fold the boundary-side
      // coefficient into the current diagonal. The north boundary also injects the moving-lid source.
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          i < nx,
          [&]() { return static_cast<int>(grid.u_index(i + 1, j)); },
          coeffs.a_e,
          &a_p);
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          i > 1,
          [&]() { return static_cast<int>(grid.u_index(i - 1, j)); },
          coeffs.a_w,
          &a_p);
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          j < ny,
          [&]() { return static_cast<int>(grid.u_index(i, j + 1)); },
          coeffs.a_n,
          &a_p,
          &assembly.system.rhs(row),
          2.0 * coeffs.a_n * config.lid_velocity);
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          j > 1,
          [&]() { return static_cast<int>(grid.u_index(i, j - 1)); },
          coeffs.a_s,
          &a_p);

      assembly.diagonal(row) = a_p;
      triplets.emplace_back(row, row, a_p);
    }
  }

  // Finalize the sparse matrix from the collected five-point-stencil triplets.
  assembly.system.matrix.setFromTriplets(triplets.begin(), triplets.end());
  return assembly;
}

// Assemble the relaxed v-momentum system and expose its diagonal for SIMPLE corrections.
// This mirrors the u-equation assembly, but uses the v pressure source and does not
// apply a moving-lid contribution on the north boundary.
MomentumAssembly assemble_v_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields) {
  const int nx = grid.nx();
  const int ny = grid.ny();
  const double mu = config.viscosity();
  const double dx = grid.dx();
  const double dy = grid.dy();
  const double diff_x = mu * dy / dx;
  const double diff_y = mu * dx / dy;
  const double alpha = config.controls.alpha_v;

  MomentumAssembly assembly;
  assembly.system.matrix.resize(static_cast<int>(grid.v_unknown_count()),
                                static_cast<int>(grid.v_unknown_count()));
  assembly.system.rhs = Eigen::VectorXd::Zero(static_cast<int>(grid.v_unknown_count()));
  assembly.diagonal = Eigen::VectorXd::Zero(static_cast<int>(grid.v_unknown_count()));

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(grid.v_unknown_count() * 5);

  for (int j = 1; j <= ny; ++j) {
    for (int i = 1; i <= nx; ++i) {
      // Map the current control volume to one row in the flattened linear system.
      const int row = static_cast<int>(grid.v_index(i, j));

      // Build the local stencil from the surrounding face fluxes and diffusion terms.
      const FaceMassFluxes fluxes = compute_face_mass_fluxes(grid, config, fields, i, j);
      const MomentumCoefficients coeffs =
          compute_momentum_coefficients(fluxes, diff_x, diff_y);
      double a_p = coeffs.a_p_base / alpha;

      // Start the RHS from the pressure source and the SIMPLE under-relaxation correction.
      assembly.system.rhs(row) =
          pressure_force_v(fields, i, j, dx) +
          ((1.0 - alpha) / alpha) * coeffs.a_p_base * fields.v(i, j);

      // As in the u equation, each side either couples to an interior neighbor or folds its
      // boundary-side coefficient into the diagonal when the control volume touches a wall.
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          i < nx,
          [&]() { return static_cast<int>(grid.v_index(i + 1, j)); },
          coeffs.a_e,
          &a_p);
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          i > 1,
          [&]() { return static_cast<int>(grid.v_index(i - 1, j)); },
          coeffs.a_w,
          &a_p);
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          j < ny,
          [&]() { return static_cast<int>(grid.v_index(i, j + 1)); },
          coeffs.a_n,
          &a_p);
      couple_neighbor_or_fold_boundary(
          &triplets,
          row,
          j > 1,
          [&]() { return static_cast<int>(grid.v_index(i, j - 1)); },
          coeffs.a_s,
          &a_p);

      assembly.diagonal(row) = a_p;
      triplets.emplace_back(row, row, a_p);
    }
  }

  // Finalize the sparse matrix from the collected five-point-stencil triplets.
  assembly.system.matrix.setFromTriplets(triplets.begin(), triplets.end());
  return assembly;
}

// Assemble the pressure-correction system from predicted face-flux imbalance.
// The function first rebuilds face velocities from the predicted fields, then loops over
// all interior pressure cells, computes the local continuity defect, and writes the
// resulting Laplacian-like stencil plus RHS into the sparse system.
PressureCorrectionAssembly assemble_pressure_correction(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields) {
  const int nx = grid.nx();
  const int ny = grid.ny();
  const double rho = config.density;

  PressureCorrectionAssembly assembly;
  assembly.system.matrix.resize(static_cast<int>(grid.pressure_cell_count()),
                                static_cast<int>(grid.pressure_cell_count()));
  assembly.system.rhs = Eigen::VectorXd::Zero(static_cast<int>(grid.pressure_cell_count()));
  assembly.mass_imbalance = Eigen::VectorXd::Zero(static_cast<int>(grid.pressure_cell_count()));

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(grid.pressure_cell_count() * 5);

  Eigen::MatrixXd u_face_star = fields.u_face;
  Eigen::MatrixXd v_face_star = fields.v_face;

  // Reconstruct predictor face fluxes from u_star and v_star before evaluating continuity errors.
  update_face_velocities(
      grid,
      config,
      fields.u_star,
      fields.v_star,
      fields.pressure,
      fields.d_u,
      fields.d_v,
      &u_face_star,
      &v_face_star);

  // Fix one pressure-correction degree of freedom because only pressure gradients matter.
  for (int j = 1; j <= ny; ++j) {
    for (int i = 1; i <= nx; ++i) {
      // Map the current pressure control volume to one row in the flattened system.
      const int row = static_cast<int>(grid.pressure_index(i, j));

      if (i == 1 && j == 1) {
        triplets.emplace_back(row, row, 1.0);
        continue;
      }

      // Build the pressure-correction stencil from neighboring d_u / d_v values and use the
      // predicted face velocities to measure how much mass this cell currently gains or loses.
      const PressureCorrectionCoefficients coeffs =
          compute_pressure_correction_coefficients(grid, config, fields, i, j);
      const double imbalance =
          compute_mass_imbalance(grid, rho, u_face_star, v_face_star, i, j);

      const double a_p = coeffs.a_p();
      assembly.mass_imbalance(row) = imbalance;
      assembly.system.rhs(row) = imbalance;

      triplets.emplace_back(row, row, a_p);
      if (i < nx) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i + 1, j)), -coeffs.a_e);
      }
      if (i > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i - 1, j)), -coeffs.a_w);
      }
      if (j < ny) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i, j + 1)), -coeffs.a_n);
      }
      if (j > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i, j - 1)), -coeffs.a_s);
      }
    }
  }

  // Finalize the sparse matrix from the collected pressure-correction triplets.
  assembly.system.matrix.setFromTriplets(triplets.begin(), triplets.end());
  return assembly;
}

// Average the absolute mass imbalance over all interior control volumes.
// The function rebuilds face velocities from the current corrected fields, loops over the
// whole interior pressure grid, and reports the mean absolute continuity defect.
double compute_continuity_residual(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields) {
  Eigen::MatrixXd u_face = fields.u_face;
  Eigen::MatrixXd v_face = fields.v_face;

  // Recompute face velocities from the corrected fields so the residual reflects the final state
  // of the current SIMPLE iteration rather than any partially updated predictor fields.
  update_face_velocities(
      grid,
      config,
      fields.u,
      fields.v,
      fields.pressure,
      fields.d_u,
      fields.d_v,
      &u_face,
      &v_face);

  double accum = 0.0;
  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 1; i <= grid.nx(); ++i) {
      accum += std::abs(
          compute_mass_imbalance(grid, config.density, u_face, v_face, i, j));
    }
  }

  return accum / static_cast<double>(grid.nx() * grid.ny());
}

}  // namespace cfd
