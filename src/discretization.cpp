#include "cfd/discretization.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace cfd {
namespace {

double positive_part(double value) {
  return value > 0.0 ? value : 0.0;
}

double grad_p_x(const StructuredGrid& grid, const Eigen::MatrixXd& pressure, int i, int j) {
  return (pressure(i + 1, j) - pressure(i - 1, j)) / (2.0 * grid.dx());
}

double grad_p_y(const StructuredGrid& grid, const Eigen::MatrixXd& pressure, int i, int j) {
  return (pressure(i, j + 1) - pressure(i, j - 1)) / (2.0 * grid.dy());
}

double pressure_force_u(const FlowFields& fields, int i, int j, double area) {
  return (fields.pressure(i - 1, j) - fields.pressure(i + 1, j)) * 0.5 * area;
}

double pressure_force_v(const FlowFields& fields, int i, int j, double area) {
  return (fields.pressure(i, j - 1) - fields.pressure(i, j + 1)) * 0.5 * area;
}

}  // namespace

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

void apply_cavity_boundary_conditions(
    const StructuredGrid& grid,
    const CavityCase& config,
    FlowFields* fields) {
  const int nx = grid.nx();
  const int ny = grid.ny();

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

  for (int j = 1; j <= ny; ++j) {
    (*u_face)(0, j) = 0.0;
    (*u_face)(nx, j) = 0.0;
    for (int i = 1; i < nx; ++i) {
      const double interp = 0.5 * (u_cells(i, j) + u_cells(i + 1, j));
      const double dp_face = (pressure(i + 1, j) - pressure(i, j)) / grid.dx();
      const double grad_avg =
          0.5 * (grad_p_x(grid, pressure, i, j) + grad_p_x(grid, pressure, i + 1, j));
      const double d_face = 0.5 * (d_u(i, j) + d_u(i + 1, j));
      (*u_face)(i, j) = interp - d_face * (dp_face - grad_avg);
    }
  }

  for (int i = 1; i <= nx; ++i) {
    (*v_face)(i, 0) = 0.0;
    (*v_face)(i, ny) = 0.0;
    for (int j = 1; j < ny; ++j) {
      const double interp = 0.5 * (v_cells(i, j) + v_cells(i, j + 1));
      const double dp_face = (pressure(i, j + 1) - pressure(i, j)) / grid.dy();
      const double grad_avg =
          0.5 * (grad_p_y(grid, pressure, i, j) + grad_p_y(grid, pressure, i, j + 1));
      const double d_face = 0.5 * (d_v(i, j) + d_v(i, j + 1));
      (*v_face)(i, j) = interp - d_face * (dp_face - grad_avg);
    }
  }
}

MomentumAssembly assemble_u_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields) {
  const int nx = grid.nx();
  const int ny = grid.ny();
  const double rho = config.density;
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
      const int row = static_cast<int>(grid.u_index(i, j));

      const double fe = rho * fields.u_face(i, j) * dy;
      const double fw = rho * fields.u_face(i - 1, j) * dy;
      const double fn = rho * fields.v_face(i, j) * dx;
      const double fs = rho * fields.v_face(i, j - 1) * dx;

      const double a_e = diff_x + positive_part(-fe);
      const double a_w = diff_x + positive_part(fw);
      const double a_n = diff_y + positive_part(-fn);
      const double a_s = diff_y + positive_part(fs);
      const double a_p_base = a_e + a_w + a_n + a_s + (fe - fw + fn - fs);
      double a_p = a_p_base / alpha;

      assembly.system.rhs(row) =
          pressure_force_u(fields, i, j, dy) +
          ((1.0 - alpha) / alpha) * a_p_base * fields.u(i, j);

      if (i < nx) {
        triplets.emplace_back(row, static_cast<int>(grid.u_index(i + 1, j)), -a_e);
      } else {
        a_p += a_e;
      }
      if (i > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.u_index(i - 1, j)), -a_w);
      } else {
        a_p += a_w;
      }
      if (j < ny) {
        triplets.emplace_back(row, static_cast<int>(grid.u_index(i, j + 1)), -a_n);
      } else {
        a_p += a_n;
        assembly.system.rhs(row) += 2.0 * a_n * config.lid_velocity;
      }
      if (j > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.u_index(i, j - 1)), -a_s);
      } else {
        a_p += a_s;
      }

      assembly.diagonal(row) = a_p;
      triplets.emplace_back(row, row, a_p);
    }
  }

  assembly.system.matrix.setFromTriplets(triplets.begin(), triplets.end());
  return assembly;
}

MomentumAssembly assemble_v_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields) {
  const int nx = grid.nx();
  const int ny = grid.ny();
  const double rho = config.density;
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
      const int row = static_cast<int>(grid.v_index(i, j));

      const double fe = rho * fields.u_face(i, j) * dy;
      const double fw = rho * fields.u_face(i - 1, j) * dy;
      const double fn = rho * fields.v_face(i, j) * dx;
      const double fs = rho * fields.v_face(i, j - 1) * dx;

      const double a_e = diff_x + positive_part(-fe);
      const double a_w = diff_x + positive_part(fw);
      const double a_n = diff_y + positive_part(-fn);
      const double a_s = diff_y + positive_part(fs);
      const double a_p_base = a_e + a_w + a_n + a_s + (fe - fw + fn - fs);
      double a_p = a_p_base / alpha;

      assembly.system.rhs(row) =
          pressure_force_v(fields, i, j, dx) +
          ((1.0 - alpha) / alpha) * a_p_base * fields.v(i, j);

      if (i < nx) {
        triplets.emplace_back(row, static_cast<int>(grid.v_index(i + 1, j)), -a_e);
      } else {
        a_p += a_e;
      }
      if (i > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.v_index(i - 1, j)), -a_w);
      } else {
        a_p += a_w;
      }
      if (j < ny) {
        triplets.emplace_back(row, static_cast<int>(grid.v_index(i, j + 1)), -a_n);
      } else {
        a_p += a_n;
      }
      if (j > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.v_index(i, j - 1)), -a_s);
      } else {
        a_p += a_s;
      }

      assembly.diagonal(row) = a_p;
      triplets.emplace_back(row, row, a_p);
    }
  }

  assembly.system.matrix.setFromTriplets(triplets.begin(), triplets.end());
  return assembly;
}

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

  for (int j = 1; j <= ny; ++j) {
    for (int i = 1; i <= nx; ++i) {
      const int row = static_cast<int>(grid.pressure_index(i, j));

      if (i == 1 && j == 1) {
        triplets.emplace_back(row, row, 1.0);
        continue;
      }

      const double a_e =
          (i < nx) ? rho * grid.dy() * 0.5 * (fields.d_u(i, j) + fields.d_u(i + 1, j)) / grid.dx() : 0.0;
      const double a_w =
          (i > 1) ? rho * grid.dy() * 0.5 * (fields.d_u(i, j) + fields.d_u(i - 1, j)) / grid.dx() : 0.0;
      const double a_n =
          (j < ny) ? rho * grid.dx() * 0.5 * (fields.d_v(i, j) + fields.d_v(i, j + 1)) / grid.dy() : 0.0;
      const double a_s =
          (j > 1) ? rho * grid.dx() * 0.5 * (fields.d_v(i, j) + fields.d_v(i, j - 1)) / grid.dy() : 0.0;

      const double imbalance =
          rho * grid.dy() * (u_face_star(i - 1, j) - u_face_star(i, j)) +
          rho * grid.dx() * (v_face_star(i, j - 1) - v_face_star(i, j));

      const double a_p = a_e + a_w + a_n + a_s;
      assembly.mass_imbalance(row) = imbalance;
      assembly.system.rhs(row) = imbalance;

      triplets.emplace_back(row, row, a_p);
      if (i < nx) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i + 1, j)), -a_e);
      }
      if (i > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i - 1, j)), -a_w);
      }
      if (j < ny) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i, j + 1)), -a_n);
      }
      if (j > 1) {
        triplets.emplace_back(row, static_cast<int>(grid.pressure_index(i, j - 1)), -a_s);
      }
    }
  }

  assembly.system.matrix.setFromTriplets(triplets.begin(), triplets.end());
  return assembly;
}

double compute_continuity_residual(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields) {
  Eigen::MatrixXd u_face = fields.u_face;
  Eigen::MatrixXd v_face = fields.v_face;
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
      const double imbalance =
          config.density * grid.dy() * (u_face(i - 1, j) - u_face(i, j)) +
          config.density * grid.dx() * (v_face(i, j - 1) - v_face(i, j));
      accum += std::abs(imbalance);
    }
  }
  return accum / static_cast<double>(grid.nx() * grid.ny());
}

}  // namespace cfd
