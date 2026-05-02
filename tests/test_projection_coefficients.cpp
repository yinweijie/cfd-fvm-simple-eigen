#include "cfd/case.hpp"
#include "cfd/discretization.hpp"

#include <cassert>
#include <cmath>

namespace {

void expect_near(double actual, double expected, double tolerance = 1e-12) {
  assert(std::abs(actual - expected) <= tolerance);
}

double sparse_coeff(
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix,
    int row,
    int col) {
  return matrix.coeff(row, col);
}

}  // namespace

int main() {
  cfd::CavityCase config;
  config.mesh_spec = {4, 4, 1.0, 1.0};
  config.controls.projection_dt = 0.1;

  cfd::StructuredGrid grid(config.mesh_spec);
  cfd::FlowFields fields(grid);
  cfd::apply_cavity_boundary_conditions(grid, config, &fields);
  cfd::update_face_velocities(
      grid,
      config,
      fields.u,
      fields.v,
      fields.pressure,
      fields.d_u,
      fields.d_v,
      &fields.u_face,
      &fields.v_face);

  const cfd::MomentumAssembly u_assembly =
      cfd::assemble_u_projection_momentum(grid, config, fields);
  const cfd::MomentumAssembly v_assembly =
      cfd::assemble_v_projection_momentum(grid, config, fields);

  const double diff_x = config.viscosity() * grid.dy() / grid.dx();
  const double diff_y = config.viscosity() * grid.dx() / grid.dy();
  const double transient = config.density * grid.dx() * grid.dy() / config.controls.projection_dt;
  const double top_left_diag = transient + 3.0 * diff_x + 3.0 * diff_y;

  const int u_top_left = static_cast<int>(grid.u_index(1, grid.ny()));
  const int u_top_left_east = static_cast<int>(grid.u_index(2, grid.ny()));
  const int u_top_left_south = static_cast<int>(grid.u_index(1, grid.ny() - 1));
  expect_near(u_assembly.diagonal(u_top_left), top_left_diag);
  expect_near(sparse_coeff(u_assembly.system.matrix, u_top_left, u_top_left), top_left_diag);
  expect_near(sparse_coeff(u_assembly.system.matrix, u_top_left, u_top_left_east), -diff_x);
  expect_near(sparse_coeff(u_assembly.system.matrix, u_top_left, u_top_left_south), -diff_y);
  expect_near(u_assembly.system.rhs(u_top_left), 2.0 * diff_y * config.lid_velocity);

  const int v_top_left = static_cast<int>(grid.v_index(1, grid.ny()));
  const int v_top_left_east = static_cast<int>(grid.v_index(2, grid.ny()));
  const int v_top_left_south = static_cast<int>(grid.v_index(1, grid.ny() - 1));
  expect_near(v_assembly.diagonal(v_top_left), top_left_diag);
  expect_near(sparse_coeff(v_assembly.system.matrix, v_top_left, v_top_left), top_left_diag);
  expect_near(sparse_coeff(v_assembly.system.matrix, v_top_left, v_top_left_east), -diff_x);
  expect_near(sparse_coeff(v_assembly.system.matrix, v_top_left, v_top_left_south), -diff_y);
  expect_near(v_assembly.system.rhs(v_top_left), 0.0);

  fields.reset();
  fields.d_u.setConstant(config.controls.projection_dt / config.density);
  fields.d_v.setConstant(config.controls.projection_dt / config.density);
  fields.u_star.setZero();
  fields.v_star.setZero();
  const cfd::PressureCorrectionAssembly p_assembly =
      cfd::assemble_pressure_correction(grid, config, fields);

  assert(p_assembly.system.matrix.rows() == static_cast<int>(grid.pressure_cell_count()));
  assert(p_assembly.system.matrix.cols() == static_cast<int>(grid.pressure_cell_count()));
  expect_near(p_assembly.system.matrix.coeff(0, 0), 1.0);
  expect_near(p_assembly.system.rhs(0), 0.0);
  expect_near(p_assembly.mass_imbalance(0), 0.0);

  const int center = static_cast<int>(grid.pressure_index(2, 2));
  const double pressure_coeff = config.controls.projection_dt;
  expect_near(sparse_coeff(p_assembly.system.matrix, center, center), 4.0 * pressure_coeff);
  expect_near(sparse_coeff(p_assembly.system.matrix, center, static_cast<int>(grid.pressure_index(3, 2))),
              -pressure_coeff);
  expect_near(sparse_coeff(p_assembly.system.matrix, center, static_cast<int>(grid.pressure_index(1, 2))),
              -pressure_coeff);
  expect_near(sparse_coeff(p_assembly.system.matrix, center, static_cast<int>(grid.pressure_index(2, 3))),
              -pressure_coeff);
  expect_near(sparse_coeff(p_assembly.system.matrix, center, static_cast<int>(grid.pressure_index(2, 1))),
              -pressure_coeff);
  for (int row = 0; row < p_assembly.system.rhs.size(); ++row) {
    expect_near(p_assembly.system.rhs(row), 0.0);
  }

  return 0;
}
