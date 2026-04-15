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
  config.controls.alpha_u = 0.5;
  config.controls.alpha_v = 0.5;
  config.controls.alpha_p = 0.3;

  cfd::StructuredGrid grid(config.mesh_spec);
  cfd::FlowFields fields(grid);

  fields.u(1, 2) = 0.3;
  fields.u(grid.nx(), 3) = -0.4;
  fields.u_star(2, 1) = -0.2;
  fields.u_star(3, grid.ny()) = 0.15;
  fields.v(1, 3) = -0.7;
  fields.v(grid.nx(), 2) = 0.25;
  fields.v_star(2, 1) = 0.6;
  fields.v_star(3, grid.ny()) = -0.35;
  fields.pressure(1, 2) = 1.3;
  fields.pressure(grid.nx(), 3) = -0.9;
  fields.pressure_correction(2, 1) = 0.75;
  fields.pressure_correction(3, grid.ny()) = -0.45;

  cfd::apply_cavity_boundary_conditions(grid, config, &fields);

  expect_near(fields.u(0, 2), -fields.u(1, 2));
  expect_near(fields.u(grid.nx() + 1, 3), -fields.u(grid.nx(), 3));
  expect_near(fields.u_star(0, 1), -fields.u_star(1, 1));
  expect_near(fields.u_star(grid.nx() + 1, grid.ny()), -fields.u_star(grid.nx(), grid.ny()));

  expect_near(fields.v(0, 3), -fields.v(1, 3));
  expect_near(fields.v(grid.nx() + 1, 2), -fields.v(grid.nx(), 2));
  expect_near(fields.v_star(0, 1), -fields.v_star(1, 1));
  expect_near(fields.v_star(grid.nx() + 1, grid.ny()), -fields.v_star(grid.nx(), grid.ny()));

  expect_near(fields.u(2, 0), -fields.u(2, 1));
  expect_near(fields.u(3, grid.ny() + 1), 2.0 * config.lid_velocity - fields.u(3, grid.ny()));
  expect_near(fields.u_star(2, 0), -fields.u_star(2, 1));
  expect_near(
      fields.u_star(3, grid.ny() + 1), 2.0 * config.lid_velocity - fields.u_star(3, grid.ny()));

  expect_near(fields.v(2, 0), -fields.v(2, 1));
  expect_near(fields.v(3, grid.ny() + 1), -fields.v(3, grid.ny()));
  expect_near(fields.v_star(2, 0), -fields.v_star(2, 1));
  expect_near(fields.v_star(3, grid.ny() + 1), -fields.v_star(3, grid.ny()));

  expect_near(fields.pressure(0, 2), fields.pressure(1, 2));
  expect_near(fields.pressure(grid.nx() + 1, 3), fields.pressure(grid.nx(), 3));
  expect_near(fields.pressure(2, 0), fields.pressure(2, 1));
  expect_near(fields.pressure(3, grid.ny() + 1), fields.pressure(3, grid.ny()));

  expect_near(fields.pressure_correction(0, 1), fields.pressure_correction(1, 1));
  expect_near(
      fields.pressure_correction(grid.nx() + 1, grid.ny()),
      fields.pressure_correction(grid.nx(), grid.ny()));
  expect_near(fields.pressure_correction(2, 0), fields.pressure_correction(2, 1));
  expect_near(
      fields.pressure_correction(3, grid.ny() + 1), fields.pressure_correction(3, grid.ny()));

  fields.reset();
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

  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 0; i <= grid.nx(); ++i) {
      expect_near(fields.u_face(i, j), 0.0);
    }
  }
  for (int i = 1; i <= grid.nx(); ++i) {
    for (int j = 0; j <= grid.ny(); ++j) {
      expect_near(fields.v_face(i, j), 0.0);
    }
  }

  const cfd::MomentumAssembly zero_u_assembly = cfd::assemble_u_momentum(grid, config, fields);
  const cfd::MomentumAssembly zero_v_assembly = cfd::assemble_v_momentum(grid, config, fields);

  const double diff_x = config.viscosity() * grid.dy() / grid.dx();
  const double diff_y = config.viscosity() * grid.dx() / grid.dy();
  const double u_interior_diag = (2.0 * diff_x + 2.0 * diff_y) / config.controls.alpha_u;
  const double v_interior_diag = (2.0 * diff_x + 2.0 * diff_y) / config.controls.alpha_v;

  const int u_top_left = static_cast<int>(grid.u_index(1, grid.ny()));
  const int u_top_left_east = static_cast<int>(grid.u_index(2, grid.ny()));
  const int u_top_left_south = static_cast<int>(grid.u_index(1, grid.ny() - 1));
  expect_near(zero_u_assembly.diagonal(u_top_left), u_interior_diag + diff_x + diff_y);
  expect_near(sparse_coeff(zero_u_assembly.system.matrix, u_top_left, u_top_left), u_interior_diag + diff_x + diff_y);
  expect_near(sparse_coeff(zero_u_assembly.system.matrix, u_top_left, u_top_left_east), -diff_x);
  expect_near(sparse_coeff(zero_u_assembly.system.matrix, u_top_left, u_top_left_south), -diff_y);
  expect_near(zero_u_assembly.system.rhs(u_top_left), 2.0 * diff_y * config.lid_velocity);

  const int v_top_left = static_cast<int>(grid.v_index(1, grid.ny()));
  const int v_top_left_east = static_cast<int>(grid.v_index(2, grid.ny()));
  const int v_top_left_south = static_cast<int>(grid.v_index(1, grid.ny() - 1));
  expect_near(zero_v_assembly.diagonal(v_top_left), v_interior_diag + diff_x + diff_y);
  expect_near(sparse_coeff(zero_v_assembly.system.matrix, v_top_left, v_top_left), v_interior_diag + diff_x + diff_y);
  expect_near(sparse_coeff(zero_v_assembly.system.matrix, v_top_left, v_top_left_east), -diff_x);
  expect_near(sparse_coeff(zero_v_assembly.system.matrix, v_top_left, v_top_left_south), -diff_y);
  expect_near(zero_v_assembly.system.rhs(v_top_left), 0.0);

  for (int j = 0; j <= grid.ny() + 1; ++j) {
    for (int i = 0; i <= grid.nx() + 1; ++i) {
      fields.pressure(i, j) = static_cast<double>(i * i * i + j * j * j);
    }
  }
  fields.d_u.setConstant(0.25);
  fields.d_v.setConstant(0.5);

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

  expect_near(fields.u_face(2, 2), 1.5);
  expect_near(fields.v_face(2, 2), 3.0);

  for (int j = 0; j <= grid.ny() + 1; ++j) {
    for (int i = 0; i <= grid.nx() + 1; ++i) {
      fields.pressure(i, j) = 3.0 * static_cast<double>(i) - 2.0 * static_cast<double>(j);
    }
  }
  fields.u.setConstant(1.25);
  fields.v.setConstant(-0.75);
  fields.d_u.setConstant(0.8);
  fields.d_v.setConstant(0.6);

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

  expect_near(fields.u_face(2, 2), 1.25);
  expect_near(fields.v_face(2, 2), -0.75);

  const cfd::MomentumAssembly u_assembly = cfd::assemble_u_momentum(grid, config, fields);
  const cfd::MomentumAssembly v_assembly = cfd::assemble_v_momentum(grid, config, fields);

  assert(u_assembly.system.matrix.rows() == static_cast<int>(grid.u_unknown_count()));
  assert(u_assembly.system.matrix.cols() == static_cast<int>(grid.u_unknown_count()));
  assert(v_assembly.system.matrix.rows() == static_cast<int>(grid.v_unknown_count()));
  assert(v_assembly.system.matrix.cols() == static_cast<int>(grid.v_unknown_count()));

  for (int row = 0; row < u_assembly.system.matrix.outerSize(); ++row) {
    double diag = 0.0;
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(u_assembly.system.matrix, row);
         it; ++it) {
      if (it.col() == row) {
        diag = it.value();
      } else {
        assert(it.value() <= 1e-12);
      }
    }
    assert(diag > 0.0);
  }

  for (int row = 0; row < v_assembly.system.matrix.outerSize(); ++row) {
    double diag = 0.0;
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(v_assembly.system.matrix, row);
         it; ++it) {
      if (it.col() == row) {
        diag = it.value();
      } else {
        assert(it.value() <= 1e-12);
      }
    }
    assert(diag > 0.0);
  }

  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 1; i <= grid.nx(); ++i) {
      fields.d_u(i, j) =
          (grid.dx() * grid.dy()) / u_assembly.diagonal(static_cast<int>(grid.u_index(i, j)));
    }
  }
  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 1; i <= grid.nx(); ++i) {
      fields.d_v(i, j) =
          (grid.dx() * grid.dy()) / v_assembly.diagonal(static_cast<int>(grid.v_index(i, j)));
    }
  }

  fields.u_star = fields.u;
  fields.v_star = fields.v;
  cfd::update_face_velocities(
      grid,
      config,
      fields.u_star,
      fields.v_star,
      fields.pressure,
      fields.d_u,
      fields.d_v,
      &fields.u_face,
      &fields.v_face);
  const cfd::PressureCorrectionAssembly p_assembly =
      cfd::assemble_pressure_correction(grid, config, fields);

  assert(p_assembly.system.matrix.rows() == static_cast<int>(grid.pressure_cell_count()));
  assert(p_assembly.system.matrix.cols() == static_cast<int>(grid.pressure_cell_count()));
  assert(p_assembly.system.matrix.coeff(0, 0) == 1.0);
  assert(p_assembly.system.rhs(0) == 0.0);
  assert(p_assembly.mass_imbalance(0) == 0.0);

  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 1; i <= grid.nx(); ++i) {
      const int row = static_cast<int>(grid.pressure_index(i, j));
      if (row == 0) {
        continue;
      }
      expect_near(p_assembly.system.rhs(row), p_assembly.mass_imbalance(row));
    }
  }

  fields.reset();
  expect_near(cfd::compute_continuity_residual(grid, config, fields), 0.0);

  return 0;
}
