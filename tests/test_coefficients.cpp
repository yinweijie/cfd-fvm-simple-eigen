#include "cfd/case.hpp"
#include "cfd/discretization.hpp"

#include <cassert>

int main() {
  cfd::CavityCase config;
  config.mesh_spec = {4, 4, 1.0, 1.0};
  config.controls.alpha_u = 0.5;
  config.controls.alpha_v = 0.5;
  config.controls.alpha_p = 0.3;

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

  return 0;
}
