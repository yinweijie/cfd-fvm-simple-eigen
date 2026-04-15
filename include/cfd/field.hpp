#pragma once

#include <Eigen/Core>

#include "cfd/mesh.hpp"

namespace cfd {

// Collocated-grid SIMPLE fields, including ghost cells and reconstructed face values.
struct FlowFields {
  // Allocate all cell-centered and face-centered arrays for the given grid.
  explicit FlowFields(const StructuredGrid& grid);

  // Reset pressure, velocity, correction, and face fields to zero.
  void reset();

  // Cell-centered pressure and its SIMPLE correction increment.
  Eigen::MatrixXd pressure;
  Eigen::MatrixXd pressure_correction;

  // Cell-centered corrected velocity and momentum-predicted velocity.
  Eigen::MatrixXd u;
  Eigen::MatrixXd v;
  Eigen::MatrixXd u_star;
  Eigen::MatrixXd v_star;

  // SIMPLE velocity-correction factors and Rhie-Chow face velocities.
  // Cell fields use (nx + 2, ny + 2) storage with one ghost-cell layer.
  // Face fields store the physical-face direction plus one ghost-aligned index in the transverse
  // direction: u_face on vertical faces, v_face on horizontal faces.
  Eigen::MatrixXd d_u;
  Eigen::MatrixXd d_v;
  Eigen::MatrixXd u_face;
  Eigen::MatrixXd v_face;
};

}  // namespace cfd
