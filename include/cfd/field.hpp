#pragma once

#include <Eigen/Core>

#include "cfd/mesh.hpp"

namespace cfd {

struct FlowFields {
  explicit FlowFields(const StructuredGrid& grid);

  void reset();

  Eigen::MatrixXd pressure;
  Eigen::MatrixXd pressure_correction;

  Eigen::MatrixXd u;
  Eigen::MatrixXd v;
  Eigen::MatrixXd u_star;
  Eigen::MatrixXd v_star;

  Eigen::MatrixXd d_u;
  Eigen::MatrixXd d_v;
  Eigen::MatrixXd u_face;
  Eigen::MatrixXd v_face;
};

}  // namespace cfd
