#pragma once

#include "cfd/case.hpp"
#include "cfd/field.hpp"
#include "cfd/linear_system.hpp"

namespace cfd {

struct MomentumAssembly {
  LinearSystem system;
  Eigen::VectorXd diagonal;
};

struct PressureCorrectionAssembly {
  LinearSystem system;
  Eigen::VectorXd mass_imbalance;
};

void apply_cavity_boundary_conditions(
    const StructuredGrid& grid,
    const CavityCase& config,
    FlowFields* fields);

void update_face_velocities(
    const StructuredGrid& grid,
    const CavityCase& config,
    const Eigen::MatrixXd& u_cells,
    const Eigen::MatrixXd& v_cells,
    const Eigen::MatrixXd& pressure,
    const Eigen::MatrixXd& d_u,
    const Eigen::MatrixXd& d_v,
    Eigen::MatrixXd* u_face,
    Eigen::MatrixXd* v_face);

MomentumAssembly assemble_u_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

MomentumAssembly assemble_v_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

PressureCorrectionAssembly assemble_pressure_correction(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

double compute_continuity_residual(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

}  // namespace cfd
