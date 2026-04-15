#pragma once

#include "cfd/case.hpp"
#include "cfd/field.hpp"
#include "cfd/linear_system.hpp"

namespace cfd {

// Sparse momentum system together with the relaxed diagonal used in SIMPLE corrections.
struct MomentumAssembly {
  LinearSystem system;
  Eigen::VectorXd diagonal;
};

// Pressure-correction system and the cell-wise continuity defect that drives its RHS.
struct PressureCorrectionAssembly {
  LinearSystem system;
  Eigen::VectorXd mass_imbalance;
};

// Apply cavity-wall ghost-cell updates for pressure, pressure_correction, and velocity fields.
void apply_cavity_boundary_conditions(
    const StructuredGrid& grid,
    const CavityCase& config,
    FlowFields* fields);

// Reconstruct face velocities from cell-centered fields using a Rhie-Chow-style pressure correction.
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

// Assemble the cell-centered u-momentum linear system.
MomentumAssembly assemble_u_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

// Assemble the cell-centered v-momentum linear system.
MomentumAssembly assemble_v_momentum(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

// Assemble the SIMPLE pressure-correction equation from predicted face mass imbalance.
PressureCorrectionAssembly assemble_pressure_correction(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

// Compute the average absolute cell-wise continuity defect of the current corrected fields.
double compute_continuity_residual(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields);

}  // namespace cfd
