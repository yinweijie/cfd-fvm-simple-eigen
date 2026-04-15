#pragma once

#include <string>
#include <utility>
#include <vector>

#include "cfd/simple_solver.hpp"

namespace cfd {

// Write interior-cell fields, centerlines, and convergence history to the output directory.
void write_results(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields,
    const SolveSummary& summary,
    const std::string& output_dir);

// Interpolate the cell-centered u field to the cavity vertical centerline x = 0.5.
std::vector<std::pair<double, double>> extract_u_centerline(
    const StructuredGrid& grid,
    const FlowFields& fields);

// Interpolate the cell-centered v field to the cavity horizontal centerline y = 0.5.
std::vector<std::pair<double, double>> extract_v_centerline(
    const StructuredGrid& grid,
    const FlowFields& fields);

}  // namespace cfd
