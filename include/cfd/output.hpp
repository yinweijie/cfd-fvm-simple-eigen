#pragma once

#include <string>
#include <utility>
#include <vector>

#include "cfd/simple_solver.hpp"

namespace cfd {

void write_results(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields,
    const SolveSummary& summary,
    const std::string& output_dir);

std::vector<std::pair<double, double>> extract_u_centerline(
    const StructuredGrid& grid,
    const FlowFields& fields);

std::vector<std::pair<double, double>> extract_v_centerline(
    const StructuredGrid& grid,
    const FlowFields& fields);

}  // namespace cfd
