#include "cfd/output.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>

namespace cfd {
namespace {

// Write a 2-column CSV series with fixed headers.
void write_pairs(
    const std::filesystem::path& path,
    const char* header_x,
    const char* header_y,
    const std::vector<std::pair<double, double>>& values) {
  std::ofstream out(path);
  out << header_x << "," << header_y << "\n";
  out << std::setprecision(10);
  for (const auto& [x, y] : values) {
    out << x << "," << y << "\n";
  }
}

// Write a dense matrix to CSV using the repo's column-major viewing convention.
void write_matrix_csv(
    const std::filesystem::path& path,
    const Eigen::MatrixXd& matrix) {
  std::ofstream out(path);
  out << std::setprecision(10);
  for (int j = 0; j < matrix.cols(); ++j) {
    for (int i = 0; i < matrix.rows(); ++i) {
      if (i > 0) {
        out << ",";
      }
      out << matrix(i, j);
    }
    out << "\n";
  }
}

// Extract the interior pressure field and drop the ghost cells.
Eigen::MatrixXd pressure_cells(const StructuredGrid& grid, const FlowFields& fields) {
  // Only export interior control-volume values; ghost cells only support the stencil.
  Eigen::MatrixXd cells(grid.nx(), grid.ny());
  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 1; i <= grid.nx(); ++i) {
      cells(i - 1, j - 1) = fields.pressure(i, j);
    }
  }
  return cells;
}

// Extract the interior corrected u field and drop the ghost cells.
Eigen::MatrixXd u_cells(const StructuredGrid& grid, const FlowFields& fields) {
  Eigen::MatrixXd cells(grid.nx(), grid.ny());
  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 1; i <= grid.nx(); ++i) {
      cells(i - 1, j - 1) = fields.u(i, j);
    }
  }
  return cells;
}

// Extract the interior corrected v field and drop the ghost cells.
Eigen::MatrixXd v_cells(const StructuredGrid& grid, const FlowFields& fields) {
  Eigen::MatrixXd cells(grid.nx(), grid.ny());
  for (int j = 1; j <= grid.ny(); ++j) {
    for (int i = 1; i <= grid.nx(); ++i) {
      cells(i - 1, j - 1) = fields.v(i, j);
    }
  }
  return cells;
}

}  // namespace

// Interpolate the corrected u field to the cavity's vertical centerline.
std::vector<std::pair<double, double>> extract_u_centerline(
    const StructuredGrid& grid,
    const FlowFields& fields) {
  std::vector<std::pair<double, double>> values;
  values.reserve(static_cast<std::size_t>(grid.ny()));

  const int left_i = grid.nx() / 2;
  const int right_i = left_i + 1;
  const double left_x = grid.cell_center_x(left_i);
  const double right_x = grid.cell_center_x(right_i);
  // Benchmark samples lie at x = 0.5, which is between two cell centers on an even grid.
  const double weight = (0.5 - left_x) / (right_x - left_x);

  for (int j = 1; j <= grid.ny(); ++j) {
    const double y = grid.cell_center_y(j);
    const double left_u = fields.u(left_i, j);
    const double right_u = fields.u(right_i, j);
    const double value = left_u + weight * (right_u - left_u);
    values.emplace_back(y, value);
  }

  return values;
}

// Interpolate the corrected v field to the cavity's horizontal centerline.
std::vector<std::pair<double, double>> extract_v_centerline(
    const StructuredGrid& grid,
    const FlowFields& fields) {
  std::vector<std::pair<double, double>> values;
  values.reserve(static_cast<std::size_t>(grid.nx()));

  const int bottom_j = grid.ny() / 2;
  const int top_j = bottom_j + 1;
  const double bottom_y = grid.cell_center_y(bottom_j);
  const double top_y = grid.cell_center_y(top_j);
  // Benchmark samples lie at y = 0.5, which is between two cell centers on an even grid.
  const double weight = (0.5 - bottom_y) / (top_y - bottom_y);

  for (int i = 1; i <= grid.nx(); ++i) {
    const double x = grid.cell_center_x(i);
    const double bottom_v = fields.v(i, bottom_j);
    const double top_v = fields.v(i, top_j);
    const double value = bottom_v + weight * (top_v - bottom_v);
    values.emplace_back(x, value);
  }

  return values;
}

// Write field snapshots, centerlines, and residual history for one finished run.
void write_results(
    const StructuredGrid& grid,
    const CavityCase& config,
    const FlowFields& fields,
    const SolveSummary& summary,
    const std::string& output_dir) {
  const std::filesystem::path out_dir(output_dir);
  std::filesystem::create_directories(out_dir);

  write_matrix_csv(out_dir / "u.csv", u_cells(grid, fields));
  write_matrix_csv(out_dir / "v.csv", v_cells(grid, fields));
  write_matrix_csv(out_dir / "p.csv", pressure_cells(grid, fields));

  write_pairs(out_dir / "centerline_u.csv", "y", "u", extract_u_centerline(grid, fields));
  write_pairs(out_dir / "centerline_v.csv", "x", "v", extract_v_centerline(grid, fields));

  std::ofstream residuals(out_dir / "residuals.csv");
  residuals << "iteration,continuity,u_momentum,v_momentum,pressure_correction,max_velocity_correction\n";
  residuals << std::setprecision(10);
  for (const IterationMetrics& item : summary.residual_history) {
    residuals << item.iteration << ","
              << item.continuity_residual << ","
              << item.u_momentum_residual << ","
              << item.v_momentum_residual << ","
              << item.pressure_correction_residual << ","
              << item.max_velocity_correction << "\n";
  }

  std::ofstream summary_out(out_dir / "summary.txt");
  summary_out << "converged=" << (summary.converged ? "true" : "false") << "\n";
  summary_out << "iterations=" << summary.iterations << "\n";
  summary_out << "re=" << config.reynolds << "\n";
  summary_out << "nx=" << config.mesh_spec.nx << "\n";
  summary_out << "ny=" << config.mesh_spec.ny << "\n";
}

}  // namespace cfd
