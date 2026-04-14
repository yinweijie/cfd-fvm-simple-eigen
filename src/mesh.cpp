#include "cfd/mesh.hpp"

#include <stdexcept>

namespace cfd {

StructuredGrid::StructuredGrid(MeshSpec spec) : spec_(spec), dx_(0.0), dy_(0.0) {
  if (spec_.nx <= 0 || spec_.ny <= 0) {
    throw std::invalid_argument("StructuredGrid requires nx > 0 and ny > 0");
  }
  if (spec_.lx <= 0.0 || spec_.ly <= 0.0) {
    throw std::invalid_argument("StructuredGrid requires lx > 0 and ly > 0");
  }

  dx_ = spec_.lx / static_cast<double>(spec_.nx);
  dy_ = spec_.ly / static_cast<double>(spec_.ny);
}

std::size_t StructuredGrid::pressure_cell_count() const noexcept {
  return static_cast<std::size_t>(spec_.nx) * static_cast<std::size_t>(spec_.ny);
}

std::size_t StructuredGrid::u_unknown_count() const noexcept {
  return pressure_cell_count();
}

std::size_t StructuredGrid::v_unknown_count() const noexcept {
  return pressure_cell_count();
}

std::size_t StructuredGrid::pressure_index(int i, int j) const {
  validate_pressure_ij(i, j);
  return static_cast<std::size_t>(j - 1) * static_cast<std::size_t>(spec_.nx) +
         static_cast<std::size_t>(i - 1);
}

std::size_t StructuredGrid::u_index(int i, int j) const {
  return pressure_index(i, j);
}

std::size_t StructuredGrid::v_index(int i, int j) const {
  return pressure_index(i, j);
}

bool StructuredGrid::is_boundary_cell(int i, int j) const noexcept {
  return i == 1 || i == spec_.nx || j == 1 || j == spec_.ny;
}

double StructuredGrid::cell_center_x(int i) const {
  validate_pressure_ij(i, 1);
  return (static_cast<double>(i) - 0.5) * dx_;
}

double StructuredGrid::cell_center_y(int j) const {
  validate_pressure_ij(1, j);
  return (static_cast<double>(j) - 0.5) * dy_;
}

void StructuredGrid::validate_pressure_ij(int i, int j) const {
  if (i < 1 || i > spec_.nx || j < 1 || j > spec_.ny) {
    throw std::out_of_range("pressure index out of range");
  }
}

}  // namespace cfd
