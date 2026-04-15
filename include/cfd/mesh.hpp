#pragma once

#include <cstddef>
#include <stdexcept>

namespace cfd {

// Uniform structured mesh specification for a rectangular cavity domain.
struct MeshSpec {
  int nx;
  int ny;
  double lx;
  double ly;
};

// Store uniform cell geometry and index mappings for the collocated cavity grid.
class StructuredGrid {
 public:
  // Validate the mesh specification and precompute uniform cell sizes.
  explicit StructuredGrid(MeshSpec spec);

  const MeshSpec& spec() const noexcept { return spec_; }

  int nx() const noexcept { return spec_.nx; }
  int ny() const noexcept { return spec_.ny; }
  double dx() const noexcept { return dx_; }
  double dy() const noexcept { return dy_; }

  // Return the number of interior pressure control volumes.
  std::size_t pressure_cell_count() const noexcept;
  // Return the number of collocated u-momentum unknowns.
  std::size_t u_unknown_count() const noexcept;
  // Return the number of collocated v-momentum unknowns.
  std::size_t v_unknown_count() const noexcept;

  // Map one interior pressure cell to its flattened linear-system row.
  std::size_t pressure_index(int i, int j) const;
  // Map one interior u-momentum cell to its flattened linear-system row.
  std::size_t u_index(int i, int j) const;
  // Map one interior v-momentum cell to its flattened linear-system row.
  std::size_t v_index(int i, int j) const;

  // Report whether an interior cell touches a physical cavity wall.
  bool is_boundary_cell(int i, int j) const noexcept;

  // Return the x coordinate of interior cell center i.
  double cell_center_x(int i) const;
  // Return the y coordinate of interior cell center j.
  double cell_center_y(int j) const;

 private:
  MeshSpec spec_;
  double dx_;
  double dy_;

  // Reject indices that do not refer to an interior pressure control volume.
  void validate_pressure_ij(int i, int j) const;
};

}  // namespace cfd
