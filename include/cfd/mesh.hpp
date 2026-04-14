#pragma once

#include <cstddef>
#include <stdexcept>

namespace cfd {

struct MeshSpec {
  int nx;
  int ny;
  double lx;
  double ly;
};

class StructuredGrid {
 public:
  explicit StructuredGrid(MeshSpec spec);

  const MeshSpec& spec() const noexcept { return spec_; }

  int nx() const noexcept { return spec_.nx; }
  int ny() const noexcept { return spec_.ny; }
  double dx() const noexcept { return dx_; }
  double dy() const noexcept { return dy_; }

  std::size_t pressure_cell_count() const noexcept;
  std::size_t u_unknown_count() const noexcept;
  std::size_t v_unknown_count() const noexcept;

  std::size_t pressure_index(int i, int j) const;
  std::size_t u_index(int i, int j) const;
  std::size_t v_index(int i, int j) const;

  bool is_boundary_cell(int i, int j) const noexcept;

  double cell_center_x(int i) const;
  double cell_center_y(int j) const;

 private:
  MeshSpec spec_;
  double dx_;
  double dy_;

  void validate_pressure_ij(int i, int j) const;
};

}  // namespace cfd
