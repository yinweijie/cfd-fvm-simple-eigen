#include "cfd/mesh.hpp"

#include <cassert>
#include <cmath>

int main() {
  cfd::StructuredGrid grid({4, 3, 2.0, 1.5});

  assert(grid.nx() == 4);
  assert(grid.ny() == 3);
  assert(std::abs(grid.dx() - 0.5) < 1e-12);
  assert(std::abs(grid.dy() - 0.5) < 1e-12);
  assert(grid.pressure_cell_count() == 12);
  assert(grid.u_unknown_count() == 12);
  assert(grid.v_unknown_count() == 12);
  assert(grid.pressure_index(1, 1) == 0);
  assert(grid.pressure_index(4, 3) == 11);
  assert(grid.u_index(1, 1) == 0);
  assert(grid.u_index(4, 3) == 11);
  assert(grid.v_index(1, 1) == 0);
  assert(grid.v_index(4, 3) == 11);
  assert(grid.is_boundary_cell(1, 2));
  assert(!grid.is_boundary_cell(2, 2));
  assert(std::abs(grid.cell_center_x(1) - 0.25) < 1e-12);
  assert(std::abs(grid.cell_center_y(3) - 1.25) < 1e-12);

  return 0;
}
