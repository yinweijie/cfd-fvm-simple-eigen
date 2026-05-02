cmake_minimum_required(VERSION 3.24)

project(cfd_fvm_simple_eigen VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(ResolveEigen)

cfd_resolve_eigen(CFD_EIGEN_INCLUDE_DIR)

add_library(cfd_eigen INTERFACE)
add_library(cfd::eigen ALIAS cfd_eigen)
target_include_directories(cfd_eigen SYSTEM INTERFACE "${CFD_EIGEN_INCLUDE_DIR}")

add_library(cfd_core STATIC
  src/mesh.cpp
  src/discretization.cpp
  src/simple_solver.cpp
  src/projection_solver.cpp
  src/output.cpp
)
add_library(cfd::core ALIAS cfd_core)

target_include_directories(cfd_core
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(cfd_core PUBLIC cfd::eigen)
target_compile_features(cfd_core PUBLIC cxx_std_17)

add_executable(cfd_solver
  src/main.cpp
)
target_link_libraries(cfd_solver PRIVATE cfd::core)

enable_testing()

add_executable(test_mesh
  tests/test_mesh.cpp
)
target_link_libraries(test_mesh PRIVATE cfd::core)
add_test(NAME test_mesh COMMAND test_mesh)

add_executable(test_coefficients
  tests/test_coefficients.cpp
)
target_link_libraries(test_coefficients PRIVATE cfd::core)
add_test(NAME test_coefficients COMMAND test_coefficients)

add_executable(test_simple_smoke
  tests/test_simple_smoke.cpp
)
target_link_libraries(test_simple_smoke PRIVATE cfd::core)
add_test(NAME test_simple_smoke COMMAND test_simple_smoke)

add_executable(test_projection_coefficients
  tests/test_projection_coefficients.cpp
)
target_link_libraries(test_projection_coefficients PRIVATE cfd::core)
add_test(NAME test_projection_coefficients COMMAND test_projection_coefficients)

add_executable(test_projection_smoke
  tests/test_projection_smoke.cpp
)
target_link_libraries(test_projection_smoke PRIVATE cfd::core)
add_test(NAME test_projection_smoke COMMAND test_projection_smoke)

add_executable(test_cli_smoke
  tests/test_cli_smoke.cpp
)
target_link_libraries(test_cli_smoke PRIVATE cfd::core)
add_test(
  NAME test_cli_smoke
  COMMAND test_cli_smoke
    $<TARGET_FILE:cfd_solver>
    ${CMAKE_CURRENT_BINARY_DIR}/test_cli_smoke_output
)

find_program(PYTHON3_EXECUTABLE python3)
if(PYTHON3_EXECUTABLE)
  add_test(
    NAME test_projection_validation_run
    COMMAND $<TARGET_FILE:cfd_solver>
      --case cavity
      --solver projection
      --nx 32
      --ny 32
      --re 100
      --max-iters 1000
      --min-iters 100
      --projection-dt 0.1
      --output-dir ${CMAKE_CURRENT_BINARY_DIR}/test_projection_validation_output
  )
  add_test(
    NAME test_projection_validation
    COMMAND ${PYTHON3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/scripts/validate_cavity.py
      --results ${CMAKE_CURRENT_BINARY_DIR}/test_projection_validation_output
  )
  set_tests_properties(test_projection_validation PROPERTIES
    DEPENDS test_projection_validation_run
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
endif()
