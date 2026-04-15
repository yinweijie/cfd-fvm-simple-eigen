#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include <stdexcept>
#include <string>

namespace cfd {

// Sparse linear system Ax = b assembled for Eigen-based iterative solvers.
struct LinearSystem {
  Eigen::SparseMatrix<double, Eigen::RowMajor> matrix;
  Eigen::VectorXd rhs;
};

// Iterative solver choices used by the momentum and pressure-correction systems.
enum class SolverKind {
  kBiCGSTAB,
  kConjugateGradient,
};

// Linear solve output including the solution vector and convergence diagnostics.
struct LinearSolveResult {
  Eigen::VectorXd solution;
  double relative_residual = 0.0;
  int iterations = 0;
};

// Solve a sparse linear system with the requested iterative method and tolerances.
inline LinearSolveResult solve_linear_system(
    const LinearSystem& system,
    SolverKind kind,
    double tolerance,
    int max_iterations) {
  LinearSolveResult result;

  // Use BiCGSTAB for the generally nonsymmetric momentum systems.
  if (kind == SolverKind::kBiCGSTAB) {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.setTolerance(tolerance);
    solver.setMaxIterations(max_iterations);
    solver.compute(system.matrix);
    result.solution = solver.solve(system.rhs);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("BiCGSTAB failed to solve the linear system");
    }
    result.relative_residual = solver.error();
    result.iterations = solver.iterations();
    return result;
  }

  // Use ConjugateGradient for the symmetric pressure-correction system.
  Eigen::ConjugateGradient<
      Eigen::SparseMatrix<double, Eigen::RowMajor>,
      Eigen::Lower | Eigen::Upper>
      solver;
  solver.setTolerance(tolerance);
  solver.setMaxIterations(max_iterations);
  solver.compute(system.matrix);
  result.solution = solver.solve(system.rhs);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("ConjugateGradient failed to solve the linear system");
  }
  result.relative_residual = solver.error();
  result.iterations = solver.iterations();
  return result;
}

}  // namespace cfd
