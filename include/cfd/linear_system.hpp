#pragma once

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include <stdexcept>
#include <string>

namespace cfd {

struct LinearSystem {
  Eigen::SparseMatrix<double, Eigen::RowMajor> matrix;
  Eigen::VectorXd rhs;
};

enum class SolverKind {
  kBiCGSTAB,
  kConjugateGradient,
};

struct LinearSolveResult {
  Eigen::VectorXd solution;
  double relative_residual = 0.0;
  int iterations = 0;
};

inline LinearSolveResult solve_linear_system(
    const LinearSystem& system,
    SolverKind kind,
    double tolerance,
    int max_iterations) {
  LinearSolveResult result;

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
