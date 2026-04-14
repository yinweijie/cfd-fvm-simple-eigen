#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "cfd/case.hpp"
#include "cfd/mesh.hpp"
#include "cfd/output.hpp"
#include "cfd/simple_solver.hpp"

namespace {

double parse_double(const std::string& value, const char* option_name) {
  try {
    return std::stod(value);
  } catch (const std::exception&) {
    throw std::invalid_argument(std::string("Invalid value for ") + option_name);
  }
}

int parse_int(const std::string& value, const char* option_name) {
  try {
    return std::stoi(value);
  } catch (const std::exception&) {
    throw std::invalid_argument(std::string("Invalid value for ") + option_name);
  }
}

void print_usage() {
  std::cout
      << "Usage: cfd_solver [--case cavity] [--nx N] [--ny N] [--re Re] "
      << "[--max-iters N] [--min-iters N] [--alpha-u A] [--alpha-v A] [--alpha-p A] "
      << "[--output-dir DIR]\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    cfd::CavityCase config;
    std::string case_name = "cavity";

    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--help") {
        print_usage();
        return 0;
      }
      if (i + 1 >= argc) {
        throw std::invalid_argument("Missing value for " + arg);
      }
      const std::string value = argv[++i];

      if (arg == "--case") {
        case_name = value;
      } else if (arg == "--nx") {
        config.mesh_spec.nx = parse_int(value, "--nx");
      } else if (arg == "--ny") {
        config.mesh_spec.ny = parse_int(value, "--ny");
      } else if (arg == "--re") {
        config.reynolds = parse_double(value, "--re");
      } else if (arg == "--max-iters") {
        config.controls.max_iterations = parse_int(value, "--max-iters");
      } else if (arg == "--min-iters") {
        config.controls.min_iterations = parse_int(value, "--min-iters");
      } else if (arg == "--alpha-u") {
        config.controls.alpha_u = parse_double(value, "--alpha-u");
      } else if (arg == "--alpha-v") {
        config.controls.alpha_v = parse_double(value, "--alpha-v");
      } else if (arg == "--alpha-p") {
        config.controls.alpha_p = parse_double(value, "--alpha-p");
      } else if (arg == "--output-dir") {
        config.output_dir = value;
      } else {
        throw std::invalid_argument("Unknown option: " + arg);
      }
    }

    if (case_name != "cavity") {
      throw std::invalid_argument("Only --case cavity is currently supported");
    }

    if (config.output_dir.empty()) {
      config.output_dir = "results/cavity_re" + std::to_string(static_cast<int>(config.reynolds));
    }

    cfd::SimpleSolver solver(config);
    cfd::SolveSummary summary = solver.run();
    cfd::write_results(solver.grid(), solver.config(), solver.fields(), summary, config.output_dir);

    const cfd::IterationMetrics& final_metrics = summary.residual_history.back();
    std::cout << "Case: cavity\n";
    std::cout << "Grid: " << config.mesh_spec.nx << "x" << config.mesh_spec.ny << "\n";
    std::cout << "Re: " << config.reynolds << "\n";
    std::cout << "Converged: " << (summary.converged ? "true" : "false") << "\n";
    std::cout << "Iterations: " << summary.iterations << "\n";
    std::cout << "Continuity residual: " << final_metrics.continuity_residual << "\n";
    std::cout << "U residual: " << final_metrics.u_momentum_residual << "\n";
    std::cout << "V residual: " << final_metrics.v_momentum_residual << "\n";
    std::cout << "Output: " << std::filesystem::absolute(config.output_dir).string() << "\n";

    return summary.converged ? 0 : 1;
  } catch (const std::exception& ex) {
    std::cerr << "cfd_solver error: " << ex.what() << "\n";
    print_usage();
    return 2;
  }
}
