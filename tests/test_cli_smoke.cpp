#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

bool has_solver_marker(const std::filesystem::path& summary_path, const std::string& expected) {
  std::ifstream in(summary_path);
  std::string line;
  while (std::getline(in, line)) {
    if (line == "solver=" + expected) {
      return true;
    }
  }
  return false;
}

bool has_standard_outputs(const std::filesystem::path& output_dir) {
  return std::filesystem::exists(output_dir / "u.csv") &&
         std::filesystem::exists(output_dir / "v.csv") &&
         std::filesystem::exists(output_dir / "p.csv") &&
         std::filesystem::exists(output_dir / "residuals.csv") &&
         std::filesystem::exists(output_dir / "summary.txt");
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Expected solver path and output dir\n";
    return 2;
  }

  const std::filesystem::path solver_path = argv[1];
  const std::filesystem::path output_dir = argv[2];
  std::filesystem::remove_all(output_dir);

  const std::filesystem::path default_output = output_dir / "default";
  const std::string default_command =
      solver_path.string() +
      " --case cavity --nx 8 --ny 8 --re 100 --max-iters 250 "
      "--min-iters 10 --projection-dt 0.15 --output-dir " +
      default_output.string();

  if (std::system(default_command.c_str()) != 0) {
    return 1;
  }
  if (!has_standard_outputs(default_output) ||
      !has_solver_marker(default_output / "summary.txt", "projection")) {
    return 1;
  }

  const std::filesystem::path simple_output = output_dir / "simple";
  const std::string simple_command =
      solver_path.string() +
      " --case cavity --solver simple --nx 8 --ny 8 --re 100 --max-iters 100 "
      "--alpha-u 0.5 --alpha-v 0.5 --alpha-p 0.3 --output-dir " +
      simple_output.string();
  if (std::system(simple_command.c_str()) != 0) {
    return 1;
  }
  if (!has_standard_outputs(simple_output) ||
      !has_solver_marker(simple_output / "summary.txt", "simple")) {
    return 1;
  }

  const std::filesystem::path projection_output = output_dir / "projection";
  const std::string projection_command =
      solver_path.string() +
      " --case cavity --solver projection --nx 8 --ny 8 --re 100 --max-iters 250 "
      "--min-iters 10 --projection-dt 0.15 --output-dir " +
      projection_output.string();
  if (std::system(projection_command.c_str()) != 0) {
    return 1;
  }
  if (!has_standard_outputs(projection_output) ||
      !has_solver_marker(projection_output / "summary.txt", "projection")) {
    return 1;
  }

  const std::string invalid_command =
      solver_path.string() +
      " --case cavity --solver invalid --nx 8 --ny 8 --output-dir " +
      (output_dir / "invalid").string();
  if (std::system(invalid_command.c_str()) == 0) {
    return 1;
  }

  return 0;
}
