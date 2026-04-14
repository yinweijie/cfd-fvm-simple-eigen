#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Expected solver path and output dir\n";
    return 2;
  }

  const std::filesystem::path solver_path = argv[1];
  const std::filesystem::path output_dir = argv[2];
  std::filesystem::remove_all(output_dir);

  const std::string command =
      solver_path.string() +
      " --case cavity --nx 8 --ny 8 --re 100 --max-iters 100 "
      "--alpha-u 0.5 --alpha-v 0.5 --alpha-p 0.3 --output-dir " +
      output_dir.string();

  const int rc = std::system(command.c_str());
  if (rc != 0) {
    return 1;
  }

  if (!std::filesystem::exists(output_dir / "u.csv") ||
      !std::filesystem::exists(output_dir / "v.csv") ||
      !std::filesystem::exists(output_dir / "p.csv") ||
      !std::filesystem::exists(output_dir / "residuals.csv")) {
    return 1;
  }

  return 0;
}
