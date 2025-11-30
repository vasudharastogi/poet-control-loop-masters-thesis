// #include "IO/StatsIO.hpp"
#include "Control/ControlModule.hpp"
#include "IO/HDF5Functions.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

int write_metrics(const std::vector<poet::CellMetrics> &metrics_history,
                  const std::vector<std::string> &species_names,
                  const std::string &dir_path, const std::string &file_name) {

  if (!fs::exists(dir_path)) {
    std::cerr << "Directory does not exist: " << dir_path << std::endl;
    return -1;
  }
  fs::path file_path = fs::path(dir_path) / file_name;

  // Use a std::string path to avoid filesystem path conversion issues.
  H5Easy::File file(file_path.string(), H5Easy::File::Truncate);

  for (size_t idx = 0; idx < metrics_history.size(); ++idx) {
    const auto &metrics = metrics_history[idx];
    /*
    std::string grp = "entry_" + std::to_string(idx) + "_iter_" +
                      std::to_string(metrics.iteration) + "_rb_" +
                      std::to_string(metrics.rollback_count);
    */
    std::string grp = "iter_" + std::to_string(metrics.iteration) + "_rb_" +
                      std::to_string(metrics.rb_count);

    size_t n_cells = metrics.id.size();
    // Use provided species_names as the source of truth to avoid OOB when mape is empty.
    size_t n_species = species_names.size();

    // Create a scalar dataset "meta" and attach attributes to it (no explicit groups).
    H5Easy::dump(file, grp + "/meta", 0, H5Easy::DumpMode::Overwrite);

    H5Easy::dumpAttribute(file, grp + "/meta", "species_names", species_names,
                          H5Easy::DumpMode::Overwrite);
    H5Easy::dumpAttribute(file, grp + "/meta", "iteration", metrics.iteration,
                          H5Easy::DumpMode::Overwrite);
    H5Easy::dumpAttribute(file, grp + "/meta", "rollback_count",
                          metrics.rb_count, H5Easy::DumpMode::Overwrite);
    H5Easy::dumpAttribute(file, grp + "/meta", "n_cells", n_cells,
                          H5Easy::DumpMode::Overwrite);
    H5Easy::dumpAttribute(file, grp + "/meta", "n_species", n_species,
                          H5Easy::DumpMode::Overwrite);

    // Dump only if data is non-empty to avoid corrupting the file on failures.
    if (!metrics.mape.empty()) {
      H5Easy::dump(file, grp + "/mape", metrics.mape,
                   H5Easy::DumpMode::Overwrite);
    }
    if (!metrics.rrmse.empty()) {
      H5Easy::dump(file, grp + "/rrmse", metrics.rrmse,
                   H5Easy::DumpMode::Overwrite);
    }
  }

  // Ensure all buffers are flushed to disk.
  file.flush();

  return 0;
}

void writeCellStatsToCSV(const std::vector<poet::CellMetrics> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &out_dir,
                         const std::string &filename) {
  std::ofstream out(std::filesystem::path(out_dir) / filename);
  if (!out.is_open()) {
    std::cerr << "Could not open " << filename << " !" << std::endl;
    return;
  }

  // Header
  out << std::left << std::setw(15) << "CellID" << std::setw(15) << "Iteration"
      << std::setw(15) << "Rollback" << std::setw(15) << "Species"
      << std::setw(15) << "MAPE" << std::setw(15) << "RRMSE"
      << "\n"
      << std::string(90, '-') << "\n";

  // Data rows (fix column ordering: include rb_count before Species)
  for (const auto &metrics : all_stats) {
    for (size_t cell_idx = 0; cell_idx < metrics.id.size(); ++cell_idx) {
      for (size_t sp_idx = 0; sp_idx < species_names.size(); ++sp_idx) {
        out << std::left << std::setw(15) << metrics.id[cell_idx]
            << std::setw(15) << metrics.iteration
            << std::setw(15) << metrics.rb_count
            << std::setw(15) << species_names[sp_idx]
            << std::setw(15) << metrics.mape[cell_idx][sp_idx]
            << std::setw(15) << metrics.rrmse[cell_idx][sp_idx] << "\n";
      }
    }
    out << "\n";
  }
  out.close();
  std::cout << "Cell error metrics written to " << out_dir << "/" << filename
            << "\n";
}

void writeSpeciesStatsToCSV(
    const std::vector<poet::SpeciesMetrics> &all_stats,
    const std::vector<std::string> &species_names, const std::string &out_dir,
    const std::string &filename) {
  std::ofstream out(std::filesystem::path(out_dir) / filename);
  if (!out.is_open()) {
    std::cerr << "Could not open " << filename << " !" << std::endl;
    return;
  }

  // Header
  out << std::left << std::setw(15) << "Iteration" << std::setw(15)
      << "Rollback" << std::setw(15) << "Species" << std::setw(15) << "MAPE"
      << std::setw(15) << "RRMSE"
      << "\n"
      << std::string(75, '-') << "\n";

  // Data rows
  for (const auto &metrics : all_stats) {
    for (size_t sp_idx = 0; sp_idx < species_names.size(); ++sp_idx) {
      out << std::left << std::setw(15) << metrics.iteration << std::setw(15)
          << metrics.rb_count << std::setw(15) << species_names[sp_idx]
          << std::setw(15) << metrics.mape[sp_idx] << std::setw(15)
          << metrics.rrmse[sp_idx] << "\n";
    }
    out << "\n";
  }
  out.close();
  std::cout << "Species error metrics written to " << out_dir << "/" << filename
            << "\n";
}
