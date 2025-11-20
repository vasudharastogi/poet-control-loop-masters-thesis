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

int write_metrics(const std::vector<poet::SpeciesErrorMetrics> &metrics_history,
                  const std::vector<std::string> &species_names,
                  const std::string &dir_path, const std::string &file_name) {

  if (!fs::exists(dir_path)) {
    std::cerr << "Directory does not exist: " << dir_path << std::endl;
    return -1;
  }
  fs::path file_path = fs::path(dir_path) / file_name;

  H5Easy::File file(file_path, H5Easy::File::Truncate);

  for (size_t idx = 0; idx < metrics_history.size(); ++idx) {
    const auto &metrics = metrics_history[idx];
    std::string grp = "iter_" + std::to_string(metrics.iteration) + "_" +
                      std::to_string(metrics.rollback_count);

    size_t n_cells = metrics.id.size();
    size_t n_species = metrics.mape[0].size();

    H5Easy::dump(file, grp + "/meta", 0);

    // Attach attributes
    H5Easy::dumpAttribute(file, grp + "/meta", "species_names", species_names);
    H5Easy::dumpAttribute(file, grp + "/meta", "iteration", metrics.iteration);
    H5Easy::dumpAttribute(file, grp + "/meta", "rollback_count",
                          metrics.rollback_count);
    H5Easy::dumpAttribute(file, grp + "/meta", "n_cells", n_cells);
    H5Easy::dumpAttribute(file, grp + "/meta", "n_species", n_species);

    // ─────────────────────────────────────────────
    // 2. Real datasets
    // ─────────────────────────────────────────────
    H5Easy::dump(file, grp + "/cell_id", metrics.id);
    H5Easy::dump(file, grp + "/mape", metrics.mape);
    H5Easy::dump(file, grp + "/rrmse", metrics.rrmse);
  }

  return 0;
}

void writeStatsToCSV(const std::vector<poet::SpeciesErrorMetrics> &all_stats,
                     const std::vector<std::string> &species_names,
                     const std::string &out_dir, const std::string &filename) {
  std::filesystem::path full_path = std::filesystem::path(out_dir) / filename;

  std::ofstream out(full_path);
  if (!out.is_open()) {
    std::cerr << "Could not open " << filename << " !" << std::endl;
    return;
  }

  // header: CellID, Iteration, Rollback, Species, MAPE, RRMSE
  out << std::left << std::setw(15) << "CellID" << std::setw(15) << "Iteration"
      << std::setw(15) << "Rollback" << std::setw(15) << "Species"
      << std::setw(15) << "MAPE" << std::setw(15) << "RRMSE"
      << "\n";

  out << std::string(90, '-') << "\n";

  // data rows: iterate over iterations
  for (size_t iter_idx = 0; iter_idx < all_stats.size(); ++iter_idx) {
    const auto &metrics = all_stats[iter_idx];

    // Iterate over cells
    for (size_t cell_idx = 0; cell_idx < metrics.id.size(); ++cell_idx) {
      // Iterate over species for this cell
      for (size_t species_idx = 0; species_idx < species_names.size();
           ++species_idx) {
        out << std::left << std::setw(15) << metrics.id[cell_idx]
            << std::setw(15) << metrics.iteration << std::setw(15)
            << metrics.rollback_count << std::setw(15)
            << species_names[species_idx] << std::setw(15)
            << metrics.mape[cell_idx][species_idx] << std::setw(15)
            << metrics.rrmse[cell_idx][species_idx] << "\n";
      }
    }
    out << "\n";
  }

  out.close();
  std::cout << "Error metrics written to " << out_dir << "/" << filename
            << "\n";
}
