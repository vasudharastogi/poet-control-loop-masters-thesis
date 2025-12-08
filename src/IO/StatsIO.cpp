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
#include <limits>
#include <string>
#include <unordered_map>
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

  const std::size_t n_species = species_names.size();
  if (n_species == 0) {
    std::cerr << "Species names list is empty; skipping HDF5 metrics dump.\n";
    return -1;
  }

  // Always start from a clean file: no stale datasets, no mismatched shapes.
  HighFive::File file(file_path.string(), HighFive::File::Truncate);

  // cell_id -> all rows over time
  using Row = std::vector<double>;
  using Matrix = std::vector<Row>;
  std::unordered_map<std::uint32_t, Matrix> mape_per_cell;
  std::unordered_map<std::uint32_t, Matrix> rrmse_per_cell;
  std::unordered_map<std::uint32_t, Matrix> conc_per_cell;

  // Build per-cell time series from the *entire* history.
  for (const auto &metrics : metrics_history) {
    const std::uint32_t iter = metrics.iteration;

    for (std::size_t cell_idx = 0; cell_idx < metrics.id.size(); ++cell_idx) {
      const auto cell_id = metrics.id[cell_idx];

      // --- MAPE row: [iteration, mape...]
      Row mape_row;
      mape_row.reserve(1 + n_species);
      mape_row.push_back(static_cast<double>(iter));

      const auto &src_m = metrics.mape[cell_idx];
      const std::size_t len_m = std::min(src_m.size(), n_species);
      mape_row.insert(mape_row.end(), src_m.begin(), src_m.begin() + len_m);
      if (len_m < n_species) {
        mape_row.resize(1 + n_species, std::numeric_limits<double>::quiet_NaN());
      }

      // --- RRMSE row: [iteration, rrmse...]
      Row rrmse_row;
      rrmse_row.reserve(1 + n_species);
      rrmse_row.push_back(static_cast<double>(iter));

      const auto &src_r = metrics.rrmse[cell_idx];
      const std::size_t len_r = std::min(src_r.size(), n_species);
      rrmse_row.insert(rrmse_row.end(), src_r.begin(), src_r.begin() + len_r);
      if (len_r < n_species) {
        rrmse_row.resize(1 + n_species, std::numeric_limits<double>::quiet_NaN());
      }

      Row conc_row;
      conc_row.reserve(1 + n_species);
      conc_row.push_back(static_cast<double>(iter));

      const auto &src_c = metrics.conc[cell_idx];
      const std::size_t len_c = std::min(src_c.size(), n_species);
      conc_row.insert(conc_row.end(), src_c.begin(), src_c.begin() + len_c);
      if (len_c < n_species) {
        conc_row.resize(1 + n_species, std::numeric_limits<double>::quiet_NaN());
      }

      // Append to per-cell matrices
      mape_per_cell[cell_id].push_back(std::move(mape_row));
      rrmse_per_cell[cell_id].push_back(std::move(rrmse_row));
      conc_per_cell[cell_id].push_back(std::move(conc_row));
    }
  }

  // Now dump each cell’s full time series
  for (const auto &kv : mape_per_cell) {
    const auto cell_id = kv.first;
    const auto &mape_matrix = kv.second;
    const auto &rrmse_matrix = rrmse_per_cell[cell_id];
    const auto &conc_matrix = conc_per_cell[cell_id];

    const std::string cell_grp = "cells/" + std::to_string(cell_id);

    H5Easy::dump(file, cell_grp + "/mape", mape_matrix, H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, cell_grp + "/rrmse", rrmse_matrix, H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, cell_grp + "/conc", conc_matrix, H5Easy::DumpMode::Overwrite);
  }

  return 0;
}

void writeCellStatsToCSV(const std::vector<poet::CellMetrics> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &out_dir, const std::string &filename) {
  std::ofstream out(std::filesystem::path(out_dir) / filename);
  if (!out.is_open()) {
    std::cerr << "Could not open " << filename << " !" << std::endl;
    return;
  }

  // Header
  out << std::left << std::setw(15) << "CellID" << std::setw(15) << "Iteration"
      << std::setw(15) << "Rollback" << std::setw(15) << "Species" << std::setw(15)
      << "MAPE" << std::setw(15) << "RRMSE"
      << "\n"
      << std::string(90, '-') << "\n";

  // Data rows (fix column ordering: include rb_count before Species)
  for (const auto &metrics : all_stats) {
    for (size_t cell_idx = 0; cell_idx < metrics.id.size(); ++cell_idx) {
      for (size_t sp_idx = 0; sp_idx < species_names.size(); ++sp_idx) {
        out << std::left << std::setw(15) << metrics.id[cell_idx] << std::setw(15)
            << metrics.iteration << std::setw(15) << metrics.rb_count << std::setw(15)
            << species_names[sp_idx] << std::setw(15) << metrics.mape[cell_idx][sp_idx]
            << std::setw(15) << metrics.rrmse[cell_idx][sp_idx] << "\n";
      }
    }
    out << "\n";
  }
  out.close();
  std::cout << "Cell error metrics written to " << out_dir << "/" << filename << "\n";
}

void writeSpeciesStatsToCSV(const std::vector<poet::SpeciesMetrics> &all_stats,
                            const std::vector<std::string> &species_names,
                            const std::string &out_dir, const std::string &filename) {
  std::ofstream out(std::filesystem::path(out_dir) / filename);
  if (!out.is_open()) {
    std::cerr << "Could not open " << filename << " !" << std::endl;
    return;
  }

  // Header
  out << std::left << std::setw(15) << "Iteration" << std::setw(15) << "Rollback"
      << std::setw(15) << "Species" << std::setw(15) << "MAPE" << std::setw(15) << "RRMSE"
      << "\n"
      << std::string(75, '-') << "\n";

  // Data rows
  for (const auto &metrics : all_stats) {
    for (size_t sp_idx = 0; sp_idx < species_names.size(); ++sp_idx) {
      out << std::left << std::setw(15) << metrics.iteration << std::setw(15)
          << metrics.rb_count << std::setw(15) << species_names[sp_idx] << std::setw(15)
          << metrics.mape[sp_idx] << std::setw(15) << metrics.rrmse[sp_idx] << "\n";
    }
    out << "\n";
  }
  out.close();
  std::cout << "Species error metrics written to " << out_dir << "/" << filename << "\n";
}
