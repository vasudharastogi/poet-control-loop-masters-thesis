#include "IO/StatsIO.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace poet {
void writeStatsToCSV(const std::vector<SpeciesMetrics> &all_stats,
                     const std::vector<std::string> &species_names,
                     const std::string &out_dir, const std::string &filename,
                     const ControlConfig &config) {
  std::filesystem::path full_path = std::filesystem::path(out_dir) / filename;

  std::ofstream out(full_path);
  if (!out.is_open()) {
    std::cerr << "Could not open " << filename << " !" << std::endl;
    return;
  }

  out << "Configuration Parameters:\n";
  out << "ctrl_interval: " << config.ctrl_interval << std::endl;
  out << "rb_limit: " << config.rb_limit << std::endl;
  out << "mape_threshold: " << config.mape_threshold[0] << std::endl;

  out << "\n\n";

  // header
  out << std::left << std::setw(15) << "Iteration" << std::setw(15) << "Rollback"
      << std::setw(15) << "Species" << std::setw(15) << "MAPE" << std::setw(15) << "RRSME"
      << "\n";

  out << std::string(75, '-') << "\n";

  // data rows
  for (size_t i = 0; i < all_stats.size(); ++i) {
    for (size_t j = 0; j < species_names.size(); ++j) {
      out << std::left << std::setw(15) << all_stats[i].iteration << std::setw(15)
          << all_stats[i].rb_count << std::setw(15) << species_names[j] << std::setw(15)
          << all_stats[i].mape[j] << std::setw(15) << all_stats[i].rrmse[j] << "\n";
    }
    out << "\n";
  }

  out.close();
  std::cout << "Error metrics written to " << out_dir << "/" << filename << "\n";
}
} // namespace poet
  // namespace poet