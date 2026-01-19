#pragma once
#include "Control/ControlModule.hpp"
#include <string>
#include <vector>


int write_metrics(const std::vector<poet::CellMetrics> &metrics_history,
                  const std::vector<std::string> &species_names,
                  const std::string &dir_path, const std::string &file_name);

void writeCellStatsToCSV(const std::vector<poet::CellMetrics> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &out_dir, const std::string &filename);

void writeSpeciesStatsToCSV(const std::vector<poet::SpeciesMetrics> &all_stats,
                            const std::vector<std::string> &species_names,
                            const std::string &out_dir, const std::string &filename,
                            const poet::ControlConfig &config);