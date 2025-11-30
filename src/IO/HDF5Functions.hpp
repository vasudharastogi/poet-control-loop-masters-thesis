#pragma once

#include "Datatypes.hpp"
#include <string>
#include <vector>

namespace poet {
  struct SpeciesMetrics;
  struct CellMetrics;
}

int write_checkpoint(const std::string &dir_path, const std::string &file_name,
                     struct Checkpoint_s &&checkpoint);

int read_checkpoint(const std::string &dir_path, const std::string &file_name,
                    struct Checkpoint_s &checkpoint);

int write_metrics(const std::vector<poet::CellMetrics> &metrics_history,
                  const std::vector<std::string> &species_names,
                  const std::string &dir_path, const std::string &file_name);