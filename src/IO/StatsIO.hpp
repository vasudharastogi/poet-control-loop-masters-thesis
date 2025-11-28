#include "Control/ControlModule.hpp"
#include <string>

namespace poet {
void writeStatsToCSV(const std::vector<SpeciesMetrics> &all_stats, const std::vector<std::string> &species_names,
                     const std::string &out_dir, const std::string &filename);
} // namespace poet
