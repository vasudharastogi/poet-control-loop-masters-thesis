#include <string>
#include "Control/ControlModule.hpp"

namespace poet
{
    void writeStatsToCSV(const std::vector<ControlModule::SpeciesErrorMetrics> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &out_dir,
                         const std::string &filename);
} // namespace poet
