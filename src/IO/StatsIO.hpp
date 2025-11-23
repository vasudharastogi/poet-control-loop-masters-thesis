#include <string>
#include <vector>

void writeSpeciesStatsToCSV(
    const std::vector<poet::SpeciesErrorMetrics> &all_stats,
    const std::vector<std::string> &species_names, const std::string &out_dir,
    const std::string &filename);
    
void writeCellStatsToCSV(const std::vector<poet::CellErrorMetrics> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &out_dir,
                         const std::string &filename);
// namespace poet
