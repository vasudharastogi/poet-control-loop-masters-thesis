#include "IO/StatsIO.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip> 
#include <filesystem>

namespace poet
{
    void writeStatsToCSV(const std::vector<ControlModule::SpeciesErrorMetrics> &all_stats,
                         const std::vector<std::string> &species_names,
                         const std::string &out_dir,
                         const std::string &filename)
    {
        std::filesystem::path full_path = std::filesystem::path(out_dir) / filename;

        std::ofstream out(full_path);
        if (!out.is_open())
        {
            std::cerr << "Could not open " << filename << " !" << std::endl;
            return;
        }

        // header: CellID, Iteration, Rollback, Species, MAPE, RRMSE
        out << std::left << std::setw(15) << "CellID"
            << std::setw(15) << "Iteration"
            << std::setw(15) << "Rollback"
            << std::setw(15) << "Species"
            << std::setw(15) << "MAPE"
            << std::setw(15) << "RRMSE" << "\n";

        out << std::string(90, '-') << "\n"; 

        // data rows: iterate over iterations
        for (size_t iter_idx = 0; iter_idx < all_stats.size(); ++iter_idx)
        {
            const auto &metrics = all_stats[iter_idx];
            
            // Iterate over cells
            for (size_t cell_idx = 0; cell_idx < metrics.id.size(); ++cell_idx)
            {
                // Iterate over species for this cell
                for (size_t species_idx = 0; species_idx < species_names.size(); ++species_idx)
                {
                    out << std::left
                        << std::setw(15) << metrics.id[cell_idx]
                        << std::setw(15) << metrics.iteration
                        << std::setw(15) << metrics.rollback_count
                        << std::setw(15) << species_names[species_idx]
                        << std::setw(15) << metrics.mape[cell_idx][species_idx]
                        << std::setw(15) << metrics.rrmse[cell_idx][species_idx]
                        << "\n";
                }
            }
            out << "\n"; 
        }

        out.close();
        std::cout << "Error metrics written to " << out_dir << "/" << filename << "\n";
    }
}
 // namespace poet