#ifndef CONTROLMODULE_H_
#define CONTROLMODULE_H_

#include "Base/Macros.hpp"
#include "Chemistry/ChemistryModule.hpp"
#include "Transport/DiffusionModule.hpp"
#include "IO/HDF5Functions.hpp"
#include "poet.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace poet {

class ChemistryModule;
class DiffusionModule;

struct ControlConfig {
  uint32_t stab_interval = 0;
  uint32_t checkpoint_interval = 0;
  double zero_abs = 0.0;
  std::vector<double> mape_threshold;
};

struct SpeciesErrorMetrics {
  std::vector<std::uint32_t> id;
  std::vector<std::vector<double>> mape;
  std::vector<std::vector<double>> rrmse;
  uint32_t iteration = 0;
  uint32_t rollback_count = 0;

  SpeciesErrorMetrics(uint32_t n_cells, uint32_t n_species, uint32_t iter,
                      uint32_t rb_count)
      : mape(n_cells, std::vector<double>(n_species, 0.0)),
        rrmse(n_cells, std::vector<double>(n_species, 0.0)), iteration(iter),
        rollback_count(rb_count) {}
};



class ControlModule {

public:
  explicit ControlModule(const ControlConfig &config);

  void beginIteration(ChemistryModule &chem, const uint32_t &iter,
                      const bool &dht_enabled, const bool &interp_enaled);

  void writeErrorMetrics(const std::string &out_dir,
                         const std::vector<std::string> &species);

  std::optional<uint32_t> getRollbackTarget();

  void computeErrorMetrics(std::vector<std::vector<double>> &reference_values,
                           std::vector<std::vector<double>> &surrogate_values,
                           const std::vector<std::string> &species);

  void processCheckpoint(DiffusionModule &diffusion, uint32_t &current_iter,
                         const std::string &out_dir,
                         const std::vector<std::string> &species);

  std::optional<uint32_t>
  getRollbackTarget(const std::vector<std::string> &species);

   
  bool shouldBcastFlags();

  bool getFlushRequest() const { return flush_request; }
  void clearFlushRequest() { flush_request = false; }

  auto getGlobalIteration() const noexcept { return global_iteration; }

  // void setChemistryModule(poet::ChemistryModule *c) { chem = c; }

  std::vector<double> getMapeThreshold() const {
    return this->config.mape_threshold;
  }

  std::vector<uint32_t> getCtrlCellIds() const { return this->ctrl_cell_ids; }

  /* Profiling getters */
  auto getUpdateCtrlLogicTime() const { return prep_t; }
  auto getWriteCheckpointTime() const { return w_check_t; }
  auto getReadCheckpointTime() const { return r_check_t; }
  auto getWriteMetricsTime() const { return stats_t; }

private:
  void updateStabilizationPhase(ChemistryModule &chem, bool dht_enabled,
                                bool interp_enabled);

  void readCheckpoint(DiffusionModule &diffusion, uint32_t &current_iter,
                      uint32_t rollback_iter, const std::string &out_dir);
  void writeCheckpoint(DiffusionModule &diffusion, uint32_t &iter,
                       const std::string &out_dir);

  uint32_t getRollbackIter();

  ControlConfig config;

  std::uint32_t global_iteration = 0;
  std::uint32_t rollback_count = 0;
  std::uint32_t disable_surr_counter = 0;
  std::vector<uint32_t> ctrl_cell_ids;
  std::uint32_t last_checkpoint_written = 0;
  std::uint32_t penalty_interval = 0;

  bool rollback_enabled = false;
  bool flush_request = false;
  bool stab_phase_ended = false;

  bool bcast_flags = false;

  std::vector<SpeciesErrorMetrics> metrics_history;

  double prep_t = 0.;
  double r_check_t = 0.;
  double w_check_t = 0.;
  double stats_t = 0.;
};

} // namespace poet

#endif // CONTROLMODULE_H_