#ifndef CONTROLMODULE_H_
#define CONTROLMODULE_H_

#include "Base/Macros.hpp"
#include "Chemistry/ChemistryModule.hpp"
#include "poet.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace poet {

class ChemistryModule;

class ControlModule {

public:
  /* Control configuration*/

  // std::uint32_t global_iter = 0;
  // std::uint32_t sur_disabled_counter = 0;
  // std::uint32_t rollback_counter = 0;

  void updateControlIteration(const uint32_t &iter, const bool &dht_enabled,
                              const bool &interp_enaled);

  void initiateWarmupPhase(bool dht_enabled, bool interp_enabled);

  bool checkAndRollback(ChemistryModule &chem, uint32_t &iter);

  struct SpeciesErrorMetrics {
    std::vector<double> mape;
    std::vector<double> rrmse;
    uint32_t iteration; // iterations in simulation after rollbacks
    uint32_t rollback_count;

    SpeciesErrorMetrics(uint32_t species_count, uint32_t iter, uint32_t counter)
        : mape(species_count, 0.0), rrmse(species_count, 0.0), iteration(iter),
          rollback_count(counter) {}
  };

  void computeSpeciesErrorMetrics(const std::vector<double> &reference_values,
                                  const std::vector<double> &surrogate_values,
                                  const uint32_t size_per_prop);

  std::vector<SpeciesErrorMetrics> metricsHistory;

  struct ControlSetup {
    std::string out_dir;
    std::uint32_t checkpoint_interval;
    std::uint32_t control_interval;
    std::vector<std::string> species_names;
    std::vector<double> mape_threshold;
  };

  void enableControlLogic(const ControlSetup &setup) {
    this->out_dir = setup.out_dir;
    this->checkpoint_interval = setup.checkpoint_interval;
    this->control_interval = setup.control_interval;
    this->species_names = setup.species_names;
    this->mape_threshold = setup.mape_threshold;
  }

  bool getControlIntervalEnabled() const {
    return this->control_interval_enabled;
  }

  void applyControlLogic(ChemistryModule &chem, uint32_t &iter);

  void writeCheckpointAndMetrics(ChemistryModule &chem, uint32_t iter);

  auto getGlobalIteration() const noexcept { return global_iteration; }

  void setChemistryModule(poet::ChemistryModule *c) { chem = c; }

  auto getControlInterval() const { return this->control_interval; }

  std::vector<double> getMapeThreshold() const { return this->mape_threshold; }

  /* Profiling getters */

  auto getUpdateCtrlLogicTime() const { return this->prep_t; }

  auto getWriteCheckpointTime() const { return this->w_check_t; }

  auto getReadCheckpointTime() const { return this->r_check_t; }

  auto getWriteMetricsTime() const { return this->stats_t; }

private:
  bool rollback_enabled = false;
  bool control_interval_enabled = false;

  poet::ChemistryModule *chem = nullptr;

  std::uint32_t checkpoint_interval = 0;
  std::uint32_t control_interval = 0;
  std::uint32_t global_iteration = 0;
  std::uint32_t rollback_count = 0;
  std::uint32_t sur_disabled_counter = 0;
  std::vector<double> mape_threshold;

  std::vector<std::string> species_names;
  std::string out_dir;

  double prep_t = 0.;
  double r_check_t = 0.;
  double w_check_t = 0.;
  double stats_t = 0.;

  /* Buffer for shuffled surrogate data */
  std::vector<double> sur_shuffled;
};

} // namespace poet

#endif // CONTROLMODULE_H_