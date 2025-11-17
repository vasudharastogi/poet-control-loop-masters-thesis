#ifndef CONTROLMODULE_H_
#define CONTROLMODULE_H_

#include "Base/Macros.hpp"
#include "Chemistry/ChemistryModule.hpp"
#include "poet.hpp"
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace poet {

class ChemistryModule;

struct ControlConfig {
  uint32_t control_interval = 0;
  uint32_t checkpoint_interval = 0;
  double zero_abs = 0.0;
  std::vector<double> mape_threshold;
};

struct SpeciesErrorMetrics {
  std::vector<double> mape;
  std::vector<double> rrmse;
  uint32_t iteration = 0;
  uint32_t rollback_count = 0;

  SpeciesErrorMetrics(uint32_t n_species, uint32_t iter, uint32_t rb_count)
      : mape(n_species, 0.0), rrmse(n_species, 0.0), iteration(iter),
        rollback_count(rb_count) {}
};

class ControlModule {
public:
  explicit ControlModule(const ControlConfig &config, ChemistryModule *chem);

  void beginIteration(const uint32_t &iter, const bool &dht_enabled,
                      const bool &interp_enaled);

  void writeCheckpoint(uint32_t &iter, const std::string &out_dir);

  void writeErrorMetrics(const std::string &out_dir,
                         const std::vector<std::string> &species);

  std::optional<uint32_t> getRollbackTarget();

  void computeErrorMetrics(const std::vector<double> &reference_values,
                           const std::vector<double> &surrogate_values,
                           const uint32_t size_per_prop,
                           const std::vector<std::string> &species);

  void processCheckpoint(uint32_t &current_iter,
                         const std::string &out_dir,
                         const std::vector<std::string> &species);

  std::optional<uint32_t>
  getRollbackTarget(const std::vector<std::string> &species);

  bool shouldBcastFlags() const;
  bool getControlIntervalEnabled() const {
    return this->control_interval_enabled;
  }

  bool getFlushRequest() const { return flush_request; }
  void clearFlushRequest() { flush_request = false; }

  /* Profiling getters */

  double getUpdateCtrlLogicTime() const { return prep_t; }
  double getWriteCheckpointTime() const { return w_check_t; }
  double getReadCheckpointTime() const { return r_check_t; }
  double getWriteMetricsTime() const { return stats_t; }

private:
  void updateStabilizationPhase(bool dht_enabled, bool interp_enabled);

  void readCheckpoint(uint32_t &current_iter,
                      uint32_t rollback_iter, const std::string &out_dir);

  uint32_t getRollbackIter();

  ControlConfig config;
  ChemistryModule *chem = nullptr;

  std::uint32_t global_iteration = 0;
  std::uint32_t rollback_count = 0;
  std::uint32_t disable_surr_counter = 0;
  std::uint32_t last_checkpoint_written = 0;

  bool rollback_enabled = false;
  bool control_interval_enabled = false;
  bool flush_request = false;

  bool bcast_flags = false;

  std::vector<SpeciesErrorMetrics> metrics_history;

  double prep_t = 0.;
  double r_check_t = 0.;
  double w_check_t = 0.;
  double stats_t = 0.;
};

} // namespace poet

#endif // CONTROLMODULE_H_