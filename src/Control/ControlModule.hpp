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
  uint32_t ctrl_interval = 0;
  uint32_t chkpt_interval = 0;
  uint32_t rb_limit = 0;
  uint32_t rb_aging_limit = 0;
  double zero_abs = 0.0;
  std::vector<double> mape_threshold;
};

struct SpeciesMetrics {
  std::vector<double> mape;
  std::vector<double> rrmse;
  uint32_t iteration = 0;
  uint32_t rb_count = 0;

  SpeciesMetrics(uint32_t n_species, uint32_t iter, uint32_t count)
      : mape(n_species, 0.0), rrmse(n_species, 0.0), iteration(iter), rb_count(count) {}
};

struct SurrState {
  bool stab_enabled;
  bool dht_enabled;
  bool interp_enabled;
};

enum { WARMUP, SUR_ENABLED, STAB, SUR_DISABLED };

class ControlModule {
public:
  explicit ControlModule(const ControlConfig &config, ChemistryModule *chem);

  void beginIteration(const uint32_t &iter, const bool &dht_enabled,
                      const bool &interp_enaled);

  void writeCheckpoint(uint32_t &iter, const std::string &out_dir);

  void writeMetrics(const std::string &out_dir, const std::vector<std::string> &species);

  std::optional<uint32_t> findRbTarget();

  void computeMetrics(const std::vector<double> &reference_values,
                      const std::vector<double> &surrogate_values,
                      const uint32_t size_per_prop,
                      const std::vector<std::string> &species);

  void processCheckpoint(uint32_t &current_iter, const std::string &out_dir,
                         const std::vector<std::string> &species);

  std::optional<uint32_t> findRbTarget(const std::vector<std::string> &species);

  bool needsFlagBcast() const;
  bool isCtrlIntervalActive() const { return this->ctrl_active; }

  bool getFlushRequest() const { return flush_request; }
  void clearFlushRequest() { flush_request = false; }

  /* Profiling getters */

  double getCtrlLogicTime() const { return prep_t; }
  double getChkptWriteTime() const { return w_check_t; }
  double getChkptReadTime() const { return r_check_t; }
  double getMetricsWriteTime() const { return stats_t; }

private:
  void updateSurrState(bool dht_enabled, bool interp_enabled);

  void readCheckpoint(uint32_t &current_iter, uint32_t rollback_iter,
                      const std::string &out_dir);

  uint32_t calcRbIter();

  inline bool rbLimitReached() const;

  bool inWarmup() const;

  SurrState getCurrPhase(bool dht_enabled, bool interp_enabled);
  void setSurrState(const SurrState &state);

  void trackStabPhase();

  void trackSurrUptime();

  void triggerRb(uint32_t &curr_iter, const std::string &out_dir);

  ControlConfig config;
  ChemistryModule *chem = nullptr;

  std::uint32_t global_iter = 0;
  std::uint32_t rb_count = 0;
  std::uint32_t stab_countdown = 0;
  std::uint32_t surr_active = 0;
  std::uint32_t last_chkpt_written = 0;

  bool rb_enabled = false;
  bool ctrl_active = false;
  bool flush_request = false;

  std::vector<SpeciesMetrics> s_history;

  double prep_t = 0.;
  double r_check_t = 0.;
  double w_check_t = 0.;
  double stats_t = 0.;
};

} // namespace poet

#endif // CONTROLMODULE_H_