#ifndef CONTROLMODULE_H_
#define CONTROLMODULE_H_

#include "Base/Macros.hpp"
#include "Chemistry/ChemistryModule.hpp"
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
  uint32_t chkpt_interval = 0;
  uint32_t rb_limit = 0;
  uint32_t rb_aging_limit = 0;
  double zero_abs = 0.0;
  std::vector<double> mape_threshold;
};

struct CellMetrics {
  std::vector<std::uint32_t> id;
  std::vector<std::vector<double>> mape;
  std::vector<std::vector<double>> rrmse;
  std::vector<std::vector<double>> conc;
  uint32_t iteration = 0;
  uint32_t rb_count = 0;

  CellMetrics(uint32_t n_cells, uint32_t n_species, uint32_t iter, uint32_t rb_count)
      : id(n_cells, 0), mape(n_cells, std::vector<double>(n_species, 0.0)),
        rrmse(n_cells, std::vector<double>(n_species, 0.0)),
        conc(n_cells, std::vector<double>(n_species, 0.0)), iteration(iter),
        rb_count(rb_count) {}
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

class ControlModule {

public:
  explicit ControlModule(const ControlConfig &config, ChemistryModule *chem);

  /* store global iteration, dht and pht settings */
  void beginIteration(uint32_t &iter, const bool &dht_enabled, const bool &interp_enaled);

  void writeMetrics(uint32_t &iter, const std::string &out_dir,
                    const std::vector<std::string> &species);

  /* Computes accuracy metrics (MAPE, RRMSE) by comparing reference and surrogate values */
  void computeMetrics(std::vector<std::vector<double>> &reference_values,
                      std::vector<std::vector<double>> &surrogate_values,
                      const std::vector<std::string> &species,
                      const uint32_t size_per_prop, const std::string &out_dir);

  void processCheckpoint(uint32_t &current_iter, const std::string &out_dir,
                         const std::vector<std::string> &species);

  /* Returns rollback target or nullopt and sets flush_request when threshold is exceeded.
  It reports the cells which exceed the threshold. */
  std::optional<uint32_t> findRbTarget(const std::vector<std::string> &species);

  bool getFlushRequest() const { return flush_request; }
  void clearFlushRequest() { flush_request = false; }

  auto getGlobalIteration() const noexcept { return global_iter; }

  std::vector<double> getMapeThreshold() const { return this->config.mape_threshold; }

  std::vector<uint32_t> getCtrlCellIds() const { return this->ctrl_cell_ids; }

  /* not being used currently */
  bool needsFlagBcast() const;

  /* profiling getters */
  auto getCtrlLogicTime() const { return prep_t; }
  auto getChkptWriteTime() const { return w_check_t; }
  auto getChkptReadTime() const { return r_check_t; }
  auto getMetricsWriteTime() const { return stats_t; }

private:
  /* updates DHT and interpolation settings based on the control logic */
  void updateSurrState(bool dht_enabled, bool interp_enabled);

  void readCheckpoint(uint32_t &current_iter, uint32_t rollback_iter,
                      const std::string &out_dir);
  void writeCheckpoint(uint32_t &iter, const std::string &out_dir);

  uint32_t calcRbIter();

  /* computes species-level metrics (MAPE, RRMSE) across all cells*/
  void computeSpeciesMetrics(const std::vector<std::vector<double>> &ref_values,
                             const std::vector<std::vector<double>> &sur_values,
                             const std::vector<std::string> &species,
                             uint32_t size_per_prop, SpeciesMetrics &s_metrics);

  /* computes species-level metrics (MAPE, RRMSE) across all cells */
  void computeCellMetrics(const std::vector<std::vector<double>> &ref_values,
                          const std::vector<std::vector<double>> &sur_values,
                          CellMetrics &c_metrics);

  /* computes alpha using the adapted formula */
  inline double computeAlpha(double ref, double sur) const;

  static void printExceedingCells(const CellMetrics &c_hist, size_t sp_idx,
                                  double sp_thr);

  /* checks if rollback limit has been reached */
  inline bool rbLimitReached() const;

  bool inWarmup() const;

  void setSurrState(const SurrState &state);

  /* tracks stabilization phase via counters */
  void trackStabPhase();
  
  /* tracks active interpolation uptime via counters */
  void trackSurrUptime();

  ControlConfig config;
  ChemistryModule *chem = nullptr;

  std::uint32_t global_iter = 0;
  std::uint32_t rb_count = 0;
  std::uint32_t rb_limit = 0;
  std::uint32_t stab_countdown = 0;
  std::uint32_t surr_active = 0;
  std::vector<uint32_t> ctrl_cell_ids;

  bool rb_enabled = false;
  /* signals a checkpoint needs to be restored*/
  bool flush_request = false;

  std::vector<CellMetrics> c_history;
  std::vector<SpeciesMetrics> s_history;

  double prep_t = 0.;
  double r_check_t = 0.;
  double w_check_t = 0.;
  double stats_t = 0.;
};

} // namespace poet

#endif // CONTROLMODULE_H_