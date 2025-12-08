#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>
#include <numeric>

poet::ControlModule::ControlModule(const ControlConfig &config_, ChemistryModule &chem)
    : config(config_), chem_(chem) {}

void poet::ControlModule::beginIteration(uint32_t &iter, const bool &dht_enabled,
                                         const bool &interp_enabled) {
  global_iter = iter;
  double prep_a, prep_b;
  prep_a = MPI_Wtime();
  updateSurrState(dht_enabled, interp_enabled);
  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

void poet::ControlModule::updateSurrState(bool dht_enabled, bool interp_enabled) {
  bool in_warmup = (global_iter <= config.stab_interval);
  bool rb_limit_reached = rbLimitReached();

  if (rb_enabled && stab_countdown > 0) {
    --stab_countdown;
    std::cout << "Rollback counter: " << std::to_string(stab_countdown) << std::endl;
    if (stab_countdown == 0) {
      rb_enabled = false;
    }
  }
  /* disables surrogate during warmup, active rollback or after rollback limit is reached
   */
  if (in_warmup || rb_enabled || rb_limit_reached) {
    chem_.SetStabEnabled(!rb_limit_reached);
    chem_.SetDhtEnabled(false);
    chem_.SetInterpEnabled(false);

    if (rb_limit_reached) {
      std::cout << "Interpolation completely disabled." << std::endl;
    } else {
      std::cout << "In stabilization phase." << std::endl;
    }
    return;
  }

  if (rb_count > 0 && !rb_enabled && !in_warmup) {
    surr_active++;
    if (surr_active > config.rb_interval_limit) {
      surr_active = 0;
      rb_count -= 1;
      std::cout << "Rollback count reset to: " << rb_count << "." << std::endl;
    }
  }

  /* enable user-requested surrogates */
  chem_.SetStabEnabled(false);
  chem_.SetDhtEnabled(dht_enabled);
  chem_.SetInterpEnabled(interp_enabled);
  std::cout << "Interpolating." << std::endl;
}

void poet::ControlModule::writeCheckpoint(uint32_t &iter, const std::string &out_dir) {
  if (global_iter % config.chkpt_interval == 0) {
    double w_check_a = MPI_Wtime();
    write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                     {.field = chem_.getField(), .iteration = iter});
    double w_check_b = MPI_Wtime();
    this->w_check_t += w_check_b - w_check_a;
  }
}

void poet::ControlModule::readCheckpoint(uint32_t &current_iter, uint32_t rollback_iter,
                                         const std::string &out_dir) {
  double r_check_a, r_check_b;

  r_check_a = MPI_Wtime();
  Checkpoint_s checkpoint_read{.field = chem_.getField()};
  read_checkpoint(out_dir, "checkpoint" + std::to_string(rollback_iter) + ".hdf5",
                  checkpoint_read);
  current_iter = checkpoint_read.iteration;
  r_check_b = MPI_Wtime();
  r_check_t += r_check_b - r_check_a;
}

void poet::ControlModule::writeMetrics(uint32_t &iter, const std::string &out_dir,
                                       const std::vector<std::string> &species) {
  double stats_a = MPI_Wtime();
  writeSpeciesStatsToCSV(s_history, species, out_dir, "species_overview.csv");
  if (global_iter % (config.chkpt_interval / 2) == 0) {
    write_metrics(c_history, species, out_dir, "metrics_overview.hdf5");
  }
  double stats_b = MPI_Wtime();
  this->stats_t += stats_b - stats_a;
}

uint32_t poet::ControlModule::calcRbIter() {
  uint32_t last_iter =
      ((global_iter - 1) / config.chkpt_interval) * config.chkpt_interval;
  return last_iter;
}

std::optional<uint32_t>
poet::ControlModule::findRbTarget(const std::vector<std::string> &species) {

  /* Skip threshold checking if already in stabilization phase*/
  if (s_history.empty() || rb_enabled) {
    return std::nullopt;
  }
  const auto &s_hist = s_history.back();

  /* skipping cell_id and id */
  for (size_t sp_idx = 2; sp_idx < species.size(); sp_idx++) {
    /* skip Charge */
    if (sp_idx == 4) {
      continue;
    }
    if (sp_idx >= config.mape_threshold.size()) {
      std::cerr << "Warning: No threshold defined for species " << species[sp_idx]
                << " at index " << std::to_string(sp_idx) << std::endl;
      continue;
    }
    const double sp_mape = s_hist.mape[sp_idx];
    const double sp_thr = config.mape_threshold[sp_idx];

    if (sp_mape > sp_thr) {
      flush_request = true;
      std::cout << "Species " << species[sp_idx]
                << " exceeded threshold (MAPE=" << sp_mape << " > thr=" << sp_thr
                << "). Offending cells:" << std::endl;

      if (!c_history.empty()) {
        const auto &c_hist = c_history.back();
        printExceedingCells(c_hist, sp_idx, sp_thr);
      }
      return calcRbIter();
    }
  }
  return std::nullopt;
}

void poet::ControlModule::printExceedingCells(const CellMetrics &c_hist, size_t sp_idx,
                                              double sp_thr) {
  for (size_t cell_idx = 0; cell_idx < c_hist.mape.size(); ++cell_idx) {
    const double mape = c_hist.mape[cell_idx][sp_idx];
    if (mape > sp_thr) {
      const uint32_t id = c_hist.id[cell_idx];
      std::cout << "  cell_id=" << id << " mape=" << mape << std::endl;
    }
  }
}

void poet::ControlModule::computeMetrics(std::vector<std::vector<double>> &ref_values,
                                         std::vector<std::vector<double>> &sur_values,
                                         const std::vector<std::string> &species,
                                         const uint32_t size_per_prop,
                                         const std::string &out_dir) {

  std::cout << "DEBUG: computeMetrics called at iteration " << global_iter
            << ", rb_count=" << rb_count << "/" << config.rb_limit << std::endl;

  const uint32_t n_cells = ref_values.size();
  const uint32_t n_species = species.size();

  CellMetrics c_metrics(n_cells, n_species, global_iter, rb_count);
  SpeciesMetrics s_metrics(n_species, global_iter, rb_count);

  /* compute cell-level metrics */
  computeCellMetrics(ref_values, sur_values, c_metrics);

  /* compute species-level metrics */
  computeSpeciesMetrics(ref_values, sur_values, species, size_per_prop, s_metrics);

  s_history.push_back(s_metrics);
  c_history.push_back(c_metrics);

  writeMetrics(global_iter, out_dir, species);
}

inline double poet::ControlModule::computeAlpha(double ref, double sur) const {
  const double zero_abs = config.zero_abs;
  if (std::isnan(ref) || std::isnan(sur)) {
    return 0.0;
  }
  if (std::abs(ref) < zero_abs) {
    return (std::abs(sur) < zero_abs) ? 0.0 : 1.0;
  }
  return 1.0 - (sur / ref);
}

void poet::ControlModule::computeSpeciesMetrics(
    const std::vector<std::vector<double>> &ref_values,
    const std::vector<std::vector<double>> &sur_values,
    const std::vector<std::string> &species, uint32_t size_per_prop,
    SpeciesMetrics &s_metrics) {

  const size_t n_cells = ref_values.size();
  const size_t n_species = species.size();
  std::vector<double> err_sum(n_species, 0.0);
  std::vector<double> sqr_sum(n_species, 0.0);

  for (size_t cell_idx = 0; cell_idx < n_cells; cell_idx++) {
    for (size_t sp_idx = 2; sp_idx < n_species; sp_idx++) {
      const double alpha =
          computeAlpha(ref_values[cell_idx][sp_idx], sur_values[cell_idx][sp_idx]);
      err_sum[sp_idx] += std::abs(alpha);
      sqr_sum[sp_idx] += alpha * alpha;
    }
  }

  for (size_t sp_idx = 2; sp_idx < n_species; sp_idx++) {
    s_metrics.mape[sp_idx] =
        100.0 * (err_sum[sp_idx] / static_cast<double>(size_per_prop));
    s_metrics.rrmse[sp_idx] =
        std::sqrt(sqr_sum[sp_idx] / static_cast<double>(size_per_prop));
  }
}

void poet::ControlModule::computeCellMetrics(
    const std::vector<std::vector<double>> &ref_values,
    const std::vector<std::vector<double>> &sur_values, CellMetrics &c_metrics) {

  const size_t n_cells = ref_values.size();
  // Use ref_values to get n_species instead of accessing potentially uninitialized mape
  const size_t n_species = (n_cells > 0) ? ref_values[0].size() : 0;

  if (n_cells == 0 || n_species == 0) {
    std::cerr << "Warning: computeCellMetrics called with empty data" << std::endl;
    return;
  }

  for (size_t cell_idx = 0; cell_idx < n_cells; ++cell_idx) {
    // Assign the per-cell id correctly
    c_metrics.id[cell_idx] = static_cast<uint32_t>(ref_values[cell_idx][0]);

    for (size_t sp_idx = 2; sp_idx < n_species; ++sp_idx) {

      const double alpha =
          computeAlpha(ref_values[cell_idx][sp_idx], sur_values[cell_idx][sp_idx]);
      c_metrics.mape[cell_idx][sp_idx] = 100.0 * std::abs(alpha);
      c_metrics.rrmse[cell_idx][sp_idx] = alpha * alpha;
      c_metrics.conc[cell_idx][sp_idx] = ref_values[cell_idx][sp_idx];
    }
  }

  /* cell_ID and ID columns */
}

void poet::ControlModule::processCheckpoint(uint32_t &current_iter,
                                            const std::string &out_dir,
                                            const std::vector<std::string> &species) {
  if (rbLimitReached()) {
    return;
  }
  if (flush_request) {
    uint32_t target = calcRbIter();
    readCheckpoint(current_iter, target, out_dir);

    rb_enabled = true;
    rb_count++;
    stab_countdown = config.stab_interval;
    flush_request = false;

    std::cout << "Restored checkpoint " << std::to_string(target)
              << ", surrogate disabled for " << std::to_string(config.stab_interval)
              << std::endl;
  } else {
    writeCheckpoint(global_iter, out_dir);
  }
}

bool poet::ControlModule::needsFlagBcast() const {
  return (config.rb_limit > 0) && !rbLimitReached();
}

inline bool poet::ControlModule::rbLimitReached() const {
  /* rollback is completly disabled */
  if (config.rb_limit == 0)
    return false;
  return rb_count >= config.rb_limit;
}