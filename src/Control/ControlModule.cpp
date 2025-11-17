#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

poet::ControlModule::ControlModule(const ControlConfig &config_,
                                   ChemistryModule *chem_)
    : config(config_), chem(chem_) {
  assert(chem && "ChemistryModule pointer must not be null");
}

void poet::ControlModule::beginIteration(const uint32_t &iter,
                                         const bool &dht_enabled,
                                         const bool &interp_enabled) {

  /* dht_enabled and inter_enabled are user settings set before startig the
   * simulation*/
  double prep_a, prep_b;

  prep_a = MPI_Wtime();
  if (config.control_interval == 0) {
    control_interval_enabled = false;
    return;
  }
  global_iteration = iter;

  updateStabilizationPhase(dht_enabled, interp_enabled);

  control_interval_enabled =
      (config.control_interval > 0 && (iter % config.control_interval == 0));

  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

void poet::ControlModule::updateStabilizationPhase(bool dht_enabled,
                                                   bool interp_enabled) {
  if (rollback_enabled) {
    if (disable_surr_counter > 0) {
      --disable_surr_counter;
      MSG("Rollback counter: " + std::to_string(disable_surr_counter));
    } else {
      rollback_enabled = false;
    }
  }
  // user requested DHT/INTEP? keep them disabled but enable warmup-phase so
  if (global_iteration <= config.control_interval || rollback_enabled) {
    chem->SetStabEnabled(true);
    chem->SetDhtEnabled(false);
    chem->SetInterpEnabled(false);
    return;
  }
  chem->SetStabEnabled(false);
  chem->SetDhtEnabled(dht_enabled);
  chem->SetInterpEnabled(interp_enabled);
}

void poet::ControlModule::writeCheckpoint(uint32_t &iter,
                                          const std::string &out_dir) {
  double w_check_a, w_check_b;
  w_check_a = MPI_Wtime();
  write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                   {.field = chem->getField(), .iteration = iter});
  w_check_b = MPI_Wtime();
  this->w_check_t += w_check_b - w_check_a;

  last_checkpoint_written = iter;
}

void poet::ControlModule::readCheckpoint(uint32_t &current_iter,
                                         uint32_t rollback_iter,
                                         const std::string &out_dir) {
  double r_check_a, r_check_b;
  r_check_a = MPI_Wtime();
  Checkpoint_s checkpoint_read{.field = chem->getField()};
  read_checkpoint(out_dir,
                  "checkpoint" + std::to_string(rollback_iter) + ".hdf5",
                  checkpoint_read);
  current_iter = checkpoint_read.iteration;
  r_check_b = MPI_Wtime();
  r_check_t += r_check_b - r_check_a;
}

void poet::ControlModule::writeErrorMetrics(
    const std::string &out_dir, const std::vector<std::string> &species) {
  double stats_a, stats_b;

  stats_a = MPI_Wtime();
  writeStatsToCSV(metrics_history, species, out_dir, "metrics_overview");
  stats_b = MPI_Wtime();

  this->stats_t += stats_b - stats_a;
}

uint32_t poet::ControlModule::getRollbackIter() {

  uint32_t last_iter = ((global_iteration - 1) / config.checkpoint_interval) *
                       config.checkpoint_interval;

  uint32_t rollback_iter = (last_iter <= last_checkpoint_written)
                               ? last_iter
                               : last_checkpoint_written;
  return rollback_iter;
}

std::optional<uint32_t> poet::ControlModule::getRollbackTarget(
    const std::vector<std::string> &species) {
  double r_check_a, r_check_b;

  if (metrics_history.empty()) {
    MSG("No error history yet, skipping rollback check.");
    flush_request = false;
    return std::nullopt;
  }
  if (rollback_count > 3) {
    MSG("Rollback limit reached, skipping rollback.");
    flush_request = false;
    return std::nullopt;
  }

  const auto &mape = metrics_history.back().mape;
  for (uint32_t i = 0; i < species.size(); ++i) {
    // skip Charge
    if (i == 4 || mape[i] == 0) {
      continue;
    }
    if (mape[i] > config.mape_threshold[i]) {
      if (last_checkpoint_written == 0) {
        MSG(" Threshold exceeded but no checkpoint exists yet.");
        return std::nullopt;
      }

      flush_request = true;

      MSG("T hreshold exceeded " + species[i] +
          " has MAPE = " + std::to_string(mape[i]) +
          " exceeding threshold = " + std::to_string(config.mape_threshold[i]));
      return getRollbackIter();
    }
  }
  MSG("All species are within their MAPE thresholds.");
  flush_request = false;
  return std::nullopt;
}

void poet::ControlModule::computeErrorMetrics(
    const std::vector<double> &reference_values,
    const std::vector<double> &surrogate_values, const uint32_t size_per_prop,
    const std::vector<std::string> &species) {

  SpeciesErrorMetrics metrics(species.size(), global_iteration, rollback_count);

  for (uint32_t i = 0; i < species.size(); ++i) {
    double err_sum = 0.0;
    double sqr_err_sum = 0.0;
    uint32_t base_idx = i * size_per_prop;

    for (uint32_t j = 0; j < size_per_prop; ++j) {
      const double ref_value = reference_values[base_idx + j];
      const double sur_value = surrogate_values[base_idx + j];
      const double ZERO_ABS = config.zero_abs;

      if (std::isnan(ref_value) || std::isnan(sur_value)) {
        continue;
      }
      if (std::abs(ref_value) < ZERO_ABS) {
        if (std::abs(sur_value) >= ZERO_ABS) {
          err_sum += 1.0;
          sqr_err_sum += 1.0;
        }
      }
      // Both zero: skip
      else {
        double alpha = 1.0 - (sur_value / ref_value);
        err_sum += std::abs(alpha);
        sqr_err_sum += alpha * alpha;
      }
    }
    metrics.mape[i] = 100.0 * (err_sum / size_per_prop);
    metrics.rrmse[i] = std::sqrt(sqr_err_sum / size_per_prop);
  }
  metrics_history.push_back(metrics);
}

void poet::ControlModule::processCheckpoint(
    uint32_t &current_iter, const std::string &out_dir,
    const std::vector<std::string> &species) {

  if (!control_interval_enabled)
    return;

  if (flush_request && rollback_count < 3) {
    uint32_t target = getRollbackIter();
    readCheckpoint(current_iter, target, out_dir);

    rollback_enabled = true;
    rollback_count++;
    disable_surr_counter = config.control_interval;

    MSG("Restored checkpoint " + std::to_string(target) +
        ", surrogates disabled for " + std::to_string(config.control_interval));
  } else {
    writeCheckpoint(global_iteration, out_dir);
  }
}

bool poet::ControlModule::shouldBcastFlags() const {
  if (global_iteration == 1 ||
      global_iteration % config.control_interval == 1) {
    return true;
  }
  return false;
}
