#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

poet::ControlModule::ControlModule(const ControlConfig &config_)
    : config(config_) {}

void poet::ControlModule::beginIteration(ChemistryModule &chem,
                                         const uint32_t &iter,
                                         const bool &dht_enabled,
                                         const bool &interp_enabled) {

  /* dht_enabled and inter_enabled are user settings set before startig the
   * simulation*/
  double prep_a, prep_b;

  prep_a = MPI_Wtime();

  global_iteration = iter;
  updateStabilizationPhase(chem, dht_enabled, interp_enabled);
  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

void poet::ControlModule::updateStabilizationPhase(ChemistryModule &chem,
                                                   bool dht_enabled,
                                                   bool interp_enabled) {
  if (rollback_enabled) {
    if (disable_surr_counter > 0) {
      --disable_surr_counter;
      flush_request = false;
      MSG("Rollback counter: " + std::to_string(disable_surr_counter));
    } else {
      rollback_enabled = false;
    }
  }
  
  bool prev_stab_state = chem.GetStabEnabled();
  
  // user requested DHT/INTEP? keep them disabled but enable warmup-phase so
  if (global_iteration <= config.stab_interval || rollback_enabled) {
    chem.SetStabEnabled(true);
    chem.SetDhtEnabled(false);
    chem.SetInterpEnabled(false);
    return;
  }
  
  chem.SetStabEnabled(false);
  chem.SetDhtEnabled(dht_enabled);
  chem.SetInterpEnabled(interp_enabled);
  
  // Mark that we need to broadcast flags if stab phase just ended
  if (prev_stab_state && !chem.GetStabEnabled()) {
    stab_phase_ended = true;
  }
}

void poet::ControlModule::writeCheckpoint(DiffusionModule &diffusion,
                                          uint32_t &iter,
                                          const std::string &out_dir) {
  double w_check_a, w_check_b;

  w_check_a = MPI_Wtime();
  write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                   {.field = diffusion.getField(), .iteration = iter});
  w_check_b = MPI_Wtime();
  this->w_check_t += w_check_b - w_check_a;

  last_checkpoint_written = iter;
}

void poet::ControlModule::readCheckpoint(DiffusionModule &diffusion,
                                         uint32_t &current_iter,
                                         uint32_t rollback_iter,
                                         const std::string &out_dir) {
  double r_check_a, r_check_b;

  r_check_a = MPI_Wtime();
  Checkpoint_s checkpoint_read{.field = diffusion.getField()};
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
    rollback_enabled = false;
    // flush_request = false;
    return std::nullopt;
  }

  const auto &mape = metrics_history.back().mape;
  for (size_t row = 0; row < mape.size(); row++) {
    for (size_t col = 0; col < species.size() && col < mape[row].size();
         col++) {

      if (mape[row][col] == 0) {
        continue;
      }

      if (mape[row][col] > config.mape_threshold[col]) {

        if (last_checkpoint_written == 0) {
          MSG(" Threshold exceeded but no checkpoint exists yet.");
          return std::nullopt;
        }
        rollback_enabled = true;
        flush_request = true;

        MSG("Threshold exceeded " + species[col] + " has MAPE = " +
            std::to_string(mape[row][col]) + " exceeding threshold = " +
            std::to_string(config.mape_threshold[col]));

        return getRollbackIter();
      }
    }
  }
  rollback_enabled = false;
  // flush_request = false;
  return std::nullopt;
}

void poet::ControlModule::computeErrorMetrics(
    std::vector<std::vector<double>> &reference_values,
    std::vector<std::vector<double>> &surrogate_values,
    const std::vector<std::string> &species) {

  const uint32_t n_cells = reference_values.size();

  SpeciesErrorMetrics metrics(n_cells, species.size(), global_iteration,
                              rollback_count);

  for (size_t cell_i = 0; cell_i < n_cells; cell_i++) {

    metrics.id.push_back(reference_values[cell_i][0]);

    for (size_t sp_i = 0; sp_i < species.size(); sp_i++) {
      const double ref_value = reference_values[cell_i][sp_i];
      const double sur_value = surrogate_values[cell_i][sp_i];
      const double ZERO_ABS = config.zero_abs;

      if (std::isnan(ref_value) || std::isnan(sur_value)) {
        metrics.mape[cell_i][sp_i] = 0.0;
        metrics.rrmse[cell_i][sp_i] = 0.0;
        continue;
      }

      if (std::abs(ref_value) < ZERO_ABS) {
        if (std::abs(sur_value) >= ZERO_ABS) {
          metrics.mape[cell_i][sp_i] = 1.0;
          metrics.rrmse[cell_i][sp_i] = 1.0;
        } else {
          metrics.mape[cell_i][sp_i] = 0.0;
          metrics.rrmse[cell_i][sp_i] = 0.0;
        }
      } else {
        double alpha = 1.0 - (sur_value / ref_value);
        metrics.mape[cell_i][sp_i] = 100.0 * std::abs(alpha);
        metrics.rrmse[cell_i][sp_i] = alpha * alpha;
      }
    }
  }

  std::cout << "[DEBUG] metrics.id.size()=" << metrics.id.size() << std::endl;
  metrics_history.push_back(metrics);
  std::cout << "[DEBUG] metricsHistory.size()=" << metrics_history.size()
            << std::endl;
}

void poet::ControlModule::processCheckpoint(
    DiffusionModule &diffusion, uint32_t &current_iter,
    const std::string &out_dir, const std::vector<std::string> &species) {

  if (flush_request) {
    uint32_t target = getRollbackIter();
    readCheckpoint(diffusion, current_iter, target, out_dir);

    rollback_enabled = true;
    rollback_count++;
    disable_surr_counter = config.stab_interval;

    MSG("Restored checkpoint " + std::to_string(target) +
        ", surrogates disabled for " + std::to_string(config.stab_interval));
  } else {
    writeCheckpoint(diffusion, global_iteration, out_dir);
  }
}

bool poet::ControlModule::shouldBcastFlags() {
  if (global_iteration == 1) {
    return true;
  }
  
  if (stab_phase_ended) {
    stab_phase_ended = false;  
    return true;
  }
  
  if (flush_request) {
    return true;
  }
  
  return false;
}