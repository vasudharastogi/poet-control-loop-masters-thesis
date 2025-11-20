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
      MSG("Rollback stabilization complete, re-enabling surrogates");
    }
  }

  bool prev_stab_state = chem.GetStabEnabled();

  // user requested DHT/INTEP? keep them disabled but enable warmup-phase so
  if (global_iteration <= config.stab_interval || rollback_enabled) {
    chem.SetStabEnabled(true);
    chem.SetDhtEnabled(false);
    chem.SetInterpEnabled(false);
  } else {
    chem.SetStabEnabled(false);
    chem.SetDhtEnabled(dht_enabled);
    chem.SetInterpEnabled(interp_enabled);
  }

  // Mark that we need to broadcast flags if stab state changed
  if (prev_stab_state != chem.GetStabEnabled()) {
    stab_phase_ended = true;
  }
}

void poet::ControlModule::writeCheckpoint(DiffusionModule &diffusion,
                                          uint32_t &iter,
                                          const std::string &out_dir) {
  double w_check_a, w_check_b;

  w_check_a = MPI_Wtime();
  if (global_iteration % config.checkpoint_interval == 0) {
    write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                     {.field = diffusion.getField(), .iteration = iter});
  }
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
  writeSpeciesStatsToCSV(metrics_history, species, out_dir,
                         "species_overview.csv");
  writeCellStatsToCSV(cell_metrics_history, species, out_dir,
                      "cell_overview.csv");
  write_metrics(cell_metrics_history, species, out_dir,
                "metrics_overview.hdf5");
  stats_b = MPI_Wtime();

  this->stats_t += stats_b - stats_a;
}

uint32_t poet::ControlModule::getRollbackIter() {

  uint32_t last_iter = ((global_iteration - 1) / config.checkpoint_interval) *
                       config.checkpoint_interval;

  uint32_t rollback_iter = (last_iter <= last_checkpoint_written)
                               ? last_iter
                               : last_checkpoint_written;

  MSG("getRollbackIter: global_iteration=" + std::to_string(global_iteration) +
      ", checkpoint_interval=" + std::to_string(config.checkpoint_interval) +
      ", last_checkpoint_written=" + std::to_string(last_checkpoint_written) +
      ", returning=" + std::to_string(last_checkpoint_written));

  return last_checkpoint_written;
}

std::optional<uint32_t> poet::ControlModule::getRollbackTarget(
    const std::vector<std::string> &species) {
  double r_check_a, r_check_b;

  if (metrics_history.empty()) {
    MSG("No error history yet, skipping rollback check.");
    rollback_enabled = false;
    return std::nullopt;
  }
  // Skip threshold checking if already in rollback/stabilization phase
  if (rollback_enabled) {
    return std::nullopt;
  }

  const auto &s_mape = metrics_history.back().mape;

  for (size_t sp_i = 2; sp_i < species.size(); sp_i++) {

    if (s_mape[sp_i] == 0) {
      continue;
    }
    if (s_mape[sp_i] > config.mape_threshold[sp_i]) {
      if (last_checkpoint_written == 0) {
        MSG(" Threshold exceeded but no checkpoint exists yet.");
        return std::nullopt;
      }

      const auto &c_mape = cell_metrics_history.back().mape;
      const auto &c_id = cell_metrics_history.back().id;

      auto max_it = std::max_element(
          c_mape.begin(), c_mape.end(),
          [sp_i](const auto &a, const auto &b) { return a[sp_i] < b[sp_i]; });

      size_t max_idx = std::distance(c_mape.begin(), max_it);
      uint32_t cell_id = c_id[max_idx];
      double cell_mape = (*max_it)[sp_i];

      rollback_enabled = true;
      flush_request = true;

      MSG("Threshold exceeded for " + species[sp_i] +
          " with species-level MAPE = " + std::to_string(s_mape[sp_i]) +
          " exceeding threshold = " +
          std::to_string(config.mape_threshold[sp_i]) + ". Worst cell: ID=" +
          std::to_string(cell_id) + " with MAPE=" + std::to_string(cell_mape));

      return getRollbackIter();
    }
  }
  rollback_enabled = false;
  flush_request = false;
  return std::nullopt;
}

void poet::ControlModule::computeErrorMetrics(
    std::vector<std::vector<double>> &reference_values,
    std::vector<std::vector<double>> &surrogate_values,
    const std::vector<std::string> &species, const uint32_t size_per_prop) {

  // Skip metric computation if already in rollback/stabilization phase
  if (rollback_enabled) {
    return;
  }

  const uint32_t n_cells = reference_values.size();
  const uint32_t n_species = species.size();
  const double ZERO_ABS = config.zero_abs;

  CellErrorMetrics c_metrics(n_cells, n_species, global_iteration,
                             rollback_count);
  SpeciesErrorMetrics s_metrics(n_species, global_iteration, rollback_count);

  std::vector<double> species_err_sum(n_species, 0.0);
  std::vector<double> species_sqr_sum(n_species, 0.0);

  for (size_t cell_i = 0; cell_i < n_cells; cell_i++) {

    c_metrics.id.push_back(reference_values[cell_i][0]);

    for (size_t sp_i = 2; sp_i < n_species; sp_i++) {
      const double ref_value = reference_values[cell_i][sp_i];
      const double sur_value = surrogate_values[cell_i][sp_i];

      if (std::isnan(ref_value) || std::isnan(sur_value)) {
        continue;
      }
      if (std::abs(ref_value) < ZERO_ABS) {
        if (std::abs(sur_value) >= ZERO_ABS) {

          species_err_sum[sp_i] += 1.0;
          species_sqr_sum[sp_i] += 1.0;

          c_metrics.mape[cell_i][sp_i] = 100.0;
          c_metrics.rrmse[cell_i][sp_i] = 1.0;
        }
      } else {
        double alpha = 1.0 - (sur_value / ref_value);

        species_err_sum[sp_i] += std::abs(alpha);
        species_sqr_sum[sp_i] += alpha * alpha;

        c_metrics.mape[cell_i][sp_i] = 100.0 * std::abs(alpha);
        c_metrics.rrmse[cell_i][sp_i] = alpha * alpha;
      }
    }
  }
  for (size_t sp_i = 2; sp_i < n_species; sp_i++) {
    s_metrics.mape[sp_i] = 100.0 * (species_err_sum[sp_i] / size_per_prop);
    s_metrics.rrmse[sp_i] = std::sqrt(species_sqr_sum[sp_i] / size_per_prop);
  }

  metrics_history.push_back(s_metrics);
  cell_metrics_history.push_back(c_metrics);
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
    return true;
  }

  return false;
}