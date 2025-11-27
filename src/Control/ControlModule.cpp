#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

poet::ControlModule::ControlModule(const ControlConfig &config_)
    : config(config_) {}

void poet::ControlModule::beginIteration(ChemistryModule &chem, uint32_t &iter,
                                         const bool &dht_enabled,
                                         const bool &interp_enabled) {

  global_iteration = iter;
  double prep_a, prep_b;

  prep_a = MPI_Wtime();

  updateStabilizationPhase(chem, dht_enabled, interp_enabled);
  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

/* Disables dht and/or interp during stabilzation phase */
void poet::ControlModule::updateStabilizationPhase(ChemistryModule &chem,
                                                   bool dht_enabled,
                                                   bool interp_enabled) {
  if (rollback_enabled) {
    if (disable_surr_counter > 0) {
      --disable_surr_counter;
      flush_request = false;
      std::cout << "Rollback counter: " << std::to_string(disable_surr_counter)
                << std::endl;
    } else {
      rollback_enabled = false;
      std::cout << "Rollback stabilization complete, re-enabling surrogate."
                << std::endl;
    }
  }

  if (global_iteration <= config.stab_interval || rollback_enabled) {
    chem.SetStabEnabled(true);
    chem.SetDhtEnabled(false);
    chem.SetInterpEnabled(false);
  } else {
    chem.SetStabEnabled(false);
    chem.SetDhtEnabled(dht_enabled);
    chem.SetInterpEnabled(interp_enabled);
  }
}

void poet::ControlModule::writeCheckpoint(DiffusionModule &diffusion,
                                          uint32_t &iter,
                                          const std::string &out_dir) {
  if (global_iteration % config.checkpoint_interval == 0) {
    double w_check_a = MPI_Wtime();
    write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                     {.field = diffusion.getField(), .iteration = iter});
    double w_check_b = MPI_Wtime();
    this->w_check_t += w_check_b - w_check_a;
    last_checkpoint_written = iter;
  }
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
    uint32_t &iter, const std::string &out_dir,
    const std::vector<std::string> &species) {

  double stats_a = MPI_Wtime();
  writeSpeciesStatsToCSV(metrics_history, species, out_dir,
                         "species_overview.csv");
  write_metrics(cell_metrics_history, species, out_dir,
                "metrics_overview.hdf5");
  double stats_b = MPI_Wtime();
  this->stats_t += stats_b - stats_a;
}

uint32_t poet::ControlModule::getRollbackIter() {

  uint32_t last_iter = ((global_iteration - 1) / config.checkpoint_interval) *
                       config.checkpoint_interval;
  return last_iter;
}

std::optional<uint32_t> poet::ControlModule::getRollbackTarget(
    const std::vector<std::string> &species) {
  double r_check_a, r_check_b;

  /* Skip threshold checking if already in stabilization phase*/
  if (metrics_history.empty() || rollback_enabled) {
    return std::nullopt;
  }

  const auto &s_hist = metrics_history.back();

  /* skipping cell_id and id */
  for (size_t sp_i = 2; sp_i < species.size(); sp_i++) {

    /* check bounds of threshold vector*/
    if (sp_i >= config.mape_threshold.size()) {
      std::cerr << "Warning: No threshold defined for species " << species[sp_i]
                << " at index " << std::to_string(sp_i) << std::endl;
      continue;
    }
    if (s_hist.mape[sp_i] > config.mape_threshold[sp_i]) {

      const auto &c_hist = cell_metrics_history.back();
      auto max_it = std::max_element(
          c_hist.mape.begin(), c_hist.mape.end(),
          [sp_i](const auto &a, const auto &b) { return a[sp_i] < b[sp_i]; });

      size_t max_idx = std::distance(c_hist.mape.begin(), max_it);
      uint32_t cell_id = c_hist.id[max_idx];
      double cell_mape = (*max_it)[sp_i];

      rollback_enabled = true;
      flush_request = true;

      std::cout << "Threshold exceeded for " << species[sp_i]
                << " with species-level MAPE = " << std::to_string(s_hist.mape[sp_i])
                << " exceeding threshold = "
                << std::to_string(config.mape_threshold[sp_i])
                << ". Worst cell: ID=" << std::to_string(cell_id)
                << ", species= " << species[sp_i] 
                << " with MAPE=" << std::to_string(cell_mape) << std::endl;

      return getRollbackIter();
    }
  }
  return std::nullopt;
}

void poet::ControlModule::computeErrorMetrics(
    std::vector<std::vector<double>> &reference_values,
    std::vector<std::vector<double>> &surrogate_values,
    const std::vector<std::string> &species, const uint32_t size_per_prop) {

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
    c_metrics.mape[cell_i][0] = reference_values[cell_i][0];
    c_metrics.rrmse[cell_i][0] = reference_values[cell_i][0];

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
        // Log extreme MAPE values for debugging
        if (c_metrics.mape[cell_i][sp_i] > 100.0) {
          std::cout << "WARNING: High MAPE detected - Cell="
                    << c_metrics.id[cell_i] << ", Species=" << species[sp_i]
                    << ", MAPE=" << c_metrics.mape[cell_i][sp_i]
                    << "%, Ref=" << ref_value << ", Sur=" << sur_value
                    << ", Alpha=" << alpha << std::endl;
        }
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

  // Use max_rollbacks from config, default to 3 if not set
  // uint32_t max_rollbacks =
  //  (config.max_rollbacks > 0) ? config.max_rollbacks : 3;

  if (flush_request /* && rollback_count < 3 */) {
    uint32_t target = getRollbackIter();
    readCheckpoint(diffusion, current_iter, target, out_dir);

    rollback_enabled = true;
    rollback_count++;
    disable_surr_counter = config.stab_interval;

    std::cout << "Restored checkpoint " << std::to_string(target)
              << ", surrogate disabled for "
              << std::to_string(config.stab_interval) << std::endl;
  } else {
    writeCheckpoint(diffusion, global_iteration, out_dir);
  }
}