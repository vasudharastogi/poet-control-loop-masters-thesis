#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

void poet::ControlModule::updateControlIteration(const uint32_t &iter,
                                                 const bool &dht_enabled,
                                                 const bool &interp_enabled) {

  /* dht_enabled and inter_enabled are user settings set before startig the
   * simulation*/
  double prep_a, prep_b;

  prep_a = MPI_Wtime();

  /*
  if (control_interval == 0) {
    control_interval_enabled = false;
    return;
  }
    */
  global_iteration = iter;
  initiateWarmupPhase(dht_enabled, interp_enabled);

  /*
  control_interval_enabled =
      (control_interval > 0 && iter % control_interval == 0);


  if (control_interval_enabled) {
    MSG("[Control] Control interval enabled at iteration " +
        std::to_string(iter));
  }
  */
  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

void poet::ControlModule::initiateWarmupPhase(bool dht_enabled,
                                              bool interp_enabled) {

  // user requested DHT/INTEP? keep them disabled but enable warmup-phase so
  if (global_iteration < stabilization_interval || rollback_enabled) {
    chem->SetWarmupEnabled(true);
    chem->SetDhtEnabled(false);
    chem->SetInterpEnabled(false);

    MSG("Stabilization enabled until next control interval at iteration " +
        std::to_string(stabilization_interval) + ".");

    if (sur_disabled_counter > 0) {
      --sur_disabled_counter;
      MSG("Rollback counter: " + std::to_string(sur_disabled_counter));
    } else {
      rollback_enabled = false;
    }

    return;
  }

  chem->SetWarmupEnabled(false);
  chem->SetDhtEnabled(dht_enabled);
  chem->SetInterpEnabled(interp_enabled);
}

void poet::ControlModule::applyControlLogic(DiffusionModule &diffusion,
                                            uint32_t &iter) {

  /*
  if (!control_interval_enabled) {
    return;
  }
  */
  writeCheckpointAndMetrics(diffusion, iter);

  if (checkAndRollback(diffusion, iter)) {
    rollback_enabled = true;
    rollback_count++;
    sur_disabled_counter = stabilization_interval;

    MSG("Interpolation disabled for the next " +
        std::to_string(stabilization_interval) + ".");
  }
}

void poet::ControlModule::writeCheckpointAndMetrics(DiffusionModule &diffusion,
                                                    uint32_t iter) {

  double w_check_a, w_check_b, stats_a, stats_b;
  MSG("Writing checkpoint of iteration " + std::to_string(iter));

  w_check_a = MPI_Wtime();
  write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                   {.field = diffusion.getField(), .iteration = iter});
  w_check_b = MPI_Wtime();
  this->w_check_t += w_check_b - w_check_a;

  stats_a = MPI_Wtime();
  writeStatsToCSV(metricsHistory, species_names, out_dir, "stats_overview");
  stats_b = MPI_Wtime();

  this->stats_t += stats_b - stats_a;
}

bool poet::ControlModule::checkAndRollback(DiffusionModule &diffusion,
                                           uint32_t &iter) {
  double r_check_a, r_check_b;

  if (global_iteration < stabilization_interval) {
    return false;
  }

  if (metricsHistory.empty()) {
    MSG("No error history yet; skipping rollback check.");
    return false;
  }

  const auto &mape = metricsHistory.back().mape;

  for (size_t row = 0; row < mape.size(); row++) {
    for (size_t col = 0; col < species_names.size() && col < mape[row].size(); col++) {
      if (mape[row][col] == 0) {
        continue;
      }

      if (mape[row][col] > mape_threshold[col]) {
        uint32_t rollback_iter =
            ((iter - 1) / checkpoint_interval) * checkpoint_interval;

        MSG("[THRESHOLD EXCEEDED] " + species_names[col] +
            " has MAPE = " + std::to_string(mape[row][col]) +
            " exceeding threshold = " + std::to_string(mape_threshold[col]) +
            ", rolling back to iteration " + std::to_string(rollback_iter));

        r_check_a = MPI_Wtime();
        Checkpoint_s checkpoint_read{.field = diffusion.getField()};
        read_checkpoint(out_dir,
                        "checkpoint" + std::to_string(rollback_iter) + ".hdf5",
                        checkpoint_read);
        iter = checkpoint_read.iteration;
        r_check_b = MPI_Wtime();
        r_check_t += r_check_b - r_check_a;
        return true;
      }
    }
  }
  MSG("All species are within their MAPE thresholds.");

  return false;
}

void poet::ControlModule::computeSpeciesErrorMetrics(
    std::vector<std::vector<double>> &reference_values,
    std::vector<std::vector<double>> &surrogate_values) {

  const uint32_t num_cells = reference_values.size();
  const uint32_t species_count = this->species_names.size();

  std::cout << "[DEBUG] computeSpeciesErrorMetrics: num_cells=" << num_cells 
            << ", species_count=" << species_count << std::endl;

  SpeciesErrorMetrics metrics(num_cells, species_count, global_iteration,
                              rollback_count);

  if (reference_values.size() != surrogate_values.size()) {
    MSG(" Reference and surrogate vectors differ in size: " +
        std::to_string(reference_values.size()) + " vs " +
        std::to_string(surrogate_values.size()));
    return;
  }

  for (size_t cell_i = 0; cell_i < num_cells; cell_i++) {

    metrics.id.push_back(reference_values[cell_i][0]);

    for (size_t sp_i = 0; sp_i < reference_values[cell_i].size(); sp_i++) {
      const double ref_value = reference_values[cell_i][sp_i];
      const double sur_value = surrogate_values[cell_i][sp_i];
      const double ZERO_ABS = 1e-13;

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
  metricsHistory.push_back(metrics);
  std::cout << "[DEBUG] metricsHistory.size()=" << metricsHistory.size() << std::endl;
}
