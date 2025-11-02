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
  if (control_interval == 0) {
    control_interval_enabled = false;
    return;
  }
  global_iteration = iter;
  initiateWarmupPhase(dht_enabled, interp_enabled);

  control_interval_enabled =
      (control_interval > 0 && iter % control_interval == 0);

  if (control_interval_enabled) {
    MSG("[Control] Control interval enabled at iteration " +
        std::to_string(iter));
  }
  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

void poet::ControlModule::initiateWarmupPhase(bool dht_enabled,
                                              bool interp_enabled) {

  // user requested DHT/INTEP? keep them disabled but enable warmup-phase so
  if (global_iteration <= control_interval || rollback_enabled) {
    chem->SetWarmupEnabled(true);
    chem->SetDhtEnabled(false);
    chem->SetInterpEnabled(false);
   // MSG("Warmup enabled until next control interval at iteration " +
     //   std::to_string(control_interval) + ".");

    if (rollback_enabled) {
      if (sur_disabled_counter > 0) {
        --sur_disabled_counter;
        //MSG("Rollback counter: " + std::to_string(sur_disabled_counter));
      } else {
        rollback_enabled = false;
      }
    }
    return;
  }

  chem->SetWarmupEnabled(false);
  chem->SetDhtEnabled(dht_enabled);
  chem->SetInterpEnabled(interp_enabled);
}

void poet::ControlModule::applyControlLogic(ChemistryModule &chem,
                                            uint32_t &iter) {
  if (!control_interval_enabled) {
    return;
  }
  writeCheckpointAndMetrics(chem, iter);

  if (checkAndRollback(chem, iter) /* && rollback_count < 4*/) {
    rollback_enabled = true;
    rollback_count++;
    sur_disabled_counter = control_interval;
    /*
    MSG("Interpolation disabled for the next " +
        std::to_string(control_interval) + ".");
    */
  }
}

void poet::ControlModule::writeCheckpointAndMetrics(ChemistryModule &chem,
                                                    uint32_t iter) {

  double w_check_a, w_check_b, stats_a, stats_b;
  MSG("Writing checkpoint of iteration " + std::to_string(iter));

  w_check_a = MPI_Wtime();
  write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                   {.field = chem.getField(), .iteration = iter});
  w_check_b = MPI_Wtime();
  this->w_check_t += w_check_b - w_check_a;

  stats_a = MPI_Wtime();
  writeStatsToCSV(metricsHistory, species_names, out_dir, "stats_overview");
  stats_b = MPI_Wtime();

  this->stats_t += stats_b - stats_a;
}

bool poet::ControlModule::checkAndRollback(ChemistryModule &chem,
                                           uint32_t &iter) {
  double r_check_a, r_check_b;

  if (metricsHistory.empty()) {
    MSG("No error history yet; skipping rollback check.");
    return false;
  }

  const auto &mape = metricsHistory.back().mape;

  for (uint32_t i = 0; i < species_names.size(); ++i) {
    if (mape[i] == 0) {
      continue;
    }

    if (mape[i] > mape_threshold[i]) {
      uint32_t rollback_iter =
          ((iter - 1) / checkpoint_interval) * checkpoint_interval;

      MSG("[THRESHOLD EXCEEDED] " + species_names[i] +
          " has MAPE = " + std::to_string(mape[i]) +
          " exceeding threshold = " + std::to_string(mape_threshold[i]) +
          " → rolling back to iteration " + std::to_string(rollback_iter));

      r_check_a = MPI_Wtime();
      Checkpoint_s checkpoint_read{.field = chem.getField()};
      read_checkpoint(out_dir,
                      "checkpoint" + std::to_string(rollback_iter) + ".hdf5",
                      checkpoint_read);
      iter = checkpoint_read.iteration;
      r_check_b = MPI_Wtime();
      r_check_t += r_check_b - r_check_a;
      return true;
    }
  }
  MSG("All species are within their MAPE thresholds.");

  return false;
}

void poet::ControlModule::computeSpeciesErrorMetrics(
    const std::vector<double> &reference_values,
    const std::vector<double> &surrogate_values, const uint32_t size_per_prop) {

  SpeciesErrorMetrics metrics(this->species_names.size(), global_iteration,
                              rollback_count);

  if (reference_values.size() != surrogate_values.size()) {
    MSG(" Reference and surrogate vectors differ in size: " +
        std::to_string(reference_values.size()) + " vs " +
        std::to_string(surrogate_values.size()));
    return;
  }

  const std::size_t expected =
      static_cast<std::size_t>(this->species_names.size()) * size_per_prop;
  if (reference_values.size() < expected) {
    std::cerr << "[CTRL ERROR] input vectors too small: expected >= "
              << expected << " entries, got " << reference_values.size()
              << "\n";
    return;
  }

  for (uint32_t i = 0; i < this->species_names.size(); ++i) {
    double err_sum = 0.0;
    double sqr_err_sum = 0.0;
    uint32_t base_idx = i * size_per_prop;

    int count = 0;

    for (uint32_t j = 0; j < size_per_prop; ++j) {
      const double ref_value = reference_values[base_idx + j];
      const double sur_value = surrogate_values[base_idx + j];
      const double ZERO_ABS = 1e-13;

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
  metricsHistory.push_back(metrics);
}
