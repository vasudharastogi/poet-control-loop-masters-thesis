#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

void poet::ControlModule::UpdateControlIteration(const uint32_t &iter,
                                                 const bool &dht_enabled,
                                                 const bool &interp_enabled) {

  /* dht_enabled and inter_enabled are user settings set before startig the
   * simulation*/

  if (control_interval == 0) {
    control_interval_enabled = false;
    return;
  }
  // InitiateWarmupPhase(dht_enabled, interp_enabled);
  global_iteration = iter;

  if (global_iteration <= control_interval) {
    chem->SetWarmupEnabled(true);
    chem->SetDhtEnabled(false);
    chem->SetInterpEnabled(false);
    MSG("Warmup enabled until first control interval at iteration " +
        std::to_string(control_interval) + ".");
  } else {
    chem->SetWarmupEnabled(false);
    chem->SetDhtEnabled(true);
    chem->SetInterpEnabled(true);
  }

  control_interval_enabled =
      (control_interval > 0 && iter % control_interval == 0);

  if (control_interval_enabled) {
    MSG("[Control] Control interval enabled at iteration " +
        std::to_string(iter));
  }
}

void poet::ControlModule::InitiateWarmupPhase(bool dht_enabled,
                                              bool interp_enabled) {

  // user requested DHT/INTEP? keep them disabled but enable warmup-phase so
  // workers do prepareKeys/fillDHT/writePairs as required.
  if (global_iteration < control_interval) {
    /* warmup phase: keep dht and interp disabled,
   workers do prepareKeys/fillDHT/writePairs*/
    chem->SetWarmupEnabled(true);
    // chem->SetDhtEnabled(false);
    // chem->SetInterpEnabled(false);
    MSG("Warmup enabled until first control interval at iteration " +
        std::to_string(control_interval) + ".");
  } else {
    /* after warmup phase: restore according to user's request*/
    chem->SetWarmupEnabled(false);
    // chem->SetDhtEnabled(dht_enabled);
    // chem->SetInterpEnabled(interp_enabled);
  }
}

/*
void poet::ControlModule::beginIteration() {
  if (rollback_enabled) {
    if (sur_disabled_counter > 0) {
      sur_disabled_counter--;
      MSG("Rollback counter: " + std::to_string(sur_disabled_counter));
    } else {
      rollback_enabled = false;
    }
  }
}
*/

void poet::ControlModule::EndIteration(const uint32_t iter) {

  if (!control_interval_enabled) {
    return;
  }
  /* Writing a checkpointing */
  /* Control Logic*/
  if (!chem) {
    MSG("chem pointer is null — skipping checkpoint/stats write");
  } else {
    MSG("Writing checkpoint of iteration " + std::to_string(iter));
    write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                     {.field = chem->getField(), .iteration = iter});
    writeStatsToCSV(error_history, species_names, out_dir, "stats_overview");

    // if()
    /*

    if (triggerRollbackIfExceeded(*chem, *params, iter)) {
      rollback_enabled = true;
      rollback_counter++;
      sur_disabled_counter = control_interval;
      MSG("Interpolation disabled for the next " +
          std::to_string(control_interval) + ".");
    }

    */
  }
}

/*
void poet::ControlModule::BCastControlFlags() {
  int interp_flag = rollback_enabled ? 0 : 1;
  int dht_fill_flag = rollback_enabled ? 1 : 0;
  chem->ChemBCast(&interp_flag, 1, MPI_INT);
  chem->ChemBCast(&dht_fill_flag, 1, MPI_INT);
}

*/


bool poet::ControlModule::RollbackIfThresholdExceeded(ChemistryModule &chem) {

  /**
  if (error_history.empty()) {
    MSG("No error history yet; skipping rollback check.");
    return false;
  }

  const auto &mape = error_history.back().mape;

  for (uint32_t i = 0; i < species_names.size(); ++i) {
    if (mape[i] == 0) {
      continue;
    }

    if (mape[i] > mape_threshold[i]) {
      uint32_t rollback_iter = ((global_iteration - 1) / checkpoint_interval) *
                               checkpoint_interval;

      MSG("[THRESHOLD EXCEEDED] " + species_names[i] +
          " has MAPE = " + std::to_string(mape[i]) +
          " exceeding threshold = " + std::to_string(mape_threshold[i])
+ " → rolling back to iteration " + std::to_string(rollback_iter));

      Checkpoint_s checkpoint_read{.field = chem.getField()};
      read_checkpoint(out_dir,
                      "checkpoint" + std::to_string(rollback_iter) + ".hdf5",
                      checkpoint_read);
      global_iteration = checkpoint_read.iteration;
      return true;
    }
  }
  MSG("All species are within their MAPE thresholds.");
  return false;
  */
}

void poet::ControlModule::ComputeSpeciesErrorMetrics(
    const std::vector<double> &reference_values,
    const std::vector<double> &surrogate_values, const uint32_t size_per_prop) {

  SimulationErrorStats species_error_stats(this->species_names.size(),
                                           global_iteration,
                                           /*rollback_counter*/ 0);

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
    uint32_t nan_count = 0;
    uint32_t valid_count = 0;

    for (uint32_t j = 0; j < size_per_prop; ++j) {
      const double ref_value = reference_values[base_idx + j];
      const double sur_value = surrogate_values[base_idx + j];
      const double ZERO_ABS = 1e-13;

      if (std::isnan(ref_value) || std::isnan(sur_value)) {
        nan_count++;
        continue;
      }
      valid_count++;

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
    if (valid_count > 0) {
      species_error_stats.mape[i] = 100.0 * (err_sum / valid_count);
      species_error_stats.rrmse[i] = std::sqrt(sqr_err_sum / valid_count);
    } else {
      species_error_stats.mape[i] = 0.0;
      species_error_stats.rrmse[i] = 0.0;
      std::cerr << "[CTRL WARN] no valid samples for species " << i << " ("
                << this->species_names[i] << "), setting errors to 0\n";
    }
  }
  error_history.push_back(species_error_stats);
}
