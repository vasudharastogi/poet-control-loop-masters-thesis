#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

poet::ControlModule::ControlModule(const ControlConfig &config_, ChemistryModule *chem_)
    : config(config_), chem(chem_) {
  assert(chem && "ChemistryModule pointer must not be null");
}

void poet::ControlModule::beginIteration(const uint32_t &iter, const bool &dht_enabled,
                                         const bool &interp_enabled) {
  double prep_a, prep_b;

  prep_a = MPI_Wtime();
  if (config.ctrl_interval == 0) {
    ctrl_active = false;
    return;
  }
  global_iter = iter;

  updateSurrState(dht_enabled, interp_enabled);

  ctrl_active = (config.ctrl_interval > 0 && (iter % config.ctrl_interval == 0));

  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

/* manages the overall surrogate state, by enabling/disabling state based on
 * warmup logic and rollback conditions*/
void poet::ControlModule::updateSurrState(bool dht_enabled, bool interp_enabled) {

  bool in_warmup = (global_iter <= config.ctrl_interval);
  bool rb_limit_reached = (rb_count >= config.rb_limit);

  if (rb_enabled && stab_countdown > 0 && !rb_limit_reached) {
    --stab_countdown;
    std::cout << "Rollback counter: " << stab_countdown << std::endl;
    if (stab_countdown == 0) {
      rb_enabled = false;
    }
    flush_request = false;
  }
  /* disable surrogates during warmup, active rollback or after limit */
  if (in_warmup || rb_enabled || rb_limit_reached) {
    chem->SetStabEnabled(!rb_limit_reached);
    chem->SetDhtEnabled(false);
    chem->SetInterpEnabled(false);

    if (rb_limit_reached) {
      std::cout << "Interpolation completly disabled." << std::endl;
    } else {
      std::cout << "In stabilization phase." << std::endl;
    }

    return;
  }
  /* enable user-requested surrogates */
  chem->SetStabEnabled(false);
  chem->SetDhtEnabled(dht_enabled);
  chem->SetInterpEnabled(interp_enabled);
  std::cout << "Interpolating." << std::endl;
}

void poet::ControlModule::writeCheckpoint(uint32_t &iter, const std::string &out_dir) {
  double w_check_a, w_check_b;
  w_check_a = MPI_Wtime();
  write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                   {.field = chem->getField(), .iteration = iter});
  w_check_b = MPI_Wtime();
  this->w_check_t += w_check_b - w_check_a;

  last_chkpt_written = iter;
}

void poet::ControlModule::readCheckpoint(uint32_t &current_iter, uint32_t rollback_iter,
                                         const std::string &out_dir) {
  double r_check_a, r_check_b;
  r_check_a = MPI_Wtime();
  Checkpoint_s checkpoint_read{.field = chem->getField()};
  read_checkpoint(out_dir, "checkpoint" + std::to_string(rollback_iter) + ".hdf5", checkpoint_read);
  current_iter = checkpoint_read.iteration;
  r_check_b = MPI_Wtime();
  r_check_t += r_check_b - r_check_a;
}

void poet::ControlModule::writeMetrics(const std::string &out_dir,
                                       const std::vector<std::string> &species) {
  if (rb_count > config.rb_limit) {
    return;
  }
  double stats_a, stats_b;

  stats_a = MPI_Wtime();
  writeStatsToCSV(metrics_history, species, out_dir, "metrics_overview");
  stats_b = MPI_Wtime();

  this->stats_t += stats_b - stats_a;
}

uint32_t poet::ControlModule::calcRbIter() {

  uint32_t last_iter = ((global_iter - 1) / config.chkpt_interval) * config.chkpt_interval;

  uint32_t rb_iter = (last_iter <= last_chkpt_written) ? last_iter : last_chkpt_written;
  return rb_iter;
}

std::optional<uint32_t> poet::ControlModule::findRbTarget(const std::vector<std::string> &species) {

  if (metrics_history.empty()) {
    std::cout << "No error history yet, skipping rollback check." << std::endl;
    flush_request = false;
    return std::nullopt;
  }
  if (rb_count > config.rb_limit) {
    std::cout << "Rollback limit reached, skipping control logic." << std::endl;
    flush_request = false;
    return std::nullopt;
  }

  std::cout << "findRbTarget called at iter " << global_iter << ", rb_count=" << rb_count
            << ", rb_limit=" << config.rb_limit << std::endl;

  double r_check_a, r_check_b;
  const auto &mape = metrics_history.back().mape;
  for (uint32_t i = 0; i < species.size(); ++i) {

    if (mape[i] == 0) {
      continue;
    }
    if (mape[i] > config.mape_threshold[i]) {
      std::cout << "Species " << species[i] << " MAPE=" << mape[i]
                << " threshold=" << config.mape_threshold[i] << std::endl;

      if (last_chkpt_written == 0) {
        std::cout << " Threshold exceeded but no checkpoint exists yet." << std::endl;
        return std::nullopt;
      }
      // rb_enabled = true;
      flush_request = true;
      std::cout << "Threshold exceeded " << species[i] << " has MAPE = " << std::to_string(mape[i])
                << " exceeding threshold = " << std::to_string(config.mape_threshold[i])
                << std::endl;
      return calcRbIter();
    }
  }
  // std::cout << "All species are within their MAPE thresholds." << std::endl;
  flush_request = false;
  return std::nullopt;
}

void poet::ControlModule::computeMetrics(const std::vector<double> &reference_values,
                                         const std::vector<double> &surrogate_values,
                                         const uint32_t size_per_prop,
                                         const std::vector<std::string> &species) {

  if (rb_count > config.rb_limit) {
    return;
  }

  SpeciesMetrics metrics(species.size(), global_iter, rb_count);

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
      } else {
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

void poet::ControlModule::processCheckpoint(uint32_t &current_iter, const std::string &out_dir,
                                            const std::vector<std::string> &species) {

  if (!ctrl_active || rb_count > config.rb_limit) {
    return;
  }

  if (flush_request) {
    uint32_t target = calcRbIter();
    readCheckpoint(current_iter, target, out_dir);

    rb_enabled = true;
    rb_count++;
    stab_countdown = config.ctrl_interval;

    std::cout << "Restored checkpoint " << std::to_string(target) << ", surrogates disabled for "
              << config.ctrl_interval << std::endl;
  } else {
    writeCheckpoint(global_iter, out_dir);
  }
}

bool poet::ControlModule::needsFlagBcast() const {
  if (rb_count > config.rb_limit) {
    return false;
  }
  if (global_iter == 1 || global_iter % config.ctrl_interval == 1) {
    return true;
  }
  return false;
}
