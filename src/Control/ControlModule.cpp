#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

poet::ControlModule::ControlModule(const ControlConfig &cfg, ChemistryModule *chm)
    : config(cfg), chem(chm) {
  assert(chem && "ChemistryModule pointer must not be null");
}

bool poet::ControlModule::inWarmup() const { return global_iter <= config.ctrl_interval; }

void poet::ControlModule::beginIteration(const uint32_t &iter, const bool &dht_enabled,
                                         const bool &interp_enabled) {
  double prep_a, prep_b;

  prep_a = MPI_Wtime();
  if (config.ctrl_interval == 0) {
    ctrl_active = false;
    return;
  }
  global_iter = iter;
  ctrl_active = (config.ctrl_interval > 0 && (iter % config.ctrl_interval == 0));
  updateSurrState(dht_enabled, interp_enabled);
  prep_b = MPI_Wtime();
  this->prep_t += prep_b - prep_a;
}

void poet::ControlModule::setSurrState(const SurrState &state) {
  chem->SetStabEnabled(state.stab_enabled);
  chem->SetDhtEnabled(state.dht_enabled);
  chem->SetInterpEnabled(state.interp_enabled);
}

void poet::ControlModule::trackStabPhase() {
  if (!rb_enabled || stab_countdown == 0 || rbLimitReached()) {
    return;
  }
  --stab_countdown;
  std::cout << "Stabilization countdown: " << stab_countdown << std::endl;
  if (stab_countdown == 0) {
    rb_enabled = false;
    std::cout << "Stabilization complete, surrogate re-enabled." << std::endl;
  }
  flush_request = false;
}

/* Decrease rollback counter after stable surrogate uptime */
void poet::ControlModule::trackSurrUptime() {
  if (rb_enabled || rb_count == 0 || inWarmup())
    return;
  if (config.rb_aging_limit == 0 || config.rb_limit == 0)
    return;
  ++surr_active;
  if (surr_active >= config.rb_aging_limit) {
    surr_active = 0;
    --rb_count;
    std::cout << "Rollback count decreased by one: " << rb_count << std::endl;
  }
}

/* manages the overall surrogate state, by enabling/disabling state based on
 * warmup logic and rollback conditions*/
void poet::ControlModule::updateSurrState(bool dht_enabled, bool interp_enabled) {

  if (inWarmup()) {
    setSurrState({true, false, false});
    return;
  }
  if (rbLimitReached()) {
    std::cout << "Interpolation completely disabled." << std::endl;
    setSurrState({false, false, false});
    return;
  }
  if (rb_enabled) {
    std::cout << "In stabilization phase." << std::endl;
    setSurrState({true, false, false});
    trackStabPhase();
    return;
  }

  std::cout << "Interpolating." << std::endl;
  setSurrState({false, dht_enabled, interp_enabled});
  if (interp_enabled) trackSurrUptime();
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

void poet::ControlModule::readCheckpoint(uint32_t &curr_iter, uint32_t rollback_iter,
                                         const std::string &out_dir) {
  double r_check_a, r_check_b;
  r_check_a = MPI_Wtime();
  Checkpoint_s checkpoint_read{.field = chem->getField()};
  read_checkpoint(out_dir, "checkpoint" + std::to_string(rollback_iter) + ".hdf5",
                  checkpoint_read);
  curr_iter = checkpoint_read.iteration;
  r_check_b = MPI_Wtime();
  r_check_t += r_check_b - r_check_a;
}

void poet::ControlModule::writeMetrics(const std::string &out_dir,
                                       const std::vector<std::string> &species) {

  /*
    if (rb_count > config.rb_limit) {
    return;
  }
  */
  double stats_a, stats_b;

  stats_a = MPI_Wtime();
  writeStatsToCSV(s_history, species, out_dir, "metrics_overview", config);
  stats_b = MPI_Wtime();

  this->stats_t += stats_b - stats_a;
}

uint32_t poet::ControlModule::calcRbIter() {

  uint32_t last_iter =
      ((global_iter - 1) / config.chkpt_interval) * config.chkpt_interval;

  uint32_t rb_iter = (last_iter <= last_chkpt_written) ? last_iter : last_chkpt_written;
  return rb_iter;
}

std::optional<uint32_t>
poet::ControlModule::findRbTarget(const std::vector<std::string> &species) {

  if (s_history.empty()) {
    std::cout << "No error history yet, skipping rollback check." << std::endl;
    return std::nullopt;
  }
  if (rbLimitReached()) {
    std::cout << "Rollback limit reached, skipping control logic." << std::endl;
    return std::nullopt;
  }

  std::cout << "findRbTarget called at iter " << global_iter << ", rb_count=" << rb_count
            << ", rb_limit=" << config.rb_limit << std::endl;

  double r_check_a, r_check_b;
  const auto &mape = s_history.back().mape;
  for (uint32_t sp_idx = 0; sp_idx < species.size(); ++sp_idx) {

    if (mape[sp_idx] == 0) {
      continue;
    }
    /* skip Charge */
    if (sp_idx == 4) {
      continue;
    }

    if (mape[sp_idx] > config.mape_threshold[sp_idx]) {
      std::cout << "Species " << species[sp_idx] << " MAPE=" << mape[sp_idx]
                << " threshold=" << config.mape_threshold[sp_idx] << std::endl;

      if (last_chkpt_written == 0) {
        std::cout << " Threshold exceeded but no checkpoint exists yet." << std::endl;
        return std::nullopt;
      }
      flush_request = true;
      std::cout << "Threshold exceeded " << species[sp_idx]
                << " has MAPE = " << std::to_string(mape[sp_idx])
                << " exceeding threshold = "
                << std::to_string(config.mape_threshold[sp_idx]) << std::endl;
      return calcRbIter();
    }
  }
  return std::nullopt;
}

void poet::ControlModule::computeMetrics(const std::vector<double> &reference_values,
                                         const std::vector<double> &surrogate_values,
                                         const uint32_t size_per_prop,
                                         const std::vector<std::string> &species) {
  if (rbLimitReached()) {
    return;
  }
  SpeciesMetrics s_metrics(species.size(), global_iter, rb_count);

  for (uint32_t sp_idx = 0; sp_idx < species.size(); ++sp_idx) {
    double err_sum = 0.0;
    double sqr_err_sum = 0.0;
    uint32_t base_idx = sp_idx * size_per_prop;

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
    s_metrics.mape[sp_idx] = 100.0 * (err_sum / size_per_prop);
    s_metrics.rrmse[sp_idx] = std::sqrt(sqr_err_sum / size_per_prop);
  }
  s_history.push_back(s_metrics);
}

void poet::ControlModule::processCheckpoint(uint32_t &curr_iter,
                                            const std::string &out_dir,
                                            const std::vector<std::string> &species) {
  if (!ctrl_active || rbLimitReached()) {
    return;
  }
  if (flush_request) {
    triggerRb(curr_iter, out_dir);
  } else {
    writeCheckpoint(global_iter, out_dir);
  }
}

void poet::ControlModule::triggerRb(uint32_t &curr_iter, const std::string &out_dir) {

  uint32_t target = calcRbIter();
  readCheckpoint(curr_iter, target, out_dir);

  surr_active = 0;
  rb_enabled = true;
  rb_count++;
  stab_countdown = config.ctrl_interval;

  std::cout << "Restored checkpoint " << std::to_string(target)
            << ", surrogates disabled for " << config.ctrl_interval << std::endl;
}

bool poet::ControlModule::needsFlagBcast() const {
  if (global_iter == 1 || global_iter % config.ctrl_interval == 1) {
    return true;
  }
  return false;
}

inline bool poet::ControlModule::rbLimitReached() const {
  /* rollback is completly disabled */
  if (config.rb_limit == 0)
    return false;
  return rb_count >= config.rb_limit;
}
