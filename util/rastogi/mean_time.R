#!/usr/bin/env Rscript


#suppressPackageStartupMessages({
#  if (!requireNamespace("qs2", quietly = TRUE)) {
#    stop("Package 'qs2' not installed. Install with: remotes::install_github('qs2io/qs2')")
#  }
#})

# calculate mean timing information from multiple qs2 files

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript mean_time.R timings.qs2 [more.qs2 ...]")
}

chem_vectors <- c(
  "idle_worker", "phreeqc_time", "dht_get_time", "dht_fill_time",
  "interp_w", "interp_r", "interp_g", "interp_fc"
)

format_hms <- function(sec) {
  if (is.na(sec)) return(NA_character_)
  h <- floor(sec / 3600)
  m <- floor((sec - h * 3600) / 60)
  s <- sec - h * 3600 - m * 60
  sprintf("%02d:%02d:%06.2f", h, m, s)
}

results <- lapply(args, function(f) {
  if (!file.exists(f)) {
    warning("File not found: ", f)
    return(NULL)
  }
  obj <- qs2::qs_read(f)

  #obj<-readRDS(f)

  chem <- obj$chemistry
  ctrl <- obj$control_loop

  # Vector means (seconds)
  vec_means <- sapply(chem_vectors, function(v) {
    x <- chem[[v]]
    if (!is.numeric(x)) return(NA_real_)
    mean(x, na.rm = TRUE)
  })

  # Worker control loop mean
  worker_mean <- if (is.numeric(ctrl$worker)) mean(ctrl$worker, na.rm = TRUE) else NA_real_

  # Vector means stored only as formatted h:m:s (seconds removed per user request)
  all_means <- c(vec_means, worker_mean)
  vector_df <- data.frame(
    vector = c(chem_vectors, "worker"),
    hms_mean = vapply(all_means, format_hms, character(1)),
    stringsAsFactors = FALSE
  )

  # Scalars: these are already single values, just format them
  scalar_names_chem <- c("simtime", "loop", "sequential", "idle_master")
  scalar_vals_chem <- sapply(scalar_names_chem, function(n) if (is.numeric(chem[[n]])) chem[[n]] else NA_real_)

  scalar_names_ctrl <- c(
    "compute_metrics_master", "unshuffle_field_master", "w_checkpoint_master",
    "r_checkpoint_master", "write_stats", "ctrl_logic_master", "recv_data_master"
  )
  scalar_vals_ctrl <- sapply(scalar_names_ctrl, function(n) if (is.numeric(ctrl[[n]])) ctrl[[n]] else NA_real_)

  diffusion_simtime <- if (!is.null(obj$diffusion$simtime) && is.numeric(obj$diffusion$simtime)) obj$diffusion$simtime else NA_real_
  simtime_total <- if (!is.null(obj$simtime) && is.numeric(obj$simtime)) obj$simtime else NA_real_

  # Scalars: keep only formatted h:m:s (remove raw seconds column)
  scalar_seconds <- c(simtime_total, diffusion_simtime, scalar_vals_chem, scalar_vals_ctrl)
  scalar_df <- data.frame(
    scalar = c(
      "simtime_total", "diffusion_simtime",
      paste0("chem_", scalar_names_chem),
      paste0("ctrl_", scalar_names_ctrl)
    ),
    hms = vapply(scalar_seconds, format_hms, character(1)),
    stringsAsFactors = FALSE
  )

  list(file = basename(f), vectors = vector_df, scalars = scalar_df)
})

# Combine
vector_summary <- do.call(rbind, lapply(results, `[[`, "vectors"))
scalar_summary <- do.call(rbind, lapply(results, `[[`, "scalars"))

cat("\nMean timing (per vector) in h:m:s:\n\n")
print(vector_summary, row.names = FALSE)

cat("\nScalar timings (single values) in h:m:s:\n\n")
print(scalar_summary, row.names = FALSE)

write.csv(vector_summary, "timings_vector_mean_summary.csv", row.names = FALSE)
write.csv(scalar_summary, "timings_scalar_values.csv", row.names = FALSE)
cat("\nSaved vector mean summary to timings_vector_mean_summary.csv")
cat("\nSaved scalar values to timings_scalar_values.csv\n")
