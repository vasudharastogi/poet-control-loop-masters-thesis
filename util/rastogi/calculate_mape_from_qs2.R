#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(qs2)})
# read qs2 files and calculate MAPE, RRMSE, and MAXMAPE between reference and predicted data

# alpha_i definition exactly as in the paper (vectorized)
calc_alpha <- function(actual, predicted) {
  stopifnot(length(actual) == length(predicted))
  
  alpha <- numeric(length(actual))
  
  idx_nonzero <- actual != 0
  idx_zero    <- actual == 0
  
  # if y_i != 0: (y_i - yhat_i)/y_i
  alpha[idx_nonzero] <- (actual[idx_nonzero] - predicted[idx_nonzero]) / actual[idx_nonzero]
  
  # if y_i == 0 and yhat_i != 0: 1
  alpha[idx_zero & predicted != 0] <- 1
  
  # if y_i == 0 and yhat_i == 0: 0  (already 0 by initialization)
  alpha
}

calculate_mape <- function(actual, predicted) {
  alpha <- calc_alpha(actual, predicted)
  100 * mean(abs(alpha), na.rm = TRUE)   # = (100/N) * sum |alpha_i|
}

calculate_rrmse <- function(actual, predicted) {
  alpha <- calc_alpha(actual, predicted)
  sqrt(mean(alpha^2, na.rm = TRUE))      # = sqrt((1/N) * sum alpha_i^2)
}

calculate_maxmape <- function(actual, predicted) {
  alpha <- calc_alpha(actual, predicted)
  100 * max(abs(alpha), na.rm = TRUE)
}


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript calculate_mape_from_qs2.R <reference_file.qs2> <predicted_file.qs2>")
}

reference_file <- args[1]
predicted_file <- args[2]

# Check if files exist
if (!file.exists(reference_file)) {
  stop("Reference file not found: ", reference_file)
}
if (!file.exists(predicted_file)) {
  stop("Predicted file not found: ", predicted_file)
}

cat("Reading reference file:", reference_file, "\n")
ref_obj <- qs2::qs_read(reference_file)

cat("Reading predicted file:", predicted_file, "\n")
pred_obj <- qs2::qs_read(predicted_file)

# Extract C data frames
ref_data <- ref_obj$C
pred_data <- pred_obj$C

# Species to calculate MAPE for
species_list <- c("H", "O", "Charge", "C", "Ca", "Cl", "Mg", 
                  "Calcite", "Dolomite")

# Calculate MAPE and RRMSE for each species
results <- data.frame(
  Species = character(),
  MAPE = numeric(),
  RRMSE = numeric(),
  MaxMAPE = numeric(),
  stringsAsFactors = FALSE
)

cat("\nCalculating metrics for species:\n")

for (species in species_list) {
  if (!(species %in% names(ref_data)) || !(species %in% names(pred_data))) {
    warning("Species '", species, "' not found in data. Skipping.")
    next
  }
  
  mape <- calculate_mape(ref_data[[species]], pred_data[[species]])
  rrmse <- calculate_rrmse(ref_data[[species]], pred_data[[species]])
  maxmape <- calculate_maxmape(ref_data[[species]], pred_data[[species]])
  
  results <- rbind(results, data.frame(
    Species = species,
    MAPE = mape,
    RRMSE = rrmse,
    MaxMAPE = maxmape,
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  %-15s MAPE: %.2e  RRMSE: %.2e  MaxMAPE: %.2e\n", species, mape, rrmse, maxmape))
}

# Print summary
cat("\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 70), "\n", sep = "")
cat(sprintf("%-15s %-20s %-20s %-20s\n",
            "Species", "MAPE", "RRMSE", "MAXMAPE"))

cat("-" , rep("-", 70), "\n", sep = "")

for (i in 1:nrow(results)) {
  cat(sprintf("%-15s %-20.2e %-20.2e %-20.2e\n", 
              results$Species[i], 
              results$MAPE[i], 
              results$RRMSE[i],
              results$MaxMAPE[i]))
}

# Save results to CSV
output_file <- "mape_results.csv"
write.csv(results, output_file, row.names = FALSE)
cat("\nResults saved to:", output_file, "\n")
