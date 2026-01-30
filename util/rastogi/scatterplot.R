#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(qs2)})
library(ggplot2)

# read qs2 files and create scatterplot of predicted vs reference 

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
  stop("Usage: Rscript scatterplot.R <reference_file.qs2> <predicted_file.qs2>")
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
x_ref <- ref_obj$C[["Dolomite"]]
y_pred <- pred_obj$C[["Dolomite"]]

# Open PDF device
output_pdf <- "dolo_scatter.pdf"
pdf(output_pdf, width = 5.5, height = 5.5)

par(
  mar = c(5, 5, 4, 2),  # mehr Platz links
  las = 1,
  family = "serif"
)


df <- data.frame(x = x_ref, y = y_pred)

plot(
  x_ref, y_pred,
  asp = 1,
  pch = 4,
  cex = 1,          # Punkte größer
  lwd = 1,
  las = 1,            # y-Achse horizontal
  cex.axis = 0.9,     # Achsenwerte größer
  cex.lab  = 1,     # Achsentitel größer
  cex.main = 1.1,     # Titel größer
  family = "serif",
  xlab = "PHREEQC",
  ylab = "Prototyp 1",
  main = "ctrl_interval=100, mape_threshold=7e-3%, rb_limit=3",
  log = "xy"
)

abline(a = 0, b = 1, col = "red", lwd = 2)



mape  <- calculate_mape(x_ref, y_pred)
rrmse <- calculate_rrmse(x_ref, y_pred)

legend(
  "topleft",
  legend = sprintf(
    "MAPE  = %.2e %%\nRRMSE = %.2e",
    mape, rrmse
  ),
  bty = "n",
  text.font = 1,
  cex = 1.1
)


dev.off()

cat("Plot saved to:", output_pdf, "\n")
