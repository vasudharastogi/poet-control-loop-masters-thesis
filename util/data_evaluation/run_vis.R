# Load the new functions
source("/mnt/beegfs/home/rastogi/poet/util/data_evaluation/RFun_Eval.R")

# Set base path
base_dir <- "/mnt/beegfs/home/rastogi/poet/bin"

sim1 <- ReadRTSims(file.path(base_dir, "proto2_eps01_no_rb_v2"))
sim2 <- ReadRTSims(file.path(base_dir, "proto2_eps0035_no_rb_v2"))


# ========================================
# Compare two simulations
# ========================================
rmse_errors <- ComputeErrors(sim1, sim2, FUN = RMSE)
mape_errors <- ComputeErrors(sim1, sim2, FUN = MAPE)

# Print summary
cat("RMSE errors computed for", ncol(rmse_errors), "variables\n")
cat("MAPE errors computed for", ncol(mape_errors), "variables\n")
cat("Number of iterations compared:", nrow(rmse_errors), "\n\n")

# Set output path explicitly
output_pdf <- file.path(base_dir, "comparison_plots.pdf")
cat("Saving plots to:", output_pdf, "\n")

# Save plots to PDF
pdf(output_pdf, width = 10, height = 6)

# Plot error progression
cat("Creating ErrorProgress plot...\n")
ErrorProgress(rmse_errors, ignore = c("Charge"), metric = "RMSE")

# Scatter plot for specific iteration (if that iteration exists)
if (length(sim1) >= 10 && length(sim2) >= 10) {
  cat("Creating scatter plot for iteration 10...\n")
  PlotScatter(sim1[[10]], sim2[[10]], 
              labs = c("Proto2 Eps=0.01%", "Proto2 Eps=0.0035%"),
              which = c("Ca", "Mg", "Cl", "C"))
} else {
  cat("Not enough iterations for scatter plot\n")
}

dev.off()
cat("Plots saved successfully to:", output_pdf, "\n")

# ========================================
# DHT/PHT analysis (if you have snapshot files)
# ========================================
#dht_snaps <- ReadAllDHT("poet/bin/proto1_only_interp")
#pht_snaps <- ReadAllPHT("poet/bin/proto1_only_interp")