#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

# Plot MAPE and rollbacks of one metrics_overview file

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript plot_mape_stats_one.R <stats_overview_file>")

stats_file <- args[1]
if (!file.exists(stats_file)) stop("File not found: ", stats_file)

cat("Reading:", stats_file, "\n")

# ----------------------------
# Read + parse ONE input file
# ----------------------------
lines <- readLines(stats_file)
data_lines <- lines[!grepl("^-+$", lines) & nchar(lines) > 0]

parsed <- lapply(data_lines, function(line) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  if (length(parts) >= 5) {
    data.frame(
      Iteration = as.numeric(parts[1]),
      Species   = parts[3],
      MAPE      = as.numeric(parts[4]),
      RRMSE     = as.numeric(parts[5]),
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }
})

df <- bind_rows(parsed) %>% filter(!is.na(Iteration))

species_list <- c("H", "O", "C", "Ca", "Cl", "Mg", "Calcite", "Dolomite")

summary_df <- df %>%
  filter(Species %in% species_list) %>%
  group_by(Iteration) %>%
  summarise(
    MedianMAPE = median(MAPE, na.rm = TRUE),
    MaxMAPE    = max(MAPE, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  filter(Iteration <= 10000) %>%
  filter(is.finite(MedianMAPE) & MedianMAPE > 0) %>%
  filter(is.finite(MaxMAPE) & MaxMAPE > 0) %>%
  arrange(Iteration)

cat("\nData summary:\n")
print(head(summary_df))

# ----------------------------
# Helper: log-axis tick marks
# ----------------------------
log_breaks <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(numeric())
  10^seq(floor(log10(min(x))), ceiling(log10(max(x))), by = 1)
}

median_log_breaks <- log_breaks(summary_df$MedianMAPE)
max_log_breaks    <- log_breaks(summary_df$MaxMAPE)

# ----------------------------
# Output directory
# ----------------------------
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- getwd()

# ----------------------------
# Base R plot: Median MAPE
# ----------------------------
median_pdf <- file.path(script_dir, "median_mape.pdf")
pdf(median_pdf, width = 6, height = 6)

xlim <- c(0, 5000)
ylim <- c(1e-14, 1e-3)

plot(summary_df$Iteration, summary_df$MedianMAPE,
     type="n", log="y", xlim=xlim, ylim=ylim,
     xlab="Iteration", ylab="Median MAPE",
     main="Median MAPE ...",
     family="serif", cex.axis=0.75, cex.lab=0.85)

# NOW shade using those exact limits (this always works)
rect(0,     ylim[1], 100,  ylim[2],
     col=adjustcolor("green", 0.40), border=NA)

rect(4200,  ylim[1], xlim[2], ylim[2],
     col=adjustcolor("orange", 0.40), border=NA)

lines(summary_df$Iteration, summary_df$MedianMAPE, lwd=2)
points(summary_df$Iteration, summary_df$MedianMAPE, pch=19, cex=0.7)
box()
dev.off()
# ----------------------------
# Base R plot: Max MAPE (with shaded phases)
# ----------------------------
max_pdf <- file.path(script_dir, "max_mape.pdf")
pdf(max_pdf, width = 7, height = 5)


xlim <- c(0, 4700)

# constrain ylim via log breaks
ylim <- range(max_log_breaks, finite = TRUE)
ylim[2] <- ylim[2] * 1.05

plot(summary_df$Iteration, summary_df$MaxMAPE,
     type="n", log="y", xlim=xlim, ylim=ylim,
     main="Maximaler MAPE über die Iterationen",
     xlab="Iteration", ylab="MAPE",
     family="serif", cex.axis=0.9,
     yaxt="n", xaxt="n")

# custom axes
axis(1, at = seq(xlim[1], xlim[2], by = 1000), family="serif", cex.axis=0.9)
axis(2, at = max_log_breaks,
     labels = format(max_log_breaks, scientific = TRUE),
     family = "serif", cex.axis = 0.9)

# shaded phases
rect(0,    ylim[1], 100,  ylim[2],
     col=adjustcolor("darkgreen", 0.40), border=NA)
rect(4200, ylim[1], xlim[2], ylim[2],
     col=adjustcolor("darkorange", 0.40), border=NA)

rollback_iters <- c(2300, 2500, 2800, 4000, 4100, 4200)

# vertical lines
abline(v = rollback_iters,
       col = "red",
       lty = 6,
       lwd = 1.5)

# data
lines(summary_df$Iteration, summary_df$MaxMAPE, lwd=1.5)
points(summary_df$Iteration, summary_df$MaxMAPE, pch=16, cex=0.7)

abline(h=3.5e-3, col="red", lwd=2)

legend(
  x = "bottom",
  bty = "n",
  text.font = 1,
  cex = 0.8,
  fill = c(
    adjustcolor("darkgreen",  alpha.f = 0.4),
    adjustcolor("darkorange", alpha.f = 0.4)
  ),
  legend = c("Aufwärmphase", "Deaktiviertes Ersatzmodell")
)

#box(lwd=0.3)
dev.off()




cat("\nPlots saved:\n")
cat("  -", median_pdf, "\n")
cat("  -", max_pdf, "\n")

# ----------------------------
# Save data
# ----------------------------
out_csv <- file.path(script_dir, "mape_summary.csv")
write.csv(summary_df, out_csv, row.names = FALSE)
cat("  -", out_csv, "\n")
