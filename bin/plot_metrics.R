#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript plot_mape_stats.R <stats_overview_file1> [stats_overview_file2] ...")

cat("Reading", length(args), "stats file(s)...\n")

# Process all input files
all_data <- lapply(args, function(stats_file) {
  if (!file.exists(stats_file)) {
    warning("File not found: ", stats_file)
    return(NULL)
  }
  
  cat("  -", stats_file, "\n")
  
  lines <- readLines(stats_file)
  data_lines <- lines[!grepl("^-+$", lines) & nchar(lines) > 0]
  
  parsed <- lapply(data_lines, function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) >= 5) {
      data.frame(
        Iteration = as.numeric(parts[1]),
        Rollback = as.numeric(parts[2]),
        Species = parts[3],
        MAPE = as.numeric(parts[4]),
        RRMSE = as.numeric(parts[5]),
        stringsAsFactors = FALSE
      )
    }
  })
  
  df <- bind_rows(parsed) %>% filter(!is.na(Iteration))
  
  species_list <- c("H", "O", "C", "Ca", "Cl", "Mg", "Calcite", "Dolomite")
  #species_list <- "Dolomite"
  df_filtered <- df %>% 
    filter(Species %in% species_list) %>%
    group_by(Iteration) %>%
    summarise(
      MedianMAPE = median(MAPE, na.rm = TRUE),
      MaxMAPE = max(MAPE, na.rm = TRUE),
      Rollback = first(Rollback),
      .groups = "drop"
    ) %>%
    filter(Iteration %% 100 == 0) %>%
    mutate(Folder = basename(dirname(stats_file)))
  
  # Detect rollback changes
  df_filtered <- df_filtered %>%
    arrange(Iteration) %>%
    mutate(RollbackChange = Rollback != lag(Rollback, default = first(Rollback)))
  
  return(df_filtered)
})

combined_data <- bind_rows(all_data) %>%
  filter(Iteration <= 3000) %>%
  filter(is.finite(MedianMAPE) & MedianMAPE > 0) %>%
  filter(is.finite(MaxMAPE) & MaxMAPE > 0)

# Identify rollback transitions for each folder
rollback_points <- combined_data %>%
  filter(RollbackChange == TRUE) %>%
  select(Folder, Iteration, Rollback)

cat("\nData summary:\n")
print(head(combined_data))
cat("\nLegend:", unique(combined_data$Folder), "\n")
cat("\nRollback transitions detected:\n")
print(rollback_points)

# A consistent style for both plots
pretty_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_line(color = "grey85"),
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks = element_line(colour = "black")
  )

# Determine nice log-scale breaks (1e-1, 1e-2, 1e-3, etc.)
log_breaks <- 10^seq(max(-6, floor(log10(min(combined_data$MedianMAPE, combined_data$MaxMAPE, na.rm = TRUE)))),
                     ceiling(log10(max(combined_data$MedianMAPE, combined_data$MaxMAPE, na.rm = TRUE))),
                     by = 1)

# Common log label formatter
log_labels <- function(x) sprintf("1e%d", log10(x))

# Plot Median MAPE
p1 <- ggplot(combined_data, aes(x = Iteration, y = MedianMAPE, color = Folder)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(data = rollback_points, aes(xintercept = Iteration, color = Folder), 
             linetype = "dashed", alpha = 0.6, linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, max(combined_data$Iteration), by = 1000)) +
  scale_y_log10(breaks = log_breaks, labels = log_labels) +
  labs(
    title = "Median MAPE Across H, O, C, Ca, Cl, Mg, Calcite, Dolomite",
    x = "Iteration",
    y = "Median MAPE",
    color = "Legend"
  ) +
  pretty_theme

# Plot Max MAPE
p2 <- ggplot(combined_data, aes(x = Iteration, y = MaxMAPE, color = Folder)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(data = rollback_points, aes(xintercept = Iteration, color = Folder), 
             linetype = "dashed", alpha = 0.6, linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, max(combined_data$Iteration), by = 1000)) +
  scale_y_log10(breaks = log_breaks, labels = log_labels, limits = c(1e-5, NA)) +
  labs(
    title = "Max MAPE Across H, O, C, Ca, Cl, Mg, Calcite, Dolomite",
    x = "Iteration",
    y = "Max MAPE",
    color = "Legend"
  ) +
  pretty_theme


# Save plots
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") script_dir <- getwd()
ggsave(file.path(script_dir, "median_mape.pdf"), p1, width = 10, height = 6)
ggsave(file.path(script_dir, "max_mape.pdf"), p2, width = 10, height = 6)

cat("\nPlots saved:\n")
cat("  -", file.path(script_dir, "median_mape.pdf"), "\n")
cat("  -", file.path(script_dir, "max_mape.pdf"), "\n")

# Also save data
write.csv(combined_data, file.path(script_dir, "mape_summary.csv"), row.names = FALSE)
cat("  -", file.path(script_dir, "mape_summary.csv"), "\n")
