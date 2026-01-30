#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(qs2)
})
# Compute MAXMAPE time series across multiple predicted qs2 directories

# ----------------------------
# MAPE (paper-style alpha_i)
# ----------------------------
calc_alpha <- function(actual, predicted) {
  stopifnot(length(actual) == length(predicted))

  alpha <- numeric(length(actual))
  idx_nonzero <- actual != 0
  idx_zero    <- actual == 0

  alpha[idx_nonzero] <- (actual[idx_nonzero] - predicted[idx_nonzero]) / actual[idx_nonzero]
  alpha[idx_zero & predicted != 0] <- 1
  alpha
}

calculate_mape <- function(actual, predicted) {
  alpha <- calc_alpha(actual, predicted)
  100 * mean(abs(alpha), na.rm = TRUE)
}

# ----------------------------
# Helpers
# ----------------------------
extract_iter <- function(path) {
  bn <- basename(path)
  m <- regexpr("[0-9]+", bn)
  if (m[1] == -1) return(NA_integer_)
  as.integer(regmatches(bn, m))
}

list_qs2_files <- function(x) {
  if (dir.exists(x)) {
    list.files(x, pattern = "\\.qs2$", full.names = TRUE)
  } else {
    Sys.glob(x)
  }
}

safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

# Compute MAXMAPE series for one predicted input vs reference
compute_series_maxmape <- function(ref_in, pred_in, species_list) {
  ref_files  <- list_qs2_files(ref_in)
  pred_files <- list_qs2_files(pred_in)

  if (length(ref_files) == 0)  stop("No reference .qs2 files found for: ", ref_in)
  if (length(pred_files) == 0) stop("No predicted .qs2 files found for: ", pred_in)

  ref_tbl  <- data.frame(file_ref  = ref_files,  iter = vapply(ref_files,  extract_iter, integer(1)))
  pred_tbl <- data.frame(file_pred = pred_files, iter = vapply(pred_files, extract_iter, integer(1)))

  ref_tbl  <- ref_tbl[!is.na(ref_tbl$iter), ]
  pred_tbl <- pred_tbl[!is.na(pred_tbl$iter), ]

  pairs <- merge(ref_tbl, pred_tbl, by = "iter")
  pairs <- pairs[order(pairs$iter), ]
  if (nrow(pairs) == 0) stop("No matching iterations between ref and pred: ", pred_in)

  series <- data.frame(iter = integer(), maxmape = numeric(), stringsAsFactors = FALSE)

  for (i in seq_len(nrow(pairs))) {
    iter <- pairs$iter[i]
    ref_obj  <- qs2::qs_read(pairs$file_ref[i])
    pred_obj <- qs2::qs_read(pairs$file_pred[i])

    ref_C  <- ref_obj$C
    pred_C <- pred_obj$C

    mape_vec <- vapply(species_list, function(sp) {
      if (!(sp %in% names(ref_C)) || !(sp %in% names(pred_C))) return(NA_real_)
      calculate_mape(ref_C[[sp]], pred_C[[sp]])
    }, numeric(1))

    series <- rbind(series, data.frame(iter = iter, maxmape = safe_max(mape_vec)))
  }

  series
}

# ----------------------------
# Args
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "  Rscript multi_maxmape_timeseries.R <phreeqc_ref_dir_or_glob> <reference_pred_dir_or_glob> <pred1_dir_or_glob> [pred2 ...]",
    "",
    "Meaning:",
    "  1) phreeqc_ref      = PHREEQC reference (denominator for MAXMAPE for ALL lines)",
    "  2) reference_pred   = special prediction to plot in RED (compared vs phreeqc_ref)",
    "  3+) preds           = other predictions to plot in DARKGREEN (compared vs phreeqc_ref)",
    "",
    "Example:",
    "  Rscript multi_maxmape_timeseries.R pqc m_ref m_p1_100_1_3_1 m_p1_100_2_3_1 m_p1_100_3_3_1",
    sep = "\n"
  ))
}

ref_in      <- args[1]          # PHREEQC reference
ref_pred_in <- args[2]          # special "reference prediction" (red line)
pred_list   <- args[-c(1, 2)]   # all other predictions (green lines)

species_list <- c("H","O","C","Ca","Cl","Mg","Calcite","Dolomite")

# ----------------------------
# User mappings (optional)
# ----------------------------
# keys must match basename(<dir_or_glob>)
rename_map <- c(
  "dolo_fgcs_3_interp_skip" = "Interpolation (3, 3)", 
  
  "m_p1_200_1_3_1" = "(200, 1.75e-3)",
  "m_p1_200_2_3_1" = "(200, 3.5e-3)",
  "m_p1_200_3_3_1" = "(200, 7e-3)"
  
  #"m_p1_100_1_3_1" = "(100, 1.75e-3)",
  #"m_p1_100_2_3_1" = "(100, 3.5e-3)",
  #"m_p1_100_3_3_1" = "(100, 7e-3)"

  #"m_p1_50_1_3_1" = "(50, 1.75e-3)",
  #"m_p1_50_2_3_1" = "(50, 3.5e-3)",
  #"m_p1_50_3_3_1" = "(50, 7e-3)"

  #"m_p2_50_1_3" = "(50, 3.5e-3)",
  #"m_p2_50_2_3" = "(50, 3.5e-2)",
  #"m_p2_50_3_3" = "(50, 3.5e-1)"

  #"m_p2_100_1_3" = "(100, 3.5e-3)",
  #"m_p2_100_2_3" = "(100, 3.5e-2)",
  #"m_p2_100_3_3" = "(100, 3.5e-1)"

  #"m_p2_100_1_3" = "(100, 3.5e-3)",
  #"m_p2_100_2_3" = "(100, 3.5e-2)",
  #"m_p2_100_3_3" = "(100, 3.5e-1)"

  #"m_p2_200_1_3" = "(200, 3.5e-3)",
  #"m_p2_200_2_3" = "(200, 3.5e-2)",
  #"m_p2_200_3_3" = "(200, 3.5e-1)"
)

# line types (base R lty: 1 solid, 2 dashed, 3 dotted, 4 dotdash, 5 longdash, 6 twodash)
lty_map <- c(
  # "m_ref"           = 1,  # (will be forced to 1 anyway)
  "m_p1_200_1_3_1" = 1,
  "m_p1_200_2_3_1" = 2,
  "m_p1_200_3_3_1" = 3

  #"m_p1_50_1_3_1" = 1,
  #"m_p1_50_2_3_1" = 2,
  #"m_p1_50_3_3_1" = 3

  #"m_p1_100_1_3_1" = 1,
  #"m_p1_100_2_3_1" = 2,
  #"m_p1_100_3_3_1" = 3

  #"m_p2_50_1_3" = 1,
  #"m_p2_50_2_3" = 2,
  #"m_p2_50_3_3" = 3

  #"m_p2_100_1_3" = 1,
  #"m_p2_100_2_3" = 2,
  #"m_p2_100_3_3" = 3

  #"m_p2_100_1_3" = 1,
  #"m_p2_100_2_3" = 2,   
  #"m_p2_100_3_3" = 3

  #"m_p2_200_1_3" = 1,
  #"m_p2_200_2_3" = 2,   
  #"m_p2_200_3_3" = 3
)

# ----------------------------
# Compute series
# ----------------------------
all_series <- list()

# Special reference curve (RED): ref_pred_in vs PHREEQC ref_in
cat("\nComputing REFERENCE:", ref_pred_in, "vs", ref_in, "\n")
s_ref <- compute_series_maxmape(ref_in, ref_pred_in, species_list)
key_ref <- basename(ref_pred_in)
s_ref$key   <- key_ref
s_ref$label <- if (key_ref %in% names(rename_map)) rename_map[[key_ref]] else key_ref
s_ref$is_ref <- TRUE
all_series[[key_ref]] <- s_ref
cat("  matched iterations:", nrow(s_ref), "range:", min(s_ref$iter), "to", max(s_ref$iter), "\n")

# All other curves (DARKGREEN): pred_in vs PHREEQC ref_in
for (pred_in in pred_list) {
  cat("\nComputing:", pred_in, "vs", ref_in, "\n")
  s <- compute_series_maxmape(ref_in, pred_in, species_list)

  key <- basename(pred_in)
  s$key   <- key
  s$label <- if (key %in% names(rename_map)) rename_map[[key]] else key
  s$is_ref <- FALSE

  all_series[[key]] <- s
  cat("  matched iterations:", nrow(s), "range:", min(s$iter), "to", max(s$iter), "\n")
}

combined <- do.call(rbind, all_series)
combined <- combined[order(combined$iter), ]

# Start plots at iteration 200 (skip warm-up)
combined <- combined[combined$iter >= 0, ]

# Determine plot limits
xlim <- range(combined$iter, finite = TRUE)
ylim <- c(1e-3, 1e+0)

# ----------------------------
# Plot (base R)
# ----------------------------
output_pdf <- "multi_maxmape.pdf"
pdf(output_pdf, width = 5, height = 4)

par(family = "serif",
mar = c(4, 4, 2, 2))  
plot(NA,
     xlim = xlim,
     ylim = ylim,          # <-- REQUIRED
     las  = 1,
     log  = "y",
     cex.axis = 0.9,
     cex.lab  = 1,
     xlab = "Iteration",
     ylab = "MAPE (%)",
     yaxt = "n")

# --- decade ticks: 1e-3, 1e-2, 1e-1, 1e0, 1e1, ...
y_ticks <- 10^seq(
  floor(log10(ylim[1])),
  ceiling(log10(ylim[2]))
)

axis(2,
     at = y_ticks,
     labels = formatC(y_ticks, format = "e", digits = 0),
     las = 1,
     cex.axis = 0.9)

keys <- unique(combined$key)

# stable mapping key -> label
key_to_label <- tapply(combined$label, combined$key, function(x) x[1])

# COLORS: reference always red, all others same dark green
#cols <- setNames(rep("orange", length(keys)), keys)
cols <- grDevices::hcl.colors(length(keys), palette = "Orange", rev = TRUE) 
names(cols) <- keys 

if (key_ref %in% keys) cols[key_ref] <- "red"

# LINE TYPES: from lty_map, default solid; reference forced solid
default_lty <- 1
lty_for_key <- setNames(rep(default_lty, length(keys)), keys)
for (k in keys) if (k %in% names(lty_map)) lty_for_key[k] <- lty_map[[k]]
if (key_ref %in% keys) lty_for_key[key_ref] <- 1

# draw lines (reference first)
plot_order <- c(intersect(key_ref, keys), setdiff(keys, key_ref))

for (k in plot_order) {
  df <- combined[combined$key == k, ]
  df <- df[order(df$iter), ]
  y <- df$maxmape
  y[!is.finite(y) | y <= 0] <- NA_real_
  lines(df$iter, y, lwd = 2, col = cols[k], lty = lty_for_key[k])
}

legend("topleft",
       legend = unname(key_to_label[plot_order]),
       col    = cols[plot_order],
       lwd    = 2,
       lty    = lty_for_key[plot_order],
       bty    = "n",
       cex    = 0.8)

legend("topright",
       legend = "Max MAPE",
       bty = "n",
       cex = 1)


box()
dev.off()
cat("\nSaved plot to:", output_pdf, "\n")
