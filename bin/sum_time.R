#!/usr/bin/env Rscript

# Summarize timing vectors from timings.qs2 files.
# Usage: Rscript summarize_timings.R file1.qs2 [file2.qs2 ...]

if (!requireNamespace("qs2", quietly = TRUE)) {
  stop("Package 'qs2' not installed. Install with: remotes::install_github('qs2io/qs2')")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript summarize_timings.R timings.qs2 [more.qs2 ...]")
}

chem_vectors <- c(
  "idle_worker", "phreeqc_time", "dht_get_time", "dht_fill_time",
  "interp_w", "interp_r", "interp_g", "interp_fc"
)

summaries <- lapply(args, function(f) {

  if (!file.exists(f)) {
    warning("File not found: ", f)
    return(NULL)
  }

  obj <- qs2::qs_read(f)
  chem <- obj$chemistry
  ctrl <- obj$control_loop

  # ---- sum chemistry vectors, round 2 digits ----
  chem_sums <- sapply(chem_vectors, function(v) {
    x <- chem[[v]]
    if (!is.numeric(x)) return(NA_real_)
    round(sum(x, na.rm = TRUE), 2)
  })

  # ---- sum worker ----
  worker_sum <- {
    x <- ctrl$worker
    if (!is.numeric(x)) NA_real_ else round(sum(x, na.rm = TRUE), 2)
  }

  # ---- assemble a long-format table ----
  data.frame(
    file = basename(f),
    vector = c(chem_vectors, "worker"),
    sum = c(chem_sums, worker_sum),
    stringsAsFactors = FALSE
  )
})

summaries <- do.call(rbind, summaries)

# ---- print result ----
cat("\nTiming sums (rounded to 2 digits):\n\n")
print(summaries, row.names = FALSE)

# ---- save CSV ----
write.csv(summaries, "timings_summary.csv", row.names = FALSE)
cat("\nSaved summary to timings_summary.csv\n")
