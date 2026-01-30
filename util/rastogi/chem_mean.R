#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(qs2)
})
#compute mean chemistry ($C) across multiple qs2 files and write to output qs2 file

# ----------------------------
# Helpers
# ----------------------------
extract_iter <- function(path) {
  bn <- basename(path)
  m <- regexpr("[0-9]+", bn)
  if (m[1] == -1) return(NA_integer_)
  as.integer(regmatches(bn, m))
}

list_qs2_files <- function(dir_or_glob) {
  if (dir.exists(dir_or_glob)) {
    list.files(dir_or_glob, pattern = "\\.qs2$", full.names = TRUE)
  } else {
    Sys.glob(dir_or_glob)
  }
}

mean_chemistry_files <- function(files, out_file) {
  stopifnot(length(files) >= 1)

  ref <- qs2::qs_read(files[1])
  if (!("C" %in% names(ref)))
    stop("No $C chemistry component in: ", files[1])

  chem_names <- names(ref$C)

  # accumulator
  sum_C <- ref$C
  for (v in chem_names) sum_C[[v]] <- as.numeric(sum_C[[v]])

  if (length(files) > 1) {
    for (f in files[-1]) {
      obj <- qs2::qs_read(f)
      if (!("C" %in% names(obj)))
        stop("No $C chemistry component in: ", f)

      if (!identical(names(obj$C), chem_names))
        stop("Chemistry variables differ in: ", f)

      for (v in chem_names) sum_C[[v]] <- sum_C[[v]] + obj$C[[v]]
    }
  }

  n <- length(files)
  mean_C <- sum_C
  for (v in chem_names) mean_C[[v]] <- mean_C[[v]] / n

  out <- ref
  out$C <- mean_C

  # optional: keep iteration metadata from template
  out$ai_surrogate_info <- NULL

  qs2::qs_save(out, out_file)
}

# ----------------------------
# Main: mean per iteration across folders
# ----------------------------
mean_by_iteration_across_folders <- function(input_dirs, output_dir, pattern = "\\.qs2$") {
  if (length(input_dirs) < 2)
    stop("Need at least 2 input folders to compute a mean.")

  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Build: for each folder -> map(iter -> filepath)
  maps <- vector("list", length(input_dirs))
  names(maps) <- input_dirs

iter_pattern <- "^iter_.*\\.qs2$"

for (i in seq_along(input_dirs)) {
  d <- input_dirs[i]
  files <- list_qs2_files(d)

  files <- files[grepl(iter_pattern, basename(files))]
  if (length(files) == 0)
    stop("No iter_*.qs2 files found in: ", d)

  iters <- vapply(files, extract_iter, integer(1))
  ok <- !is.na(iters)
  files <- files[ok]
  iters <- iters[ok]

  if (any(duplicated(iters))) {
    dup <- unique(iters[duplicated(iters)])
    stop("Duplicate iteration numbers in folder ", d, ": ", paste(dup, collapse = ", "))
  }

  maps[[i]] <- setNames(files, iters)
  message("Found ", length(files), " iteration files in: ", d)
}


  # intersection of iterations present in ALL folders
  iter_sets <- lapply(maps, names)
  common_iters <- Reduce(intersect, iter_sets)
  common_iters <- as.integer(common_iters)
  common_iters <- sort(common_iters)

  if (length(common_iters) == 0)
    stop("No common iteration numbers across all input folders.")

  message("Common iterations across all folders: ", length(common_iters))

  # Compute mean for each iteration
  for (it in common_iters) {
    files_for_it <- vapply(maps, function(m) m[[as.character(it)]], character(1))
    out_file <- file.path(output_dir, sprintf("iter_%d.qs2", it))

    message("Averaging iter ", it, " -> ", basename(out_file))
    mean_chemistry_files(files_for_it, out_file)
  }

  message("Done. Output written to: ", output_dir)
}

# ----------------------------
# CLI
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(paste(
    "Usage:",
    "  Rscript mean_chem.R <output_dir> <input_dir1> <input_dir2> [input_dir3 ...]",
    "",
    "What it does:",
    "  For each iteration number that exists in ALL input dirs, compute the mean chemistry ($C)",
    "  per cell and write iter_<N>.qs2 into <output_dir>.",
    "",
    "Example:",
    "  Rscript mean_chem.R mean_out runA runB runC",
    sep = "\n"
  ))
}

output_dir <- args[1]
input_dirs <- args[-1]

mean_by_iteration_across_folders(input_dirs, output_dir)
