#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1L]] else paste0(".", ""
)
script_path <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
))
workflow_root <- dirname(dirname(script_path))

suppressPackageStartupMessages(library(yaml))

expected <- list(
  Adig = c(1.4, 1.5, 1.6, 2.3, 2.4, 2.5, 2.9, 3.0, 3.1),
  Ahya = c(1.7, 1.8, 1.9, 2.7, 2.8, 2.9, 3.0, 3.1),
  Amur = c(1.5, 1.6, 1.7, 1.8, 3.1, 3.2, 3.3, 3.4),
  Fusp = c(1.4, 1.5, 1.6, 1.7, 2.2, 2.3, 2.4, 2.8, 2.9, 3.0),
  Gfas = c(1.4, 1.5, 1.8, 1.9, 2.0, 2.1, 3.0, 3.1, 3.2),
  Hcoe = c(1.5, 1.6, 1.7, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.9),
  Pdam = seq(1.2, 2.0, 0.1),
  Plut = c(seq(1.5, 2.0, 0.1), 2.5, seq(2.7, 3.0, 0.1)),
  Xesp = c(seq(2.2, 2.5, 0.1), seq(2.8, 3.2, 0.1)),
  Pcru = c(seq(1.4, 1.8, 0.1), seq(2.3, 2.6, 0.1))
)

resolution_column <- function(value) {
  label <- sub("_0$", "", gsub(".", "_", sprintf("%.1f", value), fixed = TRUE))
  paste0("resolution_", label)
}

rows <- lapply(names(expected), function(sample) {
  config_path <- file.path(output_root, sample, "config", "config.yaml")
  config <- yaml::read_yaml(config_path)
  configured <- sort(as.numeric(unlist(config$resampling$formal_resolutions)))
  intended <- sort(expected[[sample]])
  grid <- as.numeric(unlist(config$clustering$resolution_grid))
  pilot <- as.numeric(unlist(config$resampling$pilot_resolutions))
  run_dir <- file.path(config$run_root, config$run_id)
  reference_path <- file.path(
    run_dir, "02_full_resolution_scan", "tables", "full_resolution_assignments.tsv.gz"
  )
  reference_header <- utils::read.delim(
    gzfile(reference_path), nrows = 1L, check.names = FALSE
  )
  audit <- utils::read.csv(
    file.path(run_dir, "01_input_audit", "tables", "canonical_preprocessing_audit.csv"),
    stringsAsFactors = FALSE
  )
  selected_pcs <- as.integer(audit$value[audit$item == "selected_pcs"])
  checkpoint_files <- list.files(
    file.path(run_dir, "03_resampling", "checkpoints"),
    recursive = TRUE, full.names = TRUE
  )
  formal_complete <- sum(
    grepl("/formal__", checkpoint_files, fixed = TRUE) &
      grepl("^iteration_[0-9]+_complete[.]txt$", basename(checkpoint_files))
  )
  exact_match <- length(configured) == length(intended) &&
    max(abs(configured - intended)) < 1e-8
  in_grid <- all(vapply(configured, function(x) any(abs(x - grid) < 1e-8), logical(1)))
  in_pilot <- all(vapply(configured, function(x) any(abs(x - pilot) < 1e-8), logical(1)))
  reference_ok <- all(vapply(
    configured, function(x) resolution_column(x) %in% names(reference_header), logical(1)
  ))
  parameter_ok <-
    as.integer(config$resampling$formal_iterations) == 100L &&
    abs(as.numeric(config$resampling$fraction) - 0.8) < 1e-8 &&
    as.integer(config$resampling$k) == 20L &&
    as.integer(config$compute$workers) <= 6L &&
    as.integer(config$compute$max_workers) <= 6L &&
    isTRUE(config$resampling$save_assignments) &&
    !isTRUE(config$resampling$save_seurat_objects)
  passed <- exact_match && in_grid && in_pilot && reference_ok && parameter_ok &&
    length(selected_pcs) == 1L
  data.frame(
    species = sample,
    requested_name = if (sample == "Pcru") "Plut_n" else sample,
    config_path = normalizePath(config_path),
    config_md5 = unname(tools::md5sum(config_path)),
    run_id = config$run_id,
    n_resolutions = length(configured),
    resolution_min = min(configured),
    resolution_max = max(configured),
    resolutions = paste(format(configured, scientific = FALSE, trim = TRUE), collapse = ","),
    formal_iterations = as.integer(config$resampling$formal_iterations),
    iteration_resolution_combinations = length(configured) *
      as.integer(config$resampling$formal_iterations),
    subsample_fraction = as.numeric(config$resampling$fraction),
    selected_pcs = selected_pcs,
    resampling_k = as.integer(config$resampling$k),
    workers = as.integer(config$compute$workers),
    all_resolutions_in_scan = in_grid && reference_ok,
    all_resolutions_in_pilot = in_pilot,
    save_assignments = isTRUE(config$resampling$save_assignments),
    save_seurat_objects = isTRUE(config$resampling$save_seurat_objects),
    existing_formal_complete_checkpoints = max(0L, formal_complete),
    validation = if (passed) "passed" else "failed",
    stringsAsFactors = FALSE
  )
})
preflight <- do.call(rbind, rows)
row.names(preflight) <- NULL

summary_path <- file.path(output_root, "summary", "formal_subsample_preflight.csv")
utils::write.csv(preflight, summary_path, row.names = FALSE, na = "")

manifest <- preflight[c("species", "config_path", "config_md5")]
utils::write.csv(
  manifest,
  file.path(output_root, "manifests", "formal_subsample_config_md5.csv"),
  row.names = FALSE
)

code_files <- c(
  script_path,
  file.path(workflow_root, "scripts", "12_run_260811_formal_subsample.sh"),
  file.path(workflow_root, "run_cluster_stability.R"),
  file.path(workflow_root, "scripts", "02_run_resampling.R"),
  file.path(workflow_root, "R", "metrics.R")
)
code_manifest <- data.frame(
  path = code_files,
  md5 = unname(tools::md5sum(code_files)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  code_manifest,
  file.path(output_root, "manifests", "formal_subsample_code_md5.csv"),
  row.names = FALSE
)

if (any(preflight$validation != "passed")) {
  print(preflight[preflight$validation != "passed", ], row.names = FALSE)
  stop("One or more formal subsample configurations failed validation.")
}
message(sprintf(
  "Formal preflight passed: %d samples, %d resolution points, %d iteration-resolution combinations.",
  nrow(preflight), sum(preflight$n_resolutions),
  sum(preflight$iteration_resolution_combinations)
))
