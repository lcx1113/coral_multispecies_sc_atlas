#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1L]] else paste0(".", ""
)
script_path <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
))
workflow_root <- dirname(dirname(script_path))
sample <- "Pcru"
tag <- "pilot__raw_integer_full_pipeline"

suppressPackageStartupMessages(library(yaml))
source(file.path(workflow_root, "R", "common.R"))

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probability, names = FALSE))
}

config_path <- file.path(output_root, sample, "config", "config.yaml")
config <- yaml::read_yaml(config_path)
threshold <- as.numeric(config$resampling$stable_jaccard_threshold)
min_sampled <- as.integer(config$resampling$min_sampled_cells_per_cluster)
expected_iterations <- as.integer(config$resampling$pilot_iterations)
expected_resolutions <- sort(unique(as.numeric(unlist(config$resampling$pilot_resolutions))))
run_dir <- file.path(config$run_root, config$run_id)
resampling_dir <- file.path(run_dir, "03_resampling")
table_dir <- file.path(resampling_dir, "tables")
# Match the checkpoint directory created by 02_run_resampling.R for this Pilot set.
checkpoint_dir <- file.path(
  resampling_dir, "checkpoints",
  paste(tag, resolution_set_id(expected_resolutions), sep = "__")
)
preprocessing <- utils::read.csv(
  file.path(run_dir, "01_input_audit", "tables", "canonical_preprocessing_audit.csv"),
  stringsAsFactors = FALSE
)
selected_pcs <- as.integer(preprocessing$value[preprocessing$item == "selected_pcs"])
if (length(selected_pcs) != 1L || is.na(selected_pcs) || selected_pcs < 1L) {
  stop("The scan preprocessing audit must contain one valid selected_pcs value: ", sample)
}

paths <- list(
  jaccard = file.path(table_dir, paste0("cluster_jaccard__", tag, ".csv")),
  iteration = file.path(table_dir, paste0("iteration_summary__", tag, ".csv")),
  resolution = file.path(table_dir, paste0("resolution_summary__", tag, ".csv")),
  cluster = file.path(table_dir, paste0("cluster_stability__", tag, ".csv")),
  assignments = file.path(table_dir, paste0("cluster_assignments__", tag, ".tsv.gz")),
  reference = file.path(
    run_dir, "02_full_resolution_scan", "tables", "full_resolution_assignments.tsv.gz"
  ),
  manifest = file.path(resampling_dir, paste0("stage_manifest__", tag, ".csv")),
  checksums = file.path(resampling_dir, paste0("file_checksums__", tag, ".csv"))
)
stopifnot(all(vapply(paths, file.exists, logical(1))))

jaccard <- utils::read.csv(paths$jaccard, stringsAsFactors = FALSE)
iteration <- utils::read.csv(paths$iteration, stringsAsFactors = FALSE)
resolution <- utils::read.csv(paths$resolution, stringsAsFactors = FALSE)
cluster <- utils::read.csv(paths$cluster, stringsAsFactors = FALSE)
reference <- utils::read.delim(gzfile(paths$reference), check.names = FALSE)
manifest <- utils::read.csv(paths$manifest, stringsAsFactors = FALSE)
checksums <- utils::read.csv(paths$checksums, stringsAsFactors = FALSE)

jaccard$eligible <- jaccard$reference_cluster_sampled_n >= min_sampled
jaccard$stable <- jaccard$eligible & jaccard$max_jaccard >= threshold
jaccard$unstable <- jaccard$eligible & !jaccard$stable

iteration_exact <- lapply(
  split(jaccard, interaction(jaccard$resolution, jaccard$iteration, drop = TRUE)),
  function(x) {
    eligible <- x$eligible
    stable <- x$stable
    unstable <- x$unstable
    data.frame(
      resolution = x$resolution[[1L]],
      iteration = x$iteration[[1L]],
      eligible_cluster_count = sum(eligible),
      stable_cluster_count = sum(stable),
      unstable_cluster_count = sum(unstable),
      eligible_sampled_cell_count = sum(x$reference_cluster_sampled_n[eligible]),
      unstable_sampled_cell_count = sum(x$reference_cluster_sampled_n[unstable]),
      unstable_sampled_cell_fraction = if (sum(x$reference_cluster_sampled_n[eligible])) {
        sum(x$reference_cluster_sampled_n[unstable]) /
          sum(x$reference_cluster_sampled_n[eligible])
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  }
)
iteration_exact <- do.call(rbind, iteration_exact)
row.names(iteration_exact) <- NULL
iteration <- merge(
  iteration, iteration_exact,
  by = c("resolution", "iteration"), all.x = TRUE, sort = FALSE
)
iteration <- iteration[order(iteration$resolution, iteration$iteration), ]

resolution_label <- function(value) {
  sub("_0$", "", gsub(".", "_", sprintf("%.1f", value), fixed = TRUE))
}
size_rows <- lapply(sort(unique(cluster$resolution)), function(value) {
  column <- paste0("resolution_", resolution_label(value))
  sizes <- table(as.character(reference[[column]]))
  data.frame(
    resolution = value,
    reference_cluster = names(sizes),
    full_cluster_size = as.integer(sizes),
    stringsAsFactors = FALSE
  )
})
sizes <- do.call(rbind, size_rows)
cluster$reference_cluster <- as.character(cluster$reference_cluster)
cluster <- merge(
  cluster, sizes,
  by = c("resolution", "reference_cluster"), all.x = TRUE, sort = FALSE
)
cluster$resample_supported <- is.finite(cluster$median_jaccard) &
  cluster$median_jaccard >= threshold
cluster$resample_status <- ifelse(
  !is.finite(cluster$median_jaccard), "insufficient",
  ifelse(cluster$resample_supported, "supported", "unsupported")
)
cluster <- cluster[order(cluster$resolution, as.numeric(cluster$reference_cluster)), ]

iteration_by_resolution <- lapply(split(iteration, iteration$resolution), function(x) {
  data.frame(
    resolution = x$resolution[[1L]],
    median_eligible_cluster_count = stats::median(x$eligible_cluster_count),
    median_unstable_cluster_count = stats::median(x$unstable_cluster_count),
    q90_unstable_cluster_count = safe_quantile(x$unstable_cluster_count, 0.9),
    max_unstable_cluster_count = max(x$unstable_cluster_count),
    median_unstable_sampled_cell_count = stats::median(x$unstable_sampled_cell_count),
    q90_unstable_sampled_cell_count = safe_quantile(x$unstable_sampled_cell_count, 0.9),
    median_unstable_sampled_cell_fraction = stats::median(
      x$unstable_sampled_cell_fraction, na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
})
iteration_by_resolution <- do.call(rbind, iteration_by_resolution)

cluster_by_resolution <- lapply(split(cluster, cluster$resolution), function(x) {
  available <- x$resample_status != "insufficient"
  unsupported <- x$resample_status == "unsupported"
  total_cells <- sum(x$full_cluster_size, na.rm = TRUE)
  data.frame(
    resolution = x$resolution[[1L]],
    reference_cluster_count = nrow(x),
    supported_reference_cluster_count = sum(x$resample_status == "supported"),
    unsupported_reference_cluster_count = sum(unsupported),
    insufficient_reference_cluster_count = sum(!available),
    supported_reference_cluster_fraction = if (sum(available)) {
      sum(x$resample_status == "supported") / sum(available)
    } else NA_real_,
    unsupported_full_cell_count = sum(x$full_cluster_size[unsupported], na.rm = TRUE),
    unsupported_full_cell_fraction = if (total_cells) {
      sum(x$full_cluster_size[unsupported], na.rm = TRUE) / total_cells
    } else NA_real_,
    min_cluster_median_jaccard = min(x$median_jaccard, na.rm = TRUE),
    q10_cluster_level_median_jaccard = safe_quantile(x$median_jaccard, 0.1),
    unsupported_reference_clusters = paste(x$reference_cluster[unsupported], collapse = ";"),
    all_reference_clusters_supported = !any(unsupported) && all(available),
    stringsAsFactors = FALSE
  )
})
cluster_by_resolution <- do.call(rbind, cluster_by_resolution)
resolution <- merge(
  resolution, iteration_by_resolution, by = "resolution", all.x = TRUE, sort = FALSE
)
resolution <- merge(
  resolution, cluster_by_resolution, by = "resolution", all.x = TRUE, sort = FALSE
)
resolution <- resolution[order(resolution$resolution), ]

complete_count <- length(list.files(
  checkpoint_dir, pattern = "iteration_[0-9]+_complete[.]txt$", recursive = FALSE
))
jaccard_count <- length(list.files(
  checkpoint_dir, pattern = "iteration_[0-9]+_jaccard[.]csv$", recursive = FALSE
))
assignment_count <- length(list.files(
  checkpoint_dir, pattern = "iteration_[0-9]+_assignments[.]tsv[.]gz$", recursive = FALSE
))
summary_count <- length(list.files(
  checkpoint_dir, pattern = "iteration_[0-9]+_summary[.]csv$", recursive = FALSE
))
checksum_ok <- vapply(seq_len(nrow(checksums)), function(i) {
  identical(unname(tools::md5sum(checksums$path[[i]])), checksums$md5[[i]])
}, logical(1))
pair_outputs <- list.files(file.path(run_dir, "04_pair_consensus"), recursive = TRUE)
sensitivity_outputs <- list.files(file.path(run_dir, "05_parameter_sensitivity"), recursive = TRUE)
rds_outputs <- list.files(run_dir, pattern = "[.]rds$", recursive = TRUE)

integrity_pass <-
  nrow(iteration) == length(expected_resolutions) * expected_iterations &&
  nrow(resolution) == length(expected_resolutions) &&
  all(abs(sort(resolution$resolution) - sort(expected_resolutions)) < 1e-8) &&
  all(resolution$n_iterations == expected_iterations) &&
  length(unique(iteration$iteration)) == expected_iterations &&
  complete_count == expected_iterations &&
  jaccard_count == expected_iterations &&
  assignment_count == expected_iterations &&
  summary_count == expected_iterations &&
  all(checksum_ok) && !length(pair_outputs) && !length(sensitivity_outputs) &&
  !length(rds_outputs)

status <- data.frame(
  species = sample,
  run_id = config$run_id,
  stage_completed_at = manifest$completed_at[[1L]],
  selected_pcs = selected_pcs,
  workers = as.integer(config$compute$workers),
  resolution_min = min(expected_resolutions),
  resolution_max = max(expected_resolutions),
  n_resolutions = length(expected_resolutions),
  expected_iterations = expected_iterations,
  complete_checkpoints = complete_count,
  jaccard_checkpoints = jaccard_count,
  assignment_checkpoints = assignment_count,
  summary_checkpoints = summary_count,
  iteration_summary_rows = nrow(iteration),
  cluster_stability_rows = nrow(cluster),
  checksum_passed = sum(checksum_ok),
  checksum_total = length(checksum_ok),
  marker_mode = "none",
  pair_consensus = FALSE,
  sensitivity = FALSE,
  full_seurat_saved = FALSE,
  integrity_status = if (integrity_pass) "passed" else "failed",
  stringsAsFactors = FALSE
)

summary_dir <- file.path(output_root, "summary")
utils::write.csv(
  status, file.path(summary_dir, "Pcru_structural_pilot_run_status.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  iteration, file.path(summary_dir, "Pcru_structural_pilot_iteration_evidence.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  cluster, file.path(summary_dir, "Pcru_structural_pilot_cluster_evidence.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  resolution, file.path(summary_dir, "Pcru_structural_pilot_resolution_evidence.csv"),
  row.names = FALSE, na = ""
)

all_supported <- resolution$resolution[resolution$all_reference_clusters_supported]
min_unsupported <- min(resolution$unsupported_reference_cluster_count)
min_unsupported_resolutions <- resolution$resolution[
  resolution$unsupported_reference_cluster_count == min_unsupported
]
min_cell_fraction <- min(resolution$unsupported_full_cell_fraction)
min_cell_resolutions <- resolution$resolution[
  abs(resolution$unsupported_full_cell_fraction - min_cell_fraction) < 1e-12
]
fmt_values <- function(x) {
  if (!length(x)) return("none")
  paste(format(x, scientific = FALSE, trim = TRUE), collapse = ", ")
}

record <- c(
  "# Pcru structural Pilot results",
  "",
  "## Current output integrity",
  "",
  sprintf("- Integrity status: %s.", status$integrity_status[[1L]]),
  sprintf("- Pilot stage manifest completion time: %s.", status$stage_completed_at[[1L]]),
  sprintf("- Parameters: %d candidate resolutions (%s to %s), %d Pilot iterations, %d selected PCs (scan audit), subsampling fraction %s, k=%d, Leiden algorithm %d, %d workers (YAML).",
          length(expected_resolutions), min(expected_resolutions), max(expected_resolutions),
          expected_iterations, selected_pcs, config$resampling$fraction,
          config$resampling$k, config$clustering$algorithm, config$compute$workers),
  "- Batch start time and elapsed time are not inferred. This summary describes current scan/Pilot outputs and does not add Pcru to the nine-species display.",
  "",
  "## Integrity checks",
  "",
  "Only the matching Pilot/count-mode/candidate-set checkpoint directory is checked; Formal checkpoints are excluded.",
  sprintf("- Complete/Jaccard/summary/assignment checkpoints: %d/%d/%d/%d (expected %d each).", complete_count, jaccard_count, summary_count, assignment_count, expected_iterations),
  sprintf("- Iteration-resolution rows: %d (expected %d).", nrow(iteration), length(expected_resolutions) * expected_iterations),
  sprintf("- Resolution summary: %d rows; expected iterations per resolution: %d.", nrow(resolution), expected_iterations),
  sprintf("- Cluster-level rows: %d.", nrow(cluster)),
  sprintf("- File checksums: %d/%d MD5 matches.", sum(checksum_ok), length(checksum_ok)),
  sprintf("- Pair-consensus/sensitivity/RDS output counts: %d/%d/%d (expected zero for this structural-only workflow).", length(pair_outputs), length(sensitivity_outputs), length(rds_outputs)),
  "",
  "## Interpretation",
  "",
  "Each subsample independently rebuilds SCT, PCA, SNN and Leiden partitions. Each full-data reference cluster independently matches the resampled cluster with the largest Jaccard: |C intersect D| / |C union D|.",
  sprintf("At least %d sampled cells are required for a cluster to enter the statistics. A best-match Jaccard >= %s is stable in that iteration. Cluster-level support uses median_jaccard >= %s across %d Pilot iterations. This describes recovery of the current full-data partition, not a permanent property of a resolution.", min_sampled, threshold, threshold, expected_iterations),
  "",
  "## Resolution coverage",
  "",
  sprintf("- Among %d candidate resolutions, those with all reference clusters supported: %s.", length(expected_resolutions), fmt_values(all_supported)),
  sprintf("- Minimum unsupported reference-cluster count: %d, at resolutions %s.", min_unsupported, fmt_values(min_unsupported_resolutions)),
  sprintf("- Minimum unsupported full-cell fraction: %.4f, at resolutions %s.", min_cell_fraction, fmt_values(min_cell_resolutions)),
  "",
  "These summaries do not automatically select one resolution. Manual selection should consider contiguous ranges together with unsupported cluster counts, unsupported full-cell fraction, lower-tail Jaccard, stable-cell fraction and ARI/NMI.",
  "",
  "## Outputs",
  "",
  "- summary/Pcru_structural_pilot_run_status.csv: current output integrity and scan/Pilot records.",
  "- summary/Pcru_structural_pilot_resolution_evidence.csv: joint structural evidence per resolution.",
  "- summary/Pcru_structural_pilot_cluster_evidence.csv: full-data reference-cluster support across Pilot iterations.",
  "- summary/Pcru_structural_pilot_iteration_evidence.csv: per-iteration metrics and unstable cluster/cell counts.",
  sprintf("- Workflow tables: %s.", table_dir),
  sprintf("- Workflow figures: %s.", file.path(resampling_dir, "figures"))
)
writeLines(record, file.path(output_root, sample, "docs", "STRUCTURAL_PILOT_RUN_RECORD.md"))

code_files <- c(
  script_path,
  file.path(workflow_root, "scripts", "02_run_resampling.R"),
  file.path(workflow_root, "R", "common.R"),
  file.path(workflow_root, "R", "metrics.R"),
  config_path
)
code_manifest <- data.frame(
  path = code_files,
  md5 = unname(tools::md5sum(code_files)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  code_manifest,
  file.path(output_root, sample, "docs", "STRUCTURAL_PILOT_CODE_MD5.csv"),
  row.names = FALSE
)

if (!integrity_pass) stop("Pcru structural Pilot integrity audit failed.")
message("Pcru structural Pilot summary and integrity audit completed.")
