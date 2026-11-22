#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1L]] else paste0(".", ""
)
script_path <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
))
workflow_root <- dirname(dirname(script_path))
species <- c("Adig", "Ahya", "Amur", "Fusp", "Gfas", "Hcoe", "Pdam", "Plut", "Xesp")
tag <- "pilot__raw_integer_full_pipeline"

suppressPackageStartupMessages(library(yaml))
source(file.path(workflow_root, "R", "common.R"))

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probability, names = FALSE))
}

iteration_rows <- list()
cluster_rows <- list()
resolution_rows <- list()
audit_rows <- list()

for (sample in species) {
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
  iteration_path <- file.path(table_dir, paste0("iteration_summary__", tag, ".csv"))
  resolution_path <- file.path(table_dir, paste0("resolution_summary__", tag, ".csv"))
  cluster_path <- file.path(table_dir, paste0("cluster_stability__", tag, ".csv"))
  assignment_path <- file.path(
    run_dir, "02_full_resolution_scan", "tables", "full_resolution_assignments.tsv.gz"
  )
  manifest_path <- file.path(resampling_dir, paste0("stage_manifest__", tag, ".csv"))

  iteration <- utils::read.csv(iteration_path, stringsAsFactors = FALSE)
  resolution <- utils::read.csv(resolution_path, stringsAsFactors = FALSE)
  cluster <- utils::read.csv(cluster_path, stringsAsFactors = FALSE)
  reference <- utils::read.delim(gzfile(assignment_path), check.names = FALSE)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)

  iteration$species <- sample
  iteration$stable_cluster_count <- as.integer(round(
    iteration$n_sufficient_clusters * iteration$stable_cluster_fraction
  ))
  iteration$unstable_cluster_count <-
    iteration$n_sufficient_clusters - iteration$stable_cluster_count
  iteration_rows[[sample]] <- iteration[c(
    "species", "resolution", "iteration", "seed", "n_sufficient_clusters",
    "stable_cluster_count", "unstable_cluster_count", "stable_cluster_fraction",
    "stable_cell_fraction", "median_cluster_jaccard", "q10_cluster_jaccard",
    "partition_ARI", "partition_NMI"
  )]

  size_rows <- lapply(sort(unique(cluster$resolution)), function(value) {
    column <- paste0("resolution_", gsub("\\.", "_", format(value, scientific = FALSE, trim = TRUE)))
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
  cluster$species <- sample
  cluster$insufficient_iteration_count <-
    cluster$n_iterations - cluster$n_sufficient_iterations
  cluster$stable_sufficient_iteration_count <- as.integer(round(
    cluster$n_sufficient_iterations * cluster$stable_iteration_fraction
  ))
  cluster$unstable_sufficient_iteration_count <-
    cluster$n_sufficient_iterations - cluster$stable_sufficient_iteration_count
  cluster$resample_cluster_support <- ifelse(
    is.finite(cluster$median_jaccard),
    cluster$median_jaccard >= threshold,
    NA
  )
  cluster$resample_status <- ifelse(
    is.na(cluster$resample_cluster_support), "insufficient",
    ifelse(cluster$resample_cluster_support, "supported", "unsupported")
  )
  cluster_rows[[sample]] <- cluster[c(
    "species", "resolution", "reference_cluster", "full_cluster_size",
    "n_iterations", "n_sufficient_iterations", "insufficient_iteration_count",
    "stable_sufficient_iteration_count", "unstable_sufficient_iteration_count",
    "median_sampled_cells", "median_jaccard", "q10_jaccard",
    "stable_iteration_fraction", "median_reference_recall",
    "median_test_precision", "median_match_margin", "resample_cluster_support",
    "resample_status"
  )]

  iteration_by_resolution <- lapply(split(iteration, iteration$resolution), function(x) {
    data.frame(
      resolution = x$resolution[[1L]],
      median_sufficient_cluster_count = stats::median(x$n_sufficient_clusters),
      median_stable_cluster_count = stats::median(x$stable_cluster_count),
      median_unstable_cluster_count = stats::median(x$unstable_cluster_count),
      q90_unstable_cluster_count = safe_quantile(x$unstable_cluster_count, 0.9),
      max_unstable_cluster_count = max(x$unstable_cluster_count),
      stringsAsFactors = FALSE
    )
  })
  iteration_by_resolution <- do.call(rbind, iteration_by_resolution)
  cluster_by_resolution <- lapply(split(cluster, cluster$resolution), function(x) {
    available <- !is.na(x$resample_cluster_support)
    unsupported <- available & !x$resample_cluster_support
    supported <- available & x$resample_cluster_support
    data.frame(
      resolution = x$resolution[[1L]],
      reference_cluster_count = nrow(x),
      resample_supported_cluster_count = sum(supported),
      resample_unsupported_cluster_count = sum(unsupported),
      insufficient_cluster_count = sum(!available),
      resample_supported_cluster_fraction = if (sum(available)) {
        sum(supported) / sum(available)
      } else NA_real_,
      unsupported_full_cell_fraction = if (sum(x$full_cluster_size, na.rm = TRUE) > 0) {
        sum(x$full_cluster_size[unsupported], na.rm = TRUE) /
          sum(x$full_cluster_size, na.rm = TRUE)
      } else NA_real_,
      min_cluster_median_jaccard = min(x$median_jaccard, na.rm = TRUE),
      q10_cluster_level_median_jaccard = safe_quantile(x$median_jaccard, 0.1),
      median_cluster_stable_iteration_fraction = stats::median(
        x$stable_iteration_fraction, na.rm = TRUE
      ),
      min_cluster_stable_iteration_fraction = min(
        x$stable_iteration_fraction, na.rm = TRUE
      ),
      unsupported_reference_clusters = paste(
        x$reference_cluster[unsupported], collapse = ";"
      ),
      all_reference_clusters_resample_supported =
        sum(unsupported) == 0L && sum(!available) == 0L,
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
  resolution$species <- sample
  resolution_rows[[sample]] <- resolution[c(
    "species", "resolution", "n_iterations", "reference_cluster_count",
    "median_sufficient_cluster_count", "median_stable_cluster_count",
    "median_unstable_cluster_count", "q90_unstable_cluster_count",
    "max_unstable_cluster_count", "median_cluster_jaccard",
    "q10_cluster_jaccard", "median_stable_cluster_fraction",
    "median_stable_cell_fraction", "median_partition_ARI", "q10_partition_ARI",
    "median_partition_NMI", "resample_supported_cluster_count",
    "resample_unsupported_cluster_count", "insufficient_cluster_count",
    "resample_supported_cluster_fraction", "unsupported_full_cell_fraction",
    "min_cluster_median_jaccard", "q10_cluster_level_median_jaccard",
    "median_cluster_stable_iteration_fraction",
    "min_cluster_stable_iteration_fraction", "unsupported_reference_clusters",
    "all_reference_clusters_resample_supported"
  )]

  complete_count <- length(list.files(
    checkpoint_dir, pattern = "iteration_[0-9]+_complete\\.txt$", recursive = FALSE
  ))
  jaccard_count <- length(list.files(
    checkpoint_dir, pattern = "iteration_[0-9]+_jaccard\\.csv$", recursive = FALSE
  ))
  assignment_count <- length(list.files(
    checkpoint_dir, pattern = "iteration_[0-9]+_assignments\\.tsv\\.gz$", recursive = FALSE
  ))
  expected_resolution_count <- length(expected_resolutions)
  expected_iteration_rows <- expected_resolution_count * expected_iterations
  pair_outputs <- list.files(resampling_dir, pattern = "pair_consensus|pair_", recursive = TRUE)
  sensitivity_outputs <- list.files(resampling_dir, pattern = "sensitivity", recursive = TRUE)
  rds_outputs <- list.files(run_dir, pattern = "\\.rds$", recursive = TRUE)
  passed <- complete_count == expected_iterations &&
    jaccard_count == expected_iterations && assignment_count == expected_iterations &&
    nrow(iteration) == expected_iteration_rows &&
    nrow(resolution) == expected_resolution_count &&
    all(abs(sort(resolution$resolution) - expected_resolutions) < 1e-8) &&
    all(resolution$n_iterations == expected_iterations) &&
    !length(pair_outputs) && !length(sensitivity_outputs) && !length(rds_outputs)
  audit_rows[[sample]] <- data.frame(
    species = sample,
    run_id = config$run_id,
    selected_pcs = selected_pcs,
    n_resolutions = expected_resolution_count,
    resolution_min = min(resolution$resolution),
    resolution_max = max(resolution$resolution),
    expected_iterations = expected_iterations,
    complete_checkpoints = complete_count,
    jaccard_checkpoints = jaccard_count,
    assignment_checkpoints = assignment_count,
    iteration_summary_rows = nrow(iteration),
    cluster_stability_rows = nrow(cluster),
    stage_completed_at = manifest$completed_at[[1L]],
    marker_mode = "none",
    pair_consensus = FALSE,
    sensitivity = FALSE,
    full_seurat_saved = FALSE,
    integrity_status = if (passed) "passed" else "failed",
    stringsAsFactors = FALSE
  )
}

iteration_all <- do.call(rbind, iteration_rows)
cluster_all <- do.call(rbind, cluster_rows)
resolution_all <- do.call(rbind, resolution_rows)
audit_all <- do.call(rbind, audit_rows)
row.names(iteration_all) <- NULL
row.names(cluster_all) <- NULL
row.names(resolution_all) <- NULL
row.names(audit_all) <- NULL

summary_dir <- file.path(output_root, "summary")
utils::write.csv(
  audit_all, file.path(summary_dir, "structural_pilot_run_status_audit.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  iteration_all, file.path(summary_dir, "structural_pilot_iteration_evidence.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  cluster_all, file.path(summary_dir, "structural_pilot_cluster_evidence.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  resolution_all, file.path(summary_dir, "structural_pilot_resolution_evidence.csv"),
  row.names = FALSE, na = ""
)

status_lines <- c(
  "| Species | Resolutions | PCs | Pilot checkpoints | Integrity |",
  "|---|---:|---:|---:|---|",
  vapply(seq_len(nrow(audit_all)), function(i) {
    sprintf(
      "| %s | %d | %d | %d/%d | %s |",
      audit_all$species[[i]], audit_all$n_resolutions[[i]], audit_all$selected_pcs[[i]],
      audit_all$complete_checkpoints[[i]], audit_all$expected_iterations[[i]],
      audit_all$integrity_status[[i]]
    )
  }, character(1))
)
coverage_lines <- c(
  "| Species | Resolution points | All clusters supported | Unsupported clusters, min/median/max | Unsupported cell fraction, min/max |",
  "|---|---:|---:|---:|---:|",
  vapply(species, function(sample) {
    x <- resolution_all[resolution_all$species == sample, , drop = FALSE]
    sprintf(
      "| %s | %d | %d | %d / %.1f / %d | %.3f / %.3f |",
      sample, nrow(x), sum(x$all_reference_clusters_resample_supported),
      min(x$resample_unsupported_cluster_count),
      stats::median(x$resample_unsupported_cluster_count),
      max(x$resample_unsupported_cluster_count),
      min(x$unsupported_full_cell_fraction), max(x$unsupported_full_cell_fraction)
    )
  }, character(1))
)

record <- c(
  "# Structural Pilot results",
  "",
  "## Current output integrity",
  "",
  status_lines,
  "",
  "PC counts come from each scan preprocessing audit. Candidate sets and expected iteration counts come from the current YAML. Checkpoints are counted only in the matching Pilot/count-mode/candidate-set directory. The audit CSV records the Pilot stage manifest completion time; batch start times, elapsed times and exit codes are not inferred.",
  "",
  "This summary covers nine species. Pcru has a separate summary and is not included here. Integrity status describes the checked files, not a historical batch execution status.",
  "",
  "## Resolution coverage",
  "",
  coverage_lines,
  "",
  "All clusters supported counts resolutions where every reference cluster has a finite median Jaccard at or above the configured stable_jaccard_threshold. This describes coverage and does not select a resolution.",
  "",
  "## Interpretation",
  "",
  "For each subsample, each full-data reference cluster independently matches the resampled cluster with the largest Jaccard: |C intersect D| / |C union D|. Only clusters meeting min_sampled_cells_per_cluster enter the stability statistics. Matching is not constrained to be one-to-one.",
  "",
  "stable_cluster_fraction is the fraction of eligible clusters meeting the Jaccard threshold; stable_cell_fraction weights these indicators by sampled reference-cluster cell counts. Median and lower-tail Jaccard summarize eligible reference clusters, while partition ARI/NMI summarize the sampled-cell partitions. unstable_cluster_count = n_sufficient_clusters - stable_cluster_count.",
  "",
  "Cluster-level evidence summarizes each full-data reference cluster across Pilot iterations. Finite median_jaccard values at or above the configured threshold are supported; values below it are unsupported; unavailable values are insufficient. These are resampling support labels, not the broader fragile-cluster classification that can also use pair consensus and sensitivity.",
  "",
  "Low recall can indicate splitting or loss; low precision can indicate merging; a low match margin indicates ambiguous matches. Pilot provides candidate-range evidence rather than a permanent property of a resolution or biological cluster.",
  "",
  "## Output tables",
  "",
  "- summary/structural_pilot_resolution_evidence.csv: global metrics, unstable cluster counts and unsupported full-cell fractions per resolution.",
  "- summary/structural_pilot_cluster_evidence.csv: cluster-level support and full-data cluster sizes; input to script 08.",
  "- summary/structural_pilot_iteration_evidence.csv: per-subsample stable/unstable counts and ARI/NMI.",
  "- summary/structural_pilot_run_status_audit.csv: current output integrity, selected PCs and stage completion records.",
  "",
  "Manual candidate selection should consider contiguous ranges with fewer unsupported clusters and smaller unsupported full-cell fractions together with lower-tail Jaccard, stable-cell fraction and ARI. Do not rank solely by a single mean or a metric that inherently favours low resolutions."
)
writeLines(record, file.path(output_root, "docs", "03_structural_pilot_run_record.md"))

review_code_files <- c(
  script_path,
  file.path(workflow_root, "R", "common.R"),
  file.path(workflow_root, "R", "metrics.R"),
  file.path(workflow_root, "scripts", "02_run_resampling.R")
)
review_code_manifest <- data.frame(
  relative_path = substring(review_code_files, nchar(workflow_root) + 2L),
  md5 = unname(tools::md5sum(review_code_files)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  review_code_manifest,
  file.path(output_root, "manifests", "structural_pilot_review_code_md5.csv"),
  row.names = FALSE
)

if (any(audit_all$integrity_status != "passed")) {
  stop("One or more structural Pilot outputs failed integrity checks.")
}
message("Structural Pilot audit passed for all nine samples.")
