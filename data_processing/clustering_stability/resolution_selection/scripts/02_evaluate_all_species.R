#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1L) args[[1L]] else paste0(".", ""
)
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(
  base_dir, "resolution_selection"
)
stable_threshold <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 0.6
min_sampled_cells <- if (length(args) >= 4L) as.integer(args[[4L]]) else 10L
species <- c(
  "Adig", "Ahya", "Amur", "Fusp", "Gfas",
  "Hcoe", "Pdam", "Plut", "Xesp", "Pcru"
)

required_packages <- c("yaml", "data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

script_dir_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_dir_arg[[1L]]))
script_dir <- dirname(script_path)
single_species_script <- file.path(
  script_dir, "01_evaluate_species_cluster_jaccard.R"
)
if (!file.exists(single_species_script)) {
  stop("Missing single-species script: ", single_species_script)
}

summary_dir <- file.path(output_dir, "03_cross_species_summary")
summary_table_dir <- file.path(summary_dir, "tables")
summary_doc_dir <- file.path(summary_dir, "docs")
summary_log_dir <- file.path(summary_dir, "logs")
reference_dir <- file.path(summary_dir, "reference_partitions")
for (path in c(
  summary_table_dir, summary_doc_dir, summary_log_dir, reference_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

find_run_dir <- function(sample_name) {
  candidates <- list.dirs(
    file.path(base_dir, sample_name, "runs"),
    recursive = FALSE,
    full.names = TRUE
  )
  candidates <- candidates[
    file.exists(file.path(
      candidates, "03_resampling", "tables",
      "cluster_stability__formal__raw_integer_full_pipeline.csv"
    ))
  ]
  if (length(candidates) != 1L) {
    stop(
      "Expected exactly one completed run for ", sample_name,
      "; found ", length(candidates)
    )
  }
  candidates[[1L]]
}

value_from_audit <- function(audit, item, default = NA_character_) {
  hit <- audit$value[audit$item == item]
  if (!length(hit)) default else as.character(hit[[1L]])
}

resolution_column <- function(x) {
  value <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  value <- sub("\\.?0+$", "", value)
  paste0("resolution_", gsub("\\.", "_", value))
}

read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

package_value <- function(package_table, package) {
  hit <- package_table$version[package_table$package == package]
  if (!length(hit)) NA_character_ else as.character(hit[[1L]])
}

combined_rows <- vector("list", length(species))
package_rows <- vector("list", length(species))
run_dirs <- setNames(character(length(species)), species)

for (i in seq_along(species)) {
  sample_name <- species[[i]]
  run_dir <- find_run_dir(sample_name)
  run_dirs[[sample_name]] <- run_dir
  jaccard_path <- file.path(
    run_dir, "03_resampling", "tables",
    "cluster_jaccard__formal__raw_integer_full_pipeline.csv"
  )
  stability_path <- file.path(
    run_dir, "03_resampling", "tables",
    "cluster_stability__formal__raw_integer_full_pipeline.csv"
  )
  assignment_path <- file.path(
    run_dir, "02_full_resolution_scan", "tables",
    "full_resolution_assignments.tsv.gz"
  )
  sample_output_dir <- file.path(output_dir, "species", sample_name)
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)
  sample_log <- file.path(summary_log_dir, paste0(sample_name, "_evaluation.log"))

  command_args <- c(
    single_species_script,
    jaccard_path,
    stability_path,
    sample_output_dir,
    stable_threshold,
    min_sampled_cells,
    sample_name,
    assignment_path
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = command_args,
    stdout = sample_log,
    stderr = sample_log
  )
  if (!identical(status, 0L)) {
    stop("Cluster-Jaccard evaluation failed for ", sample_name, ": ", sample_log)
  }

  selected_path <- file.path(
    sample_output_dir, "02_selected_resolution", "tables",
    paste0(sample_name, "_selected_resolution.csv")
  )
  selected <- read_csv(selected_path)
  selected_resolution <- selected$resolution[[1L]]

  config_path <- file.path(run_dir, "00_run_manifest", "config_frozen.yaml")
  audit_path <- file.path(
    run_dir, "01_input_audit", "tables",
    "canonical_preprocessing_audit.csv"
  )
  run_manifest_path <- file.path(
    run_dir, "00_run_manifest", "run_manifest.csv"
  )
  stage_manifest_path <- file.path(
    run_dir, "02_full_resolution_scan", "stage_manifest.csv"
  )
  package_path <- file.path(
    run_dir, "01_input_audit", "tables", "package_versions.csv"
  )
  assignment_path <- file.path(
    run_dir, "02_full_resolution_scan", "tables",
    "full_resolution_assignments.tsv.gz"
  )

  config <- yaml::read_yaml(config_path)
  audit <- read_csv(audit_path)
  run_manifest <- read_csv(run_manifest_path)
  stage_manifest <- read_csv(stage_manifest_path)
  package_table <- read_csv(package_path)
  assignment_column <- resolution_column(selected_resolution)

  reference <- data.table::fread(
    assignment_path,
    select = c("cell_id", assignment_column),
    showProgress = FALSE
  )
  data.table::setnames(reference, assignment_column, "reference_cluster")
  reference[, resolution := selected_resolution]
  data.table::setcolorder(
    reference, c("cell_id", "resolution", "reference_cluster")
  )
  verified_clusters <- data.table::uniqueN(reference$reference_cluster)
  if (verified_clusters != selected$n_clusters[[1L]]) {
    stop(
      sample_name, " selected cluster count mismatch: stability=",
      selected$n_clusters[[1L]], ", assignment=", verified_clusters
    )
  }

  reference_name <- paste0(
    sample_name, "_r", selected_resolution, "_reference_partition.tsv.gz"
  )
  reference_output <- file.path(reference_dir, reference_name)
  data.table::fwrite(
    reference, reference_output, sep = "\t", quote = FALSE, compress = "gzip"
  )

  stage_detail <- as.character(stage_manifest$detail[[1L]])
  selected_pcs <- as.integer(value_from_audit(audit, "selected_pcs"))
  package_rows[[i]] <- transform(package_table, species = sample_name)[
    , c("species", "package", "available", "version")
  ]
  combined_rows[[i]] <- data.frame(
    species = sample_name,
    run_id = config$run_id,
    selected_resolution = selected_resolution,
    n_clusters = selected$n_clusters[[1L]],
    stable_cluster_count = selected$stable_cluster_count[[1L]],
    unstable_cluster_count = selected$unstable_cluster_count[[1L]],
    stable_cluster_fraction = selected$stable_cluster_fraction[[1L]],
    stable_cluster_cell_fraction =
      selected$stable_cluster_cell_fraction[[1L]],
    minimum_cluster_median_jaccard =
      selected$minimum_cluster_median_jaccard[[1L]],
    q10_cluster_median_jaccard =
      selected$q10_cluster_median_jaccard[[1L]],
    unstable_cluster_mean_median_jaccard =
      selected$unstable_cluster_mean_median_jaccard[[1L]],
    stable_jaccard_threshold = stable_threshold,
    selection_rule = paste(
      "max stable cluster fraction;",
      "then max unstable-cluster mean median Jaccard;",
      "then max q10 cluster median Jaccard; then lower resolution"
    ),
    input_rds = run_manifest$input_rds[[1L]],
    input_rds_md5 = run_manifest$input_rds_md5[[1L]],
    count_mode = run_manifest$count_mode[[1L]],
    integerization = run_manifest$integerization[[1L]],
    count_assay = config$input$count_assay,
    count_layer = config$input$count_layer,
    sct_assay = config$input$sct_assay,
    rebuild_sct = isTRUE(config$canonical_preprocessing$rebuild_sct),
    sct_variable_features_n =
      as.integer(config$canonical_preprocessing$variable_features_n),
    sct_do_scale = TRUE,
    sct_do_center = TRUE,
    sct_return_only_var_genes = TRUE,
    pca_npcs_requested = as.integer(config$clustering$pca_npcs),
    pca_selection_method =
      as.character(config$canonical_preprocessing$pc_selection_method),
    selected_pcs_actual = selected_pcs,
    pca_reduction = config$input$reduction,
    pca_reduction_key = config$canonical_preprocessing$reduction_key,
    neighbor_dims = paste0("1:", selected_pcs),
    neighbor_k = as.integer(config$clustering$primary_k),
    neighbor_graph = "stability_scan_snn",
    clustering_algorithm = as.integer(config$clustering$algorithm),
    full_pipeline_seed = as.integer(config$compute$seed),
    findclusters_random_seed = as.integer(config$compute$seed),
    leiden_seed_mode = "shared_across_resolutions",
    resample_fraction = as.numeric(config$resampling$fraction),
    formal_subsample_iterations =
      as.integer(config$resampling$formal_iterations),
    formal_subsample_seed_rule = paste0(
      config$compute$seed, " + iteration (1..",
      config$resampling$formal_iterations, ")"
    ),
    reference_assignment_column = assignment_column,
    reference_assignment_source = assignment_path,
    reference_assignment_source_md5 =
      unname(tools::md5sum(assignment_path)),
    frozen_reference_partition = reference_output,
    frozen_reference_partition_md5 =
      unname(tools::md5sum(reference_output)),
    reference_cluster_count_verified = verified_clusters,
    scan_stage_detail = stage_detail,
    Seurat_version = package_value(package_table, "Seurat"),
    SeuratObject_version = package_value(package_table, "SeuratObject"),
    Matrix_version = package_value(package_table, "Matrix"),
    config_frozen = config_path,
    config_md5 = unname(tools::md5sum(config_path)),
    session_info = file.path(
      run_dir, "01_input_audit", "logs", "session_info.txt"
    ),
    stringsAsFactors = FALSE
  )
}

combined <- do.call(rbind, combined_rows)
package_versions <- do.call(rbind, package_rows)
write_csv(
  combined,
  file.path(
    summary_table_dir,
    "all_species_selected_resolution_and_reclustering_parameters.csv"
  )
)
write_csv(
  package_versions,
  file.path(
    summary_table_dir,
    "all_species_archived_package_versions.csv"
  )
)

compact <- combined[, c(
  "species", "selected_resolution", "n_clusters",
  "stable_cluster_count", "unstable_cluster_count",
  "stable_cluster_fraction", "stable_cluster_cell_fraction",
  "unstable_cluster_mean_median_jaccard",
  "q10_cluster_median_jaccard", "selected_pcs_actual", "neighbor_k",
  "clustering_algorithm", "findclusters_random_seed",
  "frozen_reference_partition"
)]
write_csv(
  compact,
  file.path(
    summary_table_dir,
    "all_species_selected_resolution_compact_summary.csv"
  )
)

lines <- c(
  "# Cross-species cluster-level Jaccard selection",
  "",
  paste0(
    "All species were evaluated with the same median-Jaccard threshold of ",
    stable_threshold, " using their completed 100-iteration formal runs."
  ),
  "",
  "## Selected partitions",
  "",
  "| Species | Resolution | Clusters | Stable clusters | Stable fraction | Stable-cell fraction | Unstable mean median Jaccard | Cluster-median q10 | PCs | k | Algorithm | Seed |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(compact))) {
  lines <- c(
    lines,
    paste0(
      "| ", compact$species[[i]],
      " | ", compact$selected_resolution[[i]],
      " | ", compact$n_clusters[[i]],
      " | ", compact$stable_cluster_count[[i]], "/",
      compact$n_clusters[[i]],
      " | ", sprintf("%.3f", compact$stable_cluster_fraction[[i]]),
      " | ", sprintf("%.3f", compact$stable_cluster_cell_fraction[[i]]),
      " | ", ifelse(
        is.na(compact$unstable_cluster_mean_median_jaccard[[i]]),
        "NA", sprintf("%.3f", compact$unstable_cluster_mean_median_jaccard[[i]])
      ),
      " | ", sprintf("%.3f", compact$q10_cluster_median_jaccard[[i]]),
      " | ", compact$selected_pcs_actual[[i]],
      " | ", compact$neighbor_k[[i]],
      " | ", compact$clustering_algorithm[[i]],
      " | ", compact$findclusters_random_seed[[i]],
      " |"
    )
  )
}
lines <- c(
  lines, "",
  "## Reproduction contract",
  "",
  "Reproducing the selected partition requires the exact input RDS and cell order,",
  "the frozen SCT/PCA workflow, the sample-specific selected PC count, k=20,",
  paste0(
    "Leiden algorithm 4, random seed ",
    paste(sort(unique(compact$findclusters_random_seed)), collapse = ", "),
    ", and the archived package environment."
  ),
  "The full parameter and checksum contract is stored in",
  "`all_species_selected_resolution_and_reclustering_parameters.csv`.",
  "",
  "Cluster numbers are labels, not identities. A rerun must be validated against",
  "the frozen cell-to-cluster partition using label-invariant ARI/NMI and a",
  "cluster-overlap mapping. Exact numeric labels are not required for equivalence.",
  "",
  "Even with matching explicit parameters, changes in package versions, BLAS,",
  "threading, nearest-neighbor defaults, or input ordering can alter graph ties.",
  "Therefore the table documents a reproducibility target rather than assuming",
  "that resolution and seed alone guarantee identical clusters."
)
writeLines(
  lines,
  file.path(summary_doc_dir, "all_species_selection_and_reproduction.md")
)

manifest <- c(
  paste0("run_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("base_dir: ", normalizePath(base_dir)),
  paste0("output_dir: ", normalizePath(output_dir)),
  paste0("species: ", paste(species, collapse = ",")),
  paste0("stable_jaccard_threshold: ", stable_threshold),
  paste0("min_sampled_cells_per_cluster: ", min_sampled_cells)
)
writeLines(manifest, file.path(summary_log_dir, "batch_run_manifest.txt"))
writeLines(
  capture.output(sessionInfo()),
  file.path(summary_log_dir, "batch_sessionInfo.txt")
)

message("Completed species: ", paste(species, collapse = ", "))
message("Combined outputs: ", normalizePath(summary_dir))
