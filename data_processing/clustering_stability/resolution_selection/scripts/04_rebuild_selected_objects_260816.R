#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(yaml)
  library(data.table)
  library(mclust)
  library(aricode)
})

args <- commandArgs(trailingOnly = TRUE)
analysis_root <- paste0(".", ""
)
parameter_table <- if (length(args) >= 1L) args[[1L]] else file.path(
  analysis_root, "resolution_selection", "03_cross_species_summary", "tables",
  "all_species_selected_resolution_and_reclustering_parameters.csv"
)
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(
  analysis_root, "selected_objects"
)
workflow_root <- if (length(args) >= 3L) args[[3L]] else file.path(
  analysis_root, "workflow"
)
species_filter <- if (length(args) >= 4L && nzchar(args[[4L]])) {
  trimws(strsplit(args[[4L]], ",", fixed = TRUE)[[1L]])
} else {
  character()
}

common_path <- file.path(workflow_root, "R", "common.R")
if (!file.exists(common_path)) stop("Missing workflow helper: ", common_path)
source(common_path)

if (!file.exists(parameter_table)) stop("Missing parameter table: ", parameter_table)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
verification_dir <- file.path(output_dir, "verification")
mapping_dir <- file.path(verification_dir, "cluster_mappings")
log_dir <- file.path(output_dir, "logs")
manifest_dir <- file.path(output_dir, "manifests")
dir.create(verification_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mapping_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

align_rows <- function(rows) {
  if (!length(rows)) return(data.frame())
  all_columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned <- lapply(rows, function(x) {
    for (column in setdiff(all_columns, names(x))) x[[column]] <- NA
    x[, all_columns, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  row.names(out) <- NULL
  out
}

safe_ari <- function(reference, test) {
  if (length(unique(reference)) < 2L || length(unique(test)) < 2L) {
    return(as.numeric(identical(as.character(reference), as.character(test))))
  }
  as.numeric(mclust::adjustedRandIndex(reference, test))
}

safe_nmi <- function(reference, test) {
  if (length(unique(reference)) < 2L || length(unique(test)) < 2L) {
    return(as.numeric(identical(as.character(reference), as.character(test))))
  }
  as.numeric(aricode::NMI(reference, test))
}

compare_partitions <- function(reference, test_labels) {
  reference <- data.table::as.data.table(reference)
  required <- c("cell_id", "reference_cluster")
  if (!all(required %in% names(reference))) {
    stop("Reference partition requires columns: ", paste(required, collapse = ", "))
  }
  reference[, cell_id := as.character(cell_id)]
  reference[, reference_cluster := as.character(reference_cluster)]
  if (anyDuplicated(reference$cell_id)) stop("Reference partition has duplicated cell IDs")

  test <- data.table::data.table(
    cell_id = names(test_labels),
    rebuilt_cluster = as.character(test_labels)
  )
  common_cells <- intersect(reference$cell_id, test$cell_id)
  comparison <- merge(
    reference[cell_id %in% common_cells, .(cell_id, reference_cluster)],
    test[cell_id %in% common_cells],
    by = "cell_id", all = FALSE, sort = FALSE
  )
  edges <- comparison[, .(intersection_n = .N),
                      by = .(reference_cluster, rebuilt_cluster)]
  reference_sizes <- comparison[, .(reference_cluster_n = .N),
                                by = reference_cluster]
  rebuilt_sizes <- comparison[, .(rebuilt_cluster_n = .N),
                              by = rebuilt_cluster]
  edges <- merge(edges, reference_sizes, by = "reference_cluster")
  edges <- merge(edges, rebuilt_sizes, by = "rebuilt_cluster")
  edges[, reference_recall := intersection_n / reference_cluster_n]
  edges[, rebuilt_precision := intersection_n / rebuilt_cluster_n]
  edges[, jaccard := intersection_n /
          (reference_cluster_n + rebuilt_cluster_n - intersection_n)]
  data.table::setorder(edges, reference_cluster, -intersection_n)

  one_to_one <- nrow(edges) > 0L &&
    all(edges[, .N, by = reference_cluster]$N == 1L) &&
    all(edges[, .N, by = rebuilt_cluster]$N == 1L)
  cell_sets_equal <- setequal(reference$cell_id, test$cell_id)
  ari <- safe_ari(comparison$reference_cluster, comparison$rebuilt_cluster)
  nmi <- safe_nmi(comparison$reference_cluster, comparison$rebuilt_cluster)
  reference_n_clusters <- data.table::uniqueN(reference$reference_cluster)
  rebuilt_n_clusters <- data.table::uniqueN(test$rebuilt_cluster)
  exact <- cell_sets_equal && one_to_one &&
    reference_n_clusters == rebuilt_n_clusters &&
    isTRUE(all.equal(ari, 1, tolerance = 1e-12)) &&
    isTRUE(all.equal(nmi, 1, tolerance = 1e-12)) &&
    all(abs(edges$jaccard - 1) < 1e-12)

  list(
    edges = edges,
    metrics = list(
      reference_n_cells = data.table::uniqueN(reference$cell_id),
      rebuilt_n_cells = data.table::uniqueN(test$cell_id),
      common_n_cells = length(common_cells),
      reference_only_n_cells = sum(!reference$cell_id %in% test$cell_id),
      rebuilt_only_n_cells = sum(!test$cell_id %in% reference$cell_id),
      reference_n_clusters = reference_n_clusters,
      rebuilt_n_clusters = rebuilt_n_clusters,
      cell_sets_equal = cell_sets_equal,
      one_to_one_mapping = one_to_one,
      ARI = ari,
      NMI = nmi,
      exact_partition_match = exact
    )
  )
}

parameters <- utils::read.csv(
  parameter_table, stringsAsFactors = FALSE, check.names = FALSE
)
required_columns <- c(
  "species", "selected_resolution", "n_clusters", "input_rds",
  "input_rds_md5", "count_assay", "count_layer", "sct_assay",
  "sct_variable_features_n", "pca_npcs_requested", "selected_pcs_actual",
  "pca_selection_method", "pca_reduction", "pca_reduction_key",
  "neighbor_k", "neighbor_graph",
  "clustering_algorithm", "findclusters_random_seed",
  "frozen_reference_partition", "frozen_reference_partition_md5",
  "config_frozen", "config_md5", "Seurat_version",
  "SeuratObject_version", "Matrix_version"
)
missing_columns <- setdiff(required_columns, names(parameters))
if (length(missing_columns)) {
  stop("Missing parameter columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(parameters$species)) stop("Parameter table has duplicated species")
if (length(species_filter)) {
  unknown <- setdiff(species_filter, parameters$species)
  if (length(unknown)) stop("Unknown species filter: ", paste(unknown, collapse = ", "))
  parameters <- parameters[match(species_filter, parameters$species), , drop = FALSE]
}

# Fail the whole batch before expensive work if any frozen input is no longer valid.
package_versions <- c(
  Seurat = as.character(utils::packageVersion("Seurat")),
  SeuratObject = as.character(utils::packageVersion("SeuratObject")),
  Matrix = as.character(utils::packageVersion("Matrix"))
)
preflight_rows <- vector("list", nrow(parameters))
for (i in seq_len(nrow(parameters))) {
  row <- parameters[i, , drop = FALSE]
  input_exists <- file.exists(row$input_rds[[1L]])
  reference_exists <- file.exists(row$frozen_reference_partition[[1L]])
  config_exists <- file.exists(row$config_frozen[[1L]])
  current_input_md5 <- if (input_exists) unname(tools::md5sum(row$input_rds[[1L]])) else NA_character_
  current_reference_md5 <- if (reference_exists) unname(tools::md5sum(row$frozen_reference_partition[[1L]])) else NA_character_
  current_config_md5 <- if (config_exists) unname(tools::md5sum(row$config_frozen[[1L]])) else NA_character_
  preflight_rows[[i]] <- data.frame(
    species = row$species[[1L]],
    input_exists = input_exists,
    input_md5_expected = row$input_rds_md5[[1L]],
    input_md5_current = current_input_md5,
    input_md5_match = identical(current_input_md5, row$input_rds_md5[[1L]]),
    reference_exists = reference_exists,
    reference_md5_match = identical(
      current_reference_md5, row$frozen_reference_partition_md5[[1L]]
    ),
    config_exists = config_exists,
    config_md5_match = identical(current_config_md5, row$config_md5[[1L]]),
    R_version_current = as.character(getRversion()),
    R_version_match = identical(as.character(getRversion()), "4.3.3"),
    Seurat_version_match = identical(
      package_versions[["Seurat"]], row$Seurat_version[[1L]]
    ),
    SeuratObject_version_match = identical(
      package_versions[["SeuratObject"]], row$SeuratObject_version[[1L]]
    ),
    Matrix_version_match = identical(
      package_versions[["Matrix"]], row$Matrix_version[[1L]]
    ),
    stringsAsFactors = FALSE
  )
}
preflight <- do.call(rbind, preflight_rows)
preflight$preflight_pass <- apply(
  preflight[c(
    "input_exists", "input_md5_match", "reference_exists",
    "reference_md5_match", "config_exists", "config_md5_match",
    "R_version_match", "Seurat_version_match",
    "SeuratObject_version_match", "Matrix_version_match"
  )], 1L, function(x) all(as.logical(x))
)
preflight_path <- file.path(verification_dir, "batch_preflight.csv")
write_csv(preflight, preflight_path)
if (!all(preflight$preflight_pass)) {
  stop("Preflight failed; no RDS was rebuilt. See ", preflight_path)
}

summary_path <- file.path(
  verification_dir, "all_species_reclustering_consistency.csv"
)
results <- vector("list", nrow(parameters))
names(results) <- parameters$species

for (i in seq_len(nrow(parameters))) {
  row <- parameters[i, , drop = FALSE]
  species <- row$species[[1L]]
  log_path <- file.path(log_dir, paste0(species, "_rebuild.log"))
  log_connection <- file(log_path, open = "wt")
  sink(log_connection, type = "output")
  sink(log_connection, type = "message")

  result <- tryCatch({
    started <- Sys.time()
    cat("Species:", species, "\n")
    cat("Started:", format(started), "\n")
    cat("Input:", row$input_rds[[1L]], "\n")
    cat("Selected resolution:", row$selected_resolution[[1L]], "\n")
    cat("Frozen PCs:", row$selected_pcs_actual[[1L]], "\n")

    config <- yaml::read_yaml(row$config_frozen[[1L]])
    config$input$seurat_rds <- row$input_rds[[1L]]
    checks <- c(
      identical(config$input$count_assay, row$count_assay[[1L]]),
      identical(config$input$count_layer, row$count_layer[[1L]]),
      identical(config$input$sct_assay, row$sct_assay[[1L]]),
      identical(config$input$reduction, row$pca_reduction[[1L]]),
      identical(
        config$canonical_preprocessing$reduction_key,
        row$pca_reduction_key[[1L]]
      ),
      identical(
        config$canonical_preprocessing$pc_selection_method,
        row$pca_selection_method[[1L]]
      ),
      as.integer(config$canonical_preprocessing$variable_features_n) ==
        as.integer(row$sct_variable_features_n[[1L]]),
      as.integer(config$clustering$pca_npcs) ==
        as.integer(row$pca_npcs_requested[[1L]]),
      as.integer(config$clustering$primary_k) ==
        as.integer(row$neighbor_k[[1L]]),
      as.integer(config$clustering$algorithm) ==
        as.integer(row$clustering_algorithm[[1L]]),
      identical(row$neighbor_graph[[1L]], "stability_scan_snn"),
      as.integer(config$compute$seed) ==
        as.integer(row$findclusters_random_seed[[1L]])
    )
    if (!all(checks)) stop(species, " frozen config differs from parameter table")

    object <- readRDS(row$input_rds[[1L]])
    if (!inherits(object, "Seurat")) stop(species, " input is not a Seurat object")
    counts <- get_assay_matrix(
      object, config$input$count_assay, config$input$count_layer
    )
    if (!identical(colnames(counts), colnames(object))) {
      stop(species, " count matrix cell order differs from object")
    }

    rebuilt <- rebuild_canonical_sct_pca(object, config)
    object <- rebuilt$object
    auto_selected_pcs <- as.integer(rebuilt$record$selected_pcs)
    frozen_pcs <- as.integer(row$selected_pcs_actual[[1L]])
    if (auto_selected_pcs != frozen_pcs) {
      stop(
        species, " auto-selected PCs changed: current=", auto_selected_pcs,
        ", frozen=", frozen_pcs
      )
    }
    cat("Auto-selected PCs reproduced:", auto_selected_pcs, "\n")

    set.seed(as.integer(row$findclusters_random_seed[[1L]]))
    Seurat::DefaultAssay(object) <- config$input$sct_assay
    object <- Seurat::FindNeighbors(
      object,
      reduction = config$input$reduction,
      dims = seq_len(frozen_pcs),
      k.param = as.integer(row$neighbor_k[[1L]]),
      graph.name = c("stability_scan_nn", "stability_scan_snn"),
      return.neighbor = FALSE,
      verbose = FALSE
    )
    cluster_column <- "stability_selected_cluster"
    object <- Seurat::FindClusters(
      object,
      graph.name = "stability_scan_snn",
      resolution = as.numeric(row$selected_resolution[[1L]]),
      algorithm = as.integer(row$clustering_algorithm[[1L]]),
      cluster.name = cluster_column,
      random.seed = as.integer(row$findclusters_random_seed[[1L]]),
      verbose = FALSE
    )
    rebuilt_labels <- as.character(object[[cluster_column, drop = TRUE]])
    names(rebuilt_labels) <- colnames(object)

    reference <- data.table::fread(
      row$frozen_reference_partition[[1L]], showProgress = FALSE
    )
    comparison <- compare_partitions(reference, rebuilt_labels)
    metrics <- comparison$metrics
    expected_clusters <- as.integer(row$n_clusters[[1L]])
    if (metrics$reference_n_clusters != expected_clusters ||
        metrics$rebuilt_n_clusters != expected_clusters) {
      metrics$exact_partition_match <- FALSE
    }
    mapping_path <- file.path(
      mapping_dir, paste0(species, "_cluster_mapping.csv")
    )
    write_csv(as.data.frame(comparison$edges), mapping_path)
    cat("Rebuilt clusters:", metrics$rebuilt_n_clusters, "\n")
    cat("Cell sets equal:", metrics$cell_sets_equal, "\n")
    cat("One-to-one mapping:", metrics$one_to_one_mapping, "\n")
    cat("ARI:", sprintf("%.15f", metrics$ARI), "\n")
    cat("NMI:", sprintf("%.15f", metrics$NMI), "\n")
    cat("Exact partition match:", metrics$exact_partition_match, "\n")

    if (!isTRUE(metrics$exact_partition_match)) {
      stop(species, " rebuilt partition does not exactly match frozen reference")
    }

    frozen_labels <- setNames(
      as.character(reference$reference_cluster), as.character(reference$cell_id)
    )
    object$stability_frozen_reference_cluster <- unname(
      frozen_labels[colnames(object)]
    )
    object$stability_selected_resolution <- as.numeric(
      row$selected_resolution[[1L]]
    )
    object$stability_selected_pcs <- frozen_pcs
    object$stability_neighbor_k <- as.integer(row$neighbor_k[[1L]])
    object$stability_leiden_algorithm <- as.integer(
      row$clustering_algorithm[[1L]]
    )
    object$stability_leiden_seed <- as.integer(
      row$findclusters_random_seed[[1L]]
    )
    Seurat::Idents(object) <- cluster_column

    object@misc$selected_resolution_rebuild <- list(
      species = species,
      rebuilt_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      input_rds = row$input_rds[[1L]],
      input_rds_md5 = row$input_rds_md5[[1L]],
      parameter_table = parameter_table,
      parameter_table_md5 = unname(tools::md5sum(parameter_table)),
      workflow_common = common_path,
      workflow_common_md5 = unname(tools::md5sum(common_path)),
      selected_resolution = as.numeric(row$selected_resolution[[1L]]),
      pca_npcs_requested = as.integer(row$pca_npcs_requested[[1L]]),
      auto_selected_pcs = auto_selected_pcs,
      neighbor_k = as.integer(row$neighbor_k[[1L]]),
      graph_name = "stability_scan_snn",
      algorithm = as.integer(row$clustering_algorithm[[1L]]),
      seed = as.integer(row$findclusters_random_seed[[1L]]),
      frozen_reference_partition = row$frozen_reference_partition[[1L]],
      frozen_reference_partition_md5 =
        row$frozen_reference_partition_md5[[1L]],
      comparison = metrics,
      package_versions = package_versions
    )

    resolution_text <- format(
      as.numeric(row$selected_resolution[[1L]]), scientific = FALSE, trim = TRUE
    )
    output_path <- file.path(
      output_dir,
      paste0(species, "_selected_r", resolution_text, "_reclustered.rds")
    )
    temporary_path <- paste0(output_path, ".tmp")
    if (file.exists(temporary_path)) unlink(temporary_path)
    saveRDS(object, temporary_path, compress = TRUE)
    if (!file.rename(temporary_path, output_path)) {
      unlink(temporary_path)
      stop("Failed to atomically move output to ", output_path)
    }
    output_md5 <- unname(tools::md5sum(output_path))
    completed <- Sys.time()
    cat("Output:", output_path, "\n")
    cat("Output MD5:", output_md5, "\n")
    cat("Completed:", format(completed), "\n")

    data.frame(
      species = species,
      selected_resolution = as.numeric(row$selected_resolution[[1L]]),
      selected_pcs = frozen_pcs,
      neighbor_k = as.integer(row$neighbor_k[[1L]]),
      clustering_algorithm = as.integer(row$clustering_algorithm[[1L]]),
      random_seed = as.integer(row$findclusters_random_seed[[1L]]),
      reference_n_cells = metrics$reference_n_cells,
      rebuilt_n_cells = metrics$rebuilt_n_cells,
      reference_n_clusters = metrics$reference_n_clusters,
      rebuilt_n_clusters = metrics$rebuilt_n_clusters,
      cell_sets_equal = metrics$cell_sets_equal,
      one_to_one_mapping = metrics$one_to_one_mapping,
      ARI = metrics$ARI,
      NMI = metrics$NMI,
      exact_partition_match = metrics$exact_partition_match,
      verification_status = "exact_match_saved",
      elapsed_minutes = as.numeric(difftime(completed, started, units = "mins")),
      output_rds = output_path,
      output_rds_md5 = output_md5,
      output_size_bytes = file.info(output_path)$size,
      mapping_table = mapping_path,
      log_file = log_path,
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    cat("ERROR:", conditionMessage(error), "\n")
    data.frame(
      species = species,
      selected_resolution = as.numeric(row$selected_resolution[[1L]]),
      verification_status = "failed_not_saved",
      error_message = conditionMessage(error),
      output_rds = NA_character_,
      log_file = log_path,
      stringsAsFactors = FALSE
    )
  })

  sink(type = "message")
  sink(type = "output")
  close(log_connection)
  results[[i]] <- result
  write_csv(align_rows(results[seq_len(i)]), summary_path)
  message(
    "[", i, "/", nrow(parameters), "] ", species, ": ",
    result$verification_status[[1L]]
  )
  rm(result)
  gc()
}

summary <- align_rows(results)
write_csv(summary, summary_path)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "rebuild_sessionInfo.txt"))

for (column in c("selected_pcs", "rebuilt_n_clusters", "ARI", "NMI")) {
  if (!column %in% names(summary)) summary[[column]] <- NA
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
manifest <- data.frame(
  file = c(parameter_table, common_path, script_path),
  md5 = unname(tools::md5sum(c(parameter_table, common_path, script_path))),
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(manifest_dir, "rebuild_code_and_parameter_md5.csv"))

report <- c(
  "# 260816 final RDS rebuild",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Canonical processing: RNA/counts -> SCTransform (3,000 variable features) -> ",
  "PCA (100 requested PCs) -> first-derivative PC selection -> FindNeighbors ",
  "(k=20) -> Leiden algorithm 4 (seed 260810).",
  "",
  "An RDS is saved only when its input/config/reference checksums match the ",
  "frozen records and the rebuilt partition has identical cells, cluster count, ",
  "one-to-one cluster mapping, ARI=1, NMI=1, and mapping Jaccard=1.",
  "",
  "| Species | Resolution | PCs | Clusters | ARI | NMI | Status |",
  "|---|---:|---:|---:|---:|---:|---|"
)
for (i in seq_len(nrow(summary))) {
  report <- c(report, paste0(
    "| ", summary$species[[i]], " | ", summary$selected_resolution[[i]],
    " | ", ifelse(is.na(summary$selected_pcs[[i]]), "", summary$selected_pcs[[i]]),
    " | ", ifelse(is.na(summary$rebuilt_n_clusters[[i]]), "", summary$rebuilt_n_clusters[[i]]),
    " | ", ifelse(is.na(summary$ARI[[i]]), "", sprintf("%.6f", summary$ARI[[i]])),
    " | ", ifelse(is.na(summary$NMI[[i]]), "", sprintf("%.6f", summary$NMI[[i]])),
    " | ", summary$verification_status[[i]], " |"
  ))
}
writeLines(report, file.path(output_dir, "README.md"))

if (any(summary$verification_status != "exact_match_saved")) {
  stop("One or more species failed strict verification; see ", summary_path)
}
message("All strictly verified RDS outputs written to: ", normalizePath(output_dir))
