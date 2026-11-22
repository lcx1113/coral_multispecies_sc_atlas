#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(mclust)
  library(aricode)
})

analysis_root <- paste0(".", ""
)
parameter_path <- file.path(
  analysis_root, "resolution_selection", "03_cross_species_summary", "tables",
  "all_species_selected_resolution_and_reclustering_parameters.csv"
)
output_dir <- file.path(analysis_root, "selected_objects")
rebuild_summary_path <- file.path(
  output_dir, "verification", "all_species_reclustering_consistency.csv"
)
validation_path <- file.path(
  output_dir, "verification", "post_save_rds_validation.csv"
)

parameters <- utils::read.csv(
  parameter_path, stringsAsFactors = FALSE, check.names = FALSE
)
rebuild_summary <- utils::read.csv(
  rebuild_summary_path, stringsAsFactors = FALSE, check.names = FALSE
)
rows <- vector("list", nrow(parameters))

for (i in seq_len(nrow(parameters))) {
  row <- parameters[i, , drop = FALSE]
  species <- row$species[[1L]]
  saved <- rebuild_summary[rebuild_summary$species == species, , drop = FALSE]
  if (nrow(saved) != 1L) stop("Missing rebuild summary row for ", species)
  path <- saved$output_rds[[1L]]
  message("[", i, "/", nrow(parameters), "] validating ", species)

  result <- tryCatch({
    current_md5 <- unname(tools::md5sum(path))
    object <- readRDS(path)
    reference <- data.table::fread(
      row$frozen_reference_partition[[1L]], showProgress = FALSE
    )
    reference_labels <- setNames(
      as.character(reference$reference_cluster), as.character(reference$cell_id)
    )
    cells <- colnames(object)
    required_metadata <- c(
      "stability_selected_cluster", "stability_frozen_reference_cluster",
      "stability_selected_resolution", "stability_selected_pcs",
      "stability_neighbor_k", "stability_leiden_algorithm",
      "stability_leiden_seed"
    )
    metadata_complete <- all(required_metadata %in% colnames(object[[]]))
    if (!metadata_complete) stop("Required stability metadata is incomplete")

    selected <- as.character(object$stability_selected_cluster)
    frozen <- as.character(object$stability_frozen_reference_cluster)
    expected <- unname(reference_labels[cells])
    names(selected) <- cells
    names(frozen) <- cells
    contingency <- table(expected, selected)
    nonzero_edges <- which(contingency > 0, arr.ind = TRUE)
    one_to_one <- nrow(nonzero_edges) == nrow(contingency) &&
      nrow(nonzero_edges) == ncol(contingency) &&
      all(rowSums(contingency > 0) == 1L) &&
      all(colSums(contingency > 0) == 1L)
    ari <- as.numeric(mclust::adjustedRandIndex(expected, selected))
    nmi <- as.numeric(aricode::NMI(expected, selected))

    assays_present <- all(c(row$count_assay[[1L]], row$sct_assay[[1L]]) %in%
                            names(object@assays))
    reduction_present <- row$pca_reduction[[1L]] %in% names(object@reductions)
    pca_dimensions <- if (reduction_present) {
      ncol(Seurat::Embeddings(object, reduction = row$pca_reduction[[1L]]))
    } else {
      NA_integer_
    }
    reduction_assay <- if (reduction_present) {
      object[[row$pca_reduction[[1L]]]]@assay.used
    } else {
      NA_character_
    }
    graphs_present <- all(c("stability_scan_nn", "stability_scan_snn") %in%
                            names(object@graphs))
    identity_matches <- identical(
      unname(as.character(Seurat::Idents(object))), unname(selected)
    )
    misc <- object@misc$selected_resolution_rebuild
    misc_complete <- !is.null(misc) &&
      isTRUE(misc$comparison$exact_partition_match)

    checks <- c(
      file.exists(path),
      identical(current_md5, saved$output_rds_md5[[1L]]),
      inherits(object, "Seurat"),
      setequal(cells, reference$cell_id),
      length(cells) == as.integer(saved$rebuilt_n_cells[[1L]]),
      identical(unname(frozen), expected),
      one_to_one,
      isTRUE(all.equal(ari, 1, tolerance = 1e-12)),
      isTRUE(all.equal(nmi, 1, tolerance = 1e-12)),
      length(unique(selected)) == as.integer(row$n_clusters[[1L]]),
      assays_present,
      reduction_present,
      pca_dimensions == as.integer(row$pca_npcs_requested[[1L]]),
      identical(reduction_assay, row$sct_assay[[1L]]),
      graphs_present,
      identity_matches,
      all(object$stability_selected_resolution ==
            as.numeric(row$selected_resolution[[1L]])),
      all(object$stability_selected_pcs ==
            as.integer(row$selected_pcs_actual[[1L]])),
      all(object$stability_neighbor_k == as.integer(row$neighbor_k[[1L]])),
      all(object$stability_leiden_algorithm ==
            as.integer(row$clustering_algorithm[[1L]])),
      all(object$stability_leiden_seed ==
            as.integer(row$findclusters_random_seed[[1L]])),
      misc_complete
    )

    data.frame(
      species = species,
      selected_resolution = as.numeric(row$selected_resolution[[1L]]),
      n_cells = length(cells),
      n_clusters = length(unique(selected)),
      selected_pcs = as.integer(row$selected_pcs_actual[[1L]]),
      pca_dimensions_saved = pca_dimensions,
      output_rds = path,
      output_rds_md5 = current_md5,
      output_md5_match = identical(current_md5, saved$output_rds_md5[[1L]]),
      cell_sets_equal = setequal(cells, reference$cell_id),
      frozen_metadata_exact = identical(unname(frozen), expected),
      one_to_one_mapping = one_to_one,
      ARI = ari,
      NMI = nmi,
      assays_present = assays_present,
      reduction_present = reduction_present,
      reduction_assay = reduction_assay,
      graphs_present = graphs_present,
      identity_matches_selected_cluster = identity_matches,
      misc_verification_complete = misc_complete,
      post_save_validation_pass = all(checks),
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    data.frame(
      species = species,
      selected_resolution = as.numeric(row$selected_resolution[[1L]]),
      output_rds = path,
      error_message = conditionMessage(error),
      post_save_validation_pass = FALSE,
      stringsAsFactors = FALSE
    )
  })
  rows[[i]] <- result

  columns <- unique(unlist(lapply(rows[seq_len(i)], names), use.names = FALSE))
  aligned <- lapply(rows[seq_len(i)], function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  utils::write.csv(do.call(rbind, aligned), validation_path,
                   row.names = FALSE, na = "")
  rm(result)
  gc()
}

validation <- utils::read.csv(
  validation_path, stringsAsFactors = FALSE, check.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "logs", "post_save_validation_sessionInfo.txt")
)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
utils::write.csv(
  data.frame(
    file = script_path,
    md5 = unname(tools::md5sum(script_path)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "manifests", "post_save_validation_code_md5.csv"),
  row.names = FALSE
)

report <- c(
  "# Post-save RDS validation",
  "",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Each saved object was read back from disk and checked for MD5 integrity,",
  "required assays/reduction/graphs, frozen parameters, and exact partition",
  "equivalence to the selected-resolution reference.",
  "",
  "| Species | Cells | Resolution | PCs | Clusters | ARI | NMI | Pass |",
  "|---|---:|---:|---:|---:|---:|---:|---|"
)
for (i in seq_len(nrow(validation))) {
  report <- c(report, paste0(
    "| ", validation$species[[i]], " | ", validation$n_cells[[i]],
    " | ", validation$selected_resolution[[i]],
    " | ", validation$selected_pcs[[i]],
    " | ", validation$n_clusters[[i]],
    " | ", sprintf("%.6f", validation$ARI[[i]]),
    " | ", sprintf("%.6f", validation$NMI[[i]]),
    " | ", validation$post_save_validation_pass[[i]], " |"
  ))
}
writeLines(
  report,
  file.path(output_dir, "verification", "post_save_rds_validation.md")
)

if (!all(validation$post_save_validation_pass)) {
  stop("One or more saved RDS files failed post-save validation")
}
message("All saved RDS files passed post-save validation")
