suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppAnnoy)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(yaml)
  library(mclust)
  library(aricode)
  library(clustree)
})

workflow_root <- workflow_root_from_script <- NULL
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
workflow_root <- dirname(dirname(script_path))
source(file.path(workflow_root, "R", "common.R"))
source(file.path(workflow_root, "R", "metrics.R"))

cli <- parse_cli(
  file.path(workflow_root, "config", "config.yaml"),
  extra = list(resolutions = NULL)
)
config <- load_config(cli$config)
initialize_run(config, cli$config)
paths <- ensure_stage_dirs(config, "01_input_audit")
scan_paths <- ensure_stage_dirs(config, "02_full_resolution_scan")
transition_paths <- ensure_stage_dirs(config, "07_resolution_transitions")

log_path <- file.path(scan_paths$logs, "scan_resolution_grid.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Started:", format(Sys.time()), "\n")
cat("Config:", normalizePath(cli$config), "\n")

resolutions <- if (is.null(cli$resolutions)) {
  config$clustering$resolution_grid
} else {
  validate_resolution_subset(as_number_vector(cli$resolutions), config$clustering$resolution_grid, "--resolutions")
}
if (!length(resolutions)) abort("No scan resolutions were supplied.")
existing_manifest <- file.path(stage_dir(config, "02_full_resolution_scan"), "stage_manifest.csv")
existing_assignment <- file.path(
  stage_dir(config, "02_full_resolution_scan", "tables"),
  "full_resolution_assignments.tsv.gz"
)
if (file.exists(existing_manifest) && file.exists(existing_assignment) && !as_flag(cli$overwrite)) {
  previous <- utils::read.csv(existing_manifest, stringsAsFactors = FALSE)
  if (!identical(previous$resolution_set_id[[1L]], resolution_set_id(resolutions))) {
    abort("The scan stage already exists with a different resolution set. Use a new run_id.")
  }
  cat("Scan stage already completed with the same resolution set; nothing was overwritten.\n")
  cat("Pass --overwrite to recompute it explicitly.\n")
  quit(save = "no", status = 0)
}

loaded <- load_and_validate_object(config)
object <- loaded$object
config$._resolved_primary_pcs <- loaded$primary_pcs
primary_pcs <- primary_pcs_from_config(config)
counts_info <- prepare_counts(loaded$counts, config)
counts_audit <- counts_info$audit
cat("Resolved primary PCs:", primary_pcs, "\n")

metadata <- read_external_metadata(config)
cell_column <- intersect(c("cell_id", "Cell", "barcode", "cell"), names(metadata))
metadata_cell_ids <- if (length(cell_column)) {
  as.character(metadata[[cell_column[[1L]]]])
} else {
  rownames(metadata)
}
metadata_matches <- length(metadata_cell_ids) == ncol(object) &&
  setequal(metadata_cell_ids, colnames(object))
input_audit <- data.frame(
  item = c(
    "seurat_class", "n_cells", "n_features_counts", "SCT_assay",
    "count_assay", "reduction", "pca_dimensions", "metadata_cell_set",
    "reduction_assay", "selected_primary_pcs", "fractional_count_fraction",
    "save_full_seurat_objects"
  ),
  value = c(
    paste(class(object), collapse = ","),
    ncol(object),
    nrow(loaded$counts),
    config$input$sct_assay,
    config$input$count_assay,
    config$input$reduction,
    ncol(loaded$embedding),
    metadata_matches,
    object[[config$input$reduction]]@assay.used,
    primary_pcs,
    counts_audit$noninteger_nonzero_fraction,
    config$resampling$save_seurat_objects
  ),
  status = c(
    rep("OK", 10), "INFO",
    ifelse(isTRUE(config$resampling$save_seurat_objects), "WARNING", "OK")
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(input_audit, file.path(paths$tables, "input_validation.csv"))
write_csv_atomic(counts_audit, file.path(paths$tables, "count_matrix_audit.csv"))
canonical_audit_path <- file.path(paths$tables, "canonical_preprocessing_audit.csv")
write_csv_atomic(canonical_audit_table(loaded$canonical), canonical_audit_path)
canonical_elbow_path <- NULL
if (!is.null(loaded$canonical)) {
  elbow_stdev <- Seurat::Stdev(object, reduction = config$input$reduction)
  elbow_data <- data.frame(
    pc = seq_len(min(50L, length(elbow_stdev))),
    standard_deviation = elbow_stdev[seq_len(min(50L, length(elbow_stdev)))]
  )
  canonical_elbow_path <- file.path(paths$figures, "canonical_pca_elbow.pdf")
  save_plot(
    ggplot(elbow_data, aes(pc, standard_deviation)) +
      geom_line(linewidth = 0.4) +
      geom_point(size = 1.2) +
      geom_vline(xintercept = primary_pcs, color = "red", linetype = 2) +
      theme_bw() +
      labs(x = "PC", y = "Standard deviation"),
    canonical_elbow_path,
    width = 7, height = 4.5
  )
}
core_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "yaml", "ggplot2",
  "mclust", "aricode", "data.table", "cluster", "digest"
)
package_table <- data.frame(
  package = core_packages,
  available = vapply(core_packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(core_packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
write_csv_atomic(package_table, file.path(paths$tables, "package_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(paths$logs, "session_info.txt"))
if (!metadata_matches) abort("Metadata cell IDs do not match the Seurat object.")

set.seed(config$compute$seed)
DefaultAssay(object) <- config$input$sct_assay
nn_name <- "stability_scan_nn"
snn_name <- "stability_scan_snn"
object <- FindNeighbors(
  object,
  reduction = config$input$reduction,
  dims = seq_len(primary_pcs),
  k.param = config$clustering$primary_k,
  graph.name = c(nn_name, snn_name),
  return.neighbor = FALSE,
  verbose = FALSE
)

assignment <- data.frame(cell_id = colnames(object), stringsAsFactors = FALSE)
size_rows <- vector("list", length(resolutions))
for (i in seq_along(resolutions)) {
  resolution <- resolutions[[i]]
  column <- resolution_column(resolution)
  cat("Clustering resolution", fmt_resolution(resolution), "\n")
  object <- FindClusters(
    object,
    graph.name = snn_name,
    resolution = resolution,
    algorithm = config$clustering$algorithm,
    cluster.name = column,
    random.seed = config$compute$seed,
    verbose = FALSE
  )
  assignment[[column]] <- as.character(object[[column, drop = TRUE]])
  current <- cluster_size_summary(
    assignment[[column]],
    small_threshold = config$clustering$small_cluster_cells,
    min_threshold = config$clustering$min_cluster_cells
  )
  current$resolution <- resolution
  current$small_cell_fraction <- sum(table(assignment[[column]])[
    table(assignment[[column]]) < config$clustering$small_cluster_cells
  ]) / nrow(assignment)
  current$largest_cluster_fraction <- max(table(assignment[[column]])) / nrow(assignment)
  size_rows[[i]] <- current
}

size_summary <- bind_rows_base(size_rows)
size_summary <- size_summary[c("resolution", setdiff(names(size_summary), "resolution"))]
adjacent <- adjacent_partition_metrics(assignment, resolutions)
equivalent <- equivalent_resolution_intervals(adjacent, resolutions)
edges <- parent_daughter_edges(assignment, resolutions)
edge_groups <- interaction(edges$upper_resolution, edges$daughter_cluster, drop = TRUE)
edges$clustree_core_edge <- ave(
  edges$daughter_fraction,
  edge_groups,
  FUN = function(value) value == max(value)
)
edges$clustree_displayed_edge <- edges$daughter_fraction >= 0.05
configured_pilot <- validate_resolution_subset(
  config$resampling$pilot_resolutions, resolutions, "resampling.pilot_resolutions"
)
pilot_candidates <- data.frame(
  resolution = configured_pilot,
  source = rep("configured_pilot", length(configured_pilot)),
  equivalent_interval_id = vapply(configured_pilot, function(value) {
    hit <- equivalent$equivalent_interval_id[
      value >= equivalent$resolution_min - 1e-10 &
        value <= equivalent$resolution_max + 1e-10
    ]
    if (length(hit)) hit[[1L]] else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
scan_driven_candidates <- select_scan_pilot_candidates(
  size_summary,
  adjacent,
  config$scan_candidate_selection
)

assignment_path <- file.path(scan_paths$tables, "full_resolution_assignments.tsv.gz")
size_path <- file.path(scan_paths$tables, "resolution_cluster_size_summary.csv")
adjacent_path <- file.path(transition_paths$tables, "adjacent_resolution_similarity.csv")
equivalent_path <- file.path(transition_paths$tables, "equivalent_resolution_intervals.csv")
edge_path <- file.path(transition_paths$tables, "parent_daughter_edges.csv")
clustree_edge_path <- file.path(transition_paths$tables, "clustree_correspondence_edges.csv")
pilot_path <- file.path(scan_paths$tables, "pilot_candidate_resolutions.csv")
scan_candidate_path <- file.path(
  scan_paths$tables,
  "scan_driven_pilot_candidates.csv"
)
write_tsv_gz_atomic(assignment, assignment_path)
write_csv_atomic(size_summary, size_path)
write_csv_atomic(adjacent, adjacent_path)
write_csv_atomic(equivalent, equivalent_path)
write_csv_atomic(edges, edge_path)
write_csv_atomic(
  edges[c(
    "lower_resolution", "upper_resolution", "parent_cluster", "daughter_cluster",
    "overlap_n", "parent_fraction", "daughter_fraction",
    "clustree_core_edge", "clustree_displayed_edge"
  )],
  clustree_edge_path
)
write_csv_atomic(pilot_candidates, pilot_path)
write_csv_atomic(scan_driven_candidates, scan_candidate_path)

clustree_data <- assignment[vapply(resolutions, resolution_column, character(1))]
colnames(clustree_data) <- paste0(
  "scan_res.", vapply(resolutions, fmt_resolution, character(1))
)
make_clustree <- function(values, filename, width, height) {
  columns <- paste0("scan_res.", vapply(values, fmt_resolution, character(1)))
  plot <- clustree::clustree(
    clustree_data[columns],
    prefix = "scan_res.",
    prop_filter = 0.05,
    count_filter = config$clustering$min_cluster_cells,
    layout = "sugiyama",
    use_core_edges = TRUE,
    highlight_core = TRUE,
    node_text_size = 2.2,
    node_size_range = c(3, 11),
    show_axis = TRUE
  )
  save_plot(plot, file.path(transition_paths$figures, filename), width, height)
}
make_clustree(resolutions, "clustree_full_resolution_scan.pdf", 20, 28)
clustree_breaks <- list(
  low = resolutions[resolutions <= 1.5],
  middle = resolutions[resolutions >= 1.4 & resolutions <= 3.0],
  high = resolutions[resolutions >= 2.9]
)
for (name in names(clustree_breaks)) {
  values <- clustree_breaks[[name]]
  if (length(values) >= 2L) {
    make_clustree(
      values,
      paste0("clustree_", name, "_resolution_scan.pdf"),
      16,
      max(10, 0.85 * length(values))
    )
  }
}

save_plot(
  ggplot(size_summary, aes(resolution, n_clusters)) +
    geom_line(linewidth = 0.4) + geom_point(size = 1.5) + theme_bw() +
    labs(x = "Resolution", y = "Number of clusters"),
  file.path(scan_paths$figures, "cluster_count_by_resolution.pdf")
)
save_plot(
  ggplot(size_summary, aes(resolution, small_cell_fraction)) +
    geom_line(linewidth = 0.4) + geom_point(size = 1.5) + theme_bw() +
    labs(x = "Resolution", y = "Fraction of cells in small clusters"),
  file.path(scan_paths$figures, "small_cell_fraction_by_resolution.pdf")
)
if (nrow(adjacent)) {
  long <- rbind(
    data.frame(upper_resolution = adjacent$upper_resolution, metric = "ARI", value = adjacent$adjacent_ARI),
    data.frame(upper_resolution = adjacent$upper_resolution, metric = "NMI", value = adjacent$adjacent_NMI)
  )
  save_plot(
    ggplot(long, aes(upper_resolution, value, color = metric)) +
      geom_line(linewidth = 0.4) + geom_point(size = 1.5) + theme_bw() +
      labs(x = "Upper resolution", y = "Adjacent partition similarity", color = NULL),
    file.path(transition_paths$figures, "adjacent_partition_similarity.pdf")
  )
}

write_stage_manifest(
  config, "01_input_audit", resolutions,
  inputs = c(config$input$seurat_rds, config$input$metadata_csv),
  outputs = c(
    file.path(paths$tables, "input_validation.csv"),
    file.path(paths$tables, "count_matrix_audit.csv"),
    file.path(paths$tables, "package_versions.csv"),
    canonical_audit_path,
    canonical_elbow_path
  ),
  details = list(
    annotation_used = FALSE,
    full_seurat_saved = FALSE,
    canonical_rebuilt = !is.null(loaded$canonical),
    selected_pcs = primary_pcs
  )
)
write_stage_manifest(
  config, "02_full_resolution_scan", resolutions,
  inputs = c(config$input$seurat_rds),
  outputs = c(assignment_path, size_path, pilot_path, scan_candidate_path),
  details = list(
    pcs = primary_pcs,
    k = config$clustering$primary_k,
    leiden_seed = config$compute$seed,
    seed_mode = "shared_across_resolutions"
  )
)
write_stage_manifest(
  config, "07_resolution_transitions", resolutions,
  inputs = c(assignment_path),
  outputs = c(
    adjacent_path, equivalent_path, edge_path, clustree_edge_path,
    file.path(transition_paths$figures, "clustree_full_resolution_scan.pdf"),
    file.path(transition_paths$figures, "clustree_low_resolution_scan.pdf"),
    file.path(transition_paths$figures, "clustree_middle_resolution_scan.pdf"),
    file.path(transition_paths$figures, "clustree_high_resolution_scan.pdf")
  ),
  details = list(equivalent_intervals = nrow(equivalent))
)
cat("Finished:", format(Sys.time()), "\n")
cat("No Seurat object was saved.\n")
