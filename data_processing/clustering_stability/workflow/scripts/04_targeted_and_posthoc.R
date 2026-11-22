suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppAnnoy)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(yaml)
  library(mclust)
  library(aricode)
})

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
workflow_root <- dirname(dirname(script_path))
source(file.path(workflow_root, "R", "common.R"))
source(file.path(workflow_root, "R", "metrics.R"))

cli <- parse_cli(
  file.path(workflow_root, "config", "config.yaml"),
  extra = list(
    mode = NULL,
    splits = NULL,
    resampling_mode = "formal",
    count_mode = NULL,
    run_scshc = NULL
  )
)
config <- load_config(cli$config)
initialize_run(config, cli$config)
mode <- match.arg(cli$mode, c("targeted-split", "annotation-recovery"))
count_mode <- cli$count_mode %||% config$count_input$primary_mode

if (identical(mode, "annotation-recovery")) {
  dirs <- ensure_stage_dirs(config, "10_posthoc_annotation", c("tables", "figures", "report", "logs"))
  freeze_path <- file.path(
    stage_dir(config, "08_stability_range", "tables"),
    "working_resolution_freeze.csv"
  )
  if (!file.exists(freeze_path)) {
    abort("No frozen resolution exists. Run 03_summarize_stability_range.R with --freeze-resolution= first.")
  }
  frozen <- utils::read.csv(freeze_path, stringsAsFactors = FALSE)
  if (nrow(frozen) != 1L) abort("working_resolution_freeze.csv must contain exactly one row.")
  if (!identical(as.character(frozen$input_rds_md5), file_md5(config$input$seurat_rds))) {
    abort("The frozen resolution belongs to a different input object.")
  }
  resolution <- frozen$resolution[[1L]]
  object <- readRDS(config$input$seurat_rds)
  annotation_column <- config$posthoc_annotation$annotation_column
  if (!annotation_column %in% colnames(object@meta.data)) {
    abort("Annotation column '%s' is absent from object metadata.", annotation_column)
  }
  annotations <- as.character(object[[annotation_column, drop = TRUE]])
  if (anyNA(annotations)) abort("Annotation column contains missing values.")
  assignment_path <- file.path(
    stage_dir(config, "02_full_resolution_scan", "tables"),
    "full_resolution_assignments.tsv.gz"
  )
  assignments <- read_table_auto(assignment_path)
  clusters <- as.character(assignments[[resolution_column(resolution)]])
  contingency <- as.data.frame(table(cluster = clusters, annotation = annotations), stringsAsFactors = FALSE)
  names(contingency)[[3L]] <- "cell_count"
  contingency <- contingency[contingency$cell_count > 0L, , drop = FALSE]
  cluster_total <- setNames(as.numeric(table(clusters)), names(table(clusters)))
  annotation_total <- setNames(as.numeric(table(annotations)), names(table(annotations)))
  contingency$fraction_within_cluster <- contingency$cell_count / cluster_total[contingency$cluster]
  contingency$fraction_within_annotation <- contingency$cell_count / annotation_total[contingency$annotation]
  cluster_metrics <- bind_rows_base(lapply(split(contingency, contingency$cluster), function(x) {
    probability <- x$fraction_within_cluster
    entropy <- -sum(probability * log(probability))
    data.frame(
      resolution = resolution,
      cluster = x$cluster[[1L]],
      cluster_size = sum(x$cell_count),
      dominant_annotation = x$annotation[[which.max(probability)]],
      purity = max(probability),
      normalized_entropy = if (nrow(x) > 1L) entropy / log(nrow(x)) else 0,
      annotation_count = nrow(x),
      stringsAsFactors = FALSE
    )
  }))
  global_metrics <- data.frame(
    resolution = resolution,
    n_clusters = length(unique(clusters)),
    n_annotations = length(unique(annotations)),
    cluster_annotation_ARI = mclust::adjustedRandIndex(clusters, annotations),
    cluster_annotation_NMI = aricode::NMI(clusters, annotations),
    cell_weighted_purity = stats::weighted.mean(cluster_metrics$purity, cluster_metrics$cluster_size),
    median_cluster_purity = safe_median(cluster_metrics$purity),
    stringsAsFactors = FALSE
  )
  composition_path <- file.path(dirs$tables, "cluster_annotation_composition.csv")
  cluster_path <- file.path(dirs$tables, "cluster_annotation_metrics.csv")
  global_path <- file.path(dirs$tables, "global_annotation_metrics.csv")
  write_csv_atomic(contingency, composition_path)
  write_csv_atomic(cluster_metrics, cluster_path)
  write_csv_atomic(global_metrics, global_path)
  save_plot(
    ggplot(contingency, aes(annotation, cluster, fill = fraction_within_cluster)) +
      geom_tile() + scale_fill_viridis_c(limits = c(0, 1)) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(x = "Original annotation", y = "Frozen cluster", fill = "Cluster fraction"),
    file.path(dirs$figures, "annotation_recovery_heatmap.pdf"),
    width = 10, height = 8
  )
  report_path <- file.path(dirs$report, "posthoc_annotation_report.md")
  writeLines(c(
    paste0("# Posthoc annotation recovery: resolution ", fmt_resolution(resolution)),
    "",
    "- Annotation was not used to select the stability range.",
    paste0("- Annotation column: `", annotation_column, "`"),
    paste0("- Cell-weighted cluster purity: ", round(global_metrics$cell_weighted_purity, 4)),
    paste0("- Cluster/annotation ARI: ", round(global_metrics$cluster_annotation_ARI, 4)),
    paste0("- Cluster/annotation NMI: ", round(global_metrics$cluster_annotation_NMI, 4))
  ), report_path)
  write_stage_manifest(
    config, "10_posthoc_annotation", resolution,
    inputs = c(freeze_path, assignment_path),
    outputs = c(composition_path, cluster_path, global_path, report_path),
    details = list(annotation_column = annotation_column, annotation_used_for_selection = FALSE)
  )
  cat("Annotation recovery completed for resolution", fmt_resolution(resolution), "\n")
  quit(save = "no", status = 0)
}

dirs <- ensure_stage_dirs(config, "09_targeted_split_validation", c("marker", "silhouette", "scshc", "logs"))
split_path <- file.path(stage_dir(config, "07_resolution_transitions", "tables"), "split_persistence.csv")
edge_path <- file.path(stage_dir(config, "07_resolution_transitions", "tables"), "parent_daughter_edges.csv")
if (!file.exists(split_path)) abort("Run 03_summarize_stability_range.R before targeted split analysis.")
splits <- read_table_auto(split_path)
edges <- read_table_auto(edge_path)
selected_ids <- if (!is.null(cli$splits)) {
  trimws(strsplit(cli$splits, ",", fixed = TRUE)[[1L]])
} else {
  as.character(unlist(config$targeted$selected_splits))
}
if (!length(selected_ids)) abort("No split IDs selected. Pass --splits=id1,id2 or set targeted.selected_splits.")
missing <- setdiff(selected_ids, splits$split_id)
if (length(missing)) abort("Unknown split IDs: %s", paste(missing, collapse = ", "))
selected <- splits[match(selected_ids, splits$split_id), , drop = FALSE]

loaded <- load_and_validate_object(config)
object <- loaded$object
config$._resolved_primary_pcs <- loaded$primary_pcs
assignments <- read_table_auto(file.path(
  stage_dir(config, "02_full_resolution_scan", "tables"),
  "full_resolution_assignments.tsv.gz"
))
expression <- get_assay_matrix(object, config$marker$assay, config$marker$layer)
marker_rows <- list()
silhouette_rows <- list()
membership_rows <- list()

for (i in seq_len(nrow(selected))) {
  split_row <- selected[i, , drop = FALSE]
  current_edges <- edges[
    abs(edges$lower_resolution - split_row$lower_resolution) < 1e-10 &
      abs(edges$upper_resolution - split_row$upper_resolution) < 1e-10 &
      edges$parent_cluster == as.character(split_row$parent_cluster) &
      edges$parent_fraction >= 0.05, , drop = FALSE
  ]
  parent_labels <- as.character(assignments[[resolution_column(split_row$lower_resolution)]])
  daughter_labels <- as.character(assignments[[resolution_column(split_row$upper_resolution)]])
  cells <- assignments$cell_id[
    parent_labels == as.character(split_row$parent_cluster) &
      daughter_labels %in% current_edges$daughter_cluster
  ]
  daughters <- daughter_labels[match(cells, assignments$cell_id)]
  membership_rows[[i]] <- data.frame(
    split_id = split_row$split_id, cell_id = cells, daughter_cluster = daughters,
    stringsAsFactors = FALSE
  )
  subset_object <- subset(object, cells = cells)
  Idents(subset_object) <- factor(daughters[match(colnames(subset_object), cells)])
  markers <- FindAllMarkers(
    subset_object,
    assay = config$marker$assay,
    slot = config$marker$layer,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = config$marker$min_pct,
    logfc.threshold = config$marker$log2fc_threshold,
    verbose = FALSE
  )
  if (nrow(markers)) {
    effect_column <- if ("avg_log2FC" %in% names(markers)) "avg_log2FC" else "avg_logFC"
    markers <- bind_rows_base(lapply(split(markers, markers$cluster), function(x) {
      head(x[order(x[[effect_column]], decreasing = TRUE), , drop = FALSE], config$marker$top_n)
    }))
    markers$split_id <- split_row$split_id
    marker_rows[[i]] <- markers
  }

  set.seed(config$compute$seed + i)
  sampled <- if (length(cells) > config$targeted$silhouette_max_cells) {
    unlist(lapply(split(cells, daughters), function(x) {
      sample(x, max(2L, round(config$targeted$silhouette_max_cells * length(x) / length(cells))))
    }), use.names = FALSE)
  } else cells
  sampled <- unique(sampled)
  sampled_labels <- daughter_labels[match(sampled, assignments$cell_id)]
  if (length(unique(sampled_labels)) >= 2L) {
    distance <- stats::dist(
      loaded$embedding[
        sampled,
        seq_len(primary_pcs_from_config(config)),
        drop = FALSE
      ]
    )
    values <- cluster::silhouette(as.integer(factor(sampled_labels)), distance)
    silhouette_rows[[i]] <- data.frame(
      split_id = split_row$split_id,
      lower_resolution = split_row$lower_resolution,
      upper_resolution = split_row$upper_resolution,
      parent_cluster = split_row$parent_cluster,
      daughter_count = length(unique(sampled_labels)),
      evaluated_cells = length(sampled),
      mean_silhouette = mean(values[, "sil_width"]),
      median_silhouette = stats::median(values[, "sil_width"]),
      negative_silhouette_fraction = mean(values[, "sil_width"] < 0),
      stringsAsFactors = FALSE
    )
  }
}
markers <- bind_rows_base(marker_rows)
silhouette <- bind_rows_base(silhouette_rows)
memberships <- bind_rows_base(membership_rows)
marker_output <- file.path(dirs$marker, "sister_marker_details.tsv.gz")
silhouette_output <- file.path(dirs$silhouette, "silhouette_summary.csv")
membership_output <- file.path(dirs$marker, "split_cell_membership.tsv.gz")
write_tsv_gz_atomic(markers, marker_output)
write_csv_atomic(silhouette, silhouette_output)
write_tsv_gz_atomic(memberships, membership_output)

run_scshc <- if (is.null(cli$run_scshc)) isTRUE(config$targeted$run_scshc) else as_flag(cli$run_scshc)
scshc_output <- NULL
if (run_scshc) {
  local_library <- normalizePath(config$targeted$scshc_library, mustWork = TRUE)
  .libPaths(c(local_library, .libPaths()))
  assert_packages(c("scSHC", "scry"))
  description <- packageDescription("scSHC", lib.loc = local_library)
  if (!identical(description$RemoteSha %||% NA_character_, config$targeted$scshc_commit)) {
    abort("Installed scSHC commit does not match config.targeted.scshc_commit.")
  }
  prepared <- prepare_counts(loaded$counts, config, "soupx_integer_full_pipeline")
  scshc_rows <- lapply(seq_len(nrow(selected)), function(i) {
    split_row <- selected[i, , drop = FALSE]
    membership <- memberships[memberships$split_id == split_row$split_id, , drop = FALSE]
    data <- prepared$counts[, membership$cell_id, drop = FALSE]
    labels <- membership$daughter_cluster
    if (min(table(labels)) <= 20L) {
      return(data.frame(split_id = split_row$split_id, tested = FALSE, reason = "daughter_le_20_cells"))
    }
    deviance <- scry::devianceFeatureSelection(data)
    genes <- rownames(data)[head(order(deviance, decreasing = TRUE), config$targeted$scshc_num_features)]
    pseudobulk <- rowsum(t(as.matrix(data[genes, , drop = FALSE])), group = labels)
    pseudobulk <- sweep(pseudobulk, 1L, rowSums(pseudobulk), "/")
    tree <- stats::hclust(stats::dist(pseudobulk), method = "ward.D2")
    groups <- stats::cutree(tree, k = 2L)
    ids_one <- which(labels %in% names(groups)[groups == 1L])
    ids_two <- which(labels %in% names(groups)[groups == 2L])
    batch <- rep("1", ncol(data))
    set.seed(config$compute$seed + i)
    raw_p <- scSHC:::test_split(
      data, ids_one, ids_two, genes,
      as.integer(config$targeted$scshc_num_pcs), batch,
      as.numeric(config$targeted$scshc_alpha),
      as.integer(min(config$targeted$scshc_cores, config$compute$max_workers)),
      posthoc = TRUE
    )
    data.frame(
      split_id = split_row$split_id,
      tested = TRUE,
      reason = "",
      root_group_one = paste(names(groups)[groups == 1L], collapse = ";"),
      root_group_two = paste(names(groups)[groups == 2L], collapse = ";"),
      raw_p = raw_p,
      significant = raw_p < config$targeted$scshc_alpha,
      stringsAsFactors = FALSE
    )
  })
  scshc_results <- bind_rows_base(scshc_rows)
  scshc_output <- file.path(dirs$scshc, "scshc_root_split_results.csv")
  write_csv_atomic(scshc_results, scshc_output)
}

outputs <- c(marker_output, silhouette_output, membership_output)
if (!is.null(scshc_output)) outputs <- c(outputs, scshc_output)
write_stage_manifest(
  config, "09_targeted_split_validation",
  unique(c(selected$lower_resolution, selected$upper_resolution)),
  inputs = c(split_path, edge_path),
  outputs = outputs,
  details = list(split_count = nrow(selected), scshc = run_scshc)
)
cat("Targeted split validation completed for", nrow(selected), "splits.\n")
