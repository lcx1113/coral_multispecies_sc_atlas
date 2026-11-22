run_marker_filter_step <- function(sample_row, soupx_obj, sample_dir, config, marker_cutoff) {
  sample_id <- sample_row$sample_id[[1]]
  marker_first_dir <- qc_step_dir(sample_dir, "marker_filter")
  if (qc_use_numbered_dirs()) {
    dir_before <- file.path(marker_first_dir, "03a_before_marker_filter")
    dir_marker <- file.path(marker_first_dir, "03b_marker_tables")
    dir_removed <- file.path(marker_first_dir, "03c_removed_cells")
    dir_after <- file.path(marker_first_dir, "03d_after_marker_filter")
  } else {
    dir_before <- file.path(marker_first_dir, "before_marker_filter")
    dir_marker <- file.path(marker_first_dir, "marker_filter")
    dir_removed <- file.path(marker_first_dir, "removed_cells")
    dir_after <- file.path(marker_first_dir, "after_marker_filter")
  }
  for (d in c(dir_before, dir_marker, dir_removed, dir_after)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  cached_status <- file.path(dir_marker, "cell_removal_status.tsv")
  cached_kept <- file.path(dir_marker, "kept_after_marker_filter_barcodes.txt")
  if (isTRUE(config$output$resume_from_existing) &&
      !isTRUE(config$output$overwrite_existing) &&
      file.exists(cached_status) && file.exists(cached_kept)) {
    log_msg(sample_id, ": reusing existing marker-first cell-removal status.")
    kept <- readLines(cached_kept)
    return(list(obj = subset(soupx_obj, cells = kept), status = read_tsv(cached_status),
                marker_first_dir = marker_first_dir))
  }

  log_msg(sample_id, ": RNA workflow before marker filtering.")
  obj <- soupx_obj
  DefaultAssay(obj) <- "RNA"
  obj <- join_assay_layers_if_needed(obj, "RNA")
  obj@reductions <- list()
  obj@graphs <- list()
  obj@neighbors <- list()
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                       scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, assay = "RNA", selection.method = "vst",
                              nfeatures = config$marker_filter$rna_nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, assay = "RNA", features = VariableFeatures(obj), verbose = FALSE)
  pre_dims <- safe_dims(obj, 20)
  obj <- RunPCA(obj, assay = "RNA", features = VariableFeatures(obj),
                npcs = max(pre_dims), verbose = FALSE)
  neighbor_dims <- pre_dims[pre_dims <= min(10, max(pre_dims))]
  obj <- FindNeighbors(obj, dims = neighbor_dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = config$marker_filter$rna_resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = pre_dims, verbose = FALSE)
  original_cluster_col <- paste0("RNA_snn_res.", config$marker_filter$rna_resolution)
  if (!original_cluster_col %in% colnames(obj@meta.data)) original_cluster_col <- "seurat_clusters"
  obj$original_cluster <- as.character(obj@meta.data[[original_cluster_col]])
  original_levels <- numeric_order(obj$original_cluster)

  rna_meta <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(original_cluster = as.character(original_cluster))
  rna_umap <- as.data.frame(Embeddings(obj, "umap"))
  colnames(rna_umap)[1:2] <- c("UMAP_1", "UMAP_2")
  rna_umap$cell <- rownames(rna_umap)
  rna_umap <- dplyr::left_join(
    rna_umap,
    rna_meta[, intersect(c("cell", "original_cluster", "nCount_RNA", "nFeature_RNA", "percent.mito"),
                         colnames(rna_meta))],
    by = "cell"
  )
  write_tsv(rna_meta, file.path(dir_before, "RNA_cell_metadata_before_marker_filter.tsv"))
  write_tsv(rna_umap, file.path(dir_before, "RNA_UMAP_coordinates_before_marker_filter.tsv"))
  p_rna <- ggplot2::ggplot(
    rna_umap,
    ggplot2::aes(UMAP_1, UMAP_2, colour = factor(original_cluster, levels = original_levels))
  ) +
    ggplot2::geom_point(size = 0.25, alpha = 0.75) +
    ggplot2::labs(title = paste0(sample_id, " RNA UMAP before marker filter"),
                  colour = "Original") +
    ggplot2::theme_classic()
  plot_and_save(p_rna, file.path(dir_before, "RNA_UMAP_by_original_cluster_before_marker_filter.pdf"),
                width = 8, height = 7)

  log_msg(sample_id, ": SCT Leiden clustering for marker-poor cluster detection.")
  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(obj, vst.flavor = config$marker_filter$sct_vst_flavor,
                     do.scale = TRUE, do.center = TRUE,
                     return.only.var.genes = TRUE, verbose = FALSE)
  obj <- RunPCA(obj, assay = "SCT", npcs = min(100, ncol(obj) - 1),
                reduction.name = "pca", verbose = FALSE)
  selected_pc <- select_sct_pcs(obj)
  selected_dims <- seq_len(selected_pc)
  write_tsv(data.frame(selected_pc = selected_pc, dims = paste(selected_dims, collapse = ",")),
            file.path(dir_before, "SCT_PCA_selected_dims.tsv"))
  grDevices::pdf(file.path(dir_before, "SCT_PCA_elbow.pdf"), width = 7, height = 5)
  plot(seq_along(obj@reductions$pca@stdev), obj@reductions$pca@stdev,
       type = "b", xlab = "PC", ylab = "Standard deviation",
       main = paste0(sample_id, " SCT PCA elbow"))
  abline(v = selected_pc, col = "red", lty = 2)
  grDevices::dev.off()

  obj <- FindNeighbors(obj, dims = selected_dims, reduction = "pca", verbose = FALSE)
  obj <- FindClusters(obj, resolution = config$marker_filter$sct_leiden_resolution,
                      cluster.name = "leiden", algorithm = 4, method = "igraph",
                      verbose = FALSE)
  if (!"leiden" %in% colnames(obj@meta.data)) obj$leiden <- as.character(obj$seurat_clusters)
  obj$leiden <- as.character(obj$leiden)
  obj <- RunUMAP(obj, reduction = "pca", dims = selected_dims, verbose = FALSE)

  pre_marker_meta <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(original_cluster = as.character(original_cluster),
                  leiden = as.character(leiden))
  write_tsv(pre_marker_meta, file.path(dir_before, "before_marker_filter_cell_metadata.tsv"))

  pre_umap <- as.data.frame(Embeddings(obj, "umap"))
  colnames(pre_umap)[1:2] <- c("UMAP_1", "UMAP_2")
  pre_umap$cell <- rownames(pre_umap)
  pre_umap <- dplyr::left_join(
    pre_umap,
    pre_marker_meta[, intersect(c("cell", "original_cluster", "leiden",
                                  "nCount_RNA", "nFeature_RNA", "percent.mito"),
                                colnames(pre_marker_meta))],
    by = "cell"
  )
  write_tsv(pre_umap, file.path(dir_before, "UMAP_coordinates_before_marker_filter.tsv"))
  leiden_levels <- numeric_order(pre_marker_meta$leiden)

  p_umap_original <- ggplot2::ggplot(
    pre_umap,
    ggplot2::aes(UMAP_1, UMAP_2, colour = factor(original_cluster, levels = original_levels))
  ) +
    ggplot2::geom_point(size = 0.25, alpha = 0.75) +
    ggplot2::labs(title = paste0(sample_id, " before marker filter: original cluster"),
                  colour = "Original") +
    ggplot2::theme_classic()
  plot_and_save(p_umap_original, file.path(dir_before, "UMAP_by_original_cluster.pdf"),
                width = 8, height = 7)
  p_umap_leiden <- ggplot2::ggplot(
    pre_umap,
    ggplot2::aes(UMAP_1, UMAP_2, colour = factor(leiden, levels = leiden_levels))
  ) +
    ggplot2::geom_point(size = 0.25, alpha = 0.75) +
    ggplot2::labs(title = paste0(sample_id, " before marker filter: Leiden res",
                                 config$marker_filter$sct_leiden_resolution),
                  colour = "Leiden") +
    ggplot2::theme_classic()
  plot_and_save(p_umap_leiden, file.path(dir_before, "UMAP_by_leiden_res4.pdf"),
                width = 9, height = 8)

  DefaultAssay(obj) <- "RNA"
  obj <- join_assay_layers_if_needed(obj, "RNA")
  original_marker_result <- run_find_all_markers(
    obj, "original_cluster",
    paste0(sample_id, " original clusters before marker filtering"),
    config
  )
  obj <- original_marker_result$obj
  original_markers_raw <- original_marker_result$raw
  original_markers_pass <- original_marker_result$pass
  if (isTRUE(config$output$save_raw_marker_gz)) {
    write_tsv_gz(original_markers_raw,
                 file.path(dir_marker, "original_cluster_markers_one_vs_rest_raw.tsv.gz"))
  }
  write_tsv(original_markers_pass,
            file.path(dir_marker, "original_cluster_markers_one_vs_rest_pass.tsv"))
  original_summary <- build_cluster_marker_summary(obj, "original_cluster",
                                                   original_markers_pass, marker_cutoff)
  write_tsv(original_summary, file.path(dir_marker, "original_cluster_marker_summary.tsv"))

  leiden_marker_result <- run_find_all_markers(
    obj, "leiden",
    paste0(sample_id, " Leiden clusters before marker filtering"),
    config
  )
  obj <- leiden_marker_result$obj
  leiden_markers_raw <- leiden_marker_result$raw
  leiden_markers_pass <- leiden_marker_result$pass
  if (isTRUE(config$output$save_raw_marker_gz)) {
    write_tsv_gz(leiden_markers_raw,
                 file.path(dir_marker, "leiden_markers_one_vs_rest_raw.tsv.gz"))
  }
  write_tsv(leiden_markers_pass,
            file.path(dir_marker, "leiden_markers_one_vs_rest_pass.tsv"))
  leiden_summary <- build_cluster_marker_summary(obj, "leiden", leiden_markers_pass, marker_cutoff)
  write_tsv(leiden_summary, file.path(dir_marker, "leiden_cluster_marker_summary.tsv"))

  flagged_leiden <- leiden_summary %>%
    dplyr::filter(flagged_marker_poor) %>%
    dplyr::mutate(flag_reason = paste0("n_markers_one_vs_rest <= ", marker_cutoff))
  flagged_leiden_ids <- as.character(flagged_leiden$cluster)
  write_tsv(flagged_leiden, file.path(dir_marker, "flagged_leiden_clusters_marker_poor.tsv"))

  cell_removal_status <- pre_marker_meta %>%
    dplyr::mutate(
      remove_by_marker = leiden %in% flagged_leiden_ids,
      marker_filter_reason = ifelse(remove_by_marker,
                                    paste0("leiden n_markers_one_vs_rest <= ", marker_cutoff),
                                    "kept")
    )
  removed_cells <- cell_removal_status$cell[cell_removal_status$remove_by_marker]
  kept_cells <- cell_removal_status$cell[!cell_removal_status$remove_by_marker]
  write_tsv(cell_removal_status, file.path(dir_marker, "cell_removal_status.tsv"))
  writeLines(removed_cells, file.path(dir_marker, "removed_by_marker_barcodes.txt"))
  writeLines(kept_cells, file.path(dir_marker, "kept_after_marker_filter_barcodes.txt"))
  writeLines(removed_cells, file.path(dir_removed, "removed_by_marker_barcodes.txt"))

  crosswalk <- cell_removal_status %>%
    dplyr::count(leiden, original_cluster, name = "n_cells") %>%
    dplyr::group_by(leiden) %>%
    dplyr::mutate(
      leiden_cluster_size = sum(n_cells),
      fraction_of_leiden_cluster = n_cells / leiden_cluster_size,
      dominant_original_cluster = original_cluster[which.max(n_cells)]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(original_cluster) %>%
    dplyr::mutate(
      original_cluster_size = sum(n_cells),
      fraction_of_original_cluster = n_cells / original_cluster_size
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(factor(leiden, levels = leiden_levels), dplyr::desc(n_cells))
  write_tsv(crosswalk, file.path(dir_marker, "leiden_vs_original_cluster_crosswalk.tsv"))

  crosswalk_integrated <- crosswalk %>%
    dplyr::left_join(
      leiden_summary %>%
        dplyr::rename(
          leiden_n_cells = n_cells,
          leiden_median_nUMI = median_nUMI,
          leiden_median_nFeature = median_nFeature,
          leiden_median_percent_mito = median_percent_mito
        ),
      by = c("leiden" = "cluster")
    )
  write_tsv(crosswalk_integrated,
            file.path(dir_marker, "leiden_marker_original_cluster_integrated.tsv"))

  original_count_change <- cell_removal_status %>%
    dplyr::count(original_cluster, remove_by_marker, name = "n_cells") %>%
    dplyr::mutate(state = ifelse(remove_by_marker, "Removed", "Kept")) %>%
    dplyr::select(-remove_by_marker) %>%
    tidyr::pivot_wider(names_from = state, values_from = n_cells, values_fill = 0)
  if (!"Kept" %in% colnames(original_count_change)) original_count_change$Kept <- 0L
  if (!"Removed" %in% colnames(original_count_change)) original_count_change$Removed <- 0L
  original_count_change <- original_count_change %>%
    dplyr::mutate(
      cells_before = Kept + Removed,
      cells_removed = Removed,
      cells_after = Kept,
      fraction_removed = cells_removed / cells_before
    ) %>%
    dplyr::select(original_cluster, cells_before, cells_removed, cells_after, fraction_removed) %>%
    dplyr::arrange(factor(original_cluster, levels = original_levels))
  write_tsv(original_count_change,
            file.path(dir_marker, "original_cluster_cell_count_change.tsv"))

  leiden_count_change <- cell_removal_status %>%
    dplyr::count(leiden, remove_by_marker, name = "n_cells") %>%
    dplyr::mutate(state = ifelse(remove_by_marker, "Removed", "Kept")) %>%
    dplyr::select(-remove_by_marker) %>%
    tidyr::pivot_wider(names_from = state, values_from = n_cells, values_fill = 0)
  if (!"Kept" %in% colnames(leiden_count_change)) leiden_count_change$Kept <- 0L
  if (!"Removed" %in% colnames(leiden_count_change)) leiden_count_change$Removed <- 0L
  leiden_count_change <- leiden_count_change %>%
    dplyr::mutate(
      cells_before = Kept + Removed,
      cells_removed = Removed,
      cells_after = Kept,
      fraction_removed = cells_removed / cells_before
    ) %>%
    dplyr::select(leiden, cells_before, cells_removed, cells_after, fraction_removed) %>%
    dplyr::arrange(factor(leiden, levels = leiden_levels))
  write_tsv(leiden_count_change,
            file.path(dir_marker, "leiden_cluster_cell_count_change.tsv"))

  pre_umap_status <- pre_umap %>%
    dplyr::left_join(cell_removal_status[, c("cell", "remove_by_marker", "marker_filter_reason")],
                     by = "cell")
  write_tsv(pre_umap_status, file.path(dir_marker, "umap_coordinates_removal_status.tsv"))
  write_tsv(pre_umap_status %>% dplyr::filter(remove_by_marker),
            file.path(dir_removed, "removed_cells_umap_coordinates.tsv"))
  p_removed <- ggplot2::ggplot(pre_umap_status, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(colour = "grey82", size = 0.25) +
    ggplot2::geom_point(
      data = pre_umap_status %>% dplyr::filter(remove_by_marker),
      ggplot2::aes(colour = factor(original_cluster, levels = original_levels)),
      size = 0.55,
      alpha = 0.9
    ) +
    ggplot2::labs(title = paste0(sample_id, " marker-filter removed cells on UMAP"),
                  colour = "Original") +
    ggplot2::theme_classic()
  plot_and_save(p_removed, file.path(dir_removed, "UMAP_removed_cells_by_original_cluster.pdf"),
                width = 8, height = 7)

  marker_panel_original <- original_markers_pass %>%
    dplyr::filter(gene %in% rownames(obj[["RNA"]])) %>%
    dplyr::mutate(cluster = factor(as.character(cluster), levels = original_levels)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(marker_log2FC), .by_group = TRUE) %>%
    dplyr::mutate(marker_rank = dplyr::row_number()) %>%
    dplyr::slice_head(n = config$marker_filter$heatmap_top_n) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, marker_rank)
  write_tsv(marker_panel_original,
            file.path(dir_before, "heatmap_marker_panel_original_cluster_top10.tsv"))
  marker_genes_original <- marker_panel_original %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(cluster, marker_rank) %>%
    dplyr::pull(gene)
  if (length(marker_genes_original) > 0) {
    DefaultAssay(obj) <- "RNA"
    obj <- ScaleData(obj, assay = "RNA", features = marker_genes_original, verbose = FALSE)
    obj$original_cluster_for_heatmap <- factor(obj$original_cluster, levels = original_levels)
    plot_heatmap_original <- function(cells, title) {
      DoHeatmap(
        object = obj,
        features = marker_genes_original,
        cells = cells,
        group.by = "original_cluster_for_heatmap",
        assay = "RNA",
        slot = "scale.data",
        disp.min = -2.5,
        disp.max = 2.5,
        raster = TRUE,
        size = 3
      ) + ggplot2::ggtitle(title) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
    }
    heatmap_height <- max(8, length(marker_genes_original) * 0.12 + 2)
    p_heatmap_before <- plot_heatmap_original(colnames(obj),
                                              paste0(sample_id, " before marker filter"))
    p_heatmap_after <- plot_heatmap_original(kept_cells,
                                             paste0(sample_id, " after marker filter"))
    plot_and_save(p_heatmap_before,
                  file.path(dir_before, "RNA_DoHeatmap_original_cluster_before_marker_filter.pdf"),
                  width = 12, height = heatmap_height, limitsize = FALSE)
    plot_and_save(p_heatmap_before + p_heatmap_after,
                  file.path(dir_removed, "RNA_DoHeatmap_before_after_original_cluster.pdf"),
                  width = 18, height = heatmap_height, limitsize = FALSE)
  }

  filtered <- subset(obj, cells = kept_cells)
  write_tsv(filtered@meta.data %>% tibble::rownames_to_column("cell"),
            file.path(dir_after, "cell_metadata_after_marker_filter_before_top_mito_doublet.tsv"))
  write_tsv(data.frame(
    sample_id = sample_id,
    n_cells_after_soupx = ncol(soupx_obj),
    n_cells_removed_by_marker_first = length(removed_cells),
    n_cells_after_marker_filter = ncol(filtered),
    n_flagged_leiden_clusters = length(flagged_leiden_ids),
    flagged_leiden_clusters = paste(flagged_leiden_ids, collapse = ","),
    marker_filter_cutoff = marker_cutoff
  ), file.path(marker_first_dir, "marker_first_filter_summary.tsv"))

  list(obj = filtered, status = cell_removal_status, marker_first_dir = marker_first_dir)
}

apply_top_mito_after_marker <- function(sample_row, marker_obj, sample_dir, config) {
  sample_id <- sample_row$sample_id[[1]]
  qc_dir <- qc_step_dir(sample_dir, "top_mito")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  object_dir <- if (qc_use_numbered_dirs()) qc_step_dir(sample_dir, "objects") else qc_dir
  dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(object_dir, "object_after_marker_top_mito_before_DoubletFinder.rds")
  status_file <- file.path(qc_dir, "top_mito_filter_cell_status.tsv")

  if (isTRUE(config$output$resume_from_existing) &&
      !isTRUE(config$output$overwrite_existing) &&
      file.exists(cache_file) && file.exists(status_file)) {
    log_msg(sample_id, ": reusing object after marker + top/mito filter.")
    return(list(obj = readRDS(cache_file), status = read_tsv(status_file), qc_dir = qc_dir))
  }

  log_msg(sample_id, ": applying top UMI and mito filter after marker filtering.")
  obj <- marker_obj
  if (!"percent.mito" %in% colnames(obj@meta.data)) obj$percent.mito <- 0

  apply_mito <- to_bool(sample_value(sample_row, "apply_mito_filter", "FALSE"))
  nUMI_high <- as.numeric(quantile(obj$nCount_RNA, config$top_mito$top_umi_quantile, na.rm = TRUE))
  remove_top <- colnames(obj)[obj$nCount_RNA >= nUMI_high]
  remove_mito <- if (apply_mito) colnames(obj)[obj$percent.mito >= config$top_mito$mito_cutoff] else character(0)
  remove_any <- union(remove_top, remove_mito)
  kept <- setdiff(colnames(obj), remove_any)

  status <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(
      remove_topUMI_after_marker = cell %in% remove_top,
      remove_mito_after_marker = cell %in% remove_mito,
      remove_top_mito_after_marker = cell %in% remove_any,
      topUMI_threshold_after_marker = nUMI_high,
      mito_cutoff = config$top_mito$mito_cutoff,
      apply_mito_filter = apply_mito
    )
  write_tsv(status, status_file)
  writeLines(remove_top, file.path(qc_dir, "removed_topUMI_barcodes_after_marker.txt"))
  writeLines(remove_mito, file.path(qc_dir, "removed_mito_barcodes_after_marker.txt"))
  writeLines(kept, file.path(qc_dir, "after_marker_top_mito_cell_barcodes.txt"))

  filtered <- subset(obj, cells = kept)
  if (isTRUE(config$output$save_before_doubletfinder_object)) saveRDS(filtered, cache_file)
  write_tsv(data.frame(
    sample_id = sample_id,
    n_cells_after_marker_filter = ncol(obj),
    top_umi_quantile = config$top_mito$top_umi_quantile,
    topUMI_threshold_after_marker = nUMI_high,
    apply_mito_filter = apply_mito,
    mito_cutoff = config$top_mito$mito_cutoff,
    n_removed_topUMI_after_marker = length(remove_top),
    n_removed_mito_after_marker = length(remove_mito),
    n_removed_union_top_mito_after_marker = length(remove_any),
    n_cells_after_marker_top_mito_filter = ncol(filtered),
    object_after_marker_top_mito_before_DoubletFinder = cache_file
  ), file.path(qc_dir, "top_mito_filter_summary.tsv"))

  list(obj = filtered, status = status, qc_dir = qc_dir)
}
