prepare_doubletfinder_base <- function(obj, sample_row, sample_dir, config) {
  sample_id <- sample_row$sample_id[[1]]
  common_dir <- qc_step_dir(sample_dir, "doubletfinder_common")
  dir.create(common_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg(sample_id, ": preparing RNA workflow before DoubletFinder.")
  DefaultAssay(obj) <- "RNA"
  obj <- join_assay_layers_if_needed(obj, "RNA")
  obj@reductions <- list()
  obj@graphs <- list()
  obj@neighbors <- list()
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = config$final_recluster$nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  initial_dims <- safe_dims(obj, config$doubletfinder$dims_max)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = max(initial_dims), verbose = FALSE)
  df_pcs <- initial_dims[initial_dims <= min(config$doubletfinder$neighbor_dims_max, max(initial_dims))]
  obj <- FindNeighbors(obj, dims = df_pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = config$doubletfinder$resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = initial_dims, verbose = FALSE)

  log_msg(sample_id, ": running DoubletFinder paramSweep.")
  sweep.res.list <- DoubletFinder::paramSweep(obj, PCs = df_pcs, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  write_tsv(bcmvn, file.path(common_dir, "DoubletFinder_find_pK_table.tsv"))
  pK_value <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric, na.rm = TRUE)][1]))
  homotypic.prop <- DoubletFinder::modelHomotypic(obj@meta.data$seurat_clusters)

  umap <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
  umap$cell <- rownames(umap)
  umap <- dplyr::left_join(umap, obj@meta.data %>% tibble::rownames_to_column("cell"),
                           by = "cell")
  write_tsv(umap, file.path(common_dir, "RNA_UMAP_coordinates_before_DoubletFinder.tsv"))
  p_umap <- ggplot2::ggplot(umap, ggplot2::aes(UMAP_1, UMAP_2, colour = seurat_clusters)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.75) +
    ggplot2::labs(title = paste0(sample_id, " RNA UMAP before DoubletFinder"),
                  colour = "Cluster") +
    ggplot2::theme_classic()
  plot_and_save(p_umap,
                file.path(common_dir, "RNA_UMAP_before_DoubletFinder_by_cluster.pdf"),
                width = 8, height = 7)

  write_tsv(data.frame(
    sample_id = sample_id,
    pK = pK_value,
    pK_source = "recomputed_current_marker_topUMI_mito_filtered_object",
    homotypic_prop = homotypic.prop,
    df_pcs = paste(df_pcs, collapse = ","),
    n_cells_before_doubletfinder = ncol(obj)
  ), file.path(common_dir, "DoubletFinder_common_pK_summary.tsv"))
  list(obj = obj, pK_value = pK_value, homotypic_prop = homotypic.prop,
       df_pcs = df_pcs, common_dir = common_dir)
}

write_final_cluster_umi_qc <- function(obj, sample_id, rate_label, dir_final, cluster_levels) {
  cluster_umi <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(
      Cluster = factor(as.character(new_cluster), levels = cluster_levels),
      UMI = nCount_RNA
    ) %>%
    dplyr::select(cell, Cluster, UMI)

  write_tsv(cluster_umi, file.path(dir_final, "final_cluster_UMI_per_cell.tsv"))

  cluster_stats <- cluster_umi %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      Mean_UMI = mean(UMI, na.rm = TRUE),
      Median_UMI = median(UMI, na.rm = TRUE),
      Min_UMI = min(UMI, na.rm = TRUE),
      Max_UMI = max(UMI, na.rm = TRUE),
      SD_UMI = stats::sd(UMI, na.rm = TRUE),
      Cell_Count = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Cluster)

  write_tsv(cluster_stats, file.path(dir_final, "final_cluster_UMI_summary.tsv"))

  p_umi <- ggplot2::ggplot(
    cluster_umi,
    ggplot2::aes(x = Cluster, y = UMI, fill = Cluster)
  ) +
    ggplot2::geom_violin(alpha = 0.7, trim = TRUE, scale = "width") +
    ggplot2::geom_jitter(size = 0.1, alpha = 0.2, width = 0.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste0(sample_id, " ", rate_label, " UMIplot in each final cluster"),
      x = "Cluster",
      y = "UMI counts"
    ) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::scale_y_log10()

  plot_and_save(
    p_umi,
    file.path(dir_final, "final_cluster_UMI_violin_log10.pdf"),
    width = max(8, length(cluster_levels) * 0.45),
    height = 6,
    limitsize = FALSE
  )

  invisible(cluster_stats)
}

run_doublet_and_final_recluster <- function(sample_row, soupx_before, marker_status,
                                            top_mito_status, after_top_mito_obj,
                                            df_common, sample_dir, config) {
  sample_id <- sample_row$sample_id[[1]]
  rate <- config$doubletfinder$doublet_rate
  rate_label <- sprintf("doublet_rate_%.2f", rate)
  rate_dir <- qc_rate_dir(sample_dir, rate_label)
  dir_df <- file.path(rate_dir, "DoubletFinder_QC")
  dir_final <- qc_final_dir(sample_dir, rate_label)
  dir_summary <- file.path(rate_dir, "summary")
  for (d in c(dir_df, dir_final, dir_summary)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  object_dir <- if (qc_use_numbered_dirs()) qc_step_dir(sample_dir, "objects") else dir_final
  dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
  final_rds <- file.path(object_dir, "final_filtered_reclustered_object.rds")
  summary_file <- file.path(dir_summary, "doublet_rate_QC_summary.tsv")
  if (isTRUE(config$output$resume_from_existing) &&
      !isTRUE(config$output$overwrite_existing) &&
      file.exists(final_rds) && file.exists(summary_file)) {
    log_msg(sample_id, " ", rate_label, ": reusing existing final object.")
    umi_files <- file.path(
      dir_final,
      c("final_cluster_UMI_per_cell.tsv",
        "final_cluster_UMI_summary.tsv",
        "final_cluster_UMI_violin_log10.pdf")
    )
    if (!all(file.exists(umi_files))) {
      log_msg(sample_id, " ", rate_label, ": writing missing final cluster UMI QC from existing final object.")
      existing_final <- readRDS(final_rds)
      if (!"new_cluster" %in% colnames(existing_final@meta.data)) {
        new_cluster_col <- paste0("RNA_snn_res.", config$final_recluster$resolution)
        if (!new_cluster_col %in% colnames(existing_final@meta.data)) new_cluster_col <- "seurat_clusters"
        existing_final$new_cluster <- as.character(existing_final@meta.data[[new_cluster_col]])
      }
      write_final_cluster_umi_qc(
        existing_final,
        sample_id,
        rate_label,
        dir_final,
        numeric_order(existing_final$new_cluster)
      )
    }
    return(read_tsv(summary_file))
  }

  log_msg(sample_id, " ", rate_label, ": running DoubletFinder.")
  obj <- df_common$obj
  pK_value <- df_common$pK_value
  df_pcs <- df_common$df_pcs
  homotypic.prop <- df_common$homotypic_prop
  nExp_poi <- round(rate * length(obj@meta.data$seurat_clusters))
  nExp_poi.adj <- if (isTRUE(config$doubletfinder$use_homotypic_adjustment)) {
    round(nExp_poi * (1 - homotypic.prop))
  } else {
    nExp_poi
  }
  pN_value <- config$doubletfinder$pN
  df_fun <- if (exists("doubletFinder", envir = asNamespace("DoubletFinder"), inherits = FALSE)) {
    get("doubletFinder", envir = asNamespace("DoubletFinder"))
  } else {
    get("doubletFinder_v3", envir = asNamespace("DoubletFinder"))
  }
  obj <- df_fun(obj, PCs = df_pcs, pN = pN_value, pK = pK_value, nExp = nExp_poi)
  pann_cols <- grep("^pANN", colnames(obj[[]]), value = TRUE)
  stopifnot(length(pann_cols) > 0)
  obj <- df_fun(obj, PCs = df_pcs, pN = pN_value, pK = pK_value,
                nExp = nExp_poi.adj, reuse.pANN = tail(pann_cols, 1))
  df_col <- tail(grep("^DF.classifications", colnames(obj@meta.data), value = TRUE), 1)
  pann_col <- tail(grep("^pANN", colnames(obj@meta.data), value = TRUE), 1)
  obj$Doublet <- obj@meta.data[[df_col]]
  df_meta <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(DoubletFinder_final_col = df_col, pANN_final_col = pann_col)
  write_tsv(df_meta, file.path(dir_df, "DoubletFinder_metadata.tsv"))
  pANN_cutoff <- suppressWarnings(min(df_meta[[pann_col]][df_meta[[df_col]] == "Doublet"], na.rm = TRUE))
  if (!is.finite(pANN_cutoff)) pANN_cutoff <- NA_real_
  write_tsv(data.frame(
    sample_id = sample_id,
    doublet_rate = rate,
    pN = pN_value,
    pK = pK_value,
    nExp_poi = nExp_poi,
    homotypic_prop = homotypic.prop,
    use_homotypic_adjustment = config$doubletfinder$use_homotypic_adjustment,
    nExp_poi_adjusted = nExp_poi.adj,
    pANN_cutoff = pANN_cutoff,
    singlets = sum(obj$Doublet == "Singlet"),
    doublets = sum(obj$Doublet != "Singlet")
  ), file.path(dir_df, "DoubletFinder_summary.tsv"))

  df_umap <- as.data.frame(Embeddings(obj, "umap"))
  colnames(df_umap)[1:2] <- c("UMAP_1", "UMAP_2")
  df_umap$cell <- rownames(df_umap)
  df_umap <- dplyr::left_join(df_umap, df_meta, by = "cell")
  write_tsv(df_umap, file.path(dir_df, "RNA_UMAP_coordinates_with_DoubletFinder.tsv"))
  p_df <- ggplot2::ggplot(df_umap, ggplot2::aes(UMAP_1, UMAP_2, colour = Doublet)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.75) +
    ggplot2::labs(title = paste0(sample_id, " ", rate_label, " DoubletFinder class")) +
    ggplot2::theme_classic()
  plot_and_save(p_df, file.path(dir_df, "RNA_UMAP_by_DoubletFinder_class.pdf"),
                width = 8, height = 7)

  filtered <- subset(obj, subset = Doublet == "Singlet")
  writeLines(colnames(filtered), file.path(dir_df, "after_DoubletFinder_singlet_barcodes.txt"))

  log_msg(sample_id, " ", rate_label, ": final RNA reclustering.")
  DefaultAssay(filtered) <- "RNA"
  filtered <- join_assay_layers_if_needed(filtered, "RNA")
  filtered@reductions <- list()
  filtered@graphs <- list()
  filtered@neighbors <- list()
  filtered <- NormalizeData(filtered, assay = "RNA", normalization.method = "LogNormalize",
                            scale.factor = 10000, verbose = FALSE)
  filtered <- FindVariableFeatures(filtered, assay = "RNA", selection.method = "vst",
                                   nfeatures = config$final_recluster$nfeatures, verbose = FALSE)
  filtered <- ScaleData(filtered, assay = "RNA", features = VariableFeatures(filtered),
                        verbose = FALSE)
  post_dims <- safe_dims(filtered, config$final_recluster$npcs_max)
  post_use_dims <- post_dims[post_dims <= min(config$final_recluster$dims_max, max(post_dims))]
  filtered <- RunPCA(filtered, assay = "RNA", features = VariableFeatures(filtered),
                     npcs = max(post_dims), verbose = FALSE)
  filtered <- FindNeighbors(filtered, reduction = "pca", dims = post_use_dims, verbose = FALSE)
  filtered <- FindClusters(filtered, resolution = config$final_recluster$resolution, verbose = FALSE)
  filtered <- RunUMAP(filtered, reduction = "pca", dims = post_use_dims, verbose = FALSE)
  new_cluster_col <- paste0("RNA_snn_res.", config$final_recluster$resolution)
  if (!new_cluster_col %in% colnames(filtered@meta.data)) new_cluster_col <- "seurat_clusters"
  filtered$new_cluster <- as.character(filtered@meta.data[[new_cluster_col]])
  Idents(filtered) <- "new_cluster"
  new_levels <- numeric_order(filtered$new_cluster)

  filtered_meta <- filtered@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(new_cluster = as.character(new_cluster))
  write_tsv(filtered_meta, file.path(dir_final, "reclustered_cell_metadata.tsv"))

  filtered_umap <- as.data.frame(Embeddings(filtered, "umap"))
  colnames(filtered_umap)[1:2] <- c("UMAP_1", "UMAP_2")
  filtered_umap$cell <- rownames(filtered_umap)
  filtered_umap <- dplyr::left_join(
    filtered_umap,
    filtered_meta[, intersect(c("cell", "new_cluster", "original_cluster", "leiden",
                                "nCount_RNA", "nFeature_RNA", "percent.mito", "Doublet"),
                              colnames(filtered_meta))],
    by = "cell"
  )
  write_tsv(filtered_umap, file.path(dir_final, "reclustered_umap_coordinates.tsv"))

  write_final_cluster_umi_qc(filtered, sample_id, rate_label, dir_final, new_levels)

  final_markers_raw <- FindAllMarkers(
    object = filtered,
    assay = "RNA",
    only.pos = TRUE,
    test.use = "wilcox",
    logfc.threshold = config$final_recluster$marker_logfc_cutoff,
    min.pct = config$final_recluster$marker_min_pct,
    return.thresh = 1,
    verbose = TRUE
  )
  final_markers_raw$cluster <- as.character(final_markers_raw$cluster)
  if (!"gene" %in% colnames(final_markers_raw)) final_markers_raw$gene <- rownames(final_markers_raw)
  final_fc_col <- get_fc_col(final_markers_raw)
  final_markers_pass <- final_markers_raw %>%
    dplyr::filter(
      p_val_adj < config$final_recluster$marker_padj_cutoff,
      .data[[final_fc_col]] > config$final_recluster$marker_logfc_cutoff,
      pct.1 >= config$final_recluster$marker_min_pct
    ) %>%
    dplyr::mutate(delta_pct = pct.1 - pct.2)
  if (isTRUE(config$output$save_raw_marker_gz)) {
    write_tsv_gz(final_markers_raw, file.path(dir_final, "new_cluster_markers_raw.tsv.gz"))
  }
  write_tsv(final_markers_pass, file.path(dir_final, "new_cluster_markers_pass.tsv"))

  new_summary <- filtered_meta %>%
    dplyr::count(new_cluster, name = "n_cells") %>%
    dplyr::left_join(final_markers_pass %>% dplyr::count(cluster, name = "n_markers_pass"),
                     by = c("new_cluster" = "cluster")) %>%
    dplyr::mutate(n_markers_pass = ifelse(is.na(n_markers_pass), 0L, n_markers_pass)) %>%
    dplyr::arrange(factor(new_cluster, levels = new_levels))
  write_tsv(new_summary, file.path(dir_final, "new_cluster_summary.tsv"))

  p_final <- ggplot2::ggplot(
    filtered_umap,
    ggplot2::aes(UMAP_1, UMAP_2, colour = factor(new_cluster, levels = new_levels))
  ) +
    ggplot2::geom_point(size = 0.25, alpha = 0.75) +
    ggplot2::labs(title = paste0(sample_id, " ", rate_label, " final UMAP"),
                  colour = "New cluster") +
    ggplot2::theme_classic()
  plot_and_save(p_final, file.path(dir_final, "reclustered_UMAP_by_new_cluster.pdf"),
                width = 8, height = 7)

  top_markers_final <- final_markers_pass %>%
    dplyr::filter(gene %in% rownames(filtered[["RNA"]])) %>%
    dplyr::mutate(marker_source_cluster = factor(as.character(cluster), levels = new_levels)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(p_val_adj, dplyr::desc(.data[[final_fc_col]]), .by_group = TRUE) %>%
    dplyr::mutate(marker_rank = dplyr::row_number()) %>%
    dplyr::slice_head(n = config$final_recluster$heatmap_top_n) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(marker_source_cluster, marker_rank)
  write_tsv(top_markers_final, file.path(dir_final, "new_cluster_top10_markers_for_heatmap.tsv"))
  top_genes_final <- top_markers_final %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(marker_source_cluster, marker_rank) %>%
    dplyr::pull(gene)
  if (length(top_genes_final) > 0) {
    filtered <- ScaleData(filtered, assay = "RNA", features = top_genes_final, verbose = FALSE)
    filtered$new_cluster_for_heatmap <- factor(filtered$new_cluster, levels = new_levels)
    p_heatmap_final <- DoHeatmap(
      object = filtered,
      features = top_genes_final,
      group.by = "new_cluster_for_heatmap",
      assay = "RNA",
      slot = "scale.data",
      disp.min = -2.5,
      disp.max = 2.5,
      raster = TRUE,
      size = 3
    ) +
      ggplot2::ggtitle(paste0(sample_id, " ", rate_label, " final RNA marker heatmap")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
    plot_and_save(p_heatmap_final,
                  file.path(dir_final, "new_cluster_RNA_DoHeatmap_top10_markers.pdf"),
                  width = 12, height = max(8, length(top_genes_final) * 0.12 + 2),
                  limitsize = FALSE)
  }

  if (isTRUE(config$output$save_final_object)) saveRDS(filtered, final_rds)
  rate_summary <- data.frame(
    sample_id = sample_id,
    doublet_rate = rate,
    rate_label = rate_label,
    qc_order = "emptyDrops_SoupX_markerFilter_top1UMI_mito_DoubletFinder",
    n_cells_after_soupx = ncol(soupx_before),
    n_cells_removed_by_marker_first = sum(marker_status$remove_by_marker),
    n_cells_after_marker_filter = sum(!marker_status$remove_by_marker),
    n_cells_removed_topUMI_after_marker = sum(top_mito_status$remove_topUMI_after_marker),
    n_cells_removed_mito_after_marker = sum(top_mito_status$remove_mito_after_marker),
    n_cells_before_doubletfinder = ncol(after_top_mito_obj),
    pK = pK_value,
    homotypic_prop = homotypic.prop,
    nExp_poi = nExp_poi,
    nExp_poi_adjusted = nExp_poi.adj,
    pANN_cutoff = pANN_cutoff,
    n_doublets = sum(df_meta$Doublet != "Singlet"),
    n_singlets = sum(df_meta$Doublet == "Singlet"),
    n_final_clusters = length(new_levels),
    final_filtered_rds = final_rds
  )
  write_tsv(rate_summary, summary_file)
  capture.output(sessionInfo(), file = file.path(dir_summary, "sessionInfo.txt"))
  rate_summary
}
