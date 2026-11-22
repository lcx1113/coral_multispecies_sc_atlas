draw_soupx_autoEstCont_fit <- function(sc) {
  if (is.null(sc$fit) || is.null(sc$fit$posterior) || is.null(sc$fit$rhoEst)) {
    stop("SoupX fit is missing; run autoEstCont before plotting.")
  }
  rho_probes <- seq(0, 1, 0.001)
  posterior <- sc$fit$posterior
  prior_rho <- sc$fit$priorRho
  prior_rho_sd <- sc$fit$priorRhoStdDev
  v2 <- (prior_rho_sd / prior_rho)^2
  k <- 1 + v2^-2 / 2 * (1 + sqrt(1 + 4 * v2))
  theta <- prior_rho / (k - 1)
  prior_density <- dgamma(rho_probes, k, scale = theta)
  rho_fwhm <- sc$fit$rhoFWHM
  plot(
    rho_probes, posterior, "l",
    xlim = c(0, 1),
    ylim = c(0, max(c(prior_density, posterior), na.rm = TRUE)),
    frame.plot = FALSE,
    xlab = "Contamination Fraction",
    ylab = "Probability Density"
  )
  lines(rho_probes, prior_density, lty = 2)
  abline(v = sc$fit$rhoEst, col = "red")
  legend(
    x = "topright",
    legend = c(
      sprintf("prior rho %g(+/-%g)", prior_rho, prior_rho_sd),
      sprintf("post rho %g(%g,%g)", sc$fit$rhoEst, rho_fwhm[1], rho_fwhm[2]),
      "rho max"
    ),
    lty = c(2, 1, 1),
    col = c("black", "black", "red"),
    bty = "n"
  )
}

run_autoEstCont_and_save_plot <- function(sc, soupx_dir, sample_id, config) {
  force_accept <- config$soupx$force_accept
  auto_estimate <- function(sc, do_plot = TRUE) {
    call_args <- list(sc = sc, doPlot = do_plot)
    # "default" preserves the historical call that omitted forceAccept.
    if (!identical(force_accept, "default")) call_args$forceAccept <- force_accept
    do.call(SoupX::autoEstCont, call_args)
  }
  if (!isTRUE(config$soupx$save_autoEstCont_plot)) {
    return(auto_estimate(sc))
  }
  contamination_pdf <- file.path(soupx_dir, "SoupX_autoEstCont_contamination_fraction.pdf")
  contamination_png <- file.path(soupx_dir, "SoupX_autoEstCont_contamination_fraction.png")

  grDevices::pdf(contamination_pdf, width = 7, height = 5)
  sc <- tryCatch(
    auto_estimate(sc, do_plot = TRUE),
    finally = {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    }
  )

  grDevices::png(contamination_png, width = 7, height = 5, units = "in", res = 300)
  tryCatch(
    draw_soupx_autoEstCont_fit(sc),
    finally = {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    }
  )

  write_tsv(data.frame(
    sample_id = sample_id,
    soupX_rhoEst = sc$fit$rhoEst,
    soupX_rhoFWHM_low = sc$fit$rhoFWHM[1],
    soupX_rhoFWHM_high = sc$fit$rhoFWHM[2],
    forceAccept = force_accept,
    contamination_fraction_pdf = contamination_pdf,
    contamination_fraction_png = contamination_png
  ), file.path(soupx_dir, "SoupX_autoEstCont_plot_summary.tsv"))
  sc
}

run_soupx_step <- function(sample_row, count_matrix, emptydrops_cells, sample_dir, config) {
  sample_id <- sample_row$sample_id[[1]]
  soupx_dir <- qc_step_dir(sample_dir, "soupx")
  dir.create(soupx_dir, recursive = TRUE, showWarnings = FALSE)
  object_dir <- if (qc_use_numbered_dirs()) qc_step_dir(sample_dir, "objects") else soupx_dir
  dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
  target_obj_file <- file.path(object_dir, "soupx_object_after_correction_before_marker_top_mito_doublet.rds")
  summary_file <- file.path(soupx_dir, "soupx_cache_summary.tsv")

  if (isTRUE(config$output$resume_from_existing) &&
      !isTRUE(config$output$overwrite_existing) &&
      file.exists(target_obj_file) && file.exists(summary_file)) {
    log_msg(sample_id, ": reusing existing SoupX object cache.")
    return(readRDS(target_obj_file))
  }

  cell_barcodes <- intersect(emptydrops_cells, colnames(count_matrix))
  if (length(cell_barcodes) < 100) stop(sample_id, ": too few emptyDrops cells for SoupX.")
  writeLines(cell_barcodes, file.path(soupx_dir, "emptyDrops_cell_barcodes_used_for_SoupX.txt"))

  log_msg(sample_id, ": running SoupX on emptyDrops-retained cells.")
  toc <- count_matrix[, cell_barcodes, drop = FALSE]
  tod <- count_matrix

  coral_temp <- CreateSeuratObject(counts = toc, min.cells = config$soupx$min_cells)
  coral_temp <- NormalizeData(coral_temp, verbose = FALSE)
  coral_temp <- FindVariableFeatures(coral_temp, nfeatures = config$soupx$temp_nfeatures, verbose = FALSE)
  coral_temp <- ScaleData(coral_temp, features = VariableFeatures(coral_temp), verbose = FALSE)
  temp_dims <- safe_dims(coral_temp, config$soupx$temp_npcs)
  coral_temp <- RunPCA(coral_temp, npcs = max(temp_dims), verbose = FALSE)
  coral_temp <- RunUMAP(coral_temp, dims = temp_dims, verbose = FALSE)
  neighbor_dims <- temp_dims[temp_dims <= min(config$soupx$temp_neighbor_dims_max, max(temp_dims))]
  coral_temp <- FindNeighbors(coral_temp, dims = neighbor_dims, verbose = FALSE)
  coral_temp <- FindClusters(coral_temp, resolution = config$soupx$cluster_resolution, verbose = FALSE)

  sc <- SoupX::SoupChannel(tod, toc)
  sc <- SoupX::setClusters(sc, coral_temp$seurat_clusters)
  sc <- SoupX::setDR(sc, Embeddings(coral_temp, "umap"))
  sc <- run_autoEstCont_and_save_plot(sc, soupx_dir, sample_id, config)
  write_tsv(data.frame(
    soupX_rhoEst = sc$fit$rhoEst,
    forceAccept = config$soupx$force_accept
  ),
            file.path(soupx_dir, "SoupX_contamination_summary.tsv"))
  write_tsv(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], 50),
            file.path(soupx_dir, "SoupX_top50_soup_profile.tsv"))

  log_msg(sample_id, ": adjusting SoupX counts with roundToInt=", config$soupx$round_to_int)
  corrected_counts <- SoupX::adjustCounts(sc, roundToInt = config$soupx$round_to_int)
  coral_corrected <- CreateSeuratObject(
    counts = corrected_counts,
    project = sample_id,
    min.cells = config$soupx$min_cells
  )
  DefaultAssay(coral_corrected) <- "RNA"

  mito_genes <- read_mito_genes(sample_row)
  mito_gene_in_data <- intersect(mito_genes, Features(coral_corrected[["RNA"]]))
  if (length(mito_genes) > 0) {
    coral_corrected <- PercentageFeatureSet(
      object = coral_corrected,
      features = mito_gene_in_data,
      col.name = "percent.mito",
      assay = "RNA"
    )
  } else {
    coral_corrected$percent.mito <- 0
  }
  writeLines(mito_gene_in_data, file.path(soupx_dir, "mitochondrial_genes_found_in_data.txt"))

  meta <- coral_corrected@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(filter_stage = "after_SoupX_before_marker_top_mito_doublet")
  write_tsv(meta, file.path(soupx_dir, "cell_metadata_after_SoupX_before_marker_top_mito_doublet.tsv"))

  p_vln <- VlnPlot(
    coral_corrected,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
    assay = "RNA",
    layer = "counts",
    combine = TRUE,
    ncol = 3,
    pt.size = 0.05
  )
  plot_and_save(p_vln, file.path(soupx_dir, "QC_violin_after_SoupX_before_marker_filter.pdf"),
                width = 10, height = 5)
  p_scatter1 <- FeatureScatter(coral_corrected, feature1 = "nCount_RNA", feature2 = "percent.mito")
  p_scatter2 <- FeatureScatter(coral_corrected, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot_and_save(p_scatter1 + p_scatter2,
                file.path(soupx_dir, "QC_feature_scatter_after_SoupX_before_marker_filter.pdf"),
                width = 10, height = 5)

  if (isTRUE(config$output$save_soupx_object)) saveRDS(coral_corrected, target_obj_file)
  write_tsv(data.frame(
    sample_id = sample_id,
    n_emptyDrops_cells_used = length(cell_barcodes),
    n_cells_after_soupx = ncol(coral_corrected),
    n_genes_after_soupx = nrow(coral_corrected),
    apply_mito_filter = to_bool(sample_value(sample_row, "apply_mito_filter", "FALSE")),
    mito_gene_file = sample_value(sample_row, "mito_gene_file", ""),
    n_mito_genes_in_data = length(mito_gene_in_data),
    soupX_rhoEst = sc$fit$rhoEst,
    soupx_object_after_correction_before_marker_top_mito_doublet = target_obj_file
  ), summary_file)
  coral_corrected
}
