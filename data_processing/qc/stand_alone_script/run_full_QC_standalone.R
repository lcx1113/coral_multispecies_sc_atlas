#!/usr/bin/env Rscript
# Post-analysis reorganized source, not the historical QC execution script.
# Layout/control restored 2026-11-16; confirmed parameter edits are in default_config.R.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
this_script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
standalone_script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
framework_root <- normalizePath(file.path(standalone_script_dir, ".."), mustWork = FALSE)

source(file.path(framework_root, "R", "io_utils.R"))
source(file.path(framework_root, "R", "marker_utils.R"))
source(file.path(framework_root, "R", "emptydrops_utils.R"))
source(file.path(framework_root, "R", "soupx_utils.R"))
source(file.path(framework_root, "R", "marker_filter_utils.R"))
source(file.path(framework_root, "R", "doubletfinder_final_utils.R"))
source(file.path(framework_root, "R", "summary_utils.R"))

options(qc_numbered_dirs = TRUE)

normalise_pipeline_step <- function(x, default = "none") {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  x <- tolower(gsub("[.-]", "_", as.character(x)))
  aliases <- c(
    "none" = "none",
    "no" = "none",
    "false" = "none",
    "full" = "emptydrops",
    "all" = "emptydrops",
    "emptydrop" = "emptydrops",
    "emptydrops" = "emptydrops",
    "emptydrops_qc" = "emptydrops",
    "soupx" = "soupx",
    "soup" = "soupx",
    "marker" = "marker_filter",
    "marker_filter" = "marker_filter",
    "marker_first" = "marker_filter",
    "top_mito" = "top_mito",
    "topumi_mito" = "top_mito",
    "topumi" = "top_mito",
    "doubletfinder" = "doubletfinder",
    "doublet" = "doubletfinder",
    "df" = "doubletfinder",
    "final" = "final_recluster",
    "final_recluster" = "final_recluster",
    "recluster" = "final_recluster"
  )
  if (!x %in% names(aliases)) {
    stop(
      "Unsupported step: ", x,
      ". Allowed values: none, emptydrops, soupx, marker_filter, top_mito, doubletfinder, final_recluster."
    )
  }
  unname(aliases[[x]])
}

step_rank <- function(step) {
  ranks <- c(
    emptydrops = 1L,
    soupx = 2L,
    marker_filter = 3L,
    top_mito = 4L,
    doubletfinder = 5L,
    final_recluster = 6L
  )
  ranks[[step]]
}

rate_label_from_config <- function(config) {
  sprintf("doublet_rate_%.2f", config$doubletfinder$doublet_rate)
}

delete_if_exists <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths) > 0) unlink(paths, recursive = TRUE, force = TRUE)
}

delete_outputs_from_step <- function(sample_dir, config, rerun_from) {
  if (rerun_from == "none") return(invisible(NULL))
  rank <- step_rank(rerun_from)
  rate_label <- rate_label_from_config(config)
  rate_dir <- qc_rate_dir(sample_dir, rate_label)
  dir_df <- file.path(rate_dir, "DoubletFinder_QC")
  dir_final <- qc_final_dir(sample_dir, rate_label)
  dir_rate_summary <- file.path(rate_dir, "summary")

  remove_paths <- character(0)
  if (rank <= step_rank("emptydrops")) {
    remove_paths <- c(remove_paths, qc_step_dir(sample_dir, "emptydrops"))
  }
  if (rank <= step_rank("soupx")) {
    remove_paths <- c(remove_paths, qc_step_dir(sample_dir, "soupx"),
                      file.path(qc_step_dir(sample_dir, "objects"), "soupx_object_after_correction_before_marker_top_mito_doublet.rds"))
  }
  if (rank <= step_rank("marker_filter")) {
    remove_paths <- c(remove_paths, qc_step_dir(sample_dir, "marker_filter"))
  }
  if (rank <= step_rank("top_mito")) {
    remove_paths <- c(remove_paths, qc_step_dir(sample_dir, "top_mito"),
                      file.path(qc_step_dir(sample_dir, "objects"), "object_after_marker_top_mito_before_DoubletFinder.rds"))
  }
  if (rank <= step_rank("doubletfinder")) {
    remove_paths <- c(remove_paths, qc_step_dir(sample_dir, "doubletfinder_common"), rate_dir,
                      qc_final_dir(sample_dir, rate_label),
                      file.path(qc_step_dir(sample_dir, "objects"), "final_filtered_reclustered_object.rds"),
                      qc_step_dir(sample_dir, "report"))
  } else if (rank <= step_rank("final_recluster")) {
    remove_paths <- c(
      remove_paths,
      dir_final,
      file.path(dir_rate_summary, "doublet_rate_QC_summary.tsv"),
      file.path(dir_rate_summary, "sessionInfo.txt"),
      file.path(qc_step_dir(sample_dir, "objects"), "final_filtered_reclustered_object.rds"),
      qc_step_dir(sample_dir, "report")
    )
  }
  delete_if_exists(unique(remove_paths))
  invisible(NULL)
}

require_existing_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing required ", label, ": ", path)
  path
}

load_soupx_object_cache <- function(sample_dir) {
  readRDS(require_existing_file(
    file.path(qc_step_dir(sample_dir, "objects"), "soupx_object_after_correction_before_marker_top_mito_doublet.rds"),
    "SoupX-corrected object"
  ))
}

load_marker_filter_state <- function(sample_row, soupx_obj, sample_dir) {
  sample_id <- sample_row$sample_id[[1]]
  dir_marker <- file.path(qc_step_dir(sample_dir, "marker_filter"), "03b_marker_tables")
  status_file <- require_existing_file(file.path(dir_marker, "cell_removal_status.tsv"), "marker filter status")
  kept_file <- require_existing_file(file.path(dir_marker, "kept_after_marker_filter_barcodes.txt"), "marker-kept barcodes")
  kept <- intersect(readLines(kept_file), colnames(soupx_obj))
  if (length(kept) == 0) stop(sample_id, ": no marker-kept cells found in SoupX object.")
  list(
    obj = subset(soupx_obj, cells = kept),
    status = read_tsv(status_file),
    marker_first_dir = qc_step_dir(sample_dir, "marker_filter")
  )
}

load_top_mito_state <- function(sample_row, marker_obj, sample_dir) {
  sample_id <- sample_row$sample_id[[1]]
  qc_dir <- qc_step_dir(sample_dir, "top_mito")
  status_file <- require_existing_file(file.path(qc_dir, "top_mito_filter_cell_status.tsv"), "topUMI/mito status")
  cache_file <- file.path(qc_step_dir(sample_dir, "objects"), "object_after_marker_top_mito_before_DoubletFinder.rds")
  if (file.exists(cache_file)) {
    obj <- readRDS(cache_file)
  } else {
    kept_file <- require_existing_file(file.path(qc_dir, "after_marker_top_mito_cell_barcodes.txt"), "after topUMI/mito barcodes")
    kept <- intersect(readLines(kept_file), colnames(marker_obj))
    if (length(kept) == 0) stop(sample_id, ": no topUMI/mito-kept cells found in marker-filtered object.")
    obj <- subset(marker_obj, cells = kept)
  }
  list(obj = obj, status = read_tsv(status_file), qc_dir = qc_dir)
}

rerun_final_recluster_from_existing_doubletfinder <- function(sample_row, soupx_before,
                                                              marker_status, top_mito_status,
                                                              after_top_mito_obj,
                                                              sample_dir, config) {
  sample_id <- sample_row$sample_id[[1]]
  rate <- config$doubletfinder$doublet_rate
  rate_label <- rate_label_from_config(config)
  rate_dir <- qc_rate_dir(sample_dir, rate_label)
  dir_df <- file.path(rate_dir, "DoubletFinder_QC")
  dir_final <- qc_final_dir(sample_dir, rate_label)
  dir_summary <- file.path(rate_dir, "summary")
  for (d in c(dir_final, dir_summary)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  df_summary <- read_tsv(require_existing_file(file.path(dir_df, "DoubletFinder_summary.tsv"), "DoubletFinder summary"))
  df_meta <- read_tsv(require_existing_file(file.path(dir_df, "DoubletFinder_metadata.tsv"), "DoubletFinder metadata"))
  singlet_file <- require_existing_file(file.path(dir_df, "after_DoubletFinder_singlet_barcodes.txt"), "DoubletFinder singlet barcodes")
  common_summary_file <- require_existing_file(
    file.path(qc_step_dir(sample_dir, "doubletfinder_common"), "DoubletFinder_common_pK_summary.tsv"),
    "DoubletFinder common pK summary"
  )
  common_summary <- read_tsv(common_summary_file)

  singlets <- intersect(readLines(singlet_file), colnames(after_top_mito_obj))
  if (length(singlets) == 0) stop(sample_id, ": no DoubletFinder singlets found in after-topUMI/mito object.")

  log_msg(sample_id, " ", rate_label, ": final RNA reclustering from existing DoubletFinder singlets.")
  filtered <- subset(after_top_mito_obj, cells = singlets)
  if ("Doublet" %in% colnames(df_meta)) {
    filtered$Doublet <- df_meta$Doublet[match(colnames(filtered), df_meta$cell)]
  } else {
    filtered$Doublet <- "Singlet"
  }

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

  final_rds <- file.path(qc_step_dir(sample_dir, "objects"), "final_filtered_reclustered_object.rds")
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
    pK = common_summary$pK[[1]],
    homotypic_prop = df_summary$homotypic_prop[[1]],
    nExp_poi = df_summary$nExp_poi[[1]],
    nExp_poi_adjusted = df_summary$nExp_poi_adjusted[[1]],
    pANN_cutoff = df_summary$pANN_cutoff[[1]],
    n_doublets = df_summary$doublets[[1]],
    n_singlets = df_summary$singlets[[1]],
    n_final_clusters = length(new_levels),
    final_filtered_rds = final_rds
  )
  write_tsv(rate_summary, file.path(dir_summary, "doublet_rate_QC_summary.tsv"))
  capture.output(sessionInfo(), file = file.path(dir_summary, "sessionInfo.txt"))
  rate_summary
}

pdf_info <- function(sample_dir, rel_path, step, title, parameter_summary = "", n_before = NA, n_after = NA, notes = "") {
  path <- file.path(sample_dir, rel_path)
  data.frame(
    step = step,
    title = title,
    source_pdf = path,
    parameter_summary = parameter_summary,
    n_cells_before = n_before,
    n_cells_after = n_after,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

read_first_row <- function(path) {
  if (file.exists(path)) read_tsv(path)[1, , drop = FALSE] else data.frame()
}

make_report_text_page <- function(path, title, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = 8.5, height = 11)
  grid::grid.newpage()
  grid::grid.text(title, x = 0.05, y = 0.95, just = c("left", "top"),
                  gp = grid::gpar(fontsize = 18, fontface = "bold"))
  y <- 0.88
  for (line in lines) {
    wrapped <- strwrap(as.character(line), width = 95)
    for (w in wrapped) {
      grid::grid.text(w, x = 0.06, y = y, just = c("left", "top"),
                      gp = grid::gpar(fontsize = 9))
      y <- y - 0.025
      if (y < 0.05) {
        grid::grid.newpage()
        y <- 0.95
      }
    }
    y <- y - 0.01
  }
  grDevices::dev.off()
  path
}

write_report_html_index <- function(manifest, html_file) {
  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  rows <- apply(manifest, 1, function(r) {
    src <- r[["source_pdf"]]
    link <- if (file.exists(src)) {
      paste0("<a href=\"", esc(src), "\">", esc(basename(src)), "</a>")
    } else {
      esc(src)
    }
    paste0(
      "<tr><td>", esc(r[["page_order"]]), "</td><td>", esc(r[["step"]]), "</td><td>",
      esc(r[["title"]]), "</td><td>", link, "</td><td>",
      esc(r[["parameter_summary"]]), "</td><td>", esc(r[["n_cells_before"]]),
      "</td><td>", esc(r[["n_cells_after"]]), "</td><td>", esc(r[["notes"]]), "</td></tr>"
    )
  })
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>QC report index</title>",
    "<style>body{font-family:Arial,sans-serif}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ddd;padding:4px;font-size:12px}th{background:#f3f3f3}</style>",
    "</head><body><h1>QC report index</h1>",
    "<table><thead><tr><th>Order</th><th>Step</th><th>Title</th><th>PDF</th><th>Parameters</th><th>Cells before</th><th>Cells after</th><th>Notes</th></tr></thead><tbody>",
    rows,
    "</tbody></table></body></html>"
  )
  writeLines(html, html_file)
}

build_qc_report <- function(sample_row, sample_dir, config, by_rank, marker_cutoff) {
  sample_id <- sample_row$sample_id[[1]]
  rate_label <- rate_label_from_config(config)
  report_dir <- qc_step_dir(sample_dir, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  step_summary <- read_first_row(file.path(qc_step_dir(sample_dir, "run_info"), "QC_step_cell_count_summary.tsv"))
  all_steps <- if (file.exists(file.path(qc_step_dir(sample_dir, "run_info"), "QC_step_cell_count_summary.tsv"))) {
    read_tsv(file.path(qc_step_dir(sample_dir, "run_info"), "QC_step_cell_count_summary.tsv"))
  } else data.frame()
  step_n <- function(step, col) {
    if (!nrow(all_steps) || !step %in% all_steps$step) return(NA)
    all_steps[match(step, all_steps$step), col]
  }

  pdfs <- list()
  add_section <- function(order, step, title, lines) {
    section_pdf <- file.path(report_dir, sprintf("%02d_%s_parameters.pdf", order, gsub("[^A-Za-z0-9]+", "_", step)))
    make_report_text_page(section_pdf, title, lines)
    pdf_info(dirname(section_pdf), basename(section_pdf), step, paste0(title, " parameters"), paste(lines, collapse = "; "))
  }

  pdfs[[length(pdfs) + 1]] <- add_section(0, "00_run_info", paste0(sample_id, " QC overview"), c(
    paste0("sample_id: ", sample_id),
    paste0("input_rds: ", sample_row$input_rds[[1]]),
    paste0("apply_mito_filter: ", sample_row$apply_mito_filter[[1]]),
    paste0("mito_gene_file: ", sample_value(sample_row, "mito_gene_file", "")),
    paste0("by.rank: ", by_rank),
    paste0("marker_filter_cutoff: ", marker_cutoff),
    paste0("doublet_rate: ", config$doubletfinder$doublet_rate),
    paste0("final_recluster_resolution: ", config$final_recluster$resolution)
  ))

  pdfs[[length(pdfs) + 1]] <- add_section(1, "01_emptyDrops", "01 emptyDrops", c(
    paste0("by.rank = ", by_rank),
    paste0("FDR cutoff = ", config$emptydrops$fdr_cutoff),
    paste0("lower = ", config$emptydrops$lower),
    paste0("niters = ", config$emptydrops$niters)
  ))
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("01_emptyDrops", "barcode_rank_plot.pdf"), "01_emptyDrops", "Barcode total UMI rank plot", paste0("by.rank=", by_rank), step_n("raw_matrix", "n_remaining"), step_n("emptyDrops", "n_remaining"))

  pdfs[[length(pdfs) + 1]] <- add_section(2, "02_SoupX", "02 SoupX", c(
    paste0("min.cells = ", config$soupx$min_cells),
    paste0("cluster_resolution = ", config$soupx$cluster_resolution),
    paste0("roundToInt = ", config$soupx$round_to_int),
    paste0("forceAccept = ", config$soupx$force_accept)
  ))
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("02_SoupX", "SoupX_autoEstCont_contamination_fraction.pdf"), "02_SoupX", "SoupX contamination fraction", paste0("forceAccept=", config$soupx$force_accept), step_n("emptyDrops", "n_remaining"), step_n("SoupX", "n_remaining"))
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("02_SoupX", "QC_violin_after_SoupX_before_marker_filter.pdf"), "02_SoupX", "QC violin after SoupX")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("02_SoupX", "QC_feature_scatter_after_SoupX_before_marker_filter.pdf"), "02_SoupX", "QC scatter after SoupX")

  pdfs[[length(pdfs) + 1]] <- add_section(3, "03_marker_filter", "03 marker filter", c(
    paste0("RNA resolution = ", config$marker_filter$rna_resolution),
    paste0("SCT Leiden resolution = ", config$marker_filter$sct_leiden_resolution),
    paste0("marker cutoff = n_markers_one_vs_rest <= ", marker_cutoff),
    paste0("padj < ", config$marker_filter$marker_padj_cutoff),
    paste0("log2FC >= ", config$marker_filter$marker_log2fc_cutoff),
    paste0("pct.1 >= ", config$marker_filter$marker_pct_in_cutoff),
    paste0("delta_pct >= ", config$marker_filter$marker_delta_pct_cutoff)
  ))
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("03_marker_filter", "03a_before_marker_filter", "RNA_UMAP_by_original_cluster_before_marker_filter.pdf"), "03_marker_filter", "RNA UMAP by original cluster before marker filter")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("03_marker_filter", "03a_before_marker_filter", "UMAP_by_original_cluster.pdf"), "03_marker_filter", "SCT UMAP by original cluster")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("03_marker_filter", "03a_before_marker_filter", "UMAP_by_leiden_res4.pdf"), "03_marker_filter", "SCT UMAP by Leiden cluster")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("03_marker_filter", "03a_before_marker_filter", "RNA_DoHeatmap_original_cluster_before_marker_filter.pdf"), "03_marker_filter", "RNA marker heatmap before marker filter")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("03_marker_filter", "03c_removed_cells", "UMAP_removed_cells_by_original_cluster.pdf"), "03_marker_filter", "Removed marker-poor cells on UMAP", paste0("marker_cutoff=", marker_cutoff), step_n("SoupX", "n_remaining"), step_n("marker_filter", "n_remaining"))
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("03_marker_filter", "03c_removed_cells", "RNA_DoHeatmap_before_after_original_cluster.pdf"), "03_marker_filter", "Before/after marker filter heatmap")

  pdfs[[length(pdfs) + 1]] <- add_section(4, "04_topUMI_mito_filter", "04 topUMI / mito filter", c(
    paste0("top_umi_quantile = ", config$top_mito$top_umi_quantile),
    paste0("mito_cutoff = ", config$top_mito$mito_cutoff),
    paste0("apply_mito_filter = ", sample_row$apply_mito_filter[[1]])
  ))

  pdfs[[length(pdfs) + 1]] <- add_section(5, "05_DoubletFinder", "05 DoubletFinder", c(
    paste0("doublet_rate = ", config$doubletfinder$doublet_rate),
    paste0("pN = ", config$doubletfinder$pN),
    paste0("dims_max = ", config$doubletfinder$dims_max),
    paste0("neighbor_dims_max = ", config$doubletfinder$neighbor_dims_max),
    paste0("homotypic adjustment = ", config$doubletfinder$use_homotypic_adjustment)
  ))
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("05_DoubletFinder", "05a_common_pK_sweep", "RNA_UMAP_before_DoubletFinder_by_cluster.pdf"), "05_DoubletFinder", "RNA UMAP before DoubletFinder")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("05_DoubletFinder", paste0("05b_", rate_label), "DoubletFinder_QC", "RNA_UMAP_by_DoubletFinder_class.pdf"), "05_DoubletFinder", "UMAP by DoubletFinder class")

  pdfs[[length(pdfs) + 1]] <- add_section(6, "06_final_recluster", "06 final recluster", c(
    paste0("workflow = ", config$final_recluster$workflow),
    paste0("nfeatures = ", config$final_recluster$nfeatures),
    paste0("dims_max = ", config$final_recluster$dims_max),
    paste0("resolution = ", config$final_recluster$resolution),
    paste0("marker padj < ", config$final_recluster$marker_padj_cutoff),
    paste0("marker logFC > ", config$final_recluster$marker_logfc_cutoff),
    paste0("marker min.pct >= ", config$final_recluster$marker_min_pct)
  ))
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("06_final_recluster", "reclustered_UMAP_by_new_cluster.pdf"), "06_final_recluster", "Final UMAP by new cluster")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("06_final_recluster", "new_cluster_RNA_DoHeatmap_top10_markers.pdf"), "06_final_recluster", "Final RNA marker heatmap")
  pdfs[[length(pdfs) + 1]] <- pdf_info(sample_dir, file.path("06_final_recluster", "final_cluster_UMI_violin_log10.pdf"), "06_final_recluster", "Final cluster UMI violin")

  manifest <- dplyr::bind_rows(pdfs)
  manifest$exists <- file.exists(manifest$source_pdf)
  manifest <- manifest[manifest$exists, , drop = FALSE]
  manifest$page_order <- seq_len(nrow(manifest))
  manifest <- manifest[, c("page_order", "step", "title", "source_pdf", "parameter_summary", "n_cells_before", "n_cells_after", "notes", "exists")]
  manifest_file <- file.path(report_dir, paste0(sample_id, "_QC_report_manifest.tsv"))
  write_tsv(manifest, manifest_file)

  html_file <- file.path(report_dir, paste0(sample_id, "_QC_report_index.html"))
  write_report_html_index(manifest, html_file)

  report_pdf <- file.path(report_dir, paste0(sample_id, "_QC_report_ordered.pdf"))
  pdfunite <- Sys.which("pdfunite")
  if (nzchar(pdfunite) && nrow(manifest) > 0) {
    status <- system2(pdfunite, args = c(manifest$source_pdf, report_pdf), stdout = TRUE, stderr = TRUE)
    writeLines(c("pdfunite output:", status), file.path(report_dir, paste0(sample_id, "_QC_report_merge.log")))
  } else {
    if (file.exists(report_pdf)) unlink(report_pdf)
    writeLines("pdfunite was not found; ordered PDF was not merged. Use manifest/index to inspect source PDFs.",
               file.path(report_dir, paste0(sample_id, "_QC_report_merge.log")))
  }
  invisible(manifest)
}

run_one_sample_qc <- function(sample_row, config, output_root, stop_after, rerun_from) {
  sample_id <- sample_row$sample_id[[1]]
  sample_dir <- file.path(output_root, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_step_dir(sample_dir, "run_info"), recursive = TRUE, showWarnings = FALSE)

  by_rank <- resolve_by_rank(sample_row, config)
  marker_cutoff <- sample_marker_cutoff(sample_row, config)
  options(qc_doublet_rate = config$doubletfinder$doublet_rate)

  if (rerun_from != "none") {
    log_msg(sample_id, ": rerun_from=", rerun_from, "; deleting this step and downstream outputs.")
    delete_outputs_from_step(sample_dir, config, rerun_from)
  }

  write_run_parameters(sample_row, config, sample_dir, by_rank, marker_cutoff)
  write_tsv(data.frame(
    sample_id = sample_id,
    rerun_from = rerun_from,
    stop_after = stop_after,
    run_started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ), file.path(qc_step_dir(sample_dir, "run_info"), "run_control.tsv"))

  log_msg("===== Processing sample ", sample_id, " =====")

  if (rerun_from == "final_recluster") {
    soupx_obj <- load_soupx_object_cache(sample_dir)
    marker_first <- load_marker_filter_state(sample_row, soupx_obj, sample_dir)
    top_mito <- load_top_mito_state(sample_row, marker_first$obj, sample_dir)
    rate_summary <- rerun_final_recluster_from_existing_doubletfinder(
      sample_row = sample_row,
      soupx_before = soupx_obj,
      marker_status = marker_first$status,
      top_mito_status = top_mito$status,
      after_top_mito_obj = top_mito$obj,
      sample_dir = sample_dir,
      config = config
    )
    step_summary <- write_step_count_summary(sample_row, sample_dir, by_rank, marker_cutoff)
    build_qc_report(sample_row, sample_dir, config, by_rank, marker_cutoff)
    capture.output(sessionInfo(), file = file.path(qc_step_dir(sample_dir, "run_info"), "sessionInfo.txt"))
    return(list(rate_summary = rate_summary, step_summary = step_summary))
  }

  log_msg(sample_id, ": reading count matrix from ", sample_row$input_rds[[1]])
  count_matrix <- normalise_count_matrix(readRDS(sample_row$input_rds[[1]]))

  empty_cells <- run_emptydrops_step(sample_row, count_matrix, sample_dir, config, by_rank)
  if (stop_after == "emptydrops") {
    step_summary <- write_step_count_summary(sample_row, sample_dir, by_rank, marker_cutoff)
    return(list(rate_summary = data.frame(), step_summary = step_summary))
  }

  soupx_obj <- run_soupx_step(sample_row, count_matrix, empty_cells, sample_dir, config)
  if (stop_after == "soupx") {
    step_summary <- write_step_count_summary(sample_row, sample_dir, by_rank, marker_cutoff)
    return(list(rate_summary = data.frame(), step_summary = step_summary))
  }

  marker_first <- run_marker_filter_step(sample_row, soupx_obj, sample_dir, config, marker_cutoff)
  if (stop_after == "marker_filter") {
    step_summary <- write_step_count_summary(sample_row, sample_dir, by_rank, marker_cutoff)
    return(list(rate_summary = data.frame(), step_summary = step_summary))
  }

  top_mito <- apply_top_mito_after_marker(sample_row, marker_first$obj, sample_dir, config)
  if (stop_after == "top_mito") {
    step_summary <- write_step_count_summary(sample_row, sample_dir, by_rank, marker_cutoff)
    return(list(rate_summary = data.frame(), step_summary = step_summary))
  }

  df_common <- prepare_doubletfinder_base(top_mito$obj, sample_row, sample_dir, config)
  if (stop_after == "doubletfinder") {
    step_summary <- write_step_count_summary(sample_row, sample_dir, by_rank, marker_cutoff)
    return(list(rate_summary = data.frame(), step_summary = step_summary))
  }

  rate_summary <- run_doublet_and_final_recluster(
    sample_row = sample_row,
    soupx_before = soupx_obj,
    marker_status = marker_first$status,
    top_mito_status = top_mito$status,
    after_top_mito_obj = top_mito$obj,
    df_common = df_common,
    sample_dir = sample_dir,
    config = config
  )
  step_summary <- write_step_count_summary(sample_row, sample_dir, by_rank, marker_cutoff)
  build_qc_report(sample_row, sample_dir, config, by_rank, marker_cutoff)
  capture.output(sessionInfo(), file = file.path(qc_step_dir(sample_dir, "run_info"), "sessionInfo.txt"))
  list(rate_summary = rate_summary, step_summary = step_summary)
}

parsed <- parse_args(args)
sample_tsv <- first_nonempty(parsed$sample_tsv, file.path(framework_root, "config", "sample.tsv"))
config_file <- first_nonempty(parsed$config, file.path(standalone_script_dir, "default_config.R"))
output_root <- first_nonempty(parsed$output_root, file.path(standalone_script_dir, "output"))
only_sample <- first_nonempty(parsed$sample_id, default = NA_character_)
stop_after <- normalise_pipeline_step(first_nonempty(parsed$stop_after, Sys.getenv("QC_STOP_AFTER"), default = "final"), default = "final_recluster")
rerun_from <- normalise_pipeline_step(first_nonempty(parsed$rerun_from, Sys.getenv("QC_RERUN_FROM"), default = "none"), default = "none")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_global_summary_dir(output_root), recursive = TRUE, showWarnings = FALSE)
progress_log <- file.path(qc_global_summary_dir(output_root), "progress.log")
options(qc_progress_log = progress_log)

cat(paste0("\n===== Standalone QC run started at ",
           format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " =====\n"),
    file = progress_log, append = TRUE)

log_msg("Loading config: ", config_file)
config <- load_config(config_file)
if (rerun_from != "none") {
  config$output$resume_from_existing <- TRUE
  config$output$overwrite_existing <- FALSE
}
set.seed(config$random_seed)
log_msg("Loading required R packages.")
load_qc_packages(include_emptydrops = TRUE, include_doubletfinder = TRUE)
log_msg("Required R packages loaded.")

log_msg("Reading sample table: ", sample_tsv)
samples <- read_sample_table(sample_tsv)
if (!is.na(only_sample)) samples <- samples[samples$sample_id == only_sample, , drop = FALSE]
if (nrow(samples) == 0) stop("No samples to run. Check --sample_tsv and --sample_id.")
log_msg("Samples selected: ", paste(samples$sample_id, collapse = ", "))

all_results <- lapply(seq_len(nrow(samples)), function(i) {
  run_one_sample_qc(samples[i, , drop = FALSE], config, output_root, stop_after, rerun_from)
})

write_global_summaries(output_root)
log_msg("Standalone QC results done: ", normalizePath(output_root, mustWork = FALSE))
