get_fc_col <- function(marker_table) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_diff")
  fc_col <- intersect(fc_candidates, colnames(marker_table))
  if (length(fc_col) == 0) stop("Marker table lacks avg_log2FC/avg_logFC/avg_diff.")
  fc_col[[1]]
}

empty_marker_table <- function() {
  tibble::tibble(
    p_val = numeric(), avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(),
    p_val_adj = numeric(), cluster = character(), gene = character(),
    marker_log2FC = numeric(), delta_pct = numeric()
  )
}

filter_marker_table <- function(marker_table, config) {
  if (is.null(marker_table) || nrow(marker_table) == 0) return(empty_marker_table())
  fc_col <- get_fc_col(marker_table)
  if (!"gene" %in% colnames(marker_table)) marker_table$gene <- rownames(marker_table)
  marker_table$cluster <- as.character(marker_table$cluster)
  marker_table$marker_log2FC <- marker_table[[fc_col]]
  marker_table$delta_pct <- marker_table$pct.1 - marker_table$pct.2
  marker_table %>%
    dplyr::filter(
      p_val_adj < config$marker_filter$marker_padj_cutoff,
      marker_log2FC >= config$marker_filter$marker_log2fc_cutoff,
      pct.1 >= config$marker_filter$marker_pct_in_cutoff,
      delta_pct >= config$marker_filter$marker_delta_pct_cutoff
    )
}

run_find_all_markers <- function(obj, cluster_col, label, config) {
  log_msg("FindAllMarkers for ", label, ": ", cluster_col)
  DefaultAssay(obj) <- "RNA"
  obj <- join_assay_layers_if_needed(obj, "RNA")
  Idents(obj) <- obj@meta.data[[cluster_col]]
  Idents(obj) <- droplevels(Idents(obj))
  marker_raw <- FindAllMarkers(
    object = obj,
    assay = "RNA",
    slot = "data",
    only.pos = TRUE,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.05,
    max.cells.per.ident = config$marker_filter$max_cells_per_ident,
    random.seed = config$marker_filter$random_seed,
    return.thresh = 1,
    verbose = TRUE
  )
  marker_raw$cluster <- as.character(marker_raw$cluster)
  if (!"gene" %in% colnames(marker_raw)) marker_raw$gene <- rownames(marker_raw)
  marker_pass <- filter_marker_table(marker_raw, config)
  list(obj = obj, raw = marker_raw, pass = marker_pass)
}

build_cluster_marker_summary <- function(obj, cluster_col, marker_pass, marker_cutoff) {
  meta <- obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(.cluster = as.character(.data[[cluster_col]]))
  meta$.nCount <- if ("nCount_RNA" %in% colnames(meta)) meta$nCount_RNA else NA_real_
  meta$.nFeature <- if ("nFeature_RNA" %in% colnames(meta)) meta$nFeature_RNA else NA_real_
  meta$.percent_mito <- if ("percent.mito" %in% colnames(meta)) meta$percent.mito else NA_real_

  marker_count <- marker_pass %>% dplyr::count(cluster, name = "n_markers_one_vs_rest")
  meta %>%
    dplyr::group_by(.cluster) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_nUMI = if (all(is.na(.nCount))) NA_real_ else median(.nCount, na.rm = TRUE),
      median_nFeature = if (all(is.na(.nFeature))) NA_real_ else median(.nFeature, na.rm = TRUE),
      median_percent_mito = if (all(is.na(.percent_mito))) NA_real_ else median(.percent_mito, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(cluster = .cluster) %>%
    dplyr::left_join(marker_count, by = "cluster") %>%
    dplyr::mutate(
      n_markers_one_vs_rest = ifelse(is.na(n_markers_one_vs_rest), 0L, n_markers_one_vs_rest),
      flagged_marker_poor = n_markers_one_vs_rest <= marker_cutoff
    ) %>%
    dplyr::arrange(factor(cluster, levels = numeric_order(cluster)))
}

select_sct_pcs <- function(obj, max_pcs_to_consider = 50) {
  stdev <- obj@reductions$pca@stdev
  max_i <- min(max_pcs_to_consider, length(stdev))
  selected <- 30L
  # Only the source path changes: pca_elbow.R contains the verbatim extracted function.
  source_file <- file.path(framework_root, "R", "pca_elbow.R")
  if (file.exists(source_file)) {
    try({
      source(source_file)
      if (exists("find_pca_elbow")) {
        x <- find_pca_elbow(stdev)
        if ("First derivative" %in% names(x)) selected <- as.integer(x[["First derivative"]])
      }
    }, silent = TRUE)
  }
  if (is.na(selected) || selected < 5L || selected > max_i) selected <- min(30L, max_i)
  selected
}
