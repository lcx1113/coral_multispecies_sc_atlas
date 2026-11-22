linear_row_means <- function(matrix, cells) {
  if (!length(cells)) return(rep(NA_real_, nrow(matrix)))
  value <- as(matrix[, cells, drop = FALSE], "dgCMatrix")
  value@x <- expm1(value@x)
  as.numeric(Matrix::rowMeans(value))
}

marker_effect <- function(expression, in_cells, out_cells) {
  log2((linear_row_means(expression, in_cells) + 1) /
         (linear_row_means(expression, out_cells) + 1))
}

find_reference_markers <- function(config, object, reference, resolutions) {
  assay <- config$marker$assay
  layer <- config$marker$layer
  if (!assay %in% Seurat::Assays(object)) abort("Marker assay '%s' is absent.", assay)
  expression <- get_assay_matrix(object, assay, layer)
  rows <- list()
  index <- 1L
  for (resolution in resolutions) {
    labels <- factor(reference[[resolution_column(resolution)]])
    names(labels) <- reference$cell_id
    Seurat::Idents(object) <- labels[colnames(object)]
    markers <- Seurat::FindAllMarkers(
      object,
      assay = assay,
      slot = layer,
      only.pos = TRUE,
      test.use = "wilcox",
      logfc.threshold = config$marker$log2fc_threshold,
      min.pct = config$marker$min_pct,
      verbose = FALSE
    )
    if (!nrow(markers)) next
    effect_column <- if ("avg_log2FC" %in% names(markers)) "avg_log2FC" else "avg_logFC"
    for (cluster in sort(unique(as.character(labels)))) {
      current <- markers[as.character(markers$cluster) == cluster, , drop = FALSE]
      current <- head(current[order(current[[effect_column]], decreasing = TRUE), , drop = FALSE], config$marker$top_n)
      if (!nrow(current)) next
      rows[[index]] <- data.frame(
        resolution = resolution,
        reference_cluster = cluster,
        gene = current$gene,
        rank = seq_len(nrow(current)),
        reference_avg_log2FC = current[[effect_column]],
        pct_in = current$pct.1,
        pct_out = current$pct.2,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  list(markers = bind_rows_base(rows), expression = expression)
}

run_fixed_expression_marker_recovery <- function(
    config, resolutions, object, reference, resampled, jaccard, tag) {
  assert_packages(c("Seurat", "Matrix", "ggplot2"))
  root <- ensure_dir(stage_dir(config, "06_marker_program", "fixed_expression"))
  dirs <- setNames(lapply(c("tables", "figures", "logs"), function(x) ensure_dir(file.path(root, x))),
                   c("tables", "figures", "logs"))
  reference_result <- find_reference_markers(config, object, reference, resolutions)
  reference_markers <- reference_result$markers
  expression <- reference_result$expression
  if (!nrow(reference_markers)) abort("No full-data reference markers passed the configured thresholds.")
  reference_path <- file.path(dirs$tables, paste0("reference_markers__", tag, ".csv"))
  write_csv_atomic(reference_markers, reference_path)

  rows <- list()
  index <- 1L
  for (row in seq_len(nrow(jaccard))) {
    match_row <- jaccard[row, , drop = FALSE]
    if (is.na(match_row$max_jaccard) ||
        match_row$max_jaccard < config$resampling$stable_jaccard_threshold) next
    resolution <- match_row$resolution
    iteration <- match_row$iteration
    reference_cluster <- as.character(match_row$reference_cluster)
    markers <- reference_markers[
      abs(reference_markers$resolution - resolution) < 1e-10 &
        reference_markers$reference_cluster == reference_cluster, , drop = FALSE
    ]
    if (nrow(markers) < config$marker$minimum_genes) next
    current <- resampled[
      abs(resampled$resolution - resolution) < 1e-10 &
        resampled$iteration == iteration, , drop = FALSE
    ]
    in_cells <- current$cell_id[
      as.character(current$resample_cluster) == as.character(match_row$matched_cluster)
    ]
    out_cells <- setdiff(current$cell_id, in_cells)
    if (length(in_cells) < config$clustering$min_cluster_cells || !length(out_cells)) next
    matrix <- expression[markers$gene, , drop = FALSE]
    effect <- marker_effect(matrix, in_cells, out_cells)
    pct_in <- as.numeric(Matrix::rowMeans(matrix[, in_cells, drop = FALSE] > 0))
    pct_out <- as.numeric(Matrix::rowMeans(matrix[, out_cells, drop = FALSE] > 0))
    recovered <- effect >= config$marker$log2fc_threshold &
      pct_in >= config$marker$min_pct &
      effect > 0
    rows[[index]] <- data.frame(
      resolution = resolution,
      iteration = iteration,
      reference_cluster = reference_cluster,
      matched_cluster = as.character(match_row$matched_cluster),
      match_jaccard = match_row$max_jaccard,
      matched_cluster_n = length(in_cells),
      marker_effect_recovery_fraction = mean(recovered),
      positive_effect_fraction = mean(effect > 0),
      median_marker_log2FC = safe_median(effect),
      effect_rank_correlation = suppressWarnings(stats::cor(
        markers$reference_avg_log2FC, effect, method = "spearman", use = "pairwise.complete.obs"
      )),
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }
  details <- bind_rows_base(rows)
  if (!nrow(details)) abort("No matched resamples were eligible for fixed-expression marker recovery.")
  cluster_summary <- bind_rows_base(lapply(
    split(details, interaction(details$resolution, details$reference_cluster, drop = TRUE)),
    function(x) data.frame(
      resolution = x$resolution[[1L]],
      reference_cluster = x$reference_cluster[[1L]],
      n_valid_resamples = nrow(x),
      median_match_jaccard = safe_median(x$match_jaccard),
      median_marker_effect_recovery = safe_median(x$marker_effect_recovery_fraction),
      median_positive_effect_fraction = safe_median(x$positive_effect_fraction),
      median_effect_rank_correlation = safe_median(x$effect_rank_correlation),
      stringsAsFactors = FALSE
    )
  ))
  resolution_summary <- bind_rows_base(lapply(split(cluster_summary, cluster_summary$resolution), function(x) {
    data.frame(
      resolution = x$resolution[[1L]],
      n_clusters = nrow(x),
      median_marker_effect_recovery = safe_median(x$median_marker_effect_recovery),
      median_effect_rank_correlation = safe_median(x$median_effect_rank_correlation),
      stringsAsFactors = FALSE
    )
  }))
  details_path <- file.path(dirs$tables, paste0("marker_recovery_details__", tag, ".tsv.gz"))
  cluster_path <- file.path(dirs$tables, paste0("marker_recovery_cluster_summary__", tag, ".csv"))
  resolution_path <- file.path(dirs$tables, paste0("marker_recovery_summary__", tag, ".csv"))
  write_tsv_gz_atomic(details, details_path)
  write_csv_atomic(cluster_summary, cluster_path)
  write_csv_atomic(resolution_summary, resolution_path)
  save_plot(
    ggplot2::ggplot(cluster_summary, ggplot2::aes(factor(resolution), median_marker_effect_recovery)) +
      ggplot2::geom_boxplot(outlier.size = 0.5) + ggplot2::theme_bw() +
      ggplot2::labs(x = "Resolution", y = "Fixed-expression marker recovery"),
    file.path(dirs$figures, paste0("marker_recovery__", tag, ".pdf"))
  )
  write_stage_manifest(
    config, "06_marker_program", resolutions,
    inputs = character(), outputs = c(reference_path, details_path, cluster_path, resolution_path),
    details = list(marker_mode = "fixed_expression", expression_source = paste(config$marker$assay, config$marker$layer, sep = "/")),
    manifest_tag = paste0("fixed_expression__", tag)
  )
  list(reference_markers = reference_markers, details = details, summary = resolution_summary)
}

run_resample_marker_rediscovery <- function(
    config, resolutions, counts, reference_markers, resampled, jaccard, tag) {
  assert_packages(c("Seurat", "ggplot2"))
  configured <- as.numeric(unlist(config$marker$rediscovery_resolutions))
  if (length(configured)) resolutions <- intersect(resolutions, configured)
  if (!length(resolutions)) abort("No resolutions are eligible for marker rediscovery.")
  root <- ensure_dir(stage_dir(config, "06_marker_program", "rediscovery"))
  dirs <- setNames(lapply(c("tables", "figures", "logs"), function(x) ensure_dir(file.path(root, x))),
                   c("tables", "figures", "logs"))
  rows <- list()
  index <- 1L
  for (resolution in resolutions) {
    available_iterations <- sort(unique(resampled$iteration[abs(resampled$resolution - resolution) < 1e-10]))
    available_iterations <- head(available_iterations, config$marker$rediscovery_iterations)
    for (iteration in available_iterations) {
      current <- resampled[
        abs(resampled$resolution - resolution) < 1e-10 &
          resampled$iteration == iteration, , drop = FALSE
      ]
      cells <- current$cell_id
      subset_object <- Seurat::CreateSeuratObject(
        counts = counts[, cells, drop = FALSE],
        assay = config$input$count_assay
      )
      subset_object <- Seurat::NormalizeData(subset_object, verbose = FALSE)
      Seurat::Idents(subset_object) <- factor(current$resample_cluster[match(colnames(subset_object), current$cell_id)])
      markers <- Seurat::FindAllMarkers(
        subset_object,
        assay = config$input$count_assay,
        slot = "data",
        only.pos = TRUE,
        test.use = "wilcox",
        logfc.threshold = config$marker$log2fc_threshold,
        min.pct = config$marker$min_pct,
        verbose = FALSE
      )
      if (!nrow(markers)) next
      effect_column <- if ("avg_log2FC" %in% names(markers)) "avg_log2FC" else "avg_logFC"
      matches <- jaccard[
        abs(jaccard$resolution - resolution) < 1e-10 &
          jaccard$iteration == iteration &
          jaccard$max_jaccard >= config$resampling$stable_jaccard_threshold, , drop = FALSE
      ]
      for (match_index in seq_len(nrow(matches))) {
        ref_cluster <- as.character(matches$reference_cluster[[match_index]])
        test_cluster <- as.character(matches$matched_cluster[[match_index]])
        ref <- reference_markers[
          abs(reference_markers$resolution - resolution) < 1e-10 &
            reference_markers$reference_cluster == ref_cluster, , drop = FALSE
        ]
        test <- markers[as.character(markers$cluster) == test_cluster, , drop = FALSE]
        test <- head(test[order(test[[effect_column]], decreasing = TRUE), , drop = FALSE], config$marker$top_n)
        if (nrow(ref) < config$marker$minimum_genes || nrow(test) < config$marker$minimum_genes) next
        common <- intersect(ref$gene, test$gene)
        union <- union(ref$gene, test$gene)
        ref_rank <- setNames(ref$rank, ref$gene)
        test_rank <- setNames(seq_len(nrow(test)), test$gene)
        rows[[index]] <- data.frame(
          resolution = resolution,
          iteration = iteration,
          reference_cluster = ref_cluster,
          matched_cluster = test_cluster,
          match_jaccard = matches$max_jaccard[[match_index]],
          reference_marker_n = nrow(ref),
          rediscovered_marker_n = nrow(test),
          marker_set_jaccard = length(common) / length(union),
          marker_rank_correlation = if (length(common) >= 3L) {
            suppressWarnings(stats::cor(ref_rank[common], test_rank[common], method = "spearman"))
          } else NA_real_,
          rediscovered_reference_fraction = length(common) / nrow(ref),
          stringsAsFactors = FALSE
        )
        index <- index + 1L
      }
      rm(subset_object)
      gc()
    }
  }
  details <- bind_rows_base(rows)
  if (!nrow(details)) abort("No marker rediscovery comparisons passed the eligibility rules.")
  summary <- bind_rows_base(lapply(split(details, details$resolution), function(x) {
    data.frame(
      resolution = x$resolution[[1L]],
      n_cluster_resample_comparisons = nrow(x),
      median_marker_set_jaccard = safe_median(x$marker_set_jaccard),
      median_marker_rank_correlation = safe_median(x$marker_rank_correlation),
      median_rediscovered_reference_fraction = safe_median(x$rediscovered_reference_fraction),
      stringsAsFactors = FALSE
    )
  }))
  details_path <- file.path(dirs$tables, paste0("marker_set_jaccard__", tag, ".tsv.gz"))
  summary_path <- file.path(dirs$tables, paste0("marker_rank_and_effect_summary__", tag, ".csv"))
  write_tsv_gz_atomic(details, details_path)
  write_csv_atomic(summary, summary_path)
  save_plot(
    ggplot2::ggplot(details, ggplot2::aes(factor(resolution), marker_set_jaccard)) +
      ggplot2::geom_boxplot(outlier.size = 0.5) + ggplot2::theme_bw() +
      ggplot2::labs(x = "Resolution", y = "Rediscovered top-marker Jaccard"),
    file.path(dirs$figures, paste0("marker_rediscovery__", tag, ".pdf"))
  )
  write_stage_manifest(
    config, "06_marker_program", resolutions,
    inputs = character(), outputs = c(details_path, summary_path),
    details = list(marker_mode = "rediscovery", normalization = "independent_LogNormalize_per_resample"),
    manifest_tag = paste0("rediscovery__", tag)
  )
  list(details = details, summary = summary)
}
