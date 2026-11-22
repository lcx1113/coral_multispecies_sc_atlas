run_parameter_sensitivity <- function(config, resolutions, object, reference, tag) {
  assert_packages(c("Seurat", "ggplot2"))
  dirs <- ensure_stage_dirs(config, "05_parameter_sensitivity")
  DefaultAssay(object) <- config$input$sct_assay
  max_pcs <- ncol(Seurat::Embeddings(object, config$input$reduction))
  requested_pcs <- suppressWarnings(as.integer(unlist(config$sensitivity$pcs)))
  if (anyNA(requested_pcs) || any(requested_pcs < 1L)) {
    abort("sensitivity.pcs must contain positive integers.")
  }
  if (isTRUE(config$canonical_preprocessing$include_selected_pc_in_sensitivity)) {
    requested_pcs <- c(requested_pcs, primary_pcs_from_config(config))
  }
  pcs_values <- sort(unique(requested_pcs))
  unavailable_pcs <- pcs_values[pcs_values > max_pcs]
  if (length(unavailable_pcs)) {
    abort(
      "Sensitivity PCs exceed the canonical reduction dimensions (%d): %s",
      max_pcs, paste(unavailable_pcs, collapse = ", ")
    )
  }
  k_values <- as.integer(unlist(config$sensitivity$k))
  seeds <- as.integer(unlist(config$sensitivity$seeds))
  run_rows <- list()
  jaccard_rows <- list()
  index <- 1L
  for (pcs in pcs_values) {
    for (k in k_values) {
      nn_name <- paste("sensitivity", pcs, k, "nn", sep = "_")
      snn_name <- paste("sensitivity", pcs, k, "snn", sep = "_")
      object <- Seurat::FindNeighbors(
        object,
        reduction = config$input$reduction,
        dims = seq_len(pcs),
        k.param = k,
        graph.name = c(nn_name, snn_name),
        return.neighbor = FALSE,
        verbose = FALSE
      )
      for (seed in seeds) {
        for (resolution in resolutions) {
          column <- paste("sensitivity", pcs, k, seed, resolution_column(resolution), sep = "_")
          object <- Seurat::FindClusters(
            object,
            graph.name = snn_name,
            resolution = resolution,
            algorithm = config$clustering$algorithm,
            cluster.name = column,
            random.seed = seed,
            verbose = FALSE
          )
          test <- as.character(object[[column, drop = TRUE]])
          names(test) <- colnames(object)
          ref <- as.character(reference[[resolution_column(resolution)]])
          names(ref) <- reference$cell_id
          ref <- ref[colnames(object)]
          similarity <- partition_similarity(ref, test)
          jac <- jaccard_reference_clusters(ref, test)
          stable <- jac$max_jaccard >= config$resampling$stable_jaccard_threshold
          run_id <- paste(pcs, k, seed, fmt_resolution(resolution), sep = "_")
          run_rows[[index]] <- data.frame(
            run_id = run_id,
            pcs = pcs,
            k = k,
            seed = seed,
            resolution = resolution,
            reference_n_clusters = length(unique(ref)),
            test_n_clusters = length(unique(test)),
            cluster_count_difference = length(unique(test)) - length(unique(ref)),
            reference_ARI = similarity$ARI,
            reference_NMI = similarity$NMI,
            median_cluster_jaccard = safe_median(jac$max_jaccard),
            stable_cell_fraction = stats::weighted.mean(stable, jac$reference_cluster_sampled_n),
            stringsAsFactors = FALSE
          )
          jac$run_id <- run_id
          jac$pcs <- pcs
          jac$k <- k
          jac$seed <- seed
          jac$resolution <- resolution
          jaccard_rows[[index]] <- jac
          index <- index + 1L
        }
      }
    }
  }
  runs <- bind_rows_base(run_rows)
  jaccard <- bind_rows_base(jaccard_rows)
  summary <- bind_rows_base(lapply(split(runs, runs$resolution), function(x) {
    data.frame(
      resolution = x$resolution[[1L]],
      n_runs = nrow(x),
      min_reference_ARI = min(x$reference_ARI, na.rm = TRUE),
      median_reference_ARI = safe_median(x$reference_ARI),
      min_reference_NMI = min(x$reference_NMI, na.rm = TRUE),
      median_reference_NMI = safe_median(x$reference_NMI),
      median_cluster_jaccard = safe_median(x$median_cluster_jaccard),
      min_stable_cell_fraction = min(x$stable_cell_fraction, na.rm = TRUE),
      median_stable_cell_fraction = safe_median(x$stable_cell_fraction),
      max_abs_cluster_count_difference = max(abs(x$cluster_count_difference), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  runs_path <- file.path(dirs$tables, paste0("parameter_sensitivity_runs__", tag, ".csv"))
  jaccard_path <- file.path(dirs$tables, paste0("parameter_sensitivity_cluster_jaccard__", tag, ".csv"))
  summary_path <- file.path(dirs$tables, paste0("parameter_sensitivity_summary__", tag, ".csv"))
  write_csv_atomic(runs, runs_path)
  write_csv_atomic(jaccard, jaccard_path)
  write_csv_atomic(summary, summary_path)
  save_plot(
    ggplot2::ggplot(runs, ggplot2::aes(factor(resolution), reference_ARI)) +
      ggplot2::geom_boxplot(outlier.size = 0.5) + ggplot2::theme_bw() +
      ggplot2::labs(x = "Resolution", y = "ARI versus primary graph parameters"),
    file.path(dirs$figures, paste0("parameter_sensitivity__", tag, ".pdf"))
  )
  write_stage_manifest(
    config, "05_parameter_sensitivity", resolutions,
    inputs = character(), outputs = c(runs_path, jaccard_path, summary_path),
    details = list(pcs = paste(pcs_values, collapse = ","), k = paste(k_values, collapse = ",")),
    manifest_tag = tag
  )
  list(runs = runs, summary = summary)
}
