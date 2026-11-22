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
source(file.path(workflow_root, "R", "pair_consensus.R"))
source(file.path(workflow_root, "R", "sensitivity.R"))
source(file.path(workflow_root, "R", "marker_program.R"))

cli <- parse_cli(
  file.path(workflow_root, "config", "config.yaml"),
  extra = list(
    mode = "pilot",
    iterations = NULL,
    resolutions = NULL,
    count_mode = NULL,
    marker_mode = NULL,
    pair_consensus = NULL,
    sensitivity = NULL,
    postprocess_only = "false"
  )
)
config <- load_config(cli$config)
initialize_run(config, cli$config)
require_stage(config, "02_full_resolution_scan")
dirs <- ensure_stage_dirs(config, "03_resampling", c("tables", "figures", "logs", "checkpoints"))

mode <- match.arg(cli$mode, c("pilot", "formal"))
default_resolutions <- if (mode == "pilot") {
  config$resampling$pilot_resolutions
} else {
  config$resampling$formal_resolutions
}
resolutions <- if (is.null(cli$resolutions)) default_resolutions else as_number_vector(cli$resolutions)
resolutions <- validate_resolution_subset(resolutions, config$clustering$resolution_grid, "--resolutions")
if (!length(resolutions)) abort("No %s resolutions configured. Set them in config or pass --resolutions=.", mode)
n_iterations <- if (is.null(cli$iterations)) {
  if (mode == "pilot") config$resampling$pilot_iterations else config$resampling$formal_iterations
} else {
  as.integer(cli$iterations)
}
if (is.na(n_iterations) || n_iterations < 1L) abort("--iterations must be a positive integer.")
count_mode <- cli$count_mode %||% config$count_input$primary_mode
marker_mode <- cli$marker_mode %||% config$marker$default_mode
marker_mode <- gsub("-", "_", marker_mode, fixed = TRUE)
if (!marker_mode %in% c("none", "fixed_expression", "rediscovery")) {
  abort("--marker-mode must be none, fixed_expression, or rediscovery.")
}
run_pair_consensus_stage <- if (is.null(cli$pair_consensus)) {
  isTRUE(config$pair_consensus$enabled)
} else {
  as_flag(cli$pair_consensus)
}
run_sensitivity_stage <- if (is.null(cli$sensitivity)) {
  isTRUE(config$sensitivity$enabled)
} else {
  as_flag(cli$sensitivity)
}
postprocess_only <- as_flag(cli$postprocess_only)

log_path <- file.path(
  dirs$logs,
  sprintf("%s_%s.log", if (postprocess_only) "postprocess" else "resampling", mode)
)
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Started:", format(Sys.time()), "\n")
cat("Mode:", mode, "iterations:", n_iterations, "\n")
cat("Resolutions:", paste(resolutions, collapse = ", "), "\n")
cat("Count mode:", count_mode, "\n")
cat("Marker mode requested:", marker_mode, "\n")
cat("Pair consensus enabled:", run_pair_consensus_stage, "\n")
cat("Parameter sensitivity enabled:", run_sensitivity_stage, "\n")
cat("Postprocess only:", postprocess_only, "\n")

loaded <- load_and_validate_object(config)
object <- loaded$object
config$._resolved_primary_pcs <- loaded$primary_pcs
primary_pcs <- primary_pcs_from_config(config)
resampling_pcs <- if (isTRUE(config$canonical_preprocessing$enabled)) {
  primary_pcs
} else {
  resolve_pc_setting(config$resampling$pcs, primary_pcs, "resampling.pcs")
}
cat("Resolved full-data primary PCs:", primary_pcs, "\n")
cat("Resampling PCs:", resampling_pcs, "\n")
reference_path <- file.path(stage_dir(config, "02_full_resolution_scan", "tables"), "full_resolution_assignments.tsv.gz")
reference <- read_table_auto(reference_path)
if (!identical(reference$cell_id, colnames(object))) {
  index <- match(colnames(object), reference$cell_id)
  if (anyNA(index)) abort("Reference assignment cells do not match the object.")
  reference <- reference[index, , drop = FALSE]
}
required_columns <- vapply(resolutions, resolution_column, character(1))
missing_columns <- setdiff(required_columns, names(reference))
if (length(missing_columns)) abort("Missing reference assignment columns: %s", paste(missing_columns, collapse = ", "))

tag <- paste(mode, count_mode, sep = "__")
assignment_path <- file.path(dirs$tables, paste0("cluster_assignments__", tag, ".tsv.gz"))
if (postprocess_only) {
  if (!identical(marker_mode, "none")) {
    abort("Postprocess-only mode currently requires --marker-mode=none.")
  }
  if (!file.exists(assignment_path)) {
    abort("Existing resample assignments are required for postprocess-only mode: %s", assignment_path)
  }
  assignments <- read_table_auto(assignment_path)
  missing_resolutions <- setdiff(resolutions, sort(unique(assignments$resolution)))
  if (length(missing_resolutions)) {
    abort(
      "Existing assignments do not contain requested resolutions: %s",
      paste(missing_resolutions, collapse = ", ")
    )
  }
  assignments <- assignments[assignments$resolution %in% resolutions, , drop = FALSE]
  if (run_pair_consensus_stage) {
    cat("Running sampled-pair consensus from existing assignments.\n")
    invisible(run_pair_consensus(
      config, resolutions, reference, assignments,
      loaded$embedding, tag
    ))
  }
  if (run_sensitivity_stage) {
    cat("Running PC/k/seed sensitivity.\n")
    invisible(run_parameter_sensitivity(config, resolutions, object, reference, tag))
  }
  writeLines(
    capture.output(sessionInfo()),
    file.path(dirs$logs, paste0("session_info__postprocess__", tag, ".txt"))
  )
  cat("Finished postprocess-only structural extensions:", format(Sys.time()), "\n")
  cat("No resampling was repeated and no Seurat object was saved.\n")
  quit(save = "no", status = 0)
}

prepared <- prepare_counts(loaded$counts, config, count_mode)
counts <- prepared$counts
write_csv_atomic(
  prepared$audit,
  file.path(dirs$tables, sprintf("count_matrix_audit_%s.csv", count_mode))
)

checkpoint_root <- ensure_dir(file.path(
  dirs$checkpoints,
  paste(mode, count_mode, resolution_set_id(resolutions), sep = "__")
))
signature <- paste(
  file_md5(config$input$seurat_rds),
  file_md5(reference_path),
  resolution_set_id(resolutions),
  mode,
  n_iterations,
  count_mode,
  config$resampling$fraction,
  config$clustering$pca_npcs,
  resampling_pcs,
  config$resampling$k,
  config$compute$seed,
  config$canonical_preprocessing$variable_features_n %||% 3000L,
  "shared_resolution_seed",
  sep = "|"
)
signature_path <- file.path(checkpoint_root, "checkpoint_signature.txt")
if (file.exists(signature_path)) {
  if (!identical(readLines(signature_path, warn = FALSE)[[1L]], signature)) {
    abort("Checkpoint signature mismatch. Use a new run_id or matching parameters.")
  }
} else {
  writeLines(signature, signature_path)
}

checkpoint_paths <- function(iteration) {
  prefix <- file.path(checkpoint_root, sprintf("iteration_%04d", iteration))
  list(
    jaccard = paste0(prefix, "_jaccard.csv"),
    summary = paste0(prefix, "_summary.csv"),
    assignments = paste0(prefix, "_assignments.tsv.gz"),
    complete = paste0(prefix, "_complete.txt")
  )
}

all_cells <- colnames(object)
metadata <- object@meta.data
run_iteration <- function(iteration) {
  paths <- checkpoint_paths(iteration)
  if (file.exists(paths$complete) && !as_flag(cli$overwrite)) return("skipped")
  seed <- config$compute$seed + iteration
  set.seed(seed)
  sampled_cells <- sort(sample(
    all_cells,
    floor(length(all_cells) * config$resampling$fraction),
    replace = FALSE
  ))
  subset_object <- CreateSeuratObject(
    counts = counts[, sampled_cells, drop = FALSE],
    assay = config$input$count_assay,
    meta.data = metadata[sampled_cells, , drop = FALSE]
  )
  subset_object <- SCTransform(
    subset_object,
    assay = config$input$count_assay,
    new.assay.name = config$input$sct_assay,
    do.scale = TRUE,
    do.center = TRUE,
    return.only.var.genes = TRUE,
    variable.features.n = config$canonical_preprocessing$variable_features_n %||% 3000L,
    verbose = FALSE
  )
  DefaultAssay(subset_object) <- config$input$sct_assay
  npcs <- min(config$clustering$pca_npcs, ncol(subset_object) - 1L, nrow(subset_object) - 1L)
  subset_object <- RunPCA(
    subset_object,
    assay = config$input$sct_assay,
    npcs = npcs,
    reduction.name = config$input$reduction,
    reduction.key = config$canonical_preprocessing$reduction_key %||% "PC_",
    verbose = FALSE
  )
  if (resampling_pcs > npcs) {
    abort(
      "Resolved resampling PCs=%d exceeds available subset PCs=%d in iteration %d.",
      resampling_pcs, npcs, iteration
    )
  }
  dims <- seq_len(resampling_pcs)
  subset_object <- FindNeighbors(
    subset_object,
    reduction = config$input$reduction,
    dims = dims,
    k.param = config$resampling$k,
    graph.name = c("resample_nn", "resample_snn"),
    return.neighbor = FALSE,
    verbose = FALSE
  )

  jaccard_rows <- vector("list", length(resolutions))
  summary_rows <- vector("list", length(resolutions))
  assignment_rows <- vector("list", length(resolutions))
  for (i in seq_along(resolutions)) {
    resolution <- resolutions[[i]]
    cluster_column <- paste0("resample_", resolution_column(resolution))
    resolution_seed <- seed
    subset_object <- FindClusters(
      subset_object,
      graph.name = "resample_snn",
      resolution = resolution,
      algorithm = config$clustering$algorithm,
      cluster.name = cluster_column,
      random.seed = resolution_seed,
      verbose = FALSE
    )
    test <- as.character(subset_object[[cluster_column, drop = TRUE]])
    names(test) <- colnames(subset_object)
    test <- test[sampled_cells]
    ref <- as.character(reference[[resolution_column(resolution)]])
    names(ref) <- reference$cell_id
    ref <- ref[sampled_cells]
    jac <- jaccard_reference_clusters(ref, test)
    similarity <- partition_similarity(ref, test)
    jac$mode <- mode
    jac$count_mode <- count_mode
    jac$resolution <- resolution
    jac$iteration <- iteration
    jac$seed <- resolution_seed
    jac$partition_ARI <- similarity$ARI
    jac$partition_NMI <- similarity$NMI
    jaccard_rows[[i]] <- jac
    summary_rows[[i]] <- data.frame(
      mode = mode,
      count_mode = count_mode,
      resolution = resolution,
      iteration = iteration,
      seed = resolution_seed,
      sampled_cells = length(sampled_cells),
      reference_n_clusters = length(unique(ref)),
      resample_n_clusters = length(unique(test)),
      cluster_count_difference = length(unique(test)) - length(unique(ref)),
      partition_ARI = similarity$ARI,
      partition_NMI = similarity$NMI,
      stringsAsFactors = FALSE
    )
    assignment_rows[[i]] <- data.frame(
      cell_id = sampled_cells,
      mode = mode,
      count_mode = count_mode,
      resolution = resolution,
      iteration = iteration,
      reference_cluster = ref,
      resample_cluster = test,
      stringsAsFactors = FALSE
    )
  }
  write_csv_atomic(bind_rows_base(jaccard_rows), paths$jaccard)
  write_csv_atomic(bind_rows_base(summary_rows), paths$summary)
  if (isTRUE(config$resampling$save_assignments) ||
      run_pair_consensus_stage ||
      !identical(marker_mode, "none")) {
    write_tsv_gz_atomic(bind_rows_base(assignment_rows), paths$assignments)
  }
  writeLines(format(Sys.time()), paths$complete)
  rm(subset_object)
  gc()
  "completed"
}

workers <- min(config$compute$workers, n_iterations)
cat("Workers:", workers, "(configured maximum:", config$compute$max_workers, ")\n")
iterations <- seq_len(n_iterations)
if (workers > 1L && .Platform$OS.type != "windows") {
  results <- parallel::mclapply(
    iterations, run_iteration, mc.cores = workers,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  results <- lapply(iterations, run_iteration)
}
failed <- vapply(results, inherits, logical(1), "try-error")
if (any(failed)) abort("Failed resampling iterations: %s", paste(iterations[failed], collapse = ", "))

paths_all <- lapply(iterations, checkpoint_paths)
if (!all(vapply(paths_all, function(x) file.exists(x$complete), logical(1)))) {
  abort("One or more iterations did not produce a completion checkpoint.")
}
jaccard <- bind_rows_base(lapply(paths_all, function(x) utils::read.csv(x$jaccard, check.names = FALSE)))
iteration_summary <- bind_rows_base(lapply(paths_all, function(x) utils::read.csv(x$summary, check.names = FALSE)))
summaries <- summarize_resample_jaccard(
  jaccard,
  stable_threshold = config$resampling$stable_jaccard_threshold,
  min_sampled = config$resampling$min_sampled_cells_per_cluster
)
cluster_summary <- summarize_cluster_jaccard(
  jaccard,
  stable_threshold = config$resampling$stable_jaccard_threshold,
  min_sampled = config$resampling$min_sampled_cells_per_cluster
)

jaccard_path <- file.path(dirs$tables, paste0("cluster_jaccard__", tag, ".csv"))
iteration_path <- file.path(dirs$tables, paste0("iteration_summary__", tag, ".csv"))
resolution_path <- file.path(dirs$tables, paste0("resolution_summary__", tag, ".csv"))
cluster_summary_path <- file.path(dirs$tables, paste0("cluster_stability__", tag, ".csv"))
write_csv_atomic(jaccard, jaccard_path)
write_csv_atomic(summaries$per_iteration, iteration_path)
write_csv_atomic(summaries$per_resolution, resolution_path)
write_csv_atomic(cluster_summary, cluster_summary_path)
assignments <- NULL
if (isTRUE(config$resampling$save_assignments) ||
    run_pair_consensus_stage ||
    !identical(marker_mode, "none")) {
  assignments <- bind_rows_base(lapply(paths_all, function(x) read_table_auto(x$assignments)))
}
if (isTRUE(config$resampling$save_assignments)) {
  write_tsv_gz_atomic(assignments, assignment_path)
}

save_plot(
  ggplot(summaries$per_iteration, aes(factor(resolution), median_cluster_jaccard)) +
    geom_boxplot(outlier.size = 0.5) + theme_bw() +
    labs(x = "Resolution", y = "Median matched-cluster Jaccard"),
  file.path(dirs$figures, paste0("median_jaccard__", tag, ".pdf"))
)
save_plot(
  ggplot(summaries$per_iteration, aes(factor(resolution), stable_cell_fraction)) +
    geom_boxplot(outlier.size = 0.5) + theme_bw() +
    labs(x = "Resolution", y = "Stable-cell fraction"),
  file.path(dirs$figures, paste0("stable_cell_fraction__", tag, ".pdf"))
)

outputs <- c(jaccard_path, iteration_path, resolution_path, cluster_summary_path)
if (file.exists(assignment_path)) outputs <- c(outputs, assignment_path)
write_stage_manifest(
  config, "03_resampling", resolutions,
  inputs = c(config$input$seurat_rds, reference_path),
  outputs = outputs,
  details = list(
    mode = mode,
    iterations = n_iterations,
    count_mode = count_mode,
    marker_mode = marker_mode,
    pair_consensus = run_pair_consensus_stage,
    sensitivity = run_sensitivity_stage,
    selected_primary_pcs = primary_pcs,
    resampling_pcs = resampling_pcs,
    canonical_rebuilt = !is.null(loaded$canonical),
    full_seurat_saved = FALSE
  ),
  manifest_tag = tag
)
writeLines(capture.output(sessionInfo()), file.path(dirs$logs, paste0("session_info__", tag, ".txt")))

if (run_pair_consensus_stage) {
  cat("Running sampled-pair consensus.\n")
  invisible(run_pair_consensus(
    config, resolutions, reference, assignments,
    loaded$embedding, tag
  ))
}
if (run_sensitivity_stage) {
  cat("Running PC/k/seed sensitivity.\n")
  invisible(run_parameter_sensitivity(config, resolutions, object, reference, tag))
}
if (!identical(marker_mode, "none")) {
  cat("Running fixed-expression marker recovery.\n")
  fixed_marker <- run_fixed_expression_marker_recovery(
    config, resolutions, object, reference, assignments, jaccard, tag
  )
  if (identical(marker_mode, "rediscovery")) {
    cat("Running independent per-resample marker rediscovery.\n")
    invisible(run_resample_marker_rediscovery(
      config, resolutions, counts, fixed_marker$reference_markers,
      assignments, jaccard, tag
    ))
  }
}
cat("Finished:", format(Sys.time()), "\n")
cat("No Seurat object was saved.\n")
