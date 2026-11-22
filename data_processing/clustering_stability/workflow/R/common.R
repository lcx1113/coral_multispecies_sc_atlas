`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

abort <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    abort("Missing R packages: %s", paste(missing, collapse = ", "))
  }
}

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    abort("Cannot create directory: %s", path)
  }
  normalizePath(path, mustWork = TRUE)
}

script_file <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(arg)) return(NA_character_)
  normalizePath(sub("^--file=", "", arg[[1L]]), mustWork = TRUE)
}

workflow_root_from_script <- function(levels = 2L) {
  path <- script_file()
  if (is.na(path)) abort("Cannot determine workflow root outside Rscript.")
  root <- dirname(path)
  for (i in seq_len(levels - 1L)) root <- dirname(root)
  normalizePath(root, mustWork = TRUE)
}

parse_cli <- function(default_config, extra = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- c(list(config = default_config, overwrite = FALSE), extra)
  for (arg in args) {
    if (identical(arg, "--overwrite")) {
      out$overwrite <- TRUE
      next
    }
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      abort("Unknown argument: %s", arg)
    }
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pair[[1L]], fixed = TRUE)
    value <- paste(pair[-1L], collapse = "=")
    if (!key %in% names(out)) abort("Unknown option: --%s", pair[[1L]])
    out[[key]] <- value
  }
  out
}

as_flag <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  if (!x %in% c("true", "false", "1", "0", "yes", "no")) {
    abort("Expected a boolean value, got: %s", x)
  }
  x %in% c("true", "1", "yes")
}

as_number_vector <- function(x) {
  if (is.null(x) || !length(x)) return(numeric())
  if (is.numeric(x)) return(sort(unique(as.numeric(x))))
  values <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
  values <- as.numeric(values[nzchar(values)])
  if (anyNA(values)) abort("Invalid numeric list: %s", x)
  sort(unique(values))
}

fmt_resolution <- function(x) {
  sub("\\.?0+$", "", format(as.numeric(x), scientific = FALSE, trim = TRUE))
}

resolution_column <- function(x) {
  paste0("resolution_", gsub("\\.", "_", fmt_resolution(x)))
}

resolution_set_id <- function(values) {
  paste(vapply(sort(unique(as.numeric(values))), fmt_resolution, character(1)), collapse = "__")
}

file_md5 <- function(path) {
  if (!file.exists(path)) abort("File does not exist: %s", path)
  unname(tools::md5sum(path))
}

write_csv_atomic <- function(x, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  utils::write.csv(x, tmp, row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(tmp, path)) abort("Cannot move temporary file to: %s", path)
  invisible(path)
}

write_tsv_gz_atomic <- function(x, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  connection <- gzfile(tmp, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.table(
    x, connection, sep = "\t", row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = ""
  )
  close(connection)
  on.exit(NULL, add = FALSE)
  if (!file.rename(tmp, path)) abort("Cannot move temporary file to: %s", path)
  invisible(path)
}

read_table_auto <- function(path) {
  if (grepl("\\.csv(?:\\.gz)?$", path, ignore.case = TRUE)) {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

read_external_metadata <- function(config) {
  delimiter <- config$input$metadata_delimiter %||% "comma"
  reader <- if (identical(delimiter, "tab")) utils::read.delim else utils::read.csv
  row_names <- if (identical(config$input$metadata_cell_ids %||% "", "row_names")) 1L else NULL
  reader(
    config$input$metadata_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = row_names
  )
}

bind_rows_base <- function(items) {
  items <- Filter(function(x) !is.null(x) && nrow(x) > 0L, items)
  if (!length(items)) return(data.frame())
  names_all <- unique(unlist(lapply(items, names), use.names = FALSE))
  items <- lapply(items, function(x) {
    missing <- setdiff(names_all, names(x))
    for (name in missing) x[[name]] <- NA
    x[names_all]
  })
  do.call(rbind, items)
}

safe_median <- function(x) {
  if (!length(x) || all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
}

safe_quantile <- function(x, probability) {
  if (!length(x) || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probability, na.rm = TRUE, names = FALSE))
}

stage_layout <- c(
  "00_run_manifest",
  "01_input_audit",
  "02_full_resolution_scan",
  "03_resampling",
  "04_pair_consensus",
  "05_parameter_sensitivity",
  "06_marker_program",
  "07_resolution_transitions",
  "08_stability_range",
  "09_targeted_split_validation",
  "10_posthoc_annotation"
)

run_dir_from_config <- function(config) {
  file.path(config$run_root, config$run_id)
}

analysis_contract <- function(config) {
  fields <- c(
    normalizePath(config$input$seurat_rds, mustWork = TRUE),
    file_md5(config$input$seurat_rds),
    normalizePath(config$input$metadata_csv, mustWork = TRUE),
    file_md5(config$input$metadata_csv),
    config$count_input$primary_mode,
    config$count_input$integerization,
    resolution_set_id(config$clustering$resolution_grid),
    config$clustering$primary_pcs,
    config$clustering$primary_k,
    config$clustering$algorithm
  )
  if (isTRUE(config$canonical_preprocessing$enabled)) {
    fields <- c(
      fields,
      "canonical_sct_pca",
      config$canonical_preprocessing$rebuild_sct,
      config$canonical_preprocessing$variable_features_n,
      config$clustering$pca_npcs,
      config$canonical_preprocessing$pc_selection_method,
      config$input$reduction,
      config$compute$seed
    )
  }
  paste(fields, collapse = "|")
}

stage_dir <- function(config, stage, leaf = NULL) {
  if (!stage %in% stage_layout) abort("Unknown stage: %s", stage)
  path <- file.path(run_dir_from_config(config), stage)
  if (!is.null(leaf)) path <- file.path(path, leaf)
  path
}

ensure_stage_dirs <- function(config, stage, leaves = c("tables", "figures", "logs")) {
  paths <- setNames(lapply(leaves, function(x) ensure_dir(stage_dir(config, stage, x))), leaves)
  invisible(paths)
}

initialize_run <- function(config, config_path) {
  run_dir <- ensure_dir(run_dir_from_config(config))
  for (stage in stage_layout) ensure_dir(file.path(run_dir, stage))
  manifest_dir <- ensure_dir(stage_dir(config, "00_run_manifest"))
  manifest_path <- file.path(manifest_dir, "run_manifest.csv")
  current <- data.frame(
    workflow_name = config$workflow_name,
    run_id = config$run_id,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    config_path = normalizePath(config_path, mustWork = TRUE),
    config_md5 = file_md5(config_path),
    input_rds = normalizePath(config$input$seurat_rds, mustWork = TRUE),
    input_rds_md5 = file_md5(config$input$seurat_rds),
    metadata_csv = normalizePath(config$input$metadata_csv, mustWork = TRUE),
    metadata_csv_md5 = file_md5(config$input$metadata_csv),
    count_mode = config$count_input$primary_mode,
    integerization = config$count_input$integerization,
    analysis_contract = analysis_contract(config),
    stringsAsFactors = FALSE
  )
  if (file.exists(manifest_path)) {
    previous <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
    locked <- c("input_rds_md5", "metadata_csv_md5", "count_mode", "integerization")
    if ("analysis_contract" %in% names(previous)) locked <- c(locked, "analysis_contract")
    mismatch <- locked[vapply(locked, function(x) !identical(as.character(previous[[x]][1]), as.character(current[[x]][1])), logical(1))]
    if (length(mismatch)) {
      abort(
        "Run '%s' is already locked with different values for: %s. Use a new run_id.",
        config$run_id, paste(mismatch, collapse = ", ")
      )
    }
  } else {
    write_csv_atomic(current, manifest_path)
    file.copy(config_path, file.path(manifest_dir, "config_frozen.yaml"), overwrite = FALSE)
  }
  history_path <- file.path(manifest_dir, "config_history.csv")
  history_row <- current[c("created_at", "config_path", "config_md5")]
  if (file.exists(history_path)) {
    history <- utils::read.csv(history_path, stringsAsFactors = FALSE)
    if (!current$config_md5 %in% history$config_md5) {
      write_csv_atomic(rbind(history, history_row), history_path)
    }
  } else {
    write_csv_atomic(history_row, history_path)
  }
  snapshot_path <- file.path(manifest_dir, paste0("config_", current$config_md5, ".yaml"))
  if (!file.exists(snapshot_path)) file.copy(config_path, snapshot_path, overwrite = FALSE)
  invisible(run_dir)
}

write_stage_manifest <- function(
    config, stage, resolution_values, inputs, outputs, details = list(), manifest_tag = NULL) {
  path <- stage_dir(config, stage)
  ensure_dir(path)
  suffix <- if (is.null(manifest_tag) || !nzchar(manifest_tag)) "" else paste0("__", manifest_tag)
  files <- c(inputs, outputs)
  files <- files[file.exists(files)]
  file_table <- data.frame(
    role = c(rep("input", length(inputs)), rep("output", length(outputs))),
    path = c(inputs, outputs),
    stringsAsFactors = FALSE
  )
  file_table <- file_table[file.exists(file_table$path), , drop = FALSE]
  file_table$md5 <- vapply(file_table$path, file_md5, character(1))
  write_csv_atomic(file_table, file.path(path, paste0("file_checksums", suffix, ".csv")))
  manifest <- data.frame(
    stage = stage,
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    run_id = config$run_id,
    config_md5 = file_md5(config$._config_path),
    resolution_set_id = resolution_set_id(resolution_values),
    resolutions = paste(vapply(resolution_values, fmt_resolution, character(1)), collapse = ","),
    detail = paste(sprintf("%s=%s", names(details), unlist(details)), collapse = ";"),
    stringsAsFactors = FALSE
  )
  tagged_manifest <- file.path(path, paste0("stage_manifest", suffix, ".csv"))
  write_csv_atomic(manifest, tagged_manifest)
  if (nzchar(suffix)) write_csv_atomic(manifest, file.path(path, "stage_manifest.csv"))
}

require_stage <- function(config, stage) {
  path <- file.path(stage_dir(config, stage), "stage_manifest.csv")
  if (!file.exists(path)) abort("Required stage has not completed: %s", stage)
  invisible(utils::read.csv(path, stringsAsFactors = FALSE))
}

load_config <- function(path) {
  assert_packages("yaml")
  config <- yaml::read_yaml(path)
  required <- c("workflow_name", "run_id", "run_root", "input", "clustering", "resampling")
  missing <- setdiff(required, names(config))
  if (length(missing)) abort("Config is missing sections: %s", paste(missing, collapse = ", "))
  config$compute$workers <- as.integer(config$compute$workers %||% 1L)
  config$compute$max_workers <- as.integer(config$compute$max_workers %||% 30L)
  if (config$compute$workers < 1L || config$compute$workers > config$compute$max_workers) {
    abort("compute.workers must be between 1 and %d", config$compute$max_workers)
  }
  config$clustering$resolution_grid <- sort(unique(as.numeric(unlist(config$clustering$resolution_grid))))
  config$resampling$pilot_resolutions <- sort(unique(as.numeric(unlist(config$resampling$pilot_resolutions))))
  config$resampling$formal_resolutions <- sort(unique(as.numeric(unlist(config$resampling$formal_resolutions))))
  if (isTRUE(config$canonical_preprocessing$enabled)) {
    canonical <- config$canonical_preprocessing
    canonical$rebuild_sct <- canonical$rebuild_sct %||% TRUE
    canonical$variable_features_n <- as.integer(canonical$variable_features_n %||% 3000L)
    canonical$pc_selection_method <- canonical$pc_selection_method %||% "first_derivative"
    canonical$reduction_key <- canonical$reduction_key %||% "stPC_"
    canonical$include_selected_pc_in_sensitivity <-
      canonical$include_selected_pc_in_sensitivity %||% TRUE
    if (!identical(canonical$pc_selection_method, "first_derivative")) {
      abort("canonical_preprocessing.pc_selection_method must be first_derivative.")
    }
    if (is.na(canonical$variable_features_n) || canonical$variable_features_n < 1L) {
      abort("canonical_preprocessing.variable_features_n must be a positive integer.")
    }
    if (!isTRUE(canonical$rebuild_sct)) {
      abort("Canonical mode currently requires canonical_preprocessing.rebuild_sct: true.")
    }
    config$canonical_preprocessing <- canonical
  }
  config$._config_path <- normalizePath(path, mustWork = TRUE)
  config
}

get_assay_matrix <- function(object, assay, layer) {
  assert_packages(c("Seurat", "SeuratObject"))
  tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = layer)
  )
}

find_pca_elbow <- function(sdev) {
  sdev <- as.numeric(sdev)
  if (length(sdev) < 5L || any(!is.finite(sdev))) {
    abort("PCA elbow selection requires at least five finite standard deviations.")
  }
  x <- seq_along(sdev)
  elbow_data <- data.frame(x = x, sdev = sdev)

  sse <- vapply(x, function(k) {
    indicator <- ifelse(x <= k, 0, 1)
    x2 <- (x - k) * indicator
    sum(stats::lm(sdev ~ x + x2)$residuals^2)
  }, numeric(1))
  dim_piecewise <- which.min(sse)

  first_difference <- diff(sdev, differences = 1)
  first_density <- stats::density(first_difference)
  first_runs <- rle(diff(first_density$y) > 0)$lengths
  if (length(first_runs) <= 2L) {
    candidates <- which(first_difference < mean(first_difference))
  } else {
    cutoff <- sum(first_runs[-((length(first_runs) - 1L):length(first_runs))])
    candidates <- which(first_difference < first_density$x[[cutoff]])
  }
  if (!length(candidates)) abort("First-derivative PCA elbow selection returned no candidate.")
  dim_first <- max(candidates) + 1L

  second_difference <- diff(sdev, differences = 2)
  positive_second <- second_difference[second_difference > 0]
  if (!length(positive_second)) abort("Second-derivative PCA elbow selection returned no candidate.")
  second_density <- stats::density(positive_second)
  second_runs <- rle(diff(second_density$y) > 0)$lengths
  if (length(second_runs) <= 2L) {
    dim_second <- which.max(second_difference) + 1L
  } else {
    cutoff <- sum(second_runs[seq_len(min(2L, length(second_runs)))])
    candidates <- which(second_difference > second_density$x[[cutoff]])
    if (!length(candidates)) abort("Second-derivative PCA elbow selection returned no candidate.")
    dim_second <- max(candidates) + 1L
  }

  fits <- vector("list", length(sdev) - 2L)
  residuals <- numeric(length(sdev) - 2L)
  for (i in seq_len(length(sdev) - 2L)) {
    fits[[i]] <- stats::lm(sdev ~ x, data = elbow_data[(i + 1L):length(sdev), , drop = FALSE])
    residuals[[i]] <- sdev[[i]] - stats::predict(fits[[i]], newdata = data.frame(x = i))
  }
  residual_density <- stats::density(residuals)
  residual_runs <- rle(diff(residual_density$y) > 0)$lengths
  if (length(residual_runs) <= 2L) {
    candidates <- which(residuals > mean(residuals))
  } else {
    cutoff <- sum(residual_runs[seq_len(min(2L, length(residual_runs)))])
    candidates <- which(residuals > residual_density$x[[cutoff]])
  }
  if (!length(candidates)) abort("Preceding-residual PCA elbow selection returned no candidate.")
  dim_residual <- max(candidates) + 1L

  point_a <- c(1, sdev[[1L]])
  point_b <- c(length(sdev), sdev[[length(sdev)]])
  distance <- rep(NA_real_, length(sdev))
  for (i in 2:(length(sdev) - 1L)) {
    point_c <- c(i, sdev[[i]])
    determinant <- cbind(rbind(point_a, point_b, point_c), rep(1, 3))
    area <- 0.5 * abs(det(determinant))
    distance[[i]] <- 2 * area / stats::dist(rbind(point_a, point_b))
  }
  dim_perpendicular <- which.max(distance)

  set.seed(2022)
  dim_kmeans <- min(stats::kmeans(sdev, 2)$size) + 1L

  c(
    "Piecewise linear model" = dim_piecewise,
    "First derivative" = dim_first,
    "Second derivative" = dim_second,
    "Preceding residual" = dim_residual,
    "Perpendicular line" = dim_perpendicular,
    "K-means clustering" = dim_kmeans
  )
}

resolve_pc_setting <- function(value, selected_pcs, label) {
  if (length(value) == 1L && is.character(value)) {
    if (!value %in% c("auto", "auto_primary", "auto_first_derivative")) {
      abort("%s has unsupported automatic value: %s", label, value)
    }
    if (length(selected_pcs) != 1L || is.na(selected_pcs)) {
      abort("%s uses an automatic setting before canonical PCs were resolved.", label)
    }
    return(as.integer(selected_pcs))
  }
  resolved <- suppressWarnings(as.integer(value))
  if (length(resolved) != 1L || is.na(resolved) || resolved < 1L) {
    abort("%s must be a positive integer or an automatic PC setting.", label)
  }
  resolved
}

primary_pcs_from_config <- function(config) {
  if (!is.null(config$._resolved_primary_pcs)) {
    return(as.integer(config$._resolved_primary_pcs))
  }
  resolve_pc_setting(config$clustering$primary_pcs, NA_integer_, "clustering.primary_pcs")
}

rebuild_canonical_sct_pca <- function(object, config) {
  assert_packages(c("Seurat", "SeuratObject"))
  canonical <- config$canonical_preprocessing
  count_assay <- config$input$count_assay
  sct_assay <- config$input$sct_assay
  reduction <- config$input$reduction
  original_reduction_assay <- if (reduction %in% names(object@reductions)) {
    object[[reduction]]@assay.used
  } else {
    NA_character_
  }
  original_reduction_dimensions <- if (reduction %in% names(object@reductions)) {
    ncol(Seurat::Embeddings(object, reduction = reduction))
  } else {
    0L
  }

  Seurat::DefaultAssay(object) <- count_assay
  if (sct_assay %in% names(object@assays)) object[[sct_assay]] <- NULL
  set.seed(config$compute$seed)
  object <- Seurat::SCTransform(
    object,
    assay = count_assay,
    new.assay.name = sct_assay,
    do.scale = TRUE,
    do.center = TRUE,
    return.only.var.genes = TRUE,
    variable.features.n = canonical$variable_features_n,
    verbose = FALSE
  )
  Seurat::DefaultAssay(object) <- sct_assay
  npcs <- min(
    as.integer(config$clustering$pca_npcs),
    ncol(object) - 1L,
    nrow(object[[sct_assay]]) - 1L
  )
  if (npcs < 5L) abort("Canonical PCA can calculate only %d PCs.", npcs)
  set.seed(config$compute$seed)
  object <- Seurat::RunPCA(
    object,
    assay = sct_assay,
    npcs = npcs,
    reduction.name = reduction,
    reduction.key = canonical$reduction_key,
    verbose = FALSE
  )
  elbow_candidates <- find_pca_elbow(Seurat::Stdev(object, reduction = reduction))
  selected_pcs <- as.integer(elbow_candidates[["First derivative"]])
  if (is.na(selected_pcs) || selected_pcs < 1L || selected_pcs > npcs) {
    abort("Canonical first-derivative PC selection returned invalid value: %s", selected_pcs)
  }

  canonical_record <- list(
    enabled = TRUE,
    rebuild_sct = TRUE,
    count_assay = count_assay,
    sct_assay = sct_assay,
    reduction = reduction,
    reduction_key = canonical$reduction_key,
    variable_features_n = canonical$variable_features_n,
    pca_npcs = npcs,
    pc_selection_method = canonical$pc_selection_method,
    selected_pcs = selected_pcs,
    elbow_candidates = elbow_candidates,
    original_reduction_assay = original_reduction_assay,
    original_reduction_dimensions = original_reduction_dimensions
  )
  object@misc$cluster_stability_canonical <- canonical_record
  list(object = object, record = canonical_record)
}

canonical_audit_table <- function(record) {
  if (is.null(record) || !isTRUE(record$enabled)) {
    return(data.frame(item = "canonical_preprocessing_enabled", value = FALSE))
  }
  elbow <- record$elbow_candidates
  data.frame(
    item = c(
      "canonical_preprocessing_enabled",
      "rebuild_sct",
      "count_assay",
      "sct_assay",
      "reduction",
      "original_reduction_assay",
      "original_reduction_dimensions",
      "variable_features_n",
      "pca_npcs",
      "pc_selection_method",
      "selected_pcs",
      paste0("elbow_", gsub("[^a-z0-9]+", "_", tolower(names(elbow))))
    ),
    value = c(
      TRUE,
      record$rebuild_sct,
      record$count_assay,
      record$sct_assay,
      record$reduction,
      record$original_reduction_assay,
      record$original_reduction_dimensions,
      record$variable_features_n,
      record$pca_npcs,
      record$pc_selection_method,
      record$selected_pcs,
      unname(elbow)
    ),
    stringsAsFactors = FALSE
  )
}

load_and_validate_object <- function(config) {
  assert_packages(c("Seurat", "SeuratObject", "Matrix"))
  object <- readRDS(config$input$seurat_rds)
  if (!inherits(object, "Seurat")) abort("Input RDS is not a Seurat object.")
  if (!config$input$count_assay %in% names(object@assays)) {
    abort("Count assay '%s' is absent.", config$input$count_assay)
  }
  counts <- get_assay_matrix(object, config$input$count_assay, config$input$count_layer)
  if (!identical(colnames(counts), colnames(object))) abort("Count matrix cells do not match object cells.")

  canonical_record <- NULL
  if (isTRUE(config$canonical_preprocessing$enabled)) {
    rebuilt <- rebuild_canonical_sct_pca(object, config)
    object <- rebuilt$object
    canonical_record <- rebuilt$record
  } else if (!config$input$sct_assay %in% names(object@assays)) {
    abort("SCT assay '%s' is absent.", config$input$sct_assay)
  }
  if (!config$input$reduction %in% names(object@reductions)) {
    abort("Reduction '%s' is absent.", config$input$reduction)
  }
  embedding <- Seurat::Embeddings(object, reduction = config$input$reduction)
  selected_pcs <- if (isTRUE(config$canonical_preprocessing$enabled)) {
    canonical_record$selected_pcs
  } else {
    resolve_pc_setting(config$clustering$primary_pcs, NA_integer_, "clustering.primary_pcs")
  }
  if (ncol(embedding) < selected_pcs) {
    abort("Reduction has %d dimensions, fewer than selected primary PCs=%d.", ncol(embedding), selected_pcs)
  }
  if (!identical(object[[config$input$reduction]]@assay.used, config$input$sct_assay)) {
    abort(
      "Reduction '%s' uses assay '%s', expected canonical assay '%s'.",
      config$input$reduction,
      object[[config$input$reduction]]@assay.used,
      config$input$sct_assay
    )
  }
  list(
    object = object,
    counts = counts,
    embedding = embedding,
    primary_pcs = selected_pcs,
    canonical = canonical_record
  )
}

prepare_counts <- function(counts, config, mode = config$count_input$primary_mode) {
  noninteger <- abs(counts@x - round(counts@x)) > 1e-8
  audit <- data.frame(
    mode = mode,
    n_features = nrow(counts),
    n_cells = ncol(counts),
    nonzero_entries = length(counts@x),
    noninteger_nonzero_entries = sum(noninteger),
    noninteger_nonzero_fraction = mean(noninteger),
    matrix_sum_before = sum(counts),
    matrix_sum_after = sum(counts),
    positive_entries_rounded_to_zero = 0,
    integerization = "none",
    stringsAsFactors = FALSE
  )
  if (identical(mode, "soupx_fractional_full_pipeline")) {
    if (!isTRUE(config$count_input$allow_fractional_counts)) {
      abort("Fractional count mode is disabled by count_input.allow_fractional_counts.")
    }
    return(list(counts = counts, audit = audit))
  }
  if (identical(mode, "raw_integer_full_pipeline")) {
    if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
      abort("Raw integer counts contain non-finite or negative values.")
    }
    if (any(noninteger)) {
      abort("raw_integer_full_pipeline requires integer-valued counts.")
    }
    return(list(counts = counts, audit = audit))
  }
  if (!identical(mode, "soupx_integer_full_pipeline")) {
    abort("Unsupported count mode: %s", mode)
  }
  method <- config$count_input$integerization_fallback %||% "round_nonnegative"
  if (!identical(method, "round_nonnegative")) {
    abort("Unsupported integerization fallback: %s", method)
  }
  before <- counts@x
  counts@x <- round(pmax(counts@x, 0))
  counts <- Matrix::drop0(counts)
  audit$matrix_sum_after <- sum(counts)
  audit$positive_entries_rounded_to_zero <- sum(before > 0 & round(before) == 0)
  audit$integerization <- method
  list(counts = counts, audit = audit)
}

validate_resolution_subset <- function(values, grid, label) {
  values <- sort(unique(as.numeric(values)))
  invalid <- values[!vapply(values, function(x) any(abs(grid - x) < 1e-10), logical(1))]
  if (length(invalid)) {
    abort("%s contains resolutions absent from the scan grid: %s", label, paste(invalid, collapse = ", "))
  }
  values
}

save_plot <- function(plot, path, width = 8, height = 5) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot = plot, width = width, height = height, units = "in", dpi = 160)
}
