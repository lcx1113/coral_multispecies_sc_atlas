options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit) == 0) return(getwd())
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = FALSE))
}

repo_root <- function() {
  normalizePath(file.path(script_dir(), ".."), mustWork = FALSE)
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

load_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

load_qc_packages <- function(include_emptydrops = TRUE, include_doubletfinder = TRUE) {
  pkgs <- c("ggplot2", "patchwork", "Seurat", "SeuratObject", "dplyr",
            "tidyr", "tibble", "SoupX", "Matrix")
  if (include_emptydrops) pkgs <- c(pkgs, "DropletUtils", "SingleCellExperiment")
  if (include_doubletfinder) pkgs <- c(pkgs, "DoubletFinder")
  invisible(lapply(unique(pkgs), load_package))
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "\"")
}

to_bool <- function(x, default = FALSE) {
  if (length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

first_nonempty <- function(..., default = NA_character_) {
  vals <- list(...)
  for (x in vals) {
    if (length(x) > 0 && !is.na(x) && nzchar(as.character(x))) return(x)
  }
  default
}

log_msg <- function(..., progress_log = getOption("qc_progress_log", NULL)) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(line)
  if (!is.null(progress_log)) cat(line, "\n", file = progress_log, append = TRUE)
}

plot_and_save <- function(plot, filename, width = 7, height = 5, limitsize = TRUE) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, limitsize = limitsize)
}

numeric_order <- function(x) {
  x <- unique(as.character(x))
  x[order(suppressWarnings(as.numeric(x)), x)]
}

safe_dims <- function(obj, max_dims = 20) {
  n_cells <- ncol(obj)
  n_features <- nrow(obj)
  max_possible <- max(1L, min(as.integer(max_dims), n_cells - 1L, n_features - 1L))
  seq_len(max_possible)
}

normalise_count_matrix <- function(x) {
  if (is.list(x) && !is.null(x$cm)) x <- x$cm
  if (inherits(x, "Seurat")) {
    assay <- DefaultAssay(x)
    return(GetAssayData(x, assay = assay, layer = "counts"))
  }
  if (inherits(x, "SummarizedExperiment")) return(SummarizedExperiment::assay(x))
  if (!inherits(x, "Matrix")) x <- as(x, "dgCMatrix")
  if (is.null(colnames(x))) colnames(x) <- paste0("barcode_", seq_len(ncol(x)))
  x
}

join_assay_layers_if_needed <- function(obj, assay = "RNA") {
  if (inherits(obj[[assay]], "Assay5")) {
    layer_names <- SeuratObject::Layers(obj[[assay]])
    if (any(grepl("^(counts|data)\\.", layer_names))) {
      obj <- SeuratObject::JoinLayers(obj, assay = assay)
    }
  }
  obj
}

qc_use_numbered_dirs <- function() {
  isTRUE(getOption("qc_numbered_dirs", FALSE))
}

qc_step_dir <- function(sample_dir, step) {
  if (!qc_use_numbered_dirs()) {
    legacy <- c(
      run_info = "summary",
      emptydrops = "emptyDrops_QC",
      soupx = "soupx_cache",
      marker_filter = "marker_first_filter",
      top_mito = "after_marker_top_mito_filter",
      doubletfinder_common = "doubletfinder_common",
      report = "report",
      objects = "objects"
    )
    return(file.path(sample_dir, legacy[[step]] %||% step))
  }
  numbered <- c(
    run_info = "00_run_info",
    emptydrops = "01_emptyDrops",
    soupx = "02_SoupX",
    marker_filter = "03_marker_filter",
    top_mito = "04_topUMI_mito_filter",
    doubletfinder_common = "05_DoubletFinder/05a_common_pK_sweep",
    report = "98_report",
    objects = "90_intermediate_objects"
  )
  file.path(sample_dir, numbered[[step]] %||% step)
}

qc_rate_dir <- function(sample_dir, rate_label) {
  if (qc_use_numbered_dirs()) {
    file.path(sample_dir, "05_DoubletFinder", paste0("05b_", rate_label))
  } else {
    file.path(sample_dir, rate_label)
  }
}

qc_final_dir <- function(sample_dir, rate_label) {
  if (qc_use_numbered_dirs()) {
    file.path(sample_dir, "06_final_recluster")
  } else {
    file.path(qc_rate_dir(sample_dir, rate_label), "final_after_marker_top_mito_doubletfinder_recluster")
  }
}

qc_global_summary_dir <- function(output_root) {
  if (qc_use_numbered_dirs()) file.path(output_root, "99_summary") else file.path(output_root, "summary")
}

load_config <- function(config_file) {
  env <- new.env(parent = baseenv())
  sys.source(config_file, envir = env)
  if (!exists("qc_config", envir = env)) stop("Config file must define qc_config: ", config_file)
  get("qc_config", envir = env)
}

read_sample_table <- function(sample_tsv) {
  x <- read_tsv(sample_tsv)
  required <- c("sample_id", "input_rds", "apply_mito_filter")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) stop("sample.tsv missing required columns: ", paste(missing, collapse = ", "))
  optional <- c("mito_gene_file", "expected_loaded_cells", "emptydrops_rank_multiplier",
                "emptydrops_by_rank", "marker_filter_cutoff", "notes")
  for (nm in optional) if (!nm %in% colnames(x)) x[[nm]] <- NA_character_
  x
}

sample_value <- function(sample_row, column, default = NA_character_) {
  if (!column %in% colnames(sample_row)) return(default)
  val <- sample_row[[column]][[1]]
  if (length(val) == 0 || is.na(val) || !nzchar(as.character(val))) default else val
}

resolve_by_rank <- function(sample_row, config) {
  direct <- sample_value(sample_row, "emptydrops_by_rank")
  if (!is.na(direct)) return(as.integer(as.numeric(direct)))
  config_direct <- config$emptydrops$by_rank
  if (!is.null(config_direct) && length(config_direct) > 0 && !is.na(config_direct)) {
    return(as.integer(config_direct))
  }
  expected <- sample_value(sample_row, "expected_loaded_cells")
  multiplier <- sample_value(sample_row, "emptydrops_rank_multiplier")
  if (is.na(multiplier) && !is.null(config$emptydrops$rank_multiplier)) {
    multiplier <- config$emptydrops$rank_multiplier
  }
  if (!is.na(expected) && !is.na(multiplier)) {
    return(as.integer(round(as.numeric(expected) * as.numeric(multiplier))))
  }
  stop("Cannot determine emptyDrops by.rank for sample ", sample_row$sample_id[[1]],
       ". Provide emptydrops_by_rank or expected_loaded_cells + emptydrops_rank_multiplier.")
}

sample_marker_cutoff <- function(sample_row, config) {
  x <- sample_value(sample_row, "marker_filter_cutoff")
  if (!is.na(x)) return(as.integer(as.numeric(x)))
  as.integer(config$marker_filter$cutoff)
}

read_mito_genes <- function(sample_row) {
  apply_mito <- to_bool(sample_value(sample_row, "apply_mito_filter", "FALSE"))
  mito_file <- sample_value(sample_row, "mito_gene_file")
  if (!apply_mito) return(character(0))
  if (is.na(mito_file) || !file.exists(mito_file)) {
    stop("Mito filter is enabled but mito_gene_file is missing: ", mito_file)
  }
  unique(readLines(mito_file))
}

write_run_parameters <- function(sample_row, config, sample_dir, by_rank, marker_cutoff) {
  dir.create(qc_step_dir(sample_dir, "run_info"), recursive = TRUE, showWarnings = FALSE)
  params <- data.frame(
    parameter = c(
      "sample_id", "input_rds", "mito_gene_file", "apply_mito_filter",
      "emptydrops_by_rank", "emptydrops_fdr_cutoff", "soupx_roundToInt",
      "soupx_forceAccept", "soupx_cluster_resolution", "marker_filter_cutoff",
      "marker_filter_sct_leiden_resolution", "top_umi_quantile", "mito_cutoff",
      "doublet_rate", "doubletfinder_pN", "final_workflow", "final_resolution"
    ),
    value = c(
      sample_row$sample_id[[1]],
      sample_row$input_rds[[1]],
      sample_value(sample_row, "mito_gene_file", ""),
      sample_row$apply_mito_filter[[1]],
      by_rank,
      config$emptydrops$fdr_cutoff,
      config$soupx$round_to_int,
      config$soupx$force_accept,
      config$soupx$cluster_resolution,
      marker_cutoff,
      config$marker_filter$sct_leiden_resolution,
      config$top_mito$top_umi_quantile,
      config$top_mito$mito_cutoff,
      config$doubletfinder$doublet_rate,
      config$doubletfinder$pN,
      config$final_recluster$workflow,
      config$final_recluster$resolution
    )
  )
  write_tsv(params, file.path(qc_step_dir(sample_dir, "run_info"), "run_parameters.tsv"))
}
