suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) y else x
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

read_config <- function(path) {
  if (!file.exists(path)) stopf("Config does not exist: %s", path)
  cfg <- yaml::read_yaml(path)
  cfg$project$root <- normalizePath(cfg$project$root, mustWork = TRUE)
  cfg$.config_file <- normalizePath(path, mustWork = TRUE)
  cfg
}

project_path <- function(cfg, path) {
  if (grepl("^/", path)) path else file.path(cfg$project$root, path)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_species_manifest <- function(cfg, included_only = TRUE) {
  x <- fread(project_path(cfg, cfg$inputs$species_manifest), na.strings = c("NA"))
  required <- c("object_species", "orthology_species", "include", "seurat_file",
                "assay", "counts_layer", "celltype_column")
  missing <- setdiff(required, names(x))
  if (length(missing)) stopf("Species manifest lacks columns: %s", paste(missing, collapse = ", "))
  x[, include := as.logical(include)]
  if (anyDuplicated(x$object_species)) stopf("Duplicate object_species in species manifest")
  if (included_only) x <- x[include == TRUE]
  x
}

read_comparisons <- function(cfg, comparison_id = NULL, included_only = TRUE) {
  x <- fread(project_path(cfg, cfg$inputs$celltype_comparisons), na.strings = c("NA"))
  required <- c("comparison_id", "object_species", "include", "species_order", "target_labels")
  missing <- setdiff(required, names(x))
  if (length(missing)) stopf("Comparison table lacks columns: %s", paste(missing, collapse = ", "))
  x[, include := as.logical(include)]
  if (anyDuplicated(x[, .(comparison_id, object_species)])) {
    stopf("Duplicate comparison_id/object_species rows in comparison table")
  }
  if (!is.null(comparison_id)) {
    target_comparison_id <- comparison_id
    if (!target_comparison_id %in% x$comparison_id) stopf("Unknown comparison_id: %s", target_comparison_id)
    x <- x[x$comparison_id == target_comparison_id]
  }
  if (included_only) x <- x[include == TRUE]
  setorder(x, species_order)
  x
}

split_labels <- function(x) {
  ans <- trimws(strsplit(x, ";", fixed = TRUE)[[1L]])
  unique(ans[nzchar(ans)])
}

write_tsv <- function(x, path, compress = FALSE) {
  ensure_dir(dirname(path))
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA", compress = if (compress) "gzip" else "none")
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

get_counts <- function(seu, assay = "RNA", layer = "counts") {
  SeuratObject::LayerData(seu, assay = assay, layer = layer)
}

comparison_result_dir <- function(cfg, comparison_id) {
  project_path(cfg, file.path(cfg$outputs$comparisons_dir, comparison_id))
}

utc_now <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)
