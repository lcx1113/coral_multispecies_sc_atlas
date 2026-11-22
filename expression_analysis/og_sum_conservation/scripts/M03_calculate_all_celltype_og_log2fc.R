#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: M03_calculate_all_celltype_og_log2fc.R <analysis_config.yaml> <object_species>", call. = FALSE)
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
module <- cfg$log2fc_all_celltypes
if (is.null(module)) stopf("Config lacks log2fc_all_celltypes")
object_species_arg <- args[[2L]]
manifest <- read_species_manifest(cfg)
if (!object_species_arg %in% manifest$object_species) stopf("Unknown object_species: %s", object_species_arg)
m <- manifest[object_species == object_species_arg]

object_file <- project_path(cfg, file.path(module$source_object_dir, sprintf("%s.OG_SUM.seurat.rds", object_species_arg)))
if (!file.exists(object_file)) stopf("Missing OG_SUM object: %s", object_file)
output_dir <- project_path(cfg, module$output_dir)
species_dir <- file.path(output_dir, "by_species", object_species_arg)
ensure_dir(species_dir)
started <- Sys.time()

message("Reading ", object_file)
og_obj <- readRDS(object_file)
if (!inherits(og_obj, "Seurat")) stopf("Input is not a Seurat object: %s", object_file)
if (!module$assay %in% Assays(og_obj)) stopf("%s lacks assay %s", object_species_arg, module$assay)
if (!module$celltype_column %in% colnames(og_obj[[]])) stopf("%s lacks metadata column %s", object_species_arg, module$celltype_column)

counts <- SeuratObject::LayerData(og_obj, assay = module$assay, layer = module$counts_layer)
data_matrix <- SeuratObject::LayerData(og_obj, assay = module$assay, layer = module$data_layer)
if (!identical(dim(counts), dim(data_matrix)) || !identical(dimnames(counts), dimnames(data_matrix))) {
  stopf("%s counts and data layers differ in dimensions or dimnames", object_species_arg)
}
positive_total <- Matrix::rowSums(counts) > 0
og_ids <- rownames(counts)[positive_total]
if (!length(og_ids)) stopf("%s has no positive-total OG features", object_species_arg)
counts <- counts[og_ids, , drop = FALSE]
data_matrix <- data_matrix[og_ids, , drop = FALSE]

linear_normalized <- data_matrix
if (inherits(linear_normalized, "sparseMatrix")) {
  linear_normalized@x <- expm1(linear_normalized@x)
} else {
  linear_normalized <- expm1(linear_normalized)
}

celltype <- as.character(og_obj[[module$celltype_column, drop = TRUE]])
observed_labels <- sort(unique(celltype[!is.na(celltype) & nzchar(celltype)]))
excluded <- grepl(
  module$exclude_target_regex, observed_labels,
  ignore.case = isTRUE(module$exclude_target_regex_ignore_case), perl = TRUE
)
target_labels <- observed_labels[!excluded]
if (!length(target_labels)) stopf("%s has no eligible target labels", object_species_arg)
if (isTRUE(module$merge_target_labels)) stopf("This module is configured for separate targets; merge_target_labels must be false")

explicit_separate <- module$explicit_separate_labels[[object_species_arg]]
if (!is.null(explicit_separate)) {
  missing_explicit <- setdiff(unlist(explicit_separate), target_labels)
  if (length(missing_explicit)) stopf("%s explicit separate label(s) absent or excluded: %s", object_species_arg, paste(missing_explicit, collapse = ", "))
}

resolved_targets <- data.table(
  object_species = object_species_arg,
  orthology_species = m$orthology_species,
  observed_label = observed_labels,
  target_celltype = observed_labels,
  n_cells = vapply(observed_labels, function(z) sum(celltype == z, na.rm = TRUE), integer(1L)),
  include_as_target = !excluded,
  exclusion_reason = fifelse(excluded, "excluded_by_regex_but_retained_in_background", NA_character_),
  labels_merged = FALSE
)
resolved_targets[object_species == "Hcoe" & observed_label %in% c("calicoblast", "calicoblast-like"),
                 exclusion_reason := "explicitly_kept_as_separate_target"]

results <- vector("list", length(target_labels))
comparison_summaries <- vector("list", length(target_labels))
max_fc_difference <- 0

for (i in seq_along(target_labels)) {
  target_label <- target_labels[[i]]
  target_cells <- which(!is.na(celltype) & celltype == target_label)
  background_cells <- setdiff(seq_along(celltype), target_cells)
  n_target <- length(target_cells)
  n_background <- length(background_cells)
  message(sprintf("[%d/%d] %s: %s (%d target vs %d background)",
                  i, length(target_labels), object_species_arg, target_label, n_target, n_background))

  mean_target <- Matrix::rowMeans(linear_normalized[, target_cells, drop = FALSE])
  mean_background <- Matrix::rowMeans(linear_normalized[, background_cells, drop = FALSE])
  pct_target <- Matrix::rowMeans(data_matrix[, target_cells, drop = FALSE] > 0)
  pct_background <- Matrix::rowMeans(data_matrix[, background_cells, drop = FALSE] > 0)
  # Seurat FoldChange.Assay with slot="data" adds pseudocount.use to the
  # group sum before division by the number of cells.
  descriptive_fc <- log2(
    (mean_target + module$fold_change_pseudocount / n_target) /
      (mean_background + module$fold_change_pseudocount / n_background)
  )

  base_table <- data.table(
    orthogroup = og_ids,
    mean_normalized_target = as.numeric(mean_target),
    mean_normalized_background = as.numeric(mean_background),
    pct_target = as.numeric(pct_target),
    pct_background = as.numeric(pct_background),
    calculated_avg_log2FC = as.numeric(descriptive_fc)
  )

  fm <- data.table()
  if (n_target >= module$min_cells_group && n_background >= module$min_cells_group) {
    Idents(og_obj) <- factor(
      ifelse(seq_len(ncol(og_obj)) %in% target_cells, "target", "background"),
      levels = c("background", "target")
    )
    set.seed(cfg$project$seed)
    fm_raw <- FindMarkers(
      og_obj, ident.1 = "target", ident.2 = "background",
      assay = module$assay, features = og_ids,
      slot = module$slot, fc.slot = module$fc_slot,
      min.pct = module$min_pct, min.diff.pct = module$min_diff_pct,
      logfc.threshold = module$logfc_threshold, only.pos = module$only_pos,
      test.use = module$test_use,
      min.cells.feature = module$min_cells_feature,
      min.cells.group = module$min_cells_group,
      max.cells.per.ident = module$max_cells_per_ident,
      random.seed = cfg$project$seed, densify = FALSE, verbose = FALSE
    )
    fm <- as.data.table(fm_raw, keep.rownames = "orthogroup")
    if ("avg_logFC" %in% names(fm) && !"avg_log2FC" %in% names(fm)) setnames(fm, "avg_logFC", "avg_log2FC")
    keep_columns <- intersect(c("orthogroup", "p_val", "avg_log2FC", "p_val_adj"), names(fm))
    fm <- fm[, ..keep_columns]
    if ("avg_log2FC" %in% names(fm)) setnames(fm, "avg_log2FC", "seurat_avg_log2FC")
  }

  result <- merge(base_table, fm, by = "orthogroup", all.x = TRUE, sort = FALSE)
  if (!"p_val" %in% names(result)) result[, p_val := NA_real_]
  if (!"p_val_adj" %in% names(result)) result[, p_val_adj := NA_real_]
  if (!"seurat_avg_log2FC" %in% names(result)) result[, seurat_avg_log2FC := NA_real_]
  result[, avg_log2FC := fcoalesce(seurat_avg_log2FC, calculated_avg_log2FC)]
  result[, log2fc_source := fifelse(!is.na(seurat_avg_log2FC), "FindMarkers_data", "calculated_from_data_means")]
  if (n_target < module$min_cells_group) {
    result[, wilcox_test_status := "insufficient_target_cells"]
  } else if (n_background < module$min_cells_group) {
    result[, wilcox_test_status := "insufficient_background_cells"]
  } else {
    result[, wilcox_test_status := fcase(
      !is.na(p_val), "wilcox_tested",
      !is.na(seurat_avg_log2FC), "findmarkers_returned_without_p_value",
      default = "fold_change_only_not_returned_by_findmarkers"
    )]
  }
  result[, `:=`(
    object_species = object_species_arg,
    orthology_species = m$orthology_species,
    target_celltype = target_label,
    n_target_cells = n_target,
    n_background_cells = n_background,
    assay = module$assay,
    normalization = "LogNormalize",
    normalization_scale_factor = og_obj@misc$conserved_og_workflow$scale_factor %||% 10000,
    test_use = module$test_use,
    tested_slot = module$slot,
    fc_slot = module$fc_slot
  )]
  if (nrow(fm) && "seurat_avg_log2FC" %in% names(fm)) {
    difference <- max(abs(result$calculated_avg_log2FC - result$seurat_avg_log2FC), na.rm = TRUE)
    if (is.finite(difference)) max_fc_difference <- max(max_fc_difference, difference)
  }
  setcolorder(result, c(
    "object_species", "orthology_species", "target_celltype", "orthogroup",
    "n_target_cells", "n_background_cells",
    "mean_normalized_target", "mean_normalized_background",
    "pct_target", "pct_background", "avg_log2FC",
    "p_val", "p_val_adj", "wilcox_test_status", "log2fc_source",
    "seurat_avg_log2FC", "calculated_avg_log2FC",
    "assay", "normalization", "normalization_scale_factor", "test_use", "tested_slot", "fc_slot"
  ))
  results[[i]] <- result
  comparison_summaries[[i]] <- data.table(
    object_species = object_species_arg, orthology_species = m$orthology_species,
    target_celltype = target_label, n_target_cells = n_target,
    n_background_cells = n_background, n_ogs = length(og_ids),
    n_wilcox_tested = sum(result$wilcox_test_status == "wilcox_tested"),
    n_fold_change_only = sum(result$wilcox_test_status != "wilcox_tested")
  )
}

long <- rbindlist(results, use.names = TRUE, fill = TRUE)
long[, target_order := match(target_celltype, target_labels)]
setorder(long, target_order, orthogroup)
long[, target_order := NULL]
comparison_summary <- rbindlist(comparison_summaries)

if (anyDuplicated(long[, .(target_celltype, orthogroup)])) stopf("Duplicate target/OG rows for %s", object_species_arg)
if (nrow(long) != length(target_labels) * length(og_ids)) stopf("Incomplete target/OG grid for %s", object_species_arg)
if (max_fc_difference > 1e-8) stopf("Calculated and Seurat log2FC differ by %.6g for %s", max_fc_difference, object_species_arg)

long_file <- file.path(species_dir, sprintf("%s.all_celltypes_OG_log2FC.tsv.gz", object_species_arg))
matrix_file <- file.path(species_dir, sprintf("%s.all_celltypes_OG_log2FC_matrix.tsv.gz", object_species_arg))
target_file <- file.path(species_dir, sprintf("%s.resolved_celltype_targets.tsv", object_species_arg))
comparison_file <- file.path(species_dir, sprintf("%s.comparison_summary.tsv", object_species_arg))
run_file <- file.path(species_dir, sprintf("%s.run_record.tsv", object_species_arg))

write_tsv(long, long_file, compress = TRUE)
matrix_table <- dcast(long, orthogroup ~ target_celltype, value.var = "avg_log2FC")
setcolorder(matrix_table, c("orthogroup", target_labels))
write_tsv(matrix_table, matrix_file, compress = TRUE)
write_tsv(resolved_targets, target_file)
write_tsv(comparison_summary, comparison_file)

finished <- Sys.time()
run_record <- data.table(
  object_species = object_species_arg, orthology_species = m$orthology_species,
  source_object = object_file, source_object_sha256 = sha256_file(object_file),
  n_cells = ncol(og_obj), n_ogs = length(og_ids), n_observed_labels = length(observed_labels),
  n_target_celltypes = length(target_labels), n_output_rows = nrow(long),
  n_wilcox_tested_rows = sum(long$wilcox_test_status == "wilcox_tested"),
  n_fold_change_only_rows = sum(long$wilcox_test_status != "wilcox_tested"),
  max_abs_difference_calculated_vs_seurat_log2fc = max_fc_difference,
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  finished_utc = format(finished, tz = "UTC", usetz = TRUE),
  elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
  seurat_version = as.character(packageVersion("Seurat")),
  r_version = R.version.string,
  full_long_table = long_file, log2fc_matrix = matrix_file
)
write_tsv(run_record, run_file)
message(sprintf("Completed %s: %d target cell types x %d OGs = %d rows",
                object_species_arg, length(target_labels), length(og_ids), nrow(long)))
