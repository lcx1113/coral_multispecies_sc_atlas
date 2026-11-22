#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) stop("Usage: 03_find_og_sum_markers.R <analysis_config.yaml> <comparison_id> [object_species]", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
comparison_id_arg <- args[[2L]]
manifest <- read_species_manifest(cfg)
comparison_all <- read_comparisons(cfg, comparison_id_arg)
comparison <- copy(comparison_all)
if (length(args) == 3L) {
  if (!args[[3L]] %in% comparison$object_species) stopf("%s is not included in comparison %s", args[[3L]], comparison_id_arg)
  comparison <- comparison[object_species == args[[3L]]]
}
comparison <- merge(comparison, manifest, by = "object_species", all.x = TRUE, sort = FALSE)
setorder(comparison, species_order)

object_dir <- project_path(cfg, file.path(cfg$outputs$intermediate_dir, "02_og_sum_objects"))
result_dir <- comparison_result_dir(cfg, comparison_id_arg)
full_dir <- file.path(result_dir, "markers", "full")
filtered_dir <- file.path(result_dir, "markers", "filtered")
summary_dir <- file.path(result_dir, "markers", "summaries")
invisible(lapply(c(full_dir, filtered_dir, summary_dir), ensure_dir))

summaries <- list()
for (i in seq_len(nrow(comparison))) {
  x <- comparison[i]
  message(sprintf("[%d/%d] %s: %s versus all other cell types", i, nrow(comparison), x$object_species, comparison_id_arg))
  object_file <- file.path(object_dir, sprintf("%s.OG_SUM.seurat.rds", x$object_species))
  if (!file.exists(object_file)) stopf("Missing OG_SUM object: %s", object_file)
  og_obj <- readRDS(object_file)
  labels <- split_labels(x$target_labels)
  observed <- as.character(og_obj[[x$celltype_column, drop = TRUE]])
  missing_labels <- setdiff(labels, unique(observed))
  if (length(missing_labels)) stopf("%s lacks target label(s): %s", x$object_species, paste(missing_labels, collapse = ", "))
  target <- observed %in% labels
  if (!any(target) || !any(!target)) stopf("%s comparison needs both target and background cells", x$object_species)
  Idents(og_obj) <- factor(ifelse(target, "target", "background"), levels = c("background", "target"))

  set.seed(cfg$project$seed)
  markers <- FindMarkers(
    og_obj, ident.1 = "target", ident.2 = "background", assay = "OG_SUM",
    slot = cfg$markers$slot, fc.slot = cfg$markers$fc_slot,
    min.pct = cfg$markers$min_pct_for_test,
    logfc.threshold = cfg$markers$logfc_threshold_for_test,
    only.pos = cfg$markers$only_pos_for_test,
    test.use = cfg$markers$test_use,
    random.seed = cfg$project$seed,
    verbose = FALSE
  )
  markers <- as.data.table(markers, keep.rownames = "orthogroup")
  if ("avg_logFC" %in% names(markers) && !"avg_log2FC" %in% names(markers)) setnames(markers, "avg_logFC", "avg_log2FC")
  if (!all(c("avg_log2FC", "pct.1", "pct.2", "p_val_adj") %in% names(markers))) {
    stopf("Unexpected FindMarkers columns for %s", x$object_species)
  }
  setnames(markers, c("pct.1", "pct.2"), c("pct_target", "pct_background"))
  feature_meta_file <- file.path(object_dir, sprintf("%s.OG_SUM.feature_metadata.tsv.gz", x$object_species))
  feature_meta <- fread(feature_meta_file)
  markers <- merge(markers, feature_meta, by = "orthogroup", all.x = TRUE, sort = FALSE)
  markers[, `:=`(
    comparison_id = comparison_id_arg,
    object_species = x$object_species,
    orthology_species = x$orthology_species,
    target_labels = paste(labels, collapse = ";"),
    n_target_cells = sum(target),
    n_background_cells = sum(!target),
    aggregation_method = "OG_SUM",
    normalization_method = cfg$normalization$method,
    normalization_scale_factor = cfg$normalization$scale_factor,
    test_use = cfg$markers$test_use,
    tested_slot = cfg$markers$slot
  )]
  setcolorder(markers, c(
    "comparison_id", "object_species", "orthology_species", "target_labels", "orthogroup",
    "p_val", "avg_log2FC", "pct_target", "pct_background", "p_val_adj",
    setdiff(names(markers), c("comparison_id", "object_species", "orthology_species", "target_labels", "orthogroup",
                              "p_val", "avg_log2FC", "pct_target", "pct_background", "p_val_adj"))
  ))
  filtered <- markers[
    avg_log2FC > cfg$markers$filtered_avg_log2fc_gt &
      pct_target > cfg$markers$filtered_pct_target_gt &
      p_val_adj < cfg$markers$filtered_p_val_adj_lt
  ]

  full_file <- file.path(full_dir, sprintf("%s.OG_SUM.all_tested.tsv.gz", x$object_species))
  filtered_file <- file.path(filtered_dir, sprintf("%s.OG_SUM.filtered.tsv", x$object_species))
  write_tsv(markers, full_file, compress = TRUE)
  write_tsv(filtered, filtered_file)
  summaries[[i]] <- data.table(
    comparison_id = comparison_id_arg, object_species = x$object_species,
    orthology_species = x$orthology_species, target_labels = paste(labels, collapse = ";"),
    n_target_cells = sum(target), n_background_cells = sum(!target),
    n_og_features = nrow(og_obj), n_ogs_tested = nrow(markers), n_filtered_marker_ogs = nrow(filtered),
    avg_log2fc_filter = sprintf("> %s", cfg$markers$filtered_avg_log2fc_gt),
    pct_target_filter = sprintf("> %s", cfg$markers$filtered_pct_target_gt),
    p_val_adj_filter = sprintf("< %s", cfg$markers$filtered_p_val_adj_lt),
    full_file = full_file, filtered_file = filtered_file
  )
  write_tsv(summaries[[i]], file.path(summary_dir, sprintf("%s.marker_summary.tsv", x$object_species)))
  rm(og_obj, markers, filtered); gc(FALSE)
}

if (nrow(comparison) == nrow(comparison_all)) {
  write_tsv(rbindlist(summaries), file.path(summary_dir, "marker_summary.tsv"))
}
message(sprintf("Marker analysis completed for comparison: %s", comparison_id_arg))
