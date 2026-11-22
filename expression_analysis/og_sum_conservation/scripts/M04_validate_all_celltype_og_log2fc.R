#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: M04_validate_all_celltype_og_log2fc.R <analysis_config.yaml>", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
module <- cfg$log2fc_all_celltypes
manifest <- read_species_manifest(cfg)
output_dir <- project_path(cfg, module$output_dir)
ensure_dir(output_dir)
checks <- list()
summaries <- list()
targets <- list()
output_files <- list()

for (i in seq_len(nrow(manifest))) {
  sp <- manifest$object_species[i]
  species_dir <- file.path(output_dir, "by_species", sp)
  long_file <- file.path(species_dir, sprintf("%s.all_celltypes_OG_log2FC.tsv.gz", sp))
  matrix_file <- file.path(species_dir, sprintf("%s.all_celltypes_OG_log2FC_matrix.tsv.gz", sp))
  target_file <- file.path(species_dir, sprintf("%s.resolved_celltype_targets.tsv", sp))
  comparison_file <- file.path(species_dir, sprintf("%s.comparison_summary.tsv", sp))
  run_file <- file.path(species_dir, sprintf("%s.run_record.tsv", sp))
  paths <- c(long_file, matrix_file, target_file, comparison_file, run_file)
  if (!all(file.exists(paths))) stopf("Missing LOG2FC outputs for %s", sp)
  long <- fread(long_file, na.strings = c("NA", ""))
  matrix_table <- fread(matrix_file, na.strings = c("NA", ""))
  target <- fread(target_file, na.strings = c("NA", ""))
  comparison <- fread(comparison_file, na.strings = c("NA", ""))
  run <- fread(run_file, na.strings = c("NA", ""))
  included_targets <- target[include_as_target == TRUE, target_celltype]
  expected_rows <- run$n_ogs * length(included_targets)
  formula_difference <- max(abs(
    long$avg_log2FC - log2(
      (long$mean_normalized_target + module$fold_change_pseudocount / long$n_target_cells) /
        (long$mean_normalized_background + module$fold_change_pseudocount / long$n_background_cells)
    )
  ), na.rm = TRUE)
  required <- c("object_species", "target_celltype", "orthogroup", "avg_log2FC",
                "mean_normalized_target", "mean_normalized_background",
                "pct_target", "pct_background", "p_val", "p_val_adj", "wilcox_test_status")
  required_ok <- !length(setdiff(required, names(long)))
  grid_ok <- nrow(long) == expected_rows && !anyDuplicated(long[, .(target_celltype, orthogroup)])
  targets_ok <- setequal(unique(long$target_celltype), included_targets) &&
    !any(grepl(module$exclude_target_regex, included_targets,
               ignore.case = isTRUE(module$exclude_target_regex_ignore_case), perl = TRUE))
  matrix_ok <- nrow(matrix_table) == run$n_ogs &&
    setequal(setdiff(names(matrix_table), "orthogroup"), included_targets)
  hcoe_ok <- sp != "Hcoe" || all(c("calicoblast", "calicoblast-like") %in% included_targets)
  finite_ok <- all(is.finite(long$avg_log2FC)) &&
    all(long$pct_target >= 0 & long$pct_target <= 1) &&
    all(long$pct_background >= 0 & long$pct_background <= 1)
  formula_ok <- is.finite(formula_difference) && formula_difference <= 1e-10
  status <- if (all(required_ok, grid_ok, targets_ok, matrix_ok, hcoe_ok, finite_ok, formula_ok)) "PASS" else "FAIL"
  checks[[i]] <- data.table(
    object_species = sp, n_ogs = run$n_ogs, n_target_celltypes = length(included_targets),
    n_rows = nrow(long), expected_rows = expected_rows,
    required_columns_ok = required_ok, complete_unique_grid = grid_ok,
    targets_ok = targets_ok, matrix_dimensions_ok = matrix_ok,
    hcoe_calicoblast_labels_separate = hcoe_ok,
    values_finite_and_pct_in_range = finite_ok,
    max_formula_difference = formula_difference, formula_ok = formula_ok,
    status = status
  )
  summaries[[i]] <- run
  targets[[i]] <- target
  output_files[[i]] <- data.table(
    object_species = sp,
    role = c("full_long_table", "log2fc_matrix", "resolved_targets", "comparison_summary", "run_record"),
    path = paths, size_bytes = file.info(paths)$size,
    sha256 = vapply(paths, sha256_file, character(1L))
  )
}

checks <- rbindlist(checks)
summaries <- rbindlist(summaries, fill = TRUE)
targets <- rbindlist(targets, fill = TRUE)
output_files <- rbindlist(output_files)
write_tsv(checks, file.path(output_dir, "output_validation.tsv"))
write_tsv(summaries, file.path(output_dir, "all_species.run_summary.tsv"))
write_tsv(targets, file.path(output_dir, "all_species.resolved_celltype_targets.tsv"))
write_tsv(output_files, file.path(output_dir, "output_file_manifest.tsv"))
write_tsv(data.table(
  parameter = c(
    "included_in_main_pipeline", "target_mode", "exclude_target_regex",
    "excluded_labels_remain_in_background", "merge_target_labels",
    "test_use", "slot", "fc_slot", "min_pct", "min_diff_pct",
    "logfc_threshold", "only_pos", "min_cells_feature", "min_cells_group",
    "fold_change_pseudocount"
  ),
  value = as.character(c(
    module$included_in_main_pipeline, module$target_mode, module$exclude_target_regex,
    module$excluded_labels_remain_in_background, module$merge_target_labels,
    module$test_use, module$slot, module$fc_slot, module$min_pct,
    module$min_diff_pct, module$logfc_threshold, module$only_pos,
    module$min_cells_feature, module$min_cells_group, module$fold_change_pseudocount
  ))
), file.path(output_dir, "resolved_analysis_parameters.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
if (any(checks$status != "PASS")) stopf("One or more LOG2FC outputs failed validation")
message(sprintf("Validated %d species, %d target cell types, and %d OG-celltype rows",
                nrow(checks), sum(checks$n_target_celltypes), sum(checks$n_rows)))
