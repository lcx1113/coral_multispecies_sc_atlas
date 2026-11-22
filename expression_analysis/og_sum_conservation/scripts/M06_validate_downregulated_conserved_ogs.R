#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: M06_validate_downregulated_conserved_ogs.R <analysis_config.yaml> <comparison_id>", call. = FALSE)
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
module <- cfg$downregulated_og
comparison_id_arg <- args[[2L]]
comparison <- read_comparisons(cfg, comparison_id_arg)
comparison_dir <- comparison_result_dir(cfg, comparison_id_arg)
output_dir <- file.path(comparison_dir, module$output_subdir)
filtered_dir <- file.path(output_dir, "filtered")
conservation_dir <- file.path(output_dir, "conservation")
summary_dir <- file.path(output_dir, "summaries")
prefix <- paste0(comparison_id_arg, ".OG_SUM.downregulated")

species_checks <- list()
for (i in seq_len(nrow(comparison))) {
  sp <- comparison$object_species[i]
  path <- file.path(filtered_dir, sprintf("%s.OG_SUM.downregulated.tsv", sp))
  if (!file.exists(path)) stopf("Missing downregulated filtered table: %s", path)
  x <- fread(path, na.strings = c("NA", ""))
  threshold_ok <- !nrow(x) || all(
    x$avg_log2FC < module$avg_log2fc_lt &
      x$pct_background > module$pct_background_gt &
      x$p_val_adj < module$p_val_adj_lt
  )
  identity_ok <- !nrow(x) || all(x$comparison_id == comparison_id_arg & x$object_species == sp)
  duplicate_ok <- !anyDuplicated(x$orthogroup)
  species_checks[[i]] <- data.table(
    validation_scope = "species_filtered", comparison_id = comparison_id_arg,
    object_species = sp, n_rows = nrow(x), threshold_ok = threshold_ok,
    identity_ok = identity_ok, duplicate_ok = duplicate_ok,
    status = if (threshold_ok && identity_ok && duplicate_ok) "PASS" else "FAIL"
  )
}

summary_file <- file.path(conservation_dir, paste0(prefix, ".conservation_summary.tsv"))
final_file <- sub("[.]tsv$", ".with_final_annotation.tsv", summary_file)
details_file <- file.path(conservation_dir, paste0(prefix, ".species_details.tsv"))
threshold_file <- file.path(conservation_dir, paste0(prefix, ".conservation_thresholds.tsv"))
matrix_file <- file.path(conservation_dir, paste0(prefix, ".broad_log2fc_matrix.tsv"))
required_files <- c(summary_file, final_file, details_file, threshold_file, matrix_file)
if (!all(file.exists(required_files))) stopf("One or more downregulated conservation outputs are missing")

summary <- fread(summary_file, na.strings = c("NA", ""))
final <- fread(final_file, na.strings = c("NA", ""))
details <- fread(details_file, na.strings = c("NA", ""))
thresholds <- fread(threshold_file, na.strings = c("NA", ""))
matrix_table <- fread(matrix_file, na.strings = c("NA", ""))
original_names <- names(summary)
original_preserved <- nrow(summary) == nrow(final) &&
  identical(original_names, names(final)[seq_along(original_names)]) &&
  isTRUE(all.equal(summary, final[, ..original_names], check.attributes = FALSE))
has_first <- !is.na(final$first_annotation) & nzchar(trimws(final$first_annotation))
has_second <- !is.na(final$second_annotation) & nzchar(trimws(final$second_annotation))
expected_annotation <- fifelse(has_first, final$first_annotation,
                               fifelse(has_second, final$second_annotation, NA_character_))
expected_support <- fifelse(has_first, as.numeric(final$first_support),
                            fifelse(has_second, as.numeric(final$second_support), NA_real_))
annotation_ok <- identical(final$final_annotation, expected_annotation) &&
  isTRUE(all.equal(final$final_annotation_support, expected_support, check.attributes = FALSE))
summary_unique <- !anyDuplicated(summary$orthogroup)
details_grid_ok <- nrow(details) == nrow(summary) * nrow(comparison) &&
  !anyDuplicated(details[, .(orthogroup, object_species)])
broad_n <- summary[n_downregulated_species >= broad_min_species, .N]
matrix_ok <- nrow(matrix_table) == broad_n &&
  setequal(setdiff(names(matrix_table), "orthogroup"), comparison$object_species)
conservation_counts_ok <- all(
  summary$n_downregulated_species ==
    details[, sum(is_downregulated_marker), by = orthogroup][match(summary$orthogroup, orthogroup), V1]
)
overall_status <- if (all(
  original_preserved, annotation_ok, summary_unique, details_grid_ok,
  matrix_ok, conservation_counts_ok
)) "PASS" else "FAIL"
overall <- data.table(
  validation_scope = "comparison_conservation", comparison_id = comparison_id_arg,
  object_species = NA_character_, n_rows = nrow(summary),
  threshold_ok = conservation_counts_ok, identity_ok = original_preserved && annotation_ok,
  duplicate_ok = summary_unique && details_grid_ok && matrix_ok,
  status = overall_status
)
checks <- rbindlist(c(species_checks, list(overall)), fill = TRUE)
write_tsv(checks, file.path(summary_dir, "output_validation.tsv"))
if (any(checks$status != "PASS")) stopf("Downregulated OG output validation failed for %s", comparison_id_arg)
message(sprintf("Validated downregulated OG outputs for %s", comparison_id_arg))
