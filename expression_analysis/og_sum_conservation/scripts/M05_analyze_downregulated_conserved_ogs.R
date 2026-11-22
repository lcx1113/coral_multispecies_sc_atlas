#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: M05_analyze_downregulated_conserved_ogs.R <analysis_config.yaml> <comparison_id>", call. = FALSE)
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
module <- cfg$downregulated_og
if (is.null(module)) stopf("Config lacks downregulated_og")
comparison_id_arg <- args[[2L]]
manifest <- read_species_manifest(cfg)
comparison <- read_comparisons(cfg, comparison_id_arg)
comparison <- merge(
  comparison,
  manifest[, .(object_species, orthology_species)],
  by = "object_species", all.x = TRUE, sort = FALSE
)
setorder(comparison, species_order)
species_order <- comparison$object_species
n_species <- nrow(comparison)

thresholds <- data.table(
  comparison_id = comparison_id_arg,
  direction = "downregulated",
  n_species = n_species,
  core_min_species = n_species,
  broad_min_species = ceiling(cfg$conservation$broad_fraction * n_species),
  majority_min_species = ceiling(cfg$conservation$majority_fraction * n_species),
  shared_min_species = cfg$conservation$shared_min_species,
  avg_log2fc_filter = sprintf("< %s", module$avg_log2fc_lt),
  pct_background_filter = sprintf("> %s", module$pct_background_gt),
  p_val_adj_filter = sprintf("< %s", module$p_val_adj_lt)
)

broccoli <- fread(cfg$inputs$broccoli_gene_to_og, na.strings = c("NA", ""))
annotations <- fread(cfg$inputs$og_annotations, na.strings = c("NA", ""))
if (anyDuplicated(annotations$orthogroup)) stopf("OG annotation file has duplicate orthogroup rows")
og_universe <- sort(unique(broccoli$orthogroup))

comparison_dir <- comparison_result_dir(cfg, comparison_id_arg)
source_dir <- file.path(comparison_dir, module$source_marker_subdir)
output_dir <- file.path(comparison_dir, module$output_subdir)
filtered_dir <- file.path(output_dir, "filtered")
conservation_dir <- file.path(output_dir, "conservation")
summary_dir <- file.path(output_dir, "summaries")
invisible(lapply(c(filtered_dir, conservation_dir, summary_dir), ensure_dir))
object_dir <- project_path(cfg, file.path(cfg$outputs$intermediate_dir, "02_og_sum_objects"))

details_list <- list()
filter_summaries <- list()
for (i in seq_len(nrow(comparison))) {
  x <- comparison[i]
  full_file <- file.path(source_dir, sprintf("%s.OG_SUM.all_tested.tsv.gz", x$object_species))
  feature_file <- file.path(object_dir, sprintf("%s.OG_SUM.feature_metadata.tsv.gz", x$object_species))
  if (!all(file.exists(c(full_file, feature_file)))) stopf("Missing full marker or feature metadata for %s", x$object_species)
  full <- fread(full_file, na.strings = c("NA", ""))
  feature_meta <- fread(feature_file, na.strings = c("NA", ""))
  required <- c("orthogroup", "avg_log2FC", "pct_target", "pct_background", "p_val_adj")
  if (length(setdiff(required, names(full)))) stopf("%s full marker table lacks required columns", x$object_species)
  if (anyDuplicated(full$orthogroup)) stopf("Duplicated OGs in %s", full_file)
  if (!all(full$comparison_id == comparison_id_arg) || !all(full$object_species == x$object_species)) {
    stopf("Comparison/species identity mismatch in %s", full_file)
  }

  down <- full[
    avg_log2FC < module$avg_log2fc_lt &
      pct_background > module$pct_background_gt &
      p_val_adj < module$p_val_adj_lt
  ]
  down[, `:=`(
    regulation_direction = "downregulated",
    filter_avg_log2fc = sprintf("< %s", module$avg_log2fc_lt),
    filter_pct_background = sprintf("> %s", module$pct_background_gt),
    filter_p_val_adj = sprintf("< %s", module$p_val_adj_lt)
  )]
  filtered_file <- file.path(filtered_dir, sprintf("%s.OG_SUM.downregulated.tsv", x$object_species))
  write_tsv(down, filtered_file)

  d <- data.table(orthogroup = og_universe)
  d[, `:=`(
    comparison_id = comparison_id_arg,
    direction = "downregulated",
    object_species = x$object_species,
    orthology_species = x$orthology_species,
    present_in_broccoli = orthogroup %in% broccoli[species == x$orthology_species, orthogroup],
    mapped_to_expression = orthogroup %in% feature_meta$orthogroup,
    tested_as_marker = orthogroup %in% full$orthogroup,
    is_downregulated_marker = orthogroup %in% down$orthogroup
  )]
  stat_cols <- intersect(
    c("orthogroup", "p_val", "avg_log2FC", "pct_target", "pct_background", "p_val_adj"),
    names(full)
  )
  d <- merge(d, full[, ..stat_cols], by = "orthogroup", all.x = TRUE, sort = FALSE)
  d[, marker_test_status := fcase(
    is_downregulated_marker, "downregulated_marker",
    tested_as_marker, "tested_not_downregulated",
    mapped_to_expression, "mapped_not_tested",
    present_in_broccoli, "broccoli_only_unmapped",
    default = "not_in_species_broccoli"
  )]
  details_list[[i]] <- d
  filter_summaries[[i]] <- data.table(
    comparison_id = comparison_id_arg,
    object_species = x$object_species,
    orthology_species = x$orthology_species,
    n_full_marker_rows = nrow(full),
    n_downregulated_rows = nrow(down),
    filtered_file = filtered_file
  )
}

details <- rbindlist(details_list, use.names = TRUE)
details[, species_order_index := match(object_species, species_order)]
setorder(details, orthogroup, species_order_index)
details[, species_order_index := NULL]
filter_summary <- rbindlist(filter_summaries)

summary <- details[, .(
  n_species_present_in_broccoli = sum(present_in_broccoli),
  n_species_mapped_to_expression = sum(mapped_to_expression),
  n_species_tested = sum(tested_as_marker),
  n_downregulated_species = sum(is_downregulated_marker),
  downregulated_species = paste(object_species[is_downregulated_marker], collapse = ";")
), by = orthogroup]
summary[, `:=`(
  comparison_id = comparison_id_arg,
  direction = "downregulated",
  n_species_in_comparison = n_species,
  core_min_species = thresholds$core_min_species,
  broad_min_species = thresholds$broad_min_species,
  majority_min_species = thresholds$majority_min_species,
  shared_min_species = thresholds$shared_min_species
)]
summary[, conservation_class := fcase(
  n_downregulated_species >= core_min_species, "core",
  n_downregulated_species >= broad_min_species, "broad_only",
  n_downregulated_species >= majority_min_species, "majority_only",
  n_downregulated_species >= shared_min_species, "shared_only",
  n_downregulated_species == 1L, "species_specific",
  default = "not_downregulated"
)]
summary <- merge(summary, annotations, by = "orthogroup", all.x = TRUE, sort = FALSE)
setorder(summary, -n_downregulated_species, orthogroup)

core <- summary[n_downregulated_species >= core_min_species]
broad <- summary[n_downregulated_species >= broad_min_species]
shared <- summary[n_downregulated_species >= shared_min_species]
prefix <- paste0(comparison_id_arg, ".OG_SUM.downregulated")

threshold_file <- file.path(conservation_dir, paste0(prefix, ".conservation_thresholds.tsv"))
summary_file <- file.path(conservation_dir, paste0(prefix, ".conservation_summary.tsv"))
details_file <- file.path(conservation_dir, paste0(prefix, ".species_details.tsv"))
core_file <- file.path(conservation_dir, paste0(prefix, ".core_candidates.tsv"))
broad_file <- file.path(conservation_dir, paste0(prefix, ".broad_candidates.tsv"))
shared_file <- file.path(conservation_dir, paste0(prefix, ".shared_candidates.tsv"))
matrix_file <- file.path(conservation_dir, paste0(prefix, ".broad_log2fc_matrix.tsv"))

write_tsv(thresholds, threshold_file)
write_tsv(summary, summary_file)
write_tsv(details, details_file)
write_tsv(core, core_file)
write_tsv(broad, broad_file)
write_tsv(shared, shared_file)

matrix_long <- details[orthogroup %in% broad$orthogroup, .(orthogroup, object_species, avg_log2FC)]
if (nrow(matrix_long)) {
  matrix_table <- dcast(matrix_long, orthogroup ~ object_species, value.var = "avg_log2FC")
  missing_species <- setdiff(species_order, names(matrix_table))
  for (sp in missing_species) matrix_table[, (sp) := NA_real_]
  setcolorder(matrix_table, c("orthogroup", species_order))
  matrix_table <- matrix_table[match(broad$orthogroup, orthogroup)]
} else {
  matrix_table <- data.table(orthogroup = character())
  for (sp in species_order) matrix_table[, (sp) := numeric()]
}
write_tsv(matrix_table, matrix_file)
write_tsv(filter_summary, file.path(summary_dir, "filter_summary.tsv"))
write_tsv(data.table(
  comparison_id = comparison_id_arg,
  n_species = n_species,
  n_og_universe = length(og_universe),
  n_downregulated_rows_total = sum(filter_summary$n_downregulated_rows),
  n_shared_ogs = nrow(shared),
  n_broad_including_core_ogs = nrow(broad),
  n_core_ogs = nrow(core),
  source_marker_directory = source_dir,
  conservation_summary = summary_file,
  completed_at_utc = utc_now()
), file.path(summary_dir, "run_summary.tsv"))

message(sprintf("%s downregulated OGs: %d shared, %d broad-including-core, %d core",
                comparison_id_arg, nrow(shared), nrow(broad), nrow(core)))
