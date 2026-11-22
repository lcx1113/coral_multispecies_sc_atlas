#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 05_analyze_conserved_ogs.R <analysis_config.yaml> <comparison_id>", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
comparison_id_arg <- args[[2L]]
manifest <- read_species_manifest(cfg)
comparison <- read_comparisons(cfg, comparison_id_arg)
comparison <- merge(comparison, manifest[, .(object_species, orthology_species)], by = "object_species", all.x = TRUE, sort = FALSE)
setorder(comparison, species_order)
species_order <- comparison$object_species
n_species <- nrow(comparison)
thresholds <- data.table(
  comparison_id = comparison_id_arg,
  n_species = n_species,
  core_min_species = n_species,
  broad_min_species = ceiling(cfg$conservation$broad_fraction * n_species),
  majority_min_species = ceiling(cfg$conservation$majority_fraction * n_species),
  shared_min_species = cfg$conservation$shared_min_species
)

broccoli <- fread(cfg$inputs$broccoli_gene_to_og, na.strings = c("NA", ""))
annotations <- fread(cfg$inputs$og_annotations, na.strings = c("NA", ""))
if (anyDuplicated(annotations$orthogroup)) stopf("OG annotation file has duplicate orthogroup rows")
og_universe <- sort(unique(broccoli$orthogroup))

result_dir <- comparison_result_dir(cfg, comparison_id_arg)
object_dir <- project_path(cfg, file.path(cfg$outputs$intermediate_dir, "02_og_sum_objects"))
conservation_dir <- file.path(result_dir, "conservation")
ensure_dir(conservation_dir)
details_list <- list()

for (i in seq_len(nrow(comparison))) {
  x <- comparison[i]
  full_file <- file.path(result_dir, "markers", "full", sprintf("%s.OG_SUM.all_tested.tsv.gz", x$object_species))
  filtered_file <- file.path(result_dir, "markers", "filtered", sprintf("%s.OG_SUM.filtered.tsv", x$object_species))
  feature_file <- file.path(object_dir, sprintf("%s.OG_SUM.feature_metadata.tsv.gz", x$object_species))
  if (!all(file.exists(c(full_file, filtered_file, feature_file)))) stopf("Missing marker or OG feature file for %s", x$object_species)
  full <- fread(full_file)
  filtered <- fread(filtered_file)
  feature_meta <- fread(feature_file)
  if (anyDuplicated(full$orthogroup)) stopf("Duplicated OGs in %s", full_file)
  d <- data.table(orthogroup = og_universe)
  d[, `:=`(
    comparison_id = comparison_id_arg,
    object_species = x$object_species,
    orthology_species = x$orthology_species,
    present_in_broccoli = orthogroup %in% broccoli[species == x$orthology_species, orthogroup],
    mapped_to_expression = orthogroup %in% feature_meta$orthogroup,
    tested_as_marker = orthogroup %in% full$orthogroup,
    is_positive_marker = orthogroup %in% filtered$orthogroup
  )]
  stat_cols <- intersect(c("orthogroup", "p_val", "avg_log2FC", "pct_target", "pct_background", "p_val_adj"), names(full))
  d <- merge(d, full[, ..stat_cols], by = "orthogroup", all.x = TRUE, sort = FALSE)
  d[, marker_test_status := fcase(
    is_positive_marker, "positive_marker",
    tested_as_marker, "tested_not_positive",
    mapped_to_expression, "mapped_not_tested",
    present_in_broccoli, "broccoli_only_unmapped",
    default = "not_in_species_broccoli"
  )]
  details_list[[i]] <- d
}

details <- rbindlist(details_list, use.names = TRUE)
details[, species_order := match(object_species, species_order)]
setorder(details, orthogroup, species_order)
details[, species_order := NULL]

summary <- details[, .(
  n_species_present_in_broccoli = sum(present_in_broccoli),
  n_species_mapped_to_expression = sum(mapped_to_expression),
  n_species_tested = sum(tested_as_marker),
  n_positive_species = sum(is_positive_marker),
  positive_species = paste(species_order[is_positive_marker[match(species_order, object_species)]], collapse = ";")
), by = orthogroup]
summary[, `:=`(
  comparison_id = comparison_id_arg,
  n_species_in_comparison = n_species,
  core_min_species = thresholds$core_min_species,
  broad_min_species = thresholds$broad_min_species,
  majority_min_species = thresholds$majority_min_species,
  shared_min_species = thresholds$shared_min_species
)]
summary[, conservation_class := fcase(
  n_positive_species >= core_min_species, "core",
  n_positive_species >= broad_min_species, "broad_only",
  n_positive_species >= majority_min_species, "majority_only",
  n_positive_species >= shared_min_species, "shared_only",
  n_positive_species == 1L, "species_specific",
  default = "not_positive"
)]
summary <- merge(summary, annotations, by = "orthogroup", all.x = TRUE, sort = FALSE)
setorder(summary, -n_positive_species, orthogroup)

core <- summary[n_positive_species >= core_min_species]
broad <- summary[n_positive_species >= broad_min_species]
shared <- summary[n_positive_species >= shared_min_species]

prefix <- paste0(comparison_id_arg, ".OG_SUM")
write_tsv(thresholds, file.path(conservation_dir, paste0(prefix, ".conservation_thresholds.tsv")))
write_tsv(summary, file.path(conservation_dir, paste0(prefix, ".conservation_summary.tsv")))
write_tsv(details, file.path(conservation_dir, paste0(prefix, ".species_details.tsv")))
write_tsv(core, file.path(conservation_dir, paste0(prefix, ".core_candidates.tsv")))
write_tsv(broad, file.path(conservation_dir, paste0(prefix, ".broad_candidates.tsv")))
write_tsv(shared, file.path(conservation_dir, paste0(prefix, ".shared_candidates.tsv")))
message(sprintf("%s: %d core and %d broad-including-core conserved OGs", comparison_id_arg, nrow(core), nrow(broad)))
