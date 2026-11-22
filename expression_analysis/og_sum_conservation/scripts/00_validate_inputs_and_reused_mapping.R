#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 00_validate_inputs_and_reused_mapping.R <analysis_config.yaml>", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
manifest <- read_species_manifest(cfg)
comparisons <- read_comparisons(cfg, included_only = FALSE)
mapping <- fread(cfg$inputs$gene_og_mapping, na.strings = c("NA", ""))

required_mapping <- c("object_species", "orthology_species", "gene_seurat", "orthogroup")
if (length(setdiff(required_mapping, names(mapping)))) {
  stopf("Reused mapping lacks columns: %s", paste(setdiff(required_mapping, names(mapping)), collapse = ", "))
}
mapping <- unique(mapping[!is.na(gene_seurat) & !is.na(orthogroup),
                          .(object_species, orthology_species, gene_seurat, orthogroup)])

validation <- vector("list", nrow(manifest))
celltype_counts <- list()

for (i in seq_len(nrow(manifest))) {
  m <- manifest[i]
  message(sprintf("[%d/%d] Validating %s", i, nrow(manifest), m$object_species))
  if (!file.exists(m$seurat_file)) stopf("Missing Seurat object for %s: %s", m$object_species, m$seurat_file)
  seu <- readRDS(m$seurat_file)
  if (!inherits(seu, "Seurat")) stopf("Input for %s is not a Seurat object", m$object_species)
  if (!m$assay %in% Assays(seu)) stopf("%s lacks assay %s", m$object_species, m$assay)
  if (!m$celltype_column %in% colnames(seu[[]])) stopf("%s lacks metadata column %s", m$object_species, m$celltype_column)
  counts <- get_counts(seu, m$assay, m$counts_layer)
  if (anyDuplicated(rownames(counts))) stopf("%s has duplicated RNA feature names", m$object_species)
  if (anyDuplicated(colnames(counts))) stopf("%s has duplicated cell names", m$object_species)
  if (!identical(colnames(counts), rownames(seu[[]]))) stopf("%s counts/meta.data cell order differs", m$object_species)
  if (length(counts@x) && (any(counts@x < 0) || any(abs(counts@x - round(counts@x)) > 1e-8))) {
    stopf("%s RNA counts contain negative or non-integer values", m$object_species)
  }

  map_i <- mapping[object_species == m$object_species]
  missing_features <- setdiff(map_i$gene_seurat, rownames(counts))
  extra_features <- setdiff(rownames(counts), map_i$gene_seurat)
  orthology_values <- unique(map_i$orthology_species)
  orthology_ok <- identical(orthology_values, m$orthology_species)

  comp_i <- comparisons[object_species == m$object_species & include == TRUE]
  observed_labels <- as.character(seu[[m$celltype_column, drop = TRUE]])
  for (j in seq_len(nrow(comp_i))) {
    labels <- split_labels(comp_i$target_labels[j])
    missing_labels <- setdiff(labels, unique(observed_labels))
    celltype_counts[[length(celltype_counts) + 1L]] <- data.table(
      comparison_id = comp_i$comparison_id[j], object_species = m$object_species,
      target_labels = paste(labels, collapse = ";"),
      n_target_cells = sum(observed_labels %in% labels),
      n_background_cells = sum(!observed_labels %in% labels),
      missing_target_labels = if (length(missing_labels)) paste(missing_labels, collapse = ";") else NA_character_
    )
  }

  validation[[i]] <- data.table(
    object_species = m$object_species,
    orthology_species_expected = m$orthology_species,
    orthology_species_mapping = paste(orthology_values, collapse = ";"),
    orthology_species_match = orthology_ok,
    n_cells = ncol(counts), n_rna_features = nrow(counts),
    n_mapping_rows = nrow(map_i), n_mapped_features = uniqueN(map_i$gene_seurat),
    n_mapping_features_missing_from_new_object = length(missing_features),
    mapping_features_missing_fraction = length(missing_features) / uniqueN(map_i$gene_seurat),
    n_new_object_features_without_og_mapping = length(extra_features),
    missing_mapping_features_example = if (length(missing_features)) paste(head(missing_features, 20L), collapse = ";") else NA_character_,
    status = if (!orthology_ok || !nrow(map_i)) "FAIL" else if (length(missing_features)) "PASS_WITH_MISSING_FEATURES" else "PASS"
  )
  rm(seu, counts); gc(FALSE)
}

validation <- rbindlist(validation)
celltype_counts <- rbindlist(celltype_counts)
celltype_counts[, status := fifelse(n_target_cells > 0L & n_background_cells > 0L & is.na(missing_target_labels), "PASS", "FAIL")]

input_paths <- c(
  config = cfg$.config_file,
  species_manifest = project_path(cfg, cfg$inputs$species_manifest),
  celltype_comparisons = project_path(cfg, cfg$inputs$celltype_comparisons),
  gene_og_mapping = cfg$inputs$gene_og_mapping,
  broccoli_gene_to_og = cfg$inputs$broccoli_gene_to_og,
  og_annotations = cfg$inputs$og_annotations
)
input_inventory <- rbindlist(c(
  lapply(names(input_paths), function(nm) data.table(input_type = nm, object_species = NA_character_, path = input_paths[[nm]])),
  list(manifest[, .(input_type = "seurat_object", object_species, path = seurat_file)])
), use.names = TRUE)
if (any(!file.exists(input_inventory$path))) stopf("One or more configured input files do not exist")
input_inventory[, `:=`(
  size_bytes = file.info(path)$size,
  modified_time = format(file.info(path)$mtime, tz = "UTC", usetz = TRUE),
  sha256 = vapply(path, sha256_file, character(1L))
)]

outdir <- project_path(cfg, file.path(cfg$outputs$intermediate_dir, "00_input_validation"))
write_tsv(validation, file.path(outdir, "mapping_reuse_validation.tsv"))
write_tsv(celltype_counts, file.path(outdir, "comparison_target_cell_counts.tsv"))
write_tsv(input_inventory, file.path(outdir, "input_files.tsv"))

if (any(validation$status == "FAIL")) stopf("Reused mapping validation failed; inspect mapping_reuse_validation.tsv")
if (any(celltype_counts$status != "PASS")) stopf("Cell-type target validation failed; inspect comparison_target_cell_counts.tsv")
message("All configured inputs and reused mappings passed validation")
