#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L) stop("Usage: 02_build_og_sum_objects.R <analysis_config.yaml> [object_species]", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
manifest_all <- read_species_manifest(cfg)
manifest <- copy(manifest_all)
if (length(args) == 2L) {
  if (!args[[2L]] %in% manifest$object_species) stopf("Unknown object_species: %s", args[[2L]])
  manifest <- manifest[object_species == args[[2L]]]
}
mapping <- fread(cfg$inputs$gene_og_mapping, na.strings = c("NA", ""))
broccoli <- fread(cfg$inputs$broccoli_gene_to_og, na.strings = c("NA", ""))
mapping <- unique(mapping[!is.na(gene_seurat) & !is.na(orthogroup),
                          .(object_species, orthology_species, gene_seurat, orthogroup, is_orphan)])

outdir <- project_path(cfg, file.path(cfg$outputs$intermediate_dir, "02_og_sum_objects"))
ensure_dir(outdir)
summaries <- list()

for (i in seq_len(nrow(manifest))) {
  m <- manifest[i]
  message(sprintf("[%d/%d] Building OG_SUM object for %s", i, nrow(manifest), m$object_species))
  seu <- readRDS(m$seurat_file)
  counts <- get_counts(seu, m$assay, m$counts_layer)
  map_i <- unique(mapping[object_species == m$object_species,
                          .(gene_seurat, orthogroup, is_orphan)])
  missing_features <- setdiff(map_i$gene_seurat, rownames(counts))
  if (length(missing_features)) {
    message(sprintf("%s: excluding %d reused-mapping features absent from the current RNA matrix", m$object_species, length(missing_features)))
  }
  map_i <- map_i[gene_seurat %in% rownames(counts)]
  og_ids <- sort(unique(map_i$orthogroup))
  membership <- sparseMatrix(
    i = match(map_i$orthogroup, og_ids),
    j = match(map_i$gene_seurat, rownames(counts)),
    x = rep(1, nrow(map_i)),
    dims = c(length(og_ids), nrow(counts)),
    dimnames = list(og_ids, rownames(counts))
  )
  og_counts <- membership %*% counts
  gene_multiplicity <- map_i[, .(n_orthogroups = uniqueN(orthogroup)), by = gene_seurat]
  mapped_gene_totals <- Matrix::rowSums(counts[gene_multiplicity$gene_seurat, , drop = FALSE])
  expected_og_total <- sum(mapped_gene_totals * gene_multiplicity$n_orthogroups)
  observed_og_total <- sum(og_counts)
  aggregation_total_match <- isTRUE(all.equal(as.numeric(observed_og_total), as.numeric(expected_og_total), tolerance = 0))
  if (!aggregation_total_match) stopf("%s failed the full-UMI multi-OG aggregation total check", m$object_species)
  keep <- Matrix::rowSums(og_counts) > 0
  n_zero <- sum(!keep)
  og_counts <- og_counts[keep, , drop = FALSE]
  if (length(og_counts@x) && any(abs(og_counts@x - round(og_counts@x)) > 1e-8)) {
    stopf("%s OG_SUM matrix is not integer-valued", m$object_species)
  }

  og_obj <- CreateSeuratObject(
    counts = og_counts, assay = "OG_SUM", meta.data = seu[[]],
    min.cells = 0, min.features = 0, project = paste0("OG_SUM_", m$object_species)
  )
  og_obj <- NormalizeData(
    og_obj, assay = "OG_SUM",
    normalization.method = cfg$normalization$method,
    scale.factor = cfg$normalization$scale_factor,
    verbose = FALSE
  )
  check_rows <- head(rownames(og_counts), 20L)
  check_cells <- head(colnames(og_counts), 20L)
  cell_og_totals <- Matrix::colSums(og_counts)
  expected_normalized <- log1p(
    sweep(as.matrix(og_counts[check_rows, check_cells, drop = FALSE]), 2L,
          cell_og_totals[check_cells], "/") * cfg$normalization$scale_factor
  )
  observed_normalized <- as.matrix(SeuratObject::LayerData(
    og_obj, assay = "OG_SUM", layer = "data"
  )[check_rows, check_cells, drop = FALSE])
  normalization_sample_match <- isTRUE(all.equal(observed_normalized, expected_normalized, tolerance = 1e-10, check.attributes = FALSE))
  if (!normalization_sample_match) stopf("%s failed the LogNormalize numerical spot check", m$object_species)

  gene_n_og <- map_i[, .(n_orthogroups = uniqueN(orthogroup)), by = gene_seurat]
  map_meta <- merge(map_i, gene_n_og, by = "gene_seurat", all.x = TRUE)
  og_meta <- map_meta[, .(
    n_mapped_genes = uniqueN(gene_seurat),
    n_multi_og_member_genes = uniqueN(gene_seurat[n_orthogroups > 1L]),
    n_orphan_member_genes = uniqueN(gene_seurat[is_orphan %in% TRUE])
  ), by = orthogroup]
  broccoli_meta <- broccoli[species == m$orthology_species,
                            .(n_broccoli_genes = uniqueN(gene)), by = orthogroup]
  og_meta <- merge(og_meta, broccoli_meta, by = "orthogroup", all.x = TRUE)
  og_meta <- as.data.frame(og_meta[match(rownames(og_obj), orthogroup)])
  rownames(og_meta) <- og_meta$orthogroup
  og_meta$orthogroup <- NULL
  og_obj[["OG_SUM"]] <- AddMetaData(og_obj[["OG_SUM"]], metadata = og_meta)
  og_obj@misc$conserved_og_workflow <- list(
    object_species = m$object_species,
    orthology_species = m$orthology_species,
    aggregation = "binary gene-to-OG membership; full raw UMI assigned to every mapped OG",
    multi_og_gene_weight = 1,
    normalization = cfg$normalization$method,
    scale_factor = cfg$normalization$scale_factor,
    source_seurat = m$seurat_file,
    source_mapping = cfg$inputs$gene_og_mapping,
    built_at_utc = utc_now()
  )

  object_path <- file.path(outdir, sprintf("%s.OG_SUM.seurat.rds", m$object_species))
  saveRDS(og_obj, object_path, compress = TRUE)
  write_tsv(as.data.table(og_meta, keep.rownames = "orthogroup"),
            file.path(outdir, sprintf("%s.OG_SUM.feature_metadata.tsv.gz", m$object_species)), compress = TRUE)
  summaries[[i]] <- data.table(
    object_species = m$object_species, orthology_species = m$orthology_species,
    n_cells = ncol(og_obj), n_rna_features = nrow(counts),
    n_mapped_rna_features = uniqueN(map_i$gene_seurat),
    n_reused_mapping_features_absent_from_current_object = length(missing_features),
    n_gene_og_pairs = nrow(map_i), n_og_features = nrow(og_obj),
    n_mapped_ogs_removed_for_zero_counts = n_zero,
    expected_full_umi_og_total = expected_og_total,
    observed_og_sum_total = observed_og_total,
    aggregation_total_match = aggregation_total_match,
    normalization_sample_match = normalization_sample_match,
    normalization_method = cfg$normalization$method,
    scale_factor = cfg$normalization$scale_factor,
    multi_og_gene_weight = 1,
    output_file = object_path
  )
  write_tsv(summaries[[i]], file.path(outdir, sprintf("%s.OG_SUM.build_summary.tsv", m$object_species)))
  rm(seu, counts, membership, og_counts, og_obj); gc(FALSE)
}

if (nrow(manifest) == nrow(manifest_all)) {
  write_tsv(rbindlist(summaries), file.path(outdir, "og_sum_build_summary.tsv"))
}
message("OG_SUM object construction completed")
