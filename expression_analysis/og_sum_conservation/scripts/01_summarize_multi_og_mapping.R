#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 01_summarize_multi_og_mapping.R <analysis_config.yaml>", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
manifest <- read_species_manifest(cfg)
mapping <- fread(cfg$inputs$gene_og_mapping, na.strings = c("NA", ""))
required <- c("gene", "object_species", "orthology_species", "gene_seurat", "orthogroup")
if (length(setdiff(required, names(mapping)))) stopf("Mapping lacks required columns")

mapping <- unique(mapping[
  object_species %in% manifest$object_species & !is.na(gene_seurat) & !is.na(orthogroup),
  .(gene, object_species, orthology_species, gene_seurat, orthogroup)
])
gene_stats <- mapping[, .(
  broccoli_genes = paste(sort(unique(gene)), collapse = ";"),
  n_orthogroups = uniqueN(orthogroup),
  orthogroups = paste(sort(unique(orthogroup)), collapse = ";")
), by = .(object_species, orthology_species, gene_seurat)]

summary <- gene_stats[, .(
  n_mapped_genes = .N,
  n_multi_og_genes = sum(n_orthogroups > 1L),
  multi_og_gene_fraction = sum(n_orthogroups > 1L) / .N
), by = .(object_species, orthology_species)]
summary <- merge(manifest[, .(object_species, manifest_order = .I)], summary,
                 by = "object_species", all.x = TRUE, sort = FALSE)
setorder(summary, manifest_order)
summary[, manifest_order := NULL]

details <- gene_stats[n_orthogroups > 1L,
  .(object_species, orthology_species, gene = broccoli_genes,
    seurat_feature = gene_seurat, n_orthogroups, orthogroups)]
setorder(details, object_species, seurat_feature)

outdir <- project_path(cfg, file.path(cfg$outputs$intermediate_dir, "01_mapping_qc"))
write_tsv(summary, file.path(outdir, "multi_og_mapping_summary.tsv"))
write_tsv(details, file.path(outdir, "multi_og_gene_details.tsv"))
message(sprintf("Saved mapping QC for %d species and %d multi-OG genes", nrow(summary), nrow(details)))
