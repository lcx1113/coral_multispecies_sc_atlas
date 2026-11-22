#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: 06_plot_conserved_log2fc_heatmap.R <analysis_config.yaml> <comparison_id>", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))
suppressPackageStartupMessages(library(pheatmap))

cfg <- read_config(args[[1L]])
comparison_id_arg <- args[[2L]]
comparison <- read_comparisons(cfg, comparison_id_arg)
species_order <- comparison$object_species
conservation_dir <- file.path(comparison_result_dir(cfg, comparison_id_arg), "conservation")
prefix <- paste0(comparison_id_arg, ".OG_SUM")
broad_file <- file.path(conservation_dir, paste0(prefix, ".broad_candidates.tsv"))
details_file <- file.path(conservation_dir, paste0(prefix, ".species_details.tsv"))
if (!all(file.exists(c(broad_file, details_file)))) stopf("Run Step 05 before Step 06 for %s", comparison_id_arg)
broad <- fread(broad_file)
details <- fread(details_file)
if (!nrow(broad)) stopf("No broad conserved OGs to plot for %s", comparison_id_arg)

plot_data <- details[orthogroup %in% broad$orthogroup,
                     .(orthogroup, object_species, avg_log2FC)]
wide <- dcast(plot_data, orthogroup ~ object_species, value.var = "avg_log2FC")
missing_species <- setdiff(species_order, names(wide))
for (sp in missing_species) wide[, (sp) := NA_real_]
wide <- wide[match(broad$orthogroup, orthogroup), c("orthogroup", species_order), with = FALSE]
write_tsv(wide, file.path(conservation_dir, paste0(prefix, ".broad_log2fc_matrix.tsv")))

row_annotation <- broad[, c("orthogroup", "conservation_class", intersect(c("has_TF", "first_annotation", "second_annotation"), names(broad))), with = FALSE]
write_tsv(row_annotation, file.path(conservation_dir, paste0(prefix, ".broad_log2fc_row_annotations.tsv")))

mat <- as.matrix(wide[, ..species_order])
storage.mode(mat) <- "double"
rownames(mat) <- wide$orthogroup
annotation_row <- data.frame(conservation_class = broad$conservation_class, row.names = broad$orthogroup)
limit <- max(abs(mat), na.rm = TRUE)
if (!is.finite(limit) || limit == 0) limit <- 1
breaks <- seq(-limit, limit, length.out = 101L)
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100L)
cluster_rows <- isTRUE(cfg$conservation$heatmap_cluster_rows) && nrow(mat) > 1L

draw_heatmap <- function() {
  pheatmap(
    mat, color = colors, breaks = breaks, scale = "none",
    cluster_rows = cluster_rows,
    cluster_cols = isTRUE(cfg$conservation$heatmap_cluster_columns),
    annotation_row = annotation_row, na_col = "grey85",
    border_color = NA, fontsize_row = if (nrow(mat) <= 100) 6 else 3,
    fontsize_col = 9,
    main = sprintf("%s conserved OG avg_log2FC (broad includes core)", comparison_id_arg),
    silent = TRUE
  )$gtable
}

plot_width <- max(8, 2.5 + 0.55 * ncol(mat))
plot_height <- max(6, min(60, 2.5 + 0.11 * nrow(mat)))
pdf(file.path(conservation_dir, paste0(prefix, ".broad_log2fc_heatmap.pdf")), width = plot_width, height = plot_height, useDingbats = FALSE)
grid::grid.newpage(); grid::grid.draw(draw_heatmap()); dev.off()
png(file.path(conservation_dir, paste0(prefix, ".broad_log2fc_heatmap.png")),
    width = plot_width, height = plot_height, units = "in", res = 200)
grid::grid.newpage(); grid::grid.draw(draw_heatmap()); dev.off()
message(sprintf("Saved %s heatmap for %d broad conserved OGs", comparison_id_arg, nrow(mat)))
