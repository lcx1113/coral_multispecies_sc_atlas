#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || args[[1L]] %in% c("-h", "--help")) {
  cat(paste0(
    "Usage:\n",
    "  M08_add_species_marker_stats_and_rank.R ",
    "<subset.tsv> <species.filtered_marker.tsv> <species> <output.tsv>\n\n",
    "Adds species-specific avg_log2FC and p_val_adj columns, then sorts by ",
    "avg_log2FC descending, p_val_adj ascending, and orthogroup ascending.\n",
    "The input subset is never modified.\n"
  ))
  quit(status = if (length(args) == 1L) 0L else 1L)
}

subset_file <- normalizePath(args[[1L]], mustWork = TRUE)
marker_file <- normalizePath(args[[2L]], mustWork = TRUE)
species <- args[[3L]]
output_file <- normalizePath(
  file.path(normalizePath(dirname(args[[4L]]), mustWork = TRUE), basename(args[[4L]])),
  mustWork = FALSE
)

if (file.exists(output_file)) {
  stop("Refusing to overwrite existing output: ", output_file, call. = FALSE)
}

subset <- fread(subset_file, sep = "\t", na.strings = c("NA", ""))
markers <- fread(marker_file, sep = "\t", na.strings = c("NA", ""))

subset_required <- c("orthogroup", "n_positive_species", "positive_species")
marker_required <- c("orthogroup", "object_species", "avg_log2FC", "p_val_adj")
missing_subset <- setdiff(subset_required, names(subset))
missing_marker <- setdiff(marker_required, names(markers))
if (length(missing_subset)) {
  stop("Subset lacks required columns: ", paste(missing_subset, collapse = ", "), call. = FALSE)
}
if (length(missing_marker)) {
  stop("Marker table lacks required columns: ", paste(missing_marker, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(subset$orthogroup)) {
  stop("Subset contains duplicated orthogroups", call. = FALSE)
}
if (anyDuplicated(markers$orthogroup)) {
  stop("Marker table contains duplicated orthogroups", call. = FALSE)
}
marker_species <- unique(markers$object_species)
if (length(marker_species) != 1L || marker_species != species) {
  stop(
    "Marker table object_species does not uniquely equal ", species,
    "; observed: ", paste(marker_species, collapse = ", "),
    call. = FALSE
  )
}
if (any(subset$n_positive_species != 1L | subset$positive_species != species)) {
  stop("Subset is not exclusively the requested species-specific observed set", call. = FALSE)
}

marker_index <- match(subset$orthogroup, markers$orthogroup)
if (anyNA(marker_index)) {
  missing_ogs <- subset$orthogroup[is.na(marker_index)]
  stop(
    "Marker statistics are missing for ", length(missing_ogs), " OG(s): ",
    paste(head(missing_ogs, 10L), collapse = ", "),
    call. = FALSE
  )
}

avg_col <- paste0(species, "_avg_log2FC")
p_col <- paste0(species, "_p_val_adj")
rank_col <- paste0(species, "_marker_rank")
new_columns <- c(avg_col, p_col, rank_col)
if (length(intersect(new_columns, names(subset)))) {
  stop("Subset already contains one or more requested output columns", call. = FALSE)
}

out <- copy(subset)
out[, (avg_col) := as.numeric(markers$avg_log2FC[marker_index])]
out[, (p_col) := as.numeric(markers$p_val_adj[marker_index])]
setorderv(out, c(avg_col, p_col, "orthogroup"), c(-1L, 1L, 1L), na.last = TRUE)
out[, (rank_col) := seq_len(.N)]
setcolorder(
  out,
  c("orthogroup", rank_col, avg_col, p_col, setdiff(names(out), c("orthogroup", rank_col, avg_col, p_col)))
)

if (nrow(out) != nrow(subset) || !setequal(out$orthogroup, subset$orthogroup)) {
  stop("Internal validation failed: OG membership changed", call. = FALSE)
}
if (anyNA(out[[avg_col]]) || anyNA(out[[p_col]])) {
  stop("Internal validation failed: joined statistics contain NA", call. = FALSE)
}

fwrite(out, output_file, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("Rows: %d\n", nrow(out)))
cat(sprintf("Matched marker statistics: %d/%d\n", sum(!is.na(marker_index)), nrow(out)))
cat(sprintf("Sort: %s descending, %s ascending, orthogroup ascending\n", avg_col, p_col))
cat(sprintf("Output: %s\n", output_file))
