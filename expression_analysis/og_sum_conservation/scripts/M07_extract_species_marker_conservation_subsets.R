#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || args[[1L]] %in% c("-h", "--help")) {
  cat(paste0(
    "Usage:\n",
    "  M07_extract_species_marker_conservation_subsets.R ",
    "<conservation_summary.with_final_annotation.tsv> <species> <output_dir>\n\n",
    "Writes two non-overlapping subsets while preserving all input columns:\n",
    "  1) species-specific observed markers: n_positive_species == 1\n",
    "  2) shared conserved markers: the species and >=1 other species are positive\n"
  ))
  quit(status = if (length(args) == 1L) 0L else 1L)
}

input_file <- normalizePath(args[[1L]], mustWork = TRUE)
species <- args[[2L]]
output_dir <- normalizePath(args[[3L]], mustWork = TRUE)

x <- fread(input_file, sep = "\t", na.strings = c("NA", ""))
required <- c(
  "orthogroup", "n_positive_species", "positive_species",
  "conservation_class"
)
missing_columns <- setdiff(required, names(x))
if (length(missing_columns)) {
  stop("Input lacks required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(x$orthogroup)) {
  stop("Input contains duplicated orthogroup rows", call. = FALSE)
}

contains_species <- vapply(
  strsplit(fifelse(is.na(x$positive_species), "", x$positive_species), ";", fixed = TRUE),
  function(z) species %chin% z,
  logical(1L)
)

specific <- x[
  contains_species &
    n_positive_species == 1L &
    positive_species == species
]
shared <- x[
  contains_species &
    n_positive_species >= 2L
]

if (nrow(specific) && any(specific$conservation_class != "species_specific")) {
  stop("Specific subset contains a non-species_specific conservation class", call. = FALSE)
}
if (length(intersect(specific$orthogroup, shared$orthogroup))) {
  stop("Specific and shared subsets overlap", call. = FALSE)
}
expected_species_markers <- sum(contains_species)
if (nrow(specific) + nrow(shared) != expected_species_markers) {
  stop(
    "Subsets do not partition all positive markers for ", species,
    ": expected ", expected_species_markers,
    ", observed ", nrow(specific) + nrow(shared),
    call. = FALSE
  )
}

input_suffix <- ".conservation_summary.with_final_annotation.tsv"
if (!endsWith(basename(input_file), input_suffix)) {
  stop("Unexpected input filename suffix: ", basename(input_file), call. = FALSE)
}
prefix <- substr(
  basename(input_file),
  1L,
  nchar(basename(input_file)) - nchar(input_suffix)
)
specific_file <- file.path(
  output_dir,
  paste0(prefix, ".", species, "_specific_observed.tsv")
)
shared_file <- file.path(
  output_dir,
  paste0(prefix, ".", species, "_shared_conserved.tsv")
)
summary_file <- file.path(
  output_dir,
  paste0(prefix, ".", species, "_subset_summary.tsv")
)

for (path in c(specific_file, shared_file, summary_file)) {
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
}

fwrite(specific, specific_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(shared, shared_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(
  data.table(
    subset = c(
      paste0(species, "_specific_observed"),
      paste0(species, "_shared_conserved")
    ),
    definition = c(
      paste0("n_positive_species == 1 && positive_species == '", species, "'"),
      paste0("positive_species contains '", species, "' && n_positive_species >= 2")
    ),
    n_orthogroups = c(nrow(specific), nrow(shared)),
    output_file = c(specific_file, shared_file)
  ),
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(sprintf("%s-specific observed OGs: %d\n", species, nrow(specific)))
cat(sprintf("%s shared conserved OGs: %d\n", species, nrow(shared)))
cat(sprintf("Summary: %s\n", summary_file))
