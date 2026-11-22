#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
args <- args[args != "--overwrite"]
if (length(args) == 1L && args[[1L]] %in% c("-h", "--help")) {
  cat(paste0(
    "Usage:\n",
    "  M02_add_final_annotation_to_conservation_summary.R ",
    "<input.conservation_summary.tsv> [output.tsv] [--overwrite]\n\n",
    "If output.tsv is omitted, '.with_final_annotation.tsv' is appended ",
    "before the input .tsv suffix. The input file is never modified.\n"
  ))
  quit(status = 0L)
}
if (length(args) < 1L || length(args) > 2L) {
  stop(
    paste(
      "Usage: M02_add_final_annotation_to_conservation_summary.R",
      "<input.conservation_summary.tsv> [output.tsv] [--overwrite]"
    ),
    call. = FALSE
  )
}

input_file <- normalizePath(args[[1L]], mustWork = TRUE)
if (!grepl("[.]conservation_summary[.]tsv$", input_file)) {
  stop("Input must be an uncompressed *.conservation_summary.tsv file", call. = FALSE)
}
output_file <- if (length(args) == 2L) {
  args[[2L]]
} else {
  sub("[.]tsv$", ".with_final_annotation.tsv", input_file)
}
output_file <- normalizePath(
  file.path(normalizePath(dirname(output_file), mustWork = TRUE), basename(output_file)),
  mustWork = FALSE
)
if (identical(input_file, output_file)) {
  stop("Output path must differ from input path", call. = FALSE)
}
if (file.exists(output_file) && !overwrite) {
  stop("Output already exists; choose a new output path: ", output_file, call. = FALSE)
}

x <- fread(input_file, sep = "\t", na.strings = c("NA", ""))
required <- c(
  "orthogroup", "first_annotation", "first_support",
  "second_annotation", "second_support"
)
missing_columns <- setdiff(required, names(x))
if (length(missing_columns)) {
  stop("Input lacks required columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(x$orthogroup)) {
  stop("Input contains duplicated orthogroup rows", call. = FALSE)
}

original_names <- names(x)
original_nrow <- nrow(x)
original_orthogroups <- x$orthogroup
has_first <- !is.na(x$first_annotation) & nzchar(trimws(x$first_annotation))
has_second <- !is.na(x$second_annotation) & nzchar(trimws(x$second_annotation))

x[, final_annotation := fifelse(
  has_first,
  first_annotation,
  fifelse(has_second, second_annotation, NA_character_)
)]
x[, final_annotation_support := fifelse(
  has_first,
  as.numeric(first_support),
  fifelse(has_second, as.numeric(second_support), NA_real_)
)]

if (nrow(x) != original_nrow || !identical(x$orthogroup, original_orthogroups)) {
  stop("Internal validation failed: row count or OG order changed", call. = FALSE)
}
if (!identical(names(x)[seq_along(original_names)], original_names)) {
  stop("Internal validation failed: original column order changed", call. = FALSE)
}

temporary_file <- paste0(output_file, ".tmp.", Sys.getpid())
fwrite(x, temporary_file, sep = "\t", quote = FALSE, na = "NA")
if (file.exists(output_file)) {
  if (!file.copy(temporary_file, output_file, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not overwrite output: ", output_file, call. = FALSE)
  }
  unlink(temporary_file)
} else if (!file.rename(temporary_file, output_file)) {
  stop("Could not promote temporary output to: ", output_file, call. = FALSE)
}

cat(sprintf("Input rows: %d\n", nrow(x)))
cat(sprintf("final_annotation from first_annotation: %d\n", sum(has_first)))
cat(sprintf("final_annotation from second_annotation: %d\n", sum(!has_first & has_second)))
cat(sprintf("final_annotation missing: %d\n", sum(!has_first & !has_second)))
cat(sprintf("Output: %s\n", output_file))
