#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) stop("Usage: 04_validate_marker_outputs.R <analysis_config.yaml> <comparison_id> [object_species]", call. = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "common.R"))

cfg <- read_config(args[[1L]])
comparison_id_arg <- args[[2L]]
comparison <- read_comparisons(cfg, comparison_id_arg)
if (length(args) == 3L) {
  if (!args[[3L]] %in% comparison$object_species) stopf("%s is not included in comparison %s", args[[3L]], comparison_id_arg)
  comparison <- comparison[object_species == args[[3L]]]
}
result_dir <- comparison_result_dir(cfg, comparison_id_arg)
checks <- list()

for (i in seq_len(nrow(comparison))) {
  sp <- comparison$object_species[i]
  full_file <- file.path(result_dir, "markers", "full", sprintf("%s.OG_SUM.all_tested.tsv.gz", sp))
  filtered_file <- file.path(result_dir, "markers", "filtered", sprintf("%s.OG_SUM.filtered.tsv", sp))
  exists_ok <- file.exists(full_file) && file.exists(filtered_file)
  gzip_ok <- FALSE
  content_ok <- FALSE
  message_text <- NA_character_
  n_full <- n_filtered <- NA_integer_
  if (exists_ok) {
    magic <- readBin(full_file, what = "raw", n = 2L)
    gzip_ok <- identical(as.integer(magic), c(31L, 139L))
    full <- fread(full_file)
    filtered <- fread(filtered_file)
    n_full <- nrow(full); n_filtered <- nrow(filtered)
    required <- c("orthogroup", "avg_log2FC", "pct_target", "p_val_adj", "comparison_id", "object_species", "tested_slot")
    columns_ok <- !length(setdiff(required, names(full))) && !length(setdiff(required, names(filtered)))
    filter_ok <- !nrow(filtered) || all(
      filtered$avg_log2FC > cfg$markers$filtered_avg_log2fc_gt &
        filtered$pct_target > cfg$markers$filtered_pct_target_gt &
        filtered$p_val_adj < cfg$markers$filtered_p_val_adj_lt
    )
    subset_ok <- all(filtered$orthogroup %in% full$orthogroup)
    identity_ok <- all(full$comparison_id == comparison_id_arg) && all(full$object_species == sp) &&
      all(full$tested_slot == "data")
    duplicate_ok <- !anyDuplicated(full$orthogroup) && !anyDuplicated(filtered$orthogroup)
    content_ok <- columns_ok && filter_ok && subset_ok && identity_ok && duplicate_ok
    if (!content_ok) message_text <- "column/filter/subset/identity/duplicate check failed"
  } else {
    message_text <- "one or both marker files are missing"
  }
  checks[[i]] <- data.table(
    comparison_id = comparison_id_arg, object_species = sp,
    full_file_exists = file.exists(full_file), filtered_file_exists = file.exists(filtered_file),
    full_file_is_gzip = gzip_ok, n_full_rows = n_full, n_filtered_rows = n_filtered,
    content_checks_pass = content_ok,
    status = if (exists_ok && gzip_ok && content_ok) "PASS" else "FAIL",
    message = message_text
  )
}
checks <- rbindlist(checks)
out <- file.path(result_dir, "markers", "summaries", "marker_output_validation.tsv")
write_tsv(checks, out)
if (any(checks$status != "PASS")) stopf("Marker output validation failed: %s", out)
message(sprintf("All marker outputs passed validation for %s", comparison_id_arg))
