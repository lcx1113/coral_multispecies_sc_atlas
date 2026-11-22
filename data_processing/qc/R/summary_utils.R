write_step_count_summary <- function(sample_row, sample_dir, by_rank, marker_cutoff) {
  sample_id <- sample_row$sample_id[[1]]
  empty_summary <- file.path(qc_step_dir(sample_dir, "emptydrops"), "emptyDrops_summary.tsv")
  soupx_summary <- file.path(qc_step_dir(sample_dir, "soupx"), "soupx_cache_summary.tsv")
  marker_summary <- file.path(qc_step_dir(sample_dir, "marker_filter"), "marker_first_filter_summary.tsv")
  top_summary <- file.path(qc_step_dir(sample_dir, "top_mito"), "top_mito_filter_summary.tsv")
  df_summary <- file.path(qc_rate_dir(sample_dir, sprintf("doublet_rate_%.2f", getOption("qc_doublet_rate", 0.01))),
                          "summary", "doublet_rate_QC_summary.tsv")

  ed <- if (file.exists(empty_summary)) read_tsv(empty_summary) else data.frame()
  sx <- if (file.exists(soupx_summary)) read_tsv(soupx_summary) else data.frame()
  mk <- if (file.exists(marker_summary)) read_tsv(marker_summary) else data.frame()
  tm <- if (file.exists(top_summary)) read_tsv(top_summary) else data.frame()
  df <- if (file.exists(df_summary)) read_tsv(df_summary) else data.frame()

  steps <- data.frame(
    sample_id = sample_id,
    step = c("raw_matrix", "emptyDrops", "SoupX", "marker_filter", "topUMI_mito",
             "DoubletFinder", "final_recluster"),
    n_remaining = c(
      if (nrow(ed)) ed$n_raw_barcodes else NA,
      if (nrow(ed)) ed$n_emptyDrops_cells else NA,
      if (nrow(sx)) sx$n_cells_after_soupx else NA,
      if (nrow(mk)) mk$n_cells_after_marker_filter else NA,
      if (nrow(tm)) tm$n_cells_after_marker_top_mito_filter else NA,
      if (nrow(df)) df$n_singlets else NA,
      if (nrow(df)) df$n_singlets else NA
    ),
    parameter_summary = c(
      paste0("input_rds=", sample_row$input_rds[[1]]),
      paste0("by.rank=", by_rank),
      if (nrow(sx)) paste0("rho=", sx$soupX_rhoEst) else "rho=NA",
      paste0("marker_cutoff=", marker_cutoff),
      "topUMI/mito after marker filter",
      if (nrow(df)) paste0("doublet_rate=", df$doublet_rate, "; pK=", df$pK) else "doublet_rate=NA",
      "RNA workflow final resolution"
    ),
    stringsAsFactors = FALSE
  )
  steps$n_input <- c(NA, steps$n_remaining[1:6])
  steps$n_removed <- steps$n_input - steps$n_remaining
  steps <- steps[, c("sample_id", "step", "n_input", "n_removed", "n_remaining", "parameter_summary")]
  write_tsv(steps, file.path(qc_step_dir(sample_dir, "run_info"), "QC_step_cell_count_summary.tsv"))
  steps
}

write_global_summaries <- function(output_root) {
  sample_dirs <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[!basename(sample_dirs) %in% c("summary", "99_summary")]
  summary_dir <- qc_global_summary_dir(output_root)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  step_files <- file.path(vapply(sample_dirs, qc_step_dir, character(1), step = "run_info"), "QC_step_cell_count_summary.tsv")
  step_files <- step_files[file.exists(step_files)]
  if (length(step_files) > 0) {
    write_tsv(dplyr::bind_rows(lapply(step_files, read_tsv)),
              file.path(summary_dir, "QC_step_cell_count_summary_all_samples.tsv"))
  }

  rate_files <- list.files(output_root, pattern = "^doublet_rate_QC_summary.tsv$",
                           recursive = TRUE, full.names = TRUE)
  if (length(rate_files) > 0) {
    rate_summary <- dplyr::bind_rows(lapply(rate_files, read_tsv))
    write_tsv(rate_summary, file.path(summary_dir, "doublet_rate_QC_summary_all_samples.tsv"))
    write_tsv(rate_summary[, intersect(c("sample_id", "final_filtered_rds"), colnames(rate_summary)), drop = FALSE],
              file.path(summary_dir, "final_object_paths.tsv"))
  }

  marker_files <- file.path(vapply(sample_dirs, qc_step_dir, character(1), step = "marker_filter"), "marker_first_filter_summary.tsv")
  marker_files <- marker_files[file.exists(marker_files)]
  if (length(marker_files) > 0) {
    write_tsv(dplyr::bind_rows(lapply(marker_files, read_tsv)),
              file.path(summary_dir, "marker_filter_summary_all_samples.tsv"))
  }
  capture.output(sessionInfo(), file = file.path(summary_dir, "sessionInfo.txt"))
}
