run_emptydrops_step <- function(sample_row, count_matrix, sample_dir, config, by_rank) {
  empty_dir <- qc_step_dir(sample_dir, "emptydrops")
  dir.create(empty_dir, recursive = TRUE, showWarnings = FALSE)
  barcode_file <- file.path(empty_dir, "emptyDrops_cell_barcodes.txt")
  result_file <- file.path(empty_dir, "emptyDrops_results.tsv.gz")
  summary_file <- file.path(empty_dir, "emptyDrops_summary.tsv")

  if (isTRUE(config$output$resume_from_existing) &&
      !isTRUE(config$output$overwrite_existing) &&
      file.exists(barcode_file) && file.exists(summary_file)) {
    log_msg(sample_row$sample_id[[1]], ": reusing existing emptyDrops output.")
    return(readLines(barcode_file))
  }

  total_umi <- Matrix::colSums(count_matrix)
  rank_df <- data.frame(
    barcode = names(total_umi),
    total_UMI = as.numeric(total_umi),
    stringsAsFactors = FALSE
  )
  rank_df <- rank_df[order(rank_df$total_UMI, decreasing = TRUE), , drop = FALSE]
  rank_df$rank <- seq_len(nrow(rank_df))
  write_tsv(rank_df, file.path(empty_dir, "raw_barcode_total_umi_rank.tsv"))

  p_rank <- ggplot2::ggplot(rank_df, ggplot2::aes(rank, total_UMI)) +
    ggplot2::geom_line(linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = by_rank, colour = "firebrick", linetype = "dashed") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = paste0(sample_row$sample_id[[1]], " barcode total UMI rank"),
      subtitle = paste0("emptyDrops by.rank = ", by_rank),
      x = "Barcode rank",
      y = "Total UMI"
    ) +
    ggplot2::theme_classic()
  plot_and_save(p_rank, file.path(empty_dir, "barcode_rank_plot.pdf"), width = 7, height = 5)

  set.seed(config$random_seed)
  log_msg(sample_row$sample_id[[1]], ": running emptyDrops by.rank=", by_rank)
  ed <- DropletUtils::emptyDrops(
    count_matrix,
    by.rank = by_rank,
    retain = Inf,
    lower = config$emptydrops$lower,
    niters = config$emptydrops$niters,
    test.ambient = config$emptydrops$test_ambient
  )
  ed_df <- as.data.frame(ed)
  ed_df$barcode <- rownames(ed_df)
  write_tsv_gz(ed_df, result_file)

  fdr_cutoff <- config$emptydrops$fdr_cutoff
  keep <- rownames(ed)[which(ed$FDR <= fdr_cutoff)]
  keep <- intersect(keep, colnames(count_matrix))
  writeLines(keep, barcode_file)
  write_tsv(data.frame(
    sample_id = sample_row$sample_id[[1]],
    by_rank = by_rank,
    emptydrops_fdr_cutoff = fdr_cutoff,
    emptydrops_lower = config$emptydrops$lower,
    emptydrops_niters = config$emptydrops$niters,
    n_raw_barcodes = ncol(count_matrix),
    n_emptyDrops_cells = length(keep),
    emptyDrops_cell_barcodes = barcode_file,
    emptyDrops_results = result_file
  ), summary_file)
  keep
}
