#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args)) args[[1L]] else paste0(".", "/resolution_selection"
)
script_path <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
))
species <- c(
  "Adig", "Ahya", "Amur", "Fusp", "Gfas",
  "Hcoe", "Pdam", "Plut", "Xesp", "Pcru"
)

rows <- lapply(species, function(sample) {
  sample_dir <- file.path(analysis_dir, "species", sample)
  ranking_path <- file.path(
    sample_dir, "01_resolution_screen", "tables",
    paste0(sample, "_resolution_cluster_stability_ranking.csv")
  )
  selected_path <- file.path(
    sample_dir, "02_selected_resolution", "tables",
    paste0(sample, "_selected_resolution.csv")
  )
  cluster_path <- file.path(
    sample_dir, "02_selected_resolution", "tables",
    paste0(sample, "_selected_resolution_cluster_summary.csv")
  )
  jaccard_path <- file.path(
    sample_dir, "02_selected_resolution", "tables",
    paste0(sample, "_selected_resolution_jaccard_by_iteration.csv")
  )
  ranking_png <- file.path(
    sample_dir, "01_resolution_screen", "figures",
    paste0(sample, "_resolution_cluster_stability_ranking.png")
  )
  ranking_pdf <- sub("[.]png$", ".pdf", ranking_png)
  selected <- utils::read.csv(selected_path, check.names = FALSE)
  resolution <- selected$resolution[[1L]]
  boxplot_png <- file.path(
    sample_dir, "02_selected_resolution", "figures",
    paste0(sample, "_r", resolution, "_cluster_jaccard_boxplot.png")
  )
  boxplot_pdf <- sub("[.]png$", ".pdf", boxplot_png)

  ranking <- utils::read.csv(ranking_path, check.names = FALSE)
  clusters <- utils::read.csv(cluster_path, check.names = FALSE)
  jaccard <- utils::read.csv(jaccard_path, check.names = FALSE)
  expected_order <- order(
    -ranking$stable_cluster_fraction,
    -ranking$unstable_tiebreak_score,
    -ranking$q10_cluster_median_jaccard,
    ranking$resolution
  )
  figure_paths <- c(ranking_png, ranking_pdf, boxplot_png, boxplot_pdf)
  passed <-
    nrow(selected) == 1L && selected$selection_rank[[1L]] == 1L &&
    sum(ranking$selected) == 1L && identical(expected_order, seq_len(nrow(ranking))) &&
    identical(ranking$resolution[[1L]], resolution) &&
    nrow(clusters) == selected$n_clusters[[1L]] &&
    length(unique(jaccard$iteration)) == 100L &&
    all(file.exists(figure_paths)) && all(file.info(figure_paths)$size > 0)
  data.frame(
    species = sample,
    selected_resolution = resolution,
    n_tested_resolutions = nrow(ranking),
    n_clusters = selected$n_clusters[[1L]],
    stable_cluster_count = selected$stable_cluster_count[[1L]],
    unstable_cluster_count = selected$unstable_cluster_count[[1L]],
    stable_cluster_fraction = selected$stable_cluster_fraction[[1L]],
    stable_cluster_cell_fraction =
      selected$stable_cluster_cell_fraction[[1L]],
    unstable_cluster_mean_median_jaccard =
      selected$unstable_cluster_mean_median_jaccard[[1L]],
    q10_cluster_median_jaccard = selected$q10_cluster_median_jaccard[[1L]],
    selected_jaccard_iterations = length(unique(jaccard$iteration)),
    figures_present = all(file.exists(figure_paths)),
    ranking_recomputed = identical(expected_order, seq_len(nrow(ranking))),
    stable_cell_fraction_present =
      "stable_cluster_cell_fraction" %in% names(ranking) &&
      all(is.finite(ranking$stable_cluster_cell_fraction)),
    validation = if (
      passed && "stable_cluster_cell_fraction" %in% names(ranking) &&
        all(is.finite(ranking$stable_cluster_cell_fraction))
    ) "passed" else "failed",
    stringsAsFactors = FALSE
  )
})
integrity <- do.call(rbind, rows)
row.names(integrity) <- NULL

manifest_dir <- file.path(analysis_dir, "manifests")
docs_dir <- file.path(analysis_dir, "docs")
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(docs_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  integrity, file.path(manifest_dir, "output_integrity.csv"),
  row.names = FALSE, na = ""
)

code_files <- c(
  file.path(analysis_dir, "scripts", "01_evaluate_species_cluster_jaccard.R"),
  file.path(analysis_dir, "scripts", "02_evaluate_all_species.R"),
  script_path
)
utils::write.csv(
  data.frame(
    path = normalizePath(code_files),
    md5 = unname(tools::md5sum(code_files)),
    stringsAsFactors = FALSE
  ),
  file.path(manifest_dir, "analysis_code_md5.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "# 100次subsample的cluster-level resolution判断",
  "",
  "## 方法",
  "",
  "对每个full-data reference cluster，使用100次formal subsample最佳匹配Jaccard的中位数。`median_jaccard >= 0.6`定义为stable cluster。resolution依次按stable cluster fraction、unstable clusters的mean median Jaccard、cluster median Jaccard的第10百分位排序，最后以较低resolution作为确定性兜底。",
  "",
  "## 完整性",
  "",
  paste0(
    "10个物种的排名均已重新计算验证；每个所选resolution均包含100个iteration，",
    "四指标图和cluster Jaccard箱线图的PDF/PNG均存在。"
  ),
  "",
  "详细选择见`03_cross_species_summary/docs/all_species_selection_and_reproduction.md`；逐物种排名、cluster表、警告和图见`species/<sample>/`；完整性表见`manifests/output_integrity.csv`。",
  "",
  "该规则只在本次实际测试的formal resolution集合中选择，并评价当前full-data partition；不代表唯一或全局最优resolution。"
)
writeLines(summary_lines, file.path(docs_dir, "analysis_run_record.md"))

if (any(integrity$validation != "passed")) {
  print(integrity[integrity$validation != "passed", ], row.names = FALSE)
  stop("One or more resolution_selection outputs failed validation.")
}
message("Validated resolution_selection outputs for all 10 species.")
