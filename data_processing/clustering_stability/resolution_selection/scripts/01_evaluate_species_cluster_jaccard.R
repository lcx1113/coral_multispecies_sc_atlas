#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

default_run_dir <- paste0(".", "/Adig/runs/",
  "20260810_Adig_canonical_sct_leiden4_r04_r50"
)
default_output_dir <- paste0(".", "/resolution_selection"
)

jaccard_path <- if (length(args) >= 1L) args[[1L]] else file.path(
  default_run_dir, "03_resampling", "tables",
  "cluster_jaccard__formal__raw_integer_full_pipeline.csv"
)
stability_path <- if (length(args) >= 2L) args[[2L]] else file.path(
  default_run_dir, "03_resampling", "tables",
  "cluster_stability__formal__raw_integer_full_pipeline.csv"
)
output_dir <- if (length(args) >= 3L) args[[3L]] else default_output_dir
stable_threshold <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 0.6
min_sampled_cells <- if (length(args) >= 5L) as.integer(args[[5L]]) else 10L
species_name <- if (length(args) >= 6L) args[[6L]] else "Adig"
run_dir <- dirname(dirname(dirname(normalizePath(jaccard_path))))
assignment_path <- if (length(args) >= 7L) args[[7L]] else file.path(
  run_dir, "02_full_resolution_scan", "tables",
  "full_resolution_assignments.tsv.gz"
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The R package 'ggplot2' is required.")
}
if (!is.finite(stable_threshold) || stable_threshold < 0 || stable_threshold > 1) {
  stop("stable_threshold must be between 0 and 1.")
}
if (!file.exists(jaccard_path)) stop("Missing Jaccard table: ", jaccard_path)
if (!file.exists(stability_path)) stop("Missing stability table: ", stability_path)
if (!file.exists(assignment_path)) {
  stop("Missing full-resolution assignments: ", assignment_path)
}

screen_table_dir <- file.path(output_dir, "01_resolution_screen", "tables")
screen_figure_dir <- file.path(output_dir, "01_resolution_screen", "figures")
selected_table_dir <- file.path(output_dir, "02_selected_resolution", "tables")
selected_figure_dir <- file.path(output_dir, "02_selected_resolution", "figures")
docs_dir <- file.path(output_dir, "docs")
logs_dir <- file.path(output_dir, "logs")
for (path in c(
  screen_table_dir, screen_figure_dir, selected_table_dir,
  selected_figure_dir, docs_dir, logs_dir
)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_csv <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

required_stability <- c(
  "resolution", "reference_cluster", "n_iterations",
  "n_sufficient_iterations", "median_jaccard", "q10_jaccard",
  "stable_iteration_fraction"
)
required_jaccard <- c(
  "resolution", "iteration", "reference_cluster",
  "reference_cluster_sampled_n", "max_jaccard"
)

cluster_stability <- read_csv(stability_path)
jaccard <- read_csv(jaccard_path)
missing_stability <- setdiff(required_stability, names(cluster_stability))
missing_jaccard <- setdiff(required_jaccard, names(jaccard))
if (length(missing_stability)) {
  stop("Missing columns in stability table: ", paste(missing_stability, collapse = ", "))
}
if (length(missing_jaccard)) {
  stop("Missing columns in Jaccard table: ", paste(missing_jaccard, collapse = ", "))
}

cluster_stability <- cluster_stability[
  is.finite(cluster_stability$resolution) &
    is.finite(cluster_stability$median_jaccard),
  ,
  drop = FALSE
]
cluster_stability$is_stable <- cluster_stability$median_jaccard >= stable_threshold
cluster_stability$reference_cluster <- as.character(
  cluster_stability$reference_cluster
)

resolution_column <- function(value) {
  label <- format(as.numeric(value), scientific = FALSE, trim = TRUE)
  label <- sub("\\.?0+$", "", label)
  paste0("resolution_", gsub("\\.", "_", label))
}

resolution_values <- sort(unique(cluster_stability$resolution))
assignment_columns <- vapply(
  resolution_values, resolution_column, character(1)
)
reference <- utils::read.delim(
  gzfile(assignment_path), check.names = FALSE,
  colClasses = c("cell_id" = "character")
)
missing_assignment_columns <- setdiff(assignment_columns, names(reference))
if (length(missing_assignment_columns)) {
  stop(
    "Missing resolution columns in full assignments: ",
    paste(missing_assignment_columns, collapse = ", ")
  )
}
size_rows <- lapply(seq_along(resolution_values), function(i) {
  sizes <- table(as.character(reference[[assignment_columns[[i]]]]))
  data.frame(
    resolution = resolution_values[[i]],
    reference_cluster = names(sizes),
    full_cluster_size = as.integer(sizes),
    stringsAsFactors = FALSE
  )
})
cluster_sizes <- do.call(rbind, size_rows)
cluster_stability <- merge(
  cluster_stability, cluster_sizes,
  by = c("resolution", "reference_cluster"), all.x = TRUE, sort = FALSE
)
if (anyNA(cluster_stability$full_cluster_size)) {
  stop("Failed to match one or more stability clusters to full-data sizes.")
}

summarize_resolution <- function(x) {
  unstable <- !x$is_stable
  unstable_mean <- if (any(unstable)) {
    mean(x$median_jaccard[unstable], na.rm = TRUE)
  } else {
    NA_real_
  }
  data.frame(
    resolution = x$resolution[[1L]],
    n_clusters = nrow(x),
    stable_cluster_count = sum(x$is_stable),
    unstable_cluster_count = sum(unstable),
    stable_cluster_fraction = mean(x$is_stable),
    stable_cluster_cell_fraction =
      sum(x$full_cluster_size[x$is_stable]) / sum(x$full_cluster_size),
    unstable_cluster_mean_median_jaccard = unstable_mean,
    q10_cluster_median_jaccard = as.numeric(stats::quantile(
      x$median_jaccard, 0.1, names = FALSE, na.rm = TRUE
    )),
    minimum_cluster_median_jaccard = min(x$median_jaccard, na.rm = TRUE),
    mean_cluster_median_jaccard = mean(x$median_jaccard, na.rm = TRUE),
    median_cluster_median_jaccard = stats::median(x$median_jaccard, na.rm = TRUE),
    median_stable_iteration_fraction = stats::median(
      x$stable_iteration_fraction, na.rm = TRUE
    ),
    minimum_stable_iteration_fraction = min(
      x$stable_iteration_fraction, na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

resolution_groups <- split(cluster_stability, cluster_stability$resolution)
resolution_summary <- do.call(rbind, lapply(resolution_groups, summarize_resolution))
row.names(resolution_summary) <- NULL

# A resolution with no unstable clusters wins the secondary criterion by design.
resolution_summary$unstable_tiebreak_score <- ifelse(
  resolution_summary$unstable_cluster_count == 0L,
  1,
  resolution_summary$unstable_cluster_mean_median_jaccard
)
ranking_order <- order(
  -resolution_summary$stable_cluster_fraction,
  -resolution_summary$unstable_tiebreak_score,
  -resolution_summary$q10_cluster_median_jaccard,
  resolution_summary$resolution
)
resolution_summary <- resolution_summary[ranking_order, , drop = FALSE]
resolution_summary$selection_rank <- seq_len(nrow(resolution_summary))
resolution_summary$selected <- resolution_summary$selection_rank == 1L

selected_resolution <- resolution_summary$resolution[
  resolution_summary$selected
][[1L]]
selected_clusters <- cluster_stability[
  cluster_stability$resolution == selected_resolution,
  ,
  drop = FALSE
]
selected_clusters <- selected_clusters[
  order(selected_clusters$median_jaccard, selected_clusters$reference_cluster),
  ,
  drop = FALSE
]
selected_clusters$stability_status <- ifelse(
  selected_clusters$is_stable, "stable", "unstable"
)
selected_warnings <- selected_clusters[
  selected_clusters$q10_jaccard < stable_threshold |
    selected_clusters$stable_iteration_fraction < 0.8,
  ,
  drop = FALSE
]
selected_warnings$warning_reason <- vapply(
  seq_len(nrow(selected_warnings)),
  function(i) {
    reasons <- character()
    if (selected_warnings$q10_jaccard[[i]] < stable_threshold) {
      reasons <- c(reasons, "q10_jaccard_below_0.6")
    }
    if (selected_warnings$stable_iteration_fraction[[i]] < 0.8) {
      reasons <- c(reasons, "stable_iteration_fraction_below_0.8")
    }
    paste(reasons, collapse = ";")
  },
  character(1)
)

selected_jaccard <- jaccard[
  jaccard$resolution == selected_resolution &
    jaccard$reference_cluster_sampled_n >= min_sampled_cells &
    is.finite(jaccard$max_jaccard),
  ,
  drop = FALSE
]
selected_jaccard$stability_status <- ifelse(
  selected_jaccard$reference_cluster %in%
    selected_clusters$reference_cluster[selected_clusters$is_stable],
  "stable",
  "unstable"
)
cluster_order <- as.character(selected_clusters$reference_cluster)
selected_jaccard$reference_cluster_factor <- factor(
  as.character(selected_jaccard$reference_cluster),
  levels = cluster_order
)

selected_record <- resolution_summary[resolution_summary$selected, , drop = FALSE]
selected_record$stable_jaccard_threshold <- stable_threshold
selected_record$min_sampled_cells_per_cluster <- min_sampled_cells
selected_record$primary_selection_rule <- "maximum_stable_cluster_fraction"
selected_record$secondary_selection_rule <- paste0(
  "maximum_mean_of_unstable_cluster_median_jaccard;",
  "no_unstable_clusters_score_1"
)

write_csv(
  resolution_summary,
  file.path(
    screen_table_dir,
    paste0(species_name, "_resolution_cluster_stability_ranking.csv")
  )
)
write_csv(
  cluster_stability[order(
    cluster_stability$resolution,
    cluster_stability$median_jaccard
  ), ],
  file.path(
    screen_table_dir,
    paste0(species_name, "_cluster_stability_all_resolutions.csv")
  )
)
write_csv(
  selected_record,
  file.path(selected_table_dir, paste0(species_name, "_selected_resolution.csv"))
)
write_csv(
  selected_clusters,
  file.path(
    selected_table_dir,
    paste0(species_name, "_selected_resolution_cluster_summary.csv")
  )
)
write_csv(
  selected_warnings,
  file.path(
    selected_table_dir,
    paste0(species_name, "_selected_resolution_cluster_warnings.csv")
  )
)
write_csv(
  selected_jaccard,
  file.path(
    selected_table_dir,
    paste0(species_name, "_selected_resolution_jaccard_by_iteration.csv")
  )
)

resolution_summary$selected_label <- ifelse(
  resolution_summary$selected, "selected", "not selected"
)
ranking_long <- rbind(
  data.frame(
    resolution = resolution_summary$resolution,
    metric = "Stable cluster fraction",
    value = resolution_summary$stable_cluster_fraction,
    selected_label = resolution_summary$selected_label
  ),
  data.frame(
    resolution = resolution_summary$resolution,
    metric = "Cell fraction in stable clusters",
    value = resolution_summary$stable_cluster_cell_fraction,
    selected_label = resolution_summary$selected_label
  ),
  data.frame(
    resolution = resolution_summary$resolution,
    metric = "Mean median Jaccard of unstable clusters",
    value = resolution_summary$unstable_cluster_mean_median_jaccard,
    selected_label = resolution_summary$selected_label
  ),
  data.frame(
    resolution = resolution_summary$resolution,
    metric = "10th percentile of cluster median Jaccard",
    value = resolution_summary$q10_cluster_median_jaccard,
    selected_label = resolution_summary$selected_label
  )
)
ranking_long$metric <- factor(
  ranking_long$metric,
  levels = c(
    "Stable cluster fraction",
    "Cell fraction in stable clusters",
    "Mean median Jaccard of unstable clusters",
    "10th percentile of cluster median Jaccard"
  )
)
no_unstable_annotation <- data.frame(
  resolution = resolution_summary$resolution[
    resolution_summary$unstable_cluster_count == 0L
  ],
  metric = factor(
    rep(
      "Mean median Jaccard of unstable clusters",
      sum(resolution_summary$unstable_cluster_count == 0L)
    ),
    levels = levels(ranking_long$metric)
  ),
  value = rep(0.95, sum(resolution_summary$unstable_cluster_count == 0L)),
  selected_label = rep(
    "selected",
    sum(resolution_summary$unstable_cluster_count == 0L)
  ),
  label = rep(
    "No unstable\nclusters",
    sum(resolution_summary$unstable_cluster_count == 0L)
  )
)

ranking_plot <- ggplot2::ggplot(
  ranking_long,
  ggplot2::aes(
    x = resolution, y = value,
    color = selected_label, group = 1
  )
) +
  ggplot2::geom_line(color = "#707070", linewidth = 0.45, na.rm = TRUE) +
  ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
  ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y") +
  ggplot2::scale_color_manual(
    values = c("not selected" = "#347A78", "selected" = "#C74343")
  ) +
  ggplot2::scale_x_continuous(
    breaks = sort(unique(resolution_summary$resolution))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F0F0F0")
  ) +
  ggplot2::labs(
    title = paste0(
      species_name,
      " resolution screening by cluster-level resampling stability"
    ),
    subtitle = paste0(
      "Primary: stable-cluster fraction\n",
      "Secondary: unstable-cluster mean median Jaccard; tertiary: ",
      "cluster-median Jaccard q10; threshold = ", stable_threshold
    ),
    x = "Resolution",
    y = "Metric value"
  )
if (nrow(no_unstable_annotation)) {
  ranking_plot <- ranking_plot + ggplot2::geom_text(
    data = no_unstable_annotation,
    ggplot2::aes(x = resolution, y = value, label = label),
    color = "#555555", vjust = 0.5, size = 3.4,
    inherit.aes = FALSE
  )
}

boxplot <- ggplot2::ggplot(
  selected_jaccard,
  ggplot2::aes(
    x = reference_cluster_factor,
    y = max_jaccard,
    fill = stability_status
  )
) +
  ggplot2::geom_boxplot(
    width = 0.72, outlier.size = 0.55, outlier.alpha = 0.35,
    linewidth = 0.35
  ) +
  ggplot2::geom_hline(
    yintercept = stable_threshold,
    linetype = "dashed", color = "#C74343", linewidth = 0.65
  ) +
  ggplot2::geom_point(
    data = selected_clusters,
    ggplot2::aes(
      x = factor(as.character(reference_cluster), levels = cluster_order),
      y = median_jaccard
    ),
    inherit.aes = FALSE,
    shape = 23, size = 2.2, fill = "white", color = "#202020"
  ) +
  ggplot2::scale_fill_manual(
    values = c("stable" = "#4C9A91", "unstable" = "#D9785D")
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) +
  ggplot2::labs(
    title = paste0(
      species_name, " cluster Jaccard distributions at resolution ",
      selected_resolution
    ),
    subtitle = paste0(
      "100 formal subsamples; dashed line = ", stable_threshold,
      "; diamonds = cluster medians"
    ),
    x = "Reference cluster (ordered by median Jaccard)",
    y = "Best matched-cluster Jaccard",
    fill = "Median-based status"
  )

save_plot <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot = plot,
    width = width, height = height, units = "in"
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot = plot,
    width = width, height = height, units = "in", dpi = 220
  )
}

save_plot(
  ranking_plot,
  file.path(
    screen_figure_dir,
    paste0(species_name, "_resolution_cluster_stability_ranking")
  ),
  width = 10, height = 11.5
)
save_plot(
  boxplot,
  file.path(
    selected_figure_dir,
    paste0(
      species_name, "_r", selected_resolution, "_cluster_jaccard_boxplot"
    )
  ),
  width = 12, height = 6.5
)

selected_unstable <- selected_clusters[!selected_clusters$is_stable, , drop = FALSE]
report <- c(
  paste0("# ", species_name, " cluster-level Jaccard judgement"),
  "",
  "## Selection rule",
  "",
  paste0(
    "Clusters are stable when their median matched-cluster Jaccard across ",
    "100 formal subsamples is at least ", stable_threshold, "."
  ),
  "",
  "Resolutions are ranked lexicographically:",
  "",
  "1. Highest stable-cluster fraction.",
  "2. Highest mean median Jaccard among unstable clusters.",
  "3. Highest 10th percentile of cluster median Jaccard (fallback only).",
  "4. Lower resolution (final deterministic fallback only).",
  "",
  paste0(
    "A resolution with no unstable clusters receives a secondary tie-break ",
    "score of 1; the observed unstable-cluster mean remains NA."
  ),
  "",
  "## Result",
  "",
  paste0("- Selected resolution: **", selected_resolution, "**"),
  paste0("- Number of reference clusters: ", selected_record$n_clusters),
  paste0(
    "- Stable clusters: ", selected_record$stable_cluster_count, "/",
    selected_record$n_clusters, " (",
    sprintf("%.1f%%", 100 * selected_record$stable_cluster_fraction), ")"
  ),
  paste0(
    "- Full-data cell fraction in stable clusters: ",
    sprintf("%.1f%%", 100 * selected_record$stable_cluster_cell_fraction)
  ),
  paste0("- Unstable clusters: ", selected_record$unstable_cluster_count),
  paste0(
    "- Minimum cluster median Jaccard: ",
    sprintf("%.3f", selected_record$minimum_cluster_median_jaccard)
  ),
  paste0(
    "- 10th percentile of cluster median Jaccard: ",
    sprintf("%.3f", selected_record$q10_cluster_median_jaccard)
  ),
  paste0(
    "- Median stable-iteration fraction: ",
    sprintf("%.3f", selected_record$median_stable_iteration_fraction)
  ),
  paste0(
    "- Clusters with q10 Jaccard < ", stable_threshold,
    " or stable-iteration fraction < 0.8: ", nrow(selected_warnings)
  ),
  "",
  "## Interpretation",
  "",
  paste0(
    "This selection evaluates the reference clusters produced by the current ",
    species_name,
    " partition. It does not establish that every future clustering run ",
    "at the same resolution will produce the same partition."
  )
)
if (nrow(selected_warnings)) {
  report <- c(
    report, "",
    "## Distribution warnings",
    "",
    paste0(
      "The following selected-resolution clusters have a low Jaccard tail or ",
      "pass the 0.6 threshold in fewer than 80% of subsamples:"
    ),
    "",
    paste0(
      "- Cluster ", selected_warnings$reference_cluster,
      ": median = ", sprintf("%.3f", selected_warnings$median_jaccard),
      ", q10 = ", sprintf("%.3f", selected_warnings$q10_jaccard),
      ", stable fraction = ",
      sprintf("%.2f", selected_warnings$stable_iteration_fraction)
    )
  )
}
if (nrow(selected_unstable)) {
  report <- c(
    report, "",
    "Unstable reference clusters:",
    "",
    paste0(
      "- Cluster ", selected_unstable$reference_cluster,
      ": median Jaccard = ",
      sprintf("%.3f", selected_unstable$median_jaccard)
    )
  )
}
writeLines(
  report,
  file.path(
    docs_dir,
    paste0(species_name, "_cluster_jaccard_judgement.md")
  )
)

manifest <- c(
  paste0("run_time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("jaccard_input: ", normalizePath(jaccard_path)),
  paste0("stability_input: ", normalizePath(stability_path)),
  paste0("full_assignment_input: ", normalizePath(assignment_path)),
  paste0("output_dir: ", normalizePath(output_dir)),
  paste0("stable_jaccard_threshold: ", stable_threshold),
  paste0("min_sampled_cells_per_cluster: ", min_sampled_cells),
  paste0("selected_resolution: ", selected_resolution)
  ,paste0("species: ", species_name)
)
writeLines(manifest, file.path(logs_dir, "run_manifest.txt"))
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "sessionInfo.txt"))

message("Selected resolution: ", selected_resolution)
message(
  "Stable clusters: ", selected_record$stable_cluster_count, "/",
  selected_record$n_clusters
)
message("Outputs written to: ", normalizePath(output_dir))
