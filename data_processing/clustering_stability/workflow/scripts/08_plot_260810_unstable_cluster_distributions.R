#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1L]] else paste0(".", ""
)
script_path <- normalizePath(sub(
  "^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]
))
workflow_root <- dirname(dirname(script_path))
species <- c("Adig", "Ahya", "Amur", "Fusp", "Gfas", "Hcoe", "Pdam", "Plut", "Xesp")
tag <- "pilot__raw_integer_full_pipeline"

suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

result_root <- file.path(output_root, "summary", "structural_pilot_unstable_distributions")
table_dir <- file.path(result_root, "tables")
figure_dir <- file.path(result_root, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

iteration_rows <- list()
for (sample in species) {
  config <- yaml::read_yaml(file.path(output_root, sample, "config", "config.yaml"))
  threshold <- as.numeric(config$resampling$stable_jaccard_threshold)
  min_sampled <- as.integer(config$resampling$min_sampled_cells_per_cluster)
  path <- file.path(
    config$run_root, config$run_id, "03_resampling", "tables",
    paste0("cluster_jaccard__", tag, ".csv")
  )
  jaccard <- utils::read.csv(path, stringsAsFactors = FALSE)
  groups <- split(jaccard, interaction(jaccard$resolution, jaccard$iteration, drop = TRUE))
  iteration_rows[[sample]] <- do.call(rbind, lapply(groups, function(x) {
    sufficient <- x$reference_cluster_sampled_n >= min_sampled
    unstable <- sufficient & x$max_jaccard < threshold
    stable <- sufficient & x$max_jaccard >= threshold
    eligible_cells <- sum(x$reference_cluster_sampled_n[sufficient])
    unstable_cells <- sum(x$reference_cluster_sampled_n[unstable])
    data.frame(
      species = sample,
      resolution = x$resolution[[1L]],
      iteration = x$iteration[[1L]],
      seed = x$seed[[1L]],
      reference_cluster_count = nrow(x),
      eligible_reference_cluster_count = sum(sufficient),
      stable_reference_cluster_count = sum(stable),
      unstable_reference_cluster_count = sum(unstable),
      excluded_small_sampled_cluster_count = sum(!sufficient),
      eligible_sampled_cell_count = eligible_cells,
      stable_sampled_cell_count = sum(x$reference_cluster_sampled_n[stable]),
      unstable_sampled_cell_count = unstable_cells,
      unstable_sampled_cell_fraction = if (eligible_cells) {
        unstable_cells / eligible_cells
      } else NA_real_,
      partition_ARI = x$partition_ARI[[1L]],
      partition_NMI = x$partition_NMI[[1L]],
      stringsAsFactors = FALSE
    )
  }))
}
iteration_all <- do.call(rbind, iteration_rows)
row.names(iteration_all) <- NULL

distribution_summary <- do.call(rbind, lapply(
  split(iteration_all, interaction(iteration_all$species, iteration_all$resolution, drop = TRUE)),
  function(x) {
    summarize_metric <- function(values) c(
      median = stats::median(values),
      q10 = as.numeric(stats::quantile(values, 0.1, names = FALSE)),
      q90 = as.numeric(stats::quantile(values, 0.9, names = FALSE)),
      minimum = min(values),
      maximum = max(values)
    )
    cluster_values <- summarize_metric(x$unstable_reference_cluster_count)
    cell_values <- summarize_metric(x$unstable_sampled_cell_count)
    fraction_values <- summarize_metric(x$unstable_sampled_cell_fraction)
    data.frame(
      species = x$species[[1L]],
      resolution = x$resolution[[1L]],
      n_iterations = nrow(x),
      median_unstable_reference_clusters = cluster_values[["median"]],
      q10_unstable_reference_clusters = cluster_values[["q10"]],
      q90_unstable_reference_clusters = cluster_values[["q90"]],
      min_unstable_reference_clusters = cluster_values[["minimum"]],
      max_unstable_reference_clusters = cluster_values[["maximum"]],
      median_unstable_sampled_cells = cell_values[["median"]],
      q10_unstable_sampled_cells = cell_values[["q10"]],
      q90_unstable_sampled_cells = cell_values[["q90"]],
      min_unstable_sampled_cells = cell_values[["minimum"]],
      max_unstable_sampled_cells = cell_values[["maximum"]],
      median_unstable_sampled_cell_fraction = fraction_values[["median"]],
      q10_unstable_sampled_cell_fraction = fraction_values[["q10"]],
      q90_unstable_sampled_cell_fraction = fraction_values[["q90"]],
      min_unstable_sampled_cell_fraction = fraction_values[["minimum"]],
      max_unstable_sampled_cell_fraction = fraction_values[["maximum"]],
      stringsAsFactors = FALSE
    )
  }
))
row.names(distribution_summary) <- NULL

cluster_evidence <- utils::read.csv(
  file.path(output_root, "summary", "structural_pilot_cluster_evidence.csv"),
  stringsAsFactors = FALSE
)
persistent_summary <- do.call(rbind, lapply(
  split(
    cluster_evidence,
    interaction(cluster_evidence$species, cluster_evidence$resolution, drop = TRUE)
  ),
  function(x) {
    unsupported <- x$resample_status == "unsupported"
    insufficient <- x$resample_status == "insufficient"
    total_cells <- sum(x$full_cluster_size, na.rm = TRUE)
    unsupported_cells <- sum(x$full_cluster_size[unsupported], na.rm = TRUE)
    data.frame(
      species = x$species[[1L]],
      resolution = x$resolution[[1L]],
      reference_cluster_count = nrow(x),
      unsupported_reference_cluster_count = sum(unsupported),
      insufficient_reference_cluster_count = sum(insufficient),
      unsupported_full_cell_count = unsupported_cells,
      unsupported_full_cell_fraction = if (total_cells) unsupported_cells / total_cells else NA_real_,
      unsupported_reference_clusters = paste(x$reference_cluster[unsupported], collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
))
row.names(persistent_summary) <- NULL

utils::write.csv(
  iteration_all, file.path(table_dir, "unstable_reference_clusters_by_iteration.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  distribution_summary,
  file.path(table_dir, "unstable_reference_cluster_distribution_summary.csv"),
  row.names = FALSE, na = ""
)
utils::write.csv(
  persistent_summary,
  file.path(table_dir, "unsupported_full_partition_clusters_by_resolution.csv"),
  row.names = FALSE, na = ""
)

resolution_factor <- function(data) {
  factor(data$resolution, levels = sort(unique(data$resolution)))
}
base_theme <- theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey60"),
    strip.text = element_text(face = "bold")
  )

iteration_all$resolution_label <- resolution_factor(iteration_all)
persistent_summary$resolution_label <- resolution_factor(persistent_summary)

plot_box_facets <- function(y_column, y_label, filename, percent = FALSE) {
  plot <- ggplot(iteration_all, aes(x = resolution_label, y = .data[[y_column]])) +
    geom_boxplot(
      width = 0.72, outlier.size = 0.35, linewidth = 0.3,
      fill = "#70A494", color = "#24343B"
    ) +
    facet_wrap(~species, ncol = 3, scales = "free_x") +
    base_theme +
    labs(x = "Resolution", y = y_label)
  if (percent) plot <- plot + scale_y_continuous(labels = scales::label_percent(accuracy = 1))
  ggsave(file.path(figure_dir, paste0(filename, ".pdf")), plot, width = 15, height = 11)
  ggsave(file.path(figure_dir, paste0(filename, ".png")), plot, width = 15, height = 11, dpi = 180)
}

plot_box_facets(
  "unstable_reference_cluster_count",
  "Unstable full-data reference clusters per subsample",
  "all_species_unstable_reference_cluster_count_distribution"
)
plot_box_facets(
  "unstable_sampled_cell_count",
  "Sampled cells in unstable reference clusters",
  "all_species_unstable_sampled_cell_count_distribution"
)
plot_box_facets(
  "unstable_sampled_cell_fraction",
  "Fraction of eligible sampled cells in unstable reference clusters",
  "all_species_unstable_sampled_cell_fraction_distribution",
  percent = TRUE
)

plot_persistent <- function(y_column, y_label, filename, percent = FALSE) {
  plot <- ggplot(
    persistent_summary,
    aes(x = resolution_label, y = .data[[y_column]], group = 1)
  ) +
    geom_line(linewidth = 0.45, color = "#B04A4A") +
    geom_point(size = 1.25, color = "#B04A4A") +
    facet_wrap(~species, ncol = 3, scales = "free_x") +
    base_theme +
    labs(x = "Resolution", y = y_label)
  if (percent) plot <- plot + scale_y_continuous(labels = scales::label_percent(accuracy = 1))
  ggsave(file.path(figure_dir, paste0(filename, ".pdf")), plot, width = 15, height = 11)
  ggsave(file.path(figure_dir, paste0(filename, ".png")), plot, width = 15, height = 11, dpi = 180)
}

plot_persistent(
  "unsupported_reference_cluster_count",
  "Full-data reference clusters with median Jaccard < 0.6",
  "all_species_unsupported_full_partition_cluster_count"
)
plot_persistent(
  "unsupported_full_cell_count",
  "Full-data cells in clusters with median Jaccard < 0.6",
  "all_species_unsupported_full_partition_cell_count"
)
plot_persistent(
  "unsupported_full_cell_fraction",
  "Fraction of full-data cells in clusters with median Jaccard < 0.6",
  "all_species_unsupported_full_partition_cell_fraction",
  percent = TRUE
)

metric_labels <- c(
  unstable_reference_cluster_count = "Unstable reference cluster count",
  unstable_sampled_cell_count = "Sampled cells in unstable reference clusters",
  unstable_sampled_cell_fraction = "Eligible sampled-cell fraction"
)
for (sample in species) {
  x <- iteration_all[iteration_all$species == sample, , drop = FALSE]
  long <- rbind(
    data.frame(
      resolution_label = x$resolution_label,
      metric = names(metric_labels)[[1L]], value = x$unstable_reference_cluster_count
    ),
    data.frame(
      resolution_label = x$resolution_label,
      metric = names(metric_labels)[[2L]], value = x$unstable_sampled_cell_count
    ),
    data.frame(
      resolution_label = x$resolution_label,
      metric = names(metric_labels)[[3L]], value = x$unstable_sampled_cell_fraction
    )
  )
  long$metric <- factor(long$metric, levels = names(metric_labels), labels = metric_labels)
  plot <- ggplot(long, aes(resolution_label, value)) +
    geom_boxplot(
      width = 0.72, outlier.size = 0.4, linewidth = 0.3,
      fill = "#70A494", color = "#24343B"
    ) +
    facet_wrap(~metric, ncol = 1, scales = "free_y") +
    base_theme +
    labs(
      title = paste(sample, "structural Pilot"),
      subtitle = "Unstable: sampled reference cluster n >= 10 and best Jaccard < 0.6",
      x = "Resolution", y = NULL
    )
  ggsave(
    file.path(figure_dir, paste0(sample, "_unstable_cluster_and_cell_distributions.pdf")),
    plot, width = 11, height = 10
  )
}

code_files <- c(script_path, file.path(workflow_root, "R", "metrics.R"))
utils::write.csv(
  data.frame(
    relative_path = substring(code_files, nchar(workflow_root) + 2L),
    md5 = unname(tools::md5sum(code_files)),
    stringsAsFactors = FALSE
  ),
  file.path(result_root, "plot_code_md5.csv"), row.names = FALSE
)
message("Unstable-cluster distribution tables and figures written to: ", result_root)
