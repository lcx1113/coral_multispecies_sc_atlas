suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
  library(mclust)
  library(aricode)
})

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
workflow_root <- dirname(dirname(script_path))
source(file.path(workflow_root, "R", "common.R"))
source(file.path(workflow_root, "R", "metrics.R"))

cli <- parse_cli(
  file.path(workflow_root, "config", "config.yaml"),
  extra = list(
    mode = "pilot",
    count_mode = NULL,
    resolutions = NULL,
    freeze_resolution = NULL
  )
)
config <- load_config(cli$config)
initialize_run(config, cli$config)
require_stage(config, "02_full_resolution_scan")
require_stage(config, "03_resampling")
dirs <- ensure_stage_dirs(config, "08_stability_range", c("tables", "figures", "report", "logs"))

mode <- match.arg(cli$mode, c("pilot", "formal"))
count_mode <- cli$count_mode %||% config$count_input$primary_mode
resolutions <- if (!is.null(cli$resolutions)) {
  as_number_vector(cli$resolutions)
} else if (mode == "pilot") {
  config$resampling$pilot_resolutions
} else {
  config$resampling$formal_resolutions
}
resolutions <- validate_resolution_subset(resolutions, config$clustering$resolution_grid, "--resolutions")
if (!length(resolutions)) abort("No resolutions to summarize.")
tag <- paste(mode, count_mode, sep = "__")

log_path <- file.path(dirs$logs, paste0("summarize_stability_range__", tag, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Started:", format(Sys.time()), "\n")

read_required <- function(stage, leaf, filename) {
  path <- file.path(stage_dir(config, stage, leaf), filename)
  if (!file.exists(path)) abort("Required input is missing: %s", path)
  read_table_auto(path)
}
read_optional <- function(stage, leaf, filename) {
  path <- file.path(stage_dir(config, stage, leaf), filename)
  if (!file.exists(path)) return(NULL)
  read_table_auto(path)
}
merge_resolution <- function(x, y) {
  if (is.null(y) || !nrow(y)) return(x)
  merge(x, y, by = "resolution", all.x = TRUE, sort = FALSE)
}

scan_size <- read_required("02_full_resolution_scan", "tables", "resolution_cluster_size_summary.csv")
resample <- read_required("03_resampling", "tables", paste0("resolution_summary__", tag, ".csv"))
cluster_stability <- read_required("03_resampling", "tables", paste0("cluster_stability__", tag, ".csv"))
pair <- read_optional("04_pair_consensus", "tables", paste0("pair_consensus_resolution_summary__", tag, ".csv"))
pair_cluster <- read_optional("04_pair_consensus", "tables", paste0("pair_consensus_cluster_summary__", tag, ".csv"))
sensitivity <- read_optional("05_parameter_sensitivity", "tables", paste0("parameter_sensitivity_summary__", tag, ".csv"))
sensitivity_cluster_runs <- read_optional(
  "05_parameter_sensitivity", "tables",
  paste0("parameter_sensitivity_cluster_jaccard__", tag, ".csv")
)
marker <- read_optional(
  "06_marker_program",
  file.path("fixed_expression", "tables"),
  paste0("marker_recovery_summary__", tag, ".csv")
)
marker_cluster <- read_optional(
  "06_marker_program",
  file.path("fixed_expression", "tables"),
  paste0("marker_recovery_cluster_summary__", tag, ".csv")
)
assignments_path <- file.path(
  stage_dir(config, "02_full_resolution_scan", "tables"),
  "full_resolution_assignments.tsv.gz"
)
assignments <- read_table_auto(assignments_path)
scan_size <- scan_size[scan_size$resolution %in% resolutions, , drop = FALSE]
cluster_stability <- cluster_stability[cluster_stability$resolution %in% resolutions, , drop = FALSE]
if (!is.null(sensitivity) && nrow(sensitivity)) {
  sensitivity_collision_columns <- intersect(
    setdiff(names(sensitivity), "resolution"),
    names(resample)
  )
  names(sensitivity)[match(sensitivity_collision_columns, names(sensitivity))] <-
    paste0("sensitivity_", sensitivity_collision_columns)
}

evidence <- merge_resolution(scan_size, resample)
evidence <- merge_resolution(evidence, pair)
evidence <- merge_resolution(evidence, sensitivity)
evidence <- merge_resolution(evidence, marker)
evidence <- evidence[match(resolutions, evidence$resolution), , drop = FALSE]
optional_metrics <- c(
  "median_consensus_margin",
  "median_reference_ARI",
  "median_marker_effect_recovery",
  "median_effect_rank_correlation"
)
for (column in setdiff(optional_metrics, names(evidence))) evidence[[column]] <- NA_real_

threshold <- config$stability_range
evidence$size_support <- evidence$small_cluster_fraction <= threshold$maximum_small_cluster_fraction &
  evidence$n_below_min == 0
evidence$resampling_support <- evidence$median_cluster_jaccard >= threshold$minimum_median_cluster_jaccard &
  evidence$median_stable_cell_fraction >= threshold$minimum_median_stable_cell_fraction
evidence$sensitivity_available <- !is.na(evidence$median_reference_ARI)
evidence$sensitivity_support <- evidence$sensitivity_available &
  evidence$median_reference_ARI >= threshold$minimum_sensitivity_ari
evidence$pair_consensus_available <- !is.na(evidence$median_consensus_margin)
evidence$pair_consensus_support <- evidence$pair_consensus_available &
  evidence$median_consensus_margin >= threshold$minimum_pair_consensus_margin
evidence$core_evidence_complete <- evidence$sensitivity_available | !isTRUE(threshold$require_sensitivity)
evidence$core_evidence_complete <- evidence$core_evidence_complete &
  (evidence$pair_consensus_available | !isTRUE(threshold$require_pair_consensus))
evidence$core_supported <- evidence$size_support & evidence$resampling_support &
  (evidence$sensitivity_support | !isTRUE(threshold$require_sensitivity)) &
  (evidence$pair_consensus_support | !isTRUE(threshold$require_pair_consensus)) &
  evidence$core_evidence_complete
evidence$evidence_status <- ifelse(
  !evidence$core_evidence_complete,
  "insufficient",
  ifelse(evidence$core_supported, "supported", "not_supported")
)
evidence$support_count <- rowSums(cbind(
  evidence$size_support,
  evidence$resampling_support,
  evidence$sensitivity_support,
  evidence$pair_consensus_support
), na.rm = TRUE)

cluster_evidence <- cluster_stability
if (!is.null(sensitivity_cluster_runs) && nrow(sensitivity_cluster_runs)) {
  sensitivity_cluster_runs <- sensitivity_cluster_runs[
    sensitivity_cluster_runs$resolution %in% resolutions, , drop = FALSE
  ]
  sensitivity_cluster <- bind_rows_base(lapply(
    split(
      sensitivity_cluster_runs,
      interaction(
        sensitivity_cluster_runs$resolution,
        sensitivity_cluster_runs$reference_cluster,
        drop = TRUE
      )
    ),
    function(x) {
      data.frame(
        resolution = x$resolution[[1L]],
        reference_cluster = x$reference_cluster[[1L]],
        sensitivity_n_runs = nrow(x),
        sensitivity_median_jaccard = stats::median(x$max_jaccard, na.rm = TRUE),
        sensitivity_q10_jaccard = as.numeric(stats::quantile(
          x$max_jaccard, 0.1, na.rm = TRUE, names = FALSE
        )),
        sensitivity_min_jaccard = min(x$max_jaccard, na.rm = TRUE),
        sensitivity_stable_run_fraction = mean(
          x$max_jaccard >= config$resampling$stable_jaccard_threshold,
          na.rm = TRUE
        ),
        sensitivity_median_reference_recall = stats::median(
          x$reference_recall, na.rm = TRUE
        ),
        sensitivity_median_test_precision = stats::median(
          x$test_precision, na.rm = TRUE
        ),
        sensitivity_median_match_margin = stats::median(
          x$best_second_jaccard_gap, na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  ))
  cluster_evidence <- merge(
    cluster_evidence, sensitivity_cluster,
    by = c("resolution", "reference_cluster"), all.x = TRUE, sort = FALSE
  )
}
if (!is.null(pair_cluster) && nrow(pair_cluster)) {
  cluster_evidence <- merge(
    cluster_evidence, pair_cluster,
    by = c("resolution", "reference_cluster"), all.x = TRUE, sort = FALSE
  )
}
if (!is.null(marker_cluster) && nrow(marker_cluster)) {
  cluster_evidence <- merge(
    cluster_evidence, marker_cluster,
    by = c("resolution", "reference_cluster"), all.x = TRUE, sort = FALSE,
    suffixes = c("", "_marker")
  )
}

cluster_jaccard_threshold <- config$resampling$stable_jaccard_threshold
pair_margin_threshold <- threshold$minimum_pair_consensus_margin
cluster_evidence$resample_cluster_support <-
  cluster_evidence$median_jaccard >= cluster_jaccard_threshold
cluster_evidence$resample_low_precision_merge_risk <-
  cluster_evidence$median_reference_recall >= cluster_jaccard_threshold &
  cluster_evidence$median_test_precision < 0.5

pair_cluster_columns <- c("consensus_margin", "max_between_leakage")
for (column in setdiff(pair_cluster_columns, names(cluster_evidence))) {
  cluster_evidence[[column]] <- NA_real_
}
cluster_evidence$pair_cluster_available <- !is.na(cluster_evidence$consensus_margin)
cluster_evidence$pair_cluster_support <-
  cluster_evidence$pair_cluster_available &
  cluster_evidence$consensus_margin >= pair_margin_threshold
cluster_evidence$pair_leakage_risk <-
  cluster_evidence$pair_cluster_available &
  cluster_evidence$max_between_leakage > pair_margin_threshold

sensitivity_cluster_columns <- c(
  "sensitivity_median_jaccard",
  "sensitivity_q10_jaccard",
  "sensitivity_median_reference_recall",
  "sensitivity_median_test_precision"
)
for (column in setdiff(sensitivity_cluster_columns, names(cluster_evidence))) {
  cluster_evidence[[column]] <- NA_real_
}
cluster_evidence$sensitivity_cluster_available <-
  !is.na(cluster_evidence$sensitivity_median_jaccard)
cluster_evidence$sensitivity_cluster_support <-
  cluster_evidence$sensitivity_cluster_available &
  cluster_evidence$sensitivity_median_jaccard >= cluster_jaccard_threshold
cluster_evidence$sensitivity_low_precision_merge_risk <-
  cluster_evidence$sensitivity_cluster_available &
  cluster_evidence$sensitivity_median_reference_recall >= cluster_jaccard_threshold &
  cluster_evidence$sensitivity_median_test_precision < 0.5

rank_within_resolution <- function(values, resolution, decreasing = FALSE) {
  ranked_values <- if (decreasing) -values else values
  ave(
    ranked_values, resolution,
    FUN = function(x) rank(x, ties.method = "min", na.last = "keep")
  )
}
cluster_evidence$resample_jaccard_bad_rank <- rank_within_resolution(
  cluster_evidence$median_jaccard, cluster_evidence$resolution
)
cluster_evidence$resample_q10_bad_rank <- rank_within_resolution(
  cluster_evidence$q10_jaccard, cluster_evidence$resolution
)
cluster_evidence$pair_margin_bad_rank <- rank_within_resolution(
  cluster_evidence$consensus_margin, cluster_evidence$resolution
)
cluster_evidence$pair_leakage_bad_rank <- rank_within_resolution(
  cluster_evidence$max_between_leakage, cluster_evidence$resolution, decreasing = TRUE
)
cluster_evidence$sensitivity_jaccard_bad_rank <- rank_within_resolution(
  cluster_evidence$sensitivity_median_jaccard, cluster_evidence$resolution
)
cluster_evidence$sensitivity_q10_bad_rank <- rank_within_resolution(
  cluster_evidence$sensitivity_q10_jaccard, cluster_evidence$resolution
)

cluster_evidence$is_bottom3_resample <-
  !is.na(cluster_evidence$resample_jaccard_bad_rank) &
  cluster_evidence$resample_jaccard_bad_rank <= 3
cluster_evidence$is_bottom3_pair_margin <-
  !is.na(cluster_evidence$pair_margin_bad_rank) &
  cluster_evidence$pair_margin_bad_rank <= 3
cluster_evidence$is_top3_pair_leakage <-
  !is.na(cluster_evidence$pair_leakage_bad_rank) &
  cluster_evidence$pair_leakage_bad_rank <= 3
cluster_evidence$is_bottom3_sensitivity <-
  !is.na(cluster_evidence$sensitivity_jaccard_bad_rank) &
  cluster_evidence$sensitivity_jaccard_bad_rank <= 3

cluster_evidence$failed_evidence_count <-
  as.integer(!cluster_evidence$resample_cluster_support) +
  as.integer(cluster_evidence$pair_cluster_available & !cluster_evidence$pair_cluster_support) +
  as.integer(
    cluster_evidence$sensitivity_cluster_available &
      !cluster_evidence$sensitivity_cluster_support
  )
cluster_evidence$cluster_evidence_status <- ifelse(
  !cluster_evidence$pair_cluster_available |
    !cluster_evidence$sensitivity_cluster_available,
  "insufficient",
  ifelse(cluster_evidence$failed_evidence_count == 0L, "supported", "fragile")
)
cluster_evidence$fragile_reasons <- vapply(seq_len(nrow(cluster_evidence)), function(i) {
  reasons <- character()
  if (!isTRUE(cluster_evidence$resample_cluster_support[[i]])) {
    reasons <- c(reasons, "resample_jaccard_below_threshold")
  }
  if (isTRUE(cluster_evidence$resample_low_precision_merge_risk[[i]])) {
    reasons <- c(reasons, "resample_merge_risk")
  }
  if (
    isTRUE(cluster_evidence$pair_cluster_available[[i]]) &&
      !isTRUE(cluster_evidence$pair_cluster_support[[i]])
  ) {
    reasons <- c(reasons, "pair_margin_below_threshold")
  }
  if (isTRUE(cluster_evidence$pair_leakage_risk[[i]])) {
    reasons <- c(reasons, "pair_leakage_above_threshold")
  }
  if (
    isTRUE(cluster_evidence$sensitivity_cluster_available[[i]]) &&
      !isTRUE(cluster_evidence$sensitivity_cluster_support[[i]])
  ) {
    reasons <- c(reasons, "sensitivity_jaccard_below_threshold")
  }
  if (isTRUE(cluster_evidence$sensitivity_low_precision_merge_risk[[i]])) {
    reasons <- c(reasons, "sensitivity_merge_risk")
  }
  paste(unique(reasons), collapse = ";")
}, character(1))

driver_selector <- cluster_evidence$failed_evidence_count > 0L |
  cluster_evidence$resample_low_precision_merge_risk |
  cluster_evidence$pair_leakage_risk |
  cluster_evidence$sensitivity_low_precision_merge_risk |
  cluster_evidence$is_bottom3_resample |
  cluster_evidence$is_bottom3_pair_margin |
  cluster_evidence$is_top3_pair_leakage |
  cluster_evidence$is_bottom3_sensitivity
fragile_drivers <- cluster_evidence[driver_selector, , drop = FALSE]
fragile_drivers <- fragile_drivers[
  order(fragile_drivers$resolution, -fragile_drivers$failed_evidence_count),
  , drop = FALSE
]

collapse_clusters <- function(x, metric, decreasing = FALSE, n = 3L) {
  valid <- !is.na(x[[metric]])
  x <- x[valid, , drop = FALSE]
  if (!nrow(x)) return("")
  direction <- if (decreasing) -1 else 1
  x <- x[order(direction * x[[metric]], as.numeric(x$reference_cluster)), , drop = FALSE]
  x <- head(x, n)
  paste0(x$reference_cluster, " (", format(round(x[[metric]], 3), nsmall = 3), ")", collapse = "; ")
}
driver_summary <- bind_rows_base(lapply(split(cluster_evidence, cluster_evidence$resolution), function(x) {
  cluster_size <- if ("cluster_size" %in% names(x)) x$cluster_size else rep(NA_real_, nrow(x))
  fragile <- x$failed_evidence_count > 0L
  data.frame(
    resolution = x$resolution[[1L]],
    n_clusters = nrow(x),
    fragile_cluster_count = sum(fragile),
    fragile_cluster_fraction = mean(fragile),
    fragile_cell_fraction = if (all(is.na(cluster_size))) {
      NA_real_
    } else {
      sum(cluster_size[fragile], na.rm = TRUE) / sum(cluster_size, na.rm = TRUE)
    },
    resample_fail_count = sum(!x$resample_cluster_support),
    pair_fail_count = sum(x$pair_cluster_available & !x$pair_cluster_support),
    sensitivity_fail_count = sum(
      x$sensitivity_cluster_available & !x$sensitivity_cluster_support
    ),
    resample_merge_risk_count = sum(x$resample_low_precision_merge_risk),
    pair_leakage_risk_count = sum(x$pair_leakage_risk),
    sensitivity_merge_risk_count = sum(x$sensitivity_low_precision_merge_risk),
    worst_resample_clusters = collapse_clusters(x, "median_jaccard"),
    worst_pair_margin_clusters = collapse_clusters(x, "consensus_margin"),
    highest_pair_leakage_clusters = collapse_clusters(
      x, "max_between_leakage", decreasing = TRUE
    ),
    worst_sensitivity_clusters = collapse_clusters(x, "sensitivity_median_jaccard"),
    stringsAsFactors = FALSE
  )
}))
driver_summary <- driver_summary[
  match(resolutions, driver_summary$resolution), , drop = FALSE
]

candidate_adjacent <- adjacent_partition_metrics(assignments, resolutions)
candidate_adjacent$transition_support <- candidate_adjacent$adjacent_ARI >= threshold$minimum_adjacent_ari
incoming_transition <- rep(TRUE, length(resolutions))
if (length(resolutions) > 1L) incoming_transition[-1L] <- candidate_adjacent$transition_support

interval_group <- integer(length(resolutions))
group <- 0L
for (i in seq_along(resolutions)) {
  if (!isTRUE(evidence$core_supported[[i]])) {
    interval_group[[i]] <- 0L
    next
  }
  if (i == 1L || !isTRUE(evidence$core_supported[[i - 1L]]) || !isTRUE(incoming_transition[[i]])) {
    group <- group + 1L
  }
  interval_group[[i]] <- group
}
stable_indices <- which(interval_group > 0L)
intervals <- if (!length(stable_indices)) {
  data.frame(
    interval_id = character(), resolution_min = numeric(), resolution_max = numeric(),
    n_evaluated_points = integer(), resolutions = character(), status = character()
  )
} else {
  bind_rows_base(lapply(split(stable_indices, interval_group[stable_indices]), function(indices) {
    values <- resolutions[indices]
    data.frame(
      interval_id = paste0("stable_", interval_group[[indices[[1L]]]]),
      resolution_min = min(values),
      resolution_max = max(values),
      n_evaluated_points = length(values),
      resolutions = paste(vapply(values, fmt_resolution, character(1)), collapse = ","),
      status = if (mode == "formal") "formal_supported" else "pilot_provisional",
      stringsAsFactors = FALSE
    )
  }))
}
representatives <- if (!nrow(intervals)) {
  data.frame()
} else {
  bind_rows_base(lapply(seq_len(nrow(intervals)), function(i) {
    values <- as_number_vector(intervals$resolutions[[i]])
    midpoint <- mean(range(values))
    representative <- values[[which.min(abs(values - midpoint))]]
    data.frame(
      interval_id = intervals$interval_id[[i]],
      representative_resolution = representative,
      backbone_resolution = min(values),
      finest_supported_resolution = max(values),
      selection_rule = "nearest_evaluated_resolution_to_interval_midpoint",
      stringsAsFactors = FALSE
    )
  }))
}

edge_path <- file.path(stage_dir(config, "07_resolution_transitions", "tables"), "parent_daughter_edges.csv")
edges <- read_table_auto(edge_path)
event_rows <- lapply(split(edges, interaction(edges$lower_resolution, edges$upper_resolution, drop = TRUE)), function(x) {
  parent_daughters <- table(x$parent_cluster)
  daughter_parents <- table(x$daughter_cluster)
  data.frame(
    lower_resolution = x$lower_resolution[[1L]],
    upper_resolution = x$upper_resolution[[1L]],
    n_split_parents = sum(parent_daughters > 1L),
    n_merge_daughters = sum(daughter_parents > 1L),
    n_edges = nrow(x),
    stringsAsFactors = FALSE
  )
})
events <- bind_rows_base(event_rows)

split_rows <- list()
split_index <- 1L
trace_daughter <- function(start_resolution, start_cluster) {
  current_resolution <- start_resolution
  current_cluster <- as.character(start_cluster)
  steps <- 0L
  endpoint <- current_cluster
  repeat {
    next_edges <- edges[
      abs(edges$lower_resolution - current_resolution) < 1e-10 &
        edges$parent_cluster == current_cluster, , drop = FALSE
    ]
    if (!nrow(next_edges)) break
    best <- next_edges[which.max(next_edges$parent_fraction), , drop = FALSE]
    if (best$parent_fraction[[1L]] < 0.5) break
    steps <- steps + 1L
    current_resolution <- best$upper_resolution[[1L]]
    current_cluster <- as.character(best$daughter_cluster[[1L]])
    endpoint <- current_cluster
  }
  c(steps = steps, endpoint_resolution = current_resolution, endpoint_cluster = endpoint)
}
for (key in unique(paste(edges$lower_resolution, edges$upper_resolution, edges$parent_cluster, sep = "\r"))) {
  fields <- strsplit(key, "\r", fixed = TRUE)[[1L]]
  x <- edges[
    abs(edges$lower_resolution - as.numeric(fields[[1L]])) < 1e-10 &
      abs(edges$upper_resolution - as.numeric(fields[[2L]])) < 1e-10 &
      edges$parent_cluster == fields[[3L]] &
      edges$parent_fraction >= 0.05, , drop = FALSE
  ]
  if (nrow(x) < 2L) next
  traces <- lapply(seq_len(nrow(x)), function(i) {
    trace_daughter(x$upper_resolution[[i]], x$daughter_cluster[[i]])
  })
  persistence_steps <- as.numeric(vapply(traces, `[[`, character(1), "steps"))
  endpoint_resolutions <- as.numeric(vapply(traces, `[[`, character(1), "endpoint_resolution"))
  split_rows[[split_index]] <- data.frame(
    split_id = paste0(
      "r", fmt_resolution(x$lower_resolution[[1L]]), "_p", x$parent_cluster[[1L]],
      "_to_r", fmt_resolution(x$upper_resolution[[1L]])
    ),
    lower_resolution = x$lower_resolution[[1L]],
    upper_resolution = x$upper_resolution[[1L]],
    parent_cluster = x$parent_cluster[[1L]],
    n_daughters = nrow(x),
    largest_daughter_fraction = max(x$parent_fraction),
    second_daughter_fraction = sort(x$parent_fraction, decreasing = TRUE)[[2L]],
    balanced_split = sort(x$parent_fraction, decreasing = TRUE)[[2L]] >= 0.1,
    minimum_daughter_persistence_steps = min(persistence_steps),
    median_daughter_persistence_steps = stats::median(persistence_steps),
    minimum_daughter_endpoint_resolution = min(endpoint_resolutions),
    split_status = if (min(persistence_steps) >= 2L) {
      "persistent"
    } else if (max(persistence_steps) == 0L) {
      "transient"
    } else {
      "partially_persistent"
    },
    stringsAsFactors = FALSE
  )
  split_index <- split_index + 1L
}
split_persistence <- bind_rows_base(split_rows)

evidence_path <- file.path(dirs$tables, paste0("resolution_evidence_matrix__", tag, ".csv"))
cluster_evidence_path <- file.path(dirs$tables, paste0("cluster_evidence_matrix__", tag, ".csv"))
fragile_driver_path <- file.path(
  dirs$tables, paste0("fragile_cluster_drivers__", tag, ".csv")
)
driver_summary_path <- file.path(
  dirs$tables, paste0("fragile_cluster_driver_summary__", tag, ".csv")
)
interval_path <- file.path(dirs$tables, paste0("stable_resolution_intervals__", tag, ".csv"))
representative_path <- file.path(dirs$tables, paste0("interval_representatives__", tag, ".csv"))
candidate_transition_path <- file.path(dirs$tables, paste0("candidate_transition_similarity__", tag, ".csv"))
event_path <- file.path(stage_dir(config, "07_resolution_transitions", "tables"), "transition_event_summary.csv")
split_path <- file.path(stage_dir(config, "07_resolution_transitions", "tables"), "split_persistence.csv")
write_csv_atomic(evidence, evidence_path)
write_csv_atomic(cluster_evidence, cluster_evidence_path)
write_csv_atomic(fragile_drivers, fragile_driver_path)
write_csv_atomic(driver_summary, driver_summary_path)
write_csv_atomic(intervals, interval_path)
write_csv_atomic(representatives, representative_path)
write_csv_atomic(candidate_adjacent, candidate_transition_path)
write_csv_atomic(events, event_path)
write_csv_atomic(split_persistence, split_path)

plot_data <- rbind(
  data.frame(resolution = evidence$resolution, metric = "resample_jaccard", value = evidence$median_cluster_jaccard),
  data.frame(resolution = evidence$resolution, metric = "stable_cell_fraction", value = evidence$median_stable_cell_fraction),
  data.frame(resolution = evidence$resolution, metric = "sensitivity_ARI", value = evidence$median_reference_ARI),
  data.frame(resolution = evidence$resolution, metric = "pair_margin", value = evidence$median_consensus_margin)
)
save_plot(
  ggplot(plot_data, aes(resolution, value, color = metric)) +
    geom_line(linewidth = 0.45) + geom_point(size = 1.6) + facet_wrap(~metric, scales = "free_y") +
    theme_bw() + theme(legend.position = "none") +
    labs(x = "Resolution", y = "Metric value"),
  file.path(dirs$figures, paste0("stability_vs_resolution__", tag, ".pdf")),
  width = 9, height = 6
)

cluster_heatmap <- rbind(
  data.frame(
    resolution = cluster_evidence$resolution,
    reference_cluster = cluster_evidence$reference_cluster,
    metric = "resample median Jaccard",
    value = cluster_evidence$median_jaccard
  ),
  data.frame(
    resolution = cluster_evidence$resolution,
    reference_cluster = cluster_evidence$reference_cluster,
    metric = "pair consensus margin",
    value = cluster_evidence$consensus_margin
  ),
  data.frame(
    resolution = cluster_evidence$resolution,
    reference_cluster = cluster_evidence$reference_cluster,
    metric = "sensitivity median Jaccard",
    value = cluster_evidence$sensitivity_median_jaccard
  )
)
cluster_heatmap$reference_cluster <- as.numeric(as.character(cluster_heatmap$reference_cluster))
cluster_heatmap$metric <- factor(
  cluster_heatmap$metric,
  levels = c(
    "resample median Jaccard",
    "pair consensus margin",
    "sensitivity median Jaccard"
  )
)
cluster_heatmap_path <- file.path(
  dirs$figures, paste0("cluster_joint_evidence_heatmap__", tag, ".pdf")
)
save_plot(
  ggplot(cluster_heatmap, aes(metric, reference_cluster, fill = value)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_wrap(~resolution, scales = "free_y") +
    scale_fill_viridis_c(limits = c(0, 1), na.value = "grey85") +
    scale_y_continuous(breaks = function(x) pretty(x, n = 10)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(x = NULL, y = "Reference cluster", fill = "Value"),
  cluster_heatmap_path,
  width = 12, height = 8
)

driver_plot_data <- rbind(
  data.frame(
    resolution = driver_summary$resolution,
    evidence_layer = "resampling",
    failed_clusters = driver_summary$resample_fail_count
  ),
  data.frame(
    resolution = driver_summary$resolution,
    evidence_layer = "pair consensus",
    failed_clusters = driver_summary$pair_fail_count
  ),
  data.frame(
    resolution = driver_summary$resolution,
    evidence_layer = "parameter sensitivity",
    failed_clusters = driver_summary$sensitivity_fail_count
  )
)
driver_count_plot_path <- file.path(
  dirs$figures, paste0("fragile_cluster_counts__", tag, ".pdf")
)
save_plot(
  ggplot(driver_plot_data, aes(factor(resolution), failed_clusters, fill = evidence_layer)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72) +
    theme_bw() +
    labs(x = "Resolution", y = "Clusters below layer threshold", fill = "Evidence layer"),
  driver_count_plot_path,
  width = 8, height = 4.8
)

report_path <- file.path(dirs$report, paste0("stability_range_report__", tag, ".md"))
report <- c(
  paste0("# Stability range report: ", config$run_id),
  "",
  paste0("- Mode: `", mode, "`"),
  paste0("- Count mode: `", count_mode, "`"),
  paste0("- Evaluated resolutions: ", paste(resolutions, collapse = ", ")),
  paste0("- Stable intervals found: ", nrow(intervals)),
  paste0("- Cluster Jaccard threshold: ", cluster_jaccard_threshold),
  paste0("- Pair consensus margin threshold: ", pair_margin_threshold),
  "- Annotation used in structural evidence: no",
  "- Unique best resolution written: no",
  "",
  "Stable intervals are contiguous evaluated resolutions that pass all required structural",
  "criteria and whose successive partitions meet the configured adjacent-ARI threshold.",
  "Pilot intervals are provisional until formal resampling is completed.",
  "",
  "The cluster evidence matrix joins resampling, sampled-pair consensus, and PC/k/seed",
  "sensitivity by resolution and reference cluster. The fragile-cluster table contains",
  "all threshold failures plus the three worst clusters for each principal metric.",
  "Failure counts are descriptive flags; they are not combined into a weighted score."
)
writeLines(report, report_path)

freeze_path <- NULL
if (!is.null(cli$freeze_resolution)) {
  freeze <- as.numeric(cli$freeze_resolution)
  if (is.na(freeze) || !any(abs(resolutions - freeze) < 1e-10)) {
    abort("--freeze-resolution must be one of the summarized resolutions.")
  }
  if (!isTRUE(evidence$core_supported[which.min(abs(evidence$resolution - freeze))])) {
    abort("Resolution %s is not supported by the configured structural criteria.", fmt_resolution(freeze))
  }
  assert_packages("digest")
  selected <- data.frame(
    cell_id = assignments$cell_id,
    resolution = freeze,
    cluster = assignments[[resolution_column(freeze)]],
    stringsAsFactors = FALSE
  )
  selected_path <- file.path(dirs$tables, paste0("working_assignment_r", fmt_resolution(freeze), ".tsv.gz"))
  write_tsv_gz_atomic(selected, selected_path)
  loaded <- load_and_validate_object(config)
  prepared <- prepare_counts(loaded$counts, config, count_mode)
  session_path <- file.path(dirs$logs, paste0("session_info__", tag, ".txt"))
  writeLines(capture.output(sessionInfo()), session_path)
  freeze_record <- data.frame(
    run_id = config$run_id,
    frozen_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    resolution = freeze,
    mode = mode,
    count_mode = count_mode,
    input_rds_md5 = file_md5(config$input$seurat_rds),
    count_matrix_digest = digest::digest(prepared$counts, algo = "xxhash64"),
    config_md5 = file_md5(cli$config),
    resolution_set_id = resolution_set_id(resolutions),
    full_assignment_md5 = file_md5(assignments_path),
    selected_assignment_md5 = file_md5(selected_path),
    evidence_matrix_md5 = file_md5(evidence_path),
    stable_intervals_md5 = file_md5(interval_path),
    marker_mode = if (is.null(marker)) "none" else "fixed_expression",
    seed = config$compute$seed,
    session_info_md5 = file_md5(session_path),
    stringsAsFactors = FALSE
  )
  freeze_path <- file.path(dirs$tables, "working_resolution_freeze.csv")
  write_csv_atomic(freeze_record, freeze_path)
}

outputs <- c(
  evidence_path, interval_path, representative_path, candidate_transition_path,
  cluster_evidence_path, fragile_driver_path, driver_summary_path,
  cluster_heatmap_path, driver_count_plot_path,
  report_path, event_path, split_path
)
if (!is.null(freeze_path)) outputs <- c(outputs, freeze_path)
write_stage_manifest(
  config, "08_stability_range", resolutions,
  inputs = c(assignments_path),
  outputs = outputs,
  details = list(mode = mode, count_mode = count_mode, frozen = !is.null(freeze_path)),
  manifest_tag = tag
)
writeLines(capture.output(sessionInfo()), file.path(dirs$logs, paste0("session_info__", tag, ".txt")))
cat("Stable intervals:", nrow(intervals), "\n")
cat("Finished:", format(Sys.time()), "\n")
