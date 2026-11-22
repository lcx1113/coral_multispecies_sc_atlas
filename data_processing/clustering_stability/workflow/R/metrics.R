partition_similarity <- function(reference, test) {
  keep <- !is.na(reference) & !is.na(test)
  reference <- as.character(reference[keep])
  test <- as.character(test[keep])
  if (length(reference) < 2L) return(data.frame(ARI = NA_real_, NMI = NA_real_))
  data.frame(
    ARI = mclust::adjustedRandIndex(reference, test),
    NMI = aricode::NMI(reference, test)
  )
}

cluster_size_summary <- function(labels, small_threshold = 50L, min_threshold = 20L) {
  sizes <- as.integer(table(labels))
  data.frame(
    n_clusters = length(sizes),
    min_cluster_size = min(sizes),
    median_cluster_size = stats::median(sizes),
    max_cluster_size = max(sizes),
    n_below_min = sum(sizes < min_threshold),
    n_small_clusters = sum(sizes < small_threshold),
    small_cluster_fraction = mean(sizes < small_threshold)
  )
}

jaccard_reference_clusters <- function(reference, test, sampled_cells = NULL) {
  if (is.null(names(reference)) || is.null(names(test))) abort("Cluster vectors must be named by cell.")
  common <- intersect(names(reference), names(test))
  if (!is.null(sampled_cells)) common <- intersect(common, sampled_cells)
  reference <- as.character(reference[common])
  test <- as.character(test[common])
  ref_levels <- sort(unique(reference))
  test_levels <- sort(unique(test))
  rows <- lapply(ref_levels, function(ref) {
    ref_cells <- common[reference == ref]
    scores <- vapply(test_levels, function(candidate) {
      test_cells <- common[test == candidate]
      union_n <- length(union(ref_cells, test_cells))
      if (!union_n) NA_real_ else length(intersect(ref_cells, test_cells)) / union_n
    }, numeric(1))
    order_idx <- order(scores, decreasing = TRUE, na.last = TRUE)
    best <- if (length(order_idx)) order_idx[[1L]] else NA_integer_
    second <- if (length(order_idx) > 1L) order_idx[[2L]] else NA_integer_
    matched <- if (is.na(best)) NA_character_ else test_levels[[best]]
    matched_cells <- if (is.na(best)) character() else common[test == matched]
    data.frame(
      reference_cluster = ref,
      matched_cluster = matched,
      reference_cluster_sampled_n = length(ref_cells),
      matched_cluster_n = length(matched_cells),
      intersection_n = length(intersect(ref_cells, matched_cells)),
      max_jaccard = if (is.na(best)) NA_real_ else scores[[best]],
      second_jaccard = if (is.na(second)) NA_real_ else scores[[second]],
      best_second_jaccard_gap = if (is.na(best)) NA_real_ else scores[[best]] - ifelse(is.na(second), 0, scores[[second]]),
      reference_recall = if (!length(ref_cells)) NA_real_ else length(intersect(ref_cells, matched_cells)) / length(ref_cells),
      test_precision = if (!length(matched_cells)) NA_real_ else length(intersect(ref_cells, matched_cells)) / length(matched_cells),
      stringsAsFactors = FALSE
    )
  })
  bind_rows_base(rows)
}

partitions_equivalent <- function(a, b) {
  keep <- !is.na(a) & !is.na(b)
  tab <- table(as.character(a[keep]), as.character(b[keep]))
  if (!nrow(tab) || !ncol(tab) || nrow(tab) != ncol(tab)) return(FALSE)
  all(rowSums(tab > 0L) == 1L) && all(colSums(tab > 0L) == 1L)
}

adjacent_partition_metrics <- function(assignments, resolutions) {
  if (length(resolutions) < 2L) return(data.frame())
  rows <- lapply(seq_len(length(resolutions) - 1L), function(i) {
    left <- resolutions[[i]]
    right <- resolutions[[i + 1L]]
    a <- assignments[[resolution_column(left)]]
    b <- assignments[[resolution_column(right)]]
    similarity <- partition_similarity(a, b)
    tab <- table(as.character(a), as.character(b))
    parent_sizes <- rowSums(tab)
    daughters <- rowSums(tab > 0L)
    data.frame(
      lower_resolution = left,
      upper_resolution = right,
      lower_n_clusters = length(unique(a)),
      upper_n_clusters = length(unique(b)),
      delta_clusters = length(unique(b)) - length(unique(a)),
      adjacent_ARI = similarity$ARI,
      adjacent_NMI = similarity$NMI,
      equivalent_partition = partitions_equivalent(a, b),
      n_split_parents = sum(daughters > 1L),
      weighted_split_fraction = sum(parent_sizes[daughters > 1L]) / sum(parent_sizes),
      stringsAsFactors = FALSE
    )
  })
  bind_rows_base(rows)
}

equivalent_resolution_intervals <- function(adjacent, resolutions) {
  if (!length(resolutions)) return(data.frame())
  group <- integer(length(resolutions))
  group[[1L]] <- 1L
  if (length(resolutions) > 1L) {
    for (i in 2:length(resolutions)) {
      group[[i]] <- group[[i - 1L]] + as.integer(!isTRUE(adjacent$equivalent_partition[[i - 1L]]))
    }
  }
  split_res <- split(resolutions, group)
  bind_rows_base(lapply(names(split_res), function(id) {
    values <- split_res[[id]]
    data.frame(
      equivalent_interval_id = paste0("eq_", id),
      resolution_min = min(values),
      resolution_max = max(values),
      n_grid_points = length(values),
      representative_resolution = values[[ceiling(length(values) / 2)]],
      resolutions = paste(vapply(values, fmt_resolution, character(1)), collapse = ","),
      stringsAsFactors = FALSE
    )
  }))
}

parent_daughter_edges <- function(assignments, resolutions) {
  if (length(resolutions) < 2L) return(data.frame())
  edges <- lapply(seq_len(length(resolutions) - 1L), function(i) {
    left <- resolutions[[i]]
    right <- resolutions[[i + 1L]]
    parent <- as.character(assignments[[resolution_column(left)]])
    daughter <- as.character(assignments[[resolution_column(right)]])
    tab <- as.data.frame(table(parent_cluster = parent, daughter_cluster = daughter), stringsAsFactors = FALSE)
    tab <- tab[tab$Freq > 0L, , drop = FALSE]
    parent_n <- setNames(as.numeric(table(parent)), names(table(parent)))
    daughter_n <- setNames(as.numeric(table(daughter)), names(table(daughter)))
    tab$lower_resolution <- left
    tab$upper_resolution <- right
    tab$overlap_n <- tab$Freq
    tab$parent_fraction <- tab$overlap_n / parent_n[tab$parent_cluster]
    tab$daughter_fraction <- tab$overlap_n / daughter_n[tab$daughter_cluster]
    tab$Freq <- NULL
    tab[c(
      "lower_resolution", "upper_resolution", "parent_cluster", "daughter_cluster",
      "overlap_n", "parent_fraction", "daughter_fraction"
    )]
  })
  bind_rows_base(edges)
}

summarize_resample_jaccard <- function(jaccard, stable_threshold, min_sampled) {
  if (!nrow(jaccard)) return(data.frame())
  keys <- interaction(jaccard$resolution, jaccard$iteration, drop = TRUE)
  per_iteration <- bind_rows_base(lapply(split(jaccard, keys), function(x) {
    sufficient <- x$reference_cluster_sampled_n >= min_sampled
    values <- x$max_jaccard[sufficient]
    data.frame(
      mode = x$mode[[1L]],
      resolution = x$resolution[[1L]],
      iteration = x$iteration[[1L]],
      seed = x$seed[[1L]],
      n_sufficient_clusters = sum(sufficient),
      median_cluster_jaccard = safe_median(values),
      q10_cluster_jaccard = safe_quantile(values, 0.1),
      stable_cluster_fraction = if (length(values)) mean(values >= stable_threshold) else NA_real_,
      stable_cell_fraction = if (any(sufficient)) {
        stats::weighted.mean(
          x$max_jaccard[sufficient] >= stable_threshold,
          x$reference_cluster_sampled_n[sufficient]
        )
      } else NA_real_,
      partition_ARI = x$partition_ARI[[1L]],
      partition_NMI = x$partition_NMI[[1L]],
      stringsAsFactors = FALSE
    )
  }))
  per_resolution <- bind_rows_base(lapply(split(per_iteration, per_iteration$resolution), function(x) {
    data.frame(
      mode = x$mode[[1L]],
      resolution = x$resolution[[1L]],
      n_iterations = nrow(x),
      median_cluster_jaccard = safe_median(x$median_cluster_jaccard),
      q10_cluster_jaccard = safe_quantile(x$median_cluster_jaccard, 0.1),
      median_stable_cluster_fraction = safe_median(x$stable_cluster_fraction),
      median_stable_cell_fraction = safe_median(x$stable_cell_fraction),
      median_partition_ARI = safe_median(x$partition_ARI),
      q10_partition_ARI = safe_quantile(x$partition_ARI, 0.1),
      median_partition_NMI = safe_median(x$partition_NMI),
      stringsAsFactors = FALSE
    )
  }))
  list(per_iteration = per_iteration, per_resolution = per_resolution)
}

summarize_cluster_jaccard <- function(jaccard, stable_threshold, min_sampled) {
  if (!nrow(jaccard)) return(data.frame())
  bind_rows_base(lapply(
    split(jaccard, interaction(jaccard$resolution, jaccard$reference_cluster, drop = TRUE)),
    function(x) {
      sufficient <- x$reference_cluster_sampled_n >= min_sampled
      values <- x$max_jaccard[sufficient]
      data.frame(
        mode = x$mode[[1L]],
        resolution = x$resolution[[1L]],
        reference_cluster = x$reference_cluster[[1L]],
        n_iterations = nrow(x),
        n_sufficient_iterations = sum(sufficient),
        median_sampled_cells = safe_median(x$reference_cluster_sampled_n[sufficient]),
        median_jaccard = safe_median(values),
        q10_jaccard = safe_quantile(values, 0.1),
        stable_iteration_fraction = if (length(values)) mean(values >= stable_threshold) else NA_real_,
        median_reference_recall = safe_median(x$reference_recall[sufficient]),
        median_test_precision = safe_median(x$test_precision[sufficient]),
        median_match_margin = safe_median(x$best_second_jaccard_gap[sufficient]),
        stringsAsFactors = FALSE
      )
    }
  ))
}

select_scan_pilot_candidates <- function(size_summary, adjacent, settings) {
  resolutions <- sort(unique(size_summary$resolution))
  stable_edge <- adjacent$adjacent_ARI >= settings$stable_edge_min_ari &
    adjacent$adjacent_NMI >= settings$stable_edge_min_nmi &
    abs(adjacent$delta_clusters) <= settings$stable_edge_max_abs_cluster_delta

  candidate_rows <- list()
  index <- 1L
  add_candidate <- function(resolution, reason, detail, priority) {
    candidate_rows[[index]] <<- data.frame(
      resolution = resolution,
      selection_reason = reason,
      selection_detail = detail,
      priority = priority,
      stringsAsFactors = FALSE
    )
    index <<- index + 1L
  }

  runs <- rle(stable_edge)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  for (run_index in which(runs$values & runs$lengths >= settings$minimum_stable_run_edges)) {
    edge_indices <- starts[[run_index]]:ends[[run_index]]
    node_values <- c(
      adjacent$lower_resolution[edge_indices[[1L]]],
      adjacent$upper_resolution[edge_indices]
    )
    midpoint <- mean(range(node_values))
    representative <- node_values[[which.min(abs(node_values - midpoint))]]
    add_candidate(
      representative,
      "stable_run_representative",
      sprintf(
        "stable_edges=%d;range=%s-%s",
        length(edge_indices), fmt_resolution(min(node_values)), fmt_resolution(max(node_values))
      ),
      2L
    )
  }

  edge_midpoint <- (adjacent$lower_resolution + adjacent$upper_resolution) / 2
  band <- cut(
    edge_midpoint,
    breaks = seq(min(resolutions), max(resolutions), length.out = settings$transition_bands + 1L),
    include.lowest = TRUE,
    labels = FALSE
  )
  for (band_id in sort(unique(band))) {
    edge_indices <- which(band == band_id)
    order_indices <- edge_indices[order(
      adjacent$adjacent_ARI[edge_indices],
      adjacent$adjacent_NMI[edge_indices],
      -abs(adjacent$delta_clusters[edge_indices])
    )]
    selected_edges <- head(order_indices, settings$transition_edges_per_band)
    for (edge_index in selected_edges) {
      detail <- sprintf(
        "transition=%s-%s;ARI=%.3f;NMI=%.3f;delta_clusters=%+d",
        fmt_resolution(adjacent$lower_resolution[[edge_index]]),
        fmt_resolution(adjacent$upper_resolution[[edge_index]]),
        adjacent$adjacent_ARI[[edge_index]],
        adjacent$adjacent_NMI[[edge_index]],
        adjacent$delta_clusters[[edge_index]]
      )
      add_candidate(adjacent$lower_resolution[[edge_index]], "transition_lower_boundary", detail, 1L)
      add_candidate(adjacent$upper_resolution[[edge_index]], "transition_upper_boundary", detail, 1L)
    }
  }

  if (isTRUE(settings$include_range_endpoints)) {
    add_candidate(min(resolutions), "scan_range_endpoint", "lower_endpoint", 3L)
    add_candidate(max(resolutions), "scan_range_endpoint", "upper_endpoint", 3L)
  }

  candidates <- bind_rows_base(candidate_rows)
  candidates <- bind_rows_base(lapply(split(candidates, candidates$resolution), function(x) {
    data.frame(
      resolution = x$resolution[[1L]],
      priority = min(x$priority),
      selection_reason = paste(unique(x$selection_reason), collapse = ";"),
      selection_detail = paste(unique(x$selection_detail), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
  candidates <- merge(
    candidates,
    size_summary[c(
      "resolution", "n_clusters", "min_cluster_size", "n_below_min",
      "n_small_clusters", "small_cell_fraction", "largest_cluster_fraction"
    )],
    by = "resolution",
    all.x = TRUE,
    sort = FALSE
  )
  candidates[order(candidates$priority, candidates$resolution), , drop = FALSE]
}
