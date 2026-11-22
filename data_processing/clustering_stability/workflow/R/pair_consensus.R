sample_within_pairs <- function(cells, maximum) {
  cells <- as.character(cells)
  n <- length(cells)
  possible <- n * (n - 1) / 2
  target <- min(as.numeric(maximum), possible)
  if (target <= 0) return(data.frame(cell_i = character(), cell_j = character()))
  if (possible <= maximum) {
    indices <- utils::combn(seq_len(n), 2)
    return(data.frame(cell_i = cells[indices[1L, ]], cell_j = cells[indices[2L, ]]))
  }
  keys <- character()
  while (length(keys) < target) {
    batch <- max(1000L, 2L * (target - length(keys)))
    first <- sample.int(n, batch, replace = TRUE)
    second <- sample.int(n, batch, replace = TRUE)
    keep <- first != second
    lower <- pmin(first[keep], second[keep])
    upper <- pmax(first[keep], second[keep])
    keys <- unique(c(keys, paste(lower, upper, sep = ":")))
  }
  keys <- keys[seq_len(target)]
  parts <- strsplit(keys, ":", fixed = TRUE)
  first <- as.integer(vapply(parts, `[[`, character(1), 1L))
  second <- as.integer(vapply(parts, `[[`, character(1), 2L))
  data.frame(cell_i = cells[first], cell_j = cells[second])
}

sample_between_pairs <- function(cells_a, cells_b, maximum) {
  cells_a <- as.character(cells_a)
  cells_b <- as.character(cells_b)
  possible <- length(cells_a) * length(cells_b)
  target <- min(as.numeric(maximum), possible)
  if (target <= 0) return(data.frame(cell_i = character(), cell_j = character()))
  linear <- if (target == possible) seq_len(possible) else sample.int(possible, target)
  index_a <- ((linear - 1L) %/% length(cells_b)) + 1L
  index_b <- ((linear - 1L) %% length(cells_b)) + 1L
  data.frame(cell_i = cells_a[index_a], cell_j = cells_b[index_b])
}

neighbor_cluster_pairs <- function(labels, embedding, top_k = 3L) {
  labels <- as.character(labels)
  clusters <- sort(unique(labels))
  if (length(clusters) < 2L) return(data.frame(cluster_a = character(), cluster_b = character()))
  centroids <- do.call(rbind, lapply(clusters, function(cluster) {
    colMeans(embedding[labels == cluster, , drop = FALSE])
  }))
  rownames(centroids) <- clusters
  distance <- as.matrix(stats::dist(centroids))
  diag(distance) <- Inf
  pairs <- lapply(seq_along(clusters), function(i) {
    neighbors <- order(distance[i, ], decreasing = FALSE)[seq_len(min(top_k, length(clusters) - 1L))]
    data.frame(cluster_a = clusters[[i]], cluster_b = clusters[neighbors], stringsAsFactors = FALSE)
  })
  pairs <- bind_rows_base(pairs)
  lower <- pmin(pairs$cluster_a, pairs$cluster_b)
  upper <- pmax(pairs$cluster_a, pairs$cluster_b)
  pairs <- unique(data.frame(cluster_a = lower, cluster_b = upper, stringsAsFactors = FALSE))
  pairs
}

run_pair_consensus <- function(config, resolutions, reference, resampled, embedding, tag) {
  assert_packages(c("data.table", "ggplot2"))
  dirs <- ensure_stage_dirs(config, "04_pair_consensus")
  set.seed(config$pair_consensus$seed)
  pair_sets <- list()
  pair_index <- 1L
  for (resolution in resolutions) {
    labels <- as.character(reference[[resolution_column(resolution)]])
    names(labels) <- reference$cell_id
    for (cluster in sort(unique(labels))) {
      pairs <- sample_within_pairs(
        names(labels)[labels == cluster],
        config$pair_consensus$max_within_pairs_per_cluster
      )
      if (!nrow(pairs)) next
      pairs$resolution <- resolution
      pairs$pair_type <- "within"
      pairs$reference_cluster_a <- cluster
      pairs$reference_cluster_b <- cluster
      pair_sets[[pair_index]] <- pairs
      pair_index <- pair_index + 1L
    }
    neighbors <- neighbor_cluster_pairs(
      labels,
      embedding[
        reference$cell_id,
        seq_len(primary_pcs_from_config(config)),
        drop = FALSE
      ],
      config$pair_consensus$neighboring_clusters_per_cluster
    )
    for (row in seq_len(nrow(neighbors))) {
      cluster_a <- neighbors$cluster_a[[row]]
      cluster_b <- neighbors$cluster_b[[row]]
      pairs <- sample_between_pairs(
        names(labels)[labels == cluster_a],
        names(labels)[labels == cluster_b],
        config$pair_consensus$max_between_pairs_per_cluster_pair
      )
      if (!nrow(pairs)) next
      pairs$resolution <- resolution
      pairs$pair_type <- "between"
      pairs$reference_cluster_a <- cluster_a
      pairs$reference_cluster_b <- cluster_b
      pair_sets[[pair_index]] <- pairs
      pair_index <- pair_index + 1L
    }
  }
  pairs <- data.table::as.data.table(bind_rows_base(pair_sets))
  pairs[, pair_id := .I]
  pairs[, `:=`(co_sampled_count = 0L, co_clustered_count = 0L)]
  assignments <- data.table::as.data.table(resampled)
  for (resolution_value in resolutions) {
    pair_indices <- which(abs(pairs$resolution - resolution_value) < 1e-10)
    current_resolution <- assignments[abs(resolution - resolution_value) < 1e-10]
    for (iteration_value in sort(unique(current_resolution$iteration))) {
      current <- current_resolution[iteration == iteration_value]
      map <- setNames(as.character(current$resample_cluster), current$cell_id)
      first <- unname(map[pairs$cell_i[pair_indices]])
      second <- unname(map[pairs$cell_j[pair_indices]])
      observed <- !is.na(first) & !is.na(second)
      pairs$co_sampled_count[pair_indices] <- pairs$co_sampled_count[pair_indices] + as.integer(observed)
      pairs$co_clustered_count[pair_indices] <- pairs$co_clustered_count[pair_indices] +
        as.integer(observed & first == second)
    }
  }
  pairs[, co_clustering_probability := data.table::fifelse(
    co_sampled_count >= config$pair_consensus$minimum_joint_observations,
    co_clustered_count / co_sampled_count,
    NA_real_
  )]
  median_or_na <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  within <- pairs[pair_type == "within", .(
    sampled_pair_count = .N,
    valid_pair_count = sum(!is.na(co_clustering_probability)),
    valid_pair_fraction = mean(!is.na(co_clustering_probability)),
    within_consensus = median_or_na(co_clustering_probability)
  ), by = .(resolution, reference_cluster = reference_cluster_a)]
  between <- pairs[pair_type == "between", .(
    sampled_pair_count = .N,
    valid_pair_count = sum(!is.na(co_clustering_probability)),
    valid_pair_fraction = mean(!is.na(co_clustering_probability)),
    between_leakage = median_or_na(co_clustering_probability)
  ), by = .(resolution, cluster_a = reference_cluster_a, cluster_b = reference_cluster_b)]
  sizes <- data.table::rbindlist(lapply(resolutions, function(resolution) {
    tab <- table(reference[[resolution_column(resolution)]])
    data.table::data.table(
      resolution = resolution,
      reference_cluster = names(tab),
      cluster_size = as.integer(tab)
    )
  }))
  cluster_summary <- merge(within, sizes, by = c("resolution", "reference_cluster"), all.x = TRUE)
  cluster_summary[, `:=`(max_between_leakage = NA_real_, nearest_leaking_cluster = NA_character_)]
  for (row in seq_len(nrow(cluster_summary))) {
    resolution_value <- cluster_summary$resolution[[row]]
    cluster_value <- cluster_summary$reference_cluster[[row]]
    candidates <- between[
      abs(resolution - resolution_value) < 1e-10 &
        (cluster_a == cluster_value | cluster_b == cluster_value)
    ]
    if (nrow(candidates) && !all(is.na(candidates$between_leakage))) {
      best <- which.max(candidates$between_leakage)
      cluster_summary$max_between_leakage[[row]] <- candidates$between_leakage[[best]]
      cluster_summary$nearest_leaking_cluster[[row]] <- if (candidates$cluster_a[[best]] == cluster_value) {
        candidates$cluster_b[[best]]
      } else {
        candidates$cluster_a[[best]]
      }
    }
  }
  cluster_summary[, consensus_margin := within_consensus - max_between_leakage]
  resolution_summary <- cluster_summary[, .(
    median_within_consensus = median_or_na(within_consensus),
    median_max_between_leakage = median_or_na(max_between_leakage),
    median_consensus_margin = median_or_na(consensus_margin),
    median_valid_pair_fraction = median_or_na(valid_pair_fraction)
  ), by = resolution]

  pair_path <- file.path(dirs$tables, paste0("sampled_pair_probabilities__", tag, ".tsv.gz"))
  cluster_path <- file.path(dirs$tables, paste0("pair_consensus_cluster_summary__", tag, ".csv"))
  between_path <- file.path(dirs$tables, paste0("between_neighbor_leakage__", tag, ".csv"))
  resolution_path <- file.path(dirs$tables, paste0("pair_consensus_resolution_summary__", tag, ".csv"))
  write_tsv_gz_atomic(as.data.frame(pairs), pair_path)
  write_csv_atomic(as.data.frame(cluster_summary), cluster_path)
  write_csv_atomic(as.data.frame(between), between_path)
  write_csv_atomic(as.data.frame(resolution_summary), resolution_path)
  save_plot(
    ggplot2::ggplot(cluster_summary, ggplot2::aes(factor(resolution), consensus_margin)) +
      ggplot2::geom_boxplot(outlier.size = 0.5) + ggplot2::theme_bw() +
      ggplot2::labs(x = "Resolution", y = "Within consensus - maximum neighbor leakage"),
    file.path(dirs$figures, paste0("co_clustering_margin__", tag, ".pdf"))
  )
  write_stage_manifest(
    config, "04_pair_consensus", resolutions,
    inputs = character(), outputs = c(pair_path, cluster_path, between_path, resolution_path),
    details = list(pair_sampling = "within_and_neighbor_between", pair_count = nrow(pairs)),
    manifest_tag = tag
  )
  list(cluster = as.data.frame(cluster_summary), resolution = as.data.frame(resolution_summary))
}
