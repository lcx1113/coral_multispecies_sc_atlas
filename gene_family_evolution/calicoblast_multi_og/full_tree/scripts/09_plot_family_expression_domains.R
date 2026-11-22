#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(digest)
  library(ggnewscale)
  library(ggplot2)
  library(ggtree)
  library(scales)
  library(yaml)
})

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else file.path(project_root, "config", "analysis_config.yaml")
scope_label <- if (length(args) >= 2) args[[2]] else "combined_families"
track_config_path <- if (length(args) >= 3) args[[3]] else file.path(project_root, "config", "tip_tracks.yaml")
config_path <- normalizePath(config_path, mustWork = TRUE)
track_config_path <- normalizePath(track_config_path, mustWork = TRUE)
config <- yaml::read_yaml(config_path)
track_config <- yaml::read_yaml(track_config_path)
tracks <- track_config$multiple_og_tracks

stopifnot(
  identical(tracks$active, TRUE),
  identical(tracks$celltype, "calicoblast"),
  identical(tracks$reuse_precomputed_expression, TRUE),
  identical(tracks$reuse_precomputed_cdd, TRUE),
  identical(tracks$recompute_expression, FALSE),
  identical(tracks$recompute_cdd, FALSE),
  identical(config$expression_link$missing_expression_value, "NA_not_zero"),
  identical(config$visualization$show_iqtree_ufboot, TRUE),
  identical(config$visualization$filter_iqtree_ufboot, FALSE)
)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

validate_pdf <- function(path) {
  output <- suppressWarnings(system2(
    config$visualization$ghostscript_executable,
    c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", shQuote(path)),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), message = paste(output, collapse = " | "))
}

definitions <- read_tsv(file.path(project_root, config$multiple_og$family_definitions))
definitions <- definitions[toupper(definitions$active) == "TRUE", , drop = FALSE]
track_root <- file.path(project_root, tracks$output_directory)
tip_expression <- read_tsv(file.path(track_root, "tip_relative_expression.tsv"))
domain_hits <- read_tsv(file.path(track_root, "cdd_domain_hits.tsv"))
domain_status <- read_tsv(file.path(track_root, "cdd_sequence_status.tsv"))
coverage <- read_tsv(file.path(track_root, "tip_track_coverage.tsv"))
reconciliation_dir <- file.path(project_root, "results", "03_reconciliation", scope_label)
nodes <- read_tsv(file.path(reconciliation_dir, "reconciliation_nodes.tsv"))
events <- read_tsv(file.path(reconciliation_dir, "reconciliation_events.tsv"))
summaries <- read_tsv(file.path(reconciliation_dir, "reconciliation_summary.tsv"))

stopifnot(
  !anyDuplicated(paste(tip_expression$family_id, tip_expression$leaf_id, sep = "\r")),
  !anyDuplicated(paste(domain_status$family_id, domain_status$leaf_id, sep = "\r")),
  all(tip_expression$celltype == "calicoblast"),
  setequal(unique(tip_expression$family_id), summaries$orthogroup)
)

output_dir <- file.path(project_root, tracks$plot_output_directory)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

species_color <- config$visualization$node_color_speciation
duplication_color <- config$visualization$node_color_duplication
loss_color <- config$visualization$loss_color
low_color <- config$visualization$expression_log2fc_low_color
mid_color <- config$visualization$expression_log2fc_mid_color
high_color <- config$visualization$expression_log2fc_high_color
no_dataset_color <- config$visualization$expression_no_dataset_color
gene_not_tested_color <- config$visualization$expression_gene_not_tested_color
celltype_not_present_color <- config$visualization$expression_celltype_not_present_color
display_limit <- as.numeric(config$visualization$expression_log2fc_display_limit)
domain_label_font_size <- as.numeric(config$visualization$domain_label_font_size)
domain_label_angle <- as.numeric(config$visualization$domain_label_angle_degrees)
domain_lane_spacing <- as.numeric(config$visualization$domain_lane_spacing)
domain_rectangle_half_height <- as.numeric(config$visualization$domain_rectangle_half_height)
domain_label_offset <- as.numeric(config$visualization$domain_label_offset)

descendant_tip_labels <- function(tree) {
  n_tip <- length(tree$tip.label)
  n_total <- n_tip + tree$Nnode
  children <- split(tree$edge[, 2], tree$edge[, 1])
  cache <- vector("list", n_total)
  descend <- function(node) {
    if (!is.null(cache[[node]])) return(cache[[node]])
    result <- if (node <= n_tip) {
      tree$tip.label[[node]]
    } else {
      sort(unique(unlist(lapply(children[[as.character(node)]], descend), use.names = FALSE)))
    }
    cache[[node]] <<- result
    result
  }
  result <- lapply(seq_len(n_total), descend)
  names(result) <- as.character(seq_len(n_total))
  result
}

prepare_family <- function(family_id) {
  family_nodes <- nodes[nodes$orthogroup == family_id, , drop = FALSE]
  family_events <- events[events$orthogroup == family_id, , drop = FALSE]
  family_expression <- tip_expression[tip_expression$family_id == family_id, , drop = FALSE]
  family_domain_status <- domain_status[domain_status$family_id == family_id, , drop = FALSE]
  family_domain_hits <- domain_hits[domain_hits$family_id == family_id, , drop = FALSE]
  nhx_paths <- unique(family_nodes$reconciled_nhx)
  stopifnot(length(nhx_paths) == 1L, file.exists(nhx_paths))
  tree <- ape::read.tree(nhx_paths)
  descendants <- descendant_tip_labels(tree)
  descendant_keys <- vapply(descendants, function(x) paste(sort(x), collapse = ";"), character(1))
  node_index <- match(descendant_keys, family_nodes$descendant_leaves)
  if (anyNA(node_index) || anyDuplicated(node_index)) stop("Stable-node map failed for ", family_id)
  annotation <- family_nodes[node_index, , drop = FALSE]
  annotation$node <- seq_len(nrow(annotation))

  base <- ggtree::ggtree(tree, layout = "rectangular", branch.length = "none",
                         linewidth = 0.30, colour = "#525252")
  coords <- base$data[, c("node", "parent", "x", "y")]
  annotation <- merge(annotation, coords, by = "node", sort = FALSE, all.x = TRUE)
  annotation <- annotation[match(seq_len(nrow(coords)), annotation$node), , drop = FALSE]
  if (anyNA(annotation$x) || anyNA(annotation$y)) stop("Missing coordinates for ", family_id)

  tips <- annotation[annotation$node_kind == "leaf", , drop = FALSE]
  tips <- tips[match(tree$tip.label, tips$leaf_id), , drop = FALSE]
  expression_index <- match(tips$leaf_id, family_expression$leaf_id)
  status_index <- match(tips$leaf_id, family_domain_status$leaf_id)
  if (anyNA(expression_index) || anyDuplicated(expression_index) ||
      nrow(tips) != nrow(family_expression)) {
    stop("Expression does not cover the active tip set for ", family_id)
  }
  if (anyNA(status_index) || anyDuplicated(status_index) ||
      nrow(tips) != nrow(family_domain_status)) {
    stop("CDD status does not cover the active tip set for ", family_id)
  }
  expression_columns <- setdiff(names(family_expression), names(tips))
  tips <- cbind(tips, family_expression[expression_index, expression_columns, drop = FALSE])
  tips$protein_length_domain <- as.numeric(family_domain_status$protein_length[status_index])
  tips$domain_status <- family_domain_status$domain_status[status_index]
  if (anyNA(tips$protein_length_domain) || any(tips$protein_length_domain <= 0)) {
    stop("Invalid CDD protein length for ", family_id)
  }
  tips$plot_label <- paste(tips$tip_display_label, tips$original_orthogroup, sep = " | ")
  external <- tips$tree_species %in% c("Hvul", "Aaur")
  tips$plot_label[external] <- paste(tips$tree_species[external], tips$protein_id[external], sep = "|")
  tips$plot_log2fc <- NA_real_
  tips$plot_pct <- NA_real_
  observed <- tips$expression_observation_status == "observed"
  tips$plot_log2fc[observed] <- as.numeric(tips$avg_log2FC[observed])
  tips$plot_pct[observed] <- as.numeric(tips$pct_target[observed])
  if (anyNA(tips$plot_log2fc[observed]) || anyNA(tips$plot_pct[observed])) {
    stop("Observed expression contains nonnumeric values for ", family_id)
  }
  tips$fc_label <- ifelse(
    observed, sprintf("%.2f", tips$plot_log2fc),
    ifelse(tips$expression_observation_status == "no_expression_dataset", "ND",
    ifelse(tips$expression_observation_status == "gene_not_tested", "NT",
    ifelse(tips$expression_observation_status == "celltype_not_present", "CT", "NA")))
  )
  tips$marker_state <- ifelse(tips$positive_marker == "TRUE", "marker-positive", "tested-not-marker")
  tips$marker_state[!observed] <- NA_character_
  tips$marker_state <- factor(tips$marker_state, levels = c("marker-positive", "tested-not-marker"))

  internal <- annotation[annotation$node_kind == "internal", , drop = FALSE]
  internal$event_type <- factor(internal$event_type, levels = c("S", "D"))
  node_number_by_id <- setNames(annotation$node, annotation$node_id)
  x_by_node <- setNames(annotation$x, annotation$node)
  support <- internal[nzchar(internal$split_id), , drop = FALSE]
  support$parent_number <- node_number_by_id[support$parent_node_id]
  support$parent_x <- x_by_node[as.character(support$parent_number)]
  support$support_x <- (support$x + support$parent_x) / 2
  support$support_label <- ifelse(
    support$ufboot_status == "missing_support", "NA", support$iqtree_ufboot
  )

  losses <- family_events[family_events$event_type == "L", , drop = FALSE]
  loss_edges <- data.frame()
  if (nrow(losses)) {
    losses$parent_number <- node_number_by_id[losses$edge_parent_node_id]
    losses$child_number <- node_number_by_id[losses$edge_child_node_id]
    if (anyNA(losses$parent_number) || anyNA(losses$child_number)) {
      stop("Loss edge map failed for ", family_id)
    }
    edge_key <- paste(losses$edge_parent_node_id, losses$edge_child_node_id, sep = "->")
    loss_edges <- do.call(rbind, lapply(split(seq_len(nrow(losses)), edge_key), function(index) {
      part <- losses[index, , drop = FALSE]
      parent_number <- part$parent_number[[1]]
      child_number <- part$child_number[[1]]
      data.frame(
        x = (x_by_node[as.character(parent_number)] + x_by_node[as.character(child_number)]) / 2,
        y = annotation$y[match(child_number, annotation$node)],
        n_loss = nrow(part), stringsAsFactors = FALSE
      )
    }))
  }

  domains <- family_domain_hits
  domain_colors <- character()
  if (nrow(domains)) {
    domains$y <- tips$y[match(domains$leaf_id, tips$leaf_id)]
    if (anyNA(domains$y)) stop("CDD hit absent from active tips for ", family_id)
    domains$domain_from <- as.numeric(domains$domain_from)
    domains$domain_to <- as.numeric(domains$domain_to)
    domains$overlap_lane <- as.numeric(domains$overlap_lane)
    domains$domain_name_short <- substr(domains$domain_name, 1, 28)
    levels <- sort(unique(domains$domain_accession))
    palette <- grDevices::hcl.colors(max(3, length(levels)), palette = "Dark 3")
    domain_colors <- setNames(palette[seq_along(levels)], levels)
    domains$domain_color <- unname(domain_colors[domains$domain_accession])
  }
  list(base = base, tips = tips, internal = internal, support = support,
       losses = losses, loss_edges = loss_edges, domains = domains,
       domain_colors = domain_colors)
}

make_plot <- function(family_id, plot_data) {
  family_defs <- definitions[definitions$family_id == family_id, , drop = FALSE]
  family_defs <- family_defs[order(as.integer(family_defs$og_order)), , drop = FALSE]
  family_name <- unique(family_defs$family_name)
  og_colors <- c(setNames(family_defs$og_color, family_defs$original_orthogroup), setNames("#BDBDBD", "NA"))
  tips <- plot_data$tips
  internal <- plot_data$internal
  support <- plot_data$support
  max_tree_x <- max(plot_data$base$data$x)
  expression_x <- max_tree_x + 15.5
  pct_x <- expression_x + 1.6
  domain_start_x <- pct_x + 2.2
  domain_width <- 14.5
  max_protein_length <- max(tips$protein_length_domain)
  tips$domain_xend <- domain_start_x + domain_width * tips$protein_length_domain / max_protein_length
  domains <- plot_data$domains
  if (nrow(domains)) {
    domains$xmin <- domain_start_x + domain_width * domains$domain_from / max_protein_length
    domains$xmax <- domain_start_x + domain_width * domains$domain_to / max_protein_length
    domains$xmid <- (domains$xmin + domains$xmax) / 2
    lane_max <- ave(domains$overlap_lane, domains$leaf_id, FUN = max)
    domains$ycenter <- domains$y +
      (domains$overlap_lane - (lane_max + 1) / 2) * domain_lane_spacing
    domains$ymin <- domains$ycenter - domain_rectangle_half_height
    domains$ymax <- domains$ycenter + domain_rectangle_half_height
  }
  max_y <- max(tips$y)
  observed <- tips[tips$expression_observation_status == "observed", , drop = FALSE]
  no_dataset <- tips[tips$expression_observation_status == "no_expression_dataset", , drop = FALSE]
  gene_not_tested <- tips[tips$expression_observation_status == "gene_not_tested", , drop = FALSE]
  celltype_missing <- tips[tips$expression_observation_status == "celltype_not_present", , drop = FALSE]

  plot <- plot_data$base +
    ggplot2::geom_label(
      data = tips,
      aes(x = x + 0.16, y = y, label = plot_label, fill = original_orthogroup),
      inherit.aes = FALSE, hjust = 0, size = 1.92, label.size = 0.08,
      label.padding = grid::unit(0.08, "lines"), colour = "#111111", alpha = 0.82
    ) +
    ggplot2::scale_fill_manual(values = og_colors, breaks = names(og_colors), labels = ifelse(names(og_colors) == "NA", "External reference", names(og_colors)),
                               name = "Original OG / reference") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      data = internal, aes(x = x, y = y, colour = event_type),
      inherit.aes = FALSE, size = 1.75
    ) +
    ggplot2::scale_colour_manual(
      values = c(S = species_color, D = duplication_color),
      breaks = c("S", "D"), labels = c("Speciation", "Duplication"),
      name = "GeneRax event"
    ) +
    ggplot2::geom_text(
      data = support, aes(x = support_x, y = y + 0.15, label = support_label),
      inherit.aes = FALSE, size = 1.45, colour = "#636363"
    )

  if (nrow(observed)) {
    plot <- plot +
      ggplot2::geom_tile(
        data = observed, aes(x = expression_x, y = y, fill = plot_log2fc),
        inherit.aes = FALSE, width = 0.92, height = 0.78,
        colour = "#525252", linewidth = 0.15
      ) +
      ggplot2::scale_fill_gradient2(
        low = low_color, mid = mid_color, high = high_color, midpoint = 0,
        limits = c(-display_limit, display_limit), oob = scales::squish,
        name = "calicoblast\navg_log2FC",
        breaks = c(-display_limit, 0, display_limit)
      ) +
      ggplot2::geom_point(
        data = observed, aes(x = pct_x, y = y, size = plot_pct, shape = marker_state),
        inherit.aes = FALSE, colour = "#111111", fill = "#111111", stroke = 0.5
      ) +
      ggplot2::scale_size_continuous(
        name = "pct_target", limits = c(0, 1), range = c(1.4, 4.6),
        breaks = c(0.1, 0.5, 0.9)
      ) +
      ggplot2::scale_shape_manual(
        name = "Marker status",
        values = c("marker-positive" = 21, "tested-not-marker" = 1),
        drop = FALSE
      )
  }
  if (nrow(no_dataset)) {
    plot <- plot + ggplot2::geom_tile(
      data = no_dataset, aes(x = expression_x, y = y), inherit.aes = FALSE,
      width = 0.92, height = 0.78, fill = no_dataset_color,
      colour = "#737373", linewidth = 0.15
    )
  }
  if (nrow(gene_not_tested)) {
    plot <- plot + ggplot2::geom_tile(
      data = gene_not_tested, aes(x = expression_x, y = y), inherit.aes = FALSE,
      width = 0.92, height = 0.78, fill = gene_not_tested_color,
      colour = "#252525", linewidth = 0.15
    )
  }
  if (nrow(celltype_missing)) {
    plot <- plot + ggplot2::geom_tile(
      data = celltype_missing, aes(x = expression_x, y = y), inherit.aes = FALSE,
      width = 0.92, height = 0.78, fill = celltype_not_present_color,
      colour = "#737373", linewidth = 0.15
    )
  }
  plot <- plot +
    ggplot2::geom_text(
      data = tips, aes(x = expression_x, y = y, label = fc_label),
      inherit.aes = FALSE, size = 1.55, colour = "#111111", fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = data.frame(x = c(expression_x, pct_x), y = max_y + 0.75,
                        label = c("avg_log2FC", "pct_target")),
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      size = 2.25, fontface = "bold"
    ) +
    ggplot2::geom_segment(
      data = tips, aes(x = domain_start_x, xend = domain_xend, y = y, yend = y),
      inherit.aes = FALSE, colour = config$visualization$domain_no_hit_backbone_color,
      linewidth = 0.50
    ) +
    ggplot2::geom_text(
      data = data.frame(x = domain_start_x + domain_width / 2, y = max_y + 0.75,
                        label = paste0("CDD domains (protein aa; max=", max_protein_length, ")")),
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      size = 2.25, fontface = "bold"
    )
  if (nrow(domains)) {
    plot <- plot +
      ggplot2::geom_rect(
        data = domains, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, fill = domains$domain_color,
        colour = "#252525", linewidth = 0.10
      ) +
      ggplot2::geom_text(
        data = domains,
        aes(x = xmid, y = ymax + domain_label_offset, label = domain_name_short),
        inherit.aes = FALSE, size = domain_label_font_size,
        angle = domain_label_angle, hjust = 0
      )
  }
  if (nrow(plot_data$loss_edges) && isTRUE(config$visualization$multiple_og_show_loss_edge_markers)) {
    plot <- plot + ggplot2::geom_point(
      data = plot_data$loss_edges, aes(x = x, y = y), inherit.aes = FALSE,
      shape = 23, size = 1.9, stroke = 0.4, colour = loss_color, fill = "#FB6A4A"
    )
  }

  plot +
    ggtree::theme_tree2() +
    ggplot2::theme(
      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9), plot.caption = element_text(size = 7.8, hjust = 0),
      plot.margin = margin(14, 38, 12, 12)
    ) +
    coord_cartesian(
      xlim = c(-0.4, domain_start_x + domain_width + 1.8),
      ylim = c(0.4, max_y + 1.05), clip = "off"
    ) +
    labs(
      title = paste0(family_name, " SC11+Nvec all-member tree + calicoblast expression + CDD domains (", family_id, ")"),
      subtitle = paste0(
        nrow(tips), " active tips; observed=", nrow(observed),
        ", no scRNA dataset=", nrow(no_dataset),
        ", gene_not_tested=", nrow(gene_not_tested),
        "; CDD hits=", nrow(domains), "; original OGs=", nrow(family_defs)
      ),
      caption = paste(
        "Tip label: species | gene_id|orthogroup_name|Uniprot_name | original OG; background denotes original Broccoli OG.",
        "Expression is within-species calicoblast-vs-rest avg_log2FC; ND=no expression dataset, NT=gene_not_tested.",
        "Point size=pct_target; filled point=positive marker (log2FC>0.25, pct_target>0.05, adjusted P<0.05).",
        paste0("Color is display-squished to [-", display_limit, ", ", display_limit, "]; source values remain unchanged."),
        "CDD coordinates use the prior full-length-protein CDD.v3.21 results; overlapping domains retain separate lanes.",
        "Branch numbers are IQ-TREE UFBoot for unrooted bipartitions, not S/D confidence.",
        "All family members from SC11 plus Nvec are shown; Nvec has no expression dataset; GeneRax events are conditional annotations and losses are hidden."
      )
    )
}

families <- summaries$orthogroup
width_inches <- as.numeric(tracks$pdf_width_inches)
minimum_height <- as.numeric(tracks$pdf_minimum_height_inches)
height_per_tip <- as.numeric(tracks$pdf_height_per_tip_inches)
manifest <- list()

for (family_id in families) {
  plot_data <- prepare_family(family_id)
  plot <- make_plot(family_id, plot_data)
  height_inches <- max(minimum_height, nrow(plot_data$tips) * height_per_tip + 4)
  pdf_path <- file.path(output_dir, paste0(family_id, ".og_expression_domains.pdf"))
  grDevices::pdf(pdf_path, width = width_inches, height = height_inches,
                 onefile = TRUE, useDingbats = FALSE, paper = "special")
  print(plot)
  grDevices::dev.off()
  pdf_check <- validate_pdf(pdf_path)
  if (pdf_check$status != 0L) stop("Ghostscript validation failed for ", pdf_path)
  tips <- plot_data$tips
  family_defs <- definitions[definitions$family_id == family_id, , drop = FALSE]
  manifest[[length(manifest) + 1L]] <- data.frame(
    family_id = family_id, family_name = unique(family_defs$family_name),
    original_orthogroups = paste(family_defs$original_orthogroup[order(as.integer(family_defs$og_order))], collapse = ";"),
    celltype = "calicoblast", n_tips = nrow(tips),
    n_observed_tips = sum(tips$expression_observation_status == "observed"),
    n_no_expression_dataset_tips = sum(tips$expression_observation_status == "no_expression_dataset"),
    n_gene_not_tested_tips = sum(tips$expression_observation_status == "gene_not_tested"),
    n_domain_hit_tips = sum(tips$domain_status == "domain_hit"),
    n_no_domain_hit_tips = sum(tips$domain_status == "no_domain_hit"),
    n_cdd_domain_hits = nrow(plot_data$domains), n_loss_events = nrow(plot_data$losses),
    pages = 1L, width_inches = width_inches, height_inches = height_inches,
    pdf = normalizePath(pdf_path, mustWork = TRUE), sha256 = sha256_file(pdf_path),
    size_bytes = file.info(pdf_path)$size, ghostscript_status = pdf_check$status,
    status = "PASS", stringsAsFactors = FALSE
  )
}

manifest_df <- do.call(rbind, manifest)
manifest_path <- file.path(output_dir, "family_expression_domain_plot_manifest.tsv")
write.table(manifest_df, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
environment_df <- data.frame(
  component = c("R", "ape", "digest", "ggnewscale", "ggplot2", "ggtree", "scales", "yaml"),
  version = c(R.version.string, vapply(
    c("ape", "digest", "ggnewscale", "ggplot2", "ggtree", "scales", "yaml"),
    function(x) as.character(packageVersion(x)), character(1)
  )), stringsAsFactors = FALSE
)
write.table(environment_df, file.path(output_dir, "plot_environment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
input_paths <- c(
  config = config_path,
  track_config = track_config_path,
  family_definitions = file.path(project_root, config$multiple_og$family_definitions),
  tip_relative_expression = file.path(track_root, "tip_relative_expression.tsv"),
  cdd_domain_hits = file.path(track_root, "cdd_domain_hits.tsv"),
  cdd_sequence_status = file.path(track_root, "cdd_sequence_status.tsv"),
  tip_track_coverage = file.path(track_root, "tip_track_coverage.tsv"),
  reconciliation_nodes = file.path(reconciliation_dir, "reconciliation_nodes.tsv"),
  reconciliation_events = file.path(reconciliation_dir, "reconciliation_events.tsv")
)
input_df <- data.frame(
  role = names(input_paths), path = normalizePath(unname(input_paths), mustWork = TRUE),
  sha256 = vapply(unname(input_paths), sha256_file, character(1)),
  size_bytes = file.info(unname(input_paths))$size, stringsAsFactors = FALSE
)
write.table(input_df, file.path(output_dir, "plot_input_files.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Combined-family expression-domain plotting: PASS (%d single-page PDFs)\n", nrow(manifest_df)))
