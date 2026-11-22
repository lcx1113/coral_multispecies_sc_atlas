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
scope_label <- if (length(args) >= 2) args[[2]] else "sc11_all_members"
track_config_path <- if (length(args) >= 3) args[[3]] else file.path(project_root, "config", "tip_tracks.yaml")
config <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
track_config <- yaml::read_yaml(normalizePath(track_config_path, mustWork = TRUE))
tracks <- track_config$multiple_og_tracks

stopifnot(
  identical(config$multiple_og$preserve_original_broccoli_og, TRUE),
  identical(config$multiple_og$recruit_additional_ogs, FALSE),
  identical(config$visualization$multiple_og_layout, "circular"),
  identical(config$visualization$multiple_og_tip_log2fc_color, TRUE),
  identical(config$visualization$multiple_og_tip_pct_size, TRUE),
  identical(config$visualization$multiple_og_marker_outline, TRUE),
  identical(config$visualization$multiple_og_og_ring_outside_labels, TRUE),
  identical(config$visualization$filter_iqtree_ufboot, FALSE),
  identical(config$visualization$multiple_og_show_loss_edge_markers, FALSE),
  identical(config$reconciliation$method, "GeneRax_EVAL_only")
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

scientific_abbreviation <- function(x) {
  if (x %in% c("Xenia sp.", "Fungia sp.")) return(x)
  parts <- strsplit(trimws(x), "[[:space:]]+")[[1]]
  if (length(parts) < 2L) return(x)
  paste0(substr(parts[[1]], 1L, 1L), ". ", paste(parts[-1L], collapse = " "))
}

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

definitions <- read_tsv(file.path(project_root, config$multiple_og$family_definitions))
definitions <- definitions[toupper(definitions$active) == "TRUE", , drop = FALSE]
display_names <- read_tsv(file.path(project_root, config$inputs$species_display_names))
display_name_by_species <- setNames(display_names$display_name, display_names$tree_species)
track_root <- file.path(project_root, tracks$output_directory)
tip_expression <- read_tsv(file.path(track_root, "tip_relative_expression.tsv"))
reconciliation_dir <- file.path(project_root, "results", "03_reconciliation", scope_label)
nodes <- read_tsv(file.path(reconciliation_dir, "reconciliation_nodes.tsv"))
summaries <- read_tsv(file.path(reconciliation_dir, "reconciliation_summary.tsv"))

stopifnot(
  !anyDuplicated(paste(tip_expression$family_id, tip_expression$leaf_id, sep = "\r")),
  all(tip_expression$celltype == "calicoblast"),
  setequal(unique(tip_expression$family_id), summaries$orthogroup)
)

output_dir <- file.path(project_root, config$multiple_og$plot_output_directory, scope_label)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prepare_family <- function(family_id) {
  family_nodes <- nodes[nodes$orthogroup == family_id, , drop = FALSE]
  family_expression <- tip_expression[tip_expression$family_id == family_id, , drop = FALSE]
  nhx_paths <- unique(family_nodes$reconciled_nhx)
  stopifnot(length(nhx_paths) == 1L, file.exists(nhx_paths))
  tree <- ape::read.tree(nhx_paths)
  descendants <- descendant_tip_labels(tree)
  descendant_keys <- vapply(descendants, function(x) paste(sort(x), collapse = ";"), character(1))
  node_index <- match(descendant_keys, family_nodes$descendant_leaves)
  if (anyNA(node_index) || anyDuplicated(node_index)) {
    stop("Failed to map reconciled nodes for ", family_id)
  }
  annotation <- family_nodes[node_index, , drop = FALSE]
  annotation$node <- seq_len(nrow(annotation))
  base <- ggtree::ggtree(
    tree, layout = "circular", branch.length = "none",
    open.angle = as.numeric(config$visualization$multiple_og_open_angle_degrees),
    linewidth = 0.28, colour = "#525252"
  )
  branch_scale <- as.numeric(config$visualization$multiple_og_branch_length_scale)
  if (!is.finite(branch_scale) || branch_scale <= 0 || branch_scale >= 1) {
    stop("multiple_og_branch_length_scale must lie strictly between 0 and 1")
  }
  base$data$x <- base$data$x * branch_scale
  y_scale <- base$scales$get_scales("y")
  y_scale$limits <- c(0.5, length(tree$tip.label) + 0.5)
  y_scale$expand <- ggplot2::expansion(mult = 0, add = 0)
  coords <- base$data[, c("node", "parent", "x", "y", "angle")]
  annotation <- merge(annotation, coords, by = "node", sort = FALSE, all.x = TRUE)
  annotation <- annotation[match(seq_len(nrow(coords)), annotation$node), , drop = FALSE]
  if (anyNA(annotation$x) || anyNA(annotation$y)) stop("Missing coordinates for ", family_id)

  tips <- annotation[annotation$node_kind == "leaf", , drop = FALSE]
  tips <- tips[match(tree$tip.label, tips$leaf_id), , drop = FALSE]
  expression_index <- match(tips$leaf_id, family_expression$leaf_id)
  if (anyNA(expression_index) || anyDuplicated(expression_index) ||
      nrow(tips) != nrow(family_expression)) {
    stop("Expression does not cover the active tip set for ", family_id)
  }
  expression_columns <- setdiff(names(family_expression), names(tips))
  tips <- cbind(tips, family_expression[expression_index, expression_columns, drop = FALSE])
  if (anyNA(display_name_by_species[tips$leaf_species])) {
    stop("Missing species display name for ", family_id)
  }
  tips$species_fullname <- unname(display_name_by_species[tips$leaf_species])
  tips$species_scientific_abbreviation <- vapply(
    tips$species_fullname, scientific_abbreviation, character(1)
  )
  tips$scientific_abbreviation_label <- tips$species_scientific_abbreviation
  tips$original_identifier_label <- paste(
    tips$leaf_species, tips$canonical_gene_id,
    tips$original_orthogroup, sep = " | "
  )
  external <- tips$leaf_species %in% c("Hvul", "Aaur")
  tips$original_identifier_label[external] <- paste(tips$leaf_species[external], tips$protein_id[external], sep = "|")
  observed <- tips$expression_observation_status == "observed"
  tips$plot_log2fc <- NA_real_
  tips$plot_pct <- NA_real_
  tips$plot_log2fc[observed] <- as.numeric(tips$avg_log2FC[observed])
  tips$plot_pct[observed] <- as.numeric(tips$pct_target[observed])
  if (anyNA(tips$plot_log2fc[observed]) || anyNA(tips$plot_pct[observed])) {
    stop("Observed expression contains nonnumeric values for ", family_id)
  }
  if (any(tips$plot_pct[observed] < 0 | tips$plot_pct[observed] > 1)) {
    stop("pct_target lies outside [0,1] for ", family_id)
  }
  tips$positive_marker_flag <- observed & tips$positive_marker == "TRUE"
  tips$marker_legend <- factor(
    ifelse(tips$positive_marker_flag, "Marker (outer ring)", NA_character_),
    levels = "Marker (outer ring)"
  )
  n_tip <- nrow(tips)
  tips$label_hjust <- 0
  tips$label_angle <- 360 * (tips$y - 0.5) / n_tip
  reverse_text <- tips$label_angle > 90 & tips$label_angle < 270
  tips$label_angle[reverse_text] <- tips$label_angle[reverse_text] + 180
  tips$label_hjust[reverse_text] <- 1

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
  list(base = base, tips = tips, internal = internal, support = support)
}

make_plot <- function(family_id, plot_data, label_variant) {
  family_defs <- definitions[definitions$family_id == family_id, , drop = FALSE]
  family_defs <- family_defs[order(as.integer(family_defs$og_order)), , drop = FALSE]
  family_name <- unique(family_defs$family_name)
  og_colors <- c(setNames(family_defs$og_color, family_defs$original_orthogroup), setNames("#BDBDBD", "NA"))
  tips <- plot_data$tips
  tips$plot_label <- tips[[paste0(label_variant, "_label")]]
  internal <- plot_data$internal
  max_tree_x <- max(plot_data$base$data$x)
  expression_x <- max_tree_x + 0.70
  label_x <- max_tree_x + 1.30
  max_label_chars <- max(nchar(tips$plot_label))
  if (label_variant == "scientific_abbreviation") {
    label_size <- if (nrow(tips) <= 40) 4.30 else if (nrow(tips) <= 80) 3.70 else if (nrow(tips) <= 120) 3.20 else 2.70
  } else {
    label_size <- if (nrow(tips) <= 40) 3.70 else if (nrow(tips) <= 80) 3.05 else if (nrow(tips) <= 120) 2.55 else 2.15
  }
  label_extent_multiplier <- if (label_variant == "original_identifier") {
    as.numeric(config$visualization$multiple_og_original_label_extent_multiplier)
  } else {
    1.0
  }
  label_extension <- max(
    1.10,
    max_label_chars * label_size *
      as.numeric(config$visualization$multiple_og_label_extent_factor) *
      label_extent_multiplier
  )
  og_ring_x <- label_x + label_extension +
    as.numeric(config$visualization$multiple_og_ring_label_gap)
  plot_xmax <- og_ring_x + 0.55
  ring_tips <- tips
  ring_tips$og_ring_y <- ring_tips$y
  ring_tips$og_ring_height <- 1.0
  support_size <- if (nrow(tips) <= 80) 1.45 else 1.05
  observed <- tips[tips$expression_observation_status == "observed", , drop = FALSE]
  missing <- tips[tips$expression_observation_status != "observed", , drop = FALSE]
  markers <- observed[observed$positive_marker_flag, , drop = FALSE]
  n_s <- sum(internal$event_type == "S")
  n_d <- sum(internal$event_type == "D")

  plot <- plot_data$base +
    ggplot2::geom_point(
      data = internal, aes(x = x, y = y, colour = event_type),
      inherit.aes = FALSE, size = 1.6
    ) +
    ggplot2::scale_colour_manual(
      values = c(S = config$visualization$node_color_speciation,
                 D = config$visualization$node_color_duplication),
      breaks = c("S", "D"), labels = c("Speciation", "Duplication"),
      name = "GeneRax event"
    ) +
    ggplot2::geom_tile(
      data = ring_tips,
      aes(x = og_ring_x, y = og_ring_y, height = og_ring_height,
          fill = original_orthogroup),
      inherit.aes = FALSE, width = 0.48,
      colour = NA
    ) +
    ggplot2::scale_fill_manual(
      values = og_colors, breaks = names(og_colors), labels = ifelse(names(og_colors) == "NA", "External reference", names(og_colors)), name = "Original OG / reference"
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_text(
      data = plot_data$support,
      aes(x = support_x, y = y, label = support_label),
      inherit.aes = FALSE, size = support_size, colour = "#3F3FBF"
    )

  if (nrow(observed)) {
    plot <- plot +
      ggplot2::geom_point(
        data = observed,
        aes(x = expression_x, y = y, size = plot_pct, fill = plot_log2fc),
        inherit.aes = FALSE, shape = 21, colour = "#00000000", stroke = 0
      ) +
      ggplot2::scale_size_continuous(
        name = "pct_target", limits = c(0, 1),
        range = c(
          as.numeric(config$visualization$multiple_og_expression_point_size_min),
          as.numeric(config$visualization$multiple_og_expression_point_size_max)
        ),
        breaks = c(0.1, 0.5, 0.9)
      ) +
      ggplot2::scale_fill_gradient2(
        low = config$visualization$expression_log2fc_low_color,
        mid = config$visualization$expression_log2fc_mid_color,
        high = config$visualization$expression_log2fc_high_color,
        midpoint = 0,
        limits = c(-as.numeric(config$visualization$expression_log2fc_display_limit),
                   as.numeric(config$visualization$expression_log2fc_display_limit)),
        oob = scales::squish,
        name = "calicoblast\navg_log2FC"
      )
  }
  if (nrow(markers)) {
    plot <- plot +
      ggplot2::geom_point(
        data = markers,
        aes(x = expression_x, y = y, size = plot_pct, shape = marker_legend),
        inherit.aes = FALSE, fill = NA, colour = "#111111",
        stroke = as.numeric(config$visualization$multiple_og_marker_ring_stroke)
      ) +
      ggplot2::scale_shape_manual(
        name = "Marker status", values = c("Marker (outer ring)" = 21), drop = FALSE
      )
  }
  if (nrow(missing)) {
    plot <- plot + ggplot2::geom_point(
      data = missing, aes(x = expression_x, y = y), inherit.aes = FALSE,
      shape = 21,
      size = as.numeric(config$visualization$multiple_og_missing_point_size),
      fill = config$visualization$expression_gene_not_tested_color,
      colour = "#00000000", stroke = 0
    )
  }

  plot +
    ggplot2::geom_text(
      data = tips,
      aes(x = label_x, y = y, label = plot_label,
          angle = label_angle, hjust = label_hjust),
      inherit.aes = FALSE, size = label_size, colour = "#111111"
    ) +
    ggtree::theme_tree() +
    ggplot2::theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 7.6, hjust = 0),
      plot.margin = margin(12, 12, 12, 12)
    ) +
    ggplot2::xlim(0, plot_xmax) +
    labs(
      title = paste0(family_name, " SC11+Nvec all-member family tree (", family_id, ")"),
      subtitle = paste0(
        nrow(tips), " genes / ", length(unique(tips$leaf_species)), " species; ",
        nrow(family_defs), " original OGs; GeneRax S=", n_s, ", D=", n_d
      ),
      caption = paste(
        "The colored ring outside the tip labels denotes the immutable original Broccoli OG.",
        "Dot fill=calicoblast avg_log2FC; dot size=pct_target; black outer ring=positive marker",
        "(avg_log2FC>0.25, pct_target>0.05, adjusted P<0.05).",
        "Blue branch numbers are IQ-TREE UFBoot for unrooted bipartitions, not S/D confidence.",
        "All family members from the 11 scRNA-seq species plus Nvec are retained; Nvec has no expression dataset; inferred losses are not displayed."
      )
    )
}

families <- summaries$orthogroup
manifest <- list()
width_inches <- as.numeric(config$visualization$multiple_og_pdf_width_inches)
height_inches <- as.numeric(config$visualization$multiple_og_pdf_minimum_height_inches)
label_variants <- unlist(config$visualization$multiple_og_species_label_variants)
if (!setequal(label_variants, c("scientific_abbreviation", "original_identifier"))) {
  stop("Expected exactly the scientific_abbreviation and original_identifier label variants")
}

for (family_id in families) {
  plot_data <- prepare_family(family_id)
  pdf_records <- list()
  for (label_variant in label_variants) {
    plot <- make_plot(family_id, plot_data, label_variant)
    pdf_path <- file.path(output_dir, paste0(family_id, ".circular.", label_variant, ".pdf"))
    grDevices::pdf(pdf_path, width = width_inches, height = height_inches,
                   onefile = TRUE, useDingbats = FALSE, paper = "special")
    print(plot)
    grDevices::dev.off()
    pdf_check <- validate_pdf(pdf_path)
    if (pdf_check$status != 0L) stop("Ghostscript validation failed for ", pdf_path)
    pdf_records[[label_variant]] <- list(
      path = normalizePath(pdf_path, mustWork = TRUE),
      sha256 = sha256_file(pdf_path), size = file.info(pdf_path)$size,
      ghostscript_status = pdf_check$status
    )
  }
  family_defs <- definitions[definitions$family_id == family_id, , drop = FALSE]
  manifest[[length(manifest) + 1L]] <- data.frame(
    family_id = family_id,
    family_name = unique(family_defs$family_name),
    original_orthogroups = paste(
      family_defs$original_orthogroup[order(as.integer(family_defs$og_order))], collapse = ";"
    ),
    n_original_orthogroups = nrow(family_defs),
    n_tips = nrow(plot_data$tips),
    n_species = length(unique(plot_data$tips$leaf_species)),
    n_observed_tips = sum(plot_data$tips$expression_observation_status == "observed"),
    n_marker_tips = sum(plot_data$tips$positive_marker_flag),
    n_speciation_nodes = sum(plot_data$internal$event_type == "S"),
    n_duplication_nodes = sum(plot_data$internal$event_type == "D"),
    branch_length_scale = as.numeric(config$visualization$multiple_og_branch_length_scale),
    og_ring_outside_labels = TRUE,
    loss_displayed = FALSE,
    pages_per_pdf = 1L,
    width_inches = width_inches,
    height_inches = height_inches,
    scientific_abbreviation_pdf = pdf_records$scientific_abbreviation$path,
    scientific_abbreviation_sha256 = pdf_records$scientific_abbreviation$sha256,
    scientific_abbreviation_size_bytes = pdf_records$scientific_abbreviation$size,
    scientific_abbreviation_ghostscript_status = pdf_records$scientific_abbreviation$ghostscript_status,
    original_identifier_pdf = pdf_records$original_identifier$path,
    original_identifier_sha256 = pdf_records$original_identifier$sha256,
    original_identifier_size_bytes = pdf_records$original_identifier$size,
    original_identifier_ghostscript_status = pdf_records$original_identifier$ghostscript_status,
    status = "PASS", stringsAsFactors = FALSE
  )
}

manifest_df <- do.call(rbind, manifest)
write.table(manifest_df, file.path(output_dir, "family_tree_plot_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
environment_df <- data.frame(
  component = c("R", "ape", "digest", "ggnewscale", "ggplot2", "ggtree", "scales", "yaml"),
  version = c(
    R.version.string, as.character(packageVersion("ape")), as.character(packageVersion("digest")),
    as.character(packageVersion("ggnewscale")), as.character(packageVersion("ggplot2")),
    as.character(packageVersion("ggtree")), as.character(packageVersion("scales")),
    as.character(packageVersion("yaml"))
  ),
  stringsAsFactors = FALSE
)
write.table(environment_df, file.path(output_dir, "plot_environment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("SC11 circular family plotting: PASS (%d families; %d PDFs)\n",
            nrow(manifest_df), 2L * nrow(manifest_df)))
