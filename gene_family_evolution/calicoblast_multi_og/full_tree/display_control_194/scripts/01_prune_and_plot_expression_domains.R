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
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", full_args, value = TRUE)[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
config_path <- if (length(args)) args[[1]] else file.path(project_root, "config", "analysis_config.yaml")
config_path <- normalizePath(config_path, mustWork = TRUE)
config <- yaml::read_yaml(config_path)
source_root <- normalizePath(config$source$project, mustWork = TRUE)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

validate_pdf <- function(path) {
  output <- suppressWarnings(system2(
    config$visualization$ghostscript,
    c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", shQuote(path)),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), message = paste(output, collapse = " | "))
}

split_values <- function(value) {
  if (is.na(value) || !nzchar(value)) character() else strsplit(value, ";", fixed = TRUE)[[1]]
}

descendant_tip_labels <- function(tree) {
  n_tip <- length(tree$tip.label)
  n_total <- n_tip + tree$Nnode
  children <- split(tree$edge[, 2], tree$edge[, 1])
  cache <- vector("list", n_total)
  descend <- function(node) {
    if (!is.null(cache[[node]])) return(cache[[node]])
    answer <- if (node <= n_tip) tree$tip.label[[node]] else {
      sort(unique(unlist(lapply(children[[as.character(node)]], descend), use.names = FALSE)))
    }
    cache[[node]] <<- answer
    answer
  }
  answer <- lapply(seq_len(n_total), descend)
  names(answer) <- as.character(seq_len(n_total))
  answer
}

family_file <- if (is.null(config$family_definitions)) "config/family_sets.tsv" else config$family_definitions
definitions <- read_tsv(file.path(project_root, family_file))
definitions <- definitions[toupper(definitions$active) == "TRUE", , drop = FALSE]
expression_path <- file.path(source_root, config$source$expression)
domain_hits_path <- file.path(source_root, config$source$cdd_domain_hits)
domain_status_path <- file.path(source_root, config$source$cdd_sequence_status)
nodes_path <- file.path(source_root, config$source$reconciliation_nodes)
summary_path <- file.path(source_root, config$source$reconciliation_summary)
selection_path <- file.path(project_root, config$outputs$selection_dir, "tip_selection_audit.tsv")
expression <- read_tsv(expression_path)
domains_all <- read_tsv(domain_hits_path)
domain_status_all <- read_tsv(domain_status_path)
nodes_all <- read_tsv(nodes_path)
summaries <- read_tsv(summary_path)
selection <- read_tsv(selection_path)
selection <- selection[selection$retain_tip == "TRUE", , drop = FALSE]

stopifnot(
  !anyDuplicated(paste(selection$family_id, selection$leaf_id, sep = "\r")),
  !anyDuplicated(paste(expression$family_id, expression$leaf_id, sep = "\r")),
  !anyDuplicated(paste(domain_status_all$family_id, domain_status_all$leaf_id, sep = "\r")),
  setequal(unique(selection$family_id), unique(definitions$family_id))
)

tree_dir <- file.path(project_root, config$outputs$tree_dir)
dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)

map_pruned_nodes <- function(full_tree, pruned_tree, family_nodes, retained_tips) {
  full_desc <- descendant_tip_labels(full_tree)
  pruned_desc <- descendant_tip_labels(pruned_tree)
  full_key_to_row <- setNames(seq_len(nrow(family_nodes)), family_nodes$descendant_leaves)
  full_keys <- vapply(full_desc, function(x) paste(x, collapse = ";"), character(1))
  full_index <- unname(full_key_to_row[full_keys])
  if (anyNA(full_index)) stop("Could not map source tree descendant sets to reconciliation nodes")
  full_rows <- family_nodes[full_index, , drop = FALSE]
  rownames(full_rows) <- NULL
  node_desc_by_id <- setNames(lapply(seq_len(nrow(family_nodes)), function(i) split_values(family_nodes$descendant_leaves[[i]])), family_nodes$node_id)

  candidates <- list()
  full_internal <- full_rows[full_rows$node_kind == "internal", , drop = FALSE]
  for (i in seq_len(nrow(full_internal))) {
    row <- full_internal[i, , drop = FALSE]
    kept <- sort(intersect(split_values(row$descendant_leaves), retained_tips))
    if (length(kept) < 2) next
    child_ids <- split_values(row$child_node_ids)
    informative_children <- sum(vapply(child_ids, function(id) {
      length(intersect(node_desc_by_id[[id]], retained_tips)) > 0
    }, logical(1)))
    if (informative_children < 2) next
    key <- paste(kept, collapse = ";")
    candidates[[key]] <- c(candidates[[key]], row$node_id)
  }

  rows <- list()
  n_tip <- length(pruned_tree$tip.label)
  for (node_number in seq.int(n_tip + 1L, n_tip + pruned_tree$Nnode)) {
    leaves <- sort(pruned_desc[[as.character(node_number)]])
    key <- paste(leaves, collapse = ";")
    ids <- unique(candidates[[key]])
    if (is.null(ids)) ids <- character()
    status <- if (length(ids) == 1L) "mapped_unique_informative_original_node" else if (!length(ids)) "unmapped" else "ambiguous_collapsed_original_nodes"
    original <- if (length(ids) == 1L) family_nodes[match(ids, family_nodes$node_id), , drop = FALSE] else NULL
    rows[[length(rows) + 1L]] <- data.frame(
      pruned_node_number = node_number,
      n_retained_descendant_leaves = length(leaves),
      retained_descendant_leaves = key,
      mapping_status = status,
      n_original_candidates = length(ids),
      original_candidate_node_ids = paste(ids, collapse = ";"),
      inherited_node_id = if (is.null(original)) "" else original$node_id,
      inherited_event_type = if (is.null(original)) "" else original$event_type,
      inherited_species_location = if (is.null(original)) "" else original$species_location,
      inherited_split_id = if (is.null(original)) "" else original$split_id,
      inherited_iqtree_ufboot = if (is.null(original)) "" else original$iqtree_ufboot,
      inherited_ufboot_status = if (is.null(original)) "" else original$ufboot_status,
      inherited_support_scope = if (is.null(original)) "" else original$support_scope,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

prepare_family <- function(family_id) {
  family_selection <- selection[selection$family_id == family_id, , drop = FALSE]
  family_expression <- expression[expression$family_id == family_id, , drop = FALSE]
  family_status <- domain_status_all[domain_status_all$family_id == family_id, , drop = FALSE]
  family_domains <- domains_all[domains_all$family_id == family_id, , drop = FALSE]
  family_nodes <- nodes_all[nodes_all$orthogroup == family_id, , drop = FALSE]
  nhx_paths <- unique(family_nodes$reconciled_nhx)
  nhx_paths <- ifelse(grepl("^/", nhx_paths), nhx_paths, file.path(source_root, nhx_paths))
  stopifnot(length(nhx_paths) == 1L, file.exists(nhx_paths))
  full_tree <- ape::read.tree(nhx_paths)
  retained <- family_selection$leaf_id
  if (!all(retained %in% full_tree$tip.label)) stop("Selected tip absent from tree: ", family_id)
  if (length(retained) < 2L) stop("Fewer than two retained tips: ", family_id)
  pruned <- ape::keep.tip(full_tree, retained)
  mapping <- map_pruned_nodes(full_tree, pruned, family_nodes, retained)
  pruned$node.label <- mapping$inherited_node_id[match(seq.int(length(pruned$tip.label) + 1L, length(pruned$tip.label) + pruned$Nnode), mapping$pruned_node_number)]
  tree_path <- file.path(tree_dir, paste0(family_id, ".marker_or_top_log2fc.pruned.nwk"))
  ape::write.tree(pruned, file = tree_path)
  mapping$family_id <- family_id
  mapping$source_tree <- nhx_paths
  mapping$pruned_tree <- normalizePath(tree_path, mustWork = TRUE)
  mapping_path <- file.path(tree_dir, paste0(family_id, ".pruned_node_annotation_mapping.tsv"))
  write.table(mapping, mapping_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  base <- ggtree::ggtree(pruned, layout = "rectangular", branch.length = "none", linewidth = 0.38, colour = "#4D4D4D")
  coords <- base$data[, c("node", "parent", "x", "y")]
  tips <- coords[coords$node <= length(pruned$tip.label), , drop = FALSE]
  tips$leaf_id <- pruned$tip.label[tips$node]
  tips <- tips[match(pruned$tip.label, tips$leaf_id), , drop = FALSE]
  sel_idx <- match(tips$leaf_id, family_selection$leaf_id)
  exp_idx <- match(tips$leaf_id, family_expression$leaf_id)
  status_idx <- match(tips$leaf_id, family_status$leaf_id)
  if (anyNA(sel_idx) || anyNA(exp_idx) || anyNA(status_idx)) stop("Track coverage failure: ", family_id)
  tips <- cbind(tips, family_selection[sel_idx, setdiff(names(family_selection), names(tips)), drop = FALSE])
  tips$protein_length_domain <- as.numeric(family_status$protein_length[status_idx])
  tips$domain_status <- family_status$domain_status[status_idx]
  tips$plot_label <- paste(tips$tip_display_label, tips$original_orthogroup, sep = " | ")
  tips$clean_label <- paste(tips$tree_species, tips$canonical_gene_id, sep = " | ")
  external <- tips$tree_species %in% c("Hvul", "Aaur")
  tips$plot_label[external] <- paste(tips$tree_species[external], tips$protein_id[external], sep = "|")
  tips$clean_label[external] <- tips$plot_label[external]
  tips$plot_log2fc <- NA_real_
  tips$plot_pct <- NA_real_
  observed <- tips$expression_observation_status == "observed"
  tips$plot_log2fc[observed] <- as.numeric(tips$avg_log2FC[observed])
  tips$plot_pct[observed] <- as.numeric(tips$pct_target[observed])
  tips$fc_label <- ifelse(
    observed, sprintf("%.2f", tips$plot_log2fc),
    ifelse(tips$expression_observation_status == "no_expression_dataset", "ND",
    ifelse(tips$expression_observation_status == "gene_not_tested", "NT",
    ifelse(tips$expression_observation_status == "celltype_not_present", "CT", "NA")))
  )
  tips$selection_state <- ifelse(
    tips$positive_marker == "TRUE", "marker",
    ifelse(grepl("retain_all_external_reference_no_expression_filter", tips$selection_reason),
           "Hvul/Aaur retained without expression filtering",
    ifelse(grepl("species_fallback_no_observed", tips$selection_reason),
           "no observed value fallback", "species top log2FC non-marker"))
  )
  tips$selection_state <- factor(
    tips$selection_state,
    levels = c("marker", "species top log2FC non-marker", "no observed value fallback",
               "Hvul/Aaur retained without expression filtering")
  )

  mapping <- merge(mapping, coords, by.x = "pruned_node_number", by.y = "node", all.x = TRUE, sort = FALSE)
  mapping <- mapping[match(seq.int(length(pruned$tip.label) + 1L, length(pruned$tip.label) + pruned$Nnode), mapping$pruned_node_number), , drop = FALSE]
  internal <- mapping[mapping$mapping_status == "mapped_unique_informative_original_node" & mapping$inherited_event_type %in% c("S", "D"), , drop = FALSE]
  internal$event_type <- factor(internal$inherited_event_type, levels = c("S", "D"))
  x_by_node <- setNames(coords$x, coords$node)
  support <- internal[nzchar(internal$inherited_split_id), , drop = FALSE]
  if (nrow(support)) {
    support$parent_x <- x_by_node[as.character(support$parent)]
    support$support_x <- (support$x + support$parent_x) / 2
    support$support_label <- ifelse(support$inherited_ufboot_status == "missing_support", "NA", support$inherited_iqtree_ufboot)
    support <- support[!is.na(support$support_x), , drop = FALSE]
  }

  domains <- family_domains[family_domains$leaf_id %in% retained, , drop = FALSE]
  domain_colors <- character()
  if (nrow(domains)) {
    domains$y <- tips$y[match(domains$leaf_id, tips$leaf_id)]
    domains$domain_from <- as.numeric(domains$domain_from)
    domains$domain_to <- as.numeric(domains$domain_to)
    domains$overlap_lane <- as.numeric(domains$overlap_lane)
    domains$domain_name_short <- substr(domains$domain_name, 1, 28)
    levels <- sort(unique(domains$domain_accession))
    palette <- grDevices::hcl.colors(max(3, length(levels)), palette = "Dark 3")
    domain_colors <- setNames(palette[seq_along(levels)], levels)
    domains$domain_color <- unname(domain_colors[domains$domain_accession])
  }
  list(base = base, tips = tips, internal = internal, support = support, domains = domains,
       mapping = mapping, tree_path = tree_path, full_tree = full_tree, pruned = pruned)
}

make_article_plot <- function(family_id, d) {
  defs <- definitions[definitions$family_id == family_id, , drop = FALSE]
  defs <- defs[order(as.integer(defs$og_order)), , drop = FALSE]
  og_colors <- c(setNames(defs$og_color, defs$original_orthogroup), setNames("#BDBDBD", "NA"))
  branch_scale <- as.numeric(config$visualization$article_branch_x_scale)
  base <- d$base
  base$data$x <- base$data$x * branch_scale
  tips <- d$tips
  tips$x <- tips$x * branch_scale

  uniprot_label <- tips$Uniprot_name
  uniprot_label[is.na(uniprot_label) | !nzchar(uniprot_label)] <- "NA"
  gene_label <- tips$gene_id
  gene_label[is.na(gene_label) | !nzchar(gene_label)] <- "NA"
  tips$article_label <- paste(tips$tree_species, uniprot_label, gene_label, sep = "|")
  external <- tips$tree_species %in% c("Hvul", "Aaur")
  tips$article_label[external] <- paste(tips$tree_species[external], tips$protein_id[external], sep = "|")
  tips$tip_lineage <- ifelse(tips$tree_species == "Hcoe", "Hcoe",
                             ifelse(tips$tree_species == "Xesp", "Xesp", "Other coral"))
  tips$tip_lineage <- factor(tips$tip_lineage, levels = c("Hcoe", "Xesp", "Other coral"))

  pruned_desc <- descendant_tip_labels(d$pruned)
  og_by_leaf <- setNames(tips$original_orthogroup, tips$leaf_id)
  internal_og_counts <- vapply(d$mapping$pruned_node_number, function(node) {
    values <- unique(og_by_leaf[pruned_desc[[as.character(node)]]])
    real_ogs <- values[!is.na(values) & values != "NA"]
    as.integer(length(real_ogs) + as.integer("NA" %in% values && length(real_ogs) > 0))
  }, integer(1))
  cross_node_numbers <- d$mapping$pruned_node_number[internal_og_counts > 1L]
  cross_internal <- d$internal[d$internal$pruned_node_number %in% cross_node_numbers, , drop = FALSE]
  if (nrow(cross_internal)) cross_internal$x <- cross_internal$x * branch_scale
  cross_support <- d$support[d$support$pruned_node_number %in% cross_node_numbers, , drop = FALSE]
  if (nrow(cross_support)) {
    support_numeric <- suppressWarnings(as.numeric(cross_support$inherited_iqtree_ufboot))
    cross_support <- cross_support[
      !is.na(support_numeric) & support_numeric >= as.numeric(config$visualization$article_cross_og_support_minimum),
      , drop = FALSE
    ]
    cross_support$support_x <- cross_support$support_x * branch_scale
  }

  max_tree_x <- max(base$data$x)
  og_x <- max_tree_x + 0.36
  label_x <- og_x + 0.30
  expression_x <- label_x + 5.75
  domain_start_x <- expression_x + 1.15
  domain_width <- 8.4
  max_protein_length <- max(tips$protein_length_domain)
  tips$domain_xend <- domain_start_x + domain_width * tips$protein_length_domain / max_protein_length
  observed <- tips[tips$expression_observation_status == "observed", , drop = FALSE]
  missing <- tips[tips$expression_observation_status != "observed", , drop = FALSE]
  marker_tips <- observed[observed$positive_marker == "TRUE", , drop = FALSE]
  max_y <- max(tips$y)
  display_limit <- as.numeric(config$visualization$expression_display_limit)

  domains <- d$domains
  domain_colors <- character()
  if (nrow(domains)) {
    domains$xmin <- domain_start_x + domain_width * domains$domain_from / max_protein_length
    domains$xmax <- domain_start_x + domain_width * domains$domain_to / max_protein_length
    lane_max <- ave(domains$overlap_lane, domains$leaf_id, FUN = max)
    domains$ycenter <- domains$y + (domains$overlap_lane - (lane_max + 1) / 2) * 0.13
    half_height <- as.numeric(config$visualization$article_domain_height) / 2
    domains$ymin <- domains$ycenter - half_height
    domains$ymax <- domains$ycenter + half_height
    accession_name <- unique(domains[, c("domain_accession", "domain_name")])
    duplicated_name <- duplicated(accession_name$domain_name) | duplicated(accession_name$domain_name, fromLast = TRUE)
    accession_name$domain_key <- accession_name$domain_name
    accession_name$domain_key[duplicated_name] <- paste0(
      accession_name$domain_name[duplicated_name], " (", accession_name$domain_accession[duplicated_name], ")"
    )
    key_by_accession <- setNames(accession_name$domain_key, accession_name$domain_accession)
    domains$domain_key <- unname(key_by_accession[domains$domain_accession])
    color_by_accession <- setNames(domains$domain_color, domains$domain_accession)
    domain_colors <- setNames(unname(color_by_accession[accession_name$domain_accession]), accession_name$domain_key)
  }

  p <- base
  if (nrow(cross_internal)) {
    p <- p +
      geom_point(data = cross_internal, aes(x = x, y = y, colour = event_type), inherit.aes = FALSE,
                 size = as.numeric(config$visualization$article_event_point_size), stroke = 0.2) +
      scale_colour_manual(
        values = c(S = config$visualization$node_color_speciation, D = config$visualization$node_color_duplication),
        breaks = c("S", "D"), labels = c("Speciation", "Duplication"),
        name = "OG / reference boundary node"
      )
  }
  if (nrow(cross_support)) {
    p <- p + geom_text(
      data = cross_support, aes(x = support_x, y = y + 0.17, label = support_label),
      inherit.aes = FALSE, size = as.numeric(config$visualization$article_support_font_size),
      fontface = "bold", colour = "#505050"
    )
  }
  p <- p +
    ggnewscale::new_scale_colour() +
    geom_tile(data = tips, aes(x = og_x, y = y, fill = original_orthogroup), inherit.aes = FALSE,
              width = 0.22, height = 0.92, colour = "#FFFFFF", linewidth = 0.10) +
    scale_fill_manual(values = og_colors, breaks = names(og_colors), labels = ifelse(names(og_colors) == "NA", "External reference", names(og_colors)), name = "Original OG / reference",
                      guide = guide_legend(order = 1, nrow = 1, byrow = TRUE,
                                           override.aes = list(colour = NA))) +
    geom_text(data = tips, aes(x = label_x, y = y, label = article_label, colour = tip_lineage),
              inherit.aes = FALSE, hjust = 0,
              size = as.numeric(config$visualization$article_tip_label_font_size),
              fontface = ifelse(tips$tree_species %in% c("Hcoe", "Xesp"), "bold", "plain")) +
    scale_colour_manual(
      values = c(Hcoe = config$visualization$article_hcoe_label_color,
                 Xesp = config$visualization$article_xesp_label_color,
                 `Other coral` = config$visualization$article_other_label_color),
      breaks = c("Hcoe", "Xesp", "Other coral"), guide = "none"
    ) +
    ggnewscale::new_scale_fill()
  if (nrow(observed)) {
    p <- p +
      geom_point(data = observed, aes(x = expression_x, y = y, size = plot_pct, fill = plot_log2fc),
                 inherit.aes = FALSE, shape = 21, colour = "transparent", stroke = 0) +
      scale_fill_gradient2(
        low = config$visualization$expression_low_color,
        mid = config$visualization$expression_mid_color,
        high = config$visualization$expression_high_color,
        midpoint = 0, limits = c(-display_limit, display_limit), oob = scales::squish,
        name = "calicoblast avg_log2FC", breaks = c(-display_limit, 0, display_limit),
        guide = guide_colorbar(order = 3, title.position = "top", barwidth = grid::unit(2.0, "cm"))
      ) +
      scale_size_continuous(
        name = "pct_target", limits = c(0, 1),
        range = c(as.numeric(config$visualization$article_expression_point_size_min),
                  as.numeric(config$visualization$article_expression_point_size_max)),
        breaks = c(0.1, 0.5, 0.9), guide = guide_legend(order = 4, nrow = 1)
      )
  }
  if (nrow(marker_tips)) {
    p <- p + geom_point(
      data = marker_tips, aes(x = expression_x, y = y, size = plot_pct),
      inherit.aes = FALSE, shape = 21, fill = NA, colour = "#111111",
      stroke = as.numeric(config$visualization$article_marker_ring_stroke)
    )
  }
  if (nrow(missing)) {
    p <- p + geom_point(data = missing, aes(x = expression_x, y = y), inherit.aes = FALSE,
                        shape = 4, size = 4.0, stroke = 1.0, colour = "#666666")
  }
  p <- p +
    geom_segment(data = tips, aes(x = domain_start_x, xend = domain_xend, y = y, yend = y),
                 inherit.aes = FALSE, colour = "#A6A6A6", linewidth = 0.55)
  if (nrow(domains)) {
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_rect(data = domains, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = domain_key),
                inherit.aes = FALSE, colour = NA) +
      scale_fill_manual(values = domain_colors, name = "CDD domain",
                        guide = guide_legend(order = 5, nrow = if (is.null(config$visualization$article_domain_legend_max_columns)) 2 else max(2, ceiling(length(domain_colors) / config$visualization$article_domain_legend_max_columns)), byrow = TRUE,
                                             override.aes = list(colour = NA)))
  }
  header_y <- max_y + 0.72
  p <- p +
    geom_text(data = data.frame(
      x = c(og_x, expression_x, domain_start_x + domain_width / 2),
      y = header_y, label = c("OG", "Expression", "CDD domain architecture")
    ), aes(x = x, y = y, label = label), inherit.aes = FALSE,
    size = as.numeric(config$visualization$article_header_font_size), fontface = "bold")

  marker_n <- sum(tips$positive_marker == "TRUE")
  hcoe_n <- sum(tips$tree_species == "Hcoe")
  xesp_n <- sum(tips$tree_species == "Xesp")
  plot <- p + ggtree::theme_tree2() +
    theme(
      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
      legend.position = "top", legend.box = "vertical", legend.box.just = "left",
      legend.title = element_text(size = as.numeric(config$visualization$article_legend_title_size), face = "bold"),
      legend.text = element_text(size = as.numeric(config$visualization$article_legend_text_size)),
      legend.key.height = grid::unit(0.42, "cm"), legend.key.width = grid::unit(0.55, "cm"),
      legend.spacing.y = grid::unit(0.10, "cm"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.2, hjust = 0),
      plot.margin = margin(8, 18, 10, 8)
    ) +
    coord_cartesian(
      xlim = c(-0.25, domain_start_x + domain_width + 0.45),
      ylim = c(0.35, max_y + 1.0), clip = "off"
    ) +
    labs(
      title = paste0(gsub("_", " ", unique(defs$family_name)), " family: OG deployment in calicoblasts"),
      subtitle = paste0(nrow(tips), " retained tips from ", nrow(defs), " OGs; markers=", marker_n,
                        "; Hcoe tips=", hcoe_n, ", Xesp tips=", xesp_n,
                        "; annotated boundary nodes=", nrow(cross_internal)),
      caption = paste(
        "Tree branches are compressed for display. S/D and UFBoot are shown at cross-OG or reference/Anthozoa boundary nodes; pure-OG internal events are omitted.",
        "Hcoe labels are blue and Xesp labels are green. Expression point fill=avg_log2FC, size=pct_target, heavy ring=marker; x=no observed expression value.",
        "All Hvul/Aaur proteins are retained without expression filtering; Hvul/Aaur x marks indicate no scRNA-seq dataset, not zero expression.",
        "CDD hit names are represented by the legend rather than printed on each protein. Topology and reconciliation are inherited from the newly inferred complete SC11+Hvul+Aaur analysis; SC11 CDD calls are reused for identical proteins; external reference CDD calls were newly computed."
      )
    )
  label_audit <- tips[, c("family_id", "original_orthogroup", "leaf_id", "tree_species",
                          "Uniprot_name", "gene_id", "protein_id", "article_label"), drop = FALSE]
  list(plot = plot, n_cross_og_nodes = nrow(cross_internal),
       n_cross_og_support_labels = nrow(cross_support), n_domain_legend_entries = length(domain_colors),
       tip_labels = label_audit)
}

families <- summaries$orthogroup[summaries$orthogroup %in% unique(definitions$family_id)]
all_mappings <- list()
article_manifest <- list()
article_tip_labels <- list()
article_dir <- file.path(project_root, config$outputs$article_figure_dir)
dir.create(article_dir, recursive = TRUE, showWarnings = FALSE)
for (family_id in families) {
  d <- prepare_family(family_id)
  article <- make_article_plot(family_id, d)
  article$plot$labels$caption <- paste(strwrap(article$plot$labels$caption, width = 185), collapse = "\n")
  article_width <- as.numeric(config$visualization$article_pdf_width_inches)
  article_height <- max(
    as.numeric(config$visualization$article_pdf_minimum_height_inches),
    nrow(d$tips) * as.numeric(config$visualization$article_pdf_height_per_tip_inches) + 3.5
  )
  article_pdf_path <- file.path(article_dir, paste0(family_id, ".article_final.pdf"))
  grDevices::pdf(article_pdf_path, width = article_width, height = article_height,
                 onefile = TRUE, useDingbats = FALSE, paper = "special")
  print(article$plot)
  grDevices::dev.off()
  article_check <- validate_pdf(article_pdf_path)
  if (article_check$status != 0L) stop("Ghostscript validation failed: ", article_pdf_path, " ", article_check$message)
  defs <- definitions[definitions$family_id == family_id, , drop = FALSE]
  article_manifest[[length(article_manifest) + 1L]] <- data.frame(
    n_source_tips = length(d$full_tree$tip.label),
    n_species_top_nonmarker_tips = sum(d$tips$positive_marker != "TRUE" & grepl("species_top_observed", d$tips$selection_reason)),
    n_no_observed_fallback_tips = sum(grepl("species_fallback_no_observed", d$tips$selection_reason)),
    n_domain_hit_tips = sum(d$tips$domain_status == "domain_hit"),
    n_cdd_domain_hits = nrow(d$domains),
    n_unique_inherited_internal_nodes = sum(d$mapping$mapping_status == "mapped_unique_informative_original_node"),
    n_hidden_ambiguous_or_unmapped_internal_nodes = sum(d$mapping$mapping_status != "mapped_unique_informative_original_node"),
    tree = normalizePath(d$tree_path, mustWork = TRUE),
    tree_sha256 = sha256_file(d$tree_path),
    family_id = family_id,
    family_name = unique(defs$family_name),
    original_orthogroups = paste(defs$original_orthogroup[order(as.integer(defs$og_order))], collapse = ";"),
    n_retained_tips = nrow(d$tips),
    n_marker_tips = sum(d$tips$positive_marker == "TRUE"),
    n_hcoe_tips = sum(d$tips$tree_species == "Hcoe"),
    n_xesp_tips = sum(d$tips$tree_species == "Xesp"),
    n_cross_og_annotated_nodes = article$n_cross_og_nodes,
    n_cross_og_support_labels = article$n_cross_og_support_labels,
    n_domain_legend_entries = article$n_domain_legend_entries,
    n_pure_og_internal_events_displayed = 0L,
    branch_x_scale = as.numeric(config$visualization$article_branch_x_scale),
    width_inches = article_width,
    height_inches = article_height,
    pdf = normalizePath(article_pdf_path, mustWork = TRUE),
    sha256 = sha256_file(article_pdf_path),
    size_bytes = file.info(article_pdf_path)$size,
    ghostscript_status = article_check$status,
    status = "PASS",
    stringsAsFactors = FALSE
  )
  article_tip_labels[[length(article_tip_labels) + 1L]] <- article$tip_labels
  all_mappings[[length(all_mappings) + 1L]] <- d$mapping
}

mapping_df <- do.call(rbind, all_mappings)
article_manifest_df <- do.call(rbind, article_manifest)
article_tip_labels_df <- do.call(rbind, article_tip_labels)
write.table(mapping_df, file.path(tree_dir, "all_families.pruned_node_annotation_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(article_manifest_df, file.path(article_dir, "article_final_figure_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(article_tip_labels_df, file.path(article_dir, "article_final_tip_labels.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
article_contract <- data.frame(
  property = c(
    "source_tip_selection", "source_pruned_tree", "article_tip_label", "tree_branch_geometry", "displayed_internal_events",
    "displayed_support", "OG_encoding", "Hcoe_tip_label", "Xesp_tip_label",
    "expression_encoding", "marker_encoding", "missing_expression_encoding", "CDD_encoding"
  ),
  value = c(
    "current_tip_selection_audit_retained_rows", "pruned_from_complete_reconciled_tree", "species|Uniprot_name|gene_id_missing_values_shown_as_NA", "x_scaled_cladogram_display_only",
    "GeneRax_S_or_D_only_when_descendants_span_multiple_original_OGs",
    "inherited_IQ_TREE_UFBoot_only_at_cross_OG_nodes", "colored_strip",
    config$visualization$article_hcoe_label_color, config$visualization$article_xesp_label_color,
    "point_fill_avg_log2FC_and_size_pct_target", "heavy_outer_ring", "x_symbol",
    "full_length_coordinate_blocks_with_names_in_legend_only"
  ),
  stringsAsFactors = FALSE
)
write.table(article_contract, file.path(article_dir, "article_final_design_contract.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
env <- data.frame(
  component = c("R", "ape", "digest", "ggnewscale", "ggplot2", "ggtree", "scales", "yaml"),
  version = c(R.version.string, vapply(c("ape", "digest", "ggnewscale", "ggplot2", "ggtree", "scales", "yaml"),
                                     function(x) as.character(packageVersion(x)), character(1))),
  stringsAsFactors = FALSE
)
write.table(env, file.path(article_dir, "plot_environment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
inputs <- c(config = config_path, family_definitions = file.path(project_root, family_file),
            selection = selection_path, expression = expression_path, cdd_domain_hits = domain_hits_path,
            cdd_sequence_status = domain_status_path, reconciliation_nodes = nodes_path, reconciliation_summary = summary_path)
input_df <- data.frame(role = names(inputs), path = normalizePath(unname(inputs), mustWork = TRUE),
                       sha256 = vapply(unname(inputs), sha256_file, character(1)),
                       size_bytes = file.info(unname(inputs))$size, stringsAsFactors = FALSE)
write.table(input_df, file.path(article_dir, "plot_input_files.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Article-final plotting: PASS (%d compact PDFs)\n", nrow(article_manifest_df)))
