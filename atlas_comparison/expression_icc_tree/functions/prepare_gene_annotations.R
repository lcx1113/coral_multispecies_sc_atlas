# Adapted from prepare_gene_annote_260428_reanotated.Rmd.
prepare_gene_annotations <- function(species_ids, reference_dir, work_dir = ".") {
    read_tsv <- function(path) {
        read.delim(path, quote = "", stringsAsFactors = FALSE)
    }
    members <- read_tsv(file.path(reference_dir, "orthogroups", "orthogroup_conservation_final.csv"))
    annotations <- read_tsv(file.path(reference_dir, "orthogroups", "orthogroup_final_annotations.csv"))
    stopifnot(!anyDuplicated(annotations$orthogroup))

    most_frequent <- function(values) {
        if (!length(values)) return(NA_character_)
        names(sort(table(values), decreasing = TRUE))[1]
    }
    repeated <- members$gene %in% members$gene[duplicated(members$gene)]
    combined <- members[repeated, ] %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(
            species = dplyr::first(species),
            orthogroup = paste(sort(unique(orthogroup)), collapse = "/"),
            spis_ortholog = most_frequent(spis_ortholog[!is.na(spis_ortholog)]),
            Mmus_name = most_frequent(Mmus_name[!is.na(Mmus_name) & Mmus_name != ""]),
            Uniprot_gene = most_frequent(Uniprot_gene[!is.na(Uniprot_gene) & Uniprot_gene != ""]),
            is_tf = if (all(is.na(is_tf))) NA else any(is_tf, na.rm = TRUE),
            .groups = "drop"
        )
    members <- dplyr::bind_rows(members[!repeated, ], combined)
    for (spi in species_ids) {
        source_id <- if (spi == "Spin") "Spis" else spi
        genes <- read_tsv(file.path(reference_dir, "genes", sprintf("annot.%s.genes.csv", source_id)))
        in_og <- members[members$species == source_id, ] %>%
            dplyr::select(gene, spis_ortholog, orthogroup_name = Mmus_name,
                          Uniprot_name = Uniprot_gene, is_tf, orthogroup)
        other <- genes[!genes$gene %in% in_og$gene, c("gene", "spis_ortholog", "orthogroup_name", "Uniprot_name", "is_tf")]
        other$Uniprot_name <- ifelse(
            grepl("\\|", other$Uniprot_name),
            gsub(".*\\|(.+?)_.*", "\\1", other$Uniprot_name), NA_character_
        )
        # A slash-separated assignment has no single-OG annotation in the reference.
        result <- dplyr::bind_rows(in_og, other) %>%
            dplyr::left_join(annotations, by = "orthogroup")
        stopifnot(!anyDuplicated(result$gene))
        out_dir <- file.path(work_dir, sprintf("results_metacell_%s_filt", spi))
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        write.table(result, file.path(out_dir, sprintf("annot.%s.genes.csv", spi)),
                    sep = "\t", row.names = FALSE, quote = FALSE)
    }
}
