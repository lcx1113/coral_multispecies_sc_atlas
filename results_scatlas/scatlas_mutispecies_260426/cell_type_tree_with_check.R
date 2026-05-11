#libraries
setwd("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/scripts_260424/scatlas_mutispecies_260426/")
getwd()
suppressMessages(source("../coral_cs_functions/Downstream_functions.R"))
suppressMessages(source("../coral_cs_functions/Cross_species_functions.R"))
suppressMessages(source("../coral_cs_functions/helper.R"))
suppressMessages(require("umap"))
suppressMessages(require("phangorn"))
suppressMessages(require("ape"))
suppressMessages(require("adephylo"))
library(dplyr)
graphics.off()
message("=== Starting Cell Type Tree Construction ===\n")

# -------------------------
# NEW: fixed seed for reproducibility
#seed= 1237 ##seed is not set and retry should not used in final run. 
#set.seed(seed)
#message("[seed] set.seed(seed)")

# input
out_fn = "results_csps_cell_types"
dir.create(out_fn, showWarnings = FALSE)

# species
list_comparisons = list(
	Spin = list(c("Spin", "Adig", "Ahya", "Amur", "Fspp", "Gfas", "Hcoe", "Pcru", "Pdam", "Plut", "Xesp"), ref = "Spin")
)

# orthogroups
ogm <- read.table(
    "/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/final_data/broccoli/Metazoa_all_genome_260426/orthogroup_conservation_integrated.csv",
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE)
oga <- read.table(
    "/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/final_data/broccoli/Anthozoa_all_genome_260426/orthogroup_conservation.csv",
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE)
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)
ogm_gv_seu = dic_from_vecs(names = gsub("_", "-", ogm$gene), terms = ogm$orthogroup_name)

# fc threshold
fc_thr = 3.0
# completeness
ortholog_fraction = 0.7
fraction_genes = 0.7

# identify icc markers
min_ec = 0

# bootstrap
num_bs = 1000
collapse_bs = num_bs * 0.2

# NEW: bootstrap retry control
max_bootstrap_retries = 10
min_valid_bootstrap = 900

focus = "cell_type"
focid = "cts"

for (nn in 1:length(list_comparisons)) {

    set_id = names(list_comparisons)[nn]

	sps_list = list_comparisons[[nn]][[1]]
	sps_ref = list_comparisons[[nn]][["ref"]]
	sps_query = sps_list [ !sps_list %in% sps_ref ]

    unknown_gray_map = c(
        "unknown1" = "#4D4D4D",
        "unknown2" = "#6B6B6B",
        "unknown3" = "#8A8A8A",
        "unknown4" = "#A8A8A8",
        "unknown5" = "#C6C6C6"
    )

    ann_cts = c()
    for (spi in sps_list) {
		ctt_spi_fn = sprintf("samap_annotation/annot.%s.leiden.csv", spi, spi)
		ctt_spi = read.table(ctt_spi_fn, header = TRUE, comment.char = "", sep = "\t")
		ctt_spi$color[ctt_spi$cell_type %in% names(unknown_gray_map)] = unknown_gray_map[ctt_spi$cell_type[ctt_spi$cell_type %in% names(unknown_gray_map)]]
		ann_spi_v = dic_from_vecs(sprintf("%s|%s", spi, ctt_spi$cell_type), ctt_spi$color)
		ann_cts = c(ann_cts, ann_spi_v)
	}

	ref_genes_fn = sprintf("results_metacell_%s_filt/annot.%s.genes.csv", sps_ref, sps_ref)
	if (! file.exists(ref_genes_fn)) {
		warning(sprintf("File not found for %s: %s", sps_ref, ref_genes_fn))
		next
	}
	ref_gene_all = read.table(ref_genes_fn, sep = "\t", header = TRUE)
	ref_genes = ref_gene_all %>%
		mutate(sp1 = gene) %>%
		dplyr::select(sp1)

	for (spi in sps_query) {
		message(sprintf("csps %s | load %s-%s ICC data...", set_id, sps_ref, spi))
		icc_genes_fn = sprintf("results_metacell_%s_filt/csps/dat.icc.%s-%s.ec_scores.csv", sps_ref, sps_ref, spi)
		if (file.exists(icc_genes_fn)) {
			icc_genes_i = read.table(icc_genes_fn, header = TRUE)
			icc_genes_i = icc_genes_i [ icc_genes_i$ec_value >= min_ec, ]
		} else {
			message(sprintf("MISSING! -> %s", icc_genes_fn))
		}

		if (spi == sps_query[1]) {
			icc_genes_f = merge(ref_genes, icc_genes_i[,c(1,2)], by = "sp1", all.x = TRUE, all.y = FALSE)
		} else {
			icc_genes_f = merge(icc_genes_f, icc_genes_i[,c(1,2)], by = "sp1", all.x = TRUE, all.y = FALSE)
		}
	}

    colnames(icc_genes_f) = c(sps_ref, sps_query)
    icc_genes_f = icc_genes_f[rowSums(!is.na(icc_genes_f)) >= ortholog_fraction * ncol(icc_genes_f), ]

    fix_gene_name <- function(gene_names, prefix) {
        gene_names<-gsub("-","_", gene_names)
        ifelse(grepl(paste0("^", prefix), gene_names), gene_names, paste0(prefix, gene_names))
    }

	for (spi in sps_list) {
        if (spi == "Spin") {
            spi_w = "Spis"
		} else if (spi == "HcoeH" || spi == "HcoeI") {
			spi_w = "Hcoe"
		} else {
			spi_w = spi
		}

		message(sprintf("csps %s | load %s expression data...", set_id, spi))
		fps_genes_fn = sprintf("results_metacell_%s_filt/dat.%s.expression.cts_fp.rds", spi, spi)
		if (file.exists(fps_genes_fn)) {
			fps_i = readRDS(fps_genes_fn)
            rownames(fps_i)<-fix_gene_name(rownames(fps_i), prefix = paste0(spi_w, "_"))
			fps_i = fps_i[, !grepl("^unknown", colnames(fps_i)), drop = FALSE]
			fps_i_f = quantile_normalisation(fps_i)
		} else {
			message(sprintf("MISSING! -> %s", fps_genes_fn))
		}
		if (spi == sps_list[1]) {
			fps_i_f = fps_i_f [ icc_genes_f[,spi] , ]
			colnames(fps_i_f) = sprintf("%s|%s", spi, colnames(fps_i_f))
			csps_m = fps_i_f
		} else {
			fps_i_f_df = fps_i_f %>%
				as.data.frame() %>%
				tibble::rownames_to_column(var = spi)

			fps_i_f = icc_genes_f[, c(sps_ref, spi)] %>%
				left_join(fps_i_f_df, by = spi) %>%
				dplyr::select(-all_of(sps_ref), -all_of(spi))
			colnames(fps_i_f) = sprintf("%s|%s", spi, colnames(fps_i_f))
			csps_m = cbind(csps_m, fps_i_f)
		}
	}

	csps_m = csps_m [ rowSums(!is.na(csps_m)) >= ncol(csps_m) * fraction_genes , ]
    csps_m [ is.na(csps_m) ] = 0
    message(sprintf("csps %s | %s discretise fps at fc>=%.2f...", set_id, focid, fc_thr))
	csps_m_b = (csps_m >= fc_thr) * 1
    message(sprintf("csps %s | %s n=%i markers in %i %ss...", set_id, focid, nrow(csps_m_b), ncol(csps_m_b), focus))
    # stability checks
    if (nrow(csps_m_b) < 5 || ncol(csps_m_b) < 5) {
        message("[skip] csps_m_b too small, skip.")
        next
    }
    if (all(csps_m_b == 0) || all(csps_m_b == 1)) {
        message("[skip] csps_m_b all 0/1, skip.")
        next
    }
    if (any(rowSums(csps_m_b) == 0)) {
        csps_m_b = csps_m_b[rowSums(csps_m_b) > 0, ]
        message("[filter] removed zero rows.")
    }
    if (any(colSums(csps_m_b) == 0)) {
        csps_m_b = csps_m_b[, colSums(csps_m_b) > 0]
        message("[filter] removed zero columns.")
    }
    message(sprintf("csps %s | filtered %s n=%i markers in %i %ss...", set_id, focid, nrow(csps_m_b), ncol(csps_m_b), focus))
	ali_f = phangorn::as.phyDat(t(csps_m_b), type = "USER", levels = names(table(csps_m_b)), names = rownames(csps_m_b), return.index = TRUE)
	parsimony_fun = function(x) { reorder(ape::as.phylo(fastcluster::hclust(phangorn::dist.logDet(x), method = "average")), "postorder") }

	ali_f_phy = try(parsimony_fun(ali_f), silent = TRUE)
    if (inherits(ali_f_phy, "try-error")) {
        message("[skip] parsimony tree failed.")
        next
    }

	ali_tips = ali_f_phy$tip.label
	if (any(grepl("progenitors", ali_tips))) {
		ali_tips = ali_tips [ grepl("neurosecretory_progenitors", ali_tips) ]
		if (length(ali_tips) > 1) {
            ali_root = ape::getMRCA(ali_f_phy, ali_tips)
            ali_f_phy = ape::root(ali_f_phy, node = ali_root, resolve.root = TRUE)
        } else if (length(ali_tips) == 1) {
            ali_f_phy = ape::root(ali_f_phy, outgroup = ali_tips, resolve.root = TRUE)
	    }
    }

	ali_f_phy$edge.length [ is.infinite(ali_f_phy$edge.length) ] = max(ali_f_phy$edge.length [ !is.infinite(ali_f_phy$edge.length) ]) * 1.2

	# bootstrap with retry + filtering invalid trees
	message(sprintf("csps %s | %s parsimony tree, bootstrap...", set_id, focid))
	boo_list = NULL
	for (retry_i in 1:max_bootstrap_retries) {
	    message(sprintf("[bootstrap] attempt %d/%d", retry_i, max_bootstrap_retries))
	    boo_list = try(phangorn::bootstrap.phyDat(ali_f, parsimony_fun, multicore = TRUE,  bs = num_bs), silent = TRUE)

	    if (inherits(boo_list, "try-error")) {
	        message("[bootstrap] failed, retrying...")
	        set.seed(1237 + retry_i)
	        message(sprintf("[seed] reset seed to %d", 1237 + retry_i))
            next
	    }

	    valid_idx = vapply(boo_list, function(x) inherits(x, "phylo"), logical(1))
	    n_invalid = sum(!valid_idx)

	    if (n_invalid > 0) {
	        message(sprintf("[bootstrap] filtered %d invalid trees", n_invalid))
	        boo_list = boo_list[valid_idx]
	    }

	    if (length(boo_list) >= min_valid_bootstrap) {
	        message(sprintf("[bootstrap] %d valid trees retained, proceed.", length(boo_list)))
	        break
	    } else {
	        message(sprintf("[bootstrap] only %d valid trees (<%d), retrying...", length(boo_list), min_valid_bootstrap))
	        boo_list = NULL
            set.seed(1237 + retry_i)
	        message(sprintf("[seed] reset seed to %d", 1237 + retry_i))
	    }
	}

	if (is.null(boo_list) || length(boo_list) < min_valid_bootstrap) {
	    message("[skip] bootstrap failed after retries or too few valid trees.")
	    next
	}
    class(boo_list) = "multiPhylo"
    class(boo_list)

	collapse_nodes = which(prop.clades(ali_f_phy, boo_list) < collapse_bs) + length(ali_f_phy$tip.label)
	rownames(ali_f_phy$edge) = as.character(1:nrow(ali_f_phy$edge))
	collapse_edges = sort(as.numeric(rownames(ali_f_phy$edge) [ ali_f_phy$edge[,2] %in% collapse_nodes ]))
	ali_f_phy_col = ali_f_phy
	ali_f_phy_col$edge.length[ collapse_edges ] = 0
	ali_f_phy_col = ape::di2multi(ali_f_phy_col)

	ali_f_phy = ape::makeNodeLabel(ali_f_phy)

	plot_height = ceiling(ncol(csps_m_b) / 6 + 4)
	pdf(sprintf("%s/csps.%s.dendrogram.%s.UPGMA.pdf", out_fn, set_id, focid), height = plot_height, width = 12)
	phangorn::plotBS(ali_f_phy, boo_list, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP")
	ape::add.scale.bar()
	phangorn::plotBS(ali_f_phy_col, boo_list, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP", use.edge.length = FALSE)
	phangorn::plotBS(ali_f_phy, boo_list, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with TBE", fc_thr), root.edge = TRUE, method = "TBE", digits = 1)
	ape::add.scale.bar()
	ape::plot.phylo(ali_f_phy, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, underscore = TRUE, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f", fc_thr), root.edge = TRUE, use.edge.length = TRUE, show.node.label = TRUE)
	ape::add.scale.bar()
	dev.off()

	ape::write.tree(ali_f_phy, sprintf("%s/csps.%s.dendrogram.%s.UPGMA.newick", out_fn, set_id, focid))
	saveRDS(csps_m, sprintf("%s/csps.%s.cspsmatrix.%s.rds", out_fn, set_id, focid))

	# 1) 解决0或极短分支
	min_len <- 0.05
	ali_f_phy$edge.length[ali_f_phy$edge.length < min_len] <- min_len

	# 2) 拉长根第一层边
	root <- Ntip(ali_f_phy) + 1
	root_edges <- which(ali_f_phy$edge[,1] == root)
	extra <- 0.8 * max(ali_f_phy$edge.length)
	ali_f_phy$edge.length[root_edges] <- ali_f_phy$edge.length[root_edges] + extra


	pdf(sprintf("%s/csps.%s.dendrogram_radian.%s.UPGMAP.pdf", out_fn, set_id, focid), height = 20, width = 20)

	phangorn::plotBS(
	ali_f_phy, boo_list,
	tip.color = ann_cts[ali_f_phy$tip.label],
	font = 1,
	type = "fan",  # 改为圆形布局
	p = 0.7,
	digits = 1,
	bs.col ="black",
	cex = 0.85,
	label.offset = 0.05,
	main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr),
	method = "FBP"
	)
	ape::add.scale.bar()

	phangorn::plotBS(
	ali_f_phy, boo_list,
	tip.color = ann_cts[ali_f_phy$tip.label],
	font = 1,
	type = "fan",  # 改为圆形布局
	digits = 1,
	bs.col ="black",
	cex = 0.85,
	label.offset = 0.05,
	main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr),
	method = "FBP"
	)
	ape::add.scale.bar()

	bs <- ape::prop.clades(ali_f_phy, boo_list)
	bs <- bs / num_bs

	# 先画树，不显示 bootstrap（p 设得很高）
	phangorn::plotBS(
	ali_f_phy, boo_list,
	tip.color = ann_cts[ali_f_phy$tip.label],
	font = 1,
	type = "fan",
	p = 1.1,             # 不显示任何值
	digits = 2,
	bs.col = "black",
	cex = 0.85,
	label.offset = 0.05,
	main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr),
	method = "FBP"
	)

	# 只标注低于阈值的 bootstrap
	thr <- 0.7   # 例如只显示 <0.7
	ape::nodelabels(
	text = ifelse(bs < thr, round(bs, 2), NA),
	frame = "none",
	cex = 0.7,
	col = "black"
	)

	ape::add.scale.bar()
	dev.off()

}
