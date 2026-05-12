## Adapted from: https://github.com/xgrau/oculina-coral-sc-atlas
##################################################################
#' Compute cell type footprints for seurat clusters
#'
#' @param seu_object seurat object
#' @param cluster metadata field with cluster info
#' @param data either a seurat layer (RNA [default], SCT...) or a umi matrix (dgCMatrix or matrix class)
#' @param nbins num bins to split matrix into for large-scale calculations (10k)
#' @param keep_genes list of genes to keep for fp calculations (e.g. to exclude orphan peaks)
#'
#' @return per-cluster footprint matrix
#'
sca_seurat_cell_type_fp = function(seu_object, cluster = "louvain", data = "RNA", nbins = 10L, keep_genes = NULL) {
	
	sc_ct_label = seu_object@meta.data[,cluster]
	names(sc_ct_label) = rownames(seu_object@meta.data)
	
	# get UMI matrix
	if (is.character(data)) {
		umis = seu_object[[data]]$counts
	} else if ("dgCMatrix" %in% class(seu_object[[data]]$counts) | "matrix" %in% class(seu_object[[data]]$counts)) {
		umis = data
	} else {
		stop("`data` field should be either a layer in seurat object or a umi matrix")
	}
	
	# filter genes?
	if (!is.null(keep_genes)) {
		umis = umis [ keep_genes, ]
	}
	
	# geometric mean of gene expression per cell type
	message("Calculating geom mean")
	ct_geomean = tryCatch(
	  t(apply(umis, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1))),
	  error = function(e) {
	    warning(e)
	    message("Calculating geom mean for genes in ",nbins," bins")
	    umis_list <- vector("list",nbins)
	    sl <- split(1:nrow(umis),cut(1:nrow(umis),nbins))
	    for (i in seq_along(sl)) {
	      message(i," / ", nbins)
	      umisub <- umis[sl[[i]],]
	      ct_geomean_sub <- t(apply(umisub, 1,  function(x) tapply(x, sc_ct_label, function(y) exp(mean(log(1 + y))) - 1)))
	      umis_list[[i]] <- ct_geomean_sub 
	    }
	    do.call(rbind, umis_list)
	  }
	)
	
	# cell type mean size
	message("Calculating mean cluster size")
	ct_meansize = tapply(Matrix::colSums(umis), sc_ct_label, mean)
	
	# calculate footprint
	message("Calculating footprints per cluster")
	ideal_cell_size = pmin(1000,median(ct_meansize))
	g_fp = t(ideal_cell_size * t(ct_geomean) / as.vector(ct_meansize))
	fp_reg = 0.05
	g_fp_n = (fp_reg + g_fp) / apply(fp_reg + g_fp, 1, median)
	
	# for compatibility with other functions, return as a MC-like object
	return(g_fp_n)
	
}


#' Find markers along tree with Seurat and plot
#'
#' @param seu_object seurat object
#' @param seu_phylo cluster tree (`phylo` class) where tips are clusters defined in `seu_object`
#' @param layer seurat layer (RNA [default], SCT...) to use for significance testing
#' @param tip_color_dictionary named vector where names are clusters and contents are colors, to paint the tree
#' @param background_is background for DGE calculation in each node; it can be "sister" (sister node); "other" (all other tips); or a vector containing a list of tips
#' @param test test to use ("wilcox" default)
#' @param min_pct min pct cells in ingroup for DGE (`Seurat::FindMarers`)
#' @param assay seurat assay (RNA [default], SCT...) to use for significance testing
#' @param slot seurat slot ("counts" [default], "data"...) to use for significance testing
#' @param pval_thr threshold to define significantly DE genes (p-adjusted tests)
#' @param min_positive_markers_highlight min number of positively DE genes required to highlight a node in the resulting tree
#' @param output_file output file (painted tree)
#'
#' @return table with DGE markers in each node
#'
sca_seurat_tree_markers = function(
	seu_object,
	seu_phylo,
	tip_color_dictionary = NULL,
	background_is = "sister",
	assay = "RNA",
	slot = "counts",
	test = "wilcox",
	min_pct = 0.05,
	pval_thr = 1e-3,
	min_positive_markers_highlight = 100,
	output_file = NULL,
	width = NULL, height = NULL) {
	
	# generate nodes for tree of clusters
	seu_phylo = ape::makeNodeLabel(seu_phylo)
	
	# find markers in internal nodes
	descendants_from_nodes = adephylo::listTips(seu_phylo)
	list_labels = c(seu_phylo$tip.label, seu_phylo$node.label)
	nodes_to_test = list_labels[!list_labels %in% seu_phylo$node.label[1]]
	
	# init tables
	markers_per_bipartition = data.frame()
	significant_markers_per_node = data.frame(node = nodes_to_test, positive = NA, negative = NA, row.names = nodes_to_test)
	done_ancestors = c()
	
	for (focus_node in nodes_to_test) {
		
		# get background nodes
		if (background_is == "sister") {
			focus_node_ix = which(list_labels == focus_node)
			if(focus_node %in% seu_phylo$tip.label) {
				focus_tips = focus_node
			} else {
				focus_tips = as.character(sort(names(descendants_from_nodes[[focus_node]])))
			}
			ancestor_node = list_labels[phangorn::Ancestors(seu_phylo, focus_node_ix)][1]
			sister_tips   = as.character(sort(names(descendants_from_nodes[[ancestor_node]])))
			sister_tips   = sister_tips [ ! sister_tips %in% focus_tips ]
		} else if (background_is == "other") {
			focus_node_ix = which(list_labels == focus_node)
			if(focus_node %in% seu_phylo$tip.label) {
				focus_tips = focus_node
			} else {
				focus_tips = as.character(sort(names(descendants_from_nodes[[focus_node]])))
			}
			ancestor_node = list_labels[phangorn::Ancestors(seu_phylo, focus_node_ix)][1]
			sister_tips   = seu_phylo$tip.label [ ! seu_phylo$tip.label %in% focus_tips ]
		} else {
			sister_tips = seu_phylo$tip.label [ seu_phylo$tip.label %in% background_is ]
		}
		
		# find markers ingroup v sister
		# done_ancestors = c(ancestor_node, done_ancestors)
		markers_focus_node = Seurat::FindMarkers(seu_object, ident.1 = focus_tips, ident.2 = sister_tips, assay = assay, slot = slot, fc.slot = slot, min.pct = min_pct, test.use = test)
		markers_focus_node$focus_node = focus_node
		markers_focus_node$ident.1.tips  = paste(focus_tips, collapse = ",")
		markers_focus_node$ident.2.tips = paste(sister_tips, collapse = ",")
		# keep significant only
		markers_focus_node = markers_focus_node [ markers_focus_node$p_val_adj < pval_thr, ]
		# count significant
		num_sigs = table(factor(sign(markers_focus_node$avg_log2FC), levels = c(1,-1)))
		significant_markers_per_node[focus_node,"positive"] = num_sigs[1]
		significant_markers_per_node[focus_node,"negative"] = num_sigs[2]
		message(sprintf("tree markers | node %s in tree | n=%i markers (%i/%i pos/neg)", focus_node, sum(num_sigs), num_sigs[1], num_sigs[2]))
		# store for later
		markers_focus_node$gene = rownames(markers_focus_node)
		markers_per_bipartition = rbind(markers_per_bipartition, markers_focus_node)
		
	}
	
	# recompose gene names
	rownames(markers_per_bipartition) = NULL
	
	# plot tree of clusters
	if (is.null(width)) {
		width = length(seu_phylo$tip.label) / 50 + 8
	}
	if (is.null(height)) {
		height = length(seu_phylo$tip.label) / 10 + 10
	}
	plotting_function(output_file, width = width, height = height, res = NA, EXP = {

		# plot tree
		# tip names
		tip_labels = as.character(apply(significant_markers_per_node[rownames(significant_markers_per_node) %in% seu_phylo$tip.label,], 1, function(v) {sprintf("%s | +%i/-%i", v[1], as.numeric(v[2]), as.numeric(v[3])) } ))
		if (!is.null(tip_color_dictionary) & sum(names(tip_color_dictionary) %in% seu_phylo$tip.label) > 0) {
			tip_color_vector = tip_color_dictionary [ seu_phylo$tip.label ]
		} else {
			tip_color_vector = "purple3"
		}
		seu_phylo$tip.label = tip_labels
		# plot tree per se
		ape::plot.phylo(seu_phylo, font = 2, edge.color = "darkgray", show.node.label = FALSE, cex = 0.7, show.tip.label = TRUE, tip.color = tip_color_vector, underscore = TRUE)
		# node labels	
		node_labels = c("",as.character(apply(significant_markers_per_node[rownames(significant_markers_per_node) %in% seu_phylo$node.label,], 1, function(v) {sprintf("%s | +%i/-%i", v[1], as.numeric(v[2]), as.numeric(v[3])) } )))
		node_labels_color = c("blue4",scales::alpha("slateblue2",0.4)) [ factor(c(0,significant_markers_per_node[rownames(significant_markers_per_node) %in% seu_phylo$node.label,"positive"]) >= min_positive_markers_highlight, levels = c(TRUE, FALSE)) ]
		ape::nodelabels(text = node_labels, cex = 0.7, col = node_labels_color, frame = "none")
		# ape::tiplabels(text = tip_labels, cex = 0.7, col = tip_color_vector, frame = "none", font = 2, offset = 1, adj = c(0,0.5))
	
	})
	
	# return table
	return(markers_per_bipartition)
	
}




#' Find metacell clusters in Seurat object
#'
#' @param seu_object seurat object
#' @param layer seurat layer (RNA [default], SCT...) to use for significance testing
#' @param cell_subset subset of cells to retrieve from seurat object
#' @param mc_prefix prefix for metacells? (otherwise, they're reported as indexes, as usual)
#' @param run_name run ID (default is "tmp")
#' @param tmp_dir tmp dir to store mc objects ("tmp_scdb")
#' @param tmp_dir_rm whether to remove tmp dir upon completion (FALSE)
#' @param min_mc_size min metacell size (smaller ones are reported as non-clustered, default = 10)
#' @param gset_genes_blacklisted,gset_T_tot,gset_T_top3,cgraph_K... parameters in their respective mc functions
#'
#' @return list with various objects: vector of cell-metacell assignments, and all the mc objects here used
#'
sca_seurat_make_metacells = function(
	seu_object,
	layer = "RNA",
	cell_subset = NULL,
	mc_prefix = NULL,
	run_name = "tmp",
	tmp_dir = "tmp_scdb",
	tmp_dir_rm = FALSE,
	min_mc_size = 10,
	gset_genes_blacklisted = NULL,
	gset_T_tot = 100,
	gset_T_top3 = 2,
	gset_T_szcor = -0.05,
	gset_T_niche = 0.01,
	cgraph_K = 100,
	coclust_p_resamp = 0.75,
	coclust_n_resamp = 100,
	mc_K = 30,
	mc_alpha = 2) {
	
	# create input object
	message(sprintf("metacell | %s | set up scdb...", run_name))
	dir.create(tmp_dir, showWarnings = FALSE)
	metacell::scdb_init(tmp_dir, force_reinit = TRUE)
	
	message(sprintf("metacell | %s | matrix %s, create...", run_name, layer))
	if (is.null(cell_subset)) {
		umi = seu_object[[layer]]$counts
	} else {
		umi = seu_object[[layer]]$counts[,cell_subset]
	}
	mat = metacell::scm_new_matrix(umi, seu_object@meta.data[colnames(umi),], stat_type = "umi")
	metacell::scdb_add_mat(run_name, mat)
	message(sprintf("metacell | %s | matrix %s, get %i cells and %i genes...", run_name, layer, dim(umi)[2], dim(umi)[1]))
	
	# calculate gene stats
	message(sprintf("metacell | %s | gene stats, calculate...", run_name))
	suppressMessages(metacell::mcell_add_gene_stat(
		gstat_id = run_name,
		mat_id = run_name,
		force = TRUE))
	gstat = metacell::scdb_gstat(run_name)

	# define gene set
	message(sprintf("metacell | %s | gene set, calculate...", run_name))
	suppressMessages(metacell::mcell_gset_filter_multi(
		gstat_id = run_name,
		gset_id = run_name,
		T_tot = gset_T_tot,
		T_top3 = gset_T_top3,
		T_szcor = gset_T_szcor,
		T_niche = gset_T_niche,
		force_new = TRUE,
		blacklist = gset_genes_blacklisted))
	gset = metacell::scdb_gset(run_name)
	message(sprintf("metacell | %s | gene set, use n=%i marker genes...", run_name, length(gset@gene_set)))
		
	# perform reclustering (to get a cgraph object)
	message(sprintf("metacell | %s | cgraph, calculate...", run_name))
	suppressMessages(metacell::mcell_add_cgraph_from_mat_bknn(
		mat_id = run_name,
		gset_id = run_name,
		graph_id = run_name,
		K = cgraph_K,
		dsamp = FALSE))
	cgraph = metacell::scdb_cgraph(run_name)

	# coclustering graph via resampling (to get a coclust object)
	message(sprintf("metacell | %s | coclust, calculate...", run_name))
	suppressMessages(metacell::mcell_coclust_from_graph_resamp(
		coc_id = run_name,
		graph_id = run_name,
		min_mc_size = min_mc_size,
		p_resamp = coclust_p_resamp,
		n_resamp = coclust_n_resamp))
	coclust = metacell::scdb_coclust(run_name)

	# balanced coclustering (to get a mc object)
	message(sprintf("metacell | %s | metacells, define and calculate footprints...", run_name))
	
	# if it fails, assign all cells to a single metacell (it happens if you're creating metacells for small clusters)
	tryCatch(
	  
	  suppressMessages(metacell::mcell_mc_from_coclust_balanced(
		mc_id = run_name,
		coc_id = run_name,
		mat_id = run_name,
		K = mc_K,
		min_mc_size = min_mc_size,
		alpha = mc_alpha)),
	  
	  error = function(e) {
	    warning(e)
	    message(sprintf("WARNING: mc calculations failed with n=%i cells, assigning all of them to a single metacell...", dim(umi)[2]))
		mc_v = rep(1, ncol(umi))
		names(mc_v) = colnames(umi)
		metacell::mcell_new_mc(run_name, mc = mc_v, scmat = mat, outliers = mat@ignore_cells)
	  }
	  
	)
	mc = metacell::scdb_mc(run_name)
	message(sprintf("metacell | %s | metacells, %i cells classified into %i mcs...", run_name, length(mc@mc), ncol(mc@mc_fp)))
		
	
	# metacell assignments
	message(sprintf("metacell | %s | metacells, get assignments...", run_name))
	if (is.null(mc_prefix)) {
		cell_mc_v = mc@mc
	} else {
		cell_mc_v = sprintf("%s.%04d", mc_prefix, mc@mc)
		names(cell_mc_v) = names(mc@mc)
	}
	
	# housekeeping
	if (tmp_dir_rm) {
		message(sprintf("metacell | %s | clean up...", run_name))
		unlink(tmp_dir, recursive = TRUE)
	}
	
	# output
	out = list(
		cell_mc_v = cell_mc_v,
		gstat = gstat,
		gset = gset,
		cgraph = cgraph,
		coclust = coclust,
		mc = mc,
		mat = mat)
	return(out)
	
}

#' Find balanced clusters in arbitrary matrices, using a headless approach and user-provided variable features
#'
#' @param mat seurat object
#' @param mc_prefix prefix for metacells ("mc", set to NULL to report as indexes, as usual)
#' @param min_mc_size min metacell size (smaller ones are reported as non-clustered, default = 10)
#' @param top_K K-nn parameter for coclustering (coclustering analysed from top K-nn cells)
#' @param cgraph_K,cgraph_K_expand,cgraph_K_beta... parameters from the tgstat::tgs_knn and tgstat::tgs_graph functions
#' @param coclust_p_resamp,coclust_n_resamp... parameters from the tgstat::tgs_graph_cover_resample function
#'
#' @return list with various objects: vector of cell-metacell assignments, and all the mc objects here used
#'
sca_balanced_coclustering_headless = function(
	mat,
	variable_features = NULL,
	mc_prefix = "mc.",
	min_mc_size = 10,
	top_K = 30,
	alpha_relaxation = 2,
	cgraph_K = 100,
	cgraph_K_expand = 10,
	cgraph_K_beta = 3,
	coclust_cooling = 1.05,
	coclust_burnin = 10,
	coclust_p_resamp = 0.75,
	coclust_n_resamp = 100) {
	
	message(sprintf("balanced coclustering | matrix, create..."))
	if (is.null(variable_features)) { 
		variable_features = rownames(mat)
	}
	mat_f = as.matrix(mat[ variable_features, ])

	# analogous to `metacell::mcell_coclust_from_graph_resamp`, but using user-provided variable features rather than metacell ones
	message(sprintf("balanced coclustering | cell-cell correlation using n=%i variable features...", length(variable_features)))
	r_cor = tgstat::tgs_cor(mat_f, pairwise.complete.obs = FALSE, spearman = FALSE)
	message(sprintf("balanced coclustering | get closest K=%i nn cells per cell...", cgraph_K))
	r_knn = tgstat::tgs_knn(r_cor, knn = cgraph_K, diag = FALSE)
	r_cgraph = tgstat::tgs_graph(r_knn, knn = top_K, k_expand = cgraph_K_expand, k_beta = cgraph_K_beta)
	message(sprintf("balanced coclustering | cell coclustering graph at top K=%i nn...", top_K))
	r_coclust = tgstat::tgs_graph_cover_resample(
		r_cgraph,
		knn = top_K,
		min_cluster_size = min_mc_size,
		cooling = coclust_cooling,
		burn_in = coclust_burnin,
		p_resamp = coclust_p_resamp,
		n_resamp = coclust_n_resamp,
		method = "full")
	# get edges
	r_edges = r_coclust$co_cluster
	
	# analogous to `metacell::mcell_coclust_filt_by_k_deg`
	message(sprintf("balanced coclustering | filter coclustering graph by K=%i nn...", top_K))
	deg_wgt = as.matrix(table(c(r_edges$node1, r_edges$node2), c(r_edges$cnt, r_edges$cnt)))
	deg_cum = t(apply(deg_wgt, 1, function(x) cumsum(rev(x))))
	thresh_Kr = rowSums(deg_cum > top_K)
	thresh_K = rep(NA, nlevels(r_edges$node1))
	names(thresh_K) = levels(r_edges$node1)
	if (is.na(sum(as.numeric(names(thresh_K))))) {
		thresh_K[names(thresh_Kr)] = thresh_Kr
	} else {
		thresh_K[as.numeric(names(thresh_Kr))] = thresh_Kr
	}
	# get good edges
	r_edges_bool = thresh_K[r_edges$node1] < r_edges$cnt * alpha_relaxation | thresh_K[r_edges$node2] < r_edges$cnt * alpha_relaxation
	r_edges_f = r_edges[r_edges_bool,]
	message(sprintf("balanced coclustering | filtered %i edges, keep %i edges based on co-cluster imbalance...", nrow(r_edges) - sum(r_edges_bool) , sum(r_edges_bool) ))
	
	# rename object cols
	colnames(r_edges_f) = c("col1", "col2", "weight")
	r_edges_f$weight = r_edges_f$weight/max(r_edges_f$weight)

	message(sprintf("balanced coclustering | build reciprocal coclustering graph..."))
	r_edges_f = r_edges_f[r_edges_f$col1 != r_edges_f$col2,]
	r_edges_f_r = r_edges_f[,c("col1", "col2", "weight")]
	r_edges_f_r$col1 = r_edges_f$col2
	r_edges_f_r$col2 = r_edges_f$col1
	r_edges_f = rbind(r_edges_f, r_edges_f_r)
	
	# get clusters
	message(sprintf("balanced coclustering | get metacell assignments..."))
	r_node_clust = tgstat::tgs_graph_cover(r_edges_f, min_mc_size, cooling = 1.05, burn_in = 10)
	r_node_clust_no_outliers = r_node_clust [ r_node_clust$clust != 0, ]
	cell_mc_v = r_node_clust_no_outliers$cluster
	names(cell_mc_v) = r_node_clust_no_outliers$node
	message(sprintf("balanced coclustering | %i cells classified into %i mcs (%i outlier cells)...",  length(cell_mc_v), length(table(cell_mc_v)), sum(r_node_clust$clust == 0)))
	
	# metacell assignments
	if (!is.null(mc_prefix)) {
		message(sprintf("balanced coclustering | rename metacells..."))
		cell_mc_v = sprintf("%s%04d", mc_prefix, cell_mc_v)
		names(cell_mc_v) = r_node_clust_no_outliers$node
	}
	
	# finish
	return(cell_mc_v)

	
}

#' Find balanced clusters in arbitrary matrices, using a headless approach and user-provided variable features
#'
#' @param mat seurat object
#' @param mc_prefix prefix for metacells ("mc", set to NULL to report as indexes, as usual)
#' @param min_mc_size min metacell size (smaller ones are reported as non-clustered, default = 10)
#' @param top_K K-nn parameter for coclustering (coclustering analysed from top K-nn cells)
#' @param cgraph_K,cgraph_K_expand,cgraph_K_beta... parameters from the tgstat::tgs_knn and tgstat::tgs_graph functions
#' @param coclust_p_resamp,coclust_n_resamp... parameters from the tgstat::tgs_graph_cover_resample function
#'
#' @return list with various objects: vector of cell-metacell assignments, and all the mc objects here used
#'
sca_balanced_coclustering_headless_list = function(
	mat,
	variable_features = NULL,
	mc_prefix = "mc.",
	min_mc_size = 10,
	top_K = 30,
	alpha_relaxation = 2,
	cgraph_K = 100,
	cgraph_K_expand = 10,
	cgraph_K_beta = 3,
	coclust_cooling = 1.05,
	coclust_burnin = 10,
	coclust_p_resamp = 0.75,
	coclust_n_resamp = 100) {
	
	message(sprintf("balanced coclustering | matrix, create..."))
	if (is.null(variable_features)) { 
		variable_features = rownames(mat)
	}
	mat_f = as.matrix(mat[ variable_features, ])

	# analogous to `metacell::mcell_coclust_from_graph_resamp`, but using user-provided variable features rather than metacell ones
	message(sprintf("balanced coclustering | cell-cell correlation using n=%i variable features...", length(variable_features)))
	r_cor = tgstat::tgs_cor(mat_f, pairwise.complete.obs = FALSE, spearman = FALSE)
	message(sprintf("balanced coclustering | get closest K=%i nn cells per cell...", cgraph_K))
	r_knn = tgstat::tgs_knn(r_cor, knn = cgraph_K, diag = FALSE)
	r_cgraph = tgstat::tgs_graph(r_knn, knn = top_K, k_expand = cgraph_K_expand, k_beta = cgraph_K_beta)
	message(sprintf("balanced coclustering | cell coclustering graph at top K=%i nn...", top_K))
	r_coclust = tgstat::tgs_graph_cover_resample(
		r_cgraph,
		knn = top_K,
		min_cluster_size = min_mc_size,
		cooling = coclust_cooling,
		burn_in = coclust_burnin,
		p_resamp = coclust_p_resamp,
		n_resamp = coclust_n_resamp,
		method = "full")
	# get edges
	r_edges = r_coclust$co_cluster
	
	# analogous to `metacell::mcell_coclust_filt_by_k_deg`
	message(sprintf("balanced coclustering | filter coclustering graph by K=%i nn...", top_K))
	deg_wgt = as.matrix(table(c(r_edges$node1, r_edges$node2), c(r_edges$cnt, r_edges$cnt)))
	deg_cum = t(apply(deg_wgt, 1, function(x) cumsum(rev(x))))
	thresh_Kr = rowSums(deg_cum > top_K)
	thresh_K = rep(NA, nlevels(r_edges$node1))
	names(thresh_K) = levels(r_edges$node1)
	if (is.na(sum(as.numeric(names(thresh_K))))) {
		thresh_K[names(thresh_Kr)] = thresh_Kr
	} else {
		thresh_K[as.numeric(names(thresh_Kr))] = thresh_Kr
	}
	# get good edges
	r_edges_bool = thresh_K[r_edges$node1] < r_edges$cnt * alpha_relaxation | thresh_K[r_edges$node2] < r_edges$cnt * alpha_relaxation
	r_edges_f = r_edges[r_edges_bool,]
	message(sprintf("balanced coclustering | filtered %i edges, keep %i edges based on co-cluster imbalance...", nrow(r_edges) - sum(r_edges_bool) , sum(r_edges_bool) ))
	
	# rename object cols
	colnames(r_edges_f) = c("col1", "col2", "weight")
	r_edges_f$weight = r_edges_f$weight/max(r_edges_f$weight)

	message(sprintf("balanced coclustering | build reciprocal coclustering graph..."))
	r_edges_f = r_edges_f[r_edges_f$col1 != r_edges_f$col2,]
	r_edges_f_r = r_edges_f[,c("col1", "col2", "weight")]
	r_edges_f_r$col1 = r_edges_f$col2
	r_edges_f_r$col2 = r_edges_f$col1
	r_edges_f = rbind(r_edges_f, r_edges_f_r)
	
	# get clusters
	message(sprintf("balanced coclustering | get metacell assignments..."))
	r_node_clust = tgstat::tgs_graph_cover(r_edges_f, min_mc_size, cooling = 1.05, burn_in = 10)
	r_node_clust_no_outliers = r_node_clust [ r_node_clust$clust != 0, ]
	cell_mc_v = r_node_clust_no_outliers$cluster
	names(cell_mc_v) = r_node_clust_no_outliers$node
	message(sprintf("balanced coclustering | %i cells classified into %i mcs (%i outlier cells)...",  length(cell_mc_v), length(table(cell_mc_v)), sum(r_node_clust$clust == 0)))
	
	# metacell assignments
	if (!is.null(mc_prefix)) {
		message(sprintf("balanced coclustering | rename metacells..."))
		cell_mc_v = sprintf("%s.%04d", mc_prefix, cell_mc_v)
		names(cell_mc_v) = r_node_clust_no_outliers$node
	}
	
	# finish
	return(list(assignments = cell_mc_v, coc = r_coclust, coc_edges = r_edges, coc_edges_filtered = r_edges_f, graph = r_node_clust_no_outliers))

	
}


# find PCA elbow function from findPC: Automatic selection of number of principal components
# Citation: https://doi.org/10.1093/bioinformatics/btac235
#'
#' @param sdev vector of standard deviations for each PC (ordered)
#'
#' @return list with the elbow points selected using various criteria (I [Xavi] like "Perpendicular line")
#'
find_pca_elbow = function(sdev){

	x<-1:length(sdev)
	eb<-data.frame(x,sdev)

	# piecewise linear model
	sse<-NULL
	for (k in 1:length(sdev)) {
	D=ifelse(x<=k,0,1)
	x2<-(x-k)*D
	sse[k]=sum(lm(sdev~x+x2)$residuals^2)
	}
	dim_plm<-which.min(sse)
	names(dim_plm)<-'Piecewise linear model'

	# first derivative
	df1<-diff(sdev,differences=1)
	den<-density(df1)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_fid<-max(which(df1<mean(df1)))+1
	} else {
	cutoff<-sum(rlev[-((length(rlev)-1):length(rlev))])
	dim_fid<-max(which(df1<den$x[cutoff]))+1 }
	names(dim_fid)<-'First derivative'

	# second derivative
	df2<-diff(sdev,differences=2)
	df2p<-df2[df2>0]
	den<-density(df2p)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_sed<-which.max(df2)+1
	} else {
	cutoff<-sum(rlev[1:2])
	dim_sed<-max(which(df2>den$x[cutoff]))+1 }
	names(dim_sed)<-'Second derivative'

	# preceding residual
	fit<-NULL;res<-NULL
	for (i in 1:(length(sdev)-2)) {
	fit[[i]]<-lm(sdev~x,data=eb[(i+1):length(sdev),])
	res[i]<-sdev[i]-predict(fit[[i]],newdata=data.frame(x=i))
	}
	den<-density(res)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_pr<-max(which(res>mean(res)))+1
	} else {
	cutoff<-sum(rlev[1:2])
	dim_pr<-max(which(res>den$x[cutoff]))+1 }
	names(dim_pr)<-'Preceding residual'

	# perpendicular line
	A<-c(1,sdev[1]);B<-c(length(sdev),sdev[length(sdev)]);Dist<-NULL
	for (i in 2:(length(sdev)-1)) {
	C<-c(i,sdev[i]);D<-cbind(rbind(A,B,C),rep(1,3))
	S<-1/2*abs(det(D));Dist[i]<-2*S/dist(rbind(A,B))
	}
	dim_perl<-which.max(Dist)
	names(dim_perl)<-'Perpendicular line'

	# k-means clustering
	set.seed(2022)
	dim_clu<-min(kmeans(sdev,2)$size)+1
	names(dim_clu)<-'K-means clustering'

	dim_all<-c(dim_plm,dim_fid,dim_sed,dim_pr,dim_perl,dim_clu)
	return(dim_all)
}



#' Prepare object to run differential expression on SCT assay with multiple models
#'
#' Given a merged object with multiple SCT models, this function uses minimum
#' of the median UMI (calculated using the raw UMI counts) of individual objects
#' to reverse the individual SCT regression model using minimum of median UMI
#' as the sequencing depth covariate.
#' The counts slot of the SCT assay is replaced with recorrected counts and
#' the data slot is replaced with log1p of recorrected counts.
#' @param object Seurat object with SCT assays
#' @param assay Assay name where for SCT objects are stored; Default is 'SCT'
#' @param verbose Print messages and progress
#' @importFrom Matrix Matrix
#' @importFrom SeuratObject SparseEmptyMatrix
#' @importFrom pbapply pblapply
#' @importFrom future.apply future_lapply
#' @importFrom future nbrOfWorkers
#' @importFrom sctransform correct_counts
#' @importFrom SeuratObject JoinLayers
#'
#' @return Returns a Seurat object with recorrected counts and data in the SCT assay.
#' @export
#'
#' @concept differential_expression
#' @template section-progressr
#' @template section-future
#' @examples
#' data("pbmc_small")
#' pbmc_small1 <- SCTransform(object = pbmc_small, variable.features.n = 20, vst.flavor="v1")
#' pbmc_small2 <- SCTransform(object = pbmc_small, variable.features.n = 20, vst.flavor="v1")
#' pbmc_merged <- merge(x = pbmc_small1, y = pbmc_small2)
#' pbmc_merged <- PrepSCTFindMarkers(object = pbmc_merged)
#' markers <- FindMarkers(
#'   object = pbmc_merged,
#'   ident.1 = "0",
#'   ident.2 = "1",
#'   assay = "SCT"
#' )
#' pbmc_subset <- subset(pbmc_merged, idents = c("0", "1"))
#' markers_subset <- FindMarkers(
#'   object = pbmc_subset,
#'   ident.1 = "0",
#'   ident.2 = "1",
#'   assay = "SCT",
#'   recorrect_umi = FALSE
#' )
#'
sca_PrepSCTFindMarkers_force <- function(object, assay = "SCT", verbose = TRUE) {
  if (verbose && future::nbrOfWorkers() == 1) {
    my.lapply <- pbapply::pblapply
  } else {
    my.lapply <- future.apply::future_lapply
  }
  if (length(x = levels(x = object[[assay]])) == 1) {
    if (verbose) {
      message("Only one SCT model is stored - skipping recalculating corrected counts")
    }
    return(object)
  }
  observed_median_umis <- lapply(
    X = SCTResults(object = object[[assay]], slot = "cell.attributes"),
    FUN = function(x) median(x[, "umi"])
  )
  model.list <- slot(object = object[[assay]], name = "SCTModel.list")
  median_umi.status <- lapply(X = model.list,
                              FUN = function(x) { return(tryCatch(
                                expr = slot(object = x, name = 'median_umi'),
                                error = function(...) {return(NULL)})
                              )})
  if (any(is.null(x = unlist(x = median_umi.status)))){
    # For old SCT objects  median_umi is set to median umi as calculated from obserbed UMIs
    slot(object = object[[assay]], name = "SCTModel.list") <- lapply(X = model.list,
                                                                     FUN = UpdateSlots)
    SCTResults(object = object[[assay]], slot = "median_umi") <- observed_median_umis

  }
  model_median_umis <- SCTResults(object = object[[assay]], slot = "median_umi")
  min_median_umi <- min(unlist(x = observed_median_umis), na.rm = TRUE)
  if (verbose) {
    message(paste0("Found ",
                   length(x = levels(x = object[[assay]])),
                   " SCT models.",
                   " Recorrecting SCT counts using minimum median counts: ",
                   min_median_umi))
  }
  umi.assay <- unique(
    x = unlist(
      x = SCTResults(object = object[[assay]], slot = "umi.assay")
    )
  )
  if (length(x = umi.assay) > 1) {
    stop("Multiple UMI assays are used for SCTransform: ",
         paste(umi.assay, collapse = ", ")
    )
  }
  umi.layers <- SeuratObject::Layers(object = object, assay = umi.assay, search = 'counts')
  if (length(x = umi.layers) > 1) {
    object[[umi.assay]] <- SeuratObject::JoinLayers(
      object = object[[umi.assay]],
      layers = "counts", new = "counts")
  }
  raw_umi <- SeuratObject::GetAssayData(object = object, assay = umi.assay, slot = "counts")
  corrected_counts <- Matrix::Matrix(
    nrow = nrow(x = raw_umi),
    ncol = ncol(x = raw_umi),
    data = 0,
    dimnames = dimnames(x = raw_umi),
    sparse = TRUE
  )
  cell_attr <- SCTResults(object = object[[assay]], slot = "cell.attributes")
  model_pars_fit <- lapply(
    X = SCTResults(object = object[[assay]], slot = "feature.attributes"),
    FUN = function(x) x[, c("theta", "(Intercept)", "log_umi")]
  )
  arguments <- SCTResults(object = object[[assay]], slot = "arguments")
  model_str <- SCTResults(object = object[[assay]], slot = "model")
  set_median_umi <- rep(min_median_umi, length(levels(x = object[[assay]])))
  names(set_median_umi) <- levels(x = object[[assay]])
  set_median_umi <- as.list(set_median_umi)
  all_genes <- rownames(x = object[[assay]])
  # correct counts
  my.correct_counts <- function(model_name){
    model_genes <- rownames(x = model_pars_fit[[model_name]])
      x <- list(
        model_str = model_str[[model_name]],
        arguments = arguments[[model_name]],
        model_pars_fit = as.matrix(x = model_pars_fit[[model_name]]),
        cell_attr = cell_attr[[model_name]]
      )
      cells <- rownames(x = cell_attr[[model_name]])
      umi <- raw_umi[all_genes, cells]

      umi_corrected <- sctransform::correct_counts(
        x = x,
        umi = umi,
        verbosity = 0,
        scale_factor = min_median_umi
      )
      missing_features <- setdiff(x = all_genes, y = rownames(x = umi_corrected))
      corrected_counts.list <- NULL
      gc(verbose = FALSE)
      empty <- SeuratObject::SparseEmptyMatrix(nrow = length(x = missing_features), ncol = ncol(x = umi_corrected))
      rownames(x = empty) <- missing_features
      colnames(x = umi_corrected) <- colnames(x = umi_corrected)

      umi_corrected <- rbind(umi_corrected, empty)[all_genes,]

      return(umi_corrected)
  }
  corrected_counts.list <- my.lapply(X = levels(x = object[[assay]]),
                                     FUN = my.correct_counts)
  names(x = corrected_counts.list) <- levels(x = object[[assay]])

  corrected_counts <- do.call(what = sca_MergeSparseMatrices, args = corrected_counts.list)
  corrected_counts <- as.sparse(x = corrected_counts)
  corrected_data <- log1p(x = corrected_counts)
  suppressWarnings({object <- SeuratObject::SetAssayData(object = object,
                                           assay = assay,
                                           slot = "counts",
                                           new.data = corrected_counts)})
  suppressWarnings({object <- SeuratObject::SetAssayData(object = object,
                                           assay = assay,
                                           slot = "data",
                                           new.data = corrected_data)})
  Seurat::SCTResults(object = object[[assay]], slot = "median_umi") <- set_median_umi
  return(object)
}



# Merge a list of sparse matrixes
#' @importFrom Matrix summary sparseMatrix
sca_MergeSparseMatrices <- function(...) {

  colname.new <- character()
  rowname.new <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (mat in list(...)) {
    colname.old <- colnames(x = mat)
    rowname.old <- rownames(x = mat)

    # does not check if there are overlapping cells
    colname.new <- union(x = colname.new, y = colname.old)
    rowname.new <- union(x = rowname.new, y = rowname.old)

    colindex.new <- match(x = colname.old, table = colname.new)
    rowindex.new <- match(x = rowname.old, table = rowname.new)

    ind <- summary(object = mat)
    # Expand the list of indices and x
    i <- c(i, rowindex.new[ind[,1]])
    j <- c(j, colindex.new[ind[,2]])
    x <- c(x, ind[,3])
  }

  merged.mat <- Matrix::sparseMatrix(
	i=i,
	j=j,
	x=x,
	dims=c(length(rowname.new), length(colname.new)),
	dimnames=list(rowname.new, colname.new))
  return (merged.mat)
}
