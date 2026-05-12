## Adapted from: https://github.com/xgrau/oculina-coral-sc-atlas
## With other modification to adapt to our dataset
####################################################################
#libraries
setwd("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/scripts_260424/scatlas_mutispecies_260426/")
getwd()
suppressMessages(source("../coral_cs_functions/Seurat_functions.R"))
suppressMessages(source("../coral_cs_functions/gene-set-analysis.R"))
suppressMessages(source("../coral_cs_functions/helper.R"))
suppressMessages(source("../coral_cs_functions/Downstream_functions.R"))
suppressMessages(require("Seurat"))
suppressMessages(require("SeuratWrappers"))
graphics.off()
fix_gene_name <- function(gene_names, prefix) {
        gene_names<-gsub("-","_", gene_names)
        ifelse(grepl(paste0("^", prefix), gene_names), gene_names, paste0(prefix, gene_names))
    }

# data index
sps_list = c("Adig", "Ahya", "Amur", "Fspp", "Gfas", "Hcoe", "Pcru", "Pdam", "Plut", "Spin", "Xesp")
df <- read.table("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/251219_mc_lei/id2symbol.txt", 
                header = TRUE, 
                sep = "\t", 
                stringsAsFactors = FALSE,
                comment.char = "",
                quote = "")
# 去除所有的引号字符
df$X.ID. <- gsub('"', '', df$X.ID.)
df$X.symbol. <- gsub('"', '', df$X.symbol.)
symbol_to_id<- dic_from_vecs(gsub("_","-", df$X.symbol.), gsub("_","-", df$X.ID.))

for (spi in sps_list) {
    # set working species
    if (spi == "Spin") {
        spi_w = "Spis"
	} else if (spi == "HcoeH" || spi == "HcoeI") {
        spi_w = "Hcoe"
    } else {
		spi_w = spi
	}

    # output
	out_fn = sprintf("results_metacell_%s_filt", spi)
	dir.create(out_fn, showWarnings = FALSE)

    gene_annot = read.table((sprintf("%s/annot.%s.genes.csv", out_fn, spi)),
                            sep = "\t",
                            header = TRUE,
                            comment.char = "")

    #read RDS Seurat object
    message(sprintf("metacell | %s | load Seurat...", spi))
	seu = readRDS(sprintf("%s/%s_lei_mc.rds", out_fn, spi))

    #read cell annotation
    message(sprintf("metacell | %s | load annotated clusters...", spi))
	ctt = read.table(sprintf("./samap_annotation/annot.%s.leiden.csv", spi), 
                 sep = "\t", 
                 header = TRUE, 
                 comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt_leiden    = ctt$cluster[ctt$cluster_class == "leiden"]
	ctt_lou_cts_v = dic_from_vecs(ctt$cluster, ctt$cell_type)
	ctt_lou_col_v = dic_from_vecs(ctt$cluster, ctt$color)
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)

    # update metadata
	message(sprintf("metacell | %s | update metadata...", spi))
	seu@meta.data$cell_type = ctt_lou_cts_v[as.character(seu$leiden)]
	seu@meta.data$cell_type[as.character(seu$metacell) %in% ctt$cluster] = ctt_lou_cts_v[as.character(seu$metacell)][as.character(seu$metacell) %in% ctt$cluster]
	seu@meta.data$color = ctt_cts_col_v[as.character(seu$cell_type)]
    
    #plot pie plot for cell type
    metadata <- seu@meta.data

    # 处理数据：合并所有unknown*的细胞类型
    cell_type_summary <- metadata %>%
    mutate(
        # 将所有以"unknown"开头的细胞类型合并
        cell_type_modified = ifelse(grepl("^unknown", cell_type), 
                                "unknown", 
                                as.character(cell_type)),
        # 为合并的unknown类型设置灰色
        color_modified = ifelse(grepl("^unknown", cell_type), 
                            "grey70", 
                            color)
    ) %>%
    group_by(cell_type_modified, color_modified) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(percentage = round(count / sum(count) * 100, 1)) %>%
    arrange(desc(count))

    # 创建用于图例的标签（包含百分比）
    cell_type_summary <- cell_type_summary %>%
    mutate(legend_label = paste0(cell_type_modified, " (", percentage, "%)"))

    # 创建饼状图
    pie_chart <- ggplot(cell_type_summary, aes(x = "", y = count, fill = legend_label)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = setNames(cell_type_summary$color_modified, 
                                        cell_type_summary$legend_label)) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.key.size = unit(0.8, "cm")
    ) +
    labs(
        title = "Cell Type Distribution",
        fill = "Cell Type"
    ) +
    guides(fill = guide_legend(ncol = 1))

    # 保存图形
    ggsave(sprintf("%s/seu.%s.pie.ct.pdf", out_fn, spi),
        plot = pie_chart, 
        width = 12, height = 10, 
        dpi = 300, 
        bg = "white")


	# cross-tabulation of metacells and leiden clusters
	message(sprintf("metacell | %s | cross-tabulate metacells and leiden...", spi))
	xtb = table(seu$metacell, seu$leiden)
	xtf_lei = xtb / rowSums(xtb)
	order_cols = names(sort(apply(xtf_lei, 1, function(x) which.max(x))))
	xtf_lei = xtf_lei[order_cols, ]
	top_cluster_per_mc = levels(seu$leiden)[apply(xtf_lei, 1, which.max)]
	pp = plot_complex_heatmap(
		xtf_lei,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = FALSE,
		cluster_col = FALSE,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		categories_row = top_cluster_per_mc,
		colors_row = ctt_lou_col_v,
		categories_col = colnames(xtf_lei),
		colors_col = ctt_lou_col_v
	)
	pdf(sprintf("%s/seu.%s.crosstab.lei-metacell.pdf", out_fn, spi), height = (dim(xtf_lei)[1] / 10) + 5, width = (dim(xtf_lei)[2] / 10) + 5)
	print(pp)
	dev.off()

    # cross-tabulation of metacells and cell types
	message(sprintf("metacell | %s | cross-tabulate metacells and cell types...", spi))
	xtb = table(seu$metacell, seu$cell_type)
	xtf_cts = xtb / rowSums(xtb)
	order_cols = names(sort(apply(xtf_cts, 1, function(x) which.max(x))))
	xtf_cts = xtf_cts[order_cols, ]
	top_cluster_per_mc = levels(seu$cell_type)[apply(xtf_cts, 1, which.max)]
	pp = plot_complex_heatmap(
		xtf_cts,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = FALSE,
		cluster_col = FALSE,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		categories_row = top_cluster_per_mc,
		colors_row = ctt_cts_col_v,
		categories_col = colnames(xtf_cts),
		colors_col = ctt_cts_col_v
	)
	pdf(sprintf("%s/seu.%s.crosstab.cts-metacell.pdf", out_fn, spi), height = (dim(xtf_cts)[1] / 10) + 5, width = (dim(xtf_cts)[2] / 10) + 5)
	print(pp)
	dev.off()

    # setup metacell annotation table, leiden
	message(sprintf("metacell | %s | metacell annotation table, leiden...", spi))
	mc_annot_lei_l = apply(xtf_lei, 1, function(v) { 
		vs = sort(v, decreasing = TRUE)
		vs = vs[vs > 0.05]
		vn = names(vs)
		vt = vn[1]
		vv = paste(sprintf("%s|%.2f", vn, vs), collapse = ";")
		return(list(top = vt, others = vv))
	})
	mc_annot_lei_top = unlist(lapply(mc_annot_lei_l, function(i) i$top))
	mc_annot_lei_ann = unlist(lapply(mc_annot_lei_l, function(i) i$others))

	# setup metacell annotation table, cell types
	message(sprintf("metacell | %s | metacell annotation table, cell types...", spi))
	mc_annot_cts_l = apply(xtf_cts, 1, function(v) { 
		vs = sort(v, decreasing = TRUE)
		vs = vs[vs > 0.05]
		vn = names(vs)
		vt = vn[1]
		vv = paste(sprintf("%s|%.2f", vn, vs), collapse = ";")
		return(list(top = vt, others = vv))
	})
	mc_annot_cts_top = unlist(lapply(mc_annot_cts_l, function(i) i$top))
	mc_annot_cts_ann = unlist(lapply(mc_annot_cts_l, function(i) i$others))
	
	# create metacell annotation table
	ctm = data.frame(
		metacell = names(mc_annot_cts_top),
		metacell_annotation = mc_annot_cts_top,
		color = ctt_cts_col_v[mc_annot_cts_top],
		metacell_top_cell_type = mc_annot_cts_top[names(mc_annot_cts_top)],
		metacell_ann_cell_type = mc_annot_cts_ann[names(mc_annot_cts_top)],
		metacell_top_leiden    = mc_annot_lei_top[names(mc_annot_cts_top)],
		metacell_ann_leiden = mc_annot_lei_ann[names(mc_annot_cts_top)]
	)
	# order metacells: by most abundant cell type, first
	ctm$metacell_top_cell_type = factor(ctm$metacell_top_cell_type, levels = levels(ctt$cell_type))
	ctm$metacell_top_leiden = factor(ctm$metacell_top_leiden, levels = ctt_leiden)
	# order metacells: by metacell clustering, second
	SeuratObject::Idents(seu) = seu@meta.data[, "metacell"]
	seu_num_pcs_v = find_pca_elbow(seu@reductions$pca@stdev)
	seu_num_pcs = seu_num_pcs_v["First derivative"]
	seu = Seurat::BuildClusterTree(seu, assay = "RNA", reduction = "pca", dims = 1:seu_num_pcs)
	seu_phylo = Seurat::Tool(seu, slot = "Seurat::BuildClusterTree")
	ctm_metacell_factor = factor(ctm$metacell, levels = seu_phylo$tip.label)
	# apply dual order
	ctm = ctm[order(ctm$metacell_top_cell_type, ctm$metacell_top_leiden, ctm_metacell_factor), ]
	# save
	write.table(ctm, sprintf("%s/annot.%s.mcs.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)
    
    # reorder metacells
    seu@meta.data$metacell = factor(seu@meta.data$metacell, levels = ctm$metacell)

    ## Dimensionality reduction ##	
	message(sprintf("metacell | %s | recolor umap...", spi))
    # 重新进行UMAP降维
    dims_used <- seu@commands$FindNeighbors.SCT.pca@params$dims
    seu <- RunUMAP(
        seu,
        reduction = "pca",              # 使用相同的PCA
        dims = dims_used,               # 使用相同的PC数量
        n.neighbors = 50,                   # 与cgraph_K/3相近（见下方说明）
        min.dist = 0.3,                     # 控制点的紧密程度
        metric = "euclidean",               # 距离度量
        spread = 1.0,                       # 控制嵌入的扩散程度
        verbose = FALSE,
        reduction.name = "umap",
        reduction.key = "UMAP_"
    )
    # 检查可用的降维结果
    available_reductions = names(seu@reductions)
    message(sprintf("Available reductions: %s", paste(available_reductions, collapse = ", ")))
    # 确定使用的 UMAP 名称
    if ("umap" %in% available_reductions) {
        umap_name = "umap"
    } else if ("umap_integrated" %in% available_reductions) {
        umap_name = "umap_integrated"
    } else {
        stop("No UMAP reduction found in Seurat object!")
    }
    message(sprintf("Using UMAP: %s", umap_name))

    # 获取 UMAP 坐标
    umap_coords = seu@reductions[[umap_name]]@cell.embeddings
    # 确保有两列（UMAP1 和 UMAP2）
    if (ncol(umap_coords) < 2) {
        stop("UMAP must have at least 2 dimensions!")
    }
    # 计算质心（centroids）
    # Cell type centroids
    ct_centroids_d = data.frame(
        umap1 = umap_coords[, 1], 
        umap2 = umap_coords[, 2], 
        label = seu@meta.data$cell_type,
        stringsAsFactors = FALSE
    )
    ct_centroids_a = aggregate(cbind(umap1, umap2) ~ label, ct_centroids_d, median)
    # Leiden centroids
    le_centroids_d = data.frame(
        umap1 = umap_coords[, 1], 
        umap2 = umap_coords[, 2], 
        label = seu@meta.data$leiden,
        stringsAsFactors = FALSE
    )
    le_centroids_a = aggregate(cbind(umap1, umap2) ~ label, le_centroids_d, median)
    # Metacell centroids
    mc_centroids_d = data.frame(
        umap1 = umap_coords[, 1], 
        umap2 = umap_coords[, 2], 
        label = seu@meta.data$metacell,
        stringsAsFactors = FALSE
    )
    mc_centroids_a = aggregate(cbind(umap1, umap2) ~ label, mc_centroids_d, median)

    # 绘制 UMAP 图
    pdf(sprintf("%s/seu.%s.umap.pdf", out_fn, spi), height = 12, width = 12)
    layout(matrix(1:4, nrow = 2, byrow = TRUE))
    # Plot 1: 整体背景 + metacell 质心
    plot(
        umap_coords[, 1], 
        umap_coords[, 2], 
        col = "seashell3", 
        pch = 19, 
        cex = 0.5, 
        main = sprintf("%s - UMAP by cell type", spi), 
        xlab = "UMAP1", 
        ylab = "UMAP2"
    )
    points(
        mc_centroids_a$umap1, 
        mc_centroids_a$umap2, 
        col = colorspace::darken(ctm_mcs_col_v[as.character(mc_centroids_a$label)], 0.4), 
        bg = ctm_mcs_col_v[as.character(mc_centroids_a$label)], 
        pch = 21, 
        cex = 2
    )
    # Plot 2: 图例
    plot(0, 0, xaxt = "n", xlab = "", ylab = "", yaxt = "n", col = NA, axes = FALSE)
    legend(
        "topleft", 
        names(ctt_cts_col_v), 
        col = ctt_cts_col_v, 
        pch = 19, 
        cex = 0.8, 
        title = "Cell Type", 
        bty = "n", 
        ncol = 2
    )
    # Plot 3: 按 leiden cluster 着色
    layout(1)
    plot(
        umap_coords[, 1], 
        umap_coords[, 2],
        col = seu@meta.data$color, 
        pch = 19, 
        cex = 0.5, 
        main = sprintf("%s - UMAP by leiden cluster", spi), 
        xlab = "UMAP1", 
        ylab = "UMAP2"
    )
    # 添加 leiden 标签
    text(
        le_centroids_a$umap1, 
        le_centroids_a$umap2, 
        le_centroids_a$label, 
        col = scales::alpha("gray5", 0.9), 
        cex = 0.8,
        font = 2
    )
    title(sub = sprintf("n=%i cells into %i metacells", ncol(seu), nrow(ctm)))
    # Plot 4: 按 cell type 着色 + metacell 和 cell type 标签
    layout(1)
    plot(
        umap_coords[, 1], 
        umap_coords[, 2],
        col = seu@meta.data$color, 
        pch = 19, 
        cex = 0.5, 
        main = sprintf("%s - UMAP by cell type", spi), 
        xlab = "UMAP1", 
        ylab = "UMAP2"
    )
    # 添加 metacell 标签（小号灰色）
    text(
        mc_centroids_a$umap1, 
        mc_centroids_a$umap2, 
        mc_centroids_a$label, 
        col = colorspace::darken(ctm_mcs_col_v[as.character(mc_centroids_a$label)], 0.6), 
        cex = 0.4,
        font = 1
    )
    # 添加 cell type 标签（大号黑色）
    text(
        ct_centroids_a$umap1, 
        ct_centroids_a$umap2, 
        ct_centroids_a$label, 
        col = scales::alpha("gray5", 0.9), 
        cex = 0.8,
        font = 2
    )
    title(sub = sprintf("n=%i cells into %i metacells", ncol(seu), nrow(ctm)))
    dev.off()
    message(sprintf("UMAP plots saved to: %s/seu.%s.umap.pdf", out_fn, spi))

    ## Expression tables ##
	message(sprintf("metacell | %s | orphan kept...", spi))
	keep_genes = rownames(seu@assays$RNA$counts)
	#keep_genes = keep_genes[!grepl("orphan", keep_genes)]#here removed all orphan genes from further analysis
	
	message(sprintf("metacell | %s | footprints, metacell...", spi))
	mc_fp = sca_seurat_cell_type_fp(seu_object = seu, cluster = "metacell", data = "RNA", keep_genes = keep_genes)
	mc_counts  = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes, ]), 1, function(v) { tapply(v, seu$metacell, sum) }))
	mc_fracell = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes, ]), 1, function(v) { tapply(v, seu$metacell, function(vv) { 100 * sum(vv > 0) / length(vv) }) }))
	mc_umifrac = sca_mc_gene_umifrac_noobj(mc_counts, multiplying_factor = 10000)

	message(sprintf("metacell | %s | footprints, cell type...", spi))
	ct_fp = sca_seurat_cell_type_fp(seu_object = seu, cluster = "cell_type", data = "RNA", keep_genes = keep_genes)
	ct_counts  = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes, ]), 1, function(v) { tapply(v, seu$cell_type, sum) }))
	ct_fracell = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes, ]), 1, function(v) { tapply(v, seu$cell_type, function(vv) { 100 * sum(vv > 0) / length(vv) }) }))
	ct_umifrac = sca_mc_gene_umifrac_noobj(ct_counts, multiplying_factor = 10000)
	
	message(sprintf("metacell | %s | footprints, leiden...", spi))
	le_fp = sca_seurat_cell_type_fp(seu_object = seu, cluster = "leiden", data = "RNA", keep_genes = keep_genes)
	le_counts  = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes, ]), 1, function(v) { tapply(v, seu$leiden, sum) }))
	le_fracell = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes, ]), 1, function(v) { tapply(v, seu$leiden, function(vv) { 100 * sum(vv > 0) / length(vv) }) }))
	le_umifrac = sca_mc_gene_umifrac_noobj(le_counts, multiplying_factor = 10000)
	
	# save expression
	# footprints
	saveRDS(ct_fp, sprintf("%s/dat.%s.expression.cts_fp.rds", out_fn, spi))
	saveRDS(mc_fp, sprintf("%s/dat.%s.expression.mcs_fp.rds", out_fn, spi))
	saveRDS(le_fp, sprintf("%s/dat.%s.expression.lei_fp.rds", out_fn, spi))

	# umifrac
	saveRDS(ct_umifrac, sprintf("%s/dat.%s.expression.cts_umifrac.rds", out_fn, spi))
	saveRDS(mc_umifrac, sprintf("%s/dat.%s.expression.mcs_umifrac.rds", out_fn, spi))
	saveRDS(le_umifrac, sprintf("%s/dat.%s.expression.lei_umifrac.rds", out_fn, spi))
		
	# umi counts
	saveRDS(ct_counts, sprintf("%s/dat.%s.expression.cts_umicount.rds", out_fn, spi))
	saveRDS(mc_counts, sprintf("%s/dat.%s.expression.mcs_umicount.rds", out_fn, spi))
	saveRDS(le_counts, sprintf("%s/dat.%s.expression.lei_umicount.rds", out_fn, spi))
		
	# cellfrac
	saveRDS(ct_fracell, sprintf("%s/dat.%s.expression.cts_fracell.rds", out_fn, spi))
	saveRDS(mc_fracell, sprintf("%s/dat.%s.expression.mcs_fracell.rds", out_fn, spi))
	saveRDS(le_fracell, sprintf("%s/dat.%s.expression.lei_fracell.rds", out_fn, spi))

    # fix gg_v
    prefix <- paste0(spi_w, "_")
    gg_v <- setNames(fix_gene_name(rownames(seu@assays$RNA$counts), prefix), rownames(seu@assays$RNA$counts))

    ## Find cluster markers ##
    # find markers for cell types
    SeuratObject::Idents(seu) = seu@meta.data[, "cell_type"]
    mks_cts = data.frame()
    for (loi in levels(seu@meta.data[, "cell_type"])) {
        message(sprintf("seurat | %s | markers cell_type %s...", spi, loi)) 
        mki = Seurat::FindMarkers(
            seu, 
            ident.1 = loi, 
            layer = "RNA", 
            min.pct = 0.05, 
            test.use = "wilcox"
        )
        if(spi_w == "Xesp"){
            rownames(mki) = symbol_to_id[rownames(mki)]
        }
        #mki = mki[!grepl("orphan", rownames(mki)), ]
        ct_umifrac_i_v = dic_from_vecs(rownames(ct_umifrac), ct_umifrac[, loi])
        mki$umifrac = sprintf("%.2f", ct_umifrac_i_v[rownames(mki)])
        mki$gene_seurat = rownames(mki)
        mki$focus_node = loi
        rownames(mki) = NULL   
        mks_cts = rbind(mks_cts, mki)
    }
    message(sprintf("  - Total markers: %d", nrow(mks_cts)))
    mks_cts$gene = gg_v[mks_cts$gene_seurat]
    mks_cts = mks_cts %>% left_join(gene_annot, by = "gene")
    write.table(
        mks_cts, 
        sprintf("%s/seu.%s.markers.cts.csv", out_fn, spi), 
        sep = "\t", 
        col.names = TRUE, 
        row.names = FALSE, 
        quote = FALSE
    )
        message(sprintf("  - Saved: %s/seu.%s.markers.cts.csv", out_fn, spi))

    # calculate Log2FC for all genes across cell types
    SeuratObject::Idents(seu) = seu@meta.data[, "cell_type"]
    log2fc_cts = data.frame()
    for (loi in levels(seu@meta.data[, "cell_type"])) {
        message(sprintf("seurat | %s | log2fc cell_type %s...", spi, loi)) 
        mki = Seurat::FindMarkers(
            seu, 
            ident.1 = loi,
            layer = "RNA",
            logfc.threshold = 0,
            min.pct = 0,
            test.use = "wilcox"
        ) 
        if(spi_w == "Xesp"){
            rownames(mki) = symbol_to_id[rownames(mki)]
        }
        #mki = mki[!grepl("orphan", rownames(mki)), ]   
        ct_umifrac_i_v = dic_from_vecs(rownames(ct_umifrac), ct_umifrac[, loi])
        mki$umifrac = sprintf("%.2f", ct_umifrac_i_v[rownames(mki)])
        mki$focus_node = loi
        mki$gene_seurat = rownames(mki)
        rownames(mki) = NULL
        log2fc_cts = rbind(log2fc_cts, mki)
    }
    message(sprintf("  - Total lines: %d", nrow(log2fc_cts)))
    log2fc_cts$gene = gg_v[log2fc_cts$gene_seurat]
    log2fc_cts = log2fc_cts %>% left_join(gene_annot, by = "gene")
    write.table(
        log2fc_cts, 
        sprintf("%s/seu.%s.log2fc.cts.csv", out_fn, spi), 
        sep = "\t", 
        col.names = TRUE, 
        row.names = FALSE, 
        quote = FALSE
    )
        message(sprintf("  - Saved: %s/seu.%s.log2fc.cts.csv", out_fn, spi))

    # find markers for metacells
    SeuratObject::Idents(seu) = seu@meta.data[, "metacell"]
    mks_mcs = data.frame()
    for (loi in sort(unique(seu@meta.data[, "metacell"]))) {
        message(sprintf("seurat | %s | markers metacell %s...", spi, loi))
        mki = Seurat::FindMarkers(
            seu, 
            ident.1 = loi, 
            assay = "RNA", 
            slot = "counts", 
            fc.slot = "counts", 
            min.pct = 0.05, 
            test.use = "wilcox"
        ) 
        #mki = mki[!grepl("orphan", rownames(mki)), ]  
        mc_umifrac_i_v = dic_from_vecs(rownames(mc_umifrac), mc_umifrac[, loi])
        mki$umifrac = sprintf("%.2f", mc_umifrac_i_v[rownames(mki)])
        mki$focus_node = loi
        mki$gene_seurat = rownames(mki)
        rownames(mki) = NULL   
        mks_mcs = rbind(mks_mcs, head(mki, 50))
    }
    mks_mcs$gene = gg_v[mks_mcs$gene_seurat]
    mks_mcs = mks_mcs %>% left_join(gene_annot, by = "gene")
    write.table(
        mks_mcs, 
        sprintf("%s/seu.%s.markers.mcs.csv", out_fn, spi), 
        sep = "\t", 
        col.names = TRUE, 
        row.names = FALSE, 
        quote = FALSE
    )
        message(sprintf("  - Saved: %s/seu.%s.markers.mcs.csv", out_fn, spi))


    # find markers for leidens
    SeuratObject::Idents(seu) = seu@meta.data[, "leiden"]
    mks_lds = data.frame()
    for (loi in sort(unique(seu@meta.data[, "leiden"]))) {
        message(sprintf("seurat | %s | markers leiden %s...", spi, loi))
        mki = Seurat::FindMarkers(
            seu, 
            ident.1 = loi, 
            assay = "RNA", 
            slot = "counts", 
            fc.slot = "counts", 
            min.pct = 0.05, 
            test.use = "wilcox"
        ) 
        #mki = mki[!grepl("orphan", rownames(mki)), ]  
        ld_umifrac_i_v = dic_from_vecs(rownames(le_umifrac), le_umifrac[, loi])
        mki$umifrac = sprintf("%.2f", ld_umifrac_i_v[rownames(mki)])
        mki$focus_node = loi
        mki$gene_seurat = rownames(mki)
        rownames(mki) = NULL   
        mks_lds = rbind(mks_lds, head(mki, 50))
    }
    mks_lds$gene = gg_v[mks_lds$gene_seurat]
    mks_lds = mks_lds %>% left_join(gene_annot, by = "gene")
    write.table(
        mks_lds, 
        sprintf("%s/seu.%s.markers.lds.csv", out_fn, spi), 
        sep = "\t", 
        col.names = TRUE, 
        row.names = FALSE, 
        quote = FALSE
    )
        message(sprintf("  - Saved: %s/seu.%s.markers.lds.csv", out_fn, spi))
    ## Build cell type tree ##
	# tree of cell types
	SeuratObject::Idents(seu) = seu@meta.data[, "cell_type"]
	seu = Seurat::BuildClusterTree(seu, assay = "RNA", reduction = "pca", dims = 1:10)
	seu_phylo = Seurat::Tool(seu, slot = "Seurat::BuildClusterTree")
	tip_color_dictionary = dic_from_vecs(seu$cell_type, seu$color)
	mks_bipartition_cts = sca_seurat_tree_markers(
		seu = seu, seu_phylo = seu_phylo,
		output_file = sprintf("%s/seu.%s.tree.markers.cts.pdf", out_fn, spi),
		tip_color_dictionary = tip_color_dictionary,
		assay = "RNA",
		test = "wilcox",
		pval_thr = 1e-3)  
	mks_bipartition_cts$gene = gg_v[mks_bipartition_cts$gene]
	mks_bipartition_cts = mks_bipartition_cts %>% 
        left_join(gene_annot %>% 
        dplyr::select("gene", "spis_ortholog", "orthogroup_name", "Uniprot_name", "is_tf"), 
        by = "gene")
	write.table(mks_bipartition_cts, sprintf("%s/seu.%s.tree.markers.cts.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    
    ## Celltype heatmap ##
    message(sprintf("metacell | %s | expression maps, celltype heatmaps...", spi))
    # 选择 top markers
    top_markers = select_top_markers(ct_fp, matrix_thr = 1.5, n_top_markers = 30, n_markers_rollmean = 2)
    ct_fp_f = ct_fp[top_markers, ]
    order_rows = rev(names(sort(apply(ct_fp_f, 1, function(x) which.max(x)))))
    ct_fp_f = ct_fp_f[order_rows, ]
    top_cl_per_gene = colnames(ct_fp_f)[apply(ct_fp_f, 1, which.max)]

    # integrat gene annotation
    gene_seurat = rownames(ct_fp_f)
    gene_std = gg_v[gene_seurat]
    idx = match(gene_std, gene_annot$gene)
    gene_names = gene_annot$orthogroup_name[idx]
    gene_Uniprot = gene_annot$Uniprot_name[idx]
    row_labels = paste(gene_std, gene_Uniprot, gene_names, sep = " | ")
    rownames(ct_fp_f) = stringr::str_trunc(row_labels, 80)

    pp1 = plot_complex_heatmap(
        ct_fp_f,
        name = "Footprint",
        cluster_row = FALSE,
        cluster_col = FALSE,
        color_min = 1,
        color_max = 8,
        color_mat = c("gray95", "orange", "orangered2", "#520c52"),
        cell_border = gpar(col = NA, lwd = 1, lty = 1),
        heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
        categories_row = top_cl_per_gene,
        colors_row = ctt_cts_col_v,
        categories_col = colnames(ct_fp_f),
        colors_col = ctt_cts_col_v,
        fontsize = 5
    )
    pdf(sprintf("%s/seu.%s.markers.top.cts.pdf", out_fn, spi), 
        width = 14 + (ncol(ct_fp_f) / 30), 
        height = 10 + (nrow(ct_fp_f) / 20))
    print(pp1)
    dev.off()
    message(sprintf("  Saved: %s/seu.%s.markers.top.mcs.pdf", out_fn, spi))

    ## Leiden clusters heatmap ##
    message(sprintf("metacell | %s | expression maps, Leiden heatmaps...", spi))
    # 选择 top markers
    top_markers = select_top_markers(le_fp, matrix_thr = 1.5, n_top_markers = 20, n_markers_rollmean = 3)
    le_fp_f = le_fp[top_markers, ctt$cluster[ctt$cluster_class == "leiden"]]
    order_rows = rev(names(sort(apply(le_fp_f, 1, function(x) which.max(x)))))
    le_fp_f = le_fp_f[order_rows, ]
    top_cl_per_gene = colnames(le_fp_f)[apply(le_fp_f, 1, which.max)]

    # integrat gene annotation
    gene_seurat = rownames(le_fp_f)
    gene_std = gg_v[gene_seurat]
    idx = match(gene_std, gene_annot$gene)
    gene_names = gene_annot$orthogroup_name[idx]
    gene_Uniprot = gene_annot$Uniprot_name[idx]
    row_labels = paste(gene_std, gene_Uniprot, gene_names, sep = " | ")
    rownames(le_fp_f) = stringr::str_trunc(row_labels, 80)

    # 绘制热图
    pp1 = plot_complex_heatmap(
        le_fp_f,
        name = "Footprint",
        cluster_row = FALSE,
        cluster_col = FALSE,
        color_min = 1,
        color_max = 8,
        color_mat = c("gray95", "orange", "orangered2", "#520c52"),
        cell_border = gpar(col = NA, lwd = 1, lty = 1),
        heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
        categories_row = top_cl_per_gene,
        colors_row = ctt_lou_col_v,
        categories_col = colnames(le_fp_f),
        colors_col = ctt_lou_col_v,
        fontsize = 5
    )

    pdf(sprintf("%s/seu.%s.markers.top.lei.pdf", out_fn, spi), 
        width = 10 + (ncol(le_fp_f) / 30), 
        height = 10 + (nrow(le_fp_f) / 20))
    print(pp1)
    dev.off()
    message(sprintf("  Saved: %s/seu.%s.markers.top.lei.pdf", out_fn, spi))

    ## Metacells heatmap ##
    message(sprintf("metacell | %s | expression maps, metacell heatmaps...", spi))
    # 选择 top markers
    top_markers = select_top_markers(mc_fp, matrix_thr = 1.5, n_top_markers = 20, n_markers_rollmean = 4)
    mc_fp_f = mc_fp[top_markers, ctm$metacell]
    order_rows = rev(names(sort(apply(mc_fp_f, 1, function(x) which.max(x)))))
    mc_fp_f = mc_fp_f[order_rows, ]
    top_cl_per_gene = colnames(mc_fp_f)[apply(mc_fp_f, 1, which.max)]

    # integrat gene annotation
    gene_seurat = rownames(mc_fp_f)
    gene_std = gg_v[gene_seurat]
    idx = match(gene_std, gene_annot$gene)
    gene_names = gene_annot$orthogroup_name[idx]
    gene_Uniprot = gene_annot$Uniprot_name[idx]
    row_labels = paste(gene_std, gene_Uniprot, gene_names, sep = " | ")
    rownames(mc_fp_f) = stringr::str_trunc(row_labels, 80)

    # 绘制热图
    pp1 = plot_complex_heatmap(
        mc_fp_f,
        name = "Footprint",
        cluster_row = FALSE,
        cluster_col = FALSE,
        color_min = 1,
        color_max = 8,
        color_mat = c("gray95", "orange", "orangered2", "#520c52"),
        cell_border = gpar(col = NA, lwd = 1, lty = 1),
        heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
        categories_row = top_cl_per_gene,
        colors_row = ctm_mcs_col_v,
        categories_col = colnames(mc_fp_f),
        colors_col = ctm_mcs_col_v,
        fontsize = 5
    )
    pdf(sprintf("%s/seu.%s.markers.top.mcs.pdf", out_fn, spi), 
        width = 10 + (ncol(mc_fp_f) / 30), 
        height = 10 + (nrow(mc_fp_f) / 20))
    print(pp1)
    dev.off()
    message(sprintf("  Saved: %s/seu.%s.markers.top.mcs.pdf", out_fn, spi))
    
    message(sprintf("metacell | %s | save final object...", spi))

    saveRDS(seu, sprintf("%s/dat.%s.seurat_final.rds", out_fn, spi))
    write.table(
        seu@meta.data,
        sprintf("%s/dat.%s.cell_metadata.csv", out_fn, spi),
        sep = "\t",
        col.names = TRUE,
        row.names = TRUE,
        quote = FALSE
    )
}
