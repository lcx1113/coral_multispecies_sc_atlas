## Adapted from: https://github.com/xgrau/oculina-coral-sc-atlas
## Add ortholog_fraction to select ortholog present in selection fraction of species.
## With other modification to adapt to our dataset
######################################################################################
# libraries
setwd("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/scripts_260424/scatlas_mutispecies_260426")
getwd()
suppressMessages(source("../coral_cs_functions/gene-set-analysis.R"))
suppressMessages(source("../coral_cs_functions/helper.R"))
suppressMessages(library("Seurat"))
graphics.off()
suppressMessages(library("ape"))
library(dplyr)

# output
out_fn = "results_ct_evolution/"
dir.create(out_fn, showWarnings = FALSE)

# list of comparisons to perform
comp_list = list(
  "calicoblast_t" = list(Adig = "calicoblast", Ahya = "calicoblast",
                       Amur = "calicoblast", Fspp = "calicoblast",
                       Gfas = "calicoblast", Hcoe = "calicoblast",
                       Pcru = "calicoblast", Pdam = "calicoblast",
                       Plut = "calicoblast", Spin = "calicoblast")
)

# thresholds
thr_fps = 1.1
thr_pv  = 1e-3

# load orthology
ogm = read.table("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/final_data/broccoli/Metazoa_all_genome_260426/orthogroup_conservation_integrated.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/final_data/broccoli/Anthozoa_all_genome_260426/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_gov = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)
ogm_tv = dic_from_vecs(names = ogm$transcript, terms = ogm$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)

# based on markers fc
ufs_l = list()
umi_l = list()
# load seurat and umifrac
for (spi in c("Adig", "Ahya", "Amur", "Fspp", "Gfas", "Hcoe", "Pcru", "Pdam", "Plut", "Spin", "Xesp")) {
	message(sprintf("load %s", spi))
	ufs_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_umifrac.rds", spi, spi))
	umi_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_umicount.rds", spi, spi))
}

# loop over marker sets
for (nn in 1:length(comp_list)) {

	cti = names(comp_list)[[nn]]
	ctv = comp_list[[nn]]
	message(sprintf("evo | %s | init...", cti))
	
	# tree specific for each case...
	sps_list = names(comp_list[[nn]])

	# get markers
	mks_l = list()
	mks_t = data.frame()
	for (spi in sps_list) {
		
		if (spi == "Spin") {
			spi_w = "Spis"
		} else {
			spi_w = spi
		}

		message(sprintf("evo | %s | load markers %s...", cti, spi))

		# === 新增：直接读取log2fc表 ===
		spi_expression <- read.table(
			sprintf("results_metacell_%s_filt/seu.%s.log2fc.cts.csv", spi, spi),
			header = TRUE, sep = "\t"
		)

		# 只保留目标 cell type 的记录
		mks_i <- spi_expression[ spi_expression$focus_node %in% ctv[[spi]], ]

		# 过滤缺失
		mks_i <- mks_i[ !is.na(mks_i$gene) , ]

		# === 保留原orthogroup映射 ===
		mks_i$orthogroup = oga_gov [ mks_i$gene ]
		mks_i$gene_name = ogm_gv [ mks_i$gene ]
		mks_i$cell_type = ctv[[spi]][1]

		# 计算fc与阈值字段
		mks_i$fc = as.numeric(2^mks_i$avg_log2FC)

		# umi/umifrac 仍基于原UMI表计算
		umi_i_v = rowSums(as.matrix(umi_l[[spi]][,ctv[[spi]]]))
		names(umi_i_v) = fix_gene_name(names(umi_i_v), paste0(spi_w, "_"))
		ufs_i_v = umi_i_v / sum(umi_i_v) * 1e4
		mks_i$umi = umi_i_v [ mks_i$gene ]
		mks_i$umifrac = ufs_i_v [ mks_i$gene ]

		# store
		write.table(mks_i, sprintf("%s/mks.%s.%s.all.tsv", out_fn, cti, spi), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		mks_l[[spi]] = mks_i
		
		# merge with other species
		mks_f = mks_l[[spi]] [ , c("gene","orthogroup","umi","umifrac","fc","p_val_adj","pct.1") ]
		# add nonexpressed markers
		mks_f_nonexp = data.frame(gene = names(umi_i_v) [ !names(umi_i_v) %in% mks_f$gene ], orthogroup = oga_gov [ names(umi_i_v) [ !names(umi_i_v) %in% mks_f$gene ] ])
		mks_f_nonexp$umi = 0
		mks_f_nonexp$umifrac = 0
		mks_f_nonexp$fc = 0
		mks_f_nonexp$p_val_adj = 1
		mks_f_nonexp$pct.1 = 0
		mks_f = rbind(mks_f, mks_f_nonexp)
		# is significant?
		mks_f$significant = mks_f$p_val_adj < thr_pv & mks_f$fc >= thr_fps
		# add gene name
		mks_f$gene_name = ogm_gv [ mks_f$gene ]
		mks_f = mks_f [ !is.na(mks_f$gene_name), ]
		mks_f = mks_f [ order(mks_f$orthogroup, -as.numeric(mks_f$significant), mks_f$p_val_adj, -mks_f$fc) , ]

		# add species
		mks_f$species = spi
		mks_t = rbind(mks_t, mks_f)
	
	}
	# merged table
	mks_t$species = factor(mks_t$species, levels = sps_list)

	# load set of markers for each species
	ogs_d = data.frame()
	gen_l = list()
	for (spi in sps_list) {
		message(sprintf("evo | %s | select markers %s...", cti, spi))
		# save og+gene info
		gen_i = mks_l[[spi]] [ mks_l[[spi]]$fc > thr_fps & mks_l[[spi]]$p_val_adj < thr_pv ,]
		gen_i = gen_i [ order(gen_i$orthogroup, -gen_i$fc), ]
		gen_l[[spi]] = gen_i
		write.table(gen_i, sprintf("%s/mks.%s.%s.filt.tsv", out_fn, cti, spi), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		# store og info only
		ogs_d_i = data.frame(orthogroup = oga_gov [ gen_i$gene ], gene = gen_i$gene, species = spi)
		ogs_d_i = ogs_d_i [ !is.na(ogs_d_i$orthogroup), ]
		ogs_d = rbind(ogs_d, ogs_d_i)
	}

	# save
	write.table(ogs_d, sprintf("%s/anc.%s.mat.extant.tsv", out_fn, cti), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
	# resort merged table
	mks_t_o = unique(mks_t [ , c("orthogroup","significant","species") ])
	mks_t_median_uf_per_og = aggregate(umifrac ~ orthogroup, mks_t, mean)
	mks_t_numsignif_v = table(mks_t_o$orthogroup, mks_t_o$significant)[,"TRUE"]
	mks_t_median_uf_per_og$numsignif = mks_t_numsignif_v [ mks_t_median_uf_per_og$orthogroup ]
	mks_t_median_uf_per_og = mks_t_median_uf_per_og [ order(-mks_t_median_uf_per_og$numsignif, -mks_t_median_uf_per_og$umifrac), ]
	mks_t$orthogroup = factor(mks_t$orthogroup, levels = mks_t_median_uf_per_og$orthogroup)
	# filter ogs that are in markers list, and also best gene per og
	mks_t = mks_t [ order(mks_t$orthogroup, mks_t$species), ]
	mks_t_f = mks_t [ mks_t$orthogroup %in% ogs_d$orthogroup, ] 
	mks_t_f_best = mks_t_f [ !duplicated(paste(mks_t_f$orthogroup, mks_t_f$species)), ]

	# save
	write.table(mks_t_f, sprintf("%s/mks.%s.allsps.all_genes_per_og.tsv", out_fn, cti), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	write.table(mks_t_f_best, sprintf("%s/mks.%s.allsps.top_genes_per_og.tsv", out_fn, cti), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
	# loop over pairs of species (unchanged)
	for (spi in sps_list) {
		
		pdf(sprintf("%s/mks.%s.%s.all.pdf", out_fn, cti, spi), width = 20, height = 4)
		layout(matrix(1:5, nrow = 1))
		
		# volcano plot sps i
		mks_i = mks_l[[spi]]
		gen_i = gen_l[[spi]]
		mks_i = mks_i [ order(-mks_i$p_val_adj), ]
		# cap padj
		mks_i_logpadj = -log10(mks_i$p_val_adj)
		mks_i_logpadj [ mks_i_logpadj > 12 ] = 12
		# hist(mks_i_logpadj)
		mks_i_logpadj_v = cut(mks_i_logpadj, c(0, 3, 6, 9, 12), include.lowest = TRUE)
		color_palette = colorRampPalette(interpolate="l",c("gray90", "#9ae445","#44aa2c","#1a622d"))
		mks_i_logpadj_c = color_palette(nlevels(mks_i_logpadj_v))  [ mks_i_logpadj_v ]
		# cap fc
		mks_i_logfc = mks_i$avg_log2FC
		mks_i_logfc [ mks_i_logfc > 12 ] = 12
		mks_i_logfc [ mks_i_logfc < -12 ] = -12
		# cap umi
		mks_i_umi = mks_i$umi
		mks_i_umi [ mks_i_umi > 1e4 ] = 1e4
		plot(mks_i_logfc, mks_i_umi, col = mks_i_logpadj_c, pch = 19, ylim = c(1,1e4), xlim = c(-12,12), log = "y", las = 1)
		title(main = sprintf("%s\nFC v UMI in %s", spi, cti), cex.main = 1, font.main = 1)
		abline(v = c(-log2(1.1),0,log2(1.1)), lty = 2, col = "darkred")
		
		# plot fc v pvalue
		plot(mks_i$avg_log2FC, mks_i_logpadj, col = mks_i_logpadj_c, pch = 19, las = 1, xlim = c(-12,12), ylim = c(0,12))
		title(main = sprintf("%s\nFC v UMI in %s", spi, cti), cex.main = 1, font.main = 1)
		abline(v = c(-log2(1.1),0,log2(1.1)), lty = 2, col = "darkred")
		
		layout(matrix(1:5, nrow = 1))
		for (spj in sps_list) {
			if (spi != spj) {
				message(sprintf("evo | %s | compare %s-%s...", cti, spi,spj))
				# load
				gen_i = gen_l[[spi]]
				gen_j = gen_l[[spj]]
				mks_i = mks_l[[spi]]
				mks_j = mks_l[[spj]]
				# sort
				gen_i = gen_i [ !duplicated(gen_i$orthogroup), ]
				mks_i = mks_i [ order(mks_i$orthogroup, mks_i$p_val_adj, -mks_i$fc), ]
				mks_i = mks_i [ !duplicated(mks_i$orthogroup), ]
				gen_j = gen_j [ !duplicated(gen_j$orthogroup), ]
				mks_j = mks_j [ order(mks_j$orthogroup, mks_j$p_val_adj, -mks_j$fc), ]
				mks_j = mks_j [ !duplicated(mks_j$orthogroup), ]
				# merge
				mks_m = merge(mks_i[,c("orthogroup","avg_log2FC","p_val_adj","gene")],  mks_j[,c("orthogroup","avg_log2FC","p_val_adj","gene")], by = "orthogroup", all.x = TRUE, all.y = TRUE, suffixes = c("_i","_j"))
				# OG classification
				mks_m$og_class = factor(NA, levels = c("i-only","j-only","shared"))
				mks_m$og_class [ mks_m$orthogroup %in% setdiff(gen_i$orthogroup, gen_j$orthogroup) ] = "i-only"
				mks_m$og_class [ mks_m$orthogroup %in% setdiff(gen_j$orthogroup, gen_i$orthogroup) ] = "j-only"
				mks_m$og_class [ mks_m$orthogroup %in% intersect(gen_i$orthogroup, gen_j$orthogroup) ] = "shared"
				mks_m_og_class_c = c("mediumpurple1","plum2","blue4")
				mks_m = mks_m [ order(mks_m$og_class), ]
				mks_m = mks_m [ mks_m$gene_i %in% gen_i$gene | mks_m$gene_j %in% gen_j$gene, ]
				names(mks_m_og_class_c) = levels(mks_m$og_class)				
				plot(mks_m$avg_log2FC_i, mks_m$avg_log2FC_j, col = mks_m_og_class_c[mks_m$og_class], pch = 19, xlab = spi, ylab = spj, cex = 0.8)
				abline(v = c(0,log2(1.1)), h = c(0,log2(1.1)), lty = 2, col = "darkred")
				title(main = sprintf("%s\nFC %s-%s", cti, spi, spj), cex.main = 1, font.main = 1)
			} else {
				plot(0)
				legend("topleft", c("sps i-only","sps j-only","shared"), col = c("mediumpurple1","plum2","blue4"), pch = 19, bty = "n")
			}
		}
		dev.off()
	}
}

message("All done!")
