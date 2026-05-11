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

# data index
sps_tree_s = c("((Xesp,Hcoe)N2,(((Pdam,Spis)N6,(Pcru,Fspp)N7)N4,(Plut,((Adig,(Ahya,Amur)N10)N9,Gfas)N8)N5)N3)N1;")
sps_tree = ape::read.tree(text = sps_tree_s)
# ape::write.tree(sps_tree, sprintf("%s/small_tree.newick", out_fn))
# sps_tree$edge.length = rep(1e-5, nrow(sps_tree$edge)) # add tiny edge lengths, to replicate dollo
# sps_tree_sansnodes = sps_tree
# sps_tree_sansnodes$node.label = NULL
descendants_from_nodes = adephylo::listTips(sps_tree)
nod_list = c(sps_tree$tip.label, sps_tree$node.label)
sps_list = sps_tree$tip.label

# list of comparisons to perform
comp_list = list(
  "calicoblast_t" = list(Adig = "calicoblast", Ahya = "calicoblast",
                       Amur = "calicoblast", Fspp = "calicoblast",
                       Gfas = "calicoblast", Hcoe = "calicoblast",
                       Pcru = "calicoblast", Pdam = "calicoblast",
                       Plut = "calicoblast", Spin = "calicoblast")
)

# load orthology
ogm = read.table("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/final_data/broccoli/Metazoa_all_genome_260426/orthogroup_conservation_integrated.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/final_data/broccoli/Anthozoa_all_genome_260426/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_gov = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)
ogm_tv = dic_from_vecs(names = ogm$transcript, terms = ogm$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)

# load evolutionary histories of all ogs
pr_ge_m = as.matrix(read.table(sprintf("%s/orthogroup_conservation.possvm.dollo_pres.csv", out_fn)))
ga_ge_m = as.matrix(read.table(sprintf("%s/orthogroup_conservation.possvm.dollo_gain.csv", out_fn)))
lo_ge_m = as.matrix(read.table(sprintf("%s/orthogroup_conservation.possvm.dollo_loss.csv", out_fn)))

# binarise
pr_ge_m = (pr_ge_m > 0) * 1
ga_ge_m = (ga_ge_m > 0) * 1
lo_ge_m = (lo_ge_m > 0) * 1

for (nn in 1:length(comp_list)) {

	cti = names(comp_list)[[nn]]
	ctv = comp_list[[nn]]
	message(sprintf("evo | %s | init...", cti))
	
	# species included in comparison
	sps_list_in_comp = names(comp_list[[nn]])
	
	# load set of markers for each species
	pr_ex_m = as.matrix(read.table(sprintf("%s/anc.%s.mat.possvm.dollo_pres.csv", out_fn, cti)))
	ga_ex_m = as.matrix(read.table(sprintf("%s/anc.%s.mat.possvm.dollo_gain.csv", out_fn, cti)))
	lo_ex_m = as.matrix(read.table(sprintf("%s/anc.%s.mat.possvm.dollo_loss.csv", out_fn, cti)))
	
	# output
	dir.create(sprintf("%s/tests_%s/", out_fn, cti), showWarnings = FALSE)
	
	# binarise
	pr_ex_m = (pr_ex_m > 0) * 1
	ga_ex_m = (ga_ex_m > 0) * 1
	lo_ex_m = (lo_ex_m > 0) * 1

	# restrict genome conservation info to ogs of interest
	pr_ge_m_i = pr_ge_m [ rownames(ga_ex_m), ]
	ga_ge_m_i = ga_ge_m [ rownames(ga_ex_m), ]
	lo_ge_m_i = lo_ge_m [ rownames(ga_ex_m), ]
	
	# matrix of complex events
	ev_ex_a = as.matrix(data.frame(
		expressed = colSums(pr_ex_m),
		newexp_oldgene = colSums(ga_ex_m==1 & ga_ge_m_i==0, na.rm = TRUE),
		newexp_newgene = colSums(ga_ex_m==1 & ga_ge_m_i==1, na.rm = TRUE),
		lossexp_presgene = colSums(lo_ex_m==1 & pr_ge_m_i==1, na.rm = TRUE),
		lossexp_lossgene = colSums(lo_ex_m==1 & pr_ge_m_i==0, na.rm = TRUE)
	))
	# drop gains of expression/conservation in ancestral node
	ev_ex_a [ nrow(ev_ex_a), "newexp_newgene" ] = 0
	
	# dataframes for plotting over the phylogeny
	df_ga = ev_ex_a[,c("newexp_newgene","newexp_oldgene","expressed")]
	df_ga[,"expressed"] = df_ga[,"expressed"] - df_ga[,"newexp_newgene"] - df_ga[,"newexp_oldgene"]
	df_ga = df_ga + 0.01
	co_ga = c("chartreuse4","olivedrab2","lightcyan2")
	df_lo = ev_ex_a[,c("lossexp_presgene","lossexp_lossgene")]+0.01
	co_lo = c("indianred1","indianred4")
	
	pdf(sprintf("%s/anc.%s.mat.possvm.tree.pdf", out_fn, cti), width = 7, height = 10)
	ape::plot.phylo(sps_tree, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex = 1, show.tip.label = FALSE, show.node.label = FALSE)
	ape::nodelabels(pie = df_ga[sps_tree$node.label,], piecol = co_ga, cex = sqrt(rowSums(df_ga[sps_tree$node.label,]) / 100))
	ape::tiplabels(pie = df_ga[sps_tree$tip.label,],   piecol = co_ga, cex = sqrt(rowSums(df_ga[sps_tree$tip.label,]) / 100))
	ape::nodelabels(pie = df_lo[sps_tree$node.label,], piecol = co_lo, cex = sqrt(rowSums(df_lo[sps_tree$node.label,]) / 100), adj = 0.1)
	ape::tiplabels(pie = df_lo[sps_tree$tip.label,],  piecol = co_lo, cex = sqrt(rowSums(df_lo[sps_tree$tip.label,]) / 100), adj = 0.1)
	ape::nodelabels(
		text = sprintf("%s\np=%i\n+%i(%i)\n-%i(%i)", sps_tree$node.label, ev_ex_a[sps_tree$node.label,"expressed"], colSums(ga_ex_m[,sps_tree$node.label]), floor(df_ga[sps_tree$node.label,"newexp_newgene"]), colSums(lo_ex_m[,sps_tree$node.label]), floor(df_lo[sps_tree$node.label,"lossexp_lossgene"])),
		col = scales::alpha("gray20", 0.9), frame="none", cex = 0.7)
	ape::tiplabels(
		text = sprintf("%s\np=%i\n+%i(%i)\n-%i(%i)", sps_tree$tip.label, ev_ex_a[sps_tree$tip.label,"expressed"], colSums(ga_ex_m[,sps_tree$tip.label]), floor(df_ga[sps_tree$tip.label,"newexp_newgene"]), colSums(lo_ex_m[,sps_tree$tip.label]), floor(df_lo[sps_tree$tip.label,"lossexp_lossgene"])),
		col = scales::alpha("gray20", 0.9), frame="none", cex = 0.7)
	title(main = sprintf("%s", cti))
	dev.off()

	ogs_annotation = read.table("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/scripts_260424/integrating_annotation/260428_orthogroup_annotate/orthogroup_final_annotations.csv", sep = "\t", header = TRUE, quote = "")
    for (noi in nod_list) {
		
		# ogs present, gained and lost
		ogs_pr_i = names(which(pr_ex_m[,noi] == 1))
		ogs_ga_i = names(which(ga_ex_m[,noi] == 1))
		ogs_lo_i = names(which(lo_ex_m[,noi] == 1))
		ogs_ga_genechange = names(which((ga_ex_m[, noi] == 1) & (ga_ge_m_i[, noi] == 1)))
  		ogs_lo_genechange = names(which((lo_ex_m[, noi] == 1) & (lo_ge_m_i[, noi] == 1)))

		ogs_pr_i = ogs_annotation %>% filter(orthogroup %in% ogs_pr_i)
		ogs_ga_i = ogs_annotation %>% filter(orthogroup %in% ogs_ga_i)%>%
    		mutate(gene_change = orthogroup %in% ogs_ga_genechange)
		ogs_lo_i = ogs_annotation %>% filter(orthogroup %in% ogs_lo_i)%>%
    		mutate(gene_change = orthogroup %in% ogs_lo_genechange)
		write.table(ogs_pr_i, sprintf("%s/tests_%s/tst.%s.%s.all.pres.oglist.csv", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(ogs_ga_i, sprintf("%s/tests_%s/tst.%s.%s.all.gain.oglist.csv", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(ogs_lo_i, sprintf("%s/tests_%s/tst.%s.%s.all.loss.oglist.csv", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		

    }
}



for (noi in nod_list) {
		
		# ogs present, gained and lost
		ogs_pr_i = names(which(pr_ex_m[,noi] == 1))
		ogs_ga_i = names(which(ga_ex_m[,noi] == 1))
		ogs_lo_i = names(which(lo_ex_m[,noi] == 1))
		ogs_ga_genechange = names(which((ga_ex_m[, noi] == 1) & (ga_ge_m_i[, noi] == 1)))
  		ogs_lo_genechange = names(which((lo_ex_m[, noi] == 1) & (lo_ge_m_i[, noi] == 1)))

		ogs_pr_i = ogs_annotation %>% filter(orthogroup %in% ogs_pr_i)
		ogs_ga_i = ogs_annotation %>% filter(orthogroup %in% ogs_ga_i)%>%
    		mutate(gene_change = orthogroup %in% ogs_ga_genechange)
		ogs_lo_i = ogs_annotation %>% filter(orthogroup %in% ogs_lo_i)%>%
    		mutate(gene_change = orthogroup %in% ogs_lo_genechange)
		write.table(ogs_pr_i, sprintf("%s/tests_%s/tst.%s.%s.all.pres.oglist.csv", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(ogs_ga_i, sprintf("%s/tests_%s/tst.%s.%s.all.gain.oglist.csv", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		write.table(ogs_lo_i, sprintf("%s/tests_%s/tst.%s.%s.all.loss.oglist.csv", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
		
		if (noi %in% sps_tree$tip.label) {
			spi_v = noi
		} else {
			spi_v = names(descendants_from_nodes[[noi]])
		}
		for (spi_w in spi_v) {
			go_annot_table = read.table(sprintf("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/scripts_260424/integrating_annotation/260506_get_Mmus_GO/GO_annotation/%s_orthogroup_with_Mmus_GO.tsv", spi_w),header = TRUE, sep = "\t", quote = "")
			annotation_df <- go_annot_table %>%
				select(gene, GO) %>%
				filter(!is.na(GO) & GO != "") %>%
				tidyr::separate_rows(GO, sep = ";") %>%
				distinct(gene, GO)
			go_annot <- annotation_df %>%
				group_by(gene) %>%
				summarise(GO = list(unique(GO)), .groups = "drop") %>%
				tibble::deframe()
			annotated_genes = unique(names(go_annot))

			# expand og to genes using ref sps
			gen_pr_d = oga [ oga$orthogroup %in% ogs_pr_i$orthogroup & oga$species %in% spi_w, c("gene","orthogroup") ]
			gen_ga_d = oga [ oga$orthogroup %in% ogs_ga_i$orthogroup & oga$species %in% spi_w, c("gene","orthogroup") ]
			gen_lo_d = oga [ oga$orthogroup %in% ogs_lo_i$orthogroup & oga$species %in% spi_w, c("gene","orthogroup") ]
			gen_pr_d = gen_pr_d [ gen_pr_d$gene %in% annotated_genes, ]
			gen_ga_d = gen_ga_d [ gen_ga_d$gene %in% annotated_genes, ]
			gen_lo_d = gen_lo_d [ gen_lo_d$gene %in% annotated_genes, ]


			list_tests = list(
				gain = gen_ga_d, pres = gen_pr_d, loss = gen_lo_d
			)
			
			for (ttt in 1:length(list_tests)) {
				# fetch
				ttlab = names(list_tests)[ttt]
				gen_x_d = list_tests[[ttt]]
				message(sprintf("evo | %s | test %s on %s branch, from %s...", cti, ttlab, noi, spi_w))
				
				# list of genes
				gen_x_i = gen_x_d [ !duplicated(gen_x_d$orthogroup), "gene" ]
				ttgen = gen_x_i

				if (length(ttgen) > 0) {
					# test GOs
					# message(sprintf("pathways %s | %s test %s GOs...", noi, cti, ttlab))
					go_enrichment_l = gsa_topgo_enrichment_v2(
						annotation = go_annot,
						genes_fg = ttgen,
						out_fn = sprintf("%s/tests_%s/tst.%s.%s.%s.%s.topgo.csv", out_fn, cti, cti, noi, spi_w, ttlab),
						ontologyset = c("BP","MF","CC"),
						tg_test = "fisher",
						tg_algorithm = "elim",
						top_markers = 30,
						nodesize = 10
					)
				}
			}
		}
    }
