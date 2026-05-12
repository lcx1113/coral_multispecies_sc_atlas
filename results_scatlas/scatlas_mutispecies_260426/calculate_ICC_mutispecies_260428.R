## Adapted from: https://github.com/xgrau/oculina-coral-sc-atlas
##################################################################
#libraries
setwd("/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/scripts_260424/scatlas_mutispecies_260426/")
getwd()
library(doParallel)
registerDoParallel(cores=32)
suppressMessages(source("../coral_cs_functions/Seurat_functions.R"))
suppressMessages(source("../coral_cs_functions/gene-set-analysis.R"))
suppressMessages(source("../coral_cs_functions/helper.R"))
suppressMessages(source("../coral_cs_functions/Downstream_functions.R"))
suppressMessages(source("../coral_cs_functions/Cross_species_functions.R"))
suppressMessages(require("Seurat"))
suppressMessages(require("SeuratWrappers"))
graphics.off()
fix_gene_name <- function(gene_names, prefix) {
        gene_names<-gsub("-","_", gene_names)
        ifelse(grepl(paste0("^", prefix), gene_names), gene_names, paste0(prefix, gene_names))
    }

# data index
sps_list = c("Adig", "Ahya", "Amur", "Fspp", "Gfas", "Hcoe", "Pcru", "Pdam", "Plut", "Spin", "Xesp")
sps_refl = c("Spin")

for (spi in sps_refl){
    # set working species
    if (spi == "Spin") {
        spi_w = "Spis"
	} else if (spi == "HcoeH" || spi == "HcoeI") {
        spi_w = "Hcoe"
    } else {
		spi_w = spi
	}

    # output
	out_fn = sprintf("results_metacell_%s_filt/csps/", spi)
    dir.create(out_fn, showWarnings = FALSE)

    for (spj in sps_list [ sps_list != spi ]) {
        if (spj == "Spin") {
			spj_w = "Spis"
		}else if(spj == "HcoeH" || spj == "HcoeI") {
        spj_w = "Hcoe"
        }else{
			spj_w = spj
		}

        message(sprintf("csps | %s | compare to %s, load footprints...", spi, spj))
        inp_fn_i = sprintf("results_metacell_%s_filt/", spi)
        inp_fn_j = sprintf("results_metacell_%s_filt/", spj)
        mc_fp_i = readRDS(sprintf("%s/dat.%s.expression.mcs_fp.rds", inp_fn_i, spi))
        mc_fp_j = readRDS(sprintf("%s/dat.%s.expression.mcs_fp.rds", inp_fn_j, spj))

        message(sprintf("csps | %s | compare to %s, load counts...", spi, spj))
        mc_umi_i = readRDS(sprintf("%s/dat.%s.expression.mcs_umicount.rds", inp_fn_i, spi))
        mc_umi_j = readRDS(sprintf("%s/dat.%s.expression.mcs_umicount.rds", inp_fn_j, spj))

        rownames(mc_fp_i) = fix_gene_name(rownames(mc_fp_i), paste0(spi_w,"_"))
        rownames(mc_umi_i) = fix_gene_name(rownames(mc_umi_i), paste0(spi_w,"_"))
        rownames(mc_fp_j) = fix_gene_name(rownames(mc_fp_j), paste0(spj_w,"_"))
        rownames(mc_umi_j) = fix_gene_name(rownames(mc_umi_j), paste0(spj_w,"_"))
        message(sprintf("    gene format: %s", paste(head(rownames(mc_fp_i), 3), collapse = ", ")))
        message(sprintf("    gene format: %s", paste(head(rownames(mc_fp_j), 3), collapse = ", ")))

        message(sprintf("  Loading ortholog pairs: %s - %s", spi_w, spj_w))
        og_pairs_ij = clean_og_pairs(
            og_pairs_fn = "/home/data/Chenxi/workplace/scRNAseq/cross_species_scRNAseq/final_data/broccoli/Anthozoa_all_genome_260426/orthogroup_pairs.genes.csv", 
            sp1 = spi_w,
            sp2 = spj_w,
            t2g = FALSE,
            header = TRUE
        )

        message(sprintf("  Loaded %d ortholog pairs", nrow(og_pairs_ij)))
        message(sprintf("    %s format: %s", spi_w, paste(head(og_pairs_ij[,1], 3), collapse = ", ")))
        message(sprintf("    %s format: %s", spj_w, paste(head(og_pairs_ij[,2], 3), collapse = ", ")))

        message(sprintf("csps | %s | compare to %s, drop genes with too few counts...", spi, spj))
        
        counts_per_gene_i = rowSums(mc_umi_i)
        counts_per_gene_j = rowSums(mc_umi_j)
        good_genes_i = names(counts_per_gene_i) [ counts_per_gene_i >= 10 ]
        good_genes_j = names(counts_per_gene_j) [ counts_per_gene_j >= 10 ]
        good_genes_i = good_genes_i [ good_genes_i %in% rownames(mc_fp_i) ]
        good_genes_j = good_genes_j [ good_genes_j %in% rownames(mc_fp_j) ]
        message(sprintf("  %s: %d / %d genes pass filter (>= %d UMI)", spi, length(good_genes_i), nrow(mc_fp_i), 10))
        message(sprintf("  %s: %d / %d genes pass filter (>= %d UMI)", spj, length(good_genes_j), nrow(mc_fp_j), 10))
        og_pairs_ij_filt = og_pairs_ij[og_pairs_ij[, 1] %in% good_genes_i & og_pairs_ij[, 2] %in% good_genes_j, ]
        message(sprintf("  Ortholog pairs after filtering: %d (%.1f%% retained)", nrow(og_pairs_ij_filt), 100 * nrow(og_pairs_ij_filt) / nrow(og_pairs_ij)))
        message(sprintf("csps | %s | compare to %s, run ICC...", spi, spj))
        icc = csps_markers_icc_noobj(
            mat_sp1 = mc_fp_i[good_genes_i,],
            mat_sp2 = mc_fp_j[good_genes_j,],
            og_pairs = og_pairs_ij,
            do_duplicates = TRUE,
            nthreads_icc = 32,
            nthreads_dup = 1
        )
        write.table(icc$ec_markers, sprintf("%s/dat.icc.%s-%s.ec_scores.csv", out_fn, spi, spj), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        write.table(icc$ec_duplicates, sprintf("%s/dat.icc.%s-%s.ec_scores.duplicates.csv", out_fn, spi, spj), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        saveRDS(icc, sprintf("%s/dat.icc.%s-%s.obj.rds", out_fn, spi, spj))
        # Plot EC score distribution
        pdf(sprintf("%s/dat.icc.%s-%s.ec_distribution.pdf", out_fn, spi, spj), 
            width = 10, height = 6)
        par(mfrow = c(1, 2))
        hist(
            icc$ec_markers$ec_value,
            breaks = 50,
            main = sprintf("EC Score Distribution\n%s vs %s", spi, spj),
            xlab = "EC Score",
            ylab = "Frequency",
            col = "skyblue",
            border = "white"
        )
        abline(v = 0.5, col = "red", lwd = 2, lty = 2)
        abline(v = 0.7, col = "darkred", lwd = 2, lty = 3)
        legend("topright", 
            c("EC >= 0.5", "EC >= 0.7"), 
            col = c("red", "darkred"), 
            lty = c(2, 3), 
            lwd = 2)
        plot(
            ecdf(icc$ec_markers$ec_value),
            main = sprintf("Cumulative Distribution\n%s vs %s", spi, spj),
            xlab = "EC Score",
            ylab = "Cumulative Proportion",
            col = "blue",
            lwd = 2
        )
        abline(h = c(0.5, 0.95), col = "gray", lty = 2)
        abline(v = c(0.5, 0.7), col = c("red", "darkred"), lty = 2)
        dev.off()
        message(sprintf("  ✓ %s/dat.icc.%s-%s.ec_distribution.pdf", out_fn, spi, spj))
    }

}

message("All done!")
