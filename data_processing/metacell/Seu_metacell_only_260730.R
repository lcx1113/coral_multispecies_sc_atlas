options(tgs_use.blas=TRUE)
# Usage:
#   Rscript Seu_metacell_only_260730.R /path/to/input.rds /path/to/output.rds
#
# This script accepts two command line arguments:
#   1) input RDS file (Seurat object)
#   2) output RDS file (where processed Seurat object will be saved)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  message("ERROR: Missing arguments.")
  message("Usage: Rscript Seu_metacell_only_260730.R /path/to/input.rds /path/to/output.rds")
  quit(status = 1)
}
input_rds <- args[1]
output_rds <- args[2]

if (!file.exists(input_rds)) {
  stop(sprintf("Input file does not exist: %s", input_rds))
}

library(Seurat)
library(tgstat)
library(ggplot2)
library(patchwork)
library(dplyr)
source("functions/Seurat_functions.R")
seu <- readRDS(input_rds)
seu$leiden <- seu$stability_selected_cluster
pca_mat <- t(Embeddings(seu[["stability_pca"]]))
pca_var_features <- rownames(pca_mat)[seu@commands$FindNeighbors.SCT.stability_pca@params$dims]

seu_mcs_pca_v <- sca_balanced_coclustering_headless(
  mat = pca_mat,
  variable_features = pca_var_features,
  mc_prefix = "mc.",
  min_mc_size = 10,
  top_K = 100,
  cgraph_K = 200,
  alpha_relaxation = 2,
  coclust_p_resamp = 0.75,
  coclust_n_resamp = 100
)

seu@meta.data$metacell <- NA
seu@meta.data$metacell <- seu_mcs_pca_v[rownames(seu@meta.data)]
table(seu$metacell)
print(paste("Number of metacells:", length(unique(seu$metacell[!is.na(seu$metacell)]))))

# Save to user-specified output path
saveRDS(seu, output_rds)
message("Saved processed Seurat object to: ", output_rds)