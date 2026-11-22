# Metacell construction

## Workflow

Selected Seurat objects from clustering stability → balanced coclustering using the selected PCA dimensions → metacell labels stored in the same object. The script uses `min_mc_size = 10`, `top_K = 100`, `cgraph_K = 200`, `cgraph_alpha = 2`, resampling fraction 0.75 and 100 resamples. Unassigned cells remain unassigned; no random seed is added by this entry point.

## Required files

R with Seurat, tgstat, ggplot2, patchwork, dplyr and the packages used by the supplied functions. Inputs must contain `stability_selected_cluster`, `stability_pca` and the stored `FindNeighbors.SCT.stability_pca` dimensions. [Input/output mapping](inputs_outputs.tsv) lists the adopted object names and reference counts; large RDS objects are external.

## Running the script

Run from this module directory, with a writable output directory:

```bash
mkdir -p output
Rscript Seu_metacell_only_260730.R /absolute/path/selected_objects/Adig_selected_r1.4_reclustered.rds output/Adig_ld_mc.rds
```

The same script applies to the ten query datasets. Preserve the mapping in the table: Fusp is exported as Fspp; Pcru uses the same name as its clustering input. Spis is a reused reference object and is not constructed by this command; obtain that object separately. The input/output table covers only the ten constructed query datasets.

## Final output

A Seurat RDS containing `leiden` and `metacell` assignments. Group expression matrices are produced separately in [expression/ICC/tree](../../atlas_comparison/expression_icc_tree/README.md). The full supplied function library and its existing attribution comments are retained.

## Attribution

The metacell code also draws on [Levy et al. (2025)](https://doi.org/10.1038/s41586-025-09623-6) and the [oculina-coral-sc-atlas code](https://github.com/sebepedroslab/oculina-coral-sc-atlas/) ([original repository](https://github.com/xgrau/oculina-coral-sc-atlas)).
