# Expression, ICC and cell-type tree

## Workflow

Annotated Seurat/metacell objects → group expression and footprints → pairwise iterative comparison of coexpression (ICC) against Spin → cross-species cell-type matrix and tree. ICC here is not an intraclass correlation coefficient. Within-species plots and marker calculations remain part of the expression script.

## Required files

R with Seurat/SeuratWrappers, doParallel/foreach, tidyverse, Matrix, ComplexHeatmap/circlize, ape/phangorn and the other packages loaded by `functions/`. Set the external paths for the Seurat objects in the expression script. ICC and tree construction read the supplied [Anthozoa gene pairs and Anthozoa/Metazoa membership tables](../../gene_annotation/orthogroup_inference/reference_results/).

For each species, use either the metacell object `<species>_ld_mc.rds` or the downloadable annotated object `dat.<species>.seurat_final.rds` produced by this workflow. To use the downloaded objects, keep them in a separate input directory and change the `readRDS` filename pattern in the expression script to `/path/to/downloads/dat.%s.seurat_final.rds`. The script reads the existing `leiden` and `metacell` assignments and reapplies the supplied cell-type annotations.

The expression loop uses Adig, Ahya, Amur, Fspp, Gfas, Hcoe, Pcru, Pdam, Plut, Spin and Xesp. The Spis reference object is reused from the published single-cell atlas of [Levy et al. (2025)](https://doi.org/10.1038/s41586-025-09623-6) and supplied under the `Spin` input name, retaining its original `leiden` and `metacell` fields. The other ten species use Seurat objects produced after [stability-based resolution selection](../../data_processing/clustering_stability/README.md) and subsequent [metacell clustering](../../data_processing/metacell/README.md). Final cell annotations are read from [the supplied annotation module](../celltype_annotation_evidence/final_annotations/).

The eleven prepared gene annotation tables are included in [work/](work/), one per species at `results_metacell_<species>_filt/annot.<species>.genes.csv`. At startup, the expression script regenerates them in place from the [gene and OG reference tables](../../gene_annotation/functional_annotation/reference_results/), using the Spis inputs for Spin.

## Running the scripts

Start each script from this module directory; each changes its working directory to `work`:

```bash
mkdir -p work/cell_type_tree
Rscript scripts/scatlas_mutispecies_260818_reannotated.R
Rscript scripts/calculate_ICC_mutispecies_260818.R
```

The scripts include all ten query species: Adig, Ahya, Amur, Fspp, Gfas, Hcoe, Pcru, Pdam, Plut and Xesp. Expression is prepared for these species and the Spin reference; the ICC loop compares each query against Spin. The tree uses all eleven species. Once all required expression and ICC files exist, run:

```bash
Rscript scripts/cell_type_tree_with_check.R
```

Footprints use scaled RNA counts, geometric means and a 0.05 pseudocount; count sums, fractions and footprints are distinct outputs. ICC retains the ≥10 UMI input filter and its iterative settings. The tree retains EC ≥0, 0.7 orthology/coverage filters, FP ≥2.5, within-species quantile normalization, NA-to-zero handling, LogDet distances and UPGMA. The separate `csps_similarity_ct` branch is outside this module's scope.

## Final output

`work/results_metacell_*_filt/` contains annotated objects, expression matrices and plots; the Spin `csps/` subdirectory contains ICC scores/objects. `work/cell_type_tree/results_csps_cell_types_noknown_2.5_0.7_0.7_0_radian/` receives the matrix, Newick and tree figures. The adopted 162-tip tree and three PDFs are preserved in `reference_results/tree/`.

`scripts/scatlas_mutispecies_260818_reannotated.R` uses `FindMarkers` to compare each cell type with the remaining cells within the same species. Its `seu.<species>.markers.cts.csv` and `seu.<species>.log2fc.cts.csv` tables describe gene-level up- and downregulation and carry the primary and secondary OG labels, support values and TF flag from the prepared gene table. Genes assigned to several OGs retain their slash-separated IDs and have missing single-OG labels. The log2FC tables and the cell-type UMI count matrices (`dat.<species>.expression.cts_umicount.rds`) are read directly from these species directories by downstream [Dollo enrichment analysis](../../gene_family_evolution/dollo/README.md), which selects enriched genes from these statistics.

The retained bootstrap code requests 1,000 replicates, allows reduced successful batches and retains a fixed 1,000 denominator in manual support calculations; initial sampling is not seeded. Radial display branch lengths are adjusted and differ from the inferred tree. 

## Attribution

Most helper functions in `functions/` are reused from the [oculina-coral-sc-atlas repository](https://github.com/sebepedroslab/oculina-coral-sc-atlas/). The code structure and analysis workflow also draw on this repository.
