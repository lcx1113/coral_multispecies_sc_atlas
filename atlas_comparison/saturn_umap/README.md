# Pairwise SATURN UMAP figures

## Workflow

Pairwise SATURN comparisons are used to compare our single-cell datasets with the annotated cell types in the published Spis single-cell atlas ([Levy et al. (2025)](https://doi.org/10.1038/s41586-025-09623-6)), providing reference evidence for cell-type annotation.

This module provides the adopted pairwise SATURN projections for Spis against ten query species. It is a figure collection; pairwise training and figure-generation scripts are outside the supplied scope.

## Required files

The ten PDFs and [figure index](figure_index.tsv) are included. Query species are Adig, Ahya, Amur, Fspp, Gfas, Hcoe, Pcru, Pdam, Plut and Xesp. No additional input or software beyond a PDF viewer is needed to inspect them.

## Viewing the results

Open the indexed `Spis_<query>/SATURN_leiden/umap_leiden_by_main_celltype_color.pdf`. Index counts describe the provided historical figures, not a new computation. Colours represent the supplied annotations rather than an independent cell-type prediction. The historical figure preparation used row order for annotation alignment after a barcode-name limitation; the PDFs alone cannot validate that alignment.

## Final output

The ten unchanged PDF figures. Derived proximity evidence is provided separately in [cell-type annotation evidence](../celltype_annotation_evidence/README.md); two-dimensional UMAP proximity must not be interpreted as a latent-space distance or a formal confidence score. Joint SATURN code is a separate [module](../saturn_joint/README.md), not a replacement for missing pairwise commands.
