# Cell-type annotations and supporting evidence

## Workflow

SAMap, ortholog footprint correlation and pairwise SATURN provide supporting evidence for manual cell-type annotation. This module supplies the adopted annotations and saved evidence tables; no analysis script needs to be run.

## Files

| Files | Contents |
|---|---|
| `final_annotations/annot.<species>.leiden.csv` | Adopted cluster labels and colours for ten query species and the Spin/Spis reference |
| `final_annotations/all_final_annotations.tsv` | Combined annotation table, with additional `species` and `role` columns |
| `celltype_annotation_evidence.tsv`, `per_species/` | Three-method evidence for 288 query clusters, as a combined table and ten species-specific tables |
| `SATURN_centroid_top2.tsv`, `SATURN_celltype_name_mapping.tsv` | Two nearest Spis clusters in the SATURN 2D UMAP, distances and label-name mappings |

All text tables are tab-delimited, including files ending in `.csv`. The combined annotations contain 288 query-cluster rows and 64 Spin reference rows (50 Leiden clusters and 14 metacells).

## Annotation fields and use

| Field | Meaning |
|---|---|
| `cluster` | Cluster or metacell identifier within a species |
| `cell_type` | Final adopted cell-type label |
| `color` | Display colour assigned to that label |
| `broad_cell_type` | Broader cell-type category |
| `cluster_class` | Whether the row describes a Leiden cluster or a metacell |
| `legacy_annotation_string` | Retained candidate/annotation text; use `cell_type` for the adopted label |

Use the per-species annotation files for downstream [expression and cell-type tree analysis](../expression_icc_tree/README.md). Match identifiers within each species and respect `cluster_class`; `Spin` denotes the Spis reference. The combined table adds `role` to distinguish query and reference records.

The eight-column evidence table contains `species`, `cluster`, SAMap's best type and saved similarity list, correlation's best type and saved similarity list, and SATURN's nearest and second-nearest types. Join it to query annotations by `species` and `cluster`. Evidence does not automatically replace the adopted labels. SATURN distances describe the 2D projection, and its two nearest clusters may share a cell type; `未记录` in SAMap fields means the candidate was not recorded, not a zero score.
