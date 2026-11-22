# Coral cross-species single-cell atlas

Code and selected analysis inputs and results accompanying **A comparative single-cell atlas of 11 corals reveals conserved and divergent programs in specialized cells**.

The study combines nine newly generated coral single-cell datasets with published datasets from *Xenia* sp. and *Stylophora pistillata*. It compares cell-type identities and transcriptional programs across Hexacorallia and Octocorallia, with a focus on cnidocytes, calcifying cells and alga-hosting gastrodermal cells (symbiocytes). Orthogroup-level expression, gene-family phylogenies, ancestral-state reconstruction and promoter motifs provide complementary views of conserved functions and lineage-associated differences.

## Getting started

Choose a module below and follow its README for inputs, software, working directory and execution order. Modules can start from supplied reference results or from their specified upstream inputs.

The atlas uses `Adig`, `Ahya`, `Amur`, `Fspp`, `Gfas`, `Hcoe`, `Pcru`, `Pdam`, `Plut`, `Spin` and `Xesp`. QC and clustering use `Fusp` for the dataset subsequently named `Fspp`. `Spin` denotes the *S. pistillata* reference object; orthology tables and gene prefixes use `Spis`. The ten query datasets undergo clustering and metacell construction; the published Spis reference retains its existing clusters and metacells.

## Analysis modules

| Theme | Module | Analysis and supplied materials |
|---|---|---|
| Data processing | [3′ annotation extension](data_processing/extension3p/README.md) | Strand-specific peaks, annotation extension and final counting GTF |
| Data processing | [Read processing and quantification](data_processing/read_processing_quantification/README.md) | Read preparation, STAR alignment and dropEst count matrices |
| Data processing | [QC](data_processing/qc/README.md) | Droplet filtering, ambient-RNA correction, cell filtering and doublet removal for the nine new datasets |
| Data processing | [Clustering stability](data_processing/clustering_stability/README.md) | Resolution scan, Pilot/Formal stability analysis, selection and object rebuilding |
| Data processing | [Metacells](data_processing/metacell/README.md) | Balanced coclustering of the ten query datasets |
| Gene annotation | [Orthogroups and TFs](gene_annotation/orthogroup_inference/README.md) | Anthozoa/Metazoa protein and GTF inputs, Broccoli conversion, TF inference and adopted OG/pair tables |
| Gene annotation | [Functional annotation](gene_annotation/functional_annotation/README.md) | UniProt-based gene annotation and adopted gene/OG annotation tables |
| Atlas comparison | [SAMap](atlas_comparison/samap/README.md) | Homology-map preparation, a pairwise reference-mapping notebook and an eleven-species analysis script |
| Atlas comparison | [Pairwise SATURN](atlas_comparison/saturn_umap/README.md) | Ten adopted query–Spis UMAP figures |
| Atlas comparison | [Joint SATURN](atlas_comparison/saturn_joint/README.md) | Supplied training-to-figure script for joint integration |
| Atlas comparison | [Cell-type annotation evidence](atlas_comparison/celltype_annotation_evidence/README.md) | Adopted annotations and saved SAMap, footprint-correlation and SATURN evidence |
| Atlas comparison | [Expression, ICC and cell-type tree](atlas_comparison/expression_icc_tree/README.md) | Gene expression and footprints, ten query–Spin comparisons, and the adopted cell-type tree |
| Expression analysis | [OG_SUM conservation](expression_analysis/og_sum_conservation/README.md) | OG count aggregation, conserved enrichment, eleven adopted OG_SUM objects and expression statistics |
| Expression analysis | [TF conservation](expression_analysis/tf_conservation/README.md) | TF-associated OG selection, expression statistics, species support, member contributions and heatmaps |
| Expression analysis | [Calicoblast module scores](expression_analysis/calicoblast_module_scores/README.md) | 44/266-OG scoring with shared frozen controls and weights; supplied M3/M4 plotting inputs |
| Gene-family evolution | [Calicoblast family trees](gene_family_evolution/calicoblast_multi_og/README.md) | CA, SLC4, ferlin and astacin inference from frozen members, expression/domain integration, 138-tip main and 194-tip comparison displays |
| Gene-family evolution | [Dollo reconstruction](gene_family_evolution/dollo/README.md) | Separate OG-membership and cell-type-enrichment histories, node summaries and annotated OG lists |
| Regulatory analysis | [MITF-family motifs](regulatory_analysis/mitf/README.md) | Xenia promoter extraction, motif scanning, enrichment and motif-specificity comparisons |

## Connections between modules

| From | To | Files or preparation passed between them |
|---|---|---|
| 3′ extension and read processing | QC | Final GTF for counting; unfiltered dropEst `cell.counts.rds` with its `cm` matrix for QC |
| QC → stability selection | Metacells | Post-QC objects, then selected-resolution objects carrying the clustering/PCA fields required by metacell construction |
| Orthogroup and functional annotation | Expression and OG-based analyses | Anthozoa memberships/pairs, TF-integrated Metazoa memberships and gene/OG annotations; some scripts require setting their input paths to these supplied tables |
| Metacells and adopted cell-type labels | Expression/ICC/tree | `<species>_ld_mc.rds`, per-species cell annotations and gene tables prepared by the expression script |
| OG_SUM | TF conservation | Shared objects, feature metadata, gene–OG mapping and formal marker tables for calicoblast, alga-hosting gastrodermis and cnidocyte comparisons |
| Gene-level expression statistics | Dollo and promoter analysis | Dollo uses log2FC tables and cell-type UMI summaries; MITF uses the Xesp marker table together with its genome and GTF |

The expression script combines the Spis gene and OG annotations and writes `annot.Spin.genes.csv` for the Spin reference, preserving the Spis gene identifiers. The annotation module's `cell_type` field contains the adopted label.

OG_SUM compares expression summed over OG members. Dollo enrichment instead selects OGs from qualifying **individual-gene** statistics. The family-tree module starts from frozen members and adopted expression/domain tables. These inputs should not be substituted for one another. The 44-OG and 266-OG scoring branches both reuse the same frozen controls and weights.

## Data and software

Data records listed in the manuscript are [NCBI BioProject PRJNA1462667](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1462667) for raw sequencing reads and [Mendeley Data, version 1](https://doi.org/10.17632/mzt7wmxmys.1
        
        
        
        
        
        
        
        ) for processed count matrices, cell-type annotations and metadata.

This code collection also includes selected protein/GTF inputs, orthology and annotation tables, adopted OG_SUM objects, scoring inputs, trees and figures. Prepare the additional Seurat/h5ad objects, gene statistics, genomes and reference databases specified by each module. Set the input and executable paths specified in each module. Dollo reads the gene-level log2FC tables and UMI count matrices generated by the expression module directly from its species directories.

Use the software requirements in each module README. R/Seurat analyses, Python-based integration, sequence-analysis tools and motif analysis have different dependencies.

## Attribution

Several workflows and helper libraries draw on [Levy et al. (2025)](https://doi.org/10.1038/s41586-025-09) and the [oculina-coral-sc-atlas repository](https://github.com/sebepedroslab/oculina-coral-sc-atlas/) ([original repository](https://github.com/xgrau/oculina-coral-sc-atlas)). Method-specific references, upstream attribution and any supplied third-party licenses are retained in the respective modules.
