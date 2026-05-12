# coral_multispecies_sc_atlas
## Introduction
This repository contain code and some of the data in our mutispecies coral single cell transriptome analysis. 

### Files in the repository
- `data/`: protein sequence used, with broccoli code and results and TF annotation included.
- `result_annotation/`: species annotation including gene_name, Uniprot_name, and orthogroup information is provided in `result_annotation/gene_annotation/`; GO annotation in `result_annotation/GO_annotation`; Anthozoan orthogroup functional inference in `result_annotation/260428_orthogroup_annotate`; gene gain or loss infered by Anthozoan orhogroup with dollo in `result_annotation/ct_gene_evolution/`; name conversion table for each species between gene_id, protein_id and seurat_id in `result_annotation/260429_name_conversion`.
- `result_csps/`: Code to do SAMap comparison between two specise and to all species used in `result_csps/SAMap/`, with Code to do protein phylogenic and domain analysis in `result_csps/protein_phylogenetic_tree/`.
- `result_scatlas/`: Code associated with single-cell transcriptome analysis, including raw data processing with dropest in `result_scatlas/dropEst_pipline/`; 3'extension in `result_scatlas/extension_3p/`; quality control of single cell data in `result_scatlas/QC/`; analysis with scatlas in `result_scatlas/scatlas_mutispecies_260426/` and cross-species conserved gene in `result_scatlas/Integrate_crossspecies_marker/`.


## Credits

Parts of this analysis workflow were adapted from:
https://github.com/xgrau/oculina-coral-sc-atlas
We adapted their framework of processing cross-species comparison to single cell transcriptome and modified some upstream scripts to fit our dataset. Thanks to their contribution! 
