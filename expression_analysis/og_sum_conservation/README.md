# OG_SUM expression and conservation

## Workflow

OG_SUM sums the raw counts of genes belonging to the same orthogroup (OG) within each cell. This provides an OG-level basis for comparing functional programs of corresponding cell types across species and identifying conserved transcriptional signatures. OG expression summaries and cell-type log2FC values support cross-species comparisons and expression visualization.

RNA counts + gene–OG mapping → normalized OG_SUM objects → cell-type marker statistics → conservation tables and expression heatmaps.

Positive comparisons cover `cnidocyte`, `gastrodermis_alga_hosting` and `calicoblast` across 11 species, `calicoblast_stony` across nine (excluding Hcoe and Xesp), and `epidermis` across ten (excluding Xesp). Supplied downregulated results cover the first four comparisons. Species, targets and parameters are defined in `config/`.

## Required files

Use R with Seurat/SeuratObject, Matrix, data.table, yaml, digest and pheatmap, with `Rscript` available.

| Input | Preparation |
|---|---|
| Eleven RNA Seurat objects, for rebuilding OG_SUM | Download from the public data release and set paths in `config/species_manifest.tsv`. Required: RNA counts, matching feature identifiers and `cell_type` metadata. `Spin` corresponds to `Spis` in the orthology tables. |
| Eleven adopted OG_SUM objects | Included at `intermediate/02_og_sum_objects/<species>.OG_SUM.seurat.rds` for OG-level comparisons and downstream TF statistics. |
| Seurat-gene–OG mapping | Included as `inputs/all_species.seurat_gene_og.tsv.gz`; links `gene_seurat` and `object_species` to `orthogroup`. |
| Broccoli OG membership | Included at `../../gene_annotation/orthogroup_inference/reference_results/Anthozoa/orthogroup_conservation.csv`. |
| OG functional annotations | Included at `../../gene_annotation/functional_annotation/reference_results/orthogroups/orthogroup_final_annotations.csv`. |

Table paths and analysis parameters are set in `config/analysis_config.yaml`.

## Running the scripts

Run from this module directory. To rebuild the OG_SUM objects from RNA inputs and run all five positive comparisons:

```bash
bash run_all.sh
```

Execution order: `run_shared.sh` runs steps 00–02 (input checks, mapping summary and OG_SUM construction); `run_comparison.sh` runs steps 03–06 for each comparison (markers, checks, conservation and heatmap). To use the supplied objects directly, run `bash run_comparison.sh calicoblast`, replacing `calicoblast` with the required comparison ID.

After step 05, add final annotations to each positive summary, replacing `calicoblast` with the required comparison ID:

```bash
Rscript scripts/M02_add_final_annotation_to_conservation_summary.R \
  results/comparisons/calicoblast/conservation/calicoblast.OG_SUM.conservation_summary.tsv \
  --overwrite
```

The two independent branches are:

| Command | Prerequisite and execution |
|---|---|
| `bash run_log2fc_all_celltypes.sh 2` | Requires step 02 objects; M03 calculates statistics for each observed non-unknown/undefined cell type, then M04 checks all species. |
| `bash run_downregulated_og_postprocess.sh` | Requires step 03 full marker tables and step 02 feature metadata; M05 summarizes downregulated conservation → M02 adds annotations → M06 checks outputs. Runs all five configured comparisons. |

Optional candidate inspection: M07 extracts a species' specific and shared markers from an annotated positive summary. M08 adds that species' filtered marker statistics and ranks the species-specific subset. The supplied example is Hcoe in the 11-species `calicoblast` comparison.

## Final output

| Location | Contents |
|---|---|
| `intermediate/02_og_sum_objects/` | Eleven adopted OG_SUM objects and matching feature metadata; step 02 rebuilds these files from RNA inputs. |
| `results/comparisons/<comparison>/conservation/` | Five `*.conservation_summary.with_final_annotation.tsv` tables: OG, species support, conservation class and annotation. Steps 05–06 also produce species details and broad-candidate heatmaps. |
| `results/comparisons/<comparison>/markers/full/` | Adopted full marker tables for calicoblast, alga-hosting gastrodermis and cnidocyte, reused by the TF analysis. |
| `results/comparisons/<comparison>/downregulated/conservation/` | Four supplied annotated downregulated conservation summaries. |
| `results/LOG2FC/by_species/<species>/` | Eleven `*.all_celltypes_OG_log2FC.tsv.gz` long tables: OG × cell-type expression statistics. Eleven `*.all_celltypes_OG_log2FC_matrix.tsv.gz` matrices: OG rows × cell-type columns. |
