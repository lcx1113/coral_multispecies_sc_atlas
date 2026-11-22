# Calicoblast M3/M4 module scores

## Workflow

Compare calicoblast, epidermis and gastrodermis across 11 species using the 44-OG set (M3: 16; M4: 28) and the 266-OG reference (109/157).

Both sets follow the same workflow: SCT/data → background-corrected gene expression → within-OG weighted scores → equal mean across available OGs → M3/M4 plots. They use the same species-specific frozen per-feature control lists and within-OG weights, without resampling. Unmapped OGs are omitted from each module's denominator. Both plot raw scores with median/IQR summaries and retain two species layouts.

## Required files

| Input | Location |
|---|---|
| 44-OG definitions, annotations and mapping status | `narrative44_source/inputs/` |
| 266-OG definitions, annotations and mapping status | `reference266_source/inputs/` |
| Shared frozen control lists and weights, used by both scoring scripts | `narrative44_source/inputs/shared_controls/` |
| Plotting inputs: per-species cell scores, module manifest, summaries and coverage | `narrative44_source/narrative44/scores/OG_equal/` and `reference266_source/marker308/scores/OG_equal/` |
| Seurat objects, for recalculating either set | Download externally and set `object_path` in the shared `narrative44_source/inputs/species_manifest.tsv`. Objects must contain SCT/data and the supplied `cell_type` labels; Spin denotes Spis. |

Plotting requires R, data.table, ggplot2 and patchwork, with Cairo PDF support and the Liberation Sans font. Recalculating scores additionally requires Seurat, Matrix and digest.

## Running the scripts

Run from this module directory. To recalculate scores, set the Seurat object paths and run 01 before 02 for each set:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript narrative44_source/scripts/01_score_unified.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript narrative44_source/scripts/02_plot_unified_scores.R

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript reference266_source/scripts/01_score_unified.R
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript reference266_source/scripts/02_plot_updated_scores.R
```

To redraw the supplied scores without loading Seurat objects, run only the two plotting commands. Both scoring scripts accept `--species=Adig` for one species, or `--summarize-only` after all species have been scored for that set.

## Final output

Each scoring script writes scores and summaries to its set's `scores/OG_equal/` directory listed above. Plotting writes M3/M4 figures under `narrative44_source/narrative44/figures/` and `reference266_source/marker308/figures/`, including both layouts and individual panels. Generated figures and other saved analysis outputs are not bundled.

## Attribution

The scoring strategy draws on Seurat AddModuleScore, implemented here as the supplied custom OG-equal calculation. The original plotting code and attribution comments are retained.
