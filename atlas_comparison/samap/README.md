# SAMap

## Workflow

The pairwise SAMap comparison is used to compare our single-cell data with the annotated cell types in the published Spis single-cell atlas ([Levy et al. (2025)](https://doi.org/10.1038/s41586-025-09623-6)), providing reference evidence for cell-type annotation.

Pairwise protein searches prepare homology maps. The supplied Ahya–Spis notebook and the eleven-species Python script are separate analysis examples using SAMap/SAM. Their original cell/block order, calls and parameters are retained.

## Required files

Python with SAMap, SAM, Scanpy/AnnData and their plotting dependencies; BLAST+ for homology searches. Prepare writable longest-protein FASTAs for the mapping scripts, the Ahya/Spis h5ad objects expected by the notebook, or the eleven `h5ad/dat.<species>.seurat_final.h5ad` files named by `all_species_run_samap.py`. Preserve the existing annotation columns and matrix contents expected by each input.

Replace `/absolute/path/...` in the chosen entry and mapping batch. Point `f_maps` to the prepared SAMap map directory.

## Running the scripts

For homology maps, work in `run_maps_for_SAMap/` and run `bash extended_longest_blast_all_twovstwo.sh`; it calls `extended_longest_map_genes.sh` for the species pairs. Inputs must match `*_cleaned_modified.longest.fasta`; the batch derives two-letter lowercase species IDs from filenames (ad/ah/am/fs/gf/hc/pc/pd/pl/sp/xe). Longest-isoform selection and identifier preparation are external steps. The original two-direction searches retain 24 threads and E-value 1e-6. BLAST indexes are written beside the FASTAs, so use writable input copies.

For the Ahya–Spis example, open `species_celltype_annotation_samap.ipynb` with this module as the notebook working directory and execute its cells in their original order. For the separate multi-species branch, run `python all_species_run_samap.py` from this module after preparing its `h5ad/` inputs. Both retain their successive `sm.run(...)` calls; neither is replaced by a new wrapper.

## Final output

The mapping batch writes `extend_longest_maps/`. Analysis outputs follow the original relative paths in the notebook/script: processed h5ad data, serialized SAMap objects, mapping scores and figures.

