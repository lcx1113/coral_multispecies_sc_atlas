# Cross-species SATURN — one script from training to figures

`saturn_pipeline.py` runs the whole cross-species SATURN workflow and draws the
two overview figures:

1. build the run table (one row per species: h5ad + ESM2 protein embedding),
2. train the model with SATURN's `train-saturn.py`,
3. recompute an integrated UMAP on the macrogenes and cache the coordinates,
4. plot


## Before you run

Edit the **DATA PATHS** block at the top of `saturn_pipeline.py` (or pass the
matching flags). Placeholders:

```python
DEFAULT_SATURN_ROOT    = "/path/to/SATURN"                  # SATURN repository
DEFAULT_H5AD_DIR       = "/path/to/saturn_inputs/final_h5ad" # one <species>.h5ad per species
DEFAULT_EMBEDDING_DIR  = "/path/to/protein_embedding/cross_species"
```

Each input `<species>.h5ad` must carry a cell-type column and, optionally, a
`color` column giving that label's colour:

| Flag | Default | Column used for |
|------|---------|-----------------|
| `--label-col` | `cell_type` | the annotation passed to SATURN (`--in_label_col` / `--ref_label_col`) **and** the per-species labels used to validate the output and to derive the shared colour palette |
| `--color-col` | `color` | optional per-cell colour of that label; missing colours are filled from a fallback palette |

SATURN writes the integrated annotation back as `labels2` (the unprefixed
version of `labels`); that column name is fixed by SATURN and is what all
downstream steps read. The embedding path is built from `--embedding-template`.

## Usage

```bash
# full run: run table -> train -> UMAP -> figures
python saturn_pipeline.py \
    --saturn-root /path/to/SATURN \
    --h5ad-dir /path/to/final_h5ad \
    --embedding-dir /path/to/protein_embedding/cross_species \
    --work-dir run \
    --num-macrogenes 1000 --hv-genes 4000 \
    --title "11-species SATURN"

# skip training, use an existing SATURN output
python saturn_pipeline.py --saturn-h5ad path/to/out_saturn_seed_0.h5ad \
    --h5ad-dir /path/to/final_h5ad --work-dir run

# re-draw the figures from the cached coordinates (seconds)
python saturn_pipeline.py --work-dir run --reuse-coords

# see the training command without running it
python saturn_pipeline.py ... --dry-run
```

`--work-dir` receives `cross_species_run.csv`, `saturn_results/`,
`umap_coordinates.tsv.gz`, `celltype_palette.tsv`, both PDFs, `pipeline.log`
and `_SUCCESS.json`. `saturn_results/` must exist before training; the script
creates it.

## Conventions

- **Coordinates.** The UMAP is recomputed from `.X` (the macrogenes SATURN
  actually used), PCA → kNN → UMAP with a fixed seed, so the layout is
  reproducible. `--use-macrogenes` embeds on `obsm["macrogenes"]` instead (the
  full centroid matrix, which can have more columns than were used).
- **Colours.** One colour per cell type, pooled from the per-species `color`
  columns so a cell type has the same colour everywhere; conflicts (the same
  cell type with different colours in different species) resolve to the species
  contributing the most cells, and are logged.
- **Unannotated clusters.** Labels starting with `unknown` (`unknown1`, ...)
  are drawn as one grey `unknown` category so they do not dominate the palette.
  Display-only: `umap_coordinates.tsv.gz` keeps the original labels. Pass
  `--keep-unknown-split` to keep the per-category colours.
- Figures are written as PDF only.


## Requirements

```
python >= 3.9, SATURN (checked out), a CUDA device for training
numpy, pandas, matplotlib, scanpy, anndata, h5py, torch, pytorch-lightning
```
