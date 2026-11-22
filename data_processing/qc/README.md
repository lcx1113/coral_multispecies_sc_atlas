# Single-cell QC

## Workflow

[run_full_QC_standalone.R](stand_alone_script/run_full_QC_standalone.R) runs:

`All-droplet counts → emptyDrops → SoupX → marker-poor cluster removal → top UMI/mitochondrial filtering → DoubletFinder → RNA reclustering`

The output is a QC-filtered Seurat object for downstream analysis. Final atlas annotation is a separate step.

## Required files

- One RDS per sample containing an unfiltered gene-by-barcode count matrix, or a dropEst object with its `cm` matrix. Include all droplets; do not substitute filtered or corrected counts. See [read processing and quantification](../read_processing_quantification/README.md).
- [config/sample.tsv](config/sample.tsv): the default eight-sample run table (Ahya, Adig, Amur, Fusp, Gfas, Hcoe, Pdam and Plut).
- [inputs/samples.tsv](inputs/samples.tsv): the nine-sample input list, including Pcru for the separate invocation below. Both tables use `sample_id`, `input_rds`, `mito_gene_file` and `apply_mito_filter`.
- Seven mitochondrial gene lists in [inputs/mito](inputs/mito). Hcoe and Pcru have mitochondrial filtering disabled. Fusp uses `ffun_mito_gene.txt`.

Settings are in [default_config.R](stand_alone_script/default_config.R); [sample_processing_parameters.tsv](sample_processing_parameters.tsv) summarizes adopted sample settings and is not loaded automatically. Optional sample-table overrides, such as `emptydrops_by_rank` and `marker_filter_cutoff`, take precedence over the shared configuration.

Required software: R with Seurat, SeuratObject, sctransform, glmGamPoi, DropletUtils, SingleCellExperiment, SoupX, DoubletFinder, Matrix, ggplot2, patchwork, dplyr, tidyr, tibble and igraph. Optional `pdfunite` merges report PDFs; individual PDFs and an HTML index are retained without it.

## Running the scripts

Replace `input_rds` placeholders in the table used for your run with absolute paths to the corresponding count files. Run from this QC module directory: mitochondrial paths in both tables are relative to this working directory. Choose a separate output directory for each analysis setup.

From the repository root, run one sample:

```bash
cd data_processing/qc
Rscript stand_alone_script/run_full_QC_standalone.R \
  --sample_id Ahya \
  --stop_after final --rerun_from none \
  --output_root /absolute/path/qc_Ahya
```

Alternatively, from the same QC directory, omit `--sample_id` to process all eight default samples sequentially:

```bash
Rscript stand_alone_script/run_full_QC_standalone.R \
  --stop_after final --rerun_from none \
  --output_root /absolute/path/qc_eight_species
```

For **Pcru**, prepare a working copy of the shared configuration:

```bash
mkdir -p /absolute/path/Pcru_work
cp stand_alone_script/default_config.R /absolute/path/Pcru_work/run_config.R
```

Edit these three entries in the copied `qc_config` list:

| Entry | Pcru value | Eight-sample value |
|---|---|---|
| `emptydrops$niters` | `10000L` | `100000L` |
| `marker_filter$random_seed` | `100L` | `123L` |
| `soupx$force_accept` | `TRUE` | `"default"` (omits `forceAccept`) |

Keep the workflow seed at 100, marker cutoff at 15 and Pcru mitochondrial filtering at `FALSE`. Fill the Pcru input path in `inputs/samples.tsv`, then run only that row:

```bash
Rscript stand_alone_script/run_full_QC_standalone.R \
  --sample_tsv inputs/samples.tsv --sample_id Pcru \
  --config /absolute/path/Pcru_work/run_config.R \
  --stop_after final --rerun_from none \
  --output_root /absolute/path/Pcru_work/output
```

Do not run the entire nine-row table with the eight-sample configuration.

Existing stage outputs are reused by default without checking whether inputs or parameters changed. `--rerun_from STEP` deletes that stage and downstream outputs; `final_recluster` reuses existing DoubletFinder singlets. Keep historical results outside the selected output directory. Global summaries include existing sample directories under the same output root.

`--stop_after` accepts `emptydrops`, `soupx`, `marker_filter`, `top_mito`, `doubletfinder` or `final`. Stopping at `doubletfinder` completes preparation and the pK scan only, before doublet classification/removal. `QC_STOP_AFTER` and `QC_RERUN_FROM` apply when their command-line options are omitted. There is no `--check_only` mode; unknown options are not reliably rejected. Without `--output_root`, outputs go to `stand_alone_script/output`.

## Final output

The main output is `<output_root>/<sample_id>/90_intermediate_objects/final_filtered_reclustered_object.rds`, containing SoupX-corrected, filtered RNA counts and final RNA clustering.

Each sample also has stage-specific barcode lists, filtering decisions, marker tables and plots, a `00_run_info/QC_step_cell_count_summary.tsv`, and reports in `98_report`. Cross-sample summaries and final object paths are written to `<output_root>/99_summary/`, including `final_object_paths.tsv`.
