# Read processing and quantification

## Workflow

[dropest_pipeline_change_STAR.sh](dropest_pipeline_change_STAR.sh) runs the following steps in sequence:

`droptag → cutadapt (polyA) → cutadapt (polyG) → SeqKit → STAR → dropEst → dropReport`

The pipeline starts with paired R1/R2 reads and a prepared counting reference. The final GTF from [3′ annotation extension](../extension3p/README.md) is used for dropEst; prepare the corresponding STAR index before running this module.

## Required files

- Paired R1 and R2 FASTQ files, gzip-compressed.
- Final counting GTF and a corresponding, existing STAR index.
- Bundled [split_seq.xml](config/split_seq.xml) and [barcode123_rc.txt](config/barcode123_rc.txt).

Required software: Bash and GNU command-line utilities, conda, droptag/dropEst/dropReport.Rsc, cutadapt, SeqKit, STAR, and the R dependencies required by dropReport.

## Running the scripts

Before the first run, replace the `/absolute/path` placeholders:

- In the script, set `DROP_BIN`, the conda initialization and environment paths, and the STAR executable path.
- In `config/split_seq.xml`, set `barcodes_file` to the absolute path of the bundled `config/barcode123_rc.txt`. The script locates the XML relative to its own location.

Use absolute input paths and paths without whitespace. Start in the repository root and choose a sample-specific work directory. The following example uses Gfas with 12 threads.

```bash
CODE_DIR="$(pwd)/data_processing/read_processing_quantification"
FASTQ_R1=/absolute/path/Gfas_R1.fastq.gz
FASTQ_R2=/absolute/path/Gfas_R2.fastq.gz
STAR_INDEX=/absolute/path/STAR_index_Gfas
FINAL_GTF=/absolute/path/Gfas_extended_merged.gtf
WORK_DIR=/absolute/path/counts_work_Gfas

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"
bash "$CODE_DIR/dropest_pipeline_change_STAR.sh" \
  "$FASTQ_R1" "$FASTQ_R2" "$STAR_INDEX" "$FINAL_GTF" 12 Gfas
```

The six arguments are R1, R2, STAR index, GTF, threads (default 12), and prefix (default derived from the R2 filename). Outputs are written under the work directory. The pipeline removes the `no_A_*` and `no_AG_*` intermediate FASTQ files after STAR alignment and retains the final `clean_no_AG_*` FASTQ.

## Final output

The main output is `est_out/cell.counts.rds`. Its `cm` matrix is used for downstream QC; columns represent barcodes before cell filtering, ambient-RNA correction and normalization.

```r
exp <- readRDS("/absolute/path/counts_work_Gfas/est_out/cell.counts.rds")
counts <- exp$cm
```

Additional outputs include tagged/cleaned reads and droptag parameters in the work directory, STAR BAM and logs in `alignment/`, and dropEst/dropReport outputs in `est_out/`.
