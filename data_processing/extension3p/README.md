# 3′ annotation extension

## Workflow

| Step | Script | Purpose |
|---|---|---|
| 1 | [peak3p.sh](peak3p.sh) | Align cleaned R2 reads with STAR, split alignments by strand, and call peaks with MACS2. |
| 2 | [extend3p.R](extend3p.R) | Extend annotation intervals using the peaks and add orphan features. |
| 3 | [merge_gtf_with_orphan.py](merge_gtf_with_orphan.py) | Merge the extended intervals with original exons to produce the final GTF. |

## Required files

- Species-specific genome FASTA (plain or gzip-compressed).
- Original GTF annotation for the same genome; use the same file throughout the three steps.
- Cleaned R2 FASTQ (gzip-compressed), after droptag, poly(A)/poly(G) trimming, gap removal and short-read filtering.

Required software: Bash with GNU command-line utilities, STAR, sambamba, MACS2, R with `optparse`, `GenomicRanges` and `rtracklayer`, and Python 3.

## Running the scripts

Gfas is shown below as an example. The other species processed for 3′ annotation extension used the same execution procedure and workflow, with their respective input files and species prefixes.

Start in the repository root and replace the input and work-directory placeholders. Relative output paths are resolved from the work directory. Run each step separately and inspect its outputs before continuing.

```bash
# Start in the repository root. Replace all /absolute/path placeholders.
CODE_DIR="$(pwd)/data_processing/extension3p"
GENOME=/absolute/path/genome.fa
ORIGINAL_GTF=/absolute/path/original.gtf
CLEAN_R2=/absolute/path/clean_R2.fastq.gz
WORK_DIR=/absolute/path/extension_work_Gfas
PREFIX=Gfas

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"
```

**1. Align reads and call strand-specific peaks.**

```bash
bash "$CODE_DIR/peak3p.sh" \
  "$PREFIX" "$GENOME" "$ORIGINAL_GTF" "$CLEAN_R2" annotate_3p/ 24
```

**2. Extend annotation intervals using the combined peaks.**

```bash
Rscript "$CODE_DIR/extend3p.R" \
  -g "$ORIGINAL_GTF" \
  -p "annotate_3p/pool_${PREFIX}_peaks_3p.broadPeak" \
  -o "annotate_3p/reannotate_${PREFIX}_genes.gtf" \
  -m 5000 -a -r "$PREFIX" -q 1e-3
```

**3. Merge the extended intervals with the original exons.**

```bash
python3 "$CODE_DIR/merge_gtf_with_orphan.py" \
  --extend_gtf "annotate_3p/reannotate_${PREFIX}_genes.gtf" \
  --orig_gtf "$ORIGINAL_GTF" \
  --output "${PREFIX}_extended_merged.gtf"
```

## Final output

`Gfas_extended_merged.gtf` is written to the work directory. It contains transcript and exon records, including the extended regions and orphan features. The final GTF was used as the annotation input for downstream dropEst analysis. Intermediate peak-calling and extension files are written under `annotate_3p/`.

## Attribution

The 3′-extension workflow draws on [Levy et al. (2025)](https://doi.org/10.1038/s41586-025-09623-6) and the [oculina-coral-sc-atlas code](https://github.com/sebepedroslab/oculina-coral-sc-atlas/) ([original repository](https://github.com/xgrau/oculina-coral-sc-atlas)). The exon-merging step is provided by the project's local `merge_gtf_with_orphan.py` script.
