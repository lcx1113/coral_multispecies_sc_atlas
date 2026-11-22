# MITF-family E-box analysis (Figures 1 and 2)

## 1. Overview

This module tests whether genes that are significantly enriched in the symbiont-hosting
gastrodermis (`gastrodermis_alga_hosting`) carry more predicted MITF/TFE-family E-box
elements in their promoters than the other genes in the same cross-species single-cell
marker table. It reproduces **Figure 1** (element frequency across the 16 annotated cell
types) and **Figure 2** (odds ratios across eleven PWMs).

In the published run, 150 of 413 markers (36.3%) carry at least one element versus 26.5% of
the 6,971 background genes (odds ratio 1.585, Woolf 95% CI 1.29–1.95; two-sided Fisher
*P* = 2.3 × 10⁻⁵; Benjamini–Hochberg *q* = 3.6 × 10⁻⁴). Only the five MITF/TFE-family
matrices are enriched; five unrelated matrices and a scrambled control are not.

These are **predicted elements, not binding sites**: the `CACGTG` core is shared with other
E-box-binding factors, and no chromatin-accessibility or binding data are included.

## 2. Requirements

Python 3.8.10 with NumPy 1.23.5, pandas 2.0.3, SciPy 1.10.1 and Matplotlib 3.7.5.

```bash
conda env create -f environment.yml
conda activate mitf-ebox
# or: pip install -r requirements.txt
```

No external motif-scanning software is required; the scanner is implemented in
`code/02_scan_motifs.py`.

## 3. Input files

Place the three files in `data/` (or point to them with the environment variables listed in
§4). The default names are resolved in `code/mitf_config.py`.

| File | Required content and format |
|---|---|
| `xenia.contigs.FINAL.fasta` | Genome FASTA. Lowercase bases are treated as soft-masked repeats by the scanner (scored normally by default). |
| `Xesp_long.annot.gtf` | GTF with `gene_id` and `transcript_id` attributes on `transcript` and `exon` records. The representative transcript is the one with the longest summed exon length. |
| `seu.Xesp.markers.cts.csv` | **Tab-separated** long table (despite the `.csv` suffix), one row per (gene, cell type). Required columns: `gene_seurat`, `focus_node`, `p_val_adj`, `avg_log2FC`, `pct.1`, `pct.2`, `umifrac`. |

`gene_seurat` uses the form `Xe-014349`; it is converted to the GTF form `Xe_014349` by
replacing `-` with `_`. Each `focus_node` block is a ranked gene list, **not** a
significant-marker list, so markers are recomputed explicitly (§5). The analysis drops the
ten unassigned `unknown*` clusters and retains 16 annotated cell types.

## 4. Usage

Run the whole pipeline:

```bash
bash code/run_all.sh
```

Or step by step:

```bash
python3 code/00_make_scrambled_motif.py   # scrambled control matrix (seed 20261008)
python3 code/01_build_promoters.py        # promoters.bed/.fa, gene_tss.tsv
bash    code/03_run_scans.sh              # 11 matrices x prox2kb window
python3 code/04_enrichment.py             # source tables for Figures 1 and 2
python3 code/05_figures.py                # Figure 1 and Figure 2 PDFs
```

Defaults can be overridden without editing code:

```bash
MITF_ROOT=/path/to/output \
MITF_MARKERS=/path/to/seu.Xesp.markers.cts.csv \
MITF_GENOME=/path/to/xenia.contigs.FINAL.fasta \
MITF_GTF=/path/to/Xesp_long.annot.gtf \
bash code/run_all.sh
```

Other variables: `MITF_DATA`, `MITF_MOTIFS`, `MITF_EXCLUDE_PREFIX`, `MITF_SCAN_DIR`,
`PYTHON`, `JOBS`. On 8 CPU cores the pipeline takes about 2 minutes.

## 5. Key parameters

| Parameter | Value used | Code |
|---|---|---|
| Promoter window | TSS −2,000/+500 bp, in gene orientation (reverse-complemented on the minus strand) | `01`, `WINDOWS` |
| Promoter length | Extraction drops windows that retain <50% of the target length (`MIN_FRAC_KEPT = 0.5`); gene-level enrichment then keeps only full-length windows of ≥2,000 bp (`MIN_WIN_LEN`) | `01`, `04` |
| Marker threshold | `p_val_adj < 0.05` **and** `avg_log2FC > 0.25`, applied identically to all 16 cell types | `config.PADJ_MAX`, `config.LOG2FC_MIN` |
| Element threshold | Per-position `p < 1e-4` (strictly less than; `allp < args.pval`) | `02` |
| PWM scoring | Pseudocount 0.5; log-odds against a uniform background; palindromic matrices scanned on one strand; background score distribution computed by discretised convolution of per-position distributions | `02` |
| Background set | The other genes in the same filtered marker table (not the whole genome), restricted to full-length scanned promoters | `04` |
| Statistical test | Two-sided Fisher exact test on (marker vs background) × (≥1 element vs none); odds ratio with Woolf 95% CI; Benjamini–Hochberg across the 16 cell types. Figure 2 repeats this for eleven matrices | `04` |

## 6. Outputs

| File | Content |
|---|---|
| `results/tables/cluster_mitf_enrichment.tsv` | Figure 1 source data; one row per cell type |
| `results/tables/motif_robustness.tsv` | Figure 2 source data; one row per matrix |
| `results/figures/fig1_cluster_enrichment.pdf` | Figure 1 |
| `results/figures/fig2_motif_robustness.pdf` | Figure 2 |
| `results/promoters.bed`, `promoters.fa`, `gene_tss.tsv` | Intermediate promoter windows and TSS table |
| `results/scans/<matrix>.prox2kb.{gene,hits}.tsv` | Intermediate per-gene and per-element scan results |

Column definitions for the two source tables are in [`results/README.md`](results/README.md).
Re-running the pipeline reproduces the two tables; the figure PDFs differ from the deposited
copies only in the PDF creation timestamp.
