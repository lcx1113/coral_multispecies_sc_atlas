# Dollo reconstruction of OG membership and enrichment

## Workflow

Two branches are reconstructed separately: OG membership across species, and cell-type-enriched OG membership. Their Dollo presence/gain/loss matrices are combined for node summaries and annotated OG lists on the supplied 11-species tree.

The enrichment branches are `calicoblast`, `calicoblast_like`, `cindocyte` (cnidocyte) and `alga_hosting`. The original spellings and the `Spin` object / `Spis` gene-prefix distinction are retained. Enrichment inputs use gene-level FC > 1.1 and adjusted P < 0.001; they are not OG_SUM marker calls.

## Required files

Python with the dependencies imported by `membership/possvm_reconstruction.py`; R with dplyr, Seurat, ape, scales and the packages loaded by `functions/helper.R`.

For marker preparation, first run [the expression script](../../atlas_comparison/expression_icc_tree/scripts/scatlas_mutispecies_260818_reannotated.R). Dollo reads `seu.<species>.log2fc.cts.csv` and `dat.<species>.expression.cts_umicount.rds` directly from `atlas_comparison/expression_icc_tree/work/results_metacell_<species>_filt/` in the repository. It sums the UMI counts over the target cell types and calculates the fraction per 10,000. If the expression outputs are stored elsewhere, set `expression_work` in the preparation script. The R scripts link to supplied Anthozoa/Metazoa and gene-annotation tables under `gene_annotation/`. Selected extant tables, state matrices, node OG lists and PDFs are already included.

Marker tables and node OG lists use the primary and secondary labels from the project's final OG annotation table. For genes mapped to multiple OGs, marker preparation preserves missing input annotations.

## Running the scripts

Start from this module directory. To regenerate enrichment inputs:

```bash
Rscript enrichment/prepare_dollo_260827.R
```

Run the membership branch from `membership/Anthozoa_all_genome/`:

```bash
python3 ../possvm_reconstruction.py -tree Anthozoa_all_genome.newick \
  -ort ../../../../gene_annotation/functional_annotation/reference_results/orthogroups/orthogroup_conservation_final.csv \
  -out orthogroup_conservation.possvm
```

Run the enrichment branch from `enrichment/` in Bash:

```bash
for i in results_ct_evolution/anc.*extant.tsv; do
  python3 ../membership/possvm_reconstruction.py \
    -tree ../membership/Anthozoa_all_genome/Anthozoa_all_genome.newick \
    -ort <(cut -f1,2 "${i}") -out "${i%%.extant.tsv}.possvm"
done
```

Place the membership branch's three `orthogroup_conservation.possvm.dollo_{pres,gain,loss}.csv` matrices in `enrichment/results_ct_evolution/` (the adopted copies are supplied). Then, from the module directory:

```bash
Rscript enrichment/process_dollo_260827.R
```

## Final output

`membership/Anthozoa_all_genome/` contains membership states; `enrichment/results_ct_evolution/` contains enriched states, node OG lists and tree summaries. Enrichment absence is distinct from OG-member absence, and these inferred events are conditional on the supplied sampling and Dollo model.

## Attribution

The Dollo inference workflow and code were adapted from the [oculina-coral-sc-atlas repository](https://github.com/sebepedroslab/oculina-coral-sc-atlas/). Reconstruction uses the supplied local Possvm implementation. Original helper comments and citations are retained.
