# Calicoblast multi-OG family trees

## Workflow

Frozen family members → MAFFT alignment → ClipKIT trimming → IQ-TREE inference → GeneRax reconciliation → expression/domain tracks → display selection and plotting.

The CA, FERLIN, MEP1B/astacin and SLC4 trees contain 376 tips in total: 286 members from the 11 atlas species and 90 Hvul/Aaur references. Original Broccoli OG assignments are retained. GeneRax EVAL infers rooting and duplication/loss events while retaining the IQ-TREE unrooted topology, as described in *Gene-family phylogenies and conserved-domain annotation* in the manuscript Methods.

| Directory under `full_tree/` | Role |
|---|---|
| `scripts/` | Complete-tree inference, reconciliation and expression/domain integration |
| `display_main_138/` | Main display: 104 atlas members and 34 representative references |
| `display_control_194/` | Comparison display: the same 104 atlas members and all 90 references |

The comparison selects markers and the highest observed log2FC member per family, OG and species, with a lexical fallback when expression is unavailable. The main display reuses that selection and chooses reference representatives within pure-reference clades of the complete trees. Both displays prune existing trees without repeating inference.

Each display generates only the four `article_final` PDFs, using the original article layout and parameters.

## Required files

Inputs are supplied under `full_tree/inputs/`:

- `frozen_members/`: protein sequences, membership tables and adopted reference-selection evidence. References meet the dual 70% coverage and family-specific CDD criteria recorded in `config/reference_selection_rules.yaml`.
- `species_tree/`: the named Anthozoa22 backbone, pruned to the atlas species and extended with Hvul/Aaur by the input-preparation script.
- `merged_tracks/`: adopted gene-level `avg_log2FC`, detection fractions and marker statistics; CDD domain hits and per-protein domain status. These tables are matched to tree members without rerunning differential expression or domain searches. CDD annotations derive from RPS-BLAST against CDD v3.21 on unaligned full-length proteins, followed by rpsbproc processing. Missing expression remains `NA`.

Adopted alignments, complete trees, reconciliation tables and display selections are also supplied for plotting. Earlier reference searches, longest-isoform selection and original species-tree inference precede this workflow.

Use Python with PyYAML/Biopython, MAFFT, ClipKIT, IQ-TREE, GeneRax/MPI, Ghostscript, and R with ape, ggtree, ggplot2, ggnewscale, scales, digest and yaml. Set executable paths in `full_tree/config/analysis_config.yaml` and the Rscript path in `full_tree/display_control_194/config/analysis_config.yaml`; both display launchers use the latter. The existing GeneRax 1.2.2/MPI binary-hash requirements remain in the inference configuration.

## Running the scripts

From this module directory, run the complete workflow in a working copy:

```bash
cd full_tree
python3 scripts/run_all.py
(cd display_main_138 && python3 scripts/run_display.py)
```

`run_all.py` prepares inputs, runs inference and reconciliation, assembles tracks, plots the complete trees, and generates/checks the 194-tip comparison. Run the 138-tip main display afterwards. The original stage and resume interfaces remain available.

To redraw both displays using the supplied complete-tree results, start in `full_tree/`:

```bash
(cd display_control_194 && python3 scripts/run_pipeline.py)
(cd display_main_138 && python3 scripts/run_display.py)
```

Keep this order: the main display uses the comparison's selection and domain-color baseline. Its existing snapshot check requires those parent results to remain unchanged during the main-display run.

## Final output

- `full_tree/results/trees/` and `results/03_reconciliation/`: complete trees, nodes and duplication/loss event tables.
- `full_tree/results/05_tip_tracks/`: expression/domain tracks and coverage summaries.
- Each display's `results/01_tip_selection/`, `results/02_pruned_trees/` and `results/05_article_figures/`: selected members, pruned trees and final PDFs.

The article output directory also contains the figure manifest (including pruned-tree paths and hashes), tip labels and plotting records used by the validation scripts.

The supplied final PDFs are the adopted main and comparison displays. UFBoot describes gene-tree split support, not confidence in rooting or reconciliation events; loss events remain in tables but are not drawn.
