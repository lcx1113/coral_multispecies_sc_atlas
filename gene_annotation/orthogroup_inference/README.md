# Orthogroup inference and TF integration

## Workflow

Broccoli orthology was inferred separately from protein sequences for two species sets: **Anthozoa (22 species)** and **Metazoa (34 species)**, using commands of the form shown in `run_broccoli.sh`. GTF annotations map protein/transcript identifiers to genes during result conversion.

Anthozoa orthogroups and gene pairs support expression comparisons. Metazoa orthogroups are additionally integrated with TF-family results for reference annotation. TF inference follows HMM search → DIAMOND → MCL → MAFFT → ClipKIT → alignment filtering → IQ-TREE → manual tree selection → Possvm → Metazoa integration.

## Required files

The original selected protein FASTAs and matching GTFs are included:

| Scope | Protein sequences | GTF annotations |
|---|---|---|
| Anthozoa | `inputs/Anthozoa/proteins/` — 22 files | `inputs/Anthozoa/gtf/` — 22 files |
| Metazoa | `inputs/Metazoa/proteins/` — 34 files | `inputs/Metazoa/gtf/` — 33 files |

[Species-to-file mapping](inputs/species_files.tsv) lists each input, with paths relative to this module. Metazoa Hhon has no supplied GTF; the original conversion uses transcript identifiers for species without GTFs. Original filenames and sequence/annotation contents are retained.

For reconversion, provide Broccoli step-3 `orthologous_groups.txt` and, for Anthozoa pairs, step-4 `orthologous_pairs.txt`. These raw outputs are external; the converted reference tables are supplied below. Set their paths in the relevant Rmd. TF inference additionally needs the combined Metazoa protein FASTA at `tf_prediction/Metazoa/data/seq_Metazoa_20260425.fasta` and the upstream HMM profiles; these are not bundled. Small TF parameter and reference-name/outgroup tables are included.

Dependencies: Broccoli and its tools; R with tidyverse, rtracklayer and the packages requested by `functions/helper.R`; HMMER, DIAMOND, MCL, MAFFT, ClipKIT, IQ-TREE, Python/Biopython and Possvm. Set tool/HMM paths and environment activation for your installation.

## Running the scripts

**Broccoli.** Use `run_broccoli.sh` as the LSF submission example for each species set. Set the Broccoli executable and a separate working directory containing that set's protein FASTAs. Keep the two sets separate; GTFs are used only in the subsequent conversion.

**Conversion.** Start each Rmd in a fresh R session from this module directory and execute its code chunks in order. First create the output directories:

```bash
mkdir -p work/broccoli/Anthozoa work/broccoli/Metazoa
```

| Branch | Entry | Generated output under `work/broccoli/` |
|---|---|---|
| Anthozoa orthogroups | `scripts/Anthozoa/transfer_broccoli_Anthozoa.Rmd` | `Anthozoa/orthogroup_conservation.csv` |
| Anthozoa gene pairs (independent branch) | `scripts/Anthozoa/transfer_orthogroup_pairs_Anthozoa.Rmd` | `Anthozoa/orthogroup_pairs.genes.csv` |
| Metazoa orthogroups | `scripts/Metazoa/transfer_broccoli_Metazoa.Rmd` | `Metazoa/orthogroup_conservation.csv` |
| Metazoa TF integration | `scripts/Metazoa/integrating_TF.Rmd`, after Possvm and Metazoa conversion | `Metazoa/orthogroup_conservation_integrated.csv` |

To integrate TFs using the supplied Metazoa orthogroups, copy `reference_results/Metazoa/orthogroup_conservation.csv` to `work/broccoli/Metazoa/` instead of repeating Broccoli/conversion.

**TF inference.** From `tf_prediction/Metazoa`:

```bash
cd tf_prediction/Metazoa
bash hmmer_py.sh
bash diamond.sh
bash mcl_clustering.sh
python mafft.py
python clipkit.py
python filter_clipkit.py
python iq_tree_fixed.py
# Assemble the selected family/*.treefile files in iq_tree_final/ before Possvm.
bash possvm.sh
```

The supplied `tf_parameter.csv` retains HLH inflation 1.4 and Homeodomains 1.6. The adopted tree set comprised 676 other trees, 27 HLH trees and 50 Homeodomain trees; tree files remain external. The recorded zf-C2H2 cluster 332 failure is retained.

The supplied Metazoa reference tables contain all 34 species.

## Final output

The adopted, already converted results are included for downstream use:

| Supplied file | Use |
|---|---|
| [Anthozoa orthogroups](reference_results/Anthozoa/orthogroup_conservation.csv) | Original Broccoli membership for expression comparisons |
| [Anthozoa orthopairs](reference_results/Anthozoa/orthogroup_pairs.genes.csv) | Gene-level direct pairs; not a prefiltered one-to-one set |
| [Metazoa orthogroups](reference_results/Metazoa/orthogroup_conservation.csv) | Original Broccoli membership before TF integration |
| [Metazoa TF-integrated orthogroups](reference_results/Metazoa/orthogroup_conservation_integrated.csv) | Adopted membership for reference annotation |

These `.csv` files are tab-delimited. Anthozoa, original Metazoa and TF-integrated Metazoa use independent OG numbering; matching numbers do not identify matching groups. TF integration replaces memberships and renumbers groups. New conversion outputs go to `work/broccoli/`, separately from these reference copies.

## Attribution

The TF domain-guided method and HMM resources follow Levy et al. and [oculina-coral-sc-atlas](https://github.com/sebepedroslab/oculina-coral-sc-atlas/). [Possvm](https://github.com/xgrau/possvm-orthology) is by Grau-Bové and Sebé-Pedrós ([paper](https://doi.org/10.1093/molbev/msab234)); its supplied code and [GPLv3 license](toolkit/possvm-orthology/LICENSE.md) are retained. That license does not imply a license for every other file in this module.
