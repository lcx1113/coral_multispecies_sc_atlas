# Gene and orthogroup annotation

## Workflow

BLASTP against Swiss-Prot → first qualifying hit among six reference taxa → gene annotation combining GTF, Anthozoa orthology, Metazoa TF integration and UniProt → orthogroup member table → frequency summaries and final OG labels. The supplied reference tables preserve the adopted results.

## Required files

BLAST+, R with tidyverse/rtracklayer and the packages requested by `functions/helper.R`; query protein FASTAs, the Swiss-Prot BLAST database, matching GTFs, Anthozoa membership/gene pairs, TF-integrated Metazoa membership, and `Spis_Genes_ncbi_dataset.tsv`. Set `/absolute/path/...` locations in the relevant scripts. [Reference results](reference_results/) contain 12 UniProt hit tables (including Tcoc), 11 gene tables and the three adopted OG tables.

## Running the scripts

Start from this module directory. Example search using the existing interface:

```bash
mkdir -p work/gene_annotation work/orthogroups work/uniprot
bash scripts/query_uniprot_top_6species.sh -q /absolute/path/proteins/Adig.fasta -d /absolute/path/blastdb/swiss_prot -p Adig -o work/uniprot -t 1
```

The search writes `<prefix>_top_6species.tsv`; the adopted merged tables use `_top_6species_longest.tsv`. Its final hit-count message retains a working-directory lookup that can fail when output is elsewhere; the actual output path is `$TOPHIT`.

The Rmd files are interactive, ordered chunks. Start each in a fresh R session at the module directory:

1. `scripts/integrating_annotation_260428.Rmd` enters `work/gene_annotation`, reads the supplied reference UniProt tables and writes `annot.<species>.genes.csv`. To use a new search, explicitly point its `uniprot_fn` to that search result.
2. `scripts/annotation_orthogroup.Rmd` enters `work`, reads `work/gene_annotation/annot.<species>.genes.csv`, and writes `orthogroup_conservation_final.csv`. Before its subsequent read of `work/orthogroups/orthogroup_conservation_final.csv`, manually place the intended member table there. The original handoff is preserved.

You may instead copy the adopted gene tables into `work/gene_annotation/` to start at OG annotation. The gene-generation loop uses `Spin` while the OG loop expects `Spis`: resolve this explicit file-name handoff without changing gene IDs. Do not overwrite the reference tables with new calculations.

## Final output

The tab-delimited tables under `reference_results/` include gene annotations and the 11-species OG member table (`orthogroup_conservation_final.csv`). `orthogroup_status.csv` gives candidate names and their frequencies across all member rows of each OG, with each UniProt candidate paired with its own frequency. `orthogroup_final_annotations.csv` gives the selected primary and secondary labels and their support for 19,604 OGs.
