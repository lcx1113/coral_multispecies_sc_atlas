# Data and code availability

Standalone statement for the Cell STAR Methods section **RESOURCE AVAILABILITY → Data and
code availability**. Replace every bracketed placeholder before submission.

---

## Statement

- The *Xenia* genome assembly (`xenia.contigs.FINAL.fasta`) and the long-read annotation
  (`Xesp_long.annot.gtf`) are available at [REPOSITORY / ACCESSION].
- The single-cell RNA-seq data are available at [GEO ACCESSION]. The cell-type marker table
  analysed here (`seu.Xesp.markers.cts.csv`, 2026-10-26 version) is provided as
  [DATA S1 / GEO SUPPLEMENTARY FILE].
- All original analysis code, including the code that generates Figures 1 and 2, is
  available at [GITHUB URL] and has been deposited at [ZENODO DOI]. DOIs are listed in the
  key resources table.
- Source data for Figures 1 and 2 are provided in `results/tables/` and as [DATA S1].
- Any additional information required to reanalyse the data reported in this paper is
  available from the lead contact, [NAME] ([EMAIL]), upon request.

---

## What is deposited where

| Item | Repository | Identifier |
|---|---|---|
| Genome assembly | [REPOSITORY] | [ACCESSION] |
| Genome annotation | [REPOSITORY] | [ACCESSION] |
| scRNA-seq atlas | GEO | [GSE ACCESSION] |
| Cell-type marker table | GEO / supplementary | [DATA S1] |
| Analysis code (this repository) | GitHub | [GITHUB URL] |
| Archived code release | Zenodo | [DOI] |

## Code archive

The repository contains eleven PWM matrices, the analysis and figure scripts, the two source
tables and the vector PDFs for Figures 1 and 2. To create the Zenodo archive:

1. Push the repository to GitHub at [GITHUB URL].
2. In GitHub, enable the Zenodo integration for the repository
   (<https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content>).
3. Create a release (e.g. `v1.0.0`); Zenodo mints a DOI for that release.
4. Insert the DOI above and in `KEY_RESOURCES_TABLE.csv`, and add it to the manuscript.

## Before submission — checklist

- [ ] Genome assembly accession and version.
- [ ] GEO accession for the single-cell data.
- [ ] Supplementary file number for the marker table.
- [ ] Public GitHub URL.
- [ ] Zenodo DOI for the tagged release.
- [ ] Lead-contact name and e-mail.
- [ ] Author-approved license (the repository ships with an MIT `LICENSE` placeholder).
- [ ] Confirm that the marker table can be redistributed; if not, point to GEO only.
