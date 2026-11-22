#!/usr/bin/env python3
"""Shared configuration for the MITF promoter-motif analysis (Figures 1 and 2).

All paths are resolved relative to the repository so that a fresh clone runs without
editing code. They can be overridden with environment variables:

    MITF_ROOT=...          output root (default: <repo>/results)
    MITF_DATA=...          directory holding the three input files
    MITF_MARKERS=...       cell-type marker table
    MITF_GENOME=...        genome FASTA
    MITF_GTF=...           annotation GTF
    MITF_MOTIFS=...        JASPAR matrices

Directory layout expected by the pipeline::

    <repo>/
      code/            these scripts
      motifs/          JASPAR PWM files (*.jaspar)
      data/            input genome, GTF and marker table; NOT distributed
      results/         default output root
        tables/
        figures/
        scans/         PWM scan outputs
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(_HERE)                      # repository root

# ---------------------------------------------------------------- output
ROOT = os.environ.get("MITF_ROOT", os.path.join(BASE, "results"))
TAB = os.path.join(ROOT, "tables")
FIG = os.path.join(ROOT, "figures")
SCAN_DIR = os.environ.get("MITF_SCAN_DIR", os.path.join(ROOT, "scans"))

# ---------------------------------------------------------------- input
DATA = os.environ.get("MITF_DATA", os.path.join(BASE, "data"))
MARKERS = os.environ.get(
    "MITF_MARKERS", os.path.join(DATA, "seu.Xesp.markers.cts.csv"))
GENOME = os.environ.get(
    "MITF_GENOME", os.path.join(DATA, "xenia.contigs.FINAL.fasta"))
GTF = os.environ.get(
    "MITF_GTF", os.path.join(DATA, "Xesp_long.annot.gtf"))
MOTIF_DIR = os.environ.get("MITF_MOTIFS", os.path.join(BASE, "motifs"))

# cell-type nodes to drop from the analysis (matched as a prefix). The published run
# drops the ten unassigned ``unknown*`` clusters and retains 16 annotated populations.
EXCLUDE_PREFIX = os.environ.get("MITF_EXCLUDE_PREFIX", "unknown")

# ---------------------------------------------------------------- analysis constants
TARGET_NODE = "gastrodermis_alga_hosting"
PRIMARY_MOTIF = "MA0620.4"
PRIMARY_WIN = "prox2kb"

# marker definition (applied identically to every cell type)
PADJ_MAX = 0.05
LOG2FC_MIN = 0.25

# motifs shown in Figure 2 (ten JASPAR matrices + one scrambled control)
MOTIFS = ["MA0620.4", "MA1899.1", "MA2257.1", "MA0831.3", "MA0692.2",
          "SCD_MA0620.4", "MA0036.3", "MA0143.5", "MA0901.3", "MA1870.2",
          "MA1109.2"]

# motif id -> display label (Figure 2)
MOTIF_LABEL = {
    "MA0620.4": "MITF (human)", "MA1899.1": "Mitf (Ciona)",
    "MA2257.1": "Mitf (fly)", "MA0831.3": "TFE3", "MA0692.2": "TFEB",
    "SCD_MA0620.4": "scrambled MITF", "MA0036.3": "GATA2", "MA0143.5": "SOX2",
    "MA0901.3": "HOXB13", "MA1870.2": "FoxL2", "MA1109.2": "NEUROD1",
}
MITF_FAMILY = {"MA0620.4", "MA1899.1", "MA2257.1", "MA0831.3", "MA0692.2"}

for _d in (TAB, FIG):
    os.makedirs(_d, exist_ok=True)


def load_markers(drop_excluded=True):
    """Read the marker table and (by default) drop excluded cell-type nodes."""
    import pandas as pd
    m = pd.read_csv(MARKERS, sep="\t")
    m["gene_id"] = m["gene_seurat"].str.replace("-", "_", regex=False)
    if drop_excluded and EXCLUDE_PREFIX:
        keep = ~m["focus_node"].str.startswith(EXCLUDE_PREFIX)
        n_all = m["focus_node"].nunique()
        m = m[keep].copy()
        assert len(m) > 0, "node filter removed every row"
        print(f"[config] dropped nodes starting with '{EXCLUDE_PREFIX}': "
              f"{n_all - m['focus_node'].nunique()} of {n_all} nodes removed, "
              f"{m['focus_node'].nunique()} retained", flush=True)
    return m
