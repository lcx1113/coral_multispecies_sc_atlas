#!/usr/bin/env python3
"""Step 0: generate the scrambled MITF matrix used as a composition control.

Base labels are permuted independently within each column of MA0620.4. This preserves
the four base counts (and therefore the information content) of every column while
destroying the motif, so that any enrichment driven purely by per-position base
composition would survive in the scrambled control as well.

Seed 20261008 reproduces ``motifs/SCD_MA0620.4.jaspar`` exactly.

Usage:
  python3 00_make_scrambled_motif.py [--motif-dir DIR] [--seed 20261008]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mitf_config as C

BASES = "ACGT"
SOURCE = "MA0620.4"


def read_counts(path):
    rows = {}
    for line in open(path):
        if line.startswith(">") or not line.strip():
            continue
        rows[line.split()[0]] = [int(x) for x in
                                 line.split("[")[1].rstrip().rstrip("]").split()]
    return np.array([rows[b] for b in BASES], dtype=int)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--motif-dir", default=C.MOTIF_DIR)
    ap.add_argument("--seed", type=int, default=20261008)
    args = ap.parse_args()

    counts = read_counts(os.path.join(args.motif_dir, f"{SOURCE}.jaspar"))
    rng = np.random.default_rng(args.seed)
    out = np.empty_like(counts)
    for j in range(counts.shape[1]):
        out[:, j] = counts[rng.permutation(4), j]

    mid = f"SCD_{SOURCE}"
    path = os.path.join(args.motif_dir, f"{mid}.jaspar")
    with open(path, "w") as fh:
        fh.write(f">{mid}\tscrambled_{SOURCE}\n")
        for b, row in zip(BASES, out):
            fh.write(f"{b}  [ " + "  ".join(f"{v:6d}" for v in row) + " ]\n")
    print(f"-> {path}  (seed {args.seed}; per-column base counts preserved)")


if __name__ == "__main__":
    main()
