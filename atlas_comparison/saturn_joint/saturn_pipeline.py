#!/usr/bin/env python3
"""End-to-end cross-species SATURN pipeline: build the run table, train the
model, compute an integrated UMAP and draw the overview figures.

Stages
------
1. **Run table**  one row per species: ``path`` (its h5ad), ``species`` and
   ``embedding_path`` (the ESM2 protein-embedding ``.pt`` for that species).
2. **Train**      calls SATURN's ``train-saturn.py``; writes into
   ``<work_dir>/saturn_results/`` (the directory must exist beforehand).
3. **Embed**      loads the trained AnnData, checks that the integrated labels
   match the per-species input annotation, recomputes a UMAP on the macrogenes
   (PCA -> kNN -> UMAP, fixed seed) and caches it as
   ``umap_coordinates.tsv.gz``.
4. **Plot**       two figures, PDF:
   ``umap_by_celltype_no_labels.pdf``      all cells, coloured by cell type
   ``umap_species_facets_by_celltype.pdf`` one panel per species, its cells
                                           coloured by cell type, rest grey

Stage 3 also writes ``celltype_palette.tsv`` -- one colour per integrated cell
type, pooled from the ``color`` column of the per-species input atlases so that
a cell type keeps the same colour across every figure.

Typical use
-----------
    # edit the DATA PATHS block below (or pass the matching flags), then
    python saturn_pipeline.py --work-dir run_1127 --title "11-species SATURN"

    # figures only, reusing the cached coordinates (seconds, no h5ad needed)
    python saturn_pipeline.py --work-dir run_1127 --reuse-coords

    # skip training, reuse an existing SATURN h5ad
    python saturn_pipeline.py --saturn-h5ad path/to/output_saturn_seed_0.h5ad

Requirements: SATURN (checked out), pytorch with a GPU, and
numpy / pandas / matplotlib / scanpy / anndata.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

import anndata as ad
import scanpy as sc

# ---------------------------------------------------------------------------
# DATA PATHS -- replace with your own locations (or pass the matching flags)
# ---------------------------------------------------------------------------
DEFAULT_SATURN_ROOT = "/path/to/SATURN"                 # SATURN repository
DEFAULT_H5AD_DIR = "/path/to/saturn_inputs/final_h5ad"  # one <species>.h5ad per species
DEFAULT_EMBEDDING_DIR = "/path/to/protein_embedding/cross_species"

# per-species embedding file, filled in with --embedding-dir and the species name
EMBEDDING_TEMPLATE = (
    "{embedding_dir}/{species}_data/"
    "{species}_renamed_cleaned_longest_gene_symbol_to_embedding_ESM2.pt"
)

# Columns the input atlases must provide
LABEL_COL = "cell_type"     # per-species cell-type annotation
COLOR_COL = "color"         # optional: per-species cell-type colour

# ---------------------------------------------------------------------------
# Display convention for unannotated clusters
# ---------------------------------------------------------------------------
UNKNOWN_PREFIX = "unknown"
UNKNOWN_LABEL = "unknown"
UNKNOWN_COLOR = "#8C8C8C"

# R colour names used by the Seurat atlases -> hex
R_COLOR_TO_HEX = {
    "burlywood4": "#8B7355", "cadetblue2": "#8EE5EE", "cadetblue3": "#7AC5CD",
    "cadetblue4": "#53868B", "darkolivegreen3": "#A2CD5A", "khaki2": "#EEE685",
    "khaki3": "#CDC673", "mediumorchid3": "#B452CD", "mediumpurple1": "#AB82FF",
    "mediumpurple2": "#9F79EE", "pink3": "#CD919E", "plum1": "#FFBBFF",
    "plum2": "#EEAEEE", "plum3": "#CD96CD", "plum4": "#8B668B",
    "royalblue1": "#4876FF", "royalblue2": "#436EEE", "royalblue3": "#3A5FCD",
    "royalblue4": "#27408B", "sienna1": "#FF8247", "sienna2": "#EE7942",
    "sienna3": "#CD6839", "slategray3": "#9FB6CD", "steelblue1": "#63B8FF",
    "steelblue2": "#5CACEE", "steelblue3": "#4F94CD", "wheat1": "#FFE7BA",
    "wheat3": "#CDBA96",
}

RUN_LOG: Path | None = None


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def log(message: str) -> None:
    line = f"[{now_iso()}] {message}"
    print(line, flush=True)
    if RUN_LOG is not None:
        with RUN_LOG.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")


def natural_key(value: str) -> tuple:
    """Sort 'unknown2' before 'unknown10'."""
    return tuple(int(p) if p.isdigit() else p
                 for p in re.split(r"(\d+)", str(value)))


def color_to_hex(value: str) -> str:
    color = str(value)
    if re.fullmatch(r"#[0-9A-Fa-f]{6}", color):
        return color.upper()
    if color in R_COLOR_TO_HEX:
        return R_COLOR_TO_HEX[color]
    return ""          # unknown colour name: caller assigns one from a palette


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"wrote {stem.with_suffix('.pdf')}")


# ---------------------------------------------------------------------------
# Display convention for `unknown*`
# ---------------------------------------------------------------------------
def merge_unknown_labels(labels: np.ndarray) -> np.ndarray:
    series = pd.Series(np.asarray(labels)).astype(str)
    series = series.where(~series.str.startswith(UNKNOWN_PREFIX), UNKNOWN_LABEL)
    return series.to_numpy()


def merge_unknown_palette(palette: dict[str, str], labels: np.ndarray) -> dict[str, str]:
    """Ordered by abundance, with every `unknown*` pooled into one grey entry."""
    counts = pd.Series(np.asarray(labels)).astype(str).value_counts()
    ordered = sorted(counts.index,
                     key=lambda ct: (-int(counts[ct]), natural_key(ct)))
    missing = [ct for ct in ordered
               if ct != UNKNOWN_LABEL and ct not in palette]
    if missing:
        raise ValueError(f"cell types without a palette colour: {missing}")
    return {ct: UNKNOWN_COLOR if ct == UNKNOWN_LABEL else palette[ct]
            for ct in ordered}


# ---------------------------------------------------------------------------
# Stage 1: run table
# ---------------------------------------------------------------------------
def build_run_table(h5ad_dir: Path, embedding_dir: Path, species: list[str],
                    embedding_template: str, run_csv: Path) -> pd.DataFrame:
    rows = []
    for sp in species:
        h5ad_path = (h5ad_dir / f"{sp}.h5ad").resolve()
        embedding_path = Path(embedding_template.format(
            embedding_dir=embedding_dir, species=sp)).resolve()
        for label, path in (("h5ad", h5ad_path), ("embedding", embedding_path)):
            if not path.exists():
                raise FileNotFoundError(f"{sp}: {label} not found: {path}")
        rows.append({"path": str(h5ad_path), "species": sp,
                     "embedding_path": str(embedding_path)})
    table = pd.DataFrame(rows, columns=["path", "species", "embedding_path"])
    table.to_csv(run_csv, index=False)
    log(f"run table: {run_csv} ({len(table)} species)")
    return table


# ---------------------------------------------------------------------------
# Stage 2: training
# ---------------------------------------------------------------------------
def run_training(saturn_root: Path, run_csv: Path, work_dir: Path,
                 args: argparse.Namespace) -> None:
    script = saturn_root / "train-saturn.py"
    if not script.exists():
        raise FileNotFoundError(f"train-saturn.py not found: {script}")
    (work_dir / "saturn_results").mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable, str(script),
        f"--in_data={run_csv}",
        f"--in_label_col={args.label_col}",
        f"--ref_label_col={args.label_col}",
        f"--num_macrogenes={args.num_macrogenes}",
        f"--hv_genes={args.hv_genes}",
        f"--work_dir={work_dir}",
        f"--device_num={args.device}",
    ]
    if args.centroids_init_path:
        cmd.append(f"--centroids_init_path={args.centroids_init_path}")
    cmd.extend(args.saturn_extra_args or [])

    log("training: " + " ".join(cmd))
    if args.dry_run:
        log("--dry-run: training not executed")
        return
    start = time.time()
    subprocess.run(cmd, check=True, cwd=work_dir)
    log(f"training finished in {(time.time() - start) / 60:.1f} min")


def find_saturn_output(work_dir: Path) -> Path:
    """Newest final SATURN h5ad (ignores the _ep_<n> and _pretrain checkpoints)."""
    candidates = [p for p in (work_dir / "saturn_results").glob("*.h5ad")
                  if "_ep_" not in p.name and "_pretrain" not in p.name]
    if not candidates:
        raise FileNotFoundError(
            f"no SATURN output h5ad in {work_dir / 'saturn_results'}")
    latest = max(candidates, key=lambda p: p.stat().st_mtime)
    log(f"SATURN output: {latest}")
    return latest


# ---------------------------------------------------------------------------
# Stage 3: cell-type palette, validation, UMAP
# ---------------------------------------------------------------------------
def read_input_annotation(h5ad_dir: Path, species_order: list[str],
                          label_col: str = LABEL_COL,
                          color_col: str = COLOR_COL
                          ) -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Per-species cell_type sequence (for validation) and a long table of
    (species, cell_type, colour, n_cells)."""
    expected: dict[str, list[str]] = {}
    rows: list[dict] = []
    for sp in species_order:
        path = h5ad_dir / f"{sp}.h5ad"
        adata = ad.read_h5ad(path, backed="r")
        if label_col not in adata.obs.columns:
            raise ValueError(f"{path}: missing obs column {label_col!r}")
        cell_type = adata.obs[label_col].astype(str).to_numpy()
        expected[sp] = cell_type.tolist()
        if color_col in adata.obs.columns:
            color = adata.obs[color_col].astype(str).to_numpy()
        else:
            color = np.array([""] * adata.n_obs)
        counts = pd.crosstab(pd.Series(cell_type, name="cell_type"),
                             pd.Series(color, name="color"))
        for ct in counts.index:
            for col in counts.columns:
                n = int(counts.loc[ct, col])
                if n:
                    rows.append({"species": sp, "cell_type": ct,
                                 "color_r": col, "n_cells": n})
        adata.file.close()
        log(f"{sp}: {adata.n_obs} cells, {len(set(cell_type))} cell types")
    return expected, pd.DataFrame(rows)


def build_celltype_palette(sources: pd.DataFrame) -> dict[str, str]:
    """One colour per cell type, taken from the species with most cells for it.
    Cell types without a usable colour name get one from a fallback palette."""
    fallback = list(sc.pl.palettes.godsnot_102)
    used: set[str] = set()
    palette: dict[str, str] = {}
    totals = sources.groupby("cell_type")["n_cells"].sum()

    for ct in sorted(totals.index, key=lambda c: (-int(totals[c]), natural_key(c))):
        subset = sources[sources["cell_type"] == ct]
        best = subset.sort_values("n_cells", ascending=False)
        hexes = [color_to_hex(v) for v in best["color_r"]]
        hexes = [h for h in hexes if h]
        if len(best["color_r"].unique()) > 1:
            log(f"PALETTE {ct}: {len(best['color_r'].unique())} colours across "
                f"species; using {best.iloc[0]['color_r']} (most cells)")
        if hexes:
            palette[ct] = hexes[0]
        else:
            for candidate in fallback:
                if candidate.upper() not in used:
                    palette[ct] = candidate.upper()
                    break
            else:
                raise ValueError(f"no colour left for cell type {ct!r}")
        used.add(palette[ct].upper())
    return palette


def validate_labels(adata: ad.AnnData, expected: dict[str, list[str]],
                    species_order: list[str]) -> None:
    observed_order = list(dict.fromkeys(adata.obs["species"].astype(str)))
    if observed_order != species_order:
        log(f"WARNING species order {observed_order} != input {species_order}")
    mismatched = []
    for sp in species_order:
        observed = adata.obs.loc[
            adata.obs["species"].astype(str).eq(sp), "labels2"].astype(str).tolist()
        if observed != expected[sp]:
            mismatched.append(sp)
    if mismatched:
        log(f"WARNING labels2 does not match the input annotation for "
            f"{mismatched}; check that cell order and label column are correct")
    else:
        log("validation: labels2 matches the per-species input annotation")


def compute_umap(adata: ad.AnnData, seed: int) -> None:
    n_comps = min(50, adata.n_vars - 1, adata.n_obs - 1)
    if n_comps < 2:
        raise ValueError(f"cannot compute PCA with shape {adata.shape}")
    sc.pp.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_comps, metric="euclidean",
                    random_state=seed)
    sc.tl.umap(adata, min_dist=0.5, spread=1.0, random_state=seed)


def embed(saturn_h5ad: Path, h5ad_dir: Path, output_dir: Path, seed: int,
          use_macrogenes: bool = False, label_col: str = LABEL_COL,
          color_col: str = COLOR_COL
          ) -> tuple[pd.DataFrame, dict[str, str]]:
    log(f"loading {saturn_h5ad}")
    adata = sc.read_h5ad(saturn_h5ad)
    for column in ("species", "labels2"):
        if column not in adata.obs.columns:
            raise ValueError(f"{saturn_h5ad}: missing obs column {column!r}")

    species = adata.obs["species"].astype(str)
    species_order = list(dict.fromkeys(species))
    expected, sources = read_input_annotation(h5ad_dir, species_order,
                                              label_col=label_col,
                                              color_col=color_col)
    palette = build_celltype_palette(sources)
    sources.to_csv(output_dir / "celltype_palette_sources.tsv", sep="\t",
                   index=False)
    pd.DataFrame([{"cell_type": ct, "color_hex": palette[ct]}
                  for ct in palette]).to_csv(
        output_dir / "celltype_palette.tsv", sep="\t", index=False)

    validate_labels(adata, expected, species_order)

    # SATURN writes the integrated representation as .X (one column per
    # macrogene actually used, e.g. 256), and that is what the paper's UMAP was
    # built from.  obsm['macrogenes'] holds the full centroid matrix, which can
    # have more columns than were used; only use it on request.
    if use_macrogenes:
        if "macrogenes" not in adata.obsm:
            raise ValueError(f"{saturn_h5ad}: no obsm['macrogenes']")
        log(f"using obsm['macrogenes'] ({adata.obsm['macrogenes'].shape[1]} dims)")
        adata = ad.AnnData(X=np.asarray(adata.obsm["macrogenes"]),
                           obs=adata.obs.copy())
    log(f"UMAP on {adata.shape[1]} features (PCA -> kNN -> UMAP, seed={seed})")
    compute_umap(adata, seed)

    labels = adata.obs["labels2"].astype(str).to_numpy()
    xy = np.asarray(adata.obsm["X_umap"])
    coords = pd.DataFrame({
        "species": species.to_numpy(),
        "labels2": labels,
        "color_hex": pd.Series(labels).map(palette).to_numpy(),
        "UMAP1": xy[:, 0],
        "UMAP2": xy[:, 1],
    })
    coords.to_csv(output_dir / "umap_coordinates.tsv.gz", sep="\t",
                  index=False, compression="gzip")
    log(f"wrote {output_dir / 'umap_coordinates.tsv.gz'} "
        f"({len(coords)} cells, {coords['labels2'].nunique()} cell types)")
    return coords, palette


# ---------------------------------------------------------------------------
# Stage 4: figures
# ---------------------------------------------------------------------------
def plot_by_celltype(coords: pd.DataFrame, palette: dict[str, str],
                     title: str, output_stem: Path,
                     merge_unknown: bool = True) -> None:
    labels = coords["labels2"].astype(str).to_numpy()
    if merge_unknown:
        labels = merge_unknown_labels(labels)
        palette = merge_unknown_palette(palette, labels)
    xy = coords[["UMAP1", "UMAP2"]].to_numpy()

    obs = pd.DataFrame(index=range(len(xy)))
    obs["cell_type"] = pd.Categorical(labels, categories=list(palette),
                                      ordered=True)
    adata = sc.AnnData(obs=obs)
    adata.obsm["X_umap"] = xy

    fig, ax = plt.subplots(figsize=(11.5, 8.5))
    point_size = max(2.0, min(8.0, 120000.0 / adata.n_obs))
    sc.pl.umap(adata, color="cell_type", palette=palette, size=point_size,
               alpha=0.82, frameon=False, title=title, ax=ax, show=False,
               legend_loc="right margin")

    n_categories = len(palette)
    legend = ax.get_legend()
    if legend is not None:
        fontsize = 5.5 if n_categories > 20 else 7.0
        legend.set_title("Cell type")
        legend.get_title().set_fontsize(fontsize + 0.5)
        for text in legend.get_texts():
            text.set_fontsize(fontsize)
        try:
            legend.set_ncols(2 if n_categories > 18 else 1)
        except AttributeError:
            pass

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    for artist in ax.collections:
        artist.set_rasterized(True)         # keeps the PDF small
    fig.subplots_adjust(right=0.76)
    save_figure(fig, output_stem)


def plot_species_facets(coords: pd.DataFrame, palette: dict[str, str],
                        title: str, output_stem: Path,
                        merge_unknown: bool = True) -> None:
    labels = coords["labels2"].astype(str).to_numpy()
    if merge_unknown:
        labels = merge_unknown_labels(labels)
        palette = merge_unknown_palette(palette, labels)
    xy = coords[["UMAP1", "UMAP2"]].to_numpy()
    species = coords["species"].astype(str).to_numpy()
    species_order = list(dict.fromkeys(species))
    categories = list(palette)

    xlim = (float(xy[:, 0].min()), float(xy[:, 0].max()))
    ylim = (float(xy[:, 1].min()), float(xy[:, 1].max()))

    ncols = 4
    n = len(species_order)
    nrows = int(np.ceil(n / ncols))
    fig = plt.figure(figsize=(3.4 * ncols + 2.6, 3.2 * nrows))
    grid = GridSpec(nrows, ncols + 1, figure=fig,
                    width_ratios=[1] * ncols + [0.72],
                    wspace=0.08, hspace=0.16, top=0.94)
    fig.suptitle(title, fontsize=12, y=0.985)

    for i, sp in enumerate(species_order):
        ax = fig.add_subplot(grid[i // ncols, i % ncols])
        in_species = species == sp
        ax.scatter(xy[~in_species, 0], xy[~in_species, 1], s=1.0, c="#DDDDDD",
                   alpha=0.35, linewidths=0, rasterized=True, zorder=1)
        present = [ct for ct in categories
                   if (in_species & (labels == ct)).any()]
        present.sort(key=lambda ct: -int((in_species & (labels == ct)).sum()))
        for ct in present:
            mask = in_species & (labels == ct)
            ax.scatter(xy[mask, 0], xy[mask, 1], s=4, c=palette[ct],
                       alpha=0.95, linewidths=0, rasterized=True, zorder=2)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_title(f"{sp} (n={int(in_species.sum())}, {len(present)} types)",
                     fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        if i % ncols == 0:
            ax.set_ylabel("UMAP2", fontsize=8)
        if i >= n - ncols:
            ax.set_xlabel("UMAP1", fontsize=8)

    legend_ax = fig.add_subplot(grid[:, ncols])
    legend_ax.axis("off")
    handles = [Line2D([0], [0], marker="o", ls="", markersize=4.5,
                      markerfacecolor=palette[ct], markeredgecolor="none",
                      label=ct) for ct in categories if (labels == ct).any()]
    legend = legend_ax.legend(handles=handles, loc="center left", frameon=False,
                              fontsize=6.0, title="Cell type", title_fontsize=8,
                              labelspacing=0.28, handletextpad=0.4)
    for text in legend.get_texts():
        text.set_fontsize(6.0)
    save_figure(fig, output_stem)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cross-species SATURN: train, embed and plot.")
    # data locations
    parser.add_argument("--saturn-root", type=Path, default=Path(DEFAULT_SATURN_ROOT))
    parser.add_argument("--h5ad-dir", type=Path, default=Path(DEFAULT_H5AD_DIR))
    parser.add_argument("--embedding-dir", type=Path,
                        default=Path(DEFAULT_EMBEDDING_DIR))
    parser.add_argument("--embedding-template", default=EMBEDDING_TEMPLATE)
    parser.add_argument("--species", nargs="*", default=None,
                        help="species names, in run order (default: every "
                             "<species>.h5ad found in --h5ad-dir, sorted)")

    # training
    parser.add_argument("--work-dir", type=Path, default=Path("."))
    parser.add_argument("--label-col", default=LABEL_COL,
                        help="cell-type column in the per-species h5ad "
                             "and the label used to train SATURN")
    parser.add_argument("--color-col", default=COLOR_COL,
                        help="optional colour column in the input h5ad")
    parser.add_argument("--num-macrogenes", type=int, default=1000)
    parser.add_argument("--hv-genes", type=int, default=4000)
    parser.add_argument("--device", default="0")
    parser.add_argument("--centroids-init-path", default=None)
    parser.add_argument("--saturn-extra-args", nargs=argparse.REMAINDER,
                        default=None, help="extra arguments for train-saturn.py")

    # stage control
    parser.add_argument("--saturn-h5ad", type=Path, default=None,
                        help="skip training and use this h5ad")
    parser.add_argument("--skip-train", action="store_true",
                        help="reuse the newest h5ad in <work-dir>/saturn_results")
    parser.add_argument("--reuse-coords", action="store_true",
                        help="plot from the cached umap_coordinates.tsv.gz only")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the training command without running it")
    parser.add_argument("--use-macrogenes", action="store_true",
                        help="embed on obsm['macrogenes'] instead of .X")
    parser.add_argument("--keep-unknown-split", action="store_true",
                        help="do not merge `unknown*` into one grey category")

    # output
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--title", default="Multi-species SATURN integration")
    return parser.parse_args()


def main() -> int:
    global RUN_LOG
    args = parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    RUN_LOG = args.work_dir / "pipeline.log"
    log(f"work dir: {args.work_dir.resolve()}")

    coords_path = args.work_dir / "umap_coordinates.tsv.gz"
    palette_path = args.work_dir / "celltype_palette.tsv"

    if args.reuse_coords:
        if not coords_path.exists() or not palette_path.exists():
            raise SystemExit("--reuse-coords needs umap_coordinates.tsv.gz and "
                             "celltype_palette.tsv in --work-dir")
        coords = pd.read_csv(coords_path, sep="\t")
        palette = pd.read_csv(palette_path, sep="\t").set_index(
            "cell_type")["color_hex"].to_dict()
        log(f"reusing cached coordinates ({len(coords)} cells)")
    else:
        h5ad = args.saturn_h5ad
        if h5ad is None and not args.skip_train:
            species = args.species
            if not species:
                species = sorted(p.stem for p in args.h5ad_dir.glob("*.h5ad"))
                if not species:
                    raise SystemExit(f"no <species>.h5ad found in {args.h5ad_dir}")
                log(f"species auto-detected: {species}")
            run_csv = args.work_dir / "cross_species_run.csv"
            build_run_table(args.h5ad_dir, args.embedding_dir, species,
                            args.embedding_template, run_csv)
            run_training(args.saturn_root, run_csv, args.work_dir, args)
            if args.dry_run:
                return 0
        if h5ad is None:
            h5ad = find_saturn_output(args.work_dir)
        coords, palette = embed(h5ad, args.h5ad_dir, args.work_dir, args.seed,
                                use_macrogenes=args.use_macrogenes,
                                label_col=args.label_col,
                                color_col=args.color_col)

    merge_unknown = not args.keep_unknown_split
    plot_by_celltype(coords, palette,
                     f"{args.title} | integrated cell type",
                     args.work_dir / "umap_by_celltype_no_labels",
                     merge_unknown=merge_unknown)
    plot_species_facets(coords, palette,
                        f"{args.title} | per species, coloured by cell type",
                        args.work_dir / "umap_species_facets_by_celltype",
                        merge_unknown=merge_unknown)

    summary = {
        "status": "PASS",
        "title": args.title,
        "n_cells": int(len(coords)),
        "n_species": int(coords["species"].nunique()),
        "n_cell_types": int(coords["labels2"].nunique()),
        "species_counts": coords["species"].value_counts().astype(int).to_dict(),
        "figures": ["umap_by_celltype_no_labels.pdf",
                    "umap_species_facets_by_celltype.pdf"],
        "completed_at": now_iso(),
    }
    (args.work_dir / "_SUCCESS.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
