#!/usr/bin/env python3
"""Step 5: Figures 1 and 2.

Figure 1  fraction of significant marker genes carrying >= 1 predicted MITF-family
          E-box element, by cell type.
Figure 2  odds ratio for the symbiont-hosting gastrodermis, across eleven matrices.

Palette: validated categorical slots from the dataviz reference palette
(blue #2a78d6 = signal, neutral greys = everything else).
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mitf_config as C

TAB = C.TAB
FIG = C.FIG
os.makedirs(FIG, exist_ok=True)

BLUE, MUTED = "#2a78d6", "#b9b8b2"
INK, INK2 = "#0b0b0b", "#52514e"
SURFACE = "#fcfcfb"
TARGET = C.TARGET_NODE

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.major.size": 2.5, "ytick.major.size": 0,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})


def orr_ci(a, b, c, d):
    """Woolf 95% CI for the odds ratio."""
    a, b, c, d = a + .5, b + .5, c + .5, d + .5
    lor = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return np.exp(lor), np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se)


# ------------------------------------------------------------------ Figure 1
def fig_clusters():
    d = pd.read_csv(os.path.join(TAB, "cluster_mitf_enrichment.tsv"), sep="\t")
    d = d[d.gene_set == "sig_markers"].sort_values("frac_with_hit")
    bg = d.bg_frac_with_hit.median()

    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    y = np.arange(len(d))
    colors = [BLUE if c == TARGET else MUTED for c in d.cluster]
    ax.barh(y, d.frac_with_hit * 100, height=.62, color=colors, zorder=3)
    ax.axvline(bg * 100, color=INK2, lw=1.2, ls=(0, (4, 3)), zorder=4)
    ax.text(bg * 100 - .6, len(d) - .45, f"background {bg*100:.1f}%", color=INK2,
            fontsize=8, va="bottom", ha="right")

    ax.set_yticks(y, [c.replace("_", " ") for c in d.cluster], fontsize=8.5)
    for tick, c in zip(ax.get_yticklabels(), d.cluster):
        if c == TARGET:
            tick.set_color(BLUE)
            tick.set_fontweight("bold")
    ax.set_xlabel("% of marker genes with ≥ 1 MITF-family E-box element in the promoter")
    ax.set_xlim(0, 42)
    ax.set_ylim(-1.1, len(d) + .35)
    ax.grid(axis="x", color=MUTED, alpha=.35, lw=.6, zorder=0)
    ax.set_axisbelow(True)

    for i, (_, r) in enumerate(d.iterrows()):
        if r.cluster == TARGET:
            ax.text(r.frac_with_hit * 100 + .7, i,
                    f"{r.frac_with_hit*100:.1f}%   OR {r.odds_ratio:.2f}   "
                    f"p = {r.fisher_p:.1e}",
                    va="center", fontsize=8.5, color=BLUE, fontweight="bold")
    ax.legend(handles=[Line2D([], [], color=BLUE, lw=6,
                              label="symbiont-hosting gastrodermis"),
                       Line2D([], [], color=MUTED, lw=6,
                              label="other cell types")],
              loc="lower right", frameon=False, fontsize=8.5)
    ax.set_title("MITF-family E-box enrichment in marker genes, by cell type",
                 loc="left", fontsize=10.5, fontweight="bold", pad=12)
    fig.savefig(os.path.join(FIG, "fig1_cluster_enrichment.pdf"))
    plt.close(fig)


# ------------------------------------------------------------------ Figure 2
def fig_motifs():
    r = pd.read_csv(os.path.join(TAB, "motif_robustness.tsv"), sep="\t")
    r = r[(r.gene_set == "sig_markers") & (r.window == "prox2kb")].copy()
    CI = [orr_ci(int(row.n_genes * row.frac_with_hit),
                 row.n_genes - int(row.n_genes * row.frac_with_hit),
                 int(row.n_bg_with_hit),
                 row.n_bg - int(row.n_bg_with_hit))
          for row in r.itertuples()]
    r["odds_ratio"], r["lo"], r["hi"] = zip(*CI)
    r = r.sort_values("odds_ratio").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    y = np.arange(len(r))
    for i, row in r.iterrows():
        col = BLUE if row.motif in C.MITF_FAMILY else MUTED
        ax.plot([row.lo, row.hi], [i, i], color=col, lw=2,
                solid_capstyle="round", zorder=3)
        ax.scatter(row.odds_ratio, i, s=64, color=col, zorder=4)
    ax.axvline(1, color=INK2, lw=1.2, ls=(0, (4, 3)), zorder=2)

    ax.set_yticks(y, [C.MOTIF_LABEL.get(m, m) for m in r.motif], fontsize=8.5)
    for tick, m in zip(ax.get_yticklabels(), r.motif):
        if m in C.MITF_FAMILY:
            tick.set_color(BLUE)
            tick.set_fontweight("bold")
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    for i, (_, row) in enumerate(r.iterrows()):
        ax.text(1.04, i, f"p = {row.fisher_p:.1e}", va="center", ha="left",
                fontsize=8.2, color=INK2, transform=tr)
    ax.set_xlim(.8, 2.05)
    ax.set_xlabel("odds ratio for E-box enrichment in symbiont-hosting gastrodermis markers")
    ax.grid(axis="x", color=MUTED, alpha=.35, lw=.6, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, ms=8,
                              label="MITF / TFE family PWM"),
                       Line2D([], [], marker="o", ls="", color=MUTED, ms=8,
                              label="unrelated PWMs & scrambled control")],
              loc="upper left", bbox_to_anchor=(0, -0.13), ncol=2,
              frameon=False, fontsize=8.3)
    ax.set_title("Only MITF-family E-box motifs are enriched",
                 loc="left", fontsize=10.5, fontweight="bold", pad=12)
    fig.savefig(os.path.join(FIG, "fig2_motif_robustness.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig_clusters()
    fig_motifs()
    print("figures written to", FIG)
