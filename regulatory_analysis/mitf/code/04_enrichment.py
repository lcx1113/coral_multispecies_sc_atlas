#!/usr/bin/env python3
"""
Step 4: enrichment analysis for Figures 1 and 2.

Figure 1 -- for every annotated cell type, test whether its significant marker genes
carry more predicted MITF-family E-box elements than the other genes in the same marker
table (two-sided Fisher exact test; BH correction across cell types).

Figure 2 -- repeat the test for the symbiont-hosting gastrodermis with eleven matrices
(five MITF/TFE-family matrices, five unrelated matrices and one scrambled control).

Outputs (into $MITF_ROOT/tables)
  cluster_mitf_enrichment.tsv   one row per cell type (Figure 1)
  motif_robustness.tsv          one row per matrix (Figure 2)
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mitf_config as C

SCAN_DIR = C.SCAN_DIR       # overridden by --scan-dir
MIN_WIN_LEN = {"prox2kb": 2000, "core500": 500}


def load_markers():
    # honours MITF_EXCLUDE_PREFIX via the shared config
    return C.load_markers()


def load_gene_scan(motif, window):
    f = os.path.join(SCAN_DIR, f"{motif}.{window}.gene.tsv")
    d = pd.read_csv(f, sep="\t")
    # one row per (gene, window); keep full-length windows only
    d = d[d["seq_len"] >= MIN_WIN_LEN[window]]
    return d


def hit_matrix(d, k=1):
    d = d.copy()
    d["has_hit"] = d["n_hits"] >= k
    d["hits_per_kb"] = d["n_hits"] / (d["seq_len"] / 1000.0)
    return d


def fisher(a, b, c, dd):
    """a=target with hit, b=target without, c=bg with, dd=bg without"""
    odds, p = stats.fisher_exact([[a, b], [c, dd]])
    return odds, p


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    qs = p[o] * n / np.arange(1, n + 1)
    qs = np.minimum.accumulate(qs[::-1])[::-1]      # enforce monotonicity in rank order
    q = np.empty(n)
    q[o] = np.minimum(qs, 1.0)
    return q


def fdr(p):
    """BH, using the SciPy implementation when available."""
    if hasattr(stats, "false_discovery_control"):
        return stats.false_discovery_control(p)
    return bh(p)


# significant, positively enriched markers; the same definition for every cell type
SIG = lambda g: g[(g.p_val_adj < C.PADJ_MAX) & (g.avg_log2FC > C.LOG2FC_MIN)]


def main():
    global SCAN_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", default=SCAN_DIR)
    ap.add_argument("--motif", default=C.PRIMARY_MOTIF)
    ap.add_argument("--window", default=C.PRIMARY_WIN)
    a = ap.parse_args()
    SCAN_DIR = a.scan_dir
    os.makedirs(C.TAB, exist_ok=True)

    mk = load_markers()
    scan = hit_matrix(load_gene_scan(a.motif, a.window))
    universe = set(mk.gene_id) & set(scan.gene_id)
    print(f"marker genes: {mk.gene_id.nunique()}, scanned with full promoter: "
          f"{len(scan)}, intersection: {len(universe)}")

    sc = scan.set_index("gene_id")
    universe_scan = sc.loc[sorted(universe)]

    # ------------------------------------------------------------ Figure 1
    rows = []
    for node, g in mk.groupby("focus_node"):
        ids = sorted(set(SIG(g).gene_id) & universe)
        if len(ids) < 20:
            continue
        tgt = sc.loc[ids]
        bg = universe_scan.drop(index=ids)
        a_ = int(tgt.has_hit.sum())
        b_ = len(tgt) - a_
        c_ = int(bg.has_hit.sum())
        d_ = len(bg) - c_
        orr, p = fisher(a_, b_, c_, d_)
        try:
            _, pu = stats.mannwhitneyu(tgt.best_score, bg.best_score,
                                       alternative="greater")
        except ValueError:
            pu = np.nan
        rows.append(dict(cluster=node, gene_set="sig_markers", n_genes=len(tgt),
                         n_with_hit=a_, frac_with_hit=a_ / len(tgt),
                         bg_frac_with_hit=c_ / len(bg),
                         odds_ratio=orr, fisher_p=p,
                         median_hits_per_kb=tgt.hits_per_kb.median(),
                         bg_median_hits_per_kb=bg.hits_per_kb.median(),
                         mwu_bestscore_p=pu))

    res = pd.DataFrame(rows)
    res["fisher_q"] = res.groupby("gene_set")["fisher_p"].transform(fdr)
    res = res.sort_values(["gene_set", "fisher_p"])
    res.to_csv(os.path.join(C.TAB, "cluster_mitf_enrichment.tsv"),
               sep="\t", index=False)

    print("\n=== per-cell-type MITF promoter enrichment (Figure 1) ===")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(res[["cluster", "n_genes", "n_with_hit", "frac_with_hit",
                   "bg_frac_with_hit", "odds_ratio", "fisher_p",
                   "fisher_q"]].to_string(index=False))

    # ------------------------------------------------------------ Figure 2
    rows = []
    for motif in C.MOTIFS:
        for win in ["prox2kb", "core500"]:
            try:
                s = hit_matrix(load_gene_scan(motif, win))
            except FileNotFoundError:
                continue
            s = s[s.gene_id.isin(set(mk.gene_id))]
            u = set(s.gene_id)
            t = set(SIG(mk[mk.focus_node == C.TARGET_NODE]).gene_id) & u
            sx = s.set_index("gene_id")
            tgt = sx.loc[sorted(t)]
            bg = sx.drop(index=list(t))
            a_ = int(tgt.has_hit.sum())
            b_ = len(tgt) - a_
            c_ = int(bg.has_hit.sum())
            d_ = len(bg) - c_
            orr, p = fisher(a_, b_, c_, d_)
            rows.append(dict(motif=motif, window=win, gene_set="sig_markers",
                             n_genes=len(tgt), frac_with_hit=a_ / len(tgt),
                             bg_frac_with_hit=c_ / len(bg),
                             n_bg=len(bg), n_bg_with_hit=c_,
                             odds_ratio=orr, fisher_p=p))

    rob = pd.DataFrame(rows).sort_values(["gene_set", "window", "fisher_p"])
    rob.to_csv(os.path.join(C.TAB, "motif_robustness.tsv"), sep="\t", index=False)

    print("\n=== robustness across matrices (Figure 2, prox2kb) ===")
    r = rob[(rob.gene_set == "sig_markers") &
            (rob.window == "prox2kb")].sort_values("odds_ratio", ascending=False)
    with pd.option_context("display.width", 200):
        print(r[["motif", "n_genes", "frac_with_hit", "bg_frac_with_hit",
                 "odds_ratio", "fisher_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
