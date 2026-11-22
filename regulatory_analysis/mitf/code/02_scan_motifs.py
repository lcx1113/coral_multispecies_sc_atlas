#!/usr/bin/env python3
"""
Step 2: PWM scanner.

Scores every position of every input sequence with a JASPAR PWM (both strands,
or one strand for palindromic matrices), computes the background score
distribution by discretised convolution of the per-position score distributions,
converts scores to p-values, and reports hits below a p-value threshold together
with per-sequence BH q-values.

Usage:
  python3 02_scan_motifs.py --motif motifs/MA0620.4.jaspar --out results/mitf_MA0620.4
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mitf_config as C

BASES = "ACGT"
BASE_IDX = np.full(256, 4, dtype=np.uint8)          # case-insensitive
BASE_IDX_UPPER = np.full(256, 4, dtype=np.uint8)    # lowercase (=soft-masked repeat) -> N
for i, b in enumerate(BASES):
    BASE_IDX[ord(b)] = i
    BASE_IDX[ord(b.lower())] = i
    BASE_IDX_UPPER[ord(b)] = i

SCALE = 1000


def read_jaspar(path):
    """returns (matrix_id, name, counts[4 x L] as float arrays)"""
    mid = name = None
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].split()
                mid, name = f[0], (f[1] if len(f) > 1 else f[0])
            elif line.strip():
                f = line.split()
                rows[f[0]] = [float(x) for x in line.split("[")[1].rstrip().rstrip("]").split()]
    mat = np.array([rows[b] for b in BASES], dtype=float)
    return mid, name, mat


def ppm_from_counts(counts, pseudocount):
    c = counts + pseudocount
    return c / c.sum(axis=0, keepdims=True)


def log_odds(ppm, bg):
    return np.log2(np.maximum(ppm, 1e-9) / bg[:, None])


def score_distribution(S, bg, scale=SCALE):
    """Background distribution of PWM scores under `bg`, on a discretised grid.
    Returns (lo_score, dist) with dist[k] = P(score == lo + k/scale)."""
    L = S.shape[1]
    per = []
    for i in range(L):
        d = {}
        for b in range(4):
            off = int(round(S[b, i] * scale))
            d[off] = d.get(off, 0.0) + bg[b]
        per.append(d)
    lo = sum(min(d) for d in per)
    hi = sum(max(d) for d in per)
    size = hi - lo + 1
    dist = np.zeros(size)
    dist[-lo] = 1.0                      # score 0 sits at index -lo
    for d in per:
        new = np.zeros(size)
        for off, p in d.items():
            if p <= 0 or off >= size or -off >= size:
                continue
            if off >= 0:
                new[off:] += dist[:size - off] * p
            else:
                new[:size + off] += dist[-off:] * p
        dist = new
    s = dist.sum()
    if s > 0:
        dist /= s
    return lo, dist


def make_cum(dist):
    """cum[k] = P(score >= lo + k/scale)"""
    return np.concatenate([np.cumsum(dist[::-1])[::-1], [0.0]])


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def gene_bg_composition(fasta):
    cnt = {"A": 0, "C": 0, "G": 0, "T": 0}
    total = 0
    with open(fasta) as fh:
        for line in fh:
            if line[0] == ">":
                continue
            s = line.strip().upper()
            for b in BASES:
                cnt[b] += s.count(b)
            total += sum(s.count(b) for b in BASES)
    return np.array([cnt[b] / total for b in BASES])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--motif", required=True)
    ap.add_argument("--fasta", default=os.path.join(C.ROOT, "promoters.fa"))
    ap.add_argument("--tss", default=os.path.join(C.ROOT, "gene_tss.tsv"))
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--pval", type=float, default=1e-4)
    ap.add_argument("--pseudocount", type=float, default=0.5)
    ap.add_argument("--bg", default="uniform", help="uniform | genome")
    ap.add_argument("--windows", default="", help="comma list to restrict, e.g. prox2kb")
    ap.add_argument("--genes", default="", help="optional file with gene_id list to restrict scanning")
    ap.add_argument("--mask-lower", action="store_true",
                    help="treat soft-masked (lowercase = repeat) bases as N and ignore them")
    args = ap.parse_args()

    base_idx = BASE_IDX if not args.mask_lower else BASE_IDX_UPPER

    bg = np.full(4, 0.25)
    if args.bg == "genome":
        bg = gene_bg_composition(args.fasta)
        sys.stderr.write("background composition: " +
                         ", ".join(f"{b}={v:.4f}" for b, v in zip(BASES, bg)) + "\n")

    mid, mname, counts = read_jaspar(args.motif)
    ppm = ppm_from_counts(counts, args.pseudocount)
    S = log_odds(ppm, bg)
    Sr = log_odds(ppm[::-1, ::-1], bg)          # reverse complement = match other strand
    L = S.shape[1]
    sys.stderr.write(f"motif {mid} {mname}, length {L}\n")

    lo, dist = score_distribution(S, bg)
    cum = make_cum(dist)
    size = len(dist)
    lo2, dist2 = score_distribution(Sr, bg)
    cum2 = make_cum(dist2)
    size2 = len(dist2)

    # if the motif equals its own reverse complement, forward and reverse scans
    # report the same sites -> scan one strand only to avoid double counting
    maxdev = float(np.abs(S - Sr).max())
    palindromic = maxdev <= 0.1
    sys.stderr.write(f"max |PWM - revcomp(PWM)| = {maxdev:.4f} log2 units -> "
                     f"{'palindromic, scanning forward strand only' if palindromic else 'not palindromic, scanning both strands'}\n")

    def pvals(scores, lo_, cum_, size_):
        idx = np.clip(np.rint(scores * SCALE).astype(np.int64) - lo_, 0, size_ - 1)
        return cum_[idx]

    # sequence metadata from the tss table
    meta = {}
    for i, line in enumerate(open(args.tss)):
        if i == 0:
            continue
        f = line.rstrip("\n").split("\t")
        meta[(f[0], f[5])] = dict(tid=f[1], chrom=f[2], strand=f[3], tss0=int(f[4]),
                                  start=int(f[6]), end=int(f[7]), win_len=int(f[8]),
                                  clipped=int(f[9]), tss_off=int(f[10]))
    wanted_genes = None
    if args.genes:
        wanted_genes = {l.split()[0] for l in open(args.genes) if l.strip()}
    win_filter = set(args.windows.split(",")) if args.windows else None

    hits_fh = open(args.out + ".hits.tsv", "w")
    hits_fh.write("gene_id\twindow\tchrom\tgenomic_start\tgenomic_end\tmotif_strand\t"
                  "score\tp_value\tq_value\tseq_pos\ttss_relative\tn_scan_pos\n")
    gene_fh = open(args.out + ".gene.tsv", "w")
    gene_fh.write("gene_id\twindow\tseq_len\tn_hits\tn_hits_q05\tbest_score\tbest_p\tbest_q\t"
                  "best_genomic_start\tbest_motif_strand\tbest_tss_relative\tpromoter_len\tn_scan_pos\n")

    n_seq = n_hit = 0
    for header, seq in read_fasta(args.fasta):
        f = header.split("|")
        gid, wname = f[0], f[1]
        if win_filter and wname not in win_filter:
            continue
        if wanted_genes is not None and gid not in wanted_genes:
            continue
        m = meta.get((gid, wname))
        if m is None or len(seq) < L:
            continue
        n_seq += 1
        idx = base_idx[np.frombuffer(seq.encode(), dtype=np.uint8)]
        npos = len(seq) - L + 1
        w = np.lib.stride_tricks.sliding_window_view(idx, L)
        bad = (w == 4).any(axis=1)
        ar = np.arange(L)
        wc = np.where(w == 4, 0, w)          # N -> dummy index for the gather

        sf = S[wc, ar].sum(axis=1)
        sf[bad] = -1e6
        parts_s, parts_p = [sf], [pvals(sf, lo, cum, size)]
        parts_strand, parts_pos = [np.zeros(npos, dtype=np.uint8)], [np.arange(npos)]

        if not palindromic:
            sr = Sr[wc, ar].sum(axis=1)
            sr[bad] = -1e6
            parts_s.append(sr)
            parts_p.append(pvals(sr, lo2, cum2, size2))
            parts_strand.append(np.ones(npos, dtype=np.uint8))
            parts_pos.append(np.arange(npos))

        alls = np.concatenate(parts_s)
        allp = np.concatenate(parts_p)
        allstrand = np.concatenate(parts_strand)
        allpos = np.concatenate(parts_pos)

        best = int(np.argmax(alls))
        hits = np.where(allp < args.pval)[0]
        n_hits_seq = len(hits)
        n_tests = len(alls)
        # Benjamini-Hochberg over all scanned positions in this sequence:
        # rank all positions, enforce step-up monotonicity, then read off the
        # q values of the reported hits.
        order = np.argsort(allp)
        qs = np.minimum(1.0, allp[order] * n_tests / np.arange(1, n_tests + 1))
        qs = np.minimum.accumulate(qs[::-1])[::-1]
        q_all = np.empty_like(qs)
        q_all[order] = qs
        best_q = float(q_all[best])
        if n_hits_seq:
            hits = hits[np.argsort(allp[hits])]
            q = q_all[hits]
            n_q05 = int((q < 0.05).sum())
        else:
            q = np.empty(0)
            n_q05 = 0

        gene_fh.write(
            f"{gid}\t{wname}\t{len(seq)}\t{n_hits_seq}\t{n_q05}\t"
            f"{alls[best]:.4f}\t{allp[best]:.3e}\t"
            f"{best_q:.3e}\t"
            f"{_gstart(m, int(allpos[best]))}\t"
            f"{_mstrand(m, int(allstrand[best]))}\t"
            f"{int(allpos[best]) + L // 2 - m['tss_off']}\t{m['win_len']}\t{n_tests}\n")

        if n_hits_seq:
            for k, h in enumerate(hits):
                pos, st = int(allpos[h]), int(allstrand[h])
                gs = _gstart(m, pos)
                hits_fh.write(
                    f"{gid}\t{wname}\t{m['chrom']}\t{gs}\t{gs + L}\t{_mstrand(m, st)}\t"
                    f"{alls[h]:.4f}\t{allp[h]:.3e}\t{q[k]:.3e}\t{pos}\t"
                    f"{pos + L // 2 - m['tss_off']}\t{n_tests}\n")
            n_hit += n_hits_seq

    hits_fh.close()
    gene_fh.close()
    sys.stderr.write(f"scanned {n_seq} sequences, {n_hit} hits at p<{args.pval}\n")
    sys.stderr.write(f"-> {args.out}.hits.tsv\n-> {args.out}.gene.tsv\n")


def _gstart(m, seq_pos):
    return m["start"] + seq_pos if m["strand"] == "+" else m["end"] - 1 - seq_pos


def _mstrand(m, seq_hit_strand):
    fwd = (m["strand"] == "+")          # sequence is in gene orientation
    on_plus = (seq_hit_strand == 0)
    return "+" if (on_plus == fwd) else "-"


if __name__ == "__main__":
    main()
