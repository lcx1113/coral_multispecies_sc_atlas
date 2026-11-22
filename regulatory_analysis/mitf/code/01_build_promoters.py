#!/usr/bin/env python3
"""
Step 1: Build the TSS-relative promoter window used in Figures 1 and 2 for every gene
in the Xenia (Xesp) annotation, and extract the sequences from the Xenia genome.

Inputs (paths from mitf_config.py, overridable by environment variables)
  MITF_GTF     gene annotation, GTF
  MITF_GENOME  soft-masked genome assembly, FASTA

Outputs (into $MITF_ROOT)
  promoters.bed     all windows, BED6 + window name
  promoters.fa      sequences, one record per gene window, header carries metadata
  gene_tss.tsv      gene-level table: chrom, tss(0-based), strand, transcript used,
                    window coords, clipping, tss offset inside the extracted sequence
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mitf_config as C

RESULTS = C.ROOT
GTF = C.GTF
FASTA = C.GENOME

# window layout: name -> (upstream_bp, downstream_bp)
WINDOWS = {
    "prox2kb": (2000, 500),   # primary window used in Figures 1 and 2
}
MIN_FRAC_KEPT = 0.5  # drop windows clipped by more than half (contig edge)

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- genome
def load_genome(fasta):
    seqs = {}
    name = None
    buf = []
    with open(fasta) as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


# ---------------------------------------------------------------- gtf
def parse_gtf(gtf):
    """gene_id -> {tid: {'chrom','strand','tx_start','tx_end','exon_len','source'}}"""
    genes = defaultdict(dict)
    with open(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "exon"):
                continue
            attr = f[8]
            try:
                gid = attr.split('gene_id "')[1].split('"')[0]
                tid = attr.split('transcript_id "')[1].split('"')[0]
            except IndexError:
                continue
            if f[2] == "transcript":
                genes[gid][tid] = {
                    "chrom": f[0], "strand": f[6], "tx_start": int(f[3]),
                    "tx_end": int(f[4]), "exon_len": 0, "source": f[1],
                }
            else:
                if tid in genes[gid]:
                    genes[gid][tid]["exon_len"] += int(f[4]) - int(f[3]) + 1
    return genes


def pick_transcript(tx):
    """longest transcript by summed exon length (ties -> longest genomic span)."""
    def key(t):
        return (t["exon_len"], t["tx_end"] - t["tx_start"])
    best = max(tx.items(), key=lambda kv: key(kv[1]))
    return best[0], best[1]


def main():
    print("loading genome ...", flush=True)
    genome = load_genome(FASTA)
    print(f"  {len(genome)} sequences, {sum(len(v) for v in genome.values())/1e6:.1f} Mb",
          flush=True)

    print("parsing gtf ...", flush=True)
    genes = parse_gtf(GTF)
    print(f"  {len(genes)} gene_ids", flush=True)

    bed_fh = open(os.path.join(RESULTS, "promoters.bed"), "w")
    fa_fh = open(os.path.join(RESULTS, "promoters.fa"), "w")
    tab_fh = open(os.path.join(RESULTS, "gene_tss.tsv"), "w")
    tab_fh.write("gene_id\ttranscript_id\tchrom\tstrand\ttss0\twindow\t"
                 "bed_start\tbed_end\twin_len\tn_clipped\ttss_offset\tseq_len\tsource\n")

    # `Xesp_orphan_*` models are single-exon transcripts unique to this annotation.
    # They are kept by default so that the scanned promoter set covers every annotated
    # gene model; pass --drop-orphans to reproduce the earlier orphan-excluded run.
    drop_orphans = "--drop-orphans" in sys.argv

    n_ok = n_drop = n_nogene = n_orphan = 0
    for gid, tx in genes.items():
        if gid.startswith("Xesp_orphan"):
            n_orphan += 1
            if drop_orphans:
                n_nogene += 1
                continue
        tid, t = pick_transcript(tx)
        chrom, strand = t["chrom"], t["strand"]
        if chrom not in genome:
            n_nogene += 1
            continue
        clen = len(genome[chrom])
        # 0-based index of the TSS base
        tss0 = t["tx_start"] - 1 if strand == "+" else t["tx_end"] - 1

        for wname, (up, down) in WINDOWS.items():
            # genomic span that covers [-up, +down] relative to the TSS base, 0-based half-open
            if strand == "+":
                start, end = tss0 - up, tss0 + down + 1
            else:
                start, end = tss0 - down, tss0 + up + 1
            raw_len = end - start
            cs, ce = max(0, start), min(clen, end)
            kept = (ce - cs) / raw_len
            if kept < MIN_FRAC_KEPT or ce <= cs:
                n_drop += 1
                continue

            sub = genome[chrom][cs:ce]
            if strand == "-":
                seq = revcomp(sub)
                # index 0 of `seq` is genomic base ce-1
                tss_off = (ce - 1) - tss0
            else:
                seq = sub
                tss_off = tss0 - cs

            bed_fh.write(f"{chrom}\t{cs}\t{ce}\t{gid}\t0\t{strand}\t{wname}\n")
            fa_fh.write(f">{gid}|{wname}|{chrom}:{cs}-{ce}({strand})|tss_off={tss_off}|"
                        f"len={len(seq)}|clipped={raw_len-len(seq)}\n")
            for i in range(0, len(seq), 60):
                fa_fh.write(seq[i:i + 60] + "\n")
            tab_fh.write(f"{gid}\t{tid}\t{chrom}\t{strand}\t{tss0}\t{wname}\t"
                         f"{cs}\t{ce}\t{len(seq)}\t{raw_len-len(seq)}\t{tss_off}\t{len(seq)}\t{t['source']}\n")
            n_ok += 1

    bed_fh.close()
    fa_fh.close()
    tab_fh.close()
    print(f"windows written: {n_ok}   dropped(edge): {n_drop}   "
          f"orphan models: {n_orphan} ({'dropped' if drop_orphans else 'kept'})   "
          f"genes skipped: {n_nogene}", flush=True)


if __name__ == "__main__":
    main()
