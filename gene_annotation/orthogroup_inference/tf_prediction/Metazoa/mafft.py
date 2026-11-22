#!/usr/bin/env python3
"""
Robust MAFFT pipeline for all TF families
"""

import os
import sys
import subprocess
from pathlib import Path
from Bio import SeqIO
from multiprocessing import Pool
import time

# 配置
SEQ_FILE = "data/seq_Metazoa_20260425.fasta"
MCL_DIR = "mcl_result"
MAFFT_DIR = "mafft_result"
N_PROCESSES = 8  # 并行进程数

def load_sequences(seq_file):
    """加载序列到内存（加速）"""
    print(f"Loading sequences from {seq_file}...")
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    print(f"  Loaded {len(seq_dict)} sequences")
    return seq_dict

def process_single_cluster(args):
    """处理单个 cluster"""
    
    family, cluster_num, gene_ids, seq_dict, output_dir = args
    
    family_dir = Path(output_dir) / family
    family_dir.mkdir(parents=True, exist_ok=True)
    
    n_genes = len(gene_ids)
    
    # 跳过单基因
    if n_genes < 2:
        return (family, cluster_num, "singleton", 0, 0)
    
    try:
        # 提取序列
        fasta_file = family_dir / f"cluster_{cluster_num}.fasta"
        
        extracted = []
        for gene_id in gene_ids: 
            if gene_id in seq_dict: 
                extracted.append(seq_dict[gene_id])
        
        if len(extracted) < 2:
            return (family, cluster_num, "failed_extraction", len(extracted), 0)
        
        # 保存序列
        SeqIO.write(extracted, fasta_file, "fasta")
        
        # 运行 MAFFT
        aln_file = family_dir / f"cluster_{cluster_num}.aln"
        
        cmd = [
            'mafft',
            '--genafpair',
            '--maxiterate', '10000',
            '--quiet',
            '--thread', '8',
            str(fasta_file)
        ]
        
        with open(aln_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.DEVNULL, check=True)
        
        # 检查结果
        with open(aln_file, 'r') as f:
            n_aligned = sum(1 for line in f if line.startswith('>'))
        
        if n_aligned >= 2:
            # 获取比对长度
            seqs = list(SeqIO.parse(aln_file, "fasta"))
            aln_len = len(seqs[0].  seq) if seqs else 0
            
            return (family, cluster_num, "success", n_aligned, aln_len)
        else:
            return (family, cluster_num, "failed_alignment", 0, 0)
    
    except Exception as e: 
        return (family, cluster_num, f"error: {e}", 0, 0)

def main():
    
    print("=" * 70)
    print("MAFFT Alignment Pipeline (Python)")
    print("=" * 70)
    print("")
    
    # 加载序列
    seq_dict = load_sequences(SEQ_FILE)
    print("")
    
    # 收集所有任务
    print("Collecting tasks...")
    tasks = []
    
    mcl_dir = Path(MCL_DIR)
    for cluster_file in sorted(mcl_dir.glob("*_clusters.txt")):
        
        family = cluster_file.stem.replace('_clusters', '')
        
        with open(cluster_file, 'r') as f:
            for cluster_num, line in enumerate(f, 1):
                # 正确分割基因 ID
                gene_ids = []
                for part in line.replace('\t', ' ').split():
                    part = part.strip()
                    if part: 
                        gene_ids.append(part)
                
                if gene_ids:
                    tasks.append((family, cluster_num, gene_ids, seq_dict, MAFFT_DIR))
    
    print(f"Total tasks: {len(tasks)}")
    print(f"Running with {N_PROCESSES} parallel processes...")
    print("")
    
    # 并行处理
    start_time = time.time()
    
    with Pool(processes=N_PROCESSES) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(process_single_cluster, tasks), 1):
            results.append(result)
            
            family, cluster_num, status, n_seqs, aln_len = result
            
            if i % 10 == 0 or status == "success":
                if status == "success":
                    print(f"[{i}/{len(tasks)}] {family} cluster {cluster_num}:  ✓ {n_seqs} seqs, len={aln_len}")
                elif i % 50 == 0:
                    print(f"[{i}/{len(tasks)}] processed...")
    
    elapsed = time.time() - start_time
    
    # 统计
    print("")
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print("")
    
    family_stats = {}
    for family, cluster_num, status, n_seqs, aln_len in results: 
        if family not in family_stats:
            family_stats[family] = {'success': 0, 'singleton': 0, 'failed':  0}
        
        if status == "success":
            family_stats[family]['success'] += 1
        elif status == "singleton":
            family_stats[family]['singleton'] += 1
        else:
            family_stats[family]['failed'] += 1
    
    # 打印报告
    print(f"{'TF Family':<30} {'Alignments':<12} {'Singletons':<12} {'Failed':<12}")
    print("-" * 70)
    
    for family in sorted(family_stats.keys(), key=lambda f: family_stats[f]['success'], reverse=True):
        stats = family_stats[family]
        print(f"{family:<30} {stats['success']:<12} {stats['singleton']:<12} {stats['failed']:<12}")
    
    print("-" * 70)
    
    total_success = sum(s['success'] for s in family_stats.values())
    total_singleton = sum(s['singleton'] for s in family_stats.values())
    total_failed = sum(s['failed'] for s in family_stats.values())
    
    print(f"{'TOTAL':<30} {total_success:<12} {total_singleton:<12} {total_failed: <12}")
    print("")
    print(f"Time elapsed: {elapsed/60:.1f} minutes")
    print(f"Results saved in: {MAFFT_DIR}/")

if __name__ == '__main__':
    main()