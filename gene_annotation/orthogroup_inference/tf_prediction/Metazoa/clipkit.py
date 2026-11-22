#!/usr/bin/env python3
"""
Parallel ClipKIT trimming for all alignments
"""

import os
import sys
import subprocess
from pathlib import Path
from multiprocessing import Pool
import re

# 配置
MAFFT_DIR = "mafft_result"
CLIPKIT_DIR = "clipkit_result"
GAP_THRESHOLD = 0.7
N_PROCESSES = 10

def parse_clipkit_log(log_file):
    """解析 ClipKIT 日志文件"""
    
    if not os.path.exists(log_file):
        return None
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # 提取统计信息
        stats = {}
        
        # 原始长度
        match = re.search(r'alignment length:\s+(\d+)', content)
        if match:
            stats['original_length'] = int(match.group(1))
        
        # 修剪后长度
        match = re.search(r'trimmed alignment length:\s+(\d+)', content)
        if match:
            stats['trimmed_length'] = int(match.group(1))
        
        # 保留比例
        if 'original_length' in stats and 'trimmed_length' in stats: 
            stats['kept_percent'] = (stats['trimmed_length'] / stats['original_length']) * 100
        
        return stats
    
    except Exception as e:
        return None

def trim_single_alignment(args):
    """处理单个比对文件"""
    
    family, aln_file, output_dir, gap_threshold = args
    
    family_output = Path(output_dir) / family
    family_output.mkdir(parents=True, exist_ok=True)
    
    cluster_name = Path(aln_file).stem
    output_file = family_output / f"{cluster_name}.clipkit.aln"
    log_file = family_output / f"{cluster_name}.clipkit.log"
    
    try:
        # 运行 ClipKIT
        cmd = [
            'clipkit',
            str(aln_file),
            '-m', 'kpic-gappy',
            '-g', str(gap_threshold),
            '-o', str(output_file),
            '-l'
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        # 解析日志
        stats = parse_clipkit_log(log_file)
        
        if stats:
            return (family, cluster_name, "success", stats)
        else:
            return (family, cluster_name, "success_no_stats", None)
    
    except Exception as e:
        return (family, cluster_name, f"failed: {e}", None)

def main():
    
    print("=" * 70)
    print("ClipKIT Alignment Trimming (Parallel)")
    print("=" * 70)
    print(f"Mode: kpic-gappy")
    print(f"Gap threshold: {GAP_THRESHOLD}")
    print(f"Parallel processes: {N_PROCESSES}")
    print("")
    
    # 收集所有比对文件
    print("Collecting alignment files...")
    tasks = []
    
    mafft_dir = Path(MAFFT_DIR)
    for family_dir in mafft_dir.iterdir():
        if not family_dir.is_dir():
            continue
        
        family = family_dir.name
        
        for aln_file in family_dir.glob("*.aln"):
            tasks.append((family, str(aln_file), CLIPKIT_DIR, GAP_THRESHOLD))
    
    print(f"Total alignments to trim: {len(tasks)}")
    print("")
    
    # 并行处理
    print("Processing...")
    
    with Pool(processes=N_PROCESSES) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(trim_single_alignment, tasks), 1):
            results.append(result)
            
            family, cluster_name, status, stats = result
            
            if i % 50 == 0:
                print(f"[{i}/{len(tasks)}] processed...")
            elif status == "success" and stats and i <= 10:
                kept = stats.get('kept_percent', 0)
                print(f"[{i}/{len(tasks)}] {family} {cluster_name}: {kept:.1f}% kept")
    
    print("")
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print("")
    
    # 统计
    family_stats = {}
    
    for family, cluster_name, status, stats in results:
        if family not in family_stats:
            family_stats[family] = {
                'success': 0,
                'failed': 0,
                'total_kept': 0,
                'count_with_stats': 0
            }
        
        if status.startswith("success"):
            family_stats[family]['success'] += 1
            
            if stats and 'kept_percent' in stats:
                family_stats[family]['total_kept'] += stats['kept_percent']
                family_stats[family]['count_with_stats'] += 1
        else:
            family_stats[family]['failed'] += 1
    
    # 打印报告
    print(f"{'TF Family':<30} {'Trimmed':<10} {'Avg Kept':<12} {'Failed':<10}")
    print("-" * 70)
    
    for family in sorted(family_stats.keys(), key=lambda f: family_stats[f]['success'], reverse=True):
        stats = family_stats[family]
        
        if stats['count_with_stats'] > 0:
            avg_kept = stats['total_kept'] / stats['count_with_stats']
            avg_kept_str = f"{avg_kept:.1f}%"
        else:
            avg_kept_str = "N/A"
        
        print(f"{family:<30} {stats['success']:<10} {avg_kept_str:<12} {stats['failed']:<10}")
    
    print("-" * 70)
    
    total_success = sum(s['success'] for s in family_stats.values())
    total_failed = sum(s['failed'] for s in family_stats.values())
    
    print(f"{'TOTAL':<30} {total_success:<10} {'':<12} {total_failed: <10}")
    print("")
    print(f"Results saved in: {CLIPKIT_DIR}/")

if __name__ == '__main__':
    main()