#!/usr/bin/env python3
"""
Parallel filtering of empty sequences from ClipKIT results
"""

import os
import sys
import subprocess
from pathlib import Path
from multiprocessing import Pool
from Bio import SeqIO
import re

# 配置
CLIPKIT_DIR = "clipkit_result"
FILTERED_DIR = "filtered_result"
N_PROCESSES = 10
MIN_NON_GAP_CHARS = 1  # 最少需要的非gap字符数

def filter_empty_sequences_single(input_file, output_file, min_non_gap_chars=1):
    """
    过滤单个文件中的空序列
    """
    try:
        valid_sequences = []
        removed_sequences = []
        
        for record in SeqIO.parse(input_file, "fasta"):
            # 计算非gap字符的数量
            sequence = str(record.seq).upper()
            non_gap_chars = sum(1 for char in sequence if char not in ['-', 'N', 'X', '?', '.'])
            
            if non_gap_chars >= min_non_gap_chars:
                valid_sequences.append(record)
            else:
                removed_sequences.append(record.id)
        
        # 写入过滤后的序列
        if valid_sequences:
            SeqIO.write(valid_sequences, output_file, "fasta")
            return {
                'status': 'success',
                'original_count': len(valid_sequences) + len(removed_sequences),
                'kept_count': len(valid_sequences),
                'removed_count': len(removed_sequences),
                'removed_ids': removed_sequences
            }
        else:
            # 如果所有序列都被过滤掉了，创建空文件但记录状态
            with open(output_file, 'w') as f:
                pass
            return {
                'status': 'all_empty',
                'original_count': len(removed_sequences),
                'kept_count': 0,
                'removed_count': len(removed_sequences),
                'removed_ids': removed_sequences
            }
    
    except Exception as e:
        return {
            'status': 'error',
            'error': str(e)
        }

def filter_single_file(args):
    """处理单个ClipKIT输出文件"""
    
    family, clipkit_file, output_dir, min_non_gap_chars = args
    
    family_output = Path(output_dir) / family
    family_output.mkdir(parents=True, exist_ok=True)
    
    cluster_name = Path(clipkit_file).stem.replace('.clipkit', '')  # 移除.clipkit后缀
    output_file = family_output / f"{cluster_name}.filtered.aln"
    
    try:
        result = filter_empty_sequences_single(clipkit_file, output_file, min_non_gap_chars)
        result['family'] = family
        result['cluster_name'] = cluster_name
        return result
    
    except Exception as e:
        return {
            'status': 'error',
            'family': family,
            'cluster_name': cluster_name,
            'error': str(e)
        }

def main():
    
    print("=" * 70)
    print("Empty Sequence Filtering (Parallel)")
    print("=" * 70)
    print(f"Input directory: {CLIPKIT_DIR}")
    print(f"Output directory: {FILTERED_DIR}")
    print(f"Minimum non-gap characters: {MIN_NON_GAP_CHARS}")
    print(f"Parallel processes: {N_PROCESSES}")
    print("")
    
    # 收集所有ClipKIT输出文件
    print("Collecting ClipKIT output files...")
    tasks = []
    
    clipkit_dir = Path(CLIPKIT_DIR)
    
    if not clipkit_dir.exists():
        print(f"Error: {CLIPKIT_DIR} directory not found!")
        sys.exit(1)
    
    for family_dir in clipkit_dir.iterdir():
        if not family_dir.is_dir():
            continue
        
        family = family_dir.name
        
        # 寻找 .clipkit.aln 文件
        for clipkit_file in family_dir.glob("*.clipkit.aln"):
            tasks.append((family, str(clipkit_file), FILTERED_DIR, MIN_NON_GAP_CHARS))
    
    if not tasks:
        print("No ClipKIT output files found!")
        print("Expected pattern: {CLIPKIT_DIR}/*/cluster_*.clipkit.aln")
        sys.exit(1)
    
    print(f"Total files to filter: {len(tasks)}")
    print("")
    
    # 并行处理
    print("Processing...")
    
    with Pool(processes=N_PROCESSES) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(filter_single_file, tasks), 1):
            results.append(result)
            
            family = result.get('family', 'unknown')
            cluster_name = result.get('cluster_name', 'unknown')
            status = result.get('status', 'unknown')
            
            if i % 50 == 0:
                print(f"[{i}/{len(tasks)}] processed...")
            elif status in ['success', 'all_empty'] and i <= 10:
                kept = result.get('kept_count', 0)
                removed = result.get('removed_count', 0)
                print(f"[{i}/{len(tasks)}] {family} {cluster_name}: kept {kept}, removed {removed}")
    
    print("")
    print("=" * 70)
    print("Filtering Summary")
    print("=" * 70)
    print("")
    
    # 统计
    family_stats = {}
    detailed_removals = []
    
    for result in results:
        family = result.get('family', 'unknown')
        cluster_name = result.get('cluster_name', 'unknown')
        status = result.get('status', 'unknown')
        
        if family not in family_stats:
            family_stats[family] = {
                'success': 0,
                'all_empty': 0,
                'error': 0,
                'total_kept': 0,
                'total_removed': 0,
                'total_files': 0
            }
        
        family_stats[family]['total_files'] += 1
        
        if status == 'success':
            family_stats[family]['success'] += 1
            family_stats[family]['total_kept'] += result.get('kept_count', 0)
            family_stats[family]['total_removed'] += result.get('removed_count', 0)
            
            # 记录被移除的序列
            if result.get('removed_count', 0) > 0:
                detailed_removals.append({
                    'family': family,
                    'cluster': cluster_name,
                    'removed_ids': result.get('removed_ids', [])
                })
        
        elif status == 'all_empty':
            family_stats[family]['all_empty'] += 1
            family_stats[family]['total_removed'] += result.get('removed_count', 0)
            
            detailed_removals.append({
                'family': family,
                'cluster': cluster_name,
                'removed_ids': result.get('removed_ids', [])
            })
        
        elif status == 'error':
            family_stats[family]['error'] += 1
    
    # 打印统计报告
    print(f"{'TF Family':<30} {'Success':<10} {'All Empty':<12} {'Errors':<10} {'Kept Seqs':<12} {'Removed Seqs':<12}")
    print("-" * 100)
    
    total_success = 0
    total_all_empty = 0
    total_errors = 0
    total_kept_seqs = 0
    total_removed_seqs = 0
    
    for family in sorted(family_stats.keys()):
        stats = family_stats[family]
        
        print(f"{family:<30} {stats['success']:<10} {stats['all_empty']:<12} {stats['error']:<10} "
              f"{stats['total_kept']:<12} {stats['total_removed']:<12}")
        
        total_success += stats['success']
        total_all_empty += stats['all_empty']
        total_errors += stats['error']
        total_kept_seqs += stats['total_kept']
        total_removed_seqs += stats['total_removed']
    
    print("-" * 100)
    print(f"{'TOTAL':<30} {total_success:<10} {total_all_empty:<12} {total_errors:<10} "
          f"{total_kept_seqs:<12} {total_removed_seqs:<12}")
    
    print("")
    print(f"Results saved in: {FILTERED_DIR}/")
    
    # 输出详细的移除信息
    if detailed_removals:
        print("")
        print("=" * 70)
        print("Removed Sequences Details")
        print("=" * 70)
        
        removal_summary_file = Path(FILTERED_DIR) / "removed_sequences_summary.txt"
        with open(removal_summary_file, 'w') as f:
            f.write("Removed Sequences Summary\n")
            f.write("=" * 50 + "\n\n")
            
            for removal in detailed_removals:
                family = removal['family']
                cluster = removal['cluster']
                removed_ids = removal['removed_ids']
                
                f.write(f"Family: {family}, Cluster: {cluster}\n")
                f.write(f"Removed {len(removed_ids)} sequences:\n")
                for seq_id in removed_ids:
                    f.write(f"  - {seq_id}\n")
                f.write("\n")
        
        print(f"Detailed removal log saved to: {removal_summary_file}")
        
        # 显示前几个移除案例
        print("\nSample removed sequences:")
        for i, removal in enumerate(detailed_removals[:5]):
            family = removal['family']
            cluster = removal['cluster']
            removed_ids = removal['removed_ids']
            print(f"  {family}/{cluster}: {len(removed_ids)} sequences removed")
            if removed_ids:
                print(f"    Example: {removed_ids[0]}")
        
        if len(detailed_removals) > 5:
            print(f"  ... and {len(detailed_removals) - 5} more clusters with removals")
    
    print("")
    print("Note: Files with all sequences removed are kept as empty files")
    print("      for consistency with the directory structure.")

if __name__ == '__main__':
    main()