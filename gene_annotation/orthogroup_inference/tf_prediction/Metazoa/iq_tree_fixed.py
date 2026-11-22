#!/usr/bin/env python3
"""
Robust IQ-TREE phylogenetic inference with better error handling
"""

import os
import sys
import subprocess
from pathlib import Path
from multiprocessing import Pool
import time
import re

# 配置
CLIPKIT_DIR = "filtered_result"
IQTREE_DIR = "iqtree_result"
N_PROCESSES = 6
N_THREADS_PER_TREE = 12

def check_alignment_quality(aln_file):
    """检查比对质量，决定是否值得建树"""
    
    try:
        n_seqs = 0
        total_len = 0
        
        with open(aln_file, 'r') as f:
            seq = ""
            for line in f: 
                if line.startswith('>'):
                    if seq:
                        total_len = len(seq)
                    n_seqs += 1
                    seq = ""
                else:
                    seq += line.strip()
            
            if seq:
                total_len = len(seq)
        
        # 质量检查
        if n_seqs < 4:
            return False, f"too_few_seqs_{n_seqs}"
        
        if total_len < 50:
            return False, f"too_short_{total_len}bp"
        
        return True, "ok"
    
    except Exception as e: 
        return False, f"read_error:  {e}"

def build_tree(args):
    """为单个比对构建系统发育树"""
    
    family, aln_file, output_dir, n_threads = args
    
    family_output = Path(output_dir) / family
    family_output.mkdir(parents=True, exist_ok=True)
    
    cluster_name = Path(aln_file).stem.replace('.clipkit', '')
    prefix = str(family_output / cluster_name)
    
    treefile = Path(f"{prefix}.treefile")
    log_file = Path(f"{prefix}.log")
    
    # ===== 步骤 1: 检查是否已有有效的树 =====
    if treefile.exists() and treefile.stat().st_size > 0:
        # 验证树文件格式
        try:
            with open(treefile, 'r') as f:
                tree_content = f.read().strip()
                if tree_content and '(' in tree_content and ';' in tree_content:
                    # 检查 log 是否完整
                    if log_file.exists():
                        with open(log_file, 'r') as lf:
                            log_content = lf.read()
                            if "Total wall-clock time" in log_content:
                                stats = parse_iqtree_log(log_file)
                                return (family, cluster_name, "success", stats)
                            else:
                                return (family, cluster_name, "success_incomplete_log", None)
                    else:
                        return (family, cluster_name, "success_no_log", None)
        except Exception as e:
            pass  # 继续尝试重建
    
    # ===== 步骤 2: 检查比对质量 =====
    quality_ok, quality_msg = check_alignment_quality(aln_file)
    
    if not quality_ok: 
        return (family, cluster_name, f"skipped_{quality_msg}", None)
    
    # ===== 步骤 3: 运行 IQ-TREE =====
    try:
        cmd = [
            'iqtree3',
            '-s', str(aln_file),
            '-m', 'MFP',
            '-mset', 'LG,WAG,JTT',
            '-bb', '1000',
            '-bnni',
            '-nt', str(n_threads),
            '-pre', prefix,
            '-nm', '10000',
            '-nstop', '200',
            '-me', '0.001',
            '-quiet'
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=36000  # 1小时
        )
        
        # ===== 步骤 4: 处理结果 =====
        
        # 4a. 检查返回码
        if result.returncode != 0:
            # IQ-TREE 报错，但可能已生成树
            if treefile.exists() and treefile.stat().st_size > 0:
                return (family, cluster_name, f"success_with_error_{result.returncode}", None)
            else:
                # 尝试从 stderr 获取错误信息
                error_msg = result.stderr[: 200] if result.stderr else "unknown"
                return (family, cluster_name, f"failed_rc{result.returncode}:  {error_msg}", None)
        
        # 4b.  正常完成，检查输出
        if treefile.exists() and treefile.stat().st_size > 0:
            stats = parse_iqtree_log(log_file) if log_file.exists() else None
            return (family, cluster_name, "success", stats)
        else:
            return (family, cluster_name, "no_output", None)
    
    # ===== 步骤 5: 异常处理 =====
    
    except subprocess.TimeoutExpired:
        # 超时，检查是否已生成树
        if treefile.exists() and treefile.stat().st_size > 0:
            return (family, cluster_name, "success_timeout", None)
        else:
            return (family, cluster_name, "timeout_no_tree", None)
    
    except Exception as e:
        # 其他异常
        if treefile.exists() and treefile.stat().st_size > 0:
            return (family, cluster_name, f"success_exception:  {type(e).__name__}", None)
        else:
            return (family, cluster_name, f"error:  {type(e).__name__}: {str(e)[:100]}", None)

def parse_iqtree_log(log_file):
    """解析 IQ-TREE 日志文件"""
    
    if not os.path.exists(log_file):
        return None
    
    stats = {}
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # 提取最佳模型
        match = re.search(r'Best-fit model.*?:\s+(\S+)', content)
        if match:
            stats['model'] = match.group(1)
        
        # 提取对数似然
        match = re.search(r'Log-likelihood of the tree:\s+([-\d.]+)', content)
        if match:
            stats['log_likelihood'] = float(match.group(1))
        
        # 提取序列数
        match = re.search(r'Number of sequences:\s+(\d+)', content)
        if match:
            stats['n_seqs'] = int(match.group(1))
        
        # 提取位点数
        match = re.search(r'Number of positions:\s+(\d+)', content)
        if match:
            stats['total_sites'] = int(match.group(1))
        
        return stats
    
    except Exception: 
        return None

def main():
    
    print("=" * 70)
    print("IQ-TREE Phylogenetic Inference (Robust Version)")
    print("=" * 70)
    print(f"Parallel processes: {N_PROCESSES}")
    print(f"Threads per tree: {N_THREADS_PER_TREE}")
    print("")
    
    # 收集所有比对文件
    print("Collecting alignment files...")
    tasks = []
    
    clipkit_dir = Path(CLIPKIT_DIR)
    for family_dir in clipkit_dir.iterdir():
        if not family_dir.is_dir():
            continue
        
        family = family_dir.name
        
        for aln_file in family_dir.glob("*.filtered.aln"):
            tasks.append((family, str(aln_file), IQTREE_DIR, N_THREADS_PER_TREE))
    
    print(f"Total trees to build: {len(tasks)}")
    print("")
    
    # 并行处理
    print("Building trees...")
    print("")
    
    start_time = time.time()
    
    with Pool(processes=N_PROCESSES) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(build_tree, tasks), 1):
            results.append(result)
            
            family, cluster_name, status, stats = result
            
            if status.startswith("success"):
                if i <= 20 or i % 50 == 0:  # 只显示前20个和每50个
                    model = stats.get('model', 'N/A') if stats else 'N/A'
                    print(f"[{i}/{len(tasks)}] {family} {cluster_name}:  ✓ {status} (model: {model})")
            elif status.startswith("skipped"):
                if i <= 10:   # 只显示前10个跳过的
                    print(f"[{i}/{len(tasks)}] {family} {cluster_name}:  - {status}")
            else:
                print(f"[{i}/{len(tasks)}] {family} {cluster_name}:  ✗ {status}")
    
    elapsed = time.time() - start_time
    
    # ===== 详细统计 =====
    print("")
    print("=" * 70)
    print("Detailed Summary")
    print("=" * 70)
    print("")
    
    # 按状态分类
    status_counts = {}
    family_stats = {}
    model_counts = {}
    
    for family, cluster_name, status, stats in results: 
        # 状态统计
        status_category = status.split('_')[0]  # 取第一部分（success, failed, skipped 等）
        status_counts[status] = status_counts.get(status, 0) + 1
        
        # 家族统计
        if family not in family_stats: 
            family_stats[family] = {
                'success': 0,
                'success_partial': 0,
                'failed': 0,
                'skipped': 0
            }
        
        if status == "success":
            family_stats[family]['success'] += 1
            # 模型统计
            if stats and 'model' in stats: 
                model = stats['model'].split('+')[0]
                model_counts[model] = model_counts.get(model, 0) + 1
        elif status.startswith("success"):
            family_stats[family]['success_partial'] += 1
        elif status.startswith("skipped"):
            family_stats[family]['skipped'] += 1
        else:
            family_stats[family]['failed'] += 1
    
    # 打印状态统计
    print("Status breakdown:")
    for status in sorted(status_counts.keys()):
        count = status_counts[status]
        percent = (count / len(results)) * 100
        print(f"  {status: <40} {count:>6} ({percent: >5.1f}%)")
    
    print("")
    print("-" * 70)
    
    # 打印家族统计
    print(f"{'TF Family':<30} {'Complete':<10} {'Partial':<10} {'Skipped':<10} {'Failed':<10}")
    print("-" * 70)
    
    for family in sorted(family_stats.keys(), 
                        key=lambda f: family_stats[f]['success'] + family_stats[f]['success_partial'], 
                        reverse=True):
        stats = family_stats[family]
        print(f"{family:<30} {stats['success']:<10} {stats['success_partial']:<10} "
              f"{stats['skipped']:<10} {stats['failed']:<10}")
    
    print("-" * 70)
    
    total_success = sum(s['success'] for s in family_stats.values())
    total_partial = sum(s['success_partial'] for s in family_stats.values())
    total_skipped = sum(s['skipped'] for s in family_stats.values())
    total_failed = sum(s['failed'] for s in family_stats.values())
    
    print(f"{'TOTAL':<30} {total_success:<10} {total_partial:<10} {total_skipped:<10} {total_failed: <10}")
    print("")
    
    # 打印模型统计
    if model_counts:
        print("Best-fit model distribution:")
        for model, count in sorted(model_counts.items(), key=lambda x: x[1], reverse=True):
            percent = (count / total_success) * 100 if total_success > 0 else 0
            print(f"  {model}:  {count} ({percent:.1f}%)")
        print("")
    
    print(f"Time elapsed: {elapsed/60:.1f} minutes")
    if len(tasks) > 0:
        print(f"Average time per tree: {elapsed/len(tasks):.1f} seconds")
    print("")
    print(f"Results saved in: {IQTREE_DIR}/")
    print("")
    print("Legend:")
    print("  Complete:   树和日志都完整")
    print("  Partial:   树存在但日志不完整（仍可使用）")
    print("  Skipped:   比对质量不足，未建树")
    print("  Failed:    IQ-TREE 失败且未生成树")

if __name__ == '__main__':
    main()