#!/bin/bash

# 激活环境
source ~/anaconda3/etc/profile.d/conda.sh
conda activate tf_annotation

# 设置路径
DIAMOND_DIR="diamond_result"
MCL_DIR="mcl_result"
CSV_FILE="tf_parameter.csv"

mkdir -p $MCL_DIR

echo "========================================"
echo "MCL Clustering"
echo "========================================"
echo ""

# 使用 Python 读取 CSV 获取每个 TF 的 inflation 参数
python3 << 'PYTHON_SCRIPT'
import csv
import subprocess
import os

DIAMOND_DIR = "diamond_result"
MCL_DIR = "mcl_result"
CSV_FILE = "tf_parameter.csv"

# 读取 inflation 参数
inflation_map = {}
with open(CSV_FILE, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        family = row['TF family']. strip()
        inflation = float(row['MCL inflation parameter']. strip())
        inflation_map[family] = inflation

print(f"Loaded inflation parameters for {len(inflation_map)} TF families\n")

# 对每个有 Diamond 结果的 TF 家族运行 MCL
for abc_file in sorted(os.listdir(DIAMOND_DIR)):
    if not abc_file.endswith('_diamond.abc'):
        continue
    
    family = abc_file.replace('_diamond.abc', '')
    abc_path = f"{DIAMOND_DIR}/{abc_file}"
    
    # 获取 inflation 参数
    inflation = inflation_map. get(family, 1.1)  # 默认 1.1
    
    print(f"[{family}]")
    print(f"  Inflation: {inflation}")
    
    # 检查 ABC 文件
    if not os.path.exists(abc_path):
        print(f"  ✗ ABC file not found")
        continue
    
    # 检查文件是否为空
    if os.path.getsize(abc_path) == 0:
        print(f"  ✗ ABC file is empty")
        continue
    
    # 统计输入
    with open(abc_path, 'r') as f:
        n_edges = len(f.readlines())
    
    # 统计唯一节点数
    nodes = set()
    with open(abc_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                nodes.add(parts[0])
                nodes.add(parts[1])
    n_nodes = len(nodes)
    
    print(f"  Input:  {n_nodes} genes, {n_edges} edges")
    
    # 运行 MCL
    output_file = f"{MCL_DIR}/{family}_clusters.txt"
    cmd = [
        'mcl',
        abc_path,
        '--abc',
        '-I', str(inflation),
        '-o', output_file
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # 统计输出
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                clusters = [line.strip().split() for line in f if line.strip()]
            
            n_clusters = len(clusters)
            cluster_sizes = [len(c) for c in clusters]
            
            print(f"  ✓ {n_clusters} clusters")
            print(f"    Sizes: min={min(cluster_sizes)}, max={max(cluster_sizes)}, mean={sum(cluster_sizes)/len(cluster_sizes):.1f}")
            
            # 显示最大的几个 clusters
            if n_clusters > 0:
                top_clusters = sorted(cluster_sizes, reverse=True)[:3]
                print(f"    Top 3: {', '.join(map(str, top_clusters))}")
        else:
            print(f"  ✗ Output file not created")
    
    except subprocess.CalledProcessError as e:
        print(f"  ✗ MCL failed: {e}")
    
    print("")

print("=" * 50)
print("Summary")
print("=" * 50)
print("")

# 汇总所有结果
summary = []
for cluster_file in sorted(os.listdir(MCL_DIR)):
    if cluster_file.endswith('_clusters.txt'):
        family = cluster_file.replace('_clusters.txt', '')
        path = f"{MCL_DIR}/{cluster_file}"
        
        with open(path, 'r') as f:
            clusters = [line.strip().split() for line in f if line.strip()]
        
        n_clusters = len(clusters)
        n_genes = sum(len(c) for c in clusters)
        
        summary.append((family, n_clusters, n_genes))

# 按 cluster 数排序
summary.sort(key=lambda x: x[1], reverse=True)

print(f"{'TF Family':<30} {'Clusters':<10} {'Total Genes':<15}")
print("-" * 60)
for family, n_clusters, n_genes in summary:
    print(f"{family:<30} {n_clusters:<10} {n_genes:<15}")

PYTHON_SCRIPT

echo ""
echo "Done!  Results in $MCL_DIR/"