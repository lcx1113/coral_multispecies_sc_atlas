#!/bin/bash

# 激活环境
source ~/anaconda3/etc/profile.d/conda.sh
conda activate tf_annotation

# 设置路径
HMM_DIR="/absolute/path/hmms"
SEQ_FILE="data/seq_Metazoa_20260425.fasta"
OUT_DIR="hmmer_result"
CSV_FILE="tf_parameter.csv"

mkdir -p $OUT_DIR
esl-sfetch --index $SEQ_FILE 2>/dev/null

echo "Processing TF families..."
echo ""

# 使用 Python 正确解析 CSV
python3 << 'PYTHON_SCRIPT'
import csv
import subprocess
import os

HMM_DIR = "/absolute/path/hmms"
SEQ_FILE = "data/seq_Metazoa_20260425.fasta"
OUT_DIR = "hmmer_result"
CSV_FILE = "tf_parameter.csv"

with open(CSV_FILE, 'r') as f:
    reader = csv.DictReader(f)
    
    for row in reader:
        family = row['TF family']. strip()
        domains_str = row['Pfam domains'].strip()
        inflation = row['MCL inflation parameter'].strip()
        threshold = row['HMMsearch threshold parameter'].strip()
        
        print(f"[{family}]")
        print(f"  Domains: {domains_str}")
        
        # 分割 domains
        domains = [d.strip() for d in domains_str.split(',')]
        print(f"  Parsed:  {len(domains)} domains - {', '.join(domains)}")
        
        # 收集 HMM 文件
        hmm_files = []
        missing = []
        
        for domain in domains:
            hmm_path = f"{HMM_DIR}/{domain}.hmm"
            if os.path.exists(hmm_path):
                hmm_files.append(hmm_path)
            else:
                missing.append(domain)
        
        print(f"  Found:  {len(hmm_files)}/{len(domains)} HMM files")
        if missing:
            print(f"  Missing: {', '.join(missing)}")
        
        if not hmm_files:
            print("  ✗ No HMM files found, skipping\n")
            continue
        
        # 合并 HMM 文件
        combined_hmm = f"{OUT_DIR}/{family}_combined.hmm"
        with open(combined_hmm, 'w') as out:
            for hmm_file in hmm_files:
                with open(hmm_file, 'r') as f:
                    out.write(f.read())
        
        print(f"  Combined HMM size: {os.path.getsize(combined_hmm)} bytes")
        
        # 运行 hmmsearch
        cmd = [
            'hmmsearch',
            '--cut_ga',
            '--cpu', '8',
            '--tblout', f"{OUT_DIR}/{family}_hits.txt",
            '--domtblout', f"{OUT_DIR}/{family}_domains.txt",
            '-o', f"{OUT_DIR}/{family}_search.log",
            combined_hmm,
            SEQ_FILE
        ]
        
        print(f"  Running hmmsearch...")
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            
            # 统计结果
            hits_file = f"{OUT_DIR}/{family}_hits.txt"
            if os.path.exists(hits_file):
                with open(hits_file, 'r') as f:
                    hits = [l for l in f if not l.startswith('#')]
                    nhits = len(hits)
                
                print(f"  ✓ {nhits} hits")
                
                if nhits > 0:
                    # 按物种统计
                    species_counts = {}
                    for hit in hits:
                        gene = hit.split()[0]
                        species = gene.split('_')[0]
                        species_counts[species] = species_counts.get(species, 0) + 1
                    
                    print("  Top species:")
                    for sp, count in sorted(species_counts.items(), key=lambda x: x[1], reverse=True)[:5]:
                        print(f"    {sp}: {count}")
            else:
                print("  ✗ No results file")
        except subprocess.CalledProcessError as e:
            print(f"  ✗ hmmsearch failed: {e}")
        
        print("")

print("\n=== Summary ===")
# 统计所有结果
for fname in os.listdir(OUT_DIR):
    if fname.endswith('_hits.txt'):
        with open(f"{OUT_DIR}/{fname}", 'r') as f:
            hits = len([l for l in f if not l.startswith('#')])
            if hits > 0:
                family = fname.replace('_hits.txt', '')
                print(f"{family}: {hits} hits")

PYTHON_SCRIPT

echo ""
echo "Done!  Results in $OUT_DIR/"