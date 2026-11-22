#!/bin/bash

# 激活环境
source ~/anaconda3/etc/profile.d/conda.sh
conda activate tf_annotation

# 设置路径
SEQ_FILE="data/seq_Metazoa_20260425.fasta"
HMMER_DIR="hmmer_result"
DIAMOND_DIR="diamond_result"

mkdir -p "$DIAMOND_DIR"

echo "========================================"
echo "Diamond All-vs-All Alignment"
echo "========================================"
echo ""

# 索引序列文件
if [ ! -f "${SEQ_FILE}.ssi" ]; then
    echo "Indexing sequence file..."
    esl-sfetch --index "$SEQ_FILE"
fi

# 步骤 1: 提取每个 TF 家族的序列
echo "Step 1: Extracting sequences for each TF family"
echo "------------------------------------------------"

for hits_file in "$HMMER_DIR"/*_hits.txt; do
    
    if [ !  -f "$hits_file" ]; then
        continue
    fi
    
    family=$(basename "$hits_file" _hits.txt)
    
    # 提取序列 ID
    grep -v "^#" "$hits_file" | \
        awk '{print $1}' | \
        sort -u > "$DIAMOND_DIR/${family}_ids.txt"
    
    n_ids=$(wc -l < "$DIAMOND_DIR/${family}_ids.txt")
    
    if [ $n_ids -eq 0 ]; then
        echo "$family: 0 sequences (skipping)"
        continue
    fi
    
    # 提取序列
    esl-sfetch -f "$SEQ_FILE" "$DIAMOND_DIR/${family}_ids.txt" \
        > "$DIAMOND_DIR/${family}_sequences.fasta" 2>/dev/null
    
    echo "$family: $n_ids sequences extracted"
    
done

echo ""
echo "Step 2: Running Diamond all-vs-all alignments"
echo "------------------------------------------------"

# 步骤 2: Diamond 比对
for fasta in "$DIAMOND_DIR"/*_sequences.fasta; do
    
    if [ ! -f "$fasta" ]; then
        continue
    fi
    
    family=$(basename "$fasta" _sequences.fasta)
    n_seqs=$(grep -c "^>" "$fasta" 2>/dev/null || echo 0)
    
    echo ""
    echo "[$family]"
    echo "  Sequences: $n_seqs"
    
    if [ $n_seqs -lt 2 ]; then
        echo "  Skipping (need ≥2 sequences for alignment)"
        continue
    fi
    
    # 创建数据库
    echo "  Creating Diamond database..."
    diamond makedb \
        --in "$fasta" \
        --db "$DIAMOND_DIR/${family}_db" 2>&1 | grep -v "^diamond" || true
    
    # 运行比对
    echo "  Running Diamond blastp..."
    diamond blastp \
        --query "$fasta" \
        --db "$DIAMOND_DIR/${family}_db" \
        --out "$DIAMOND_DIR/${family}_diamond.txt" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --more-sensitive \
        --max-target-seqs 100 \
        --threads 8 \
        --evalue 1e-5 2>&1 | grep -v "^diamond" || true
    
    # 检查结果
    if [ -f "$DIAMOND_DIR/${family}_diamond.txt" ]; then
        n_alignments=$(wc -l < "$DIAMOND_DIR/${family}_diamond.txt")
        echo "  ✓ $n_alignments alignments"
        
        # 转换为 ABC 格式 (用于 MCL)
        awk '{print $1, $2, $12}' "$DIAMOND_DIR/${family}_diamond.txt" \
            > "$DIAMOND_DIR/${family}_diamond.abc"
        
        # 统计比对覆盖率
        n_unique_queries=$(awk '{print $1}' "$DIAMOND_DIR/${family}_diamond.txt" | sort -u | wc -l)
        coverage=$(awk "BEGIN {printf \"%.1f\", $n_unique_queries/$n_seqs*100}")
        echo "  Coverage: $coverage% ($n_unique_queries/$n_seqs sequences with hits)"
    else
        echo "  ✗ Diamond failed"
    fi
    
done

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo ""
echo "TF families processed:"
for abc in "$DIAMOND_DIR"/*.abc; do
    if [ -f "$abc" ]; then
        family=$(basename "$abc" _diamond.abc)
        n_seqs=$(wc -l < "$DIAMOND_DIR/${family}_ids.txt")
        n_alignments=$(wc -l < "$abc")
        printf "%-30s %6d seqs, %8d alignments\n" "$family" "$n_seqs" "$n_alignments"
    fi
done | sort -k3 -rn

echo ""
echo "Results saved in:  $DIAMOND_DIR/"
echo "  *_sequences.fasta  - Protein sequences"
echo "  *_diamond.txt      - Diamond output (tabular)"
echo "  *_diamond.abc      - ABC format for MCL clustering"
echo ""
echo "Next step: MCL clustering"
echo "  Example:  mcl ${DIAMOND_DIR}/bZIP_diamond.abc --abc -I 1. 2 -o ${DIAMOND_DIR}/bZIP_clusters.txt"