#!/usr/bin/bash
# ============================================
# Split-seq + STAR + DropEst Pipeline
# 作者: Ruijie,Chenxi
# ============================================

# 参数：
# $1 = R1.fastq.gz
# $2 = R2.fastq.gz
# $3 = STAR index 路径
# $4 = GTF 文件路径
# $5 = (可选)线程数，默认12
# $6 = (可选)输出文件前缀，默认使用R2文件名

set -e  # 出错即停止
set -o pipefail

# ---------- 参数输入 ----------
FASTQ_R1=$1
FASTQ_R2=$2
STAR_INDEX=$3
GTF_FILE=$4
THREADS=${5:-12}
CUSTOM_PREFIX=$6

# ---------- 检查输入 ----------
if [[ -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$STAR_INDEX" || -z "$GTF_FILE" ]]; then
    echo "用法: $0 R1.fastq.gz R2.fastq.gz STAR_index_path GTF_file [threads]"
    exit 1
fi

# ---------- 常用路径 ----------
DROP_BIN=/home/users/Ruijie/pipeline/dropEst/bin
CONFIG_XML=/home/users/Ruijie/scRNAseq/wjs_data/config/split_seq.xml
source /home/users/Chenxi/anaconda3/etc/profile.d/conda.sh
conda activate /home/users/Chenxi/anaconda3/envs/genome_tools

# ---------- 基于 R2 文件名生成前缀 ----------
if [[ -n "$CUSTOM_PREFIX" ]]; then
    R2_PREFIX=$CUSTOM_PREFIX
else
    R2_BASENAME=$(basename "$FASTQ_R2" .gz)
    R2_PREFIX=${R2_BASENAME%.fastq*}  # 去掉.fastq或.fastq.gz
fi
echo "样本标识：$R2_PREFIX"
echo "-------------------------------------------"

# ---------- Step 1. droptag ----------
echo "[1/6] droptag..."
$DROP_BIN/droptag -c $CONFIG_XML -S -s -p $THREADS -n "${R2_PREFIX}.tagged" "$FASTQ_R1" "$FASTQ_R2"

# 输出: ${R2_PREFIX}.tagged.1.fastq.gz 和 ${R2_PREFIX}.tagged.params.gz

# ---------- Step 2. 去掉 polyA ----------
echo "[2/6] 去掉 polyA..."
cutadapt -j $THREADS -a AAAAAAAAAA \
  -o no_A_${R2_PREFIX}.tagged.1.fastq.gz \
  ${R2_PREFIX}.tagged.1.fastq.gz

# ---------- Step 3. 去掉 polyG ----------
echo "[3/6] 去掉 polyG..."
cutadapt -j $THREADS -a GGGGGGGGGG \
  -o no_AG_${R2_PREFIX}.tagged.1.fastq.gz \
  no_A_${R2_PREFIX}.tagged.1.fastq.gz

# ---------- Step 4. 清理空序列 ----------
echo "[4/6] 清理空序列..."
seqkit seq -m 15 -g ./no_AG_${R2_PREFIX}.tagged.1.fastq.gz \
  -o clean_no_AG_${R2_PREFIX}.tagged.1.fastq.gz

# ---------- Step 5. STAR 比对 ----------
echo "[5/6] STAR 比对..."
ulimit -n 65535
mkdir -p alignment && cd alignment

/home/users/Chenxi/.local/share/mamba/envs/extend3p/bin/STAR-avx2 --outSAMtype BAM SortedByCoordinate \
     --runThreadN $THREADS \
     --genomeDir $STAR_INDEX \
     --readFilesCommand zcat \
     --readFilesIn ../clean_no_AG_${R2_PREFIX}.tagged.1.fastq.gz \
     --outFileNamePrefix clean_no_AG_${R2_PREFIX}.tagged.1. 

cd ..
rm no_AG_${R2_PREFIX}.tagged.1.fastq.gz
rm no_A_${R2_PREFIX}.tagged.1.fastq.gz

# ---------- Step 6. DropEst ----------
echo "[6/6] DropEst 计数..."
mkdir -p est_out && cd est_out
$DROP_BIN/dropest -V -m \
  -r ../${R2_PREFIX}.tagged.params.gz \
  -g $GTF_FILE \
  -c $CONFIG_XML \
  ../alignment/clean_no_AG_${R2_PREFIX}.tagged.1.Aligned.sortedByCoord.out.bam

# ---------- Step 7. DropReport ----------
echo "[Final] 生成 DropReport..."
$DROP_BIN/dropReport.Rsc cell.counts.rds

echo "处理完成：$R2_PREFIX"
