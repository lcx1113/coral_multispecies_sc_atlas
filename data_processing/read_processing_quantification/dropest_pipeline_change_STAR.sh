#!/usr/bin/bash
# ============================================
# Split-seq + STAR + DropEst Pipeline
# ============================================

# Arguments:
# $1 = R1.fastq.gz
# $2 = R2.fastq.gz
# $3 = STAR index directory
# $4 = GTF file path
# $5 = (optional) number of threads; default: 12
# $6 = (optional) output prefix; defaults to the R2 filename-derived prefix

set -e  # Exit on error
set -o pipefail

# ---------- Input arguments ----------
FASTQ_R1=$1
FASTQ_R2=$2
STAR_INDEX=$3
GTF_FILE=$4
THREADS=${5:-12}
CUSTOM_PREFIX=$6

# ---------- Check required inputs ----------
if [[ -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$STAR_INDEX" || -z "$GTF_FILE" ]]; then
    echo "用法: $0 R1.fastq.gz R2.fastq.gz STAR_index_path GTF_file [threads]"
    exit 1
fi

# ---------- Software and configuration paths ----------
DROP_BIN=/absolute/path/dropEst/bin
# Set the software paths below; the bundled XML is resolved before changing directories.
CONFIG_XML="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/config/split_seq.xml"
source /absolute/path/conda/etc/profile.d/conda.sh
conda activate /absolute/path/conda/envs/genome_tools

# ---------- Use the supplied prefix or derive it from the R2 filename ----------
if [[ -n "$CUSTOM_PREFIX" ]]; then
    R2_PREFIX=$CUSTOM_PREFIX
else
    R2_BASENAME=$(basename "$FASTQ_R2" .gz)
    R2_PREFIX=${R2_BASENAME%.fastq*}  # Remove the .fastq suffix after stripping .gz above
fi
echo "样本标识：$R2_PREFIX"
echo "-------------------------------------------"

# ---------- Step 1. droptag ----------
echo "[1/6] droptag..."
$DROP_BIN/droptag -c $CONFIG_XML -S -s -p $THREADS -n "${R2_PREFIX}.tagged" "$FASTQ_R1" "$FASTQ_R2"

# Outputs: ${R2_PREFIX}.tagged.1.fastq.gz and ${R2_PREFIX}.tagged.params.gz

# ---------- Step 2. Trim poly(A) ----------
echo "[2/6] 去掉 polyA..."
cutadapt -j $THREADS -a AAAAAAAAAA \
  -o no_A_${R2_PREFIX}.tagged.1.fastq.gz \
  ${R2_PREFIX}.tagged.1.fastq.gz

# ---------- Step 3. Trim poly(G) ----------
echo "[3/6] 去掉 polyG..."
cutadapt -j $THREADS -a GGGGGGGGGG \
  -o no_AG_${R2_PREFIX}.tagged.1.fastq.gz \
  no_A_${R2_PREFIX}.tagged.1.fastq.gz

# ---------- Step 4. Remove gaps and filter reads by minimum length ----------
echo "[4/6] 清理空序列..."
seqkit seq -m 15 -g ./no_AG_${R2_PREFIX}.tagged.1.fastq.gz \
  -o clean_no_AG_${R2_PREFIX}.tagged.1.fastq.gz

# ---------- Step 5. Align reads with STAR ----------
echo "[5/6] STAR 比对..."
ulimit -n 65535
mkdir -p alignment && cd alignment

/absolute/path/STAR --outSAMtype BAM SortedByCoordinate \
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
