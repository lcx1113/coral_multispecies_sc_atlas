#!/usr/bin/env bash
set -euo pipefail


CONDA_SH="/home/users/Minjie/miniconda3/etc/profile.d/conda.sh"
CONDA_LCX="/home/users/Chenxi/anaconda3/etc/profile.d/conda.sh"
ENV_ETE="/home/users/Minjie/miniconda3/envs/ETE3"
ENV_R="/home/users/Chenxi/anaconda3/envs/r-env"

CDD_DB="/home/data/Chenxi/data/blastdb/cdd/Cdd"                 
RPSBPROC_DATA="/home/data/Chenxi/tools/RpsbProc-x64-linux/data" 
TREE_DOMAIN_R="/home/users/Ruijie/scRNAseq/cross_species/cell_type_analysis/script/tree_domain_visulize.R"

THREADS_RPS=12


if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <input.fasta>"
  exit 1
fi

INPUT="$1"
if [[ ! -f "$INPUT" ]]; then
  echo "Error: input fasta not found: $INPUT"
  exit 1
fi


if [[ "$INPUT" == *.fasta ]]; then
  BASE="${INPUT%.fasta}"
elif [[ "$INPUT" == *.fa ]]; then
  BASE="${INPUT%.fa}"
elif [[ "$INPUT" == *.faa ]]; then
  BASE="${INPUT%.faa}"
else
  BASE="${INPUT%.*}"
fi

echo "Input fasta: $INPUT"


export QT_QPA_PLATFORM='offscreen'
source "$CONDA_SH"
set +u
conda activate "$ENV_ETE"
set -u

ulimit -Sv 30000000
ete3 build -w standard_trimmed_fasttree \
  -a "$INPUT" \
  -o "${BASE}_tree"

rpsblast -db "$CDD_DB" \
  -query "$INPUT" \
  -evalue 0.01 -num_threads "$THREADS_RPS" \
  -outfmt 11 -out "${BASE}_cdd.asn"

rpsbproc -i "${BASE}_cdd.asn" \
  -o "${BASE}_cdd.out" \
  -m rep -d "$RPSBPROC_DATA" -t doms

perl -nle 'print $_ if /^QUERY/' "${BASE}_cdd.out" > "${BASE}.query"
perl -nle 'print $_ if /^\d+\tQuery/' "${BASE}_cdd.out" > "${BASE}.domain"


awk '/^>/ {if (seqlen){print seqname"\t"seqlen}; seqname=substr($1,2); seqlen=0; next} {seqlen+=length($0)} END{print seqname"\t"seqlen}' "$INPUT" > "${BASE}.lengths"
# ==========================================================

source "$CONDA_LCX"
set +u
conda activate "$ENV_R"
set -u

TREE_NW="${BASE}_tree/clustalo_default-trimal01-none-fasttree_full/$(basename "$INPUT").final_tree.nw"

if [[ ! -f "$TREE_NW" ]]; then
  echo "Error: tree file not found: $TREE_NW"
  echo "Please check ete3 output directory structure."
  exit 1
fi


Rscript "$TREE_DOMAIN_R" "$TREE_NW" "${BASE}.query" "${BASE}.domain" "${BASE}_tree.pdf" "${BASE}.lengths"

echo "Done. Output: ${BASE}_tree.pdf"