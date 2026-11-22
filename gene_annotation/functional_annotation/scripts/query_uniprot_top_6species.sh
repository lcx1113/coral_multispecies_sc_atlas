#!/bin/bash

# 帮助信息
usage() {
  echo "Usage: $0 -q QUERY [-d DB] -p PREFIX -o OUT_DIR [-t THREADS]"
  echo "  -q QUERY      输入的蛋白查询文件路径"
  echo "  -d DB         BLAST数据库路径（可选，默认: /absolute/path/blastdb/swiss_prot）"
  echo "  -p PREFIX     输出文件前缀 (如 Adig)"
  echo "  -o OUT_DIR    输出文件（OUT和TOPHIT）的存储目录"
  echo "  -t THREADS    BLAST线程数，默认1"
  exit 1
}

# 默认参数
THREADS=1
DB_DEFAULT="/absolute/path/blastdb/swiss_prot"

# 解析参数
while getopts "q:d:p:o:t:" opt; do
  case $opt in
    q) QUERY="$OPTARG" ;;
    d) DB="$OPTARG" ;;
    p) PREFIX="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done

# 使用默认DB路径（如果未指定）
if [[ -z "$DB" ]]; then
  DB="$DB_DEFAULT"
fi

# 检查其它必需参数
if [[ -z "$QUERY" || -z "$PREFIX" || -z "$OUT_DIR" ]]; then
  usage
fi

OUT="$OUT_DIR/${PREFIX}_blast_result.tsv"
TOPHIT="$OUT_DIR/${PREFIX}_top_6species.tsv"

# 运行 BLAST
blastp -query "$QUERY" -db "$DB" -out "$OUT" -evalue 1e-5 -outfmt 6 -num_threads "$THREADS"

# 提取 top hit
awk '$2 ~ /_(MOUSE|ARATH|MIMIV|BACSU|YEAST|DICDI)/ && !seen[$1]++ {print $1 "\t" $2 "\t" $11}' "$OUT" > "$TOPHIT"

echo "$PREFIX:$(wc -l < ${PREFIX}_top_6species.tsv) hits"
echo "完成！Top hits 已保存至 $TOPHIT"