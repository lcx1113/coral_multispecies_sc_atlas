#!/usr/bin/env bash
# Scan the proximal promoter window (TSS -2,000/+500 bp) with the eleven matrices
# shown in Figure 2.
#
# Environment overrides: MITF_ROOT, MITF_MOTIFS, PYTHON, JOBS.
set -euo pipefail

BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ROOT=${MITF_ROOT:-$BASE/results}
MOTIF_DIR=${MITF_MOTIFS:-$BASE/motifs}
PYTHON=${PYTHON:-python3}
JOBS=${JOBS:-8}
SCANNER="$BASE/code/02_scan_motifs.py"

MOTIFS="MA0620.4 MA1899.1 MA2257.1 MA0831.3 MA0692.2 SCD_MA0620.4 \
MA0036.3 MA0143.5 MA0901.3 MA1870.2 MA1109.2"
WINDOW="prox2kb"

if [[ ! -s "$ROOT/promoters.fa" || ! -s "$ROOT/gene_tss.tsv" ]]; then
  echo "ERROR: $ROOT/promoters.fa and $ROOT/gene_tss.tsv are required." >&2
  echo "       Run code/01_build_promoters.py first (or set MITF_ROOT)." >&2
  exit 1
fi

mkdir -p "$ROOT/scans" "$ROOT/logs"

scan() {
  local m=$1
  "$PYTHON" "$SCANNER" \
      --motif "$MOTIF_DIR/$m.jaspar" \
      --fasta "$ROOT/promoters.fa" \
      --tss "$ROOT/gene_tss.tsv" \
      --out "$ROOT/scans/$m.$WINDOW" \
      --windows "$WINDOW" --bg uniform \
      > "$ROOT/logs/$m.$WINDOW.log" 2>&1
}
export -f scan
export ROOT MOTIF_DIR PYTHON SCANNER WINDOW

for m in $MOTIFS; do
  echo "$m"
done | xargs -P "$JOBS" -n 1 bash -c 'scan "$0"'

echo "done"
ls -l "$ROOT/scans" | head -30
