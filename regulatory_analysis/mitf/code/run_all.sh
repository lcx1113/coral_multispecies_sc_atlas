#!/usr/bin/env bash
# End-to-end reproduction of Figures 1 and 2.
#
# Usage:
#   bash code/run_all.sh
#
# Environment overrides: MITF_ROOT (output root), MITF_DATA, MITF_MARKERS,
# MITF_GENOME, MITF_GTF, MITF_MOTIFS, PYTHON, JOBS.
set -euo pipefail

BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ROOT=${MITF_ROOT:-$BASE/results}
PYTHON=${PYTHON:-python3}
mkdir -p "$ROOT/tables" "$ROOT/figures" "$ROOT/logs"

echo "== 00  scrambled control matrix"
"$PYTHON" "$BASE/code/00_make_scrambled_motif.py"

echo "== 01  promoter windows"
"$PYTHON" "$BASE/code/01_build_promoters.py"

echo "== 02/03  PWM scans"
bash "$BASE/code/03_run_scans.sh"

echo "== 04  enrichment (Figures 1 and 2 source tables)"
"$PYTHON" "$BASE/code/04_enrichment.py"

echo "== 05  figures"
"$PYTHON" "$BASE/code/05_figures.py"

echo
echo "all done -> $ROOT"
