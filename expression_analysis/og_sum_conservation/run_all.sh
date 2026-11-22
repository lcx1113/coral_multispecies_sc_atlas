#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${project_root}/run_shared.sh"

while IFS= read -r comparison_id; do
    "${project_root}/run_comparison.sh" "${comparison_id}"
done < <(awk -F '\t' 'NR > 1 && $3 == "TRUE" {print $1}' "${project_root}/config/celltype_comparisons.tsv" | sort -u)
