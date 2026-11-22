#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf 'Usage: %s <comparison_id>\n' "$0" >&2
    exit 2
fi

comparison_id="$1"
project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rscript="Rscript"
config="${project_root}/config/analysis_config.yaml"
log_dir="${project_root}/logs/${comparison_id}"
mkdir -p "${log_dir}"

run_step() {
    local step="$1"
    shift
    local log_file="${log_dir}/${step}.log"
    printf '[%s] START %s %s\n' "$(date -u +%FT%TZ)" "${comparison_id}" "${step}" | tee "${log_file}"
    "${rscript}" "$@" 2>&1 | tee -a "${log_file}"
    printf '[%s] COMPLETE %s %s\n' "$(date -u +%FT%TZ)" "${comparison_id}" "${step}" | tee -a "${log_file}"
}

run_step 03_find_markers "${project_root}/scripts/03_find_og_sum_markers.R" "${config}" "${comparison_id}"
run_step 04_validate_markers "${project_root}/scripts/04_validate_marker_outputs.R" "${config}" "${comparison_id}"
run_step 05_conservation "${project_root}/scripts/05_analyze_conserved_ogs.R" "${config}" "${comparison_id}"
run_step 06_heatmap "${project_root}/scripts/06_plot_conserved_log2fc_heatmap.R" "${config}" "${comparison_id}"
