#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rscript="Rscript"
config="${project_root}/config/analysis_config.yaml"
log_dir="${project_root}/logs"
mkdir -p "${log_dir}"

run_step() {
    local step="$1"
    shift
    local log_file="${log_dir}/${step}.log"
    printf '[%s] START %s\n' "$(date -u +%FT%TZ)" "${step}" | tee "${log_file}"
    "${rscript}" "$@" 2>&1 | tee -a "${log_file}"
    printf '[%s] COMPLETE %s\n' "$(date -u +%FT%TZ)" "${step}" | tee -a "${log_file}"
}

run_step 00_validate_inputs "${project_root}/scripts/00_validate_inputs_and_reused_mapping.R" "${config}"
run_step 01_mapping_qc "${project_root}/scripts/01_summarize_multi_og_mapping.R" "${config}"
run_step 02_build_og_sum_objects "${project_root}/scripts/02_build_og_sum_objects.R" "${config}"
