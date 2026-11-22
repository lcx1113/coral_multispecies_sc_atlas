#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rscript="Rscript"
config="${project_root}/config/analysis_config.yaml"
comparison_config="${project_root}/config/celltype_comparisons.tsv"
log_dir="${project_root}/logs/downregulated_og"
mkdir -p "${log_dir}"
timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
main_log="${log_dir}/run_all_${timestamp}.log"

mapfile -t comparisons < <(awk -F '\t' 'NR > 1 && $3 == "TRUE" {print $1}' "${comparison_config}" | sort -u)
for comparison_id in "${comparisons[@]}"; do
    printf '[%s] START %s\n' "$(date -u +%FT%TZ)" "${comparison_id}" | tee -a "${main_log}"
    "${rscript}" "${project_root}/scripts/M05_analyze_downregulated_conserved_ogs.R" \
        "${config}" "${comparison_id}" 2>&1 | tee "${log_dir}/${comparison_id}_analysis_${timestamp}.log" | tee -a "${main_log}"

    summary_file="${project_root}/results/comparisons/${comparison_id}/downregulated/conservation/${comparison_id}.OG_SUM.downregulated.conservation_summary.tsv"
    final_file="${project_root}/results/comparisons/${comparison_id}/downregulated/conservation/${comparison_id}.OG_SUM.downregulated.conservation_summary.with_final_annotation.tsv"
    "${rscript}" "${project_root}/scripts/M02_add_final_annotation_to_conservation_summary.R" \
        "${summary_file}" "${final_file}" --overwrite \
        2>&1 | tee "${log_dir}/${comparison_id}_annotation_${timestamp}.log" | tee -a "${main_log}"

    "${rscript}" "${project_root}/scripts/M06_validate_downregulated_conserved_ogs.R" \
        "${config}" "${comparison_id}" 2>&1 | tee "${log_dir}/${comparison_id}_validation_${timestamp}.log" | tee -a "${main_log}"
    printf '[%s] COMPLETE %s\n' "$(date -u +%FT%TZ)" "${comparison_id}" | tee -a "${main_log}"
done

printf 'started_utc\tfinished_utc\tcommand\tconfig\n%s\t%s\t%s\t%s\n' \
    "${timestamp}" "$(date -u +%FT%TZ)" "$0" "${config}" \
    > "${log_dir}/run_invocation_${timestamp}.tsv"
