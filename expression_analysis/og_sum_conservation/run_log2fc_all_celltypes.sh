#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rscript="Rscript"
config="${project_root}/config/analysis_config.yaml"
manifest="${project_root}/config/species_manifest.tsv"
output_dir="${project_root}/results/LOG2FC"
max_parallel="${1:-2}"

if ! [[ "${max_parallel}" =~ ^[1-9][0-9]*$ ]]; then
    printf 'max_parallel must be a positive integer\n' >&2
    exit 2
fi

mkdir -p "${output_dir}/logs"
timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
main_log="${output_dir}/logs/run_all_${timestamp}.log"
printf 'started_utc\tcommand\tconfig\tmax_parallel\n%s\t%s\t%s\t%s\n' \
    "$(date -u +%FT%TZ)" "$0 $*" "${config}" "${max_parallel}" \
    > "${output_dir}/run_invocation_${timestamp}.tsv"

mapfile -t species < <(awk -F '\t' 'NR > 1 && $3 == "TRUE" {print $1}' "${manifest}")

run_species() {
    local sp="$1"
    local log_file="${output_dir}/logs/${sp}_${timestamp}.log"
    printf '[%s] START %s\n' "$(date -u +%FT%TZ)" "${sp}" | tee -a "${main_log}"
    if "${rscript}" "${project_root}/scripts/M03_calculate_all_celltype_og_log2fc.R" \
        "${config}" "${sp}" > "${log_file}" 2>&1; then
        printf '[%s] COMPLETE %s\n' "$(date -u +%FT%TZ)" "${sp}" | tee -a "${main_log}"
    else
        printf '[%s] FAILED %s; inspect %s\n' "$(date -u +%FT%TZ)" "${sp}" "${log_file}" | tee -a "${main_log}"
        return 1
    fi
}

batch_pids=()
batch_species=()
failed=0
for sp in "${species[@]}"; do
    run_species "${sp}" &
    batch_pids+=("$!")
    batch_species+=("${sp}")
    if [[ "${#batch_pids[@]}" -ge "${max_parallel}" ]]; then
        for i in "${!batch_pids[@]}"; do
            if ! wait "${batch_pids[$i]}"; then
                failed=1
            fi
        done
        batch_pids=()
        batch_species=()
    fi
done
for i in "${!batch_pids[@]}"; do
    if ! wait "${batch_pids[$i]}"; then
        failed=1
    fi
done
if [[ "${failed}" -ne 0 ]]; then
    printf 'At least one species failed; validation was not run\n' | tee -a "${main_log}" >&2
    exit 1
fi

printf '[%s] START validation\n' "$(date -u +%FT%TZ)" | tee -a "${main_log}"
"${rscript}" "${project_root}/scripts/M04_validate_all_celltype_og_log2fc.R" \
    "${config}" 2>&1 | tee "${output_dir}/logs/validation_${timestamp}.log" | tee -a "${main_log}"
printf '[%s] COMPLETE validation\n' "$(date -u +%FT%TZ)" | tee -a "${main_log}"
