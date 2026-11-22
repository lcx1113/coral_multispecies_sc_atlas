#!/usr/bin/env bash

set -uo pipefail

ROOT="$(pwd)"
RSCRIPT="Rscript"
DRIVER="${ROOT}/workflow/run_cluster_stability.R"
VALIDATOR="${ROOT}/workflow/scripts/11_validate_260811_formal_config.R"
RUN_LOG_DIR="${ROOT}/logs/formal_subsample_20260811"
STATUS_FILE="${ROOT}/summary/formal_subsample_batch_status.csv"
SAMPLES=(Adig Amur Fusp Gfas Hcoe Pdam Plut Xesp Pcru)

"${RSCRIPT}" --vanilla "${VALIDATOR}" "${ROOT}"

mkdir -p "${RUN_LOG_DIR}"
printf 'species,started_at,completed_at,exit_status,status,log\n' > "${STATUS_FILE}"

overall_status=0
for sample in "${SAMPLES[@]}"; do
  config="${ROOT}/${sample}/config/config.yaml"
  log="${RUN_LOG_DIR}/${sample}.log"
  started_at=$(date '+%Y-%m-%d %H:%M:%S %z')

  echo "[$(date '+%F %T')] Starting ${sample} formal subsample"
  "${RSCRIPT}" --vanilla "${DRIVER}" "${config}" \
    --stage=formal \
    --marker-mode=none \
    --pair-consensus=false \
    --sensitivity=false 2>&1 | tee "${log}"
  status=${PIPESTATUS[0]}

  completed_at=$(date '+%Y-%m-%d %H:%M:%S %z')
  if [[ ${status} -eq 0 ]]; then
    state="completed"
  else
    state="failed"
    overall_status=1
  fi
  printf '"%s","%s","%s",%d,"%s","%s"\n' \
    "${sample}" "${started_at}" "${completed_at}" "${status}" "${state}" "${log}" \
    >> "${STATUS_FILE}"
  echo "[$(date '+%F %T')] Finished ${sample}: ${state}"
done

exit "${overall_status}"
