#!/usr/bin/env bash
# run_e2e.sh <sample_id> — run the workflow on a chr-prefixed BAM and matching BAI.
# eggd_chr_prefix runs separately. Requires: dx login, a project selected, and
# WORKFLOW_ID / INPUT_BAM / INPUT_BAI / PON_FILE / DESTINATION set (edit below or export).
set -euo pipefail

SAMPLE="${1:?usage: run_e2e.sh <sample_id>}"

WORKFLOW_ID="${WORKFLOW_ID:?export WORKFLOW_ID=workflow-XXXX}"
INPUT_BAM="${INPUT_BAM:?export INPUT_BAM=project-XXXX:file-CHR-PREFIXED-BAM}"
INPUT_BAI="${INPUT_BAI:?export INPUT_BAI=project-XXXX:file-CHR-PREFIXED-BAI}"
PON_FILE="${PON_FILE:?export PON_FILE=project-XXXX:file-PON}"
DESTINATION="${DESTINATION:-/atlas_cnv/${SAMPLE}/}"

echo "Running ${WORKFLOW_ID} on ${SAMPLE}..."
ANALYSIS=$(dx run "${WORKFLOW_ID}" \
  -iinput_bam="${INPUT_BAM}" \
  -iinput_bai="${INPUT_BAI}" \
  -isample_id="${SAMPLE}" \
  -iploidy_cap_purity_threshold=0.35 \
  -icnvkit_cn_reference="${PON_FILE}" \
  --destination "${DESTINATION}" --priority high --brief -y --watch)

echo "Analysis: ${ANALYSIS}"
echo "Checking igv_html output..."
dx describe "${ANALYSIS}" --json | jq '.output.igv_html // .output'
