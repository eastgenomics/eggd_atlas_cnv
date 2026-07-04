#!/usr/bin/env bash
# run_e2e.sh <sample_id> — run the workflow on one raw (Ensembl-named) BAM end-to-end.
# Stage 0 adds the chr prefix. Requires: dx login, a project selected, and
# WORKFLOW_ID / INPUT_BAM / PON_FILE / DESTINATION set (edit below or export).
set -euo pipefail

SAMPLE="${1:?usage: run_e2e.sh <sample_id>}"

WORKFLOW_ID="${WORKFLOW_ID:?export WORKFLOW_ID=workflow-XXXX}"
INPUT_BAM="${INPUT_BAM:?export INPUT_BAM=project-XXXX:file-BAM}"
PON_FILE="${PON_FILE:?export PON_FILE=project-XXXX:file-PON}"
DESTINATION="${DESTINATION:-/atlas_cnv/${SAMPLE}/}"

echo "Running ${WORKFLOW_ID} on ${SAMPLE}..."
ANALYSIS=$(dx run "${WORKFLOW_ID}" \
  -iinput_bam="${INPUT_BAM}" \
  -isample_id="${SAMPLE}" \
  -iploidy_cap_purity_threshold=0.35 \
  -icnvkit_cn_reference="${PON_FILE}" \
  --destination "${DESTINATION}" --priority high --brief -y --watch)

echo "Analysis: ${ANALYSIS}"
echo "Checking igv_html output..."
dx describe "${ANALYSIS}" --json | jq '.output.igv_html // .output'
