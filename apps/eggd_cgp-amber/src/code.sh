#!/bin/bash
# eggd_cgp-amber v1.0.0 — AMBER 4.3-beta.4 tumour-only BAF
# Converted from cgp-amber applet: metadata + timeoutPolicy + execDepends.
# Tool flags and output names are FROZEN (downstream links depend on them).
set -eo pipefail

main() {
    echo "====================================================="
    echo " eggd_cgp-amber: AMBER tumour-only BAF"
    echo " Sample  : ${sample_id}"
    echo "====================================================="

    # ── 1. Download inputs ──────────────────────────────────────────────────
    # (system deps come from execDepends — no run-time apt-get)
    java     -version  2>&1 | head -1
    samtools --version 2>&1 | head -1

    echo "[1/4] Downloading inputs..."
    dx download "${tumour_bam}"     -o tumour.bam
    dx download "${tumour_bai}"     -o tumour.bam.bai
    dx download "${amber_jar}"      -o amber.jar
    dx download "${germline_sites}" -o germline_sites.tsv.gz

    # ── 2. Run AMBER ────────────────────────────────────────────────────────
    echo "[2/4] Running AMBER..."
    mkdir -p "${sample_id}"

    java -Xmx6G -jar amber.jar \
        -tumor             "${sample_id}" \
        -tumor_bam         tumour.bam \
        -loci              germline_sites.tsv.gz \
        -tumor_min_depth   10 \
        -bam_validation    LENIENT \
        -ref_genome_version 38 \
        -threads           "$(nproc)" \
        -output_dir        "${sample_id}/"

    # ── 3. Verify outputs ──────────────────────────────────────────────────
    echo "[3/4] Verifying outputs..."
    BAF_FILE="${sample_id}/${sample_id}.amber.baf.tsv.gz"
    PCF_FILE="${sample_id}/${sample_id}.amber.baf.pcf"
    [[ -s "${BAF_FILE}" ]] || { echo "ERROR: AMBER BAF TSV missing"; exit 1; }
    [[ -s "${PCF_FILE}" ]] || { echo "ERROR: AMBER BAF PCF missing (Java PCF failed?)"; exit 1; }

    SITES=$(zcat "${BAF_FILE}" | tail -n +2 | wc -l)
    echo "BAF sites: ${SITES}"
    ls -lh "${sample_id}/"

    # ── 4. Tar and upload ──────────────────────────────────────────────────
    echo "[4/4] Uploading..."
    tar --no-same-owner -czf "${sample_id}.amber.tar.gz" "${sample_id}/"

    amber_tar=$(dx upload "${sample_id}.amber.tar.gz" --brief)
    dx-jobutil-add-output amber_tar "${amber_tar}" --class=file

    echo "====================================================="
    echo " eggd_cgp-amber DONE: ${sample_id}  sites=${SITES}"
    echo "====================================================="
}
