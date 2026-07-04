#!/bin/bash
# eggd_cgp-cobalt v1.0.0 — COBALT 3.0-beta.5 targeted tumour-only depth ratios
# Converted from cgp-cobalt applet: metadata + timeoutPolicy + execDepends.
# Tool flags and output names are FROZEN (downstream links depend on them).
set -eo pipefail

main() {
    echo "====================================================="
    echo " eggd_cgp-cobalt: COBALT targeted depth ratios"
    echo " Sample  : ${sample_id}"
    echo "====================================================="

    # ── 1. Download inputs ──────────────────────────────────────────────────
    # (system deps come from execDepends — no run-time apt-get)
    java     -version  2>&1 | head -1
    samtools --version 2>&1 | head -1

    echo "[1/4] Downloading inputs..."
    dx download "${tumour_bam}"      -o tumour.bam
    dx download "${tumour_bai}"      -o tumour.bam.bai
    dx download "${cobalt_jar}"      -o cobalt.jar
    dx download "${norm_file}"       -o target_regions_normalisation.tsv
    dx download "${diploid_regions}" -o DiploidRegions.38.bed.gz
    dx download "${gc_profile}"      -o GC_profile.1000bp.38.cnp
    dx download "${ref_fasta}"       -o ref.fa.gz
    dx download "${ref_fai}"         -o ref.fa.fai
    # COBALT needs the .fai next to the gz
    ln -sf ref.fa.fai ref.fa.gz.fai

    # ── 2. Run COBALT ───────────────────────────────────────────────────────
    echo "[2/4] Running COBALT..."
    mkdir -p "${sample_id}"

    java -Xmx12G -jar cobalt.jar \
        -tumor                   "${sample_id}" \
        -tumor_bam               tumour.bam \
        -target_region_norm_file target_regions_normalisation.tsv \
        -tumor_only_diploid_bed  DiploidRegions.38.bed.gz \
        -gc_profile              GC_profile.1000bp.38.cnp \
        -ref_genome              ref.fa.gz \
        -ref_genome_version      38 \
        -bam_validation          SILENT \
        -threads                 "$(nproc)" \
        -output_dir              "${sample_id}/"

    # ── 3. Verify outputs ──────────────────────────────────────────────────
    echo "[3/4] Verifying outputs..."
    RATIO_FILE="${sample_id}/${sample_id}.cobalt.ratio.tsv.gz"
    PCF_FILE="${sample_id}/${sample_id}.cobalt.ratio.pcf"
    [[ -s "${RATIO_FILE}" ]] || { echo "ERROR: COBALT ratio TSV missing"; exit 1; }
    [[ -s "${PCF_FILE}" ]]   || { echo "ERROR: COBALT ratio PCF missing"; exit 1; }

    WINDOWS=$(zcat "${RATIO_FILE}" | tail -n +2 | wc -l)
    echo "Ratio windows: ${WINDOWS}"
    ls -lh "${sample_id}/"

    # ── 4. Tar and upload ──────────────────────────────────────────────────
    echo "[4/4] Uploading..."
    tar --no-same-owner -czf "${sample_id}.cobalt.tar.gz" "${sample_id}/"

    cobalt_tar=$(dx upload "${sample_id}.cobalt.tar.gz" --brief)
    dx-jobutil-add-output cobalt_tar "${cobalt_tar}" --class=file

    echo "====================================================="
    echo " eggd_cgp-cobalt DONE: ${sample_id}  windows=${WINDOWS}"
    echo "====================================================="
}
