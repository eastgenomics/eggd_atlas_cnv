#!/bin/bash
# eggd_cgp-cobalt v1.0.0 — COBALT 3.0-beta.5 targeted tumour-only depth ratios
# Converted from cgp-cobalt applet: metadata + timeoutPolicy + execDepends.
# Tool flags and output names are FROZEN (downstream links depend on them).
set -euo pipefail

main() {
    case "${sample_id}" in
        *[!A-Za-z0-9._-]* | "" | .* | -* )
            echo "ERROR: unsafe sample_id '${sample_id}' (allowed: A-Za-z0-9._-, no leading '-'/'.')" >&2; exit 1 ;;
    esac
    echo "====================================================="
    echo " eggd_cgp-cobalt: COBALT targeted depth ratios"
    echo " Sample  : ${sample_id}"
    echo "====================================================="

    # ── 1. Download inputs ──────────────────────────────────────────────────
    # (system deps come from execDepends — no run-time apt-get)
    java -version  2>&1 | sed -n '1p'

    echo "[1/4] Downloading inputs..."
    pids=()
    dx download "${tumour_bam}"      -o tumour.bam                       & pids+=("$!")
    dx download "${tumour_bai}"      -o tumour.bam.bai                   & pids+=("$!")
    dx download "${cobalt_jar}"      -o cobalt.jar                       & pids+=("$!")
    dx download "${norm_file}"       -o target_regions_normalisation.tsv & pids+=("$!")
    dx download "${diploid_regions}" -o DiploidRegions.38.bed.gz         & pids+=("$!")
    dx download "${gc_profile}"      -o GC_profile.1000bp.38.cnp         & pids+=("$!")
    dx download "${ref_fasta}"       -o ref.fa.gz                        & pids+=("$!")
    dx download "${ref_fai}"         -o ref.fa.fai                       & pids+=("$!")
    for pid in "${pids[@]}"; do wait "${pid}" || { echo "ERROR: dx download failed (pid ${pid})"; exit 1; }; done
    # COBALT needs the .fai next to the gz
    ln -sf ref.fa.fai ref.fa.gz.fai

    # ── 1b. Validate chr-prefix on GRCh38 inputs ───────────────────────────
    echo "Verifying chr-prefixed contigs on GRCh38 inputs..."
    # Use capture-then-test to avoid SIGPIPE false-failures under set -o pipefail
    bam_sq=$(samtools view -H tumour.bam | awk '/^@SQ/ && /SN:chr/{print; exit}' || true)
    [[ -n "${bam_sq}" ]] \
        || { echo "ERROR: tumour BAM does not have chr-prefixed contigs" >&2; exit 1; }
    bed_chr=$(zcat DiploidRegions.38.bed.gz | grep -m1 '^[^#]' | cut -f1 || true)
    [[ "${bed_chr}" == chr* ]] \
        || { echo "ERROR: diploid_regions does not have chr-prefixed contigs" >&2; exit 1; }
    cnp_chr=$(grep -m1 '^[^#]' GC_profile.1000bp.38.cnp | cut -f1 || true)
    [[ "${cnp_chr}" == chr* ]] \
        || { echo "ERROR: gc_profile does not have chr-prefixed contigs" >&2; exit 1; }

    # ── 2. Run COBALT ───────────────────────────────────────────────────────
    echo "[2/4] Running COBALT..."
    mkdir -p "${sample_id}"

    HEAP_MB=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 - 2048 ))
    java -Xmx${HEAP_MB}m -jar cobalt.jar \
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
    [[ "${WINDOWS}" -gt 0 ]] || { echo "ERROR: COBALT ratio TSV has no data rows"; exit 1; }
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
