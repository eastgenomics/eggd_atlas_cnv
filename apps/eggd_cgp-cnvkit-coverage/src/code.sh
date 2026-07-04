#!/bin/bash
# cgp-cnvkit-coverage/src/code.sh
# Computes CNVkit on-target coverage for one sample (amplicon mode).
# No antitarget bins — backbone CNV signal is handled by PURPLE.
# Passes the gene panel BED directly to cnvkit.py coverage (no autobin).
# Output: <sample_id>.targetcoverage.cnn
set -eo pipefail

main() {
    echo "=== CGP CNVkit coverage: ${sample_id} ==="
    echo "Instance: $(hostname); CPUs: $(nproc); RAM: $(free -h | awk '/Mem:/{print $2}')"

    # ── Load CNVkit Docker image ──────────────────────────────────────────────
    # Image stored in DNAnexus; no external internet required.
    # Update CNVKIT_IMAGE_ID after running scripts/dnanexus/docker/cgp-cnvkit/build_and_upload.sh
    CNVKIT_IMAGE_ID="project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q"   # cgp-cnvkit:1.0.0 — set after upload
    CNVKIT_IMAGE_TAG="cgp-cnvkit:1.0.0"
    echo "[setup] Loading CNVkit image..."
    dx download "${CNVKIT_IMAGE_ID}" -o cnvkit-image.tar.gz
    docker load < cnvkit-image.tar.gz
    run_cnvkit() { docker run --rm -v "$(pwd)":/work -w /work "${CNVKIT_IMAGE_TAG}" cnvkit.py "$@"; }
    run_cnvkit version

    # ── Stage inputs ──────────────────────────────────────────────────────────
    echo "[inputs] Downloading files..."
    dx download "${tumour_bam}" -o tumour.bam
    dx download "${tumour_bai}" -o tumour.bam.bai
    dx download "${baits}"      -o panel.bed

    # Deduplicate BED intervals — some panel BEDs have identical coordinates
    # from multiple transcript annotations; CNVkit 0.9.13+ raises on duplicates
    awk '!seen[$1 FS $2 FS $3]++' panel.bed > panel.dedup.bed && mv panel.dedup.bed panel.bed

    NBINS=$(grep -vc "^#" panel.bed || true)
    echo "[inputs] Panel intervals: ${NBINS}"

    # ── Compute coverage ──────────────────────────────────────────────────────
    # Amplicon mode: pass panel BED directly to coverage (no autobin needed).
    # Each target interval becomes one bin. CNVkit weights bins by depth so
    # small exon intervals are down-weighted automatically.
    echo "[coverage] Computing on-target coverage..."
    run_cnvkit coverage tumour.bam panel.bed \
        --processes "$(nproc)" \
        --output "${sample_id}.targetcoverage.cnn"

    # Sanity check: verify the coverage file has sufficient intervals and depth.
    # NLINES > 1000 confirms CNVkit produced output for most panel targets
    # (the Twist CGP panel has ~20,000 intervals; a near-empty file indicates
    # a BAM/BED mismatch such as mismatched chromosome naming).
    # MEDIAN_DEPTH is logged for QC — very low values (<10x) suggest the sample
    # failed sequencing or the wrong BAM was provided.
    NLINES=$(wc -l < "${sample_id}.targetcoverage.cnn")
    echo "[coverage] Output lines (incl. header): ${NLINES}"
    [ "${NLINES}" -gt 1000 ] || { echo "ERROR: coverage file too small (${NLINES} lines)"; exit 1; }

    MEDIAN_DEPTH=$(python3 -c "
import csv
depths = []
with open('${sample_id}.targetcoverage.cnn') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        try: depths.append(float(row.get('depth', 0)))
        except: pass
depths.sort()
print(f'{depths[len(depths)//2]:.1f}' if depths else '0')
")
    echo "[coverage] Median depth: ${MEDIAN_DEPTH}"

    # ── Upload output ─────────────────────────────────────────────────────────
    coverage_cnn=$(dx upload "${sample_id}.targetcoverage.cnn" --brief)
    dx-jobutil-add-output coverage_cnn "${coverage_cnn}" --class=file

    echo "=== Done: ${sample_id} ==="
}
