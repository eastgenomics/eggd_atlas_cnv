#!/bin/bash
# eggd_cgp-sage v1.0.0 — SAGE 5.0-beta.11 somatic SNV/indel calling, tumour-only, CGP+backbone panel
# Converted from cgp-sage applet: metadata + timeoutPolicy + execDepends.
# Tool flags and output names are FROZEN (downstream links depend on them).
set -eo pipefail

main() {
    echo "====================================================="
    echo " eggd_cgp-sage: SAGE somatic calling"
    echo " Sample  : ${sample_id}"
    echo " Instance: $(hostname)"
    echo "====================================================="

    # ── 1. Verify deps (from execDepends — no run-time apt-get) ─────────────
    java -version

    # ── 2. Download inputs ─────────────────────────────────────────────────
    echo "[1/4] Downloading inputs..."
    dx download "${tumour_bam}"   -o tumour.bam
    dx download "${tumour_bai}"   -o tumour.bam.bai
    dx download "${sage_jar}"     -o sage.jar
    dx download "${ref_fasta}"    -o ref.fa.gz
    dx download "${ref_fai}"      -o ref.fa.fai
    dx download "${hotspots_vcf}" -o hotspots.vcf.gz
    dx download "${hotspots_tbi}" -o hotspots.vcf.gz.tbi
    dx download "${panel_bed}"    -o panel.bed
    dx download "${hc_bed}"       -o hc.bed.gz
    dx download "${pon_file}"     -o pon.tsv.gz
    dx download "${ensembl_data}" -o ensembl_data.tar.gz

    # ── 3. Prepare reference + Ensembl ─────────────────────────────────────
    echo "[2/4] Preparing reference and Ensembl data..."
    # SAGE (htsjdk) cannot handle bgzf-compressed FASTA — decompress to plain .fa
    echo "Decompressing reference FASTA (~3 GB, ~1 min)..."
    bgzip -d -c ref.fa.gz > ref.fa
    samtools faidx ref.fa
    samtools dict ref.fa > ref.dict   # htsjdk requires .dict for sequenceDictionary

    tar --no-same-owner -xzf ensembl_data.tar.gz
    ENSEMBL_DIR=$(find . -maxdepth 2 -name "ensembl_gene_data.csv" -printf "%h\n" | head -1)
    [[ -n "${ENSEMBL_DIR}" ]] || { echo "ERROR: ensembl_gene_data.csv not found" >&2; exit 1; }
    echo "Ensembl dir: ${ENSEMBL_DIR}"

    # ── 4. Run SAGE ────────────────────────────────────────────────────────
    echo "[3/4] Running SAGE..."
    THREADS=$(nproc)
    MEM_GB=$(( $(free -g | awk '/^Mem:/{print $2}') - 4 ))
    mkdir -p sage_out/

    java -Xmx${MEM_GB}G -jar sage.jar \
        -tumor               "${sample_id}" \
        -tumor_bam           tumour.bam \
        -ref_genome          ref.fa \
        -ref_genome_version  38 \
        -hotspots            hotspots.vcf.gz \
        -panel_bed           panel.bed \
        -high_confidence_bed hc.bed.gz \
        -pon_file            pon.tsv.gz \
        -ensembl_data_dir    "${ENSEMBL_DIR}" \
        -panel_only \
        -high_depth_mode \
        -skip_msi_jitter \
        -skip_bqr \
        -ref_sample_count    0 \
        -bam_validation      SILENT \
        -output_vcf          "sage_out/${sample_id}.sage.somatic.vcf.gz" \
        -threads             "${THREADS}"

    # ── 5. Upload outputs ──────────────────────────────────────────────────
    echo "[4/4] Uploading outputs..."
    [[ -f "sage_out/${sample_id}.sage.somatic.vcf.gz" ]] \
        || { echo "ERROR: SAGE did not produce output VCF" >&2; ls sage_out/ >&2; exit 1; }

    # SAGE 5.x should auto-index; fallback if not
    [[ -f "sage_out/${sample_id}.sage.somatic.vcf.gz.tbi" ]] \
        || tabix -p vcf "sage_out/${sample_id}.sage.somatic.vcf.gz"

    PASS_COUNT=$(zcat "sage_out/${sample_id}.sage.somatic.vcf.gz" \
        | awk '!/^#/ && $7=="PASS"' | wc -l)
    echo "PASS variants: ${PASS_COUNT}"

    somatic_vcf=$(dx upload     "sage_out/${sample_id}.sage.somatic.vcf.gz"     --brief)
    somatic_vcf_tbi=$(dx upload "sage_out/${sample_id}.sage.somatic.vcf.gz.tbi" --brief)
    dx-jobutil-add-output somatic_vcf     "${somatic_vcf}"     --class=file
    dx-jobutil-add-output somatic_vcf_tbi "${somatic_vcf_tbi}" --class=file

    echo "====================================================="
    echo " eggd_cgp-sage DONE: ${sample_id}  PASS=${PASS_COUNT}"
    echo "====================================================="
}
