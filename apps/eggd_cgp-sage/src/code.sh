#!/bin/bash
# eggd_cgp-sage v1.0.0 — SAGE 5.0-beta.11 somatic SNV/indel calling, tumour-only, CGP+backbone panel
# Converted from cgp-sage applet: metadata + timeoutPolicy + execDepends.
# Tool flags and output names are FROZEN (downstream links depend on them).
set -euxo pipefail

main() {
    case "${sample_id}" in
        *[!A-Za-z0-9._-]* | "" | .* | -* )
            echo "ERROR: unsafe sample_id '${sample_id}' (allowed: A-Za-z0-9._-, no leading '-'/'.')" >&2; exit 1 ;;
    esac
    echo "====================================================="
    echo " eggd_cgp-sage: SAGE somatic calling"
    echo " Sample  : ${sample_id}"
    echo " Instance: $(hostname)"
    echo "====================================================="

    # ── 1. Verify deps (from execDepends — no run-time apt-get) ─────────────
    java         -version  2>&1 | sed -n '1p'
    samtools     --version 2>&1 | sed -n '1p'

    # ── 2. Download inputs ─────────────────────────────────────────────────
    echo "[1/4] Downloading inputs..."
    declare -A DL_PIDS
    dx download "${tumour_bam}"   -o tumour.bam &         DL_PIDS[$!]="tumour_bam"
    dx download "${tumour_bai}"   -o tumour.bam.bai &     DL_PIDS[$!]="tumour_bai"
    dx download "${sage_jar}"     -o sage.jar &            DL_PIDS[$!]="sage_jar"
    dx download "${ref_fasta}"    -o ref.fa.gz &           DL_PIDS[$!]="ref_fasta"
    dx download "${hotspots_vcf}" -o hotspots.vcf.gz &     DL_PIDS[$!]="hotspots_vcf"
    dx download "${hotspots_tbi}" -o hotspots.vcf.gz.tbi & DL_PIDS[$!]="hotspots_tbi"
    dx download "${panel_bed}"    -o panel.bed &           DL_PIDS[$!]="panel_bed"
    dx download "${hc_bed}"       -o hc.bed.gz &           DL_PIDS[$!]="hc_bed"
    dx download "${pon_file}"     -o pon.tsv.gz &          DL_PIDS[$!]="pon_file"
    dx download "${ensembl_data}" -o ensembl_data.tar.gz & DL_PIDS[$!]="ensembl_data"
    for pid in "${!DL_PIDS[@]}"; do
        wait "$pid" || { echo "ERROR: download failed for input '${DL_PIDS[$pid]}'" >&2; exit 1; }
    done

    # ── 3. Prepare reference + Ensembl ─────────────────────────────────────
    echo "[2/4] Preparing reference and Ensembl data..."
    # SAGE (htsjdk) cannot handle bgzf-compressed FASTA — decompress to plain .fa
    echo "Decompressing reference FASTA (~3 GB, ~1 min)..."
    bgzip -d -c ref.fa.gz > ref.fa
    samtools faidx ref.fa
    samtools dict ref.fa > ref.dict   # htsjdk requires .dict for sequenceDictionary

    mkdir -p ensembl_data
    tar --no-same-owner -xzf ensembl_data.tar.gz -C ensembl_data
    ENSEMBL_DIR=$(find ensembl_data -maxdepth 2 -name "ensembl_gene_data.csv" -printf "%h\n" -quit)
    [[ -n "${ENSEMBL_DIR}" ]] || { echo "ERROR: ensembl_gene_data.csv not found" >&2; exit 1; }
    echo "Ensembl dir: ${ENSEMBL_DIR}"

    # ── 3b. Validate chr-prefix on all GRCh38 inputs ──────────────────────────
    echo "Verifying chr-prefixed contigs on all inputs..."
    # Use capture-then-test to avoid SIGPIPE false-failures under set -o pipefail
    bam_sq=$(samtools view -H tumour.bam | awk '/^@SQ/ && /SN:chr/{print; exit}' || true)
    [[ -n "${bam_sq}" ]] \
        || { echo "ERROR: tumour BAM does not have chr-prefixed contigs" >&2; exit 1; }
    fai_chr=$(head -1 ref.fa.fai | cut -f1 || true)
    [[ "${fai_chr}" == chr* ]] \
        || { echo "ERROR: ref FASTA does not have chr-prefixed contigs" >&2; exit 1; }
    vcf_chr=$(zcat hotspots.vcf.gz | grep -m1 '^[^#]' | cut -f1 || true)
    [[ "${vcf_chr}" == chr* ]] \
        || { echo "ERROR: hotspots VCF does not have chr-prefixed contigs" >&2; exit 1; }
    bed_chr=$(grep -m1 '^[^#]' panel.bed | cut -f1 || true)
    [[ "${bed_chr}" == chr* ]] \
        || { echo "ERROR: panel BED does not have chr-prefixed contigs" >&2; exit 1; }
    hc_chr=$(zcat hc.bed.gz | grep -m1 '^[^#]' | cut -f1 || true)
    [[ "${hc_chr}" == chr* ]] \
        || { echo "ERROR: hc_bed does not have chr-prefixed contigs" >&2; exit 1; }
    pon_chr=$(zcat pon.tsv.gz | grep -m1 '^[^#]' | cut -f1 || true)
    [[ "${pon_chr}" == chr* ]] \
        || { echo "ERROR: PON file does not have chr-prefixed contigs" >&2; exit 1; }

    # ── 4. Run SAGE ────────────────────────────────────────────────────────
    echo "[3/4] Running SAGE..."
    THREADS=$(nproc)
    # Derive JVM heap from available RAM, leaving ~2 GiB headroom
    HEAP_MB=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 - 2048 ))
    (( HEAP_MB > 512 )) \
        || { echo "ERROR: insufficient RAM for SAGE — HEAP_MB=${HEAP_MB} (instance needs >2.5 GiB total RAM)" >&2; exit 1; }
    mkdir -p sage_out/

    java -Xmx${HEAP_MB}m -jar sage.jar \
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

    TOTAL_COUNT=$(zcat "sage_out/${sample_id}.sage.somatic.vcf.gz" \
        | awk '!/^#/{n++} END{print n+0}')
    PASS_COUNT=$(zcat "sage_out/${sample_id}.sage.somatic.vcf.gz" \
        | awk '!/^#/ && $7=="PASS"{n++} END{print n+0}')
    echo "Variants: PASS=${PASS_COUNT}/total=${TOTAL_COUNT}"

    somatic_vcf=$(dx upload     "sage_out/${sample_id}.sage.somatic.vcf.gz"     --brief)
    somatic_vcf_tbi=$(dx upload "sage_out/${sample_id}.sage.somatic.vcf.gz.tbi" --brief)
    dx-jobutil-add-output somatic_vcf     "${somatic_vcf}"     --class=file
    dx-jobutil-add-output somatic_vcf_tbi "${somatic_vcf_tbi}" --class=file

    echo "====================================================="
    echo " eggd_cgp-sage DONE: ${sample_id}  PASS=${PASS_COUNT}/total=${TOTAL_COUNT}"
    echo "====================================================="
}
