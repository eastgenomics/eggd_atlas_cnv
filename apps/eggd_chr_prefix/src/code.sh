#!/bin/bash
# eggd_chr_prefix v1.0.0 (folded-in stage 0) — add the `chr` prefix to a BAM header.
# Adapted from the public eggd_chr_prefix for SINGLE-FILE per-sample use with PASSTHROUGH:
#   - single-file I/O (input_bam -> output_bam + output_bai), no arrays
#   - mode fixed to add_chr
#   - ALWAYS emits a BAM + index (reheadered when needed, else the input passed through);
#     it must never emit nothing, which would starve downstream workflow stages.
# Header-only reheader (samtools reheader) — read data is never altered.
set -eo pipefail

main() {
    echo "====================================================="
    echo " eggd_chr_prefix (stage 0): add_chr, single-file + passthrough"
    echo "====================================================="
    samtools --version 2>&1 | head -1

    echo "[1/4] Downloading input BAM..."
    dx download "${input_bam}" -o in.bam
    NAME=$(dx describe --name "${input_bam}")
    OUT="${NAME}"                              # preserve the original filename

    # add-chr header transform (SN: sequence names):
    #   1..22 -> chr1..chr22, X/Y -> chrX/chrY, MT -> chrM (matches chr-prefixed GRCh38)
    ADD_CHR='s/\tSN:\([0-9][0-9]*\)\t/\tSN:chr\1\t/g; s/\tSN:X\t/\tSN:chrX\t/g; s/\tSN:Y\t/\tSN:chrY\t/g; s/\tSN:MT\t/\tSN:chrM\t/g'

    echo "[2/4] Extracting + transforming header..."
    samtools view -H in.bam > orig.sam
    sed "${ADD_CHR}" orig.sam > new.sam

    echo "[3/4] Reheader or passthrough..."
    if cmp -s orig.sam new.sam; then
        echo "  Header already chr-prefixed — passthrough."
        cp in.bam "${OUT}"
    else
        echo "  Reheadering to add chr prefix."
        samtools reheader new.sam in.bam > "${OUT}"
    fi
    samtools index "${OUT}" "${OUT}.bai"

    echo "[4/4] Uploading..."
    dx-jobutil-add-output output_bam "$(dx upload "${OUT}"     --brief)" --class=file
    dx-jobutil-add-output output_bai "$(dx upload "${OUT}.bai" --brief)" --class=file

    echo "====================================================="
    echo " eggd_chr_prefix DONE: ${OUT}"
    echo "====================================================="
}
