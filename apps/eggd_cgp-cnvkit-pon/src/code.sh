#!/bin/bash
# cgp-cnvkit-pon/src/code.sh
# Builds a CNVkit PoN reference from per-sample coverage files.
# Optionally annotates intervals with GC content (pass fasta/fasta_fai/fasta_gzi)
# for GC-bias correction in downstream fix steps.
# Outputs reference.cnn + a summary stats TSV.
# shellcheck disable=SC2154,SC2086  # SC2154: vars injected by DNAnexus; SC2086: FASTA_ARG intentionally unquoted
set -eo pipefail

main() {
    # coverage_files, pon_name, fasta, fasta_fai, fasta_gzi are injected by
    # DNAnexus at runtime from dxapp.json inputSpec; ShellCheck cannot see this.
    echo "=== CGP CNVkit PoN build ==="

    # ── Load CNVkit Docker image ──────────────────────────────────────────────
    # Image stored in DNAnexus; no external internet required.
    # Update CNVKIT_IMAGE_ID after running scripts/dnanexus/docker/cgp-cnvkit/build_and_upload.sh
    CNVKIT_IMAGE_ID="project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q"   # cgp-cnvkit:1.0.0
    CNVKIT_IMAGE_TAG="cgp-cnvkit:1.0.0"
    echo "[setup] Loading CNVkit image..."
    dx download "${CNVKIT_IMAGE_ID}" -o cnvkit-image.tar.gz
    docker load < cnvkit-image.tar.gz
    run_cnvkit() { docker run --rm -v "$(pwd)":/work -w /work "${CNVKIT_IMAGE_TAG}" cnvkit.py "$@"; }
    run_cnvkit version

    # ── Download all coverage files ───────────────────────────────────────────
    echo "[inputs] Downloading coverage files..."
    mkdir -p cnn_files

    # coverage_files is an array — DNAnexus passes as bash array of file IDs
    # dx-download-all-inputs downloads to ~/in/coverage_files/N/filename
    dx-download-all-inputs 2>/dev/null || true

    # Collect from default dx-download-all-inputs path, or download manually
    if ls ~/in/coverage_files/*/*.cnn 2>/dev/null | head -1 | grep -q ".cnn"; then
        cp ~/in/coverage_files/*/*.cnn cnn_files/
    else
        # Manual download if dx-download-all-inputs not available
        for fid in "${coverage_files[@]}"; do
            fname=$(dx describe "$fid" --json | python3 -c "import json,sys; print(json.load(sys.stdin)['name'])")
            dx download "$fid" -o "cnn_files/${fname}"
        done
    fi

    N_FILES=$(ls cnn_files/*.cnn 2>/dev/null | wc -l)
    echo "[inputs] Coverage files found: ${N_FILES}"
    [ "${N_FILES}" -ge 10 ] || { echo "ERROR: too few coverage files (${N_FILES}); expect ≥10"; exit 1; }

    ls cnn_files/*.cnn | head -5
    echo "..."

    # ── Build reference ─────────────────────────────────────────────────────────────
    # If a FASTA is provided, annotate intervals with GC content and
    # RepeatMasker fraction for GC-bias correction in downstream fix steps.
    # pyfaidx (bundled with CNVkit) handles bgzf-compressed FASTA natively
    # when both the .fai and .gzi index files are present alongside.
    FASTA_ARG=""
    # Guard: if any FASTA index is provided without the FASTA itself, fail early
    if ([ -n "${fasta_fai:-}" ] || [ -n "${fasta_gzi:-}" ]) && [ -z "${fasta:-}" ]; then
        echo "ERROR: fasta_fai/fasta_gzi supplied but fasta is missing — all three must be provided together" >&2
        exit 1
    fi

    if [ -n "${fasta:-}" ]; then
        [ -n "${fasta_fai:-}" ] || { echo "ERROR: fasta supplied but fasta_fai is missing" >&2; exit 1; }
        [ -n "${fasta_gzi:-}" ] || { echo "ERROR: fasta supplied but fasta_gzi is missing" >&2; exit 1; }
        echo "[reference] Downloading FASTA for GC correction..."
        dx download "${fasta}"     -o ref.fasta.gz
        dx download "${fasta_fai}" -o ref.fasta.gz.fai
        dx download "${fasta_gzi}" -o ref.fasta.gz.gzi
        FASTA_ARG="--fasta ref.fasta.gz"
        echo "[reference] FASTA ready; GC correction enabled"
    else
        echo "[reference] No FASTA supplied; skipping GC/repeat annotation"
    fi

    echo "[reference] Building PoN from ${N_FILES} samples..."
    run_cnvkit reference cnn_files/*.cnn \
        --output "${pon_name:-cgp_cnvkit_reference}.cnn" \
        ${FASTA_ARG}

    export REF_FILE="${pon_name:-cgp_cnvkit_reference}.cnn"
    N_INTERVALS=$(wc -l < "${REF_FILE}")
    echo "[reference] Reference intervals (incl. header): ${N_INTERVALS}"
    [ "${N_INTERVALS}" -gt 1000 ] || { echo "ERROR: reference too small"; exit 1; }

    # ── Summary stats ─────────────────────────────────────────────────────────
    python3 - << 'PYEOF'
import csv, statistics, os

ref_file = os.environ["REF_FILE"]
log2_vals = []
depths = []

with open(ref_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        try:
            log2_vals.append(float(row['log2']))
            depths.append(float(row.get('depth', 0)))
        except (ValueError, KeyError):
            pass

with open('reference_stats.tsv', 'w') as out:
    out.write('metric\tvalue\n')
    out.write(f'n_intervals\t{len(log2_vals)}\n')
    out.write(f'log2_mean\t{statistics.mean(log2_vals):.4f}\n')
    out.write(f'log2_stdev\t{statistics.stdev(log2_vals):.4f}\n')
    out.write(f'log2_median\t{statistics.median(log2_vals):.4f}\n')
    out.write(f'depth_median\t{statistics.median(depths):.1f}\n')
    out.write(f'n_samples\t{len([f for f in os.listdir("cnn_files") if f.endswith(".cnn")])}\n')

print("Reference stats:")
with open('reference_stats.tsv') as f:
    print(f.read())
PYEOF

    # ── Upload outputs ────────────────────────────────────────────────────────
    cn_reference=$(dx upload "${REF_FILE}" --brief)
    reference_stats=$(dx upload reference_stats.tsv --brief)

    dx-jobutil-add-output cn_reference    "${cn_reference}"    --class=file
    dx-jobutil-add-output reference_stats "${reference_stats}" --class=file

    echo "=== PoN build complete: ${N_FILES} samples, ${N_INTERVALS} intervals ==="
}
