#!/bin/bash
# eggd_cgp-purple v1.0.0 — PURPLE 4.4-beta.12 tumour-only purity/ploidy/CN
# Extended from cgp-purple applet with:
#   - optional max-ploidy cap (STATIC / CONDITIONAL two-pass, via bundled ploidy_gate.py)
#   - typed scalar outputs purity(float)/ploidy(int)/sample_sex(string) for stage linking
#   - standalone cnv_somatic_tsv / cnv_gene_tsv file outputs (were only inside purple_tar)
# PURPLE tool flags are FROZEN; the only added flag is -max_ploidy from the cap decision.
set -eo pipefail

ATLAS=/home/dnanexus/atlas

main() {
    echo "====================================================="
    echo " eggd_cgp-purple: PURPLE purity/ploidy/CN"
    echo " Sample  : ${sample_id}"
    echo "====================================================="

    # ── 0. Sanitise sample_id before any path/rm use ────────────────────────
    case "${sample_id}" in
        *[!A-Za-z0-9._-]* | "" | .* ) echo "ERROR: unsafe sample_id" >&2; exit 1 ;;
    esac

    echo "[setup] Verifying deps (execDepends: java/samtools/tabix/jq/R/circos)..."
    java    -version  2>&1 | head -1
    samtools --version 2>&1 | head -1
    bgzip   --version 2>&1 | head -1
    jq      --version 2>&1 | head -1
    Rscript --version 2>&1 | head -1
    circos  --version 2>&1 | head -1 || echo "circos: $(which circos 2>/dev/null || echo not found)"

    echo "[1/7] Downloading inputs..."
    dx download "${purple_jar}"        -o purple.jar
    dx download "${amber_tar}"         -o amber.tar.gz
    dx download "${cobalt_tar}"        -o cobalt.tar.gz
    dx download "${gc_profile}"        -o GC_profile.1000bp.38.cnp
    dx download "${target_regions_bed}" -o target_regions.bed
    dx download "${ref_fasta}"         -o ref.fasta.gz
    dx download "${ref_fai}"           -o ref.fasta.gz.fai
    dx download "${ensembl_data}"      -o ensembl_data.tar.gz

    SOMATIC_ARG=""
    if [[ -n "${somatic_vcf:-}" ]]; then
        dx download "${somatic_vcf}"     -o somatic.vcf.gz
        dx download "${somatic_vcf_tbi}" -o somatic.vcf.gz.tbi
        SOMATIC_ARG="-somatic_vcf somatic.vcf.gz"
        echo "  Somatic VCF provided — will anchor purity"
    fi

    # ── 2. Prepare inputs ──────────────────────────────────────────────────
    echo "[2/7] Preparing inputs..."
    tar --no-same-owner -xzf amber.tar.gz
    tar --no-same-owner -xzf cobalt.tar.gz
    AMBER_DIR=$(find . -maxdepth 2 -name "*.amber.baf.tsv.gz" -printf "%h\n" | head -1)
    COBALT_DIR=$(find . -maxdepth 2 -name "*.cobalt.ratio.tsv.gz" -printf "%h\n" | head -1)
    [[ -n "${AMBER_DIR}"  ]] || { echo "ERROR: AMBER dir not found in tar"; exit 1; }
    [[ -n "${COBALT_DIR}" ]] || { echo "ERROR: COBALT dir not found in tar"; exit 1; }
    echo "  AMBER dir:  ${AMBER_DIR}"
    echo "  COBALT dir: ${COBALT_DIR}"

    echo "Decompressing reference FASTA..."
    bgzip -d -c ref.fasta.gz > ref.fasta
    samtools faidx ref.fasta

    tar --no-same-owner -xzf ensembl_data.tar.gz
    ENSEMBL_DIR=$(find . -maxdepth 2 -name "ensembl_gene_data.csv" -printf "%h\n" | head -1)
    [[ -n "${ENSEMBL_DIR}" ]] || { echo "ERROR: ensembl_gene_data.csv not found"; exit 1; }
    echo "  Ensembl dir: ${ENSEMBL_DIR}"

    WORK="${sample_id}"           # PURPLE output dir; names its files ${sample_id}.purple.*

    run_purple() {                # "$@" = extra PURPLE args (may be empty)
        java -Xmx10G -jar purple.jar \
            -tumor              "${sample_id}" \
            -amber              "${AMBER_DIR}" \
            -cobalt             "${COBALT_DIR}" \
            -target_regions_bed target_regions.bed \
            -gc_profile         GC_profile.1000bp.38.cnp \
            -ensembl_data_dir   "${ENSEMBL_DIR}" \
            -ref_genome         ref.fasta \
            -ref_genome_version 38 \
            ${SOMATIC_ARG} \
            "$@" \
            -output_dir         "${WORK}/"
    }

    # ── 3. Resolve ploidy-cap decision (shell-safe, delimiter-free) ─────────
    echo "[3/7] Resolving ploidy-cap decision..."
    EVAL=$(max_ploidy="${max_ploidy:-}" \
           ploidy_cap_purity_threshold="${ploidy_cap_purity_threshold:-}" \
           ploidy_cap_value="${ploidy_cap_value:-2}" \
           python3 - "$ATLAS" <<'PY'
import sys, os, json
sys.path.insert(0, sys.argv[1])
from ploidy_gate import decide
d = decide(
    max_ploidy=(int(os.environ["max_ploidy"]) if os.environ.get("max_ploidy") else None),
    threshold=(float(os.environ["ploidy_cap_purity_threshold"]) if os.environ.get("ploidy_cap_purity_threshold") else None),
    cap_value=int(os.environ.get("ploidy_cap_value", "2")),
)
print("FIRST_ARGS_JSON=" + json.dumps(json.dumps(d.first_pass_args)))  # quoted JSON string
print("CONDITIONAL=" + ("yes" if d.conditional else "no"))
print("MODE=" + d.mode.value)
PY
)
eval "$EVAL"                                   # sets FIRST_ARGS_JSON + CONDITIONAL + MODE
echo "  mode=${MODE}  conditional=${CONDITIONAL}  first_pass_args=${FIRST_ARGS_JSON}"
mapfile -t FIRST_ARGS < <(jq -r '.[]' <<<"$FIRST_ARGS_JSON")

# ── 4. PURPLE pass 1 ────────────────────────────────────────────────────
echo "[4/7] Running PURPLE (pass 1)..."
rm -rf "${WORK}"; mkdir -p "${WORK}"
run_purple "${FIRST_ARGS[@]}"
PTSV="${WORK}/${sample_id}.purple.purity.tsv"
RANGE_TSV="${WORK}/${sample_id}.purple.purity.range.tsv"
[[ -s "${PTSV}" ]]      || { echo "ERROR: PURPLE purity TSV missing after pass 1"; exit 1; }
[[ -s "${RANGE_TSV}" ]] || { echo "ERROR: PURPLE range TSV missing after pass 1"; exit 1; }

# ── 5. Conditional re-run (at most one extra pass) ──────────────────────
if [ "$CONDITIONAL" = "yes" ]; then
    PURITY=$(python3 -c "import sys;sys.path.insert(0,'$ATLAS');from purity import read_purity_ploidy as r;print(r('$PTSV').purity)")
    NEED=$(python3 -c "print('yes' if float('$PURITY') < float('${ploidy_cap_purity_threshold}') else 'no')")
    echo "  pass-1 purity=${PURITY} threshold=${ploidy_cap_purity_threshold} -> rerun=${NEED}"
    if [ "$NEED" = "yes" ]; then
        echo "[4b/7] Purity below threshold — re-running PURPLE with -max_ploidy ${ploidy_cap_value:-2} (pass 2)..."
        rm -rf "${WORK}"; mkdir -p "${WORK}"
        run_purple -max_ploidy "${ploidy_cap_value:-2}"
        [[ -s "${PTSV}" ]] || { echo "ERROR: PURPLE purity TSV missing after pass 2"; exit 1; }
    fi
fi

# ── 6. Verify + parse final purity.tsv ──────────────────────────────────
echo "[5/7] Final purity result:"
cat "${PTSV}"
ls -lh "${WORK}/"

eval "$(python3 - "$ATLAS" "$PTSV" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from purity import read_purity_ploidy
f = read_purity_ploidy(sys.argv[2])
print(f'FPUR={f.purity}; FPLO={f.ploidy_int()}; FSEX={f.sample_sex or "unknown"}')
PY
)"
echo "  emitting purity=${FPUR} ploidy=${FPLO} sample_sex=${FSEX}"

# ── 6b. Collect charts (if generated) ────────────────────────────────────
PLOT_FILES=$(find "${WORK}/" -name "*.png" 2>/dev/null | sort)
if [ -n "${PLOT_FILES}" ]; then
    tar --no-same-owner -czf "${sample_id}.purple.plots.tar.gz" ${PLOT_FILES}
    HAVE_PLOTS=true
else
    echo "  No chart PNGs generated (non-fatal)"
    HAVE_PLOTS=false
fi

# ── 7. Upload outputs ───────────────────────────────────────────────────
echo "[6/7] Uploading..."
tar --no-same-owner -czf "${sample_id}.purple.tar.gz" "${WORK}/"

dx-jobutil-add-output purple_tar       "$(dx upload "${sample_id}.purple.tar.gz" --brief)" --class=file
dx-jobutil-add-output purity_tsv       "$(dx upload "${PTSV}"      --brief)" --class=file
dx-jobutil-add-output purity_range_tsv "$(dx upload "${RANGE_TSV}" --brief)" --class=file
if [ "${HAVE_PLOTS}" = "true" ]; then
    dx-jobutil-add-output plots_tar "$(dx upload "${sample_id}.purple.plots.tar.gz" --brief)" --class=file
fi

# Standalone CNV call TSVs — always exist (header-only fallback if PURPLE omitted them).
# Headers below match PURPLE 4.4's *.purple.cnv.somatic.tsv / *.purple.cnv.gene.tsv columns.
CNV_SOMATIC="${WORK}/${sample_id}.purple.cnv.somatic.tsv"
CNV_GENE="${WORK}/${sample_id}.purple.cnv.gene.tsv"
[ -f "${CNV_SOMATIC}" ] || printf 'chromosome\tstart\tend\tcopyNumber\tbafCount\tobservedBAF\tbaf\tsegmentStartSupport\tsegmentEndSupport\tmethod\tdepthWindowCount\tgcContent\tminStart\tmaxStart\tminorAlleleCopyNumber\tmajorAlleleCopyNumber\n' > "${CNV_SOMATIC}"
[ -f "${CNV_GENE}" ]    || printf 'chromosome\tstart\tend\tgene\tminCopyNumber\tmaxCopyNumber\tsomaticRegions\ttranscriptId\tisCanonical\tchromosomeBand\tminRegions\tminRegionStart\tminRegionEnd\tminRegionStartSupport\tminRegionEndSupport\tminRegionMethod\tminMinorAlleleCopyNumber\tdepthWindowCount\n' > "${CNV_GENE}"
dx-jobutil-add-output cnv_somatic_tsv "$(dx upload "${CNV_SOMATIC}" --brief)" --class=file
dx-jobutil-add-output cnv_gene_tsv    "$(dx upload "${CNV_GENE}"    --brief)" --class=file

# Typed scalar outputs (the only way purity/ploidy/sex reach downstream stages).
dx-jobutil-add-output purity     "${FPUR}" --class=float
dx-jobutil-add-output ploidy     "${FPLO}" --class=int
dx-jobutil-add-output sample_sex "${FSEX}" --class=string

echo "[7/7] Done."
echo "====================================================="
echo " eggd_cgp-purple DONE: ${sample_id}"
echo "   purity=${FPUR}  ploidy=${FPLO}  sample_sex=${FSEX}"
echo "====================================================="
}
