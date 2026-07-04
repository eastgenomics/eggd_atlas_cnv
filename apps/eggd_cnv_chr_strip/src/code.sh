#!/bin/bash
# eggd_cnv_chr_strip v0.1.0 (final stage) — Ensembl-named copies of CNV text files.
# Writes *.nochr.* copies with the chr prefix stripped from the CHROMOSOME COLUMN ONLY
# (chrM/chrMT -> MT). Never modifies the originals (CNVkit's own outputs; PURPLE's inside
# purple_tar). Delegates the transform to the bundled pure-Python atlas_helpers.chr_strip.
set -euo pipefail

ATLAS=/home/dnanexus/atlas

main() {
    echo "====================================================="
    echo " eggd_cnv_chr_strip: chr -> Ensembl CNV copies"
    echo " Sample  : ${sample_id}"
    echo "====================================================="

    dx-download-all-inputs --parallel
    case "${sample_id}" in
        *[!A-Za-z0-9._-]* | "" | .* ) echo "ERROR: unsafe sample_id" >&2; exit 1 ;;
    esac

    strip() {  # $1 = input field name, $2 = output filename, $3 = dx output field
        local f
        f=$(ls ~/in/"$1"/* 2>/dev/null || true)
        [ -n "$f" ] || { echo "  (skip $1 — not supplied)"; return 0; }
        echo "  stripping $1 -> $2"
        python3 -c "import sys;sys.path.insert(0,'$ATLAS');from chr_strip import strip_file;strip_file(sys.argv[1],sys.argv[2])" "$f" "$2"
        [ -s "$2" ] || { echo "ERROR: strip produced no output for $1" >&2; exit 1; }
        dx-jobutil-add-output "$3" "$(dx upload "$2" --brief)" --class=file
    }

    strip cnvkit_cnr          "${sample_id}.nochr.cnr"                    cnvkit_cnr_nochr
    strip cnvkit_cns          "${sample_id}.nochr.cns"                    cnvkit_cns_nochr
    strip cnvkit_call_cns     "${sample_id}.nochr.call.cns"               cnvkit_call_cns_nochr
    strip cnvkit_genemetrics  "${sample_id}.nochr.genemetrics.tsv"        cnvkit_genemetrics_nochr
    strip purple_cnv_somatic  "${sample_id}.purple.cnv.somatic.nochr.tsv" purple_cnv_somatic_nochr
    strip purple_cnv_gene     "${sample_id}.purple.cnv.gene.nochr.tsv"    purple_cnv_gene_nochr

    echo "====================================================="
    echo " eggd_cnv_chr_strip DONE: ${sample_id}"
    echo "====================================================="
}
