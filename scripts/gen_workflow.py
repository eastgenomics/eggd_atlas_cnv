#!/usr/bin/env python3
"""Generate dxworkflow.json for eggd_atlas_cnv (7 stages). Run from repo root.
Fixed reference file IDs are from REFERENCE.md §2."""
import json
import os

PROJ_REF = "project-Fkb6Gkj433GVVvj73J7x8KbV"  # canonical reference project

# Stage executables. Prefer exact app IDs from scripts/app_ids.json (built this session)
# because by-name resolution (app-<name>) resolves to any PUBLISHED app of that name,
# which for the reused cnvkit apps can shadow our dev builds with a different interface.
# Fall back to app-<name> if the map is absent.
_IDS_PATH = os.path.join(os.path.dirname(__file__), "app_ids.json")
APP_ID = json.load(open(_IDS_PATH)) if os.path.exists(_IDS_PATH) else {}


def execu(app_name):
    return APP_ID.get(app_name, f"app-{app_name}")

# Fixed reference files (REFERENCE §2). Objects without an explicit project use bare id.
F = {
    "AMBER_JAR": {"id": "file-J9936z04p12z3FYjFYb448pP"},
    "GERMLINE_SITES": {"id": "file-J88xxvQ4QyV82JKY5PkqXGbv"},
    "COBALT_JAR": {"id": "file-J893p9Q470j4zY3zzpVBjP11"},
    "NORM_FILE": {"id": "file-J89183j4zxxZ5v7Gv5F9BK06"},
    "DIPLOID_REGIONS": {"id": "file-J88xxvQ4QyV7vyPx294b4Z23"},
    "GC_PROFILE": {"id": "file-J88xxvQ4QyVPb8K6VFqX1FKB"},
    "COBALT_REF_FASTA": {"project": PROJ_REF, "id": "file-GjPxXq84qv8xz3ZV1jFq6z3g"},
    "COBALT_REF_FAI": {"id": "file-GjPxp3Q4qv8Vk76xFVjfyPgJ"},
    "SAGE_JAR": {"id": "file-J8F1bFj4gPFvK81Z330Yfk3j"},
    "SAGE_REF_FASTA": {"id": "file-G5xBZvj4yPz1xZfqKKjkG8xQ"},
    "HOTSPOTS_VCF": {"id": "file-J8F1kQQ4gPFbZQZ1BG7KY62z"},
    "HOTSPOTS_TBI": {"id": "file-J8F1kV04gPFz7kZ45XfQPp4K"},
    "PANEL_BED": {"id": "file-J8F1kv84gPFfQKqYbZVYQqv8"},
    "HC_BED": {"id": "file-J8F1kVj4gPFkBK29yj6P4BgQ"},
    "SAGE_PON": {"id": "file-J8F1kX84gPFpJ2kvz13KJ5kP"},
    "ENSEMBL_DATA": {"id": "file-J892FbQ493ZPf7pZzpjX05P6"},
    "PURPLE_JAR": {"id": "file-J893gpQ470jJPF18gjJ1zg7Y"},
    "PURPLE_TARGET_BED": {"id": "file-J88gVF84Y8X123K6JX8jB8Z5"},
    "PURPLE_REF_FASTA": {"id": "file-Gb757784XGyY3FPvkPQ74K9z"},
    "PURPLE_REF_FAI": {"id": "file-Gb7578Q4XGyQ8xfvyxJkBgx7"},
    "CNVKIT_BED": {"id": "file-J8F94G845FG5ZG55xFzF1fP6"},
}


def link_file(key):
    v = F[key]
    if "project" in v:
        return {"$dnanexus_link": {"project": v["project"], "id": v["id"]}}
    return {"$dnanexus_link": v["id"]}  # bare-string form for project-less file refs


def wf(name):  # workflow-level input
    return {"$dnanexus_link": {"workflowInputField": name}}


def stg(stage, field):  # stage output link (JBOR)
    return {"$dnanexus_link": {"stage": stage, "outputField": field}}


stages = [
    {
        "id": "amber",
        "executable": execu("eggd_cgp-amber"),
        "input": {
            "tumour_bam": wf("input_bam"),
            "tumour_bai": wf("input_bai"),
            "sample_id": wf("sample_id"),
            "amber_jar": link_file("AMBER_JAR"),
            "germline_sites": link_file("GERMLINE_SITES"),
        },
    },
    {
        "id": "cobalt",
        "executable": execu("eggd_cgp-cobalt"),
        "input": {
            "tumour_bam": wf("input_bam"),
            "tumour_bai": wf("input_bai"),
            "sample_id": wf("sample_id"),
            "cobalt_jar": link_file("COBALT_JAR"),
            "norm_file": link_file("NORM_FILE"),
            "diploid_regions": link_file("DIPLOID_REGIONS"),
            "gc_profile": link_file("GC_PROFILE"),
            "ref_fasta": link_file("COBALT_REF_FASTA"),
            "ref_fai": link_file("COBALT_REF_FAI"),
        },
    },
    {
        "id": "sage",
        "executable": execu("eggd_cgp-sage"),
        "input": {
            "tumour_bam": wf("input_bam"),
            "tumour_bai": wf("input_bai"),
            "sample_id": wf("sample_id"),
            "sage_jar": link_file("SAGE_JAR"),
            "ref_fasta": link_file("SAGE_REF_FASTA"),
            "hotspots_vcf": link_file("HOTSPOTS_VCF"),
            "hotspots_tbi": link_file("HOTSPOTS_TBI"),
            "panel_bed": link_file("PANEL_BED"),
            "hc_bed": link_file("HC_BED"),
            "pon_file": link_file("SAGE_PON"),
            "ensembl_data": link_file("ENSEMBL_DATA"),
        },
    },
    {
        "id": "purple",
        "executable": execu("eggd_cgp-purple"),
        "input": {
            "sample_id": wf("sample_id"),
            "purple_jar": link_file("PURPLE_JAR"),
            "amber_tar": stg("amber", "amber_tar"),
            "cobalt_tar": stg("cobalt", "cobalt_tar"),
            "somatic_vcf": stg("sage", "somatic_vcf"),
            "somatic_vcf_tbi": stg("sage", "somatic_vcf_tbi"),
            "gc_profile": link_file("GC_PROFILE"),
            "target_regions_bed": link_file("PURPLE_TARGET_BED"),
            "ref_fasta": link_file("PURPLE_REF_FASTA"),
            "ref_fai": link_file("PURPLE_REF_FAI"),
            "ensembl_data": link_file("ENSEMBL_DATA"),
            "max_ploidy": wf("max_ploidy"),
            "ploidy_cap_purity_threshold": wf("ploidy_cap_purity_threshold"),
            "ploidy_cap_value": wf("ploidy_cap_value"),
        },
    },
    {
        "id": "cnvkit_batch",
        "executable": execu("eggd_cgp-cnvkit-batch"),
        "input": {
            "tumour_bam": wf("input_bam"),
            "tumour_bai": wf("input_bai"),
            "sample_id": wf("sample_id"),
            "cn_reference": wf("cnvkit_cn_reference"),
            "baits": link_file("CNVKIT_BED"),
            "purity": stg("purple", "purity"),
            "ploidy": stg("purple", "ploidy"),
            "sample_sex": stg("purple", "sample_sex"),
        },
    },
    {
        "id": "cnv_chr_strip",
        "executable": execu("eggd_cnv_chr_strip"),
        "input": {
            "sample_id": wf("sample_id"),
            "cnvkit_cnr": stg("cnvkit_batch", "copy_ratios"),
            "cnvkit_cns": stg("cnvkit_batch", "segments"),
            "cnvkit_call_cns": stg("cnvkit_batch", "call_segments"),
            "cnvkit_genemetrics": stg("cnvkit_batch", "genemetrics"),
            "purple_cnv_somatic": stg("purple", "cnv_somatic_tsv"),
            "purple_cnv_gene": stg("purple", "cnv_gene_tsv"),
        },
    },
    {
        "id": "purple_plotter",
        "executable": execu("eggd_purple_plotter"),
        "input": {
            "sample_id": wf("sample_id"),
            "purple_tar": stg("purple", "purple_tar"),
            "cnvkit_cnr": stg("cnvkit_batch", "copy_ratios"),
            "cnvkit_call_cns": stg("cnvkit_batch", "call_segments"),
            "cnvkit_genemetrics": stg("cnvkit_batch", "genemetrics"),
        },
    },
]

inputs = [
    # These are chr-prefixed BAM/BAI outputs of the separately-run eggd_chr_prefix app.
    {"name": "input_bam", "class": "file"},
    {"name": "input_bai", "class": "file"},
    {"name": "sample_id", "class": "string"},
    {"name": "cnvkit_cn_reference", "class": "file"},
    {"name": "max_ploidy", "class": "int", "optional": True},
    {"name": "ploidy_cap_purity_threshold", "class": "float", "optional": True},
    {"name": "ploidy_cap_value", "class": "int", "optional": True, "default": 2},
]


def out(name, cls, stage, field):
    return {"name": name, "class": cls, "outputSource": stg(stage, field)}


outputs = [
    out("igv_html", "file", "purple_plotter", "igv_html"),
    out("purple_cnv_somatic_nochr", "file", "cnv_chr_strip", "purple_cnv_somatic_nochr"),
    out("purple_cnv_gene_nochr", "file", "cnv_chr_strip", "purple_cnv_gene_nochr"),
    out("cnvkit_call_cns_nochr", "file", "cnv_chr_strip", "cnvkit_call_cns_nochr"),
    out("cnvkit_genemetrics_nochr", "file", "cnv_chr_strip", "cnvkit_genemetrics_nochr"),
    out("cnvkit_cnr_nochr", "file", "cnv_chr_strip", "cnvkit_cnr_nochr"),
    out("purple_tar", "file", "purple", "purple_tar"),
    out("purity", "float", "purple", "purity"),
    out("ploidy", "int", "purple", "ploidy"),
]

workflow = {
    "name": "eggd_atlas_cnv",
    "title": "eggd_atlas_cnv",
    "summary": "Per-sample somatic CNV workflow: chr-prefixed BAM/BAI, AMBER+COBALT+SAGE -> PURPLE, CNVkit batch, IGV plotter, chr-strip",
    "dxapi": "1.0.0",
    "version": "0.1.0",
    "inputs": inputs,
    "outputs": outputs,
    "stages": stages,
    "regionalOptions": {"aws:eu-central-1": {}},
}

with open("dxworkflow.json", "w") as fh:
    json.dump(workflow, fh, indent=2)
    fh.write("\n")
print("wrote dxworkflow.json")
